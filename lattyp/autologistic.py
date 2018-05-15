# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import norm, gamma
from scipy.special import gammaln
import random
import sys
from collections import defaultdict

from rand_utils import rand_partition_log

class CategoricalAutologistic(object):
    S_X = 1
    S_V = 2
    S_H = 3
    S_A = 4

    def __init__(self, vec, size,
                 vnet=None, hnet=None,
                 mvs=None,
                 only_alphas=False,
                 drop_vs=False,
                 drop_hs=False,
                 norm_sigma = 5.0,
                 gamma_shape = 1.0,
                 gamma_scale = 0.001,
    ):
        self.vec = vec
        self.L = self.vec.size
        self.size = size

        self.vnet = vnet
        self.hnet = hnet

        self.only_alphas = only_alphas
        self.drop_vs = drop_vs
        self.drop_hs = drop_hs
        assert(mvs is None or vec.shape == mvs.shape)
        self.mvs = mvs # missing values; (i,p) => bool (True: missing value)
        self.mv_list = []
        bincount = np.bincount(self.vec, minlength=self.size).astype(np.float32)
        if self.mvs is not None:
            bincount -= np.bincount(self.vec[self.mvs], minlength=self.size)
            for l in range(self.L):
                if self.mvs[l]:
                    self.mv_list.append(l)

        self.norm_sigma = norm_sigma
        self.gamma_shape = gamma_shape
        self.gamma_scale = gamma_scale

        a = np.log(bincount / bincount.sum() + 0.001)
        a -= a.mean()
        self.alphas = a
        # self.alphas = 0.5 * np.random.normal(loc=0.0, scale=self.norm_sigma, size=self.size)
        if not (self.only_alphas or self.drop_vs):
            self.v = 0.0001
        else:
            self.v = 0.0
        if not (self.only_alphas or self.drop_hs):
            self.h = 0.0001
        else:
            self.h = 0.0
        self.init_tasks()

    def init_tasks(self, a_repeat=1, sample_w=True):
        self.tasks = []
        for l in self.mv_list:
            self.tasks.append((self.S_X, l))
        for a in range(a_repeat):
            if not self.only_alphas:
                if not self.drop_vs:
                    self.tasks.append((self.S_V, None))
                if not self.drop_hs:
                    self.tasks.append((self.S_H, None))
            for k in range(self.size):
                self.tasks.append((self.S_A, k))

    def sample(self):
        c_x = [0, 0]
        c_v = [0, 0]
        c_h = [0, 0]
        c_a = [0, 0]
        random.shuffle(self.tasks)
        for t_type, t_val in self.tasks:
            if t_type == self.S_X:
                l = t_val
                changed = self.sample_x(l)
                c_x[changed] += 1
            elif t_type == self.S_V:
                changed = self.sample_autologistic(t_type, t_val)
                c_v[changed] += 1
            elif t_type == self.S_H:
                changed = self.sample_autologistic(t_type, t_val)
                c_h[changed] += 1
            elif t_type == self.S_A:
                changed = self.sample_autologistic(t_type, t_val)
                c_a[changed] += 1
            else:
                raise NotImplementedError
        if sum(c_x) > 0:
            sys.stderr.write("\tx_cat\t{}\n".format(float(c_x[1]) / sum(c_x)))
        if not self.only_alphas:
            if sum(c_v) > 0:
                sys.stderr.write("\tz_v\t{}\n".format(float(c_v[1]) / sum(c_v)))
            if sum(c_h) > 0:
                sys.stderr.write("\tz_h\t{}\n".format(float(c_h[1]) / sum(c_h)))
        if sum(c_a) > 0:
            sys.stderr.write("\tz_a\t{}\n".format(float(c_a[1]) / sum(c_a)))
        sys.stderr.write("\tv\t{}\th\t{}\ta\t{}\n".format(
            self.v if not self.only_alphas and not self.drop_vs else 0.0,
            self.h if not self.only_alphas and not self.drop_hs else 0.0,
            self.alphas,
        ))

    def sample_x(self, l):
        assert(self.mvs is not None and self.mvs[l])
        x_old = self.vec[l]
        logprobs = np.zeros(self.size, dtype=np.float32)
        if not self.only_alphas:
            if not self.drop_vs:
                logprobs += self.v * np.bincount(self.vec[self.vnet[l]], minlength=self.size)
            if not self.drop_hs:
                logprobs += self.h * np.bincount(self.vec[self.hnet[l]], minlength=self.size)
        logprobs += self.alphas
        self.vec[l] = rand_partition_log(logprobs)
        return False if self.vec[l] == x_old else True

    def sample_autologistic(self, t_type, k):
        logr = 0.0
        if t_type == self.S_A:
            oldval = self.alphas[k]
            newval = np.random.normal(loc=oldval, scale=0.01)
            # P(theta') / P(theta)
            logr += (oldval ** 2 - newval ** 2) / (2.0 * self.norm_sigma * self.norm_sigma)
            alphas = self.alphas.copy()
            alphas[k] = newval
            v, h, a = self.v, self.h, alphas
        else:
            assert(not self.only_alphas)
            assert(not (t_type == self.S_V and self.drop_vs))
            assert(not (t_type == self.S_H and self.drop_hs))
            if t_type == self.S_V:
                oldval = self.v
            else:
                oldval = self.h
            P_SIGMA = 0.5
            rate = np.random.lognormal(mean=0.0, sigma=P_SIGMA)
            irate = 1.0 / rate
            newval = rate * oldval
            lograte = np.log(rate)
            logirate = np.log(irate)
            # P(theta') / P(theta)
            # logr += gamma.logpdf(newval, self.gamma_shape, scale=self.gamma_scale) \
            #         - gamma.logpdf(oldval, self.gamma_shape, scale=self.gamma_scale)
            logr += (self.gamma_shape - 1.0) * (np.log(newval) - np.log(oldval)) \
                    - (newval - oldval) / self.gamma_scale
            # q(theta|theta', x) / q(theta'|theta, x)
            logr += (lograte * lograte - logirate * logirate) / (2.0 * P_SIGMA * P_SIGMA) + lograte - logirate
            if t_type == self.S_V:
                v, h, a = newval, self.h, self.alphas
                net = self.vnet
            else:
                v, h, a = self.v, newval, self.alphas
                net = self.hnet
        vec = self.vec.copy()
        llist = np.arange(self.L)
        np.random.shuffle(llist)
        for l in llist:
            logprobs = np.zeros(self.size, dtype=np.float32)
            if not self.only_alphas:
                if not self.drop_vs:
                    logprobs += v * np.bincount(self.vec[self.vnet[l]], minlength=self.size)
                if not self.drop_hs:
                    logprobs += h * np.bincount(self.vec[self.hnet[l]], minlength=self.size)
            logprobs += a
            vec[l] = rand_partition_log(logprobs)
        if t_type == self.S_A:
            logr += (oldval - newval) * ((vec == k).sum() - (self.vec == k).sum())
            if logr >= 0 or np.log(np.random.rand()) < logr:
                # accept
                self.alphas[k] = newval
                return True
            else:
                return False
        else:
            oldsum = self._neighbor_sum(self.vec, net)
            newsum = self._neighbor_sum(vec, net)
            logr += (oldval - newval) * (newsum - oldsum)
            if logr >= 0 or np.log(np.random.rand()) < logr:
                # accept
                if t_type == self.S_V:
                    self.v = newval
                else:
                    self.h = newval
                return True
            else:
                return False

    def _neighbor_sum(self, vec, net):
        s = 0
        for l in range(self.L):
            s += (vec[net[l]] == vec[l]).sum()
        assert(s % 2 == 0)
        return s / 2
