# -*- coding: utf-8 -*-
import numpy as np
from scipy.stats import norm, gamma
from scipy.special import gammaln
import random
import sys
from collections import defaultdict

from rand_utils import rand_partition_log
from hmc import hmc

class MatrixDecompositionAutologistic(object):
    S_X_CAT = 1
    S_X_BIN = 2
    S_X_CNT = 3
    S_Z = 4
    S_W_MH = 5
    S_W_HMC = 6
    S_Z_V = 7
    S_Z_H = 8
    S_Z_A = 9

    # backward compatibility
    hmc_l = 10
    hmc_epsilon = 0.05

    def __init__(self, mat, flist,
                 sigma=1.0,
                 vnet=None, hnet=None,
                 K=50, mvs=None,
                 bias=False,
                 only_alphas=False,
                 drop_vs=False,
                 drop_hs=False,
                 norm_sigma = 5.0,
                 gamma_shape = 1.0,
                 gamma_scale = 0.001,
                 hmc_l = 10,
                 hmc_epsilon = 0.05,
    ):
        self.mat = mat # X: L x P matrix
        self.flist = flist

        self.hmc_l = hmc_l
        self.hmc_epsilon = hmc_epsilon

        # L: # of langs
        # P: # of surface features
        # M: # of linearlized weight elements
        # K: # of latent parameters
        # l: current idx of langs
        # p: current idx of surface features
        # T: size of the current feature
        # j: linearlized idx of the current feature
        self.L, self.P = self.mat.shape
        self.M = sum(map(lambda x: x["size"], self.flist))
        self.K = K

        self.j2pt = np.empty((self.M, 2), dtype=np.int32)
        self.p2jT = np.empty((self.P, 2), dtype=np.int32)
        binsize = 0
        for fstruct in self.flist:
            self.p2jT[fstruct["fid"]] = [binsize, fstruct["size"]]
            for t in range(fstruct["size"]):
                self.j2pt[binsize+t] = [fstruct["fid"], t]
            binsize += fstruct["size"]
        
        self.vnet = vnet
        self.hnet = hnet

        self.bias = bias
        self.only_alphas = only_alphas
        self.drop_vs = drop_vs
        self.drop_hs = drop_hs
        assert(mvs is None or mat.shape == mvs.shape)
        self.mvs = mvs # missing values; (i,p) => bool (True: missing value)
        self.mv_list = []
        if self.mvs is not None:
            for l in range(self.L):
                for p in range(self.P):
                    if self.mvs[l,p]:
                        self.mv_list.append((l,p))
        self.norm_sigma = norm_sigma
        self.gamma_shape = gamma_shape
        self.gamma_scale = gamma_scale
        self.alphas = 0.5 * np.random.normal(loc=0.0, scale=self.norm_sigma, size=self.K)
        if not (self.only_alphas or self.drop_vs):
            self.vks = 0.0001 * np.ones(self.K, dtype=np.float32)
        else:
            self.vks = np.zeros(self.K, dtype=np.float32)
        if not (self.only_alphas or self.drop_hs):
            self.hks = 0.0001 * np.ones(self.K, dtype=np.float32)
        else:
            self.hks = np.zeros(self.K, dtype=np.float32)
        self.zmat = np.zeros((self.K, self.L), dtype=np.bool_)
        for k, alpha in enumerate(self.alphas):
            thres = 1.0 / (1.0 + np.exp(-alpha))
            self.zmat[k] = (np.random.rand(self.L) < thres)
        self.sigma = sigma # Normal
        self.wmat = 0.1 * np.random.standard_t(df=self.sigma, size=(self.K, self.M))
        self.theta_tilde = np.zeros((self.L, self.M), dtype=np.float32)
        self.theta = np.ones((self.L, self.M), dtype=np.float32)
        self.calc_theta()
        self.init_tasks()

    def init_with_clusters(self):
        # # use only K-1 binary features
        min_mu = 0.99
        # 1st feature: fully active
        self.alphas[0] = np.log(min_mu / (1.0 - min_mu))
        self.zmat[0,:] = True

        freqlist = []
        for p in range(self.P):
            freqlist.append(defaultdict(int))
        for l in range(self.L):
            for p in range(self.P):
                if self.mvs is None or self.mvs[l,p] == False:
                    freqlist[p][self.mat[l,p]] += 1
        for p in range(self.P):
            j_start, T = self.p2jT[p]
            if self.flist[p]["type"] == "cat":
                freq = np.array([freqlist[p][t] for t in range(T)], dtype=np.float32) + 0.5
                freq /= freq.sum()
                w = np.log(freq)
                w -= max(w.min(), -10.0)
                self.wmat[0,j_start:j_start+T] = w
            elif self.flist[p]["type"] == "bin":
                pi = (freqlist[p][1] + 0.5) / (sum(freqlist[p].values()) + 1.0)
                # sys.stderr.write("p\t{}\n".format(p))
                self.wmat[0,j_start] = np.log(pi / (1.0 - pi) + 1E-5)
            elif self.flist[p]["type"] == "count":
                avgl, total = 0, 0
                for k, v in freqlist[p].items():
                    total += v
                    avgl += k * v
                avgl /= total
                # sys.stderr.write("avgl\t{}\n".format(avgl))
                self.wmat[0,j_start] = np.log(np.exp(min(avgl, 10.0)) - 1.0)

        # subsequent K-1 features
        for k in range(1, self.K):
            min_mu = max(min_mu * min(np.random.beta(19.0, 1.0), 0.99), 0.01)
            self.alphas[k] = np.log(min_mu / (1.0 - min_mu))
            self.wmat[k] = np.random.normal(loc=0.0, scale=0.1, size=self.M)
            self.zmat[k,:] = (np.random.rand(self.L) < min_mu)
            sys.stderr.write("{}\n".format(min_mu))
            sys.stderr.write("{}\n".format(self.zmat[k].sum()))
        self.calc_theta()

    def calc_loglikelihood(self):
        # self.calc_theta()
        ll = 0.0
        ls = np.arange(self.L, dtype=np.int32)
        for p in range(self.P):
            j_start, T = self.p2jT[p]
            xs = self.mat[:,p]
            if self.flist[p]["type"] == "cat":
                ll += np.log(self.theta[ls,j_start+xs] + 1E-20).sum()
            elif self.flist[p]["type"] == "bin":
                ll += np.log(np.where(xs == 0, 1.0 - self.theta[:,j_start], self.theta[:,j_start]) + 1E-20).sum()
                pass
            elif self.flist[p]["type"] == "count":
                ll += (xs * np.log(self.theta[:,j_start] + 1E-20)).sum() - self.theta[:,j_start].sum() - gammaln(xs + 1).sum()
            else:
                raise NotImplementedError
        return ll

    def calc_theta(self):
        self.theta_tilde[...] = np.matmul(self.zmat.T, self.wmat) # (K x L)^T x (K x M) -> (L x M)
        for p in range(self.P):
            j_start, T = self.p2jT[p]
            if self.flist[p]["type"] == "cat":
                e_theta_tilde = np.exp(self.theta_tilde[:,j_start:j_start+T] - self.theta_tilde[:,j_start:j_start+T].max(axis=1).reshape(self.L, 1))
                self.theta[:,j_start:j_start+T] = e_theta_tilde / e_theta_tilde.sum(axis=1).reshape(self.L, 1)
            elif self.flist[p]["type"] == "bin":
                # sigmoid
                self.theta[:,j_start] = 0.5 * np.tanh(0.5 * self.theta_tilde[:,j_start]) + 0.5
            elif self.flist[p]["type"] == "count":
                # softplus
                self.theta[:,j_start] = np.fmax(self.theta_tilde[:,j_start], 0) + np.log1p(np.exp(-np.fabs(self.theta_tilde[:,j_start])))
            else:
                raise NotImplementedError

    def init_tasks(self, a_repeat=1, sample_w=True):
        self.tasks = []
        for l, p in self.mv_list:
            if self.flist[p]["type"] == "cat":
                self.tasks.append((self.S_X_CAT, (l, p)))
            elif self.flist[p]["type"] == "bin":
                self.tasks.append((self.S_X_BIN, (l, p)))
            elif self.flist[p]["type"] == "count":
                self.tasks.append((self.S_X_CNT, (l, p)))
            else:
                raise NotImplementedError
        for k in range(self.K):
            if self.bias and k == 0:
                continue
            for a in range(a_repeat):
                if not self.only_alphas:
                    if not self.drop_vs:
                        self.tasks.append((self.S_Z_V, k))
                    if not self.drop_hs:
                        self.tasks.append((self.S_Z_H, k))
                self.tasks.append((self.S_Z_A, k))
        for l in range(self.L):
            if self.bias:
                self.tasks += map(lambda k: (self.S_Z, (l, k)), range(1, self.K))
            else:
                self.tasks += map(lambda k: (self.S_Z, (l, k)), range(self.K))
        if sample_w:
            # self.tasks += map(lambda p: (self.S_W_HMC, p), range(self.P))
            self.tasks += map(lambda k: (self.S_W_HMC, k), range(self.K))

    def sample(self, _iter=0, maxanneal=0, itemp=-1):
        # inverse of temperature
        if itemp > 0:
            sys.stderr.write("\t\titemp\t{}\n".format(itemp))
        elif _iter >= maxanneal:
            itemp = 1.0
        else:
            itemp = 0.1 + 0.9 * _iter / maxanneal
            sys.stderr.write("\t\titemp\t{}\n".format(itemp))

        c_x_cat = [0, 0]
        c_x_bin = [0, 0]
        c_x_cnt = [0, 0]
        c_z = [0, 0]
        c_zx = [0, 0] # changed, total
        c_z_v = [0, 0]
        c_z_h = [0, 0]
        c_z_a = [0, 0]
        c_w_hmc = [0, 0]
        random.shuffle(self.tasks)
        for t_type, t_val in self.tasks:
            if t_type == self.S_X_CAT:
                l, p = t_val
                changed = self.sample_x_cat(l, p)
                c_x_cat[changed] += 1
            elif t_type == self.S_X_BIN:
                l, p = t_val
                changed = self.sample_x_bin(l, p)
                c_x_bin[changed] += 1
            elif t_type == self.S_X_CNT:
                l, p = t_val
                changed = self.sample_x_count(l, p)
                c_x_cnt[changed] += 1
            elif t_type == self.S_Z:
                l, k = t_val
                # changed = self.sample_z(l, k, itemp=itemp)
                changed, c, t = self.sample_zx(l, k, itemp=itemp)
                c_z[changed] += 1
                c_zx[0] += c
                c_zx[1] += t
            elif t_type == self.S_W_HMC:
                changed = self.sample_w_hmc(t_val)
                c_w_hmc[changed] += 1
            elif t_type == self.S_Z_V:
                changed = self.sample_autologistic(t_type, t_val)
                c_z_v[changed] += 1
            elif t_type == self.S_Z_H:
                changed = self.sample_autologistic(t_type, t_val)
                c_z_h[changed] += 1
            elif t_type == self.S_Z_A:
                changed = self.sample_autologistic(t_type, t_val)
                c_z_a[changed] += 1
            else:
                raise NotImplementedError
        self.calc_theta() # fix numerical errors
        if sum(c_x_cat) > 0:
            sys.stderr.write("\tx_cat\t{}\n".format(float(c_x_cat[1]) / sum(c_x_cat)))
        if sum(c_x_bin) > 0:
            sys.stderr.write("\tx_bin\t{}\n".format(float(c_x_bin[1]) / sum(c_x_bin)))
        if sum(c_x_cnt) > 0:
            sys.stderr.write("\tx_cnt\t{}\n".format(float(c_x_cnt[1]) / sum(c_x_cnt)))
        sys.stderr.write("\tz\t{}\n".format(float(c_z[1]) / sum(c_z)))
        if c_zx[1] > 0:
            sys.stderr.write("\tzx\t{}\n".format(float(c_zx[0]) / c_zx[1]))
        if sum(c_w_hmc) > 0:
            sys.stderr.write("\tw_hmc\t{}\n".format(float(c_w_hmc[1]) / sum(c_w_hmc)))
        if not self.only_alphas:
            if sum(c_z_v) > 0:
                sys.stderr.write("\tz_v\t{}\n".format(float(c_z_v[1]) / sum(c_z_v)))
            if sum(c_z_h) > 0:
                sys.stderr.write("\tz_h\t{}\n".format(float(c_z_h[1]) / sum(c_z_h)))
        if sum(c_z_a) > 0:
            sys.stderr.write("\tz_a\t{}\n".format(float(c_z_a[1]) / sum(c_z_a)))
        if not self.only_alphas:
            sys.stderr.write("\tv\tavg\t{}\tmax\t{}\n".format(self.vks.mean(), self.vks.max()))
            sys.stderr.write("\th\tavg\t{}\tmax\t{}\n".format(self.hks.mean(), self.hks.max()))
        sys.stderr.write("\ta\tavg\t{}\tvar\t{}\n".format(self.alphas.mean(), self.alphas.var()))

    def sample_x_cat(self, l, p):
        assert(self.mvs is not None and self.mvs[l,p])
        j_start, T = self.p2jT[p]
        x_old = self.mat[l,p]
        self.mat[l,p] = np.random.choice(T, p=self.theta[l,j_start:j_start+T])
        return False if x_old == self.mat[l,p] else True

    def sample_x_bin(self, l, p):
        assert(self.mvs is not None and self.mvs[l,p])
        j_start, T = self.p2jT[p]
        x_old = self.mat[l,p]
        self.mat[l,p] = np.random.random() < self.theta[l,j_start]
        return False if x_old == self.mat[l,p] else True

    def sample_x_count(self, l, p):
        assert(self.mvs is not None and self.mvs[l,p])
        j_start, T = self.p2jT[p]
        x_old = self.mat[l,p]
        self.mat[l,p] = np.random.poisson(lam=self.theta[l,j_start])
        return False if x_old == self.mat[l,p] else True

    def sample_z(self, l, k, itemp=1.0):
        assert(not self.bias or not k == 0)
        z_old = self.zmat[k,l]
        logprob0, logprob1 = 0.0, 0.0
        if not self.only_alphas:
            vcount = np.bincount(self.zmat[k,self.vnet[l]], minlength=2)
            logprob0 += self.vks[k] * vcount[0]
            logprob1 += self.vks[k] * vcount[1]
            hcount = np.bincount(self.zmat[k,self.hnet[l]], minlength=2)
            logprob0 += self.hks[k] * hcount[0]
            logprob1 += self.hks[k] * hcount[1]
        logprob1 += self.alphas[k]

        theta_new = np.empty_like(self.theta[l,:])
        theta_tilde_new = np.copy(self.theta_tilde[l,:])
        if z_old == False:
            # proposal: 1
            theta_tilde_new += self.wmat[k,:]
            logprob_old, logprob_new = logprob0, logprob1
        else:
            # proposal: 0
            theta_tilde_new -= self.wmat[k,:]
            logprob_old, logprob_new = logprob1, logprob0
        for p in range(self.P):
            j_start, T = self.p2jT[p]
            x = self.mat[l,p]
            if self.flist[p]["type"] == "cat":
                logprob_old += np.log(self.theta[l,j_start+x] + 1E-20)
                e_theta_tilde = np.exp(theta_tilde_new[j_start:j_start+T] - theta_tilde_new[j_start:j_start+T].max())
                theta_new[j_start:j_start+T] = e_theta_tilde / e_theta_tilde.sum()
                logprob_new += np.log(theta_new[j_start+x] + 1E-20)
            elif self.flist[p]["type"] == "bin":
                logprob_old += np.log((self.theta[l,j_start] if x == 1 else (1.0 - self.theta[l,j_start])) + 1E-20)
                theta_new[j_start] = 0.5 * np.tanh(0.5 * theta_tilde_new[j_start]) + 0.5
                logprob_new += np.log((theta_new[j_start] if x == 1 else (1.0 - theta_new[j_start])) + 1E-20)
            elif self.flist[p]["type"] == "count":
                logprob_old += x * np.log(self.theta[l,j_start] + 1E-20)
                theta_new[j_start] = np.fmax(theta_tilde_new[j_start], 0) + np.log1p(np.exp(-np.fabs(theta_tilde_new[j_start])))
                logprob_new += x * np.log(theta_new[j_start] + 1E-20)
            else:
                raise NotImplementedError
        if itemp != 1.0:
            logprob_old *= itemp
            logprob_new *= itemp
        accepted = np.bool_(rand_partition_log((logprob_old, logprob_new)))
        if accepted:
            if z_old == False:
                # 0 -> 1
                self.zmat[k,l] = True
            else:
                # 1 -> 0
                self.zmat[k,l] = False
            self.theta_tilde[l,:] = theta_tilde_new
            self.theta[l,:] = theta_new
            return True
        else:
            return False

    def sample_zx(self, l, k, itemp=1.0):
        assert(not self.bias or not k == 0)
        z_old = self.zmat[k,l]
        logprob0, logprob1 = 0.0, 0.0
        if not self.only_alphas:
            vcount = np.bincount(self.zmat[k,self.vnet[l]], minlength=2)
            logprob0 += self.vks[k] * vcount[0]
            logprob1 += self.vks[k] * vcount[1]
            hcount = np.bincount(self.zmat[k,self.hnet[l]], minlength=2)
            logprob0 += self.hks[k] * hcount[0]
            logprob1 += self.hks[k] * hcount[1]
        logprob1 += self.alphas[k]

        theta_new = np.empty_like(self.theta[l,:])
        theta_tilde_new = np.copy(self.theta_tilde[l,:])
        if z_old == False:
            # proposal: 1
            theta_tilde_new += self.wmat[k,:]
            logprob_old, logprob_new = logprob0, logprob1
        else:
            # proposal: 0
            theta_tilde_new -= self.wmat[k,:]
            logprob_old, logprob_new = logprob1, logprob0
        xs_new = self.mat[l,:].copy()
        for p, (x, is_missing) in enumerate(zip(self.mat[l,:], self.mvs[l,:])):
            j_start, T = self.p2jT[p]
            if self.flist[p]["type"] == "cat":
                e_theta_tilde = np.exp(theta_tilde_new[j_start:j_start+T] - theta_tilde_new[j_start:j_start+T].max())
                theta_new[j_start:j_start+T] = e_theta_tilde / e_theta_tilde.sum()
                if is_missing:
                    xs_new[p] = np.random.choice(T, p=theta_new[j_start:j_start+T])
                else:
                    logprob_old += np.log(self.theta[l,j_start+x] + 1E-20)
                    logprob_new += np.log(theta_new[j_start+xs_new[p]] + 1E-20)
            elif self.flist[p]["type"] == "bin":
                theta_new[j_start] = 0.5 * np.tanh(0.5 * theta_tilde_new[j_start]) + 0.5
                if is_missing:
                    xs_new[p] = np.random.random() < theta_new[j_start]
                else:
                    logprob_old += np.log((self.theta[l,j_start] if x == 1 else (1.0 - self.theta[l,j_start])) + 1E-20)
                    logprob_new += np.log((theta_new[j_start] if xs_new[p] == 1 else (1.0 - theta_new[j_start])) + 1E-20)
            elif self.flist[p]["type"] == "count":
                theta_new[j_start] = np.fmax(theta_tilde_new[j_start], 0) + np.log1p(np.exp(-np.fabs(theta_tilde_new[j_start])))
                if is_missing:
                    xs_new[p] = np.random.poisson(lam=theta_new[j_start])
                else:
                    logprob_old += x * np.log(self.theta[l,j_start] + 1E-20)
                    logprob_new += xs_new[p] * np.log(theta_new[j_start] + 1E-20)
            else:
                raise NotImplementedError
        if itemp != 1.0:
            logprob_old *= itemp
            logprob_new *= itemp
        accepted = np.bool_(rand_partition_log((logprob_old, logprob_new)))
        if accepted:
            if z_old == False:
                # 0 -> 1
                self.zmat[k,l] = True
            else:
                # 1 -> 0
                self.zmat[k,l] = False
            self.theta_tilde[l,:] = theta_tilde_new
            self.theta[l,:] = theta_new
            changed = (self.mat[l,:] != xs_new).sum()
            self.mat[l,:] = xs_new
            return True, changed, self.mvs[l,:].sum()
        else:
            return False, 0, self.mvs[l,:].sum()

    def sample_autologistic(self, t_type, k):
        logr = 0.0
        if t_type == self.S_Z_A:
            oldval = self.alphas[k]
            pivot = min((self.zmat[k].sum() + 0.01) / self.L, 0.99)
            pivot = np.log(pivot / (1.0 - pivot))
            oldmean = (oldval + pivot) / 2.0
            oldscale = max(abs(oldval - pivot), 0.001)
            newval = np.random.normal(loc=oldmean, scale=oldscale)
            newmean = (newval + pivot) / 2.0
            newscale = max(abs(newval - pivot), 0.001)
            # q(theta|theta', x) / q(theta'|theta, x)
            logr += -((oldval - newmean) ** 2) / (2.0 * newscale * newscale) - np.log(newscale) \
                    + ((newval - oldmean) ** 2) / (2.0 * oldscale * oldscale) + np.log(oldscale)
            # P(theta') / P(theta)
            logr += (oldval * oldval - newval * newval) / (2.0 * self.norm_sigma * self.norm_sigma)
            # skip: q(theta|theta', x) / q(theta'|theta, x) for symmetric proposal
            v, h, a = self.vks[k], self.hks[k], newval
        else:
            assert(not self.only_alphas)
            assert(not (t_type == self.S_Z_V and self.drop_vs))
            assert(not (t_type == self.S_Z_H and self.drop_hs))
            if t_type == self.S_Z_V:
                oldval = self.vks[k]
            else:
                oldval = self.hks[k]
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
            if t_type == self.S_Z_V:
                v, h, a = newval, self.hks[k], self.alphas[k]
                net = self.vnet
            else:
                v, h, a = self.vks[k], newval, self.alphas[k]
                net = self.hnet
        zvect = self.zmat[k].copy()
        llist = np.arange(self.L)
        np.random.shuffle(llist)
        for l in llist:
            logprob0, logprob1 = (0.0, 0.0)
            if not self.only_alphas:
                vcount = np.bincount(self.zmat[k,self.vnet[l]], minlength=2)
                logprob0 += v * vcount[0]
                logprob1 += v * vcount[1]
                hcount = np.bincount(self.zmat[k,self.hnet[l]], minlength=2)
                logprob0 += h * hcount[0]
                logprob1 += h * hcount[1]
            logprob1 += a
            zvect[l] = rand_partition_log([logprob0, logprob1])
        # V_oldsum = self._neighbor_sum(self.zmat[k], self.vnet)
        # H_oldsum = self._neighbor_sum(self.zmat[k], self.hnet)
        # A_oldsum = self.zmat[k].sum()
        # V_newsum = self._neighbor_sum(zvect, self.vnet)
        # H_newsum = self._neighbor_sum(zvect, self.hnet)
        # A_newsum = zvect.sum()
        # logr += (self.vks[k] * V_newsum + self.hks[k] * H_newsum + self.alphas[k] * A_newsum) \
        #         + (v * V_oldsum + h * H_oldsum + a * A_oldsum) \
        #         - (self.vks[k] * V_oldsum + self.hks[k] * H_oldsum + self.alphas[k] * A_oldsum) \
        #         - (v * V_newsum + h * H_newsum + a * A_newsum)
        # logr += (self.vks[k] * (V_newsum - V_oldsum) + self.hks[k] * (H_newsum - H_oldsum) + self.alphas[k] * (A_newsum - A_oldsum)) \
        #         - (v * (V_newsum - V_oldsum) + h * (H_newsum - H_oldsum) + a * (A_newsum - A_oldsum))
        # logr += ((self.vks[k] - v) * (V_newsum - V_oldsum) + (self.hks[k] - h) * (H_newsum - H_oldsum) + (self.alphas[k] - a) * (A_newsum - A_oldsum))
        if t_type == self.S_Z_A:
            logr += (oldval - newval) * (zvect.sum() - self.zmat[k].sum())
            if logr >= 0 or np.log(np.random.rand()) < logr:
                # accept
                self.alphas[k] = newval
                return True
            else:
                return False
        else:
            oldsum = self._neighbor_sum(self.zmat[k], net)
            newsum = self._neighbor_sum(zvect, net)
            logr += (oldval - newval) * (newsum - oldsum)
            if logr >= 0 or np.log(np.random.rand()) < logr:
                # accept
                if t_type == self.S_Z_V:
                    self.vks[k] = newval
                else:
                    self.hks[k] = newval
                return True
            else:
                return False

    def _neighbor_sum(self, zvect, net):
        s = 0
        for l in range(self.L):
            s += (zvect[net[l]] == zvect[l]).sum()
        assert(s % 2 == 0)
        return s / 2

    # def sample_w_hmc(self, p):
    #     sigma2 = (self.sigma + 1.0) / self.sigma
    #     j_start, T = self.p2jT[p]
    #     xs = self.mat[:,p]
    #     ls = np.arange(self.L, dtype=np.int32)

    #     def U(wcs):
    #         ll = 0.5 * (self.sigma + 1.0) * np.log(1.0 + (wcs * wcs) / self.sigma).sum()
    #         theta_tilde = np.matmul(self.zmat.T, wcs) # [K, L]^T x [K, T] -> [L, T]
    #         if self.flist[p]["type"] == "cat":
    #             e_theta_tilde = np.exp(theta_tilde - theta_tilde.max(axis=1).reshape(self.L, 1))
    #             theta = e_theta_tilde / e_theta_tilde.sum(axis=1).reshape(self.L, 1)
    #             ll -= np.log(self.theta[ls,xs] + 1E-20).sum()
    #         elif self.flist[p]["type"] == "bin":
    #             theta = 0.5 * np.tanh(0.5 * theta_tilde) + 0.5
    #             ll -= np.log(np.where(xs == 0, 1.0 - theta, theta) + 1E-20).sum()
    #         elif self.flist[p]["type"] == "count":
    #             theta = np.fmax(theta_tilde, 0) + np.log1p(np.exp(-np.fabs(theta_tilde)))
    #             ll -= (xs * np.log(theta + 1E-20)).sum() - theta.sum()
    #         else:
    #             raise NotImplementedError
    #         return ll
    #     def gradU(wcs):
    #         grad = sigma2 * (wcs / (1.0 + (wcs * wcs) / self.sigma))
    #         theta_tilde = np.matmul(self.zmat.T, wcs) # [K, L]^T x [K, T] -> [L, T]
    #         if self.flist[p]["type"] == "cat":
    #             e_theta_tilde = np.exp(theta_tilde - theta_tilde.max(axis=1).reshape(self.L, 1))
    #             error = -e_theta_tilde / e_theta_tilde.sum(axis=1).reshape(self.L, 1)
    #             error[ls,xs] += 1
    #         elif self.flist[p]["type"] == "bin":
    #             theta = 0.5 * np.tanh(0.5 * theta_tilde) + 0.5
    #             error = np.where(xs.reshape(self.L, 1) == 0, -theta, 1.0 - theta)
    #         elif self.flist[p]["type"] == "count":
    #             theta = np.fmax(theta_tilde, 0) + np.log1p(np.exp(-np.fabs(theta_tilde)))
    #             error = (xs.reshape(self.L, 1) / theta - 1.0) / (1.0 + np.exp(-theta_tilde))
    #         else:
    #             raise NotImplementedError
    #         grad -= np.matmul(self.zmat, error)
    #         return grad
    #     accepted, wcs = hmc(U, gradU, self.HMC_EPSILON, self.HMC_L, self.wmat[:,j_start:j_start+T])
    #     if accepted:
    #         # update theta_tilde
    #         self.wmat[:,j_start:j_start+T] = wcs
    #         self.theta_tilde[:,j_start:j_start+T] = np.matmul(self.zmat.T, self.wmat[:,j_start:j_start+T]) # [K, L]^T x [K, T] -> [L, T]
    #         if self.flist[p]["type"] == "cat":
    #             e_theta_tilde = np.exp(self.theta_tilde[:,j_start:j_start+T] - self.theta_tilde[:,j_start:j_start+T].max(axis=1).reshape(self.L, 1))
    #             self.theta[:,j_start:j_start+T] = e_theta_tilde / e_theta_tilde.sum(axis=1).reshape(self.L, 1)
    #         elif self.flist[p]["type"] == "bin":
    #             self.theta[:,j_start] = 0.5 * np.tanh(0.5 * self.theta_tilde[:,j_start]) + 0.5
    #         elif self.flist[p]["type"] == "count":
    #             self.theta[:,j_start] = np.fmax(self.theta_tilde[:,j_start], 0) + np.log1p(np.exp(-np.fabs(self.theta_tilde[:,j_start])))
    #         else:
    #             raise NotImplementedError
    #         return True
    #     else:
    #         return False

    def sample_w_hmc(self, k):
        def U(Mvect):
            # ll = -norm.logpdf(Mvect, 0.0, scale=self.sigma).sum()
            ll = 0.5 * (self.sigma + 1.0) * np.log(1.0 + (Mvect * Mvect) / self.sigma).sum()
            for l in range(self.L):
                if self.zmat[k,l] == False:
                    continue
                theta_tilde = self.theta_tilde[l] - self.wmat[k] + Mvect
                for p in range(self.P):
                    j_start, T = self.p2jT[p]
                    x = self.mat[l,p]
                    if self.flist[p]["type"] == "cat":
                        theta_tilde2 = theta_tilde[j_start:j_start+T] - theta_tilde[j_start:j_start+T].max()
                        ll -= theta_tilde2[x] - np.log(np.exp(theta_tilde2).sum())
                    elif self.flist[p]["type"] == "bin":
                        theta = 0.5 * np.tanh(0.5 * theta_tilde[j_start]) + 0.5
                        ll -= np.log((theta if x == 1 else (1.0 - theta)) + 1E-20)
                    elif self.flist[p]["type"] == "count":
                        theta = np.fmax(theta_tilde[j_start], 0) + np.log1p(np.exp(-np.fabs(theta_tilde[j_start])))
                        ll -= x * np.log(theta + 1E-20) - theta
                    else:
                        raise NotImplementedError
            return ll
        sigma2 = (self.sigma + 1.0) / self.sigma
        def gradU(Mvect):
            grad = sigma2 * (Mvect / (1.0 + (Mvect * Mvect) / self.sigma))
            for l in range(self.L):
                if self.zmat[k,l] == False:
                    continue
                theta_tilde = self.theta_tilde[l] - self.wmat[k] + Mvect
                for p in range(self.P):
                    j_start, T = self.p2jT[p]
                    x = self.mat[l,p]
                    if self.flist[p]["type"] == "cat":
                        e_theta_tilde = np.exp(theta_tilde[j_start:j_start+T] - theta_tilde[j_start:j_start+T].max())
                        theta = e_theta_tilde / e_theta_tilde.sum()
                        grad[j_start:j_start+T] += theta
                        grad[j_start + x] -= 1
                    elif self.flist[p]["type"] == "bin":
                        theta = 0.5 * np.tanh(0.5 * theta_tilde[j_start]) + 0.5
                        grad[j_start] -= (1.0 - theta) if x == 1 else (-theta)
                    elif self.flist[p]["type"] == "count":
                        theta = np.fmax(theta_tilde[j_start], 0) + np.log1p(np.exp(-np.fabs(theta_tilde[j_start])))
                        grad[j_start] -= (x / theta - 1.0) / (1.0 + np.exp(-theta_tilde[j_start]))
                    else:
                        raise NotImplementedError                    
            return grad
        accepted, Mvect = hmc(U, gradU, self.hmc_epsilon, self.hmc_l, self.wmat[k])
        if accepted:
            # update theta_tilde
            for l in range(self.L):
                if self.zmat[k,l] == False:
                    continue
                self.theta_tilde[l] += Mvect - self.wmat[k]
                for p in range(self.P):
                    j_start, T = self.p2jT[p]
                    if self.flist[p]["type"] == "cat":
                        e_theta_tilde = np.exp(self.theta_tilde[l,j_start:j_start+T] - self.theta_tilde[l,j_start:j_start+T].max())
                        self.theta[l,j_start:j_start+T] = e_theta_tilde / e_theta_tilde.sum()
                    elif self.flist[p]["type"] == "bin":
                        self.theta[l,j_start] = 0.5 * np.tanh(0.5 * self.theta_tilde[l,j_start]) + 0.5
                    elif self.flist[p]["type"] == "count":
                        self.theta[l,j_start] = np.fmax(self.theta_tilde[l,j_start], 0) + np.log1p(np.exp(-np.fabs(self.theta_tilde[l,j_start])))
                    else:
                        raise NotImplementedError
                # assert(all(self.theta_tilde[i] > 0))
            self.wmat[k] = Mvect
            return True
        else:
            return False
