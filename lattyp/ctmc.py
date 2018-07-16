#!/bin/env python
# -*- coding: utf-8 -*-
#
# continuous time Markov chain models
#
import numpy as np

class TwoStateCTMC(object):
    """
    two-state CTMC with a gamma prior
    """
    D = 2
    def __init__(self, params=None, scale=0.00005):
        self.scale = scale
        if params is None:
            params = 0.0001 + np.random.gamma(shape=1.0, scale=self.scale, size=2)
        self.set_params(params)

    def get_params(self):
        return np.array([self.alpha, self.beta])

    def set_params(self, params):
        assert(len(params) == 2)
        self.alpha, self.beta = params
        self._sum = self.alpha + self.beta

    def get_param_logprob(self):
        return -(self.alpha + self.beta) / self.scale

    def get_stationary_probs(self):
        return np.array([self.beta / self._sum, self.alpha / self._sum])

    def get_stationary_grad_logprobs(self):
        # 1st dim. alpha or beta
        # 2nd dim. 0 or 1
        c = -1.0 / self._sum
        return np.array([[c, 1.0 / self.alpha + c],
                         [1.0 / self.beta + c, c]])

    def get_transition_probs(self, time, rate=1.0):
        assert(time >= 0.0)
        probmat = np.zeros((2,2))
        ep = np.exp(-self._sum * time * rate)
        probmat[0,0] = (self.beta + self.alpha * ep) / self._sum
        probmat[1,1] = (self.alpha + self.beta * ep) / self._sum
        probmat[0,1] = 1.0 - probmat[0,0]
        probmat[1,0] = 1.0 - probmat[1,1]
        return probmat

    def get_transition_grad_logprobs(self, time, rate=1.0):
        # for HMC
        assert(time >= 0.0)
        gradmat = -np.ones((2,2,2)) / self._sum
        rtime = time * rate
        ep = np.exp(-self._sum * rtime)
        denom0x = self.beta + self.alpha * ep
        gradmat[0,0,0] += (ep * (1.0 - self.alpha * rtime)) / denom0x
        gradmat[1,0,0] += (1.0 - self.alpha * rtime * ep) / denom0x
        denom1 = 1.0 - ep + 1E-10
        grad0 = (rtime * ep) / denom1
        gradmat[0,0,1] += 1.0 / self.alpha + grad0
        gradmat[1,0,1] += grad0
        gradmat[0,1,0] += grad0
        gradmat[1,1,0] += 1.0 / self.beta + grad0
        denom1x = self.alpha + self.beta * ep
        gradmat[0,1,1] += (1.0 - self.beta * rtime * ep) / denom1x
        gradmat[1,1,1] += (ep * (1.0 - self.beta * rtime)) / denom1x
        return gradmat

class MultiStateCTMC(object):
    def __init__(self, params=None, D=-1, scale=0.00005):
        self.scale = scale
        if params is None:
            assert(D > 0)
            params = 0.0001 + np.random.gamma(shape=1.0, scale=self.scale, size=(D,D))
        self.set_params(params)

    def get_params(self):
        return self.params

    def set_params(self, params):
        assert(len(params.shape) == 2 and params.shape[0] == params.shape[1])
        self.D = params.shape[0]
        self.params = params
        for i in range(self.D):
            s = 0.0
            for j in range(self.D):
                if i != j:
                    assert(params[i,j] >= 0)
                    s += params[i,j]
            params[i,i] = -s
        if hasattr(self, "w"):
            delattr(self, "w")
        if hasattr(self, "p"):
            delattr(self, "p")

    def get_param_logprob(self):
        return -(self.params.sum() - np.diag(self.params).sum()) / self.scale
        # return -(self.alpha + self.beta) / self.scale

    def _eig(self):
        self.w, self.v = np.linalg.eig(self.params)
        self.vi = np.linalg.inv(self.v)

    def get_stationary_probs(self):
        # normalize the right eigenvector whose eigenvalue is close to zero
        if hasattr(self, "p"):
            return self.p
        if not hasattr(self, "w"):
            self._eig()
        # a = np.isclose(np.exp(self.w), 1.0, rtol=1E-6, atol=1E-6)
        # if a.sum() <= 0:
        #     raise FloatingPointError
        # for k, v in enumerate(a):
        #     if v: break
        k = np.argmin(np.absolute(self.w))
        if not np.isclose(np.exp(self.w[k]), 1.0, rtol=1E-6, atol=1E-6):
            raise FloatingPointError
        p = self.vi[k]
        if np.iscomplexobj(p): # numerical errors
            p = p.real
        ps = p.sum() # may be negative
        self.p = p / ps
        return self.p

    def get_stationary_grad_logprobs(self):
        raise NotImplementedError

    def get_transition_probs(self, time, rate=1.0):
        assert(time >= 0.0)
        if not hasattr(self, "w"):
            self._eig()
        probs = np.dot(np.dot(self.v, np.diag(np.exp(self.w * time))), self.vi)
        if np.iscomplexobj(probs):
            probs = probs.real
        return probs

    def get_transition_grad_logprobs(self, time, rate=1.0):
        raise NotImplementedError
