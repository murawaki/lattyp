#!/bin/env python
# -*- coding: utf-8 -*-
#
# CTMC-based phylogenetic tree
#
import numpy as np
from scipy.stats import norm
from scipy.misc import logsumexp

class GaussianPrior(object):
    def __init__(self, params):
        assert(len(params) == 2)
        self.mean = params[0]
        self.sd = params[1]

    def init_val(self):
        return self.mean

    def get_interval(self):
        return [-np.inf, np.inf]

    def prob(self, x):
        return norm.pdf(x, loc=self.mean, scale=self.sd)

    def logprob(self, x):
        return norm.logpdf(x, loc=self.mean, scale=self.sd)

class GaussianMixturePrior(object):
    def __init__(self, params):
        assert(len(params) % 3 == 0)
        self.dists = []
        weights = []
        for i in range(len(params) // 3):
            weights.append(params[i*3])
            self.dists.append(GaussianPrior((params[i*3+1], params[i*3+2])))
        self.weights = np.array(weights)
        self.logweights = np.log(self.weights)

    def init_val(self):
        x = 0.0
        for i, (w, d) in enumerate(zip(self.weights, self.dists)):
            x += w * d.init_val()
        return x

    def get_interval(self):
        return [-np.inf, np.inf]

    def prob(self, x):
        p = 0.0
        for i, (w, d) in enumerate(zip(self.weights, self.dists)):
            p += w * d.prob(x)
        return p

    def logprob(self, x):
        ll = []
        for i, (lw, d) in enumerate(zip(self.logweights, self.dists)):
            ll.append(lw + d.logprob(x))
        return logsumexp(ll)

class LogNormalPrior(object):
    def __init__(self, params):
        assert(len(params) in (2, 3))
        self.mean = params[0]
        self.sd = params[1]
        self.offset = params[2] if len(params) == 3 else 0

    def init_val(self):
        return self.offset + self.mean + self.sd * self.sd

    def get_interval(self):
        return [self.offset, np.inf]

    def prob(self, x):
        if x < self.offset + self.mean:
            return 0.0
        else:
            return norm.pdf(x - self.offset, loc=self.mean, scale=self.sd)

    def logprob(self, x):
        if x < self.offset + self.mean:
            return -np.inf
        else:
            return norm.logpdf(x - self.offset, loc=self.mean, scale=self.sd)

class UniformPrior(object):
    def __init__(self, params):
        assert(len(params) == 2)
        self.lb, self.ub = params

    def init_val(self):
        return (self.lb + self.ub) / 2.0

    def get_interval(self):
        return [self.lb, self.ub]

    def prob(self, x):
        if self.lb < x < self.ub:
            return 1.0
        else:
            return 0.0

    def logprob(self, x):
        if self.lb < x < self.ub:
            return 0.0
        else:
            return -np.inf

prior_dict = {
    "Gaussian": GaussianPrior,
    "GaussianMixture": GaussianMixturePrior,
    "LogNormal": LogNormalPrior,
    "Uniform": UniformPrior,
}

def create_node_priors(prior_specs):
    node_priors = {}
    for spec in prior_specs:
        if spec["dtype"] in prior_dict:
            d = prior_dict[spec["dtype"]](spec["params"])
        else:
            raise NotImplementedError
        node_priors[spec["glottocode"]] = d
    return node_priors
