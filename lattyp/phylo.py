#!/bin/env python
# -*- coding: utf-8 -*-
#
# CTMC-based phylogenetic tree
#
import sys
import random
import numpy as np
from scipy.misc import logsumexp
from scipy.stats import norm

from ctmc import TwoStateCTMC, MultiStateCTMC
from hmc import hmc
from rand_utils import rand_partition_log, slice_sampler1d

class CTMC_Sampler(object):
    S_CTMC_HMC = 1
    S_CTMC_ONE = 2
    S_CTMC_RATE = 3
    S_DATE_TREE = 4
    S_DATE_ROOT = 5
    S_DATE_INTERNAL = 6
    S_STATE_TREE = 7
    S_STATE_ROOT = 8
    S_STATE_INTERNAL = 9
    HMC_L = 10
    HMC_EPSILON = 0.01 # 0.1

    ABS_LIMIT = 30000.0 # upper bound of the limit of a root

    def __init__(self, K, states=2, ctmc_scale=0.00005):
        self.states = states # scalar(2) or list with size K
        self.K = K
        self.ctmc_scale = ctmc_scale
        self.ctmcs = []
        self.tasks = []
        self.calibrated_nodes = []

    def init_trees(self, trees, sample_dates=True, hmc_repeat=1):
        self._init_ctmc_task(trees, hmc_repeat=hmc_repeat)
        for tree in trees:
            if not tree.is_isolate:
                if sample_dates:
                    self.tasks.append((self.S_DATE_TREE, tree))
            self._init_node(tree, ub=None, sample_dates=sample_dates)

    def _init_ctmc_task(self, trees, hmc_repeat=1):
        if isinstance(self.states, int) and self.states == 2:
            for k in range(self.K):
                self.ctmcs.append(TwoStateCTMC(scale=self.ctmc_scale))
                for i in range(hmc_repeat):
                    self.tasks.append((self.S_CTMC_HMC, trees, k))
        else:
            for k, D in enumerate(self.states):
                self.ctmcs.append(MultiStateCTMC(D=D, scale=self.ctmc_scale))
                self.tasks.append((self.S_CTMC_RATE, trees, k))
                for i in range(D):
                    for j in range(D):
                        if i != j:
                            self.tasks.append((self.S_CTMC_ONE, trees, k, (i,j)))

    def _init_node(self, node, ub=None, sample_dates=True):
        # create initial state of MCMC
        # node state: random sampling based on the states of the children
        # date: between upper and lower bounds
        #
        # if sample_dates is False, sample only states
        if not node.is_date_frozen:
            if hasattr(node, "is_root") and node.is_root:
                if sample_dates:
                    self.tasks.append((self.S_DATE_ROOT, node))
                    self.tasks.append((self.S_DATE_TREE, node))
                for k in range(self.K):
                    self.tasks.append((self.S_STATE_ROOT, node, k))
                    self.tasks.append((self.S_STATE_TREE, node, k))
            else:
                if sample_dates:
                    self.tasks.append((self.S_DATE_INTERNAL, node))
                for k in range(self.K):
                    self.tasks.append((self.S_STATE_INTERNAL, node, k))
        if sample_dates:
            if hasattr(node, "prior"):
                self.calibrated_nodes.append(node)
                node.date = node.prior.init_val()
                cub = node.date
            elif hasattr(node, "date"):
                cub = node.date
            else:
                cub = ub
        else:
            cub = None # dummy
        # depth-first
        lb = 0.0
        for child in node.children:
            lb2 = self._init_node(child, ub=cub, sample_dates=sample_dates)
            if lb2 > lb:
                lb = lb2
        if not hasattr(node, "state"):
            assert(len(node.children) > 0)
            if isinstance(self.states, int) and self.states == 2:
                prob = np.zeros(self.K, dtype=np.float32)
                for child in node.children:
                    prob += child.state
                prob /= len(node.children)
                # prob1, prob2 = prob ** 0.5, (1.0 - prob) ** 0.5
                # prob = prob1 / (prob1 + prob2)
                node.state = np.random.random(self.K) < prob
            else:
                node.state = np.zeros(self.K, dtype=np.int32)
                for k, D in enumerate(self.states):
                    counts = np.zeros(D)
                    for child in node.children:
                        counts[child.state[k]] += 1.0
                    counts /= counts.sum()
                    node.state[k] = np.random.choice(D, p=counts)
        if not hasattr(node, "date"):
            assert(sample_dates)
            if ub is None:
                ub = lb + 1000.0
            assert(ub >= lb)
            node.date = lb + (ub - lb) * np.random.beta(2.0, 8.0) # np.random.uniform(low=lb, high=ub)
        else:
            assert(node.date >= lb)
        return node.date

    def sample(self, _iter=0):
        random.shuffle(self.tasks)
        c_hmc = [0, 0]
        c_co = [0, 0]
        c_cr = [0, 0]
        c_dt = [0, 0]
        c_dr = [0, 0] # 0: total count, 1: accum diff
        c_di = [0, 0] # 0: total count, 1: accum diff
        c_st = [0, 0]
        c_sr = [0, 0]
        c_si = [0, 0]
        for task in self.tasks:
            if task[0] == self.S_CTMC_HMC:
                _, trees, k = task
                _old = self.ctmcs[k].get_params()
                changed = self.sample_ctmc_hmc(trees, k)
                _new = self.ctmcs[k].get_params()
                sys.stderr.write("{}\t{} -> {}\n".format(changed, _old, _new))
                c_hmc[changed] += 1
            elif task[0] == self.S_CTMC_ONE:
                _, trees, k, pos = task
                changed = self.sample_ctmc_one(trees, k, pos)
                c_co[changed] += 1
            elif task[0] == self.S_CTMC_RATE:
                _, trees, k = task
                changed = self.sample_ctmc_rate(trees, k)
                c_cr[changed] += 1
            elif task[0] == self.S_DATE_TREE:
                _, node = task
                changed = self.sample_node_date_tree(node)
                c_dt[changed] += 1
            elif task[0] == self.S_DATE_ROOT:
                _, node = task
                diff = self.sample_node_date_root(node)
                c_dr[0] += 1
                c_dr[1] += diff
            elif task[0] == self.S_DATE_INTERNAL:
                _, node = task
                diff = self.sample_node_date_internal(node)
                c_di[0] += 1
                c_di[1] += diff
            elif task[0] == self.S_STATE_TREE:
                _, node, k = task
                # _before = self.logprob([node])
                changed, total = self.sample_node_state_tree(node, k)
                c_st[0] += total - changed
                c_st[1] += changed
                # _after = self.logprob([node])
                # sys.stdout.write("{}\t{}\n".format(_after - _before, node.glottocode))
            elif task[0] == self.S_STATE_ROOT:
                _, node, k = task
                # changed = self.sample_node_state_root(node, k)
                # c_sr[changed] += 1
            elif task[0] == self.S_STATE_INTERNAL:
                _, node, k = task
                # changed = self.sample_node_state_internal(node, k)
                # c_si[changed] += 1
            else:
                raise NotImplementedError
        for name, clist in (("dt", c_dt),
                            ("st", c_st),
                            ("sr", c_sr),
                            ("si", c_si),
                            ("hmc", c_hmc),
                            ("co", c_co),
                            ("cr", c_cr)):
            if sum(clist) > 0:
                sys.stderr.write("\t{}\t{}\n".format(name, float(clist[1]) / sum(clist)))
        for name, clist in (("dr", c_dr),
                            ("di", c_di),):
            if clist[0] > 0:
                sys.stderr.write("\t{}\t{}\n".format(name, clist[1] / clist[0]))

    def sample_node_date_tree(self, node):
        def _logprob_branch_rec(parent, node, pdate, rate=1.0):
            if node.is_date_frozen:
                mydate = node.date
            else:
                mydate = rate * node.date
                if pdate < mydate:
                    sys.stderr.write("  {}\t{}\n".format(pdate, mydate))
                    raise ValueError
            ll = 0.0
            if node in self.calibrated_nodes:
                ll += node.prior.logprob(mydate)
            ll += self._logprob_branch_each(pdate - mydate, parent, node, self.ctmcs)
            for child in node.children:
                ll += _logprob_branch_rec(node, child, mydate, rate=rate)
            return ll

        assert(node.parent is None and len(node.children) > 0)
        P_SIGMA = 1.0 # 0.5
        while True:
            rate = np.random.lognormal(mean=0.0, sigma=P_SIGMA) + 1E-10
            if node.date * rate < self.ABS_LIMIT:
                irate = 1.0 / rate
                lograte = np.log(rate)
                logirate = np.log(irate)
                break

        clp, plp = 0.0, 0.0
        if node.is_date_frozen:
            pdate = node.date
        else:
            pdate = rate * node.date
        if node in self.calibrated_nodes:
            clp += node.prior.logprob(node.date)
            plp += node.prior.logprob(pdate)
        try:
            for child in node.children:
                clp += _logprob_branch_rec(node, child, node.date, rate=1.0)
                plp += _logprob_branch_rec(node, child, pdate, rate=rate)
        except ValueError:
            sys.stderr.write("proposed tree violates constraints...skip\trate: {}\n".format(rate))
            return False
        logp = plp - clp + logirate - lograte
        if logp >= 0.0 or np.log(np.random.uniform(0.0, 1.0) + 1E-10) < logp:
            # accept
            stack = [node]
            while len(stack) > 0:
                _node = stack.pop(0)
                if not node.is_date_frozen:
                    _node.date *= rate
                for child in _node.children:
                    stack.append(child)
            return True
        else:
            # reject
            return False

    def sample_node_date_root(self, node):
        assert(node.parent is None and len(node.children) > 0)
        lb, ub = max(map(lambda x: x.date, node.children)), self.ABS_LIMIT # HACK
        if hasattr(node, "prior"):
            lb1, ub1 = node.prior.get_interval()
            if lb1 > lb:
                lb = lb1
            if ub1 < ub:
                ub = ub1
        current = node.date

        def get_logprob(pdate):
            logprob = 0.0
            if hasattr(node, "prior"):
                logprob += node.prior.logprob(pdate)
            # prob. of generating root state does not depend on the node's date
            for child in node.children:
                logprob += self._logprob_branch_each(pdate - child.date, node, child, self.ctmcs)
            return logprob

        try:
            newv = slice_sampler1d(get_logprob, current, min_x=lb, max_x=ub)
        except:
            n = node.glottocode if hasattr(node, "glottocode") else id(node)
            sys.stderr.write("date root: exceeded max. eval trials\t{}\t{}\n".format(n, node.date))
            return 0.0
        node.date = newv
        return abs(node.date - current)

    def sample_node_date_internal(self, node):
        assert(node.parent is not None and len(node.children) > 0)
        ub = node.parent.date
        lb = max(map(lambda x: x.date, node.children))
        assert(lb <= node.date <= ub)
        if hasattr(node, "prior"):
            lb1, ub1 = node.prior.get_interval()
            if lb1 > lb:
                lb = lb1
            if ub1 < ub:
                ub = ub1
        current = node.date

        def get_logprob(pdate):
            logprob = 0.0
            if hasattr(node, "prior"):
                logprob += node.prior.logprob(pdate)
            logprob += self._logprob_branch_each(node.parent.date - pdate, node.parent, node, self.ctmcs)
            for child in node.children:
                logprob += self._logprob_branch_each(pdate - child.date, node, child, self.ctmcs)
            return logprob

        try:
            newv = slice_sampler1d(get_logprob, current, min_x=lb, max_x=ub)
        except:
            n = node.glottocode if hasattr(node, "glottocode") else id(node)
            sys.stderr.write("date internal: exceeded max. eval trials\t{}\t{}\n".format(n, node.date))
            return 0.0
        node.date = newv
        return abs(node.date - current)

    def sample_node_state_tree(self, node, k):
        # forward-filtering, backward-sampling
        def _forward(node):
            Ln = np.zeros(D)
            for child in node.children:
                # Lc[s_j]
                Lc = _forward(child)
                time = node.date - child.date
                lmat = np.log(np.maximum(ctmc.get_transition_probs(time), 1E-10))
                for d in range(D):
                    Ln[d] += logsumexp(lmat[d,:] + Lc)
                tscores[child] = lmat
            if node.is_state_frozen:
                for d in range(D):
                    if node.state[k] != d:
                        Ln[d] = -np.inf
            nscores[node] = Ln
            return Ln
        def _backward(node):
            changed, total = 0, 0
            if not node.is_state_frozen:
                total += 1
                pval = node.parent.state[k].astype(np.int32)
                logprobs = nscores[node] + tscores[node][pval,:]
                newv = rand_partition_log(logprobs)
                if node.state[k] != newv:
                    changed += 1
                node.state[k] = newv
            for child in node.children:
                _changed, _total = _backward(child)
                changed += _changed
                total += _total
            return changed, total
        ctmc = self.ctmcs[k]
        D = ctmc.D
        nscores, tscores = {}, {}
        Ln = _forward(node)
        Ln += np.log(np.maximum(ctmc.get_stationary_probs(), 1E-10))
        changed, total = 0, 1
        newv = rand_partition_log(Ln)
        if node.state[k] != newv:
            changed += 1
        for child in node.children:
            _changed, _total = _backward(child)
            changed += _changed
            total += _total
        return changed, total

    def sample_node_state_root(self, node, k):
        assert(not node.is_state_frozen)
        oldv = node.state[k]
        ctmc = self.ctmcs[k]
        logprobs = np.log(np.maximum(ctmc.get_stationary_probs(), 1E-10))
        for child in node.children:
            time = node.date - child.date
            mat = ctmc.get_transition_probs(time)
            logprobs += np.log(np.maximum(mat[:,child.state[k].astype(np.int32)], 1E-10))
        newv = rand_partition_log(logprobs)
        node.state[k] = newv
        return not(oldv == newv)

    def sample_node_state_internal(self, node, k):
        assert(not node.is_state_frozen)
        oldv = node.state[k]
        ctmc = self.ctmcs[k]
        time = node.parent.date - node.date
        mat = ctmc.get_transition_probs(time)
        logprobs = np.log(np.maximum(mat[node.parent.state[k].astype(np.int32),:], 1E-10))
        for child in node.children:
            time = node.date - child.date
            mat = ctmc.get_transition_probs(time)
            logprobs += np.log(np.maximum(mat[:,child.state[k].astype(np.int32)], 1E-10))
        newv = rand_partition_log(logprobs)
        node.state[k] = newv
        return not(oldv == newv)

    def _logprob_k(self, trees, ctmc, k):
        def _logprob_branch_k(parent, node, ctmc, k):
            time = parent.date - node.date
            mat = ctmc.get_transition_probs(time)
            prob = mat[parent.state[k].astype(np.int32), node.state[k].astype(np.int32)]
            _ll = np.log(prob + 1E-10)
            for child in node.children:
                _ll += _logprob_branch_k(node, child, ctmc, k)
            return _ll

        stationary_logprobs = np.log(np.maximum(ctmc.get_stationary_probs(), 1E-10))
        ll = 0.0
        for tree in trees:
            ll += stationary_logprobs[tree.state[k].astype(np.int32)]
            for child in tree.children:
                ll += _logprob_branch_k(tree, child, ctmc, k)
        return ll

    def sample_ctmc_one(self, trees, k, pos):
        # for MultiStateCTMC
        P_SIGMA = 2.0 # 0.5
        rate = np.random.lognormal(mean=0.0, sigma=P_SIGMA) + 1E-10
        irate = 1.0 / rate
        lograte = np.log(rate)
        logirate = np.log(irate)
        ctmc1 = self.ctmcs[k]
        params = ctmc1.get_params().copy()
        clp = params[pos] / -self.ctmc_scale
        params[pos] = rate * params[pos] + 1E-20 # diag. elemtns will be fixed later
        ctmc2 = MultiStateCTMC(params=params)
        plp = params[pos] / -self.ctmc_scale
        clp += self._logprob_k(trees, ctmc1, k)
        plp += self._logprob_k(trees, ctmc2, k)
        logp = plp - clp + logirate - lograte
        if logp >= 0.0 or np.log(np.random.uniform(0.0, 1.0) + 1E-10) < logp:
            # accept
            self.ctmcs[k] = ctmc2
            return True
        else:
            # reject
            return False

    def sample_ctmc_rate(self, trees, k):
        # for MultiStateCTMC
        P_SIGMA = 0.1 # 1.0
        rate = np.random.lognormal(mean=0.0, sigma=P_SIGMA) + 1E-10
        irate = 1.0 / rate
        lograte = np.log(rate)
        logirate = np.log(irate)
        ctmc1 = self.ctmcs[k]
        params = rate * ctmc1.get_params() + 1E-20 # diag. elemtns will be fixed later
        clp = ctmc1.get_param_logprob()
        ctmc2 = MultiStateCTMC(params=params) # TODO: efficiency
        plp = ctmc2.get_param_logprob()
        clp += self._logprob_k(trees, ctmc1, k)
        plp += self._logprob_k(trees, ctmc2, k)
        logp = plp - clp + logirate - lograte
        if logp >= 0.0 or np.log(np.random.uniform(0.0, 1.0) + 1E-10) < logp:
            # accept
            self.ctmcs[k] = ctmc2
            return True
        else:
            # reject
            return False

    def sample_ctmc_hmc(self, trees, k):
        # for TwoStateCTMC
        def U(logParams):
            params = np.exp(logParams)
            ctmc = TwoStateCTMC(params)
            ll = -params.sum() / self.ctmc_scale # log Gamma(x|shape=1.0, scale)
            ll += self._logprob_k(trees, ctmc, k)
            return -ll

        def gradU(logParams):
            params = np.exp(logParams)
            ctmc = TwoStateCTMC(params)
            grad = np.ones(2) / -self.ctmc_scale

            def _grad_logprob_branch_k(parent, node, ctmc, k):
                time = parent.date - node.date
                tensor = ctmc.get_transition_grad_logprobs(time)
                _grad = tensor[:, parent.state[k].astype(np.int32), node.state[k].astype(np.int32)]
                for child in node.children:
                    _grad += _grad_logprob_branch_k(node, child, ctmc, k)
                return _grad

            stationary_grad = ctmc.get_stationary_grad_logprobs()
            for tree in trees:
                grad += stationary_grad[:,tree.state[k].astype(np.int32)]
                for child in tree.children:
                    grad += _grad_logprob_branch_k(tree, child, ctmc, k)
            return -grad * params

        params = self.ctmcs[k].get_params()
        accepted, logParams = hmc(U, gradU, self.HMC_EPSILON, self.HMC_L, np.log(params))
        if accepted:
            self.ctmcs[k].set_params(np.exp(logParams))
            return True
        else:
            return False

    def logprob(self, trees):
        ll = 0.0
        for node in self.calibrated_nodes:
            ll += node.prior.logprob(node.date)

        def _logprob_branch_rec(parent, node):
            _ll = 0.0
            time = parent.date - node.date
            _ll += self._logprob_branch_each(time, parent, node, self.ctmcs)
            for child in node.children:
                _ll += _logprob_branch_rec(node, child)
            return _ll

        stationary_probs = []
        for k in range(self.K):
            ll += self.ctmcs[k].get_param_logprob()
            stationary_probs.append(np.log(np.maximum(self.ctmcs[k].get_stationary_probs(), 1E-10)))
        for tree in trees:
            for k in range(self.K):
                # prob. of generating the root state
                ll += stationary_probs[k][tree.state[k].astype(np.int32)]
            for child in tree.children:
                ll += _logprob_branch_rec(tree, child)
        return ll

    def _logprob_branch_each(self, time, parent, node, ctmcs):
        ll = 0.0
        for k in range(self.K):
            mat = ctmcs[k].get_transition_probs(time)
            prob = mat[parent.state[k].astype(np.int32), node.state[k].astype(np.int32)]
            ll += np.log(prob + 1E-10)
        return ll
