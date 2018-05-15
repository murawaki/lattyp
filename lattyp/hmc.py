# -*- coding: utf-8 -*-
"""
Hamiltonian Monte Carlo
http://www.mcmchandbook.net/HandbookChapter5.pdf
"""
import numpy as np

def hmc(U, gradU, epsilon, L, current_q):
    q = current_q
    p = np.random.normal(size=q.shape)
    current_p = p

    # make a half step for momentum at the beginning
    p = p - epsilon * gradU(q) / 2.0
    # alternate full steps for position and momentum
    for i in range(L):
        # make a full step for the position
        q = q + epsilon * p
        # make a full step for the momentum, except at end of trajectory
        if i < L - 1:
            p = p - epsilon * gradU(q)
    # make a half step for momentum at the end
    p = p - epsilon * gradU(q) / 2.0
    # negate momentum at end of trajectory to make the proposal symmetric
    p = -p

    # evaluate the potential and kinetic energies at start and end of trajectory
    currentU = U(current_q)
    currentK = (current_p * current_p).sum() / 2.0
    proposedU = U(q)
    proposedK = (p * p).sum() / 2.0

    # print "hmc ", np.exp(currentU - proposedU + currentK - proposedK)
    # accept or reject the state at end of trajectory,
    logp = currentU - proposedU + currentK - proposedK
    if logp >= 0 or np.log(np.random.uniform(0.0, 1.0)) < logp:
        # accept
        return (True, q)
    else:
        # reject
        return (False, current_q)
