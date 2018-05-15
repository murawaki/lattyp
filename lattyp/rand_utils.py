# -*- coding: utf-8 -*-
import numpy as np

def rand_partition_log(log_list):
    base = max(log_list)
    prob_list = [np.exp(l - base) for l in log_list]
    return rand_partition(prob_list)

def rand_partition(prob_list):
    s = sum(prob_list)
    r = np.random.uniform(0, s)
    for i in range(0, len(prob_list)):
        r -= prob_list[i]
        if r <= 0.0:
            return i
    return len(prob_list) - 1
