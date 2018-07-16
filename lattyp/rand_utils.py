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


def slice_sampler1d(logfunc, x, min_x=-np.inf, max_x=np.inf, w=0.0, nsamples=1, max_nfeval=200):
    def F(x):
        if min_x < x < max_x:
            nfeval[0] += 1
            if nfeval[0] > max_nfeval:
                raise Exception
            fx = logfunc(x)
            # assert(np.isfinite(fx))
            return fx
        else:
            return -np.inf

    assert(np.isfinite(x))
    if w <= 0.0:
        if min_x > -np.inf and max_x < np.inf:
            w = (max_x - min_x) / 8.0
        else:
            w = max(abs(x) / 4.0, 0.1)
            # _w = -x0 if x0 < 0.0 else x0
            # w = max(_w / 2.0, 1e-7)

    nfeval = [0] # use list since nonlocal variables cannot be modified in Python2
    logfx = F(x)
    for sample in range(0, nsamples):
        logy = logfx + np.log(np.random.random() + 1e-100)
        assert(np.isfinite(logy))

        xl = max(min_x + 1E-10, x - w * np.random.random())
        logfxl = F(xl)
        xr = min(max_x - 1E-10, xl + w)
        logfxr = F(xr)

        left_limit, right_limit = False, False
        while logy < logfxl or logy < logfxr: # doubling
            if left_limit:
                if right_limit:
                    break
                else:
                    is_left = False
            else:
                if right_limit:
                    is_left = True
                else:
                    is_left = np.random.random() < 0.5
            if is_left:
                xl -= xr - xl
                if xl < min_x:
                    xl = min_x + 1E-10
                    left_limit = True
                logfxl = F(xl)
            else:
                xr += xr - xl
                if xr > max_x:
                    xr = max_x - 1E-10
                    right_limit = True
                logfxr = F(xr);

        xl1 = xl
        xr1 = xr
        while True: # shringking
            x1 = xl1 + np.random.random() * (xr1 - xl1);
            if logy < F(x1):
                xl2 = xl
                xr2 = xr
                d = False
                iflag = False
                acceptable = None
                while xr2 - xl2 > 1.1 * w:
                    xm = (xl2 + xr2) / 2.0
                    if (x < xm <= x1) or (x >= xm > x1):
                        d = True
                    if x1 < xm:
                        xr2 = xm
                    else:
                        xl2 = xm
                    if d and logy >= F(xl2) and logy >= F(xr2):
                        iflag = True
                        acceptable = False
                        break
                if iflag and acceptable is False:
                    break
                x = x1
                acceptable = True
                break
            else:
                acceptable = True
                break
            if acceptable:
                break
            else:
                if x1 < x:
                    xl1 = x1
                else :
                    xr1 = x1
        w = (4.0 * w + (xr1 - xl1)) / 5.0
    return x
