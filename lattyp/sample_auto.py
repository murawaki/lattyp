# -*- coding: utf-8 -*-
#
# sample x and z for languages while keeping W fixed
#
import sys
import json
import numpy as np
import random
import pickle
from argparse import ArgumentParser

from json_utils import load_json_file, load_json_stream
from mda import MatrixDecompositionAutologistic

def dumps(mda, _iter):
    return {
        "iter": _iter,
        "x": mda.mat.tolist(),
        "z": mda.zmat.tolist(),
        "v": mda.vks.tolist(),
        "h": mda.hks.tolist(),
        "a": mda.alphas.tolist(),
    }

def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", metavar="INT", type=int, default=None,
                        help="random seed")
    parser.add_argument("-i", "--iter", dest="_iter", metavar="INT", type=int, default=10,
                        help="# of iterations")
    parser.add_argument("--a_repeat", dest="a_repeat", metavar="INT", type=int, default=1)
    parser.add_argument("model", metavar="MODEL", default=None)
    parser.add_argument("output", metavar="OUTPUT", default=None)
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    spec = pickle.load(open(args.model, "rb"))
    mda = spec["model"]

    f = sys.stdout if args.output == "-" else open(args.output, "w")
    sys.stderr.write("iter 0\n")
    f.write("%s\n" % json.dumps(dumps(mda, 0)))
    mda.init_tasks(a_repeat=args.a_repeat, sample_w=False)
    for _iter in range(args._iter - 1): # already have iter 0
        sys.stderr.write("iter {}\n".format(_iter + 1))
        mda.sample(_iter=_iter)
        f.write("%s\n" % json.dumps(dumps(mda, _iter + 1)))
        f.flush()
    if not f == sys.stdout:
        f.close()

if __name__ == "__main__":
    main()
