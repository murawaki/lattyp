# -*- coding: utf-8 -*-

import sys, os
import json
import numpy as np
import random 
from collections import defaultdict
from argparse import ArgumentParser
from statsmodels.stats.contingency_tables import mcnemar

from json_utils import load_json_file, load_json_stream

def eval_mv(filled_langs1, filled_langs2, langs):
    mat = np.zeros((2, 2), dtype=np.int32)
    for lang, filled_lang1, filled_lang2 in zip(langs, filled_langs1, filled_langs2):
        for fname, v in lang["annotation"]["features"].items():
            if fname not in filled_lang1["annotation"]["features"]:
                p1, p2 = 0, 0
                if filled_lang1["annotation"]["features_filled"][fname] == v:
                    p1 = 1
                if filled_lang2["annotation"]["features_filled"][fname] == v:
                    p2 = 1
                mat[p1,p2] += 1
            else:
                assert(fname in filled_lang2["annotation"]["features"])
    return mat

def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", dest="seed", metavar="INT", type=int, default=None,
                        help="random seed")
    parser.add_argument("--random", dest="random", action="store_true", default=False)
    parser.add_argument("--freq", dest="most_frequent", action="store_true", default=False)
    parser.add_argument("--cvn", metavar="INT", type=int, default=10)
    parser.add_argument("langs", metavar="LANGS", default=None)
    parser.add_argument("f1", metavar="LANGS2 PREFIX", default=None)
    parser.add_argument("f2", metavar="LANGS2 PREFIX", default=None)
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    langs = list(load_json_stream(open(args.langs)))
    mat = np.zeros((2, 2), dtype=np.int32)
    for cvi in range(args.cvn):
        fp1 = args.f1.format(cvi)
        fp2 = args.f2.format(cvi)
        sys.stderr.write("processsing {} and {}\n".format(fp1, fp2))
        filled_langs1 = list(load_json_stream(open(fp1)))
        filled_langs2 = list(load_json_stream(open(fp2)))
        mat += eval_mv(filled_langs1, filled_langs2, langs)
    print(mat)
    bunch = mcnemar(mat, exact=False)
    print("mcnemar\t{}".format(bunch))
    
if __name__ == "__main__":
    main()
