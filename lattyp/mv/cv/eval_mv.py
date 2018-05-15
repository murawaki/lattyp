# -*- coding: utf-8 -*-

import sys, os
import json
import numpy as np
import random 
from collections import defaultdict
from argparse import ArgumentParser

from json_utils import load_json_file, load_json_stream

def eval_mv(filled_langs, langs):
    total = 0
    correct = 0
    for lang, filled_lang in zip(langs, filled_langs):
        for fname, v in lang["annotation"]["features"].items():
            if fname not in filled_lang["annotation"]["features"]:
                total += 1
                if filled_lang["annotation"]["features_filled"][fname] == v:
                    correct += 1
    return (total, correct)

def eval_most_frequent(flist, hidelist, langs):
    # NOTE: most frequent values are based on hidelist, not langs (reference)
    fsize = len(flist)
    fname2histogram = {}
    for hlang in hidelist:
        for fname, v in hlang["annotation"]["features"].items():
            if fname not in fname2histogram:
                fname2histogram[fname] = defaultdict(int)
            fname2histogram[fname][v] += 1
    fname2maxk = {}
    for fname, histogram in fname2histogram.items():
        maxk, maxv = None, -1
        for k, v in histogram.items():
            if v >= maxv:
                maxv = v
                maxk = k
        fname2maxk[fname] = maxk

    total = 0
    correct = 0
    for lang, hlang in zip(langs, hidelist):
        for fname, v in lang["annotation"]["features"].items():
            if fname not in hlang["annotation"]["features"]:
                total += 1
                if fname2maxk[fname] == v:
                    correct += 1
    return (total, correct)

def eval_random(flist, langs):
    total = 0
    correct = 0

    fname2size = {}
    for fstruct in flist:
        if fstruct["type"] == "bin":
            size = 2
        else:
            size = len(fstruct["annotation"]["vid2label"])
        fname2size[fstruct["annotation"]["name"]] = size
    for lang in langs:
        for fname, v in lang["annotation"]["features"].items():
            total += 1
            r = np.random.random_integers(0, fname2size[fstruct["annotation"]["name"]] - 1)
            if v == r:
                correct += 1
    return (total, correct)


def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", dest="seed", metavar="INT", type=int, default=None,
                        help="random seed")
    parser.add_argument("--random", dest="random", action="store_true", default=False)
    parser.add_argument("--freq", dest="most_frequent", action="store_true", default=False)
    parser.add_argument("langs", metavar="LANGS", default=None)
    parser.add_argument("f1", metavar="DUMMY_OR_LANGS_FILLED_OR_LANGS_HIDDEN", default=None)
    parser.add_argument("f2", metavar="FLIST_OR_DUMMY_OR_LANGS_HIDDEN", default=None)
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    langs = list(load_json_stream(open(args.langs)))
    if args.random:
        flist = load_json_file(args.f2)
        total, correct = eval_random(flist, langs)
    elif args.most_frequent:
        hidelist = list(load_json_stream(open(args.f1)))
        flist = load_json_file(args.f2)
        total, correct = eval_most_frequent(flist, hidelist, langs)
    else:
        filled_langs = list(load_json_stream(open(args.f1)))
        total, correct = eval_mv(filled_langs, langs)
    sys.stdout.write("%f\t%d\t%d\n" % (float(correct) / total, correct, total))

    
if __name__ == "__main__":
    main()
