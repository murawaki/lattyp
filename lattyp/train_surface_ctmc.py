#!/bin/env python
# -*- coding: utf-8 -*-
#
# train a CTMC model for a surface feature
#
import sys
import os
import codecs
import re
import six
# from six.moves.cPickle import load, dump
import pickle
from argparse import ArgumentParser
import numpy as np
import random

from json_utils import load_json_file, load_json_stream
from phylo import CTMC_Sampler
from rand_utils import rand_partition_log, slice_sampler1d
from priors import * # create_node_priors

def update_state(node, wals_ids):
    node.oldstate = node.state
    delattr(node, "state")
    if hasattr(node, "lang"):
        state = []
        for wals_id in wals_ids:
            state.append(node.lang["annotation"]["features_filled"][wals_id]) # **HACK**
        node.state = np.array(state, dtype=np.int32)
    for child in node.children:
        update_state(child, wals_ids=wals_ids)

def main():
    # sys.stderr = codecs.getwriter("utf-8")(sys.stderr)

    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", dest="seed", metavar="INT", type=int, default=None,
                        help="random seed")
    # parser.add_argument("--sidx", metavar="IDX", type=int, default=0,
    #                     help="i-th sample of leaf states (-1: last sample)")
    # parser.add_argument("--npriors", metavar="NODE_PRIORS", default=None, help="priors for nodes (json)")
    parser.add_argument("-i", "--iter", dest="_iter", metavar="INT", type=int, default=1000,
                        help="# of iterations")
    # parser.add_argument("--resume_if", action="store_true", default=False,
    #                     help="resume training if the output exists")
    parser.add_argument("model", metavar="FILE", default=None,
                        help="resume training from model dump")
    parser.add_argument("flist", metavar="FLIST", default=None)
    # parser.add_argument("trees", metavar="TREES", default=None, help="merged trees (pkl)")
    parser.add_argument("langs", metavar="LANG", default=None) # **HACK**
    # parser.add_argument("samples", metavar="SAMPLES", default=None, help="parameter states (json stream)")
    parser.add_argument("out", metavar="OUT", default=None, help="out (pkl)")
    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    # if args.resume_if:
    #     if os.path.isfile(args.out + ".current"):
    #         args.resume = args.out + ".current"
    #     elif os.path.isfile(args.out + ".best"):
    #         args.resume = args.out + ".best"
    # if args.resume:
    flist = load_json_file(args.flist)
    for fid, fnode in enumerate(flist):
        if fnode["annotation"]["fullname"] == "81A Order of Subject, Object and Verb":
            wals_id = fnode["annotation"]["name"]
            T = fnode["size"] # len(fnode["vid2label"])
            break
    # "label2vid": {"1 SOV": 0, "2 SVO": 1, "3 VSO": 2, "6 OSV": 5, "4 VOS": 3, "5 OVS": 4, "7 No dominant order": 6}
    # fval = 0
    # j_start, T = ibp.bmap(fid)

    sys.stderr.write("loading model from %s\n" % args.model)
    spec = pickle.load(open(args.model, "rb"), encoding="latin-1")
    trees = spec["trees"]

    
    # HACK
    from train_bin_ctmc import register_node
    langs = list(load_json_stream(open(args.langs, "r")))
    idx2id = {}
    for i, lang in enumerate(langs):
        if "glottocode" in lang:
            idx2id[i] = lang["glottocode"] + ":" + lang["annotation"]["name"]
    id2node = {}
    glottocode2node = {}
    for tree in trees:
        register_node(tree, id2node, glottocode2node)
    for i, lang in enumerate(langs):
        if i in idx2id:
            _id = idx2id[i]
            if _id in id2node:
                node = id2node[_id]
                node.lang = lang
    
    # sampler = spec["sampler"]
    # if "logprob" not in spec:
    #     logprob = sampler.logprob(trees)
    # else:
    #     logprob = spec["logprob"]
    # sys.stderr.write("iter {}\t{}\n".format(spec["iter"], logprob))
    # _start = spec["iter"] + 1
    # else:
    _start = 1
    #     trees = load(open(args.trees, 'rb'))
    #     # trees2 = []
    #     # for tree in trees:
    #     #     if tree.is_isolate is False:
    #     #         trees2.append(tree)
    #     # trees = trees2
    #     sys.stderr.write("{} trees\n".format(len(trees)))

    #     node_priors = None
    #     if args.npriors is not None:
    #         prior_specs = load_json_file(args.npriors)
    #         node_priors = create_node_priors(prior_specs)
    
    #     langs = list(load_json_stream(open(args.langs)))
    #     with open(args.samples, 'r') as f:
    #         for i, sample in enumerate(load_json_stream(f)):
    #             if i == args.sidx:
    #                 break
    for tree in trees:
        update_state(tree, [wals_id])

    sampler = CTMC_Sampler(1, states=[T], ctmc_scale=0.00005)
    sampler.init_trees(trees, sample_dates=False)
    sys.stderr.write("iter 0\t{}\n".format(sampler.logprob(trees)))
    for _iter in six.moves.range(_start, args._iter + 1):
        sampler.sample(_iter=_iter)
        logprob = sampler.logprob(trees)
        sys.stderr.write("iter {}\t{}\n".format(_iter, logprob))
        with open(args.out + ".current", "wb") as f:
            pickle.dump({ "sampler": sampler, "trees": trees, "iter": _iter, "logprob": logprob }, f)
    with open(args.out + ".final", "wb") as f:
        pickle.dump({ "sampler": sampler, "trees": trees, "iter": _iter, "logprob": logprob }, f)

if __name__ == "__main__":
    main()
