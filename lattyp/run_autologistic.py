# -*- coding: utf-8 -*-
#
# autologistic model for categorical features
#
import sys, os
import shutil
import math
import json
import numpy as np
import random
import pickle
from argparse import ArgumentParser

from json_utils import load_json_file, load_json_stream
from train_mda import create_vnet, create_hnet
from autologistic import CategoricalAutologistic

def create_vec(langs, fstruct, p):
    vec = np.zeros(len(langs), dtype=np.int32)
    mvs = np.zeros(len(langs), dtype=np.bool_)
    for l, lang in enumerate(langs):
        v, v2 = lang["catvect_filled"][p], lang["catvect"][p]
        if fstruct["type"] == "count":
            vec[l] = fstruct["annotation"]["label2vid"][str(v)]
        else:
            vec[l] = v
        if v2 < 0:
            mvs[l] = True
    if fstruct["type"] == "cat":
        size = fstruct["size"]
    elif fstruct["type"] == "bin":
        size = 2
    elif fstruct["type"] == "count":
        size = len(fstruct["annotation"]["vid2label"])
    else:
        raise NotImplementedError
    return vec, mvs, size

def create_cvlist(langs, fstruct, p):
    cvlist = []
    for l, lang in enumerate(langs):
        v1, v2 = lang["catvect"][p], lang["catvect_orig"][p]
        if v1 == -1 and v2 != -1:
            if fstruct["type"] == "count":
                v2 = fstruct["annotation"]["label2vid"][str(v2)]
            cvlist.append((l, v2))
    return cvlist

def eval_cvlist(al, verbose=True):
    cor = 0
    for l, v in al.cvlist:
        if al.vec[l] == v:
            cor += 1
    cv_result = { "cor": cor, "total": len(al.cvlist) }
    if verbose:
        sys.stderr.write("\tcv\t{}\t({} / {})\n".format(float(cor) / len(al.cvlist), cor, len(al.cvlist)))
    return cv_result

def get_result(al, do_cv):
    obj = {
        "vec": al.vec.copy(),
        "v": al.v,
        "h": al.h,
        "a": al.alphas.copy(),
    }
    if do_cv:
        obj["cv"] = eval_cvlist(al, verbose=False)
    return obj

def aggregate_results(results, al, fstruct, do_cv):
    vecs = np.stack(list(map(lambda x: x["vec"], results)))
    vs = np.array(list(map(lambda x: x["v"], results)))
    hs = np.array(list(map(lambda x: x["h"], results)))
    alphass = np.stack(list(map(lambda x: x["a"], results)))
    
    bincounts = []
    for i in range(al.L):
        bincounts.append(np.bincount(vecs[:,i], minlength=al.size))
    bincounts = np.stack(bincounts)
    agg_vec = np.argmax(bincounts, axis=1)
    results = {
        "fid": fstruct["fid"],
        "vec": agg_vec.tolist(),
        "v": vs.mean(),
        "h": hs.mean(),
        "a": alphass.mean(axis=0).tolist(),
        "samplenum": len(results),
        "samples": {
            "vs": vs.tolist(),
            "hs": hs.tolist(),
            "alphass": alphass.tolist(),
        },
    }
    if do_cv:
        cor = 0
        for l, v in al.cvlist:
            if agg_vec[l] == v:
                cor += 1
        results["cv_result"] = { "cor": cor, "total": len(al.cvlist), "cv": al.cvn }
    return results

def main():
    # sys.stderr = codecs.getwriter("utf-8")(sys.stderr)

    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", metavar="INT", type=int, default=None,
                        help="random seed")
    parser.add_argument("--fid", metavar="INT", type=int, default=0)
    parser.add_argument("--only_alphas", action="store_true", default=False,
                        help="autologistic: ignore v and h")
    parser.add_argument("--drop_vs", action="store_true", default=False,
                        help="autologistic: ignore h")
    parser.add_argument("--drop_hs", action="store_true", default=False,
                        help="autologistic: ignore v")
    parser.add_argument("--burnin", metavar="INT", type=int, default=1000,
                        help="# of iterations")
    parser.add_argument("--samples", metavar="INT", type=int, default=500,
                        help="save interval")
    parser.add_argument("--interval", metavar="INT", type=int, default=5,
                        help="sampling interval")
    parser.add_argument("--alpha", metavar="FLOAT", type=float, default=-1.0,
                        help="parameter alpha")
    parser.add_argument("--K", metavar="INT", type=int, default=100,
                        help="K")
    parser.add_argument('--norm_sigma', type=float, default=5.0,
                        help='standard deviation of Gaussian prior for u')
    parser.add_argument('--gamma_shape', type=float, default=1.0,
                        help='shape of Gamma prior for v and h')
    parser.add_argument('--gamma_scale', type=float, default=0.001,
                        help='scale of Gamma prior for v and h')
    parser.add_argument("--cv", action="store_true", default=False,
                        help="some features are intentionally hidden (but kept as \"catvect_orig\")")
    parser.add_argument("--output", dest="output", metavar="FILE", default=None,
                        help="save the model to the specified path")
    parser.add_argument("--resume", metavar="FILE", default=None,
                        help="resume training from model dump")
    parser.add_argument("--resume_if", action="store_true", default=False,
                        help="resume training if the output exists")
    parser.add_argument("langs", metavar="LANG", default=None)
    parser.add_argument("flist", metavar="FLIST", default=None)
    parser.add_argument("aggregated", metavar="FLIST", default=None)
    args = parser.parse_args()
    sys.stderr.write("args\t{}\n".format(args))

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    flist = load_json_file(args.flist)
    fstruct = flist[args.fid]

    # offset = 0
    # if args.resume_if:
    #     if os.path.isfile(args.output + ".current"):
    #         args.resume = args.output + ".current"
    #     elif os.path.isfile(args.output + ".best"):
    #         args.resume = args.output + ".best"
    # if args.resume:
    #     sys.stderr.write("loading model from {}\n".format(args.resume))
    #     spec = pickle.load(open(args.resume, "rb"))
    #     mda = spec["model"]
    #     sys.stderr.write("iter {}\n".format(spec["iter"] + 1))
    #     if args.cv:
    #         eval_cvlist(mda)
    #     offset = spec["iter"] + 1
    # else:
    langs = list(load_json_stream(open(args.langs)))
    vec, mvs, size = create_vec(langs, fstruct, args.fid)

    sys.stderr.write("building vnet\n")
    vnet = create_vnet(langs)
    sys.stderr.write("building hnet\n")
    hnet = create_hnet(langs)
    al = CategoricalAutologistic(vec, size,
                                 vnet=vnet, hnet=hnet,
                                 mvs=mvs,
                                 only_alphas=args.only_alphas,
                                 drop_vs=args.drop_vs,
                                 drop_hs=args.drop_hs,
                                 norm_sigma=args.norm_sigma,
                                 gamma_shape=args.gamma_shape,
                                 gamma_scale=args.gamma_scale)
    if args.cv:
        al.cvlist = create_cvlist(langs, fstruct, args.fid)
    # al.init_with_clusters()
    sys.stderr.write("iter 0\n")
    if args.cv:
        eval_cvlist(al)
        al.cvn = langs[0]["cv"]
    # ll_max = -np.inf
    offset = 0
    for _iter in range(args.burnin):
        al.sample()
        offset += 1
        # ll = mda.calc_loglikelihood()
        sys.stderr.write("iter {}\n".format(offset))
        sys.stderr.flush()
        if args.cv:
            cv_result = eval_cvlist(al)
            sys.stderr.flush()
    if args.output is not None:
        with open(args.output, "wb") as f:
            obj = { "model": al, "iter": offset }
            if args.cv:
                obj["cv_result"] = cv_result
            pickle.dump(obj, f)
    results = []
    results.append(get_result(al, args.cv))
    while len(results) < args.samples:
        for _iter in range(args.interval):
            al.sample()
            offset += 1
            sys.stderr.write("iter {}\n".format(offset))
            sys.stderr.flush()
            if args.cv:
                cv_result = eval_cvlist(al)
                sys.stderr.flush()
        results.append(get_result(al, args.cv))            
    if args.aggregated == "-":
        f = sys.stdout
    else:
        f = open(args.aggregated, "w")
    aggregated = aggregate_results(results, al, fstruct, args.cv)
    f.write("%s\n" % json.dumps(aggregated))

if __name__ == "__main__":
    main()
