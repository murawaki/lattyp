# -*- coding: utf-8 -*-
#
# latent representations of typological features
#
import sys, os
import shutil
import math
import numpy as np
import random
import pickle
from argparse import ArgumentParser

from json_utils import load_json_file, load_json_stream
from mda import MatrixDecompositionAutologistic


def create_mat(langs, flist):
    P = len(flist)
    mat = np.zeros((len(langs), P), dtype=np.int32)
    mvs = np.zeros((len(langs), P), dtype=np.bool_)
    for l, lang in enumerate(langs):
        for p, (v, v2) in enumerate(zip(lang["catvect_filled"], lang["catvect"])):
            mat[l,p] = v
            if v2 < 0:
                mvs[l,p] = True
    return mat, mvs

def create_cvlist(langs):
    cvlist = []
    for l, lang in enumerate(langs):
        for p, (v1, v2) in enumerate(zip(lang["catvect"], lang["catvect_orig"])):
            if v1 == -1 and v2 != -1:
               cvlist.append((l, p, v2))
    return cvlist

def create_vnet(langs):
    vnet = []
    vgroups = {}
    for l, lang in enumerate(langs):
        if lang["phylo"] not in vgroups:
            vgroups[lang["phylo"]] = []
        vgroups[lang["phylo"]].append(l)
    for l, lang in enumerate(langs):
        vnet.append(np.array([l2 for l2 in vgroups[lang["phylo"]] if not l2 == l], dtype=np.int32))
    return vnet

def create_hnet(langs):
    def _distance(x1, y1, x2, y2):
        """
        Calculate a distance between 2 points
        from whose longitude and latitude
        """
        A = 6378137.0
        B = 6356752.314140
        dy = math.radians(y1 - y2)
        dx = math.radians(x1 - x2)

        if(dx < -math.pi):
            dx += 2*math.pi
        if(dx > math.pi):
            dx -= 2*math.pi

        my = math.radians((y1 + y2) / 2)
        E2 = (A**2 - B**2) / A**2
        Mnum = A * (1 - E2)
        w = math.sqrt(1 - E2 * math.sin(my)**2)
        m = Mnum / w**3
        n = A / w
        return math.sqrt((dy * m)**2 + (dx * n * math.cos(my))**2)

    hnet = []
    for i, lang1 in enumerate(langs):
        hvect = []
        for j, lang2 in enumerate(langs):
            if i == j: continue
            distance = _distance(float(lang1['longitude']),
                                 float(lang1['latitude']),
                                 float(lang2['longitude']),
                                 float(lang2['latitude']))
            # in 1000km
            DISTANCE_THRESHOLD = 1000000
            if distance <= DISTANCE_THRESHOLD:
               hvect.append(j)
        hnet.append(np.array(hvect, dtype=np.int32))
    return hnet

def eval_cvlist(mda):
    cor = 0
    for l, p, v in mda.cvlist:
        if mda.mat[l,p] == v:
            cor += 1
    cv_result = { "cor": cor, "total": len(mda.cvlist) }
    sys.stderr.write("\tcv\t{}\t({} / {})\n".format(float(cor) / len(mda.cvlist), cor, len(mda.cvlist)))
    return cv_result

def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", metavar="INT", type=int, default=None,
                        help="random seed")
    parser.add_argument("--bias", action="store_true", default=False,
                        help="bias term in Z")
    parser.add_argument("--only_alphas", action="store_true", default=False,
                        help="autologistic: ignore v and h")
    parser.add_argument("--drop_vs", action="store_true", default=False,
                        help="autologistic: ignore h")
    parser.add_argument("--drop_hs", action="store_true", default=False,
                        help="autologistic: ignore v")
    parser.add_argument("-i", "--iter", dest="_iter", metavar="INT", type=int, default=1000,
                        help="# of iterations")
    parser.add_argument("--save_interval", metavar="INT", type=int, default=-1,
                        help="save interval")
    parser.add_argument("--K", metavar="INT", type=int, default=100,
                        help="K")
    parser.add_argument('--norm_sigma', type=float, default=5.0,
                        help='standard deviation of Gaussian prior for u')
    parser.add_argument('--gamma_shape', type=float, default=1.0,
                        help='shape of Gamma prior for v and h')
    parser.add_argument('--gamma_scale', type=float, default=0.001,
                        help='scale of Gamma prior for v and h')
    parser.add_argument("--maxanneal", metavar="INT", type=int, default=0)
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
    args = parser.parse_args()
    sys.stderr.write("args\t{}\n".format(args))

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    flist = load_json_file(args.flist)

    offset = 0
    if args.resume_if:
        if os.path.isfile(args.output + ".current"):
            args.resume = args.output + ".current"
        elif os.path.isfile(args.output + ".best"):
            args.resume = args.output + ".best"
    if args.resume:
        sys.stderr.write("loading model from {}\n".format(args.resume))
        spec = pickle.load(open(args.resume, "rb"))
        mda = spec["model"]
        sys.stderr.write("iter {}: {}\n".format(spec["iter"] + 1, spec["ll"]))
        if args.cv:
            eval_cvlist(mda)
        offset = spec["iter"] + 1
    else:
        langs = list(load_json_stream(open(args.langs)))
        mat, mvs = create_mat(langs, flist)

        sys.stderr.write("building vnet\n")
        vnet = create_vnet(langs)
        sys.stderr.write("building hnet\n")
        hnet = create_hnet(langs)
        mda = MatrixDecompositionAutologistic(mat, flist,
                                              vnet=vnet, hnet=hnet,
                                              K=args.K, mvs=mvs,
                                              bias=args.bias,
                                              only_alphas=args.only_alphas,
                                              drop_vs=args.drop_vs,
                                              drop_hs=args.drop_hs,
                                              norm_sigma=args.norm_sigma,
                                              gamma_shape=args.gamma_shape,
                                              gamma_scale=args.gamma_scale)
        if args.cv:
            mda.cvlist = create_cvlist(langs)
        mda.init_with_clusters()
        sys.stderr.write("iter 0: {}\n".format(mda.calc_loglikelihood()))
        if args.cv:
            eval_cvlist(mda)
    ll_max = -np.inf
    for _iter in range(offset, args._iter):
        mda.sample(_iter=_iter, maxanneal=args.maxanneal)
        ll = mda.calc_loglikelihood()
        sys.stderr.write("iter {}: {}\n".format(_iter + 1, ll))
        sys.stderr.flush()
        if args.cv:
            cv_result = eval_cvlist(mda)
            sys.stderr.flush()
        if args.save_interval >= 0 and (_iter + 1) % args.save_interval == 0:
            with open(args.output + ".{}".format(_iter), "wb") as f:
                obj = { "model": mda, "iter": _iter, "ll": ll }
                if args.cv:
                    obj["cv_result"] = cv_result
                pickle.dump(obj, f)
        if args.output is not None:
            with open(args.output + ".current", "wb") as f:
                obj = { "model": mda, "iter": _iter, "ll": ll }
                if args.cv:
                    obj["cv_result"] = cv_result
                pickle.dump(obj, f)
        if ll > ll_max:
            ll_max = ll
            shutil.copyfile(args.output + ".current", args.output + ".best")
    if args.output is not None:
        with open(args.output + ".final", "wb") as f:
            obj = { "model": mda, "iter": _iter, "ll": ll }
            if args.cv:
                obj["cv_result"] = cv_result
            pickle.dump(obj, f)

if __name__ == "__main__":
    main()
