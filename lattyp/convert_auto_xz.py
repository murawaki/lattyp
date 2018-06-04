# -*- coding: utf-8 -*-
#
# sample merger for the MDA model
#
import sys
import json
import numpy as np
from collections import Counter
from argparse import ArgumentParser

from json_utils import load_json_file, load_json_stream

def main():
    parser = ArgumentParser()
    parser.add_argument("--input", metavar="FILE", default=None)
    parser.add_argument("--burnin", metavar="INT", type=int, default=0,
                        help="# of burn-in iterations")
    parser.add_argument("--interval", metavar="INT", type=int, default=1,
                        help="pick up one per # samples")
    parser.add_argument("--update", action="store_true", default=False,
                        help="update features (for MVI)")
    parser.add_argument("langs", metavar="LANG", default=None)
    parser.add_argument("flist", metavar="FLIST", default=None)
    args = parser.parse_args()

    flist = load_json_file(args.flist)
    langs = list(load_json_stream(open(args.langs)))
    P = len(flist)
    L = len(langs)

    count = 0
    xfreq = []
    for l in range(L):
        xfreq.append([None] * P)
    zfreq = None

    if args.input is None or args.input == "-":
        f = sys.stdin
    elif args.input.endswith(".bz2"):
        import bz2
        f = bz2.open(args.input, "r")
    else:
        f = open(args.input, "r")

    for langdat in load_json_stream(f):
        sys.stderr.write("+")
        sys.stderr.flush()
        if langdat["iter"] >= args.burnin and langdat["iter"] % args.interval == 0:
            count += 1
            if zfreq is None:
                zfreq = np.zeros((L, len(langdat["z"])), dtype=np.int32)
            zfreq += np.array(langdat["z"]).T
            for l in range(L):
                for p in range(P):
                    v = langdat["x"][l][p]
                    if xfreq[l][p] is None:
                        if flist[p]["type"] == "bin":
                            xfreq[l][p] = [0, 0]
                        elif flist[p]["type"] == "cat":
                            xfreq[l][p] = [0] * flist[p]["size"]
                        elif flist[p]["type"] == "count":
                            xfreq[l][p] = Counter()
                    xfreq[l][p][v] += 1
    if args.input is not None and args.input != "-":
        f.close()
    sys.stderr.write("\n")

    for p in range(P):
        if flist[p]["type"] == "count":
            for l, lang in enumerate(langs):
                maxv = max(xfreq[l][p].keys())
                vlist = [0] * (maxv + 1)
                for k, v in xfreq[l][p].items():
                    vlist[k] = v
                xfreq[l][p] = vlist
    for l, lang in enumerate(langs):
        lang["count"] = count
        lang["xfreq"] = xfreq[l]
        lang["zfreq"] = zfreq[l].tolist()
        if args.update:
            for p, fnode in enumerate(flist):
                v = int(np.argmax(np.array(lang["xfreq"][p])))
                lang["catvect_filled"][p] = v
                lang["annotation"]["features_filled"][fnode["annotation"]["name"]] = v
        sys.stdout.write("{}\n".format(json.dumps(lang)))

if __name__ == "__main__":
    main()
