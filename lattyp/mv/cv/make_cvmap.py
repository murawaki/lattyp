# -*- coding: utf-8 -*-
#
import sys, os
import json
import random
from argparse import ArgumentParser

from json_utils import load_json_file, load_json_stream


def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--seed", dest="seed", metavar="INT", type=int, default=None,
                        help="random seed")
    parser.add_argument("src", metavar="SOURCE", default=None)
    parser.add_argument("dst", metavar="DESTINATION", default=None)
    parser.add_argument("cvn", metavar="INT", default=None)
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    src, dst = args.src, args.dst
    cvn = int(args.cvn)
    langs = list(load_json_stream(open(src)))

    filled_list = []
    for lang in langs:
        for name, v in lang["annotation"]["features"].items():
            filled_list.append((lang["annotation"]["name"], name))
    random.shuffle(filled_list)

    # N-fold cross-validation
    cell_size = len(filled_list) // cvn
    cell_size2 = len(filled_list) % cvn

    cvmap = [[] for i in range(cvn)]
    for i in range(cvn):
        cell_start = cell_size * i + min(i, cell_size2)
        cell_len = cell_size + (i < cell_size2)
        for j in range(cell_start, cell_start + cell_len):
            cvmap[i].append(filled_list[j])

    with open(dst, 'w') as f:
        f.write(json.dumps(cvmap))


if __name__ == "__main__":
    main()
