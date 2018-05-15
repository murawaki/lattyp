# -*- coding: utf-8 -*-
#
# output TSV file for imput_mca.r
#
import sys, os
import json
from collections import defaultdict

from json_utils import load_json_file, load_json_stream


def main(src, fpath, dst):
    flist = load_json_file(fpath)
    langs = list(load_json_stream(open(src, "r")))

    counts = {}
    for fstruct in flist:
        if fstruct["type"] == "count":
            counts[fstruct["fid"]] = fstruct

    with open(dst, 'w') as f:
        rv = "\t".join([fstruct["annotation"]["name"] for fstruct in flist])
        f.write(rv + "\n")

        for i, lang in enumerate(langs):
            catvect = list(lang["catvect"])
            f.write("L{}\t{}\n".format(i, "\t".join(map(lambda x: str(x) if x >= 0 else "NA", catvect))))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
