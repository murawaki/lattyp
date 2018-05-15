# -*- coding: utf-8 -*-
import sys, os
import glob
import json
from collections import Counter

from json_utils import load_json_file, load_json_stream


def main(orig, src, fpath, dst):
    langs = list(load_json_stream(open(orig)))
    flist = load_json_file(fpath)

    for lang in langs:
        lang["counted_features"] = [Counter() for feature in flist]
        lang["annotation"]["features_filled"] = {}
        for fstruct in flist:
            lang["annotation"]["features_filled"][fstruct["annotation"]["name"]] = -1

    for fpath in glob.glob(src + ".*"):
        sys.stderr.write("processing {}\n".format(fpath))
        with open(fpath) as fin:
            fin.readline() # ignore the header
            for lang, l in zip(langs, fin):
                l = l.rstrip()
                a = l.split("\t")
                label = a.pop(0)
                for fid, v in enumerate(a):
                    lang["counted_features"][fid][int(v)] += 1

    for lang in langs:
        binsize = 0
        xfreq = Counter()
        for fid, (fstruct, counts) in enumerate(zip(flist, lang["counted_features"])):
            if fstruct["type"] == "bin":
                size = 2
            else:
                size = len(fstruct["annotation"]["vid2label"])
            name = fstruct["annotation"]["name"]
            maxv, maxvv = -1, -1
            for i in range(size):
                xfreq[binsize+i] += counts[i]
                # if lang["xfreq"][binsize+i] >= maxvv:
                if counts[i] >= maxvv:
                    maxvv = counts[i]
                    maxv = i
                lang["annotation"]["features_filled"][name] = maxv
            binsize += size
        del lang["counted_features"]
        lang["xfreq"] = [xfreq[i] for i in range(binsize)]

    with open(dst, 'w') as fout:
        for lang in langs:
            fout.write("%s\n" % json.dumps(lang))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
