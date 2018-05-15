# -*- coding: utf-8 -*-
#
import sys, os
import json
import copy

from json_utils import load_json_file, load_json_stream


def main(src, fpath, dst, cvmap_file, cvn):
    langs = [lang for lang in load_json_stream(open(src))]
    flist = load_json_file(fpath)
    cvmap = load_json_file(cvmap_file)

    name2lang = {}
    for lang in langs:
        lang["annotation"]["features_orig"] = copy.copy(lang["annotation"]["features"])
        lang["catvect_orig"] = copy.copy(lang["catvect"])
        lang["cv"] = cvn
        name2lang[lang["annotation"]["name"]] = lang

    name2fstruct = {}
    for fstruct in flist:
        name2fstruct[fstruct["annotation"]["name"]] = fstruct

    for lname, fname in cvmap[cvn]:
        lang = name2lang[lname]
        fstruct = name2fstruct[fname]
        lang["catvect"][fstruct["fid"]] = -1
        del lang["annotation"]["features"][fname]

    with open(dst, 'w') as f:
        for lang in langs:
            f.write("%s\n" % json.dumps(lang))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))
