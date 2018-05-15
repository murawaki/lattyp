# -*- coding: utf-8 -*-
import sys, os
import json

from json_utils import load_json_file, load_json_stream

def main(orig, src, fpath, dst):
    flist = load_json_file(fpath)

    with open(src) as fin:
        fin.readline() # ignore the header
        with open(dst, 'w') as fout:
            for lang, line in zip(load_json_stream(open(orig)), fin):
                line = line.rstrip()
                a = line.split("\t")
                a.pop(0) # lang id
                catvect = list(map(lambda x: int(x), a))
                lang["catvect_filled"] = catvect
                lang["annotation"]["features_filled"] = {}
                for fid, v in enumerate(catvect):
                    name = flist[fid]["annotation"]["name"]
                    lang["annotation"]["features_filled"][name] = v
                    assert(name not in lang["annotation"]["features"] or lang["annotation"]["features"][name] == v)
                fout.write("%s\n" % json.dumps(lang))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
