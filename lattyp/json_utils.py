# -*- coding: utf-8 -*-

import sys
import json
from bz2 import BZ2File

def load_json_file(fname):
    if fname.endswith(".bz2"):
        reader = BZ2File(fname)
    else:
        reader = open(fname, "r")
    dat = reader.read()
    return json.loads(dat)

def load_json_stream(f, offset=0, verbose=False):
    # do not apply codecs to f; this is too slow!
    for i, line in enumerate(f):
        if i >= offset:
            yield json.loads(line)
        elif verbose:
            if i % 1 == 10:
                sys.stderr.write("#\n")
