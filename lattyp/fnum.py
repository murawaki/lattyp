# -*- coding: utf-8 -*-
import sys
import json

from json_utils import load_json_file, load_json_stream

flist = load_json_file(sys.argv[1])
print(len(flist))
