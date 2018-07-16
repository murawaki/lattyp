#!/bin/env python
# -*- coding: utf-8 -*-
#
# Merge Glottolog tree and WALS/autotyp languages
#
import sys
import os
import re
from pickle import load, dump
from argparse import ArgumentParser

from newick_tree import Node

from json_utils import load_json_file, load_json_stream


def remove_unclassified(node):
    # lower nodes take priority over attachment
    children2 = []
    for child in node.children:
        if "Unclassified" in child.name:
            sys.stderr.write("remove unclassified node: {}\n".format(child.name))
        else:
            children2.append(child)
            remove_unclassified(child)
    node.children = children2

def attach_lang(node, glotto_code2lang):
    # lower nodes take priority over attachment
    for child in node.children:
        attach_lang(child, glotto_code2lang)

    m = attach_lang.code_re.search(node.name)
    if m:
        code = m.group(1)
        node.glottocode = code
        if code in glotto_code2lang:
            if len(glotto_code2lang[code]) > 1:
                # create child nodes under the target node and attach languages to them
                sys.stderr.write("creating child nodes because glottocode {} is shared by {} langs\n".format(code, len(glotto_code2lang[code])))
                for i, lang in enumerate(glotto_code2lang[code]):
                    child = Node("%d_%d" % (node._id, i))
                    child.name = node.name
                    child.lang = lang
                    child.is_virtual_child = True
                    node.is_virtual_parent = True
                    child.parent = node
                    node.children.append(child)
                    lang["used"] = True
            else:
                lang = glotto_code2lang[code][0]
                if "used" in lang and lang["used"] == True:
                    sys.stderr.write("glottocode: {} already appeared\n".format(code))
                    node.is_virtual_parent = True
                else:
                    node.lang = lang
                    lang["used"] = True
            # if len(node.children) > 0:
            #     sys.stderr.write("NOTE: attach data to an internal node: %s\n" % node.name)
    else:
        sys.stderr.write("node without glottocode: {}\n".format(node.name))

attach_lang.code_re = re.compile(r" \[([a-z0-9]{4}[0-9]{4})\]")

def shrink_tree(node, registered_codes):
    children2 = []
    for child in node.children:
        c = shrink_tree(child, registered_codes)
        if c > 0:
            children2.append(child)
    node.children = children2
    if len(node.children) <= 0:
        if hasattr(node, "lang"):
            # becomes a leaf
            return 1
        else:
            # sys.stderr.write("node %s removed (non-leaf without children)\n" % node.name)
            return 0
    else:
        if hasattr(node, "lang"):
            sys.stderr.write("creating a node for {} because the intermediate node is attached\n".format(node.name))
            child = Node("%d_c" % (node._id))
            child.name = node.name
            child.lang = node.lang
            child.is_virtual_child = True
            node.is_virtual_parent = True
            del node.lang
            child.parent = node
            node.children.append(child)
        else:
            if len(node.children) > 1:
                # valid intermediate node
                if len(node.children) > 2:
                    sys.stderr.write(u"\tnode {} has {} children ({} nodes need to be added to bifurcate)\n".format(node.name, len(node.children), len(node.children) - 2))
            else:
                # skip intermediate node by overriding the node
                # sys.stderr.write("skipping node %s (non-leaf with one child)\n" % node.name)

                # keep glottocode specified by a prior
                if node.glottocode in registered_codes:
                    if len(node.children[0].children) <= 0:
                        sys.stderr.write("drop a prior since the node has only one child\t{}\n".format(node.glottocode))
                    else:
                        sys.stderr.write("keep glottocode associated with a prior\t{}\n".format(node.glottocode))
                        if node.children[0].glottocode in registered_codes:
                            sys.stderr.write("WARNING: parent-child associated with priors\t{}\t{}\n".format(node.glottocode, node.children[0].glottocode))
                        else:
                            # node's properties will be replaced by child's
                            node.children[0].glottocode = node.glottocode
                parent = None
                if hasattr(node, "parent"):
                    parent = node.parent
                if hasattr(node, "is_virtual_parent"):
                    delattr(node, "is_virtual_parent")
                for k, v in node.children[0].__dict__.items():
                    if k not in ("parent", "is_virtual_parent"):
                        setattr(node, k, v)
                for child in node.children:
                    child.parent = node
                if parent:
                    node.parent = parent
        return 1

def get_feature_coverage(tree):
    countlist = []
    def _get_feature_coverage(node):
        if hasattr(node, "lang"):
            catvect = node.lang["catvect"]
            if len(countlist) <= 0:
                for i in xrange(len(catvect)):
                    countlist.append(0)
            for i, v in enumerate(catvect):
                if v >= 0:
                    countlist[i] += 1
        for child in node.children:
            _get_feature_coverage(child)
    _get_feature_coverage(tree)
    c = len(filter(lambda x: x > 0, countlist))
    return float(c) / len(countlist)

def set_date(node):
    if hasattr(node, "lang"):
        assert(len(node.children) == 0)
        # TODO: generalization
        node.date = 0.0
        node.is_date_frozen = True
        node.is_state_frozen = True
    else:
        node.is_date_frozen = False
        node.is_state_frozen = False
        for child in node.children:
            set_date(child)

def main():
    parser = ArgumentParser()
    parser.add_argument("--lthres", dest="lthres", metavar="FLOAT", type=float, default=0.0,
                        help="eliminate trees with higher rate of missing values [0,1]")
    parser.add_argument("--npriors", metavar="NODE_PRIORS", default=None, help="priors for nodes (json)")
    parser.add_argument("langs", metavar="LANG", default=None)
    parser.add_argument("tree", metavar="TREE", default=None)
    parser.add_argument("out", metavar="OUTPUT", default=None)
    args = parser.parse_args()

    trees = load(open(args.tree, 'rb'))
    langs = list(load_json_stream(open(args.langs)))

    # a glottocode corresponds to one or more WALS languages
    glotto_code2lang = {}
    for lang in langs:
        if "glottocode" in lang and lang["glottocode"]:
            if lang["glottocode"] in glotto_code2lang:
                glotto_code2lang[lang["glottocode"]].append(lang)
            else:
                glotto_code2lang[lang["glottocode"]] = [lang]
        else:
            sys.stderr.write(u"dropping lang without glottocode: {}\n".format(lang["annotation"]["name"]))

    # try to keep nodes specified by priors
    registered_codes = {}
    if args.npriors is not None:
        prior_specs = load_json_file(args.npriors)
        for spec in prior_specs:
            registered_codes[spec["glottocode"]] = True

    for tree in trees:
        remove_unclassified(tree)
    for tree in trees:
        attach_lang(tree, glotto_code2lang)

    ucount, total = 0, 0
    for code, langlist in glotto_code2lang.items():
        for lang in langlist:
            total += 1
            if "used" in lang and lang["used"] == True:
                ucount += 1
                del lang["used"]
            else:
                sys.stderr.write("glottocode never appeared in trees: {}\n".format(code))
    sys.stderr.write("{} (out of {}) languages in the trees\n".format(ucount, total))
        
    trees2 = []
    for tree in trees:
        c = shrink_tree(tree, registered_codes)
        if c > 0:
            trees2.append(tree)
    sys.stderr.write("# of trees: {} -> {}\n".format(len(trees), len(trees2)))
    trees = trees2

    if args.lthres > 0.0:
        trees2 = []
        for tree in trees:
            c = get_feature_coverage(tree)
            if c >= args.lthres:
                trees2.append(tree)
        sys.stderr.write("# of trees ({} thres): {} -> {}\n".format(args.lthres, len(trees), len(trees2)))
        trees = trees2

    isolate = 0
    for tree in trees:
        tree.is_root = True
        tree.parent = None
        if len(tree.children) <= 0:
            tree.is_isolate = True
            isolate += 1
            tree.date = 0.0
            tree.is_date_frozen = True
            tree.is_state_frozen = True
        else:
            tree.is_isolate = False
            tree.is_date_frozen = False
            tree.is_state_frozen = False
            for child in tree.children:
                set_date(child)
    sys.stderr.write("# of isolate trees: {} (out of {})\n".format(isolate, len(trees)))
            
    # for tree in trees:
    #     sys.stdout.write(tree.name + "\n")
    with open(args.out, "wb") as f:
        dump(trees, f)


if __name__ == "__main__":
    main()
