#!/bin/env python
# -*- coding: utf-8 -*-
#
# simple parser of NEXUS annotated trees
#
import sys
import os
import re
import pickle

def label_clades(node):
    clade_dict = {}
    def _label_clades_main(node):
        label_list = []
        for cnode in node.children:
            label_list += _label_clades_main(cnode)
        label_list.sort()
        node.clade = ":".join(label_list)        
        # if hasattr(node, "left"):
        #     label_list = _label_clades_main(node.left) + _label_clades_main(node.right)
        #     label_list.sort()
        #     node.clade = ":".join(label_list)
        if not hasattr(node, "parent"): # root
            node.clade = "ROOT"
            label_list = [node.clade]
        elif hasattr(node, "name"):
            # named nodes including leaves
            node.clade = node.name
            label_list = [node.clade]
        clade_dict[node.clade] = node
        return label_list
    _label_clades_main(node)
    return clade_dict

class Node(object):
    def __init__(self, _id):
        self._id = _id
        self.children = []

class TreeParser(object):
    START = 1
    END = 2
    EOS = 3
    NODE = 4
    COMMA = 5
    ANNOTATION = 6
    BRANCH = 7
    taxa_re = re.compile(r"(?:(?P<uq>[A-Za-z0-9_\-\.\[\]]+)|(?:\'(?P<q>(?:\'\'|[^\'])+)\'))") # Glottolog NEXUS was buggy: 'Dzu'oasi [dzuo1238]'
    # taxa_re = re.compile(r"(?:(?P<uq>[A-Za-z0-9_\-\.\[\]]+)|(?:\'(?P<q>[^\,\(\)\:]+)\'))") # HACK
    branch_re = re.compile(r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)")

    def __init__(self, dat):
        self.dat = dat

    def parse(self):
        tokens = self._tokenize(self.dat)
        return self._parse(tokens)

    def _tokenize(self, data):
        tree_data = data[data.find('('):-1]  # skip tree-level annotation & strip the last semicolon
        idx = 0
        tokens = []
        while idx < len(tree_data):
            if tree_data[idx] in ("\n",):
                idx += 1
            elif tree_data[idx] == '(':
                tokens.append(self.START)
                idx += 1
            elif tree_data[idx] == ')':
                tokens.append(self.END)
                idx += 1
            elif tree_data[idx] == ',':
                tokens.append(self.COMMA)
                idx += 1
            elif tree_data[idx] == ';':
                tokens.append(self.EOS)
                idx += 1
            elif tree_data[idx] == '[':
                # annotation
                idx2 = tree_data.find(']', idx + 1)
                rawstr = tree_data[idx + 1:idx2]
                annotation = {}
                for kv in rawstr.split(','):
                    k, v = kv.split("=", 1)
                    annotation[k] = v
                obj = {
                    'type': self.ANNOTATION,
                    'annotation': annotation,
                }
                idx = idx2 + 1
                tokens.append(obj)
            elif tree_data[idx] == ':':
                match = self.branch_re.search(tree_data, idx + 1)
                assert(match is not None)
                obj = {
                    'type': self.BRANCH,
                    'branch': float(tree_data[match.start():match.end()]),
                }
                idx = match.end()
                tokens.append(obj)
            else:
                match = self.taxa_re.search(tree_data, idx)
                assert(match is not None)
                if match.group('uq'):
                    taxa = match.group('uq')
                else:
                    taxa = match.group('q')
                    taxa = re.sub(r"\'\'", "\'", taxa) # unescape
                obj = {
                    'type': self.NODE,
                    'taxa': taxa,
                }
                idx = match.end()
                tokens.append(obj)
        return tokens

    def _parse(self, tokens):
        count = 0
        root = Node(_id=count)
        count += 1
        node = root
        rv = []
        for token in tokens:
            if token == self.START:
                node2 = Node(_id=count)
                count += 1
                assert(len(node.children) == 0)
                node.children.append(node2)
                node2.parent = node
                node = node2
            elif token == self.END:
                if hasattr(node, "parent"):
                    node = node.parent
                else:
                    node = None
            elif token == self.EOS:
                rv.append(root)
                root = Node(_id=count)
                count += 1
                node = root
            elif token == self.COMMA:
                node2 = Node(_id=count)
                count += 1
                assert(len(node.parent.children) > 0)
                node.parent.children.append(node2)
                node2.parent = node.parent
                node = node2
            elif token['type'] == self.ANNOTATION:
                node.annotation = token['annotation']
            elif token['type'] == self.BRANCH:
                node.branch = token['branch']
            elif token['type'] == self.NODE:
                node.name = token['taxa']
        return rv

if __name__ == "__main__":
    f = open(sys.argv[1], "r")
    tp = TreeParser(f.read())
    trees = tp.parse()
    with open(sys.argv[2], "wb") as f:
        pickle.dump(trees, f)
