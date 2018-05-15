# -*- coding: utf-8 -*-
#
# convert CSV into two JSON files (language list and feature list)
#
import sys, os
import csv
import json
from argparse import ArgumentParser
from collections import defaultdict

def main():
    parser = ArgumentParser()
    parser.add_argument("--fthres", metavar="IDX", type=int, default=150)
    parser.add_argument("csv", metavar="WALS_LANG", default=None)
    parser.add_argument("langs", metavar="MERGED", default=None)
    parser.add_argument("flist", metavar="FLIST", default=None)
    args = parser.parse_args()

    # name <- wals_code (uniqueness assumed)
    legend = ("name", "iso_code", "glottocode", "fullname", "latitude", "longitude", "genus", "family", "macroarea", "countrycodes")
    with open(args.csv, newline="", encoding='utf-8') as f:
        reader = csv.reader(f)
        dat = list(reader)

    # process legend
    row = dat.pop(0)
    idx2fid = {}
    flist_raw = []
    for k in legend:
        row.pop(0)
    for i, k in enumerate(row):
        name = k[0 : k.index(" ")]
        # idx2fid[i] = len(flist_raw)
        # idx: WALS index
        # id: our own index (will be changed by filtering)
        flist_raw.append({ "size": -1, "type": "cat", "annotation": { "orig_idx": i, "name": name, "fullname": k, "label2vid": {}, "vid2label": [] }})
    sys.stderr.write("{} features (original)\n".format(len(flist_raw)))


    # sort first by suffix and then by main id
    flist_sorted = sorted(flist_raw, key=lambda x: x["annotation"]["name"][-1])
    flist_sorted = sorted(flist_sorted, key=lambda x: int(x["annotation"]["name"][:-1]))
    bad_features = {
        "10B Nasal Vowels in West Africa": "too specific",
        "39B Inclusive/Exclusive Forms in Pama-Nyungan": "too specific",
        "139A Irregular Negatives in Sign Languages": "too specific",
        "140A Question Particles in Sign Languages": "too specific",
        "141A Writing Systems": "too specific",

        "81B Languages with two Dominant Orders of Subject, Object, and Verb": "subdividing a low-coverage feature value",
        "90D Internally-headed relative clauses": "subdividing a low-coverage feature value",
        "90E Correlative relative clauses": "subdividing a low-coverage feature value",
        "90F Adjoined relative clauses": "subdividing a low-coverage feature value",
        "90G Double-headed relative clauses": "subdividing a low-coverage feature value",
        "130B Cultural Categories of Languages with Identity of 'Finger' and 'Hand'": "subdividing a low-coverage feature value",
        "143G Minor morphological means of signaling negation": "subdividing a low-coverage feature value", # None dominates
        "144C Languages with different word order in negative clauses": "subdividing a low-coverage feature value",
        "144E Multiple Negative Constructions in SVO Languages": "subdividing a low-coverage feature value",
        "144F Obligatory Double Negation in SVO languages": "subdividing a low-coverage feature value",
        "144G Optional Double Negation in SVO languages": "subdividing a low-coverage feature value",
        "144M Multiple Negative Constructions in SOV Languages": "subdividing a low-coverage feature value",
        "144N Obligatory Double Negation in SOV languages": "subdividing a low-coverage feature value",
        "144O Optional Double Negation in SOV languages": "subdividing a low-coverage feature value",
        "144T The Position of Negative Morphemes in Verb-Initial Languages": "subdividing a low-coverage feature value",
        "144U Double negation in verb-initial languages": "subdividing a low-coverage feature value",
        "144V Verb-Initial with Preverbal Negative": "subdividing a low-coverage feature value",
        "144W Verb-Initial with Negative that is Immediately Postverbal or between Subject and Object": "subdividing a low-coverage feature value",
        "144X Verb-Initial with Clause-Final Negative": "subdividing a low-coverage feature value",
        "144Y The Position of Negative Morphemes in Object-Initial Languages": "subdividing a low-coverage feature value",
    }
    flist_refined = []
    idx2fid = {}
    for fstruct in flist_sorted:
        if fstruct["annotation"]["fullname"] in bad_features:
            sys.stderr.write("drop {}\t{}\n".format(fstruct["annotation"]["fullname"], bad_features[fstruct["annotation"]["fullname"]]))
            fstruct["fid"] = -1
        else:
            fstruct["fid"] = len(flist_refined)
            idx2fid[fstruct["annotation"]["orig_idx"]] = fstruct["fid"]
            flist_refined.append(fstruct)
    sys.stderr.write("{} features after refinement\n".format(len(flist_refined)))    


    # bernoulli features [FALSE, TRUE]
    bin_features = {
        "25B Zero Marking of A and P Arguments": ["2 Non-zero marking", "1 Zero-marking"],
        "47A Intensifiers and Reflexive Pronouns": ["2 Differentiated", "1 Identical"],
        "58A Obligatory Possessive Inflection": ["2 Absent", "1 Exists"],
        "63A Noun Phrase Conjunction": ["1 'And' different from 'with'", "2 'And' identical to 'with'"],
        "65A Perfective/Imperfective Aspect": ["2 No grammatical marking", "1 Grammatical marking"],
        "67A The Future Tense": ["2 No inflectional future", "1 Inflectional future exists"],
        "73A The Optative": ["2 Inflectional optative absent", "1 Inflectional optative present"],
        "107A Passive Constructions": ["2 Absent", "1 Present"],
        "119A Nominal and Locational Predication": ["1 Different", "2 Identical"],
        "120A Zero Copula for Predicate Nominals": ["1 Impossible", "2 Possible"],
        "129A Hand and Arm": ["2 Different", "1 Identical"],
        "130A Finger and Hand": ["2 Different", "1 Identical"],
        "136B M in First Person Singular": ["1 No m in first person singular", "2 m in first person singular"],
        "137B M in Second Person Singular": ["1 No m in second person singular", "2 m in second person singular"],
    }
    for fstruct in flist_refined:
        if fstruct["annotation"]["fullname"] in bin_features:
            sys.stderr.write("bin feature {}\n".format(fstruct["annotation"]["fullname"]))
            fstruct["size"] = 1
            fstruct["type"] = "bin"
            fstruct["annotation"]["vid2label"] = list(bin_features[fstruct["annotation"]["fullname"]])
            fstruct["annotation"]["label2vid"] = { fstruct["annotation"]["vid2label"][0]: 0, fstruct["annotation"]["vid2label"][1]: 1 }
    flist = flist_refined


    # process each language
    sys.stderr.write("{} languages (original)\n".format(len(dat)))
    bad_languages = {
        "Jurchen": "Hezhen (Nanai) mixed up with historical Jurchen",
    }
    langs = []
    while len(dat) > 0:
        row = dat.pop(0)
        lang = { "annotation": { "source": "WALS", "features": {} } }
        for k in legend:
            v = row.pop(0)
            if k == "glottocode":
                lang[k] = v
            elif k in ("latitude", "longitude"):
                lang[k] = float(v)
            else:
                lang["annotation"][k] = v
        lang["phylo"] = lang["annotation"]["genus"]

        if lang["annotation"]["family"] == "other":
            # "Sign Languages" or "Creoles and Pidgins"
            sys.stderr.write("drop {}\n".format(lang["annotation"]["fullname"]))
            continue
        if lang["annotation"]["name"] in bad_languages:
            sys.stderr.write("drop bad language {}\n".format(lang["annotation"]["fullname"]))
            continue

        for idx, fid in idx2fid.items():
            v = row[idx]
            if len(v) <= 0:
                # missing value
                continue
            fstruct = flist[fid]
            name = fstruct["annotation"]["name"]
            if fstruct["type"] == "bin":
                if v not in fstruct["annotation"]["label2vid"]:
                    sys.stderr.write("corrupt bin feature?\t{}\t{}".format(fstruct["annotation"]["fullname"], v))
                    exit(1)
                lang["annotation"]["features"][name] = fstruct["annotation"]["label2vid"][v]
            elif fstruct["type"] == "cat":
                vid = int(v[0 : v.index(" ")]) - 1
                lang["annotation"]["features"][name] = vid
                fstruct["annotation"]["label2vid"][v] = vid
            else:
                sys.stderr.write("unsupported type\t{}\n".format(fstruct["type"]))
                exit(1)
        langs.append(lang)
    sys.stderr.write("{} languages (refined)\n".format(len(langs)))


    # fill vid2label
    for fstruct in flist:
        if fstruct["type"] == "cat":
            fstruct["size"] = max(fstruct["annotation"]["label2vid"].values()) + 1
            fstruct["annotation"]["vid2label"] = [None] * fstruct["size"]
            for k, vid in fstruct["annotation"]["label2vid"].items():
                fstruct["annotation"]["vid2label"][vid] = k


    # fill logically determined values
    logical_specs = [
        {
            "name": "81A",  # Order of Subject, Object and Verb
            "affects": [
                { "flist": ["144D", "144H", "144I", "144J", "144K"],
                  "svals": ["2 SVO", "7 No dominant order"] },
                # only appliy to SVO languages
                # - 144D The Position of Negative Morphemes in SVO Languages
                # - 144H NegSVO Order
                # - 144I SNegVO Order
                # - 144J SVNegO Order
                # - 144K SVONeg Order
                { "flist": ["144L", "144P", "144Q", "144R", "144S"],
                  "svals": ["1 SOV", "7 No dominant order"] },
                # only apply to SOV languages
                # - 144L The Position of Negative Morphemes in SOV Languages
                # - 144P NegSOV Order
                # - 144Q SNegOV Order
                # - 144R SONegV Order
                # - 144S SOVNeg Order
            ]
        },
        {
            "name": "90A", # Order of Relative Clause and Noun
            "affects": [
                { "flist": ["90B"],
                  "svals": ["2 Relative clause-Noun", "7 Mixed"] },
                # - 90B Prenominal relative clauses
                { "flist": ["90C"],
                  "svals": ["1 Noun-Relative clause", "7 Mixed"] },
                # - 90C Postnominal relative clauses
            ]
        }
    ]
    def get_struct_by_name(flist, name):
        for i, fstruct in enumerate(flist):
            if fstruct["annotation"]["name"] == name:
                return fstruct
        return None

    for spec in logical_specs:
        for affected in spec["affects"]:
            affected["uvid"] = []
            for name in affected["flist"]:
                fstruct = get_struct_by_name(flist, name)
                fstruct["annotation"]["label2vid"]["-1 Undefined"] = len(fstruct["annotation"]["vid2label"])
                fstruct["annotation"]["vid2label"].append("-1 Undefined")
                affected["uvid"].append(fstruct["annotation"]["label2vid"]["-1 Undefined"])
                fstruct["size"] += 1

    filled = 0
    for spec in logical_specs:
        src_name = spec["name"]
        src_fstruct = get_struct_by_name(flist, src_name)

        for lang in langs:
            if src_name not in lang["annotation"]["features"]:
                continue
            src_val = src_fstruct["annotation"]["vid2label"][lang["annotation"]["features"][src_name]]
            for affected in spec["affects"]:
                if src_val in affected["svals"]:
                    continue
                for i, name in enumerate(affected["flist"]):
                    if name in lang["annotation"]["features"]:
                        sys.stderr.write("illogical feature values: {}\t{}: {}\t{}: {}\n".format(lang["annotation"]["fullname"], src_name, src_val, name, lang["annotation"]["features"][name]))
                    else:
                        lang["annotation"]["features"][name] = affected["uvid"][i]
                        filled += 1
    sys.stderr.write("{} elements filled\n".format(filled))


    # filter out low-coverage features
    if args.fthres > 0:
        fcounts = defaultdict(int)
        for lang in langs:
            for k, v in lang["annotation"]["features"].items():
                fcounts[k] += 1
        flist2 = []
        name_list = {}
        for fstruct in flist:
            count = fcounts[fstruct["annotation"]["name"]]
            if count >= args.fthres:
                fstruct["fid"] = len(flist2)
                flist2.append(fstruct)
                name_list[fstruct["annotation"]["name"]] = True
            else:
                sys.stderr.write("remove low-coverage feature\t{}\t({} counts)\n".format(fstruct["annotation"]["fullname"], count))
        flist = flist2
        deleted = 0
        for lang in langs:
            for name in list(lang["annotation"]["features"].keys()):
                if name not in name_list:
                    del lang["annotation"]["features"][name]
                    deleted += 1
        sys.stderr.write("{} elements deleted\n".format(deleted))
        sys.stderr.write("{} features remain (binsize {})\n".format(len(flist), sum(map(lambda x: x["size"], flist))))
        filled = sum(fcounts.values()) - deleted
        total = len(langs) * len(flist)
        sys.stderr.write("coverage\t{} ({} / {})\n".format(filled / total, filled, total))


    # create catvect
    name2fid = {}
    for fstruct in flist:
        name2fid[fstruct["annotation"]["name"]] = fstruct["fid"]
    for lang in langs:
        catvect = [-1] * len(flist)
        for name, v in lang["annotation"]["features"].items():
            catvect[name2fid[name]] = v
        lang["catvect"] = catvect
        

    with open(args.flist, 'w') as f:
        f.write("%s\n" % json.dumps(flist, indent=4, sort_keys=True))

    with open(args.langs, 'w') as f:
        for lang in langs:
            f.write("%s\n" % json.dumps(lang))

if __name__ == "__main__":
    main()
