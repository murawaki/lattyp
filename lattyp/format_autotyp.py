# -*- coding: utf-8 -*-
#
# convert CSV intoto two JSON files (language list and feature list)
#
import sys, os
import csv
import json
from argparse import ArgumentParser
from collections import defaultdict

def main():
    parser = ArgumentParser()
    parser.add_argument("--binminthres", metavar="IDX", type=int, default=5)
    parser.add_argument("--fthres", metavar="IDX", type=int, default=50)
    parser.add_argument("basedir", metavar="AUTOTYP_BASEDIR", default=None)
    parser.add_argument("langs", metavar="MERGED", default=None)
    parser.add_argument("flist", metavar="FLIST", default=None)
    args = parser.parse_args()

    if not os.path.isdir(args.basedir):
        sys.stderr.write("specificy autotyp basedir\n")
        exit(1)

    # languoid_path = "/home/murawaki/research/lattyp/data/glottolog/languoid.csv"
    # glottocode2lang = {}
    # legend = ("glottocode", "family_id", "parent_glottocode", "name", "bookkeeping", "level", "status", "latitude", "longitude", "iso639", "description", "markup_description", "child_family_count", "child_language_count", "child_dialect_count", "country_ids")
    # with open(languoid_path, newline="", encoding='utf-8') as f:
    #     reader = csv.reader(f)
    #     next(reader) # skip legend
    #     for row in reader:
    #         lang = {}
    #         for k, v in zip(legend, row):
    #             if len(v) <= 0:
    #                 continue
    #             if k == "glottocode":
    #                 if v in glottocode2lang:
    #                     sys.stderr.write("non-unique glottocode\t{}\n".format(v))
    #                 glottocode2lang[v] = lang
    #             else:
    #                 lang[k] = v    

    # load languages
    # ancient written languages to be removed from the list
    bad_langs = {
        "egyp1246", # Egyptian (Ancient)
        "akka1240", # Akkadian
        "sume1241", # Sumerian
        "elam1244", # Elamite
        "hurr1240", # Hurrian
        "hitt1242", # Hittite
        "anci1244", # Hebrew (Biblical)
        "etru1241", # Etruscan
        "copt1239", # Coptic
        "clas1254", # Tibetan (Classical)
        "anci1242", # Greek (Ancient)
        "goth1244", # Gothic
        "ugar1238", # Ugaritic
        "jewi1240", # Aramaic (Biblical)
        "oldh1241", # High German (Old)
        "olds1250", # Old Saxon
        "oldn1244", # Old Norse
        "oldf1241", # Frisian (Old)
        "oldr1238", # Russian (Old)
        "bolg1250", # Bulgarian (Old)
        "pali1273", # Pali
        "lati1261", # Latin
        "aves1237", # Avestan
        "sama1313", # Samaritan
        "clas1252", # Syriac
        "oldi1245", # Irish (Old)
    }
    langs = []
    # name <- LID (uniqueness assumed)
    legend = ("name", "fullname", "stock", "lowestsubbranch", "majorbranch", "quasistock", "subbranch", "aliases", "genesis", "longitude", "latitude","nearprotohomeland", "modality", "localregion", "subsubbranch", "stock_alias", "stock_aliases_all", "aliases_all", "subsistence", "iso639", "glottocode", "origin_continent", "continent", "area", "macrocontinent")
    with open(os.path.join(args.basedir, "data", "Register.csv"), newline="", encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader) # skip legend
        dat = list(reader)
    sys.stderr.write("{} languages (original)\n".format(len(dat)))
    phylos = defaultdict(int)
    for row in dat:
        lang = { "annotation": { "source": "autotyp", "features": {} } }
        for k, v in zip(legend, row):
            if len(v) <= 0:
                continue
            if k == "glottocode":
                lang[k] = v
            if k in ("longitude", "latitude"):
                lang[k] = float(v)
            else:
                lang["annotation"][k] = v
        if lang["annotation"]["modality"] == "signed":
            sys.stderr.write("drop sign language\t{}\n".format(lang["annotation"]["name"]))
            continue
        if lang["annotation"]["genesis"] in ("creole", "mixed"):
            sys.stderr.write("drop creole mixed language\t{}\n".format(lang["annotation"]["name"]))
            continue
        if "glottocode" in lang and lang["glottocode"] in bad_langs:
            sys.stderr.write("drop ancient written language\t{} ({})\n".format(lang["annotation"]["name"], lang["glottocode"]))
            continue
        if "longitude" not in lang:
            sys.stderr.write("drop language without geo. pointer\t{}\n".format(lang["annotation"]["name"]))
        # if "glottocode" in lang and lang["glottocode"] in glottocode2lang:
        #     glang = glottocode2lang[lang["glottocode"]]
        #     if "status" in glang and glang["status"] == "extinct":
        #         sys.stderr.write("\textinct according to Glottolog\t{}\t{}\n".format(lang["annotation"]["name"], lang["glottocode"])

        # assign phylo
        #
        # NOTE: problematic 'Western Malayo-Polynesian non-clade'
        phylo = None
        for level in ("majorbranch", "stock", "subbranch", "subsubbranch", "lowestsubbranch", "quasistock"):
            if level not in lang["annotation"]:
                continue
            if level == "majorbranch" and lang["annotation"][level] in ("Malayo-Polynesian", "Indo-Iranian", "Italic-Celtic", "Finno-Ugric", "Balto-Slavic", "Daghestanian", "Eastern Mon-Khmer", "Southern Mon-Khmer"):
                level = "subbranch"
                if level not in lang["annotation"]:
                    phylo = lang["annotation"]["majorbranch"]
                    break
            phylo = lang["annotation"][level]
            break
        if phylo is None:
            sys.stderr.write("no phylogenetic information\t{}\n".format(lang["annotation"]["name"]))
            continue
        lang["phylo"] = phylo
        phylos[phylo] += 1
        # print("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(phylo, lang["annotation"].get("quasistock", "_"), lang["annotation"].get("stock", "_"), lang["annotation"].get("majorbranch", "_"), lang["annotation"].get("subbranch", "_"), lang["annotation"].get("subsubbranch", "_"), lang["annotation"].get("lowestsubbranch", "_")))
        langs.append(lang)
    sys.stderr.write("{} languages (refined)\n".format(len(langs)))
    sys.stderr.write("{} phylogenetic groups\n".format(len(phylos)))

    name2lang = {}
    for lang in langs:
        name2lang[lang["annotation"]["name"]] = lang


    specs = [
        { "path": "Agreement.csv",
          "features": [
              ("VPolyagreement.Presence.v2", "bin", None),
              ("VPolyagreement.Presence.v1", "bin", None),
          ],
        },
        { "path": "Alienability.csv",
          "features": [
              ("NPLocusDistribution", "cat", None),
              ("NPPossClassificationType", "modbin", { "classificatory": False, "semantic": True }),
              ("NPPossBoundNouns.n", "count", None),
              ("NPPossBoundNouns.Presence", None, None),     # duplicate { "none": False, "some": True }),
              ("NPPossClassesMultiple.Presence", "bin", None),
              ("NPAlienableClasses.n", "count", None),
              ("NPAlienableClasses.Presence", None, None),   # duplicate
              ("NPInalienableClasses.n", "count", None),
              ("NPInalienableClasses.Presence", None, None), # duplicate
              ("NPPossClasses.n", "count", None),
              ("NPPossClasses.Presence", None, None),        # duplicate (some incompatibilities)
              ("NPPossClassificationAny.Presence", "bin", None),
          ]
        },
        { "path": "Alignment_per_language.csv",
          "features": [
              ("AlignmentAgrAccDegree.binned3", "cat", None),
              ("AlignmentAgrAccDegree", None, None),   # ratio
              ("AlignmentAgrDominantPTG", "cat", None),
              ("AlignmentAgrDominantSAP", "cat", None),
              ("AlignmentAgrDominantSAPTG", "cat", None),
              ("AlignmentAgrErgDegree.binned3", "cat", None),
              ("AlignmentAgrErgDegree", None, None),   # ratio, duplicate
              ("AlignmentCaseAccDegree.binned3", "cat", None),
              ("AlignmentCaseAccDegree", None, None), # ratio
              ("AlignmentCaseDominantPTG", "cat", None),
              ("AlignmentCaseDominantSAP", "cat", None),
              ("AlignmentCaseDominantSAPTG", "cat", None),
              ("AlignmentCaseErgDegree.binned3", "cat", None),
              ("AlignmentCaseErgDegree", None, None), # ratio
              ("AlignmentCaseDominantLexSAP", "cat", None),
              ("AlignmentCaseDominantProSAP", "cat", None),
              ("CaseSplitByPolarity.Presence", "bin", None),
              ("AgrSplitByRef.Presence", "bin", None),
              ("CaseSplitByRef.Presence", "bin", None),
              ("AgrSplitByValClass.Presence", "bin", None),
              ("CaseSplitByValClass.Presence", "bin", None),
              ("AgrSplitByRank.Presence", "bin", None),
              ("CaseSplitByRank.Presence", "bin", None),
              ("AgrSplitByCat.Presence", "bin", None),
              ("CaseSplitByCat.Presence", "bin", None),
              ("AlignmentAgrMarkerAccDegree", None, None), # ratio
              ("AlignmentAgrMarkerDominantSAP", "cat", None),
              ("AlignmentAgrMarkerErgDegree.binned3", "cat", None),
              ("AlignmentAgrMarkerErgDegree", None, None), # ratio
          ]
        },
        { "path": "Clause_word_order.csv",
          "features": [
              ("WordOrderAPLex", "cat", None),
              ("WordOrderAPBasicLex", "modbin", { "OA": False, "AO": True }),
              ("WordOrderAVBasicLex", "modbin", { "head-dep": False, "dep-head": True }),
              ("PositionVBasicLex", "cat", None),
              ("WordOrderAPVBasicLex", "cat", None),
              ("PositionVBasicLex.binned.v4", "modbin", { "non_final": False, "final": True }),
              ("PositionVBasicLex.binned.v1", "modbin", { "other": False, "final_or_free": True }),
              ("PositionVBasicLex.binned.v3", "modbin", { "other": False, "initial": True }),
              ("PositionVBasicLex.binned.v2", "modbin", { "other": False, "medial": True }),
              ("WordOrderPVBasicLex", "modbin", { "head-dep": False, "dep-head": True }),
              ("WordOrderAPSecondaryLex", "modbin", { "OA": False, "AO": True }),
              ("WordOrderAVSecondaryLex", "modbin", { "head-dep": False, "dep-head": True }),
              ("PositionVSecondaryLex", "cat", None),
              ("WordOrderAPVSecondaryLex", "cat", None),
              ("PositionVSecondaryLex.binned.v4", "modbin", { "non_final": False, "final": True }),
              ("PositionVSecondaryLex.binned.v1", "modbin", { "other": False, "final_or_free": True }),
              ("PositionVSecondaryLex.binned.v3", "modbin", { "other": False, "initial": True }),
              ("PositionVSecondaryLex.binned.v2", "modbin", { "other": False, "medial": True }),
              ("WordOrderPVSecondaryLex", "modbin", { "head-dep": False, "dep-head": True }),
              ("VFinalLex", "cat", None),
              ("VFinalOrFreeLex", "cat", None),
              ("WordOrderAPVLexFlexible", "cat", None),
              ("VInitialLex", "cat", None),
              ("VMedialLex", "cat", None),
              ("WordOrderSVBasicLex", "cat", None),
              ("WordOrderSVSecondaryLex", None, None), # only one datapoint
          ]
        },
        { "path": "Clusivity.csv",
          "features": [
              ("InclExclAsPerson.Presence", "bin", None),
              ("InclExclAny.Presence", None, None), # duplicate
              ("InclExclType", "cat", None),
              ("InclExclAsMinAug.Presence", "bin", None),
          ]
        },
        { "path": "Gender.csv",
          "features": [
              ("Gender.n", "count", None),
              ("Gender.binned4", None, None),  # duplicate
              ("Gender.Presence", None, None), # duplicate; FALSE -> Gender.n = 0
          ],
        },
        { "path": "GR_per_language.csv",
          "features": [
              ("CoArgSensitivityAgr.Presence", "bin", None),
              ("VAgreement.Presence.v1", "bin", None),
              ("CoArgSensitivityCase.Presence", "bin", None),
          ]
        },
        { "path": "Locus_per_language.csv",
          "features": [
              ("LocusA.binned6", None, None),  # near-duplicate
              ("Locus.A.default.binned6", "cat", None),
              ("Locus.A.default.if.N.binned6", None, None),        # only one datapoint
              ("Locus.A.default.if.N_or_Pro.binned6", None, None), # only one datapoint
              ("Locus.A.default.elsewhere.binned6", None, None),
              ("Locus.A.default.if.Pro.binned6", None, None),
              ("Locus.A.default.if.SAP.binned6", None, None),      # only one datapoint
              ("Locus.A.exp.binned6", None, None),
              ("Locus.A.poss.binned6", None, None),
              ("Locus.Act.binned6", None, None),
              ("Locus.ATTR.default.binned6", "cat", None),
              ("Locus.ATTR.default.elsewhere.binned6", None, None),
              ("Locus.B.ben.binned6", None, None),                 # only one datapoint
              ("Locus.B.default.binned6", "cat", None),
              ("Locus.B.default.elsewhere.binned6", None, None),
              ("Locus.B.default.if.Pro.binned6", None, None),
              ("Locus.B.rec.binned6", None, None),
              ("Locus.G.default.binned6", None, None),
              ("Locus.I.default.binned6", None, None),
              ("LocusP.binned6", "cat", None),
              ("Locus.Pat.binned6", None, None),
              ("LocusPOSS.binned6", None, None), # near-duplicate
              ("Locus.POSS.default.binned6", "cat", None),
              ("Locus.POSS.default.if.1PERS.binned6", None, None),
              ("Locus.POSS.default.if.N.binned6", None, None),
              ("Locus.POSS.default.if.N.alien.binned6", None, None),
              ("Locus.POSS.default.if.N.inal.binned6", None, None),
              ("Locus.POSS.default.if.N.kin.binned6", None, None),
              ("Locus.POSS.default.if.N_or_Pro.binned6", None, None),
              ("Locus.POSS.default.elsewhere.binned6", None, None),
              ("Locus.POSS.default.if.Pro.binned6", None, None),
              ("LocusS.binned6", None, None),  # near-duplicate
              ("Locus.S.default.binned6", "cat", None),
              ("Locus.S.default.if.N_or_Pro.binned6", None, None),
              ("Locus.S.default.elsewhere.binned6", None, None),
              ("Locus.S.default.if.Pro.binned6", None, None),
              ("Locus.S.default.if.SAP.binned6", None, None),
              ("Locus.T.def.binned6", None, None),
              ("Locus.T.default.binned6", None, None),
              ("Locus.T.default.elsewhere.binned6", None, None),
              ("Locus.T.emp.binned6", None, None),
              ("Locus.U.def.binned6", None, None),
              ("Locus.U.default.binned6", "cat", None),
              ("Locus.U.default.if.1PERS.binned6", None, None),
              ("Locus.U.default.if.N.binned6", None, None),
              ("Locus.U.default.elsewhere.binned6", None, None),
              ("Locus.U.default.if.Pro.binned6", None, None),
              ("Locus.U.default.if.SAP.binned6", None, None),
              ("Locus.U.emp.binned6", None, None),
              ("Locus.U.ref.binned6", None, None),
              ("Locus.U.stim.binned6", None, None),
              ("LocusAandP.dm.Presence", "bin", None),
              ("LocusAandP.hm.Presence", "bin", None),
              ("LocusA.dm.Presence", "bin", None),
              ("LocusA.hm.Presence", "bin", None),
              ("LocusAorP.dm.Presence", "bin", None),
              ("LocusAorP.hm.Presence", "bin", None),
              ("LocusARG.dm.Presence", "bin", None),
              ("LocusP.dm.Presence", "bin", None),
              ("LocusP.hm.Presence", "bin", None),
              ("LocusPOSS.dm.Presence", "bin", None),
              ("LocusPOSS.hm.Presence", "bin", None),
              ("LocusSandA.dm.Presence", "bin", None),
              ("LocusSandA.hm.Presence", "bin", None),
              ("LocusS.dm.Presence", "bin", None),
              ("LocusS.hm.Presence", "bin", None),
              ("LocusA", None, None), # near-duplicate
              ("Locus.A.default", "cat", None),
              ("Locus.A.default.if.N", None, None),
              ("Locus.A.default.if.N_or_Pro", None, None),
              ("Locus.A.default.elsewhere", None, None),
              ("Locus.A.default.if.Pro", None, None),
              ("Locus.A.default.if.SAP", None, None),
              ("Locus.A.exp", None, None),
              ("Locus.A.poss", None, None),
              ("Locus.Act", None, None),
              ("LocusARG.hm.Presence", "bin", None),
              ("Locus.ATTR.default", "cat", None),
              ("Locus.ATTR.default.elsewhere", None, None),
              ("Locus.B.ben", None, None),
              ("Locus.B.default", "cat", None),
              ("Locus.B.default.elsewhere", None, None),
              ("Locus.B.default.if.Pro", None, None),
              ("Locus.B.rec", None, None),
              ("Locus.G.default", None, None),
              ("Locus.I.default", None, None),
              ("LocusP", "cat", None),
              ("Locus.Pat", None, None),
              ("LocusPOSS", None, None), # near-duplicate
              ("Locus.POSS.default", "cat", None),
              ("Locus.POSS.default.if.1PERS", None, None),
              ("Locus.POSS.default.if.N.alien", None, None),
              ("Locus.POSS.default.if.N.inal", None, None),
              ("Locus.POSS.default.if.N.kin", None, None),
              ("Locus.POSS.default.if.N", None, None),
              ("Locus.POSS.default.if.N_or_Pro", None, None),
              ("Locus.POSS.default.elsewhere", None, None),
              ("Locus.POSS.default.if.Pro", None, None),
              ("LocusS", None, None), # near-duplicate
              ("Locus.S.default", "cat", None),
              ("Locus.S.default.if.N_or_Pro", None, None),
              ("Locus.S.default.elsewhere", None, None),
              ("Locus.S.default.if.Pro", None, None),
              ("Locus.S.default.if.SAP", None, None),
              ("Locus.T.def", None, None),
              ("Locus.T.default", None, None),
              ("Locus.T.default.elsewhere", None, None),
              ("Locus.T.emp", None, None),
              ("Locus.U.def", None, None),
              ("Locus.U.default", "cat", None),
              ("Locus.U.default.if.1PERS", None, None),
              ("Locus.U.default.if.N", None, None),
              ("Locus.U.default.elsewhere", None, None),
              ("Locus.U.default.if.Pro", None, None),
              ("Locus.U.default.if.SAP", None, None),
              ("Locus.U.emp", None, None),
              ("Locus.U.ref", None, None),
              ("Locus.U.stim", None, None),
              ("LocusA.binned5", None, None), # near-duplicate
              ("Locus.A.default.binned5", "cat", None),
              ("Locus.A.default.if.N.binned5", None, None),
              ("Locus.A.default.if.N_or_Pro.binned5", None, None),
              ("Locus.A.default.elsewhere.binned5", None, None),
              ("Locus.A.default.if.Pro.binned5", None, None),
              ("Locus.A.default.if.SAP.binned5", None, None),
              ("Locus.A.exp.binned5", None, None),
              ("Locus.A.poss.binned5", None, None),
              ("Locus.Act.binned5", None, None),
              ("Locus.ATTR.default.binned5", "cat", None), # small?
              ("Locus.ATTR.default.elsewhere.binned5", None, None),
              ("Locus.B.ben.binned5", None, None),
              ("Locus.B.default.binned5", "cat", None), # small?
              ("Locus.B.default.elsewhere.binned5", None, None),
              ("Locus.B.default.if.Pro.binned5", None, None),
              ("Locus.B.rec.binned5", None, None),
              ("Locus.G.default.binned5", None, None),
              ("Locus.I.default.binned5", None, None),
              ("LocusP.binned5", "cat", None),
              ("Locus.Pat.binned5", None, None),
              ("LocusPOSS.binned5", None, None), # near-duplicate
              ("Locus.POSS.default.binned5", "cat", None),
              ("Locus.POSS.default.if.1PERS.binned5", None, None),
              ("Locus.POSS.default.if.N.binned5", None, None),
              ("Locus.POSS.default.if.N.alien.binned5", None, None),
              ("Locus.POSS.default.if.N.inal.binned5", None, None),
              ("Locus.POSS.default.if.N.kin.binned5", None, None),
              ("Locus.POSS.default.if.N_or_Pro.binned5", None, None),
              ("Locus.POSS.default.elsewhere.binned5", None, None),
              ("Locus.POSS.default.if.Pro.binned5", None, None),
              ("LocusS.binned5", None, None), # near-duplicate
              ("Locus.S.default.binned5", "cat", None),
              ("Locus.S.default.if.N_or_Pro.binned5", None, None),
              ("Locus.S.default.elsewhere.binned5", None, None),
              ("Locus.S.default.if.Pro.binned5", None, None),
              ("Locus.S.default.if.SAP.binned5", None, None),
              ("Locus.T.def.binned5", None, None),
              ("Locus.T.default.binned5", None, None),
              ("Locus.T.default.elsewhere.binned5", None, None),
              ("Locus.T.emp.binned5", None, None),
              ("Locus.U.def.binned5", None, None),
              ("Locus.U.default.binned5", "cat", None),
              ("Locus.U.default.if.1PERS.binned5", None, None),
              ("Locus.U.default.if.N.binned5", None, None),
              ("Locus.U.default.elsewhere.binned5", None, None),
              ("Locus.U.default.if.Pro.binned5", None, None),
              ("Locus.U.default.if.SAP.binned5", None, None),
              ("Locus.U.emp.binned5", None, None),
              ("Locus.U.ref.binned5", None, None),
              ("Locus.U.stim.binned5", None, None),
          ],
        },
        { "path": "Markers_per_language.csv",
          "features": [
              ("BehaviorCase.binned4", "cat", None),
              ("BehaviorNeg.binned4", "cat", None),
              ("BehaviorNounPlural.binned4", "cat", None),
              ("BehaviorTense.binned4", "modbin", { "final": False, "initial": True }),
              ("BehaviorCase.spreading.Presence", "bin", None),
              ("BehaviorNeg.spreading.Presence", None, None), # all FALSE
              ("BehaviorNounPlural.spreading.Presence", "bin", None),
              ("BehaviorTense.spreading.Presence", None, None), # all FALSE
              ("BehaviorCase", "cat", None),
              ("BehaviorNeg", "cat", None),
              ("BehaviorNounPlural", "cat", None),
              ("BehaviorTense", "cat", None),
              ("PolyexponentialityCase.Presence", "bin", None),
              ("PolyexponentialityNegation.Presence", "bin", None),
              ("PolyexponentialityNounPlural.Presence", "bin", None),
              ("PolyexponentialityTense.Presence", "bin", None),
              ("ExponentialityCase.n", "count", None),
              ("ExponentialityNegation.n", "count", None),
              ("ExponentialityNounPlural.n", "count", None),
              ("ExponentialityTense.n", "count", None),
              ("ExponentialityCase", "cat", None),
              ("ExponentialityNegation", "cat", None),
              ("ExponentialityNounPlural", "cat", None),
              ("ExponentialityTense", None, None), # 35 types for 167 datapoints
              ("FlexivityFormativeCase", "cat", None),
              ("FlexivityFormativeNeg", "cat", None),
              ("FlexivityFormativeNounPlural", "cat", None),
              ("FlexivityFormativeTense", "cat", None),
              ("FusionCase.binned.binned6", "cat", None),
              ("FusionNegation.binned.binned6", "cat", None),
              ("FusionNounPlural.binned.binned6", "cat", None),
              ("FusionTense.binned.binned6", "cat", None),
              ("FusionCase", "cat", None),
              ("FusionNegation", "cat", None),
              ("FusionNounPlural", "cat", None),
              ("FusionTense", "cat", None),
              ("FusionCase.isolating.Presence", "bin", None),
              ("FusionNegation.isolating.Presence", "bin", None),
              ("FusionNounPlural.isolating.Presence", "bin", None),
              ("FusionTense.isolating.Presence", "bin", None),
              ("FusionCase.nonlinear.Presence", "modbin", { "linear": False, "nonlinear": True }),
              ("FusionNegation.nonlinear.Presence", "modbin", { "linear": False, "nonlinear": True }),
              ("FusionNounPlural.nonlinear.Presence", "modbin", { "linear": False, "nonlinear": True }),
              ("FusionTense.nonlinear.Presence", "modbin", { "linear": False, "nonlinear": True }),
              ("FusionCase.reduplicative.Presence", None, None), # only one value 'non-reduplicative'
              ("FusionNegation.reduplicative.Presence", None, None), # only one value 'non-reduplicative'
              ("FusionNounPlural.reduplicative.Presence", "modbin", { "non-reduplicative": False, "reduplicative": True }),
              ("FusionTense.reduplicative.Presence", None, None), # only one exception
              ("FusionCase.tonal.Presence", "bin", None),
              ("FusionNegation.tonal.Presence", None, None), # all FALSE
              ("FusionNounPlural.tonal.Presence", None, None), # all FALSE
              ("FusionTense.tonal.Presence", "bin", None),
              ("MorphemeTypeCase", "cat", None),
              ("MorphemeTypeNeg", "modbin", { "formative": False, "PoSWd": True }),
              ("PositionCase.post.Presence", "bin", None),
              ("PositionNeg.post.Presence", "bin", None),
              ("PositionNounPlural.post.Presence", "bin", None),
              ("PositionTense.post.Presence", "bin", None),
              ("PositionCase.binned4", "cat", None),
              ("PositionNeg.binned4", "cat", None),
              ("PositionNounPlural.binned4", "cat", None),
              ("PositionTense.binned4", "cat", None),
              ("PositionCase.binned5", "cat", None),
              ("PositionNeg.binned5", "cat", None),
              ("PositionNounPlural.binned5", "cat", None),
              ("PositionTense.binned5", "cat", None),
              ("PositionCase", "cat", None),
              ("PositionNeg", "cat", None),
              ("PositionNounPlural", "cat", None),
              ("PositionTense", "cat", None),
              ("PositionCase.pre.Presence.v1", "bin", None),
              ("PositionNeg.pre.Presence.v1", "bin", None),
              ("PositionNounPlural.pre.Presence.v1", "bin", None),
              ("PositionTense.pre.Presence.v1", "bin", None),
              ("PositionCase.pre.Presence.v2", "bin", None),
              ("PositionNeg.pre.Presence.v2", "bin", None),
              ("PositionNounPlural.pre.Presence.v2", "bin", None),
              ("PositionTense.pre.Presence.v2", "bin", None),
              ("PositionCase.simul.Presence.v1", None, None), # all FALSE
              ("PositionNeg.simul.Presence.v1", "bin", None),
              ("PositionNounPlural.simul.Presence.v1", "bin", None),
              ("PositionTense.simul.Presence.v1", "bin", None),
              ("PositionCase.simul.Presence.v2", None, None), # all FALSE
              ("PositionNeg.simul.Presence.v2", "bin", None),
              ("PositionNounPlural.simul.Presence.v2", "bin", None),
              ("PositionTense.simul.Presence.v2", "bin", None),
              ("HostRestrictionsCase", "modbin", { "restricted": False, "semi-restricted": True }),
              ("HostRestrictionsNegation", "cat", None),
              ("HostRestrictionsNounPlural", "modbin", { "restricted": False, "semi-restricted": True }),
              ("HostRestrictionsTense", "modbin", { "restricted": False, "semi-restricted": True }),
              ("FlexivityStemCase", "cat", None),
              ("FlexivityStemNeg", "cat", None),
              ("FlexivityStemNounPlural", "cat", None),
              ("FlexivityStemTense", "cat", None),
          ],
        },
        { "path": "Morpheme_types.csv",
          "features": [
              ("Enclitic.Presence", "bin", None),
              ("Endoclitic.Presence", "bin", None),
              ("Infix.Presence", "bin", None),
              ("MorphemeType.n", "pcount", None),
              ("PreMorpheme.n", "count", None),
              ("PostMorpheme.n", "count", None),
              ("Prefix.Presence", "bin", None),
              ("Proclitic.Presence", "bin", None),
              ("Suffix.Presence", "bin", None),
          ],
        },
        { "path": "Morphology_per_language.csv",
          "features": [
              ("VAgreement.Presence.v2", "bin", None),
              ("Flexivity.Presence", "bin", None),
              ("FlexivityLexical.Presence", "bin", None),
              ("Polyexponence.Presence", "bin", None),
              ("FlexivityNoun.Presence", "bin", None),
              ("FlexivityVerb.Presence", "bin", None),
              ("FlexivityLexicalNoun.Presence", "bin", None),
              ("FlexivityLexicalVerb.Presence", "bin", None),
              ("PolyexponenceNoun.Presence", "bin", None),
              ("PolyexponenceVerb.Presence", "bin", None),
              ("SlotsWithZerosPost.n", "count", None),
              ("SlotsWithZerosPrae.n", "count", None),
          ],
        },
        { "path": "NP_per_language.csv",
          "features": [
              ("AdjAttrAgr.Presence", "bin", None),
              ("AdjAttrConstr.Presence", "bin", None),
              ("AdjAttrGvt.Presence", "bin", None),
              ("AdjAttrMarking.overt.Presence", "bin", None),
              ("NPAgr.Presence", "bin", None),
              ("NPConstr.Presence", "bin", None),
              ("NPGvt.Presence", "bin", None),
              ("NPMarking.overt.Presence", "bin", None),
          ],
        },
        { "path": "NP_word_order.csv",
          "features": [
              ("WordOrderNPBasic", "cat", None),
          ],
        },
        { "path": "Numeral_classifiers.csv",
          "features": [
              ("NumClass.n", "count", None),
              ("NumClass.Presence", "bin", None), # TRUE -> NumClass.n > 0
          ],
        },
        { "path": "Rhythm_per_language.csv",
          "features": [
              ("RhythmType", "cat", None),
          ],
        },
        { "path": "Synthesis.csv",
          "features": [
              ("VInlfCatSurveyComplete", None, None), # FALSE -> kill entries??
              ("VBipartiteStem.Presence", "bin", None),
              ("VInflCat", None, None),             # 264 types!
              ("VInflMacrocategories", None, None), # 97 types
              ("VInflCatAndAgrMax.n", "count", None),
              ("VInflCatMax.n", "count", None),
              ("VInflCatMax.binned3", None, None), # duplicate
              ("VInflCatAndAgrMax.binned3", None, None), # duplicate
              ("VInflCatAndAgrFmtvMax.n", "count", None),
              ("VInflCatFmtvMax.n", "count", None),
              ("VInflCatandAgrFmtvMax.binned3", None, None), # duplicate
              ("VAnyIncorporation.Presence", "bin", None),
              ("VInflExponence", "cat", None),
              ("VNounIncorporation.Presence", "bin", None),
              ("VPhonCoherence.Presence", "bin", None),
              ("VAgrAllPre.Presence", "bin", None),
              ("VAgrAnyPre.Presence", "bin", None),
              ("VProsCoherence.Presence", "bin", None),
              ("VAgrRoles", None, None), # 28 types
              ("VInflMax.n", "count", None),
              ("VSynCoherence.Presence", "bin", None),
              ("VInflMax.binned3", None, None), # duplicate
              ("VIncorporationVorN.Presence", "bin", None),
              ("VIncorporationV.Presence", "bin", None),
          ],
        },
        { "path": "Valence_classes_per_language.csv",
          "features": [
              ("ValClassesWithVerbsOfType.ability.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.achieve.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.activity.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.appearance_disappearance.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.aspectual.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.assuming_position.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.attempt.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.bodily_emission.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.caused_motion.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.change_of_possession.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.change_of_state.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.cognition.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.combine_attach.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.communication.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.creation_transformation.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.desideration.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.destroy.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.diseases_and_bodily_states.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.emotion.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.existence.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.fail.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.fight.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.grooming.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.happening.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.help.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.hide.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.hit.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.hold_keep.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.hurt.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.influence_cause.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.ingestion.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.judge.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.kill.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.lose.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.laugh_wink.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_contact.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_effective_action.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_emotion.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_motion.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_other.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_perception_cognition.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_pursuit.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.M_sensation.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.measure.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.modality.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.motion.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.name.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.necessity.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.obligation.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.perception.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.persuit.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.possession.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.push_pull.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.put_pour_spray_load_fill.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.quality.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.remove_wipe_clear.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.scent_emission.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.seem.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.send_carry_bring_take.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.sensation.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.separate_disassemble.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.sleep.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.social_interaction.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.sound_emission.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.spatial_position.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.spread_apply.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.touch.Presence", "bin", None),
              ("ValClassesWithVerbsOfType.use.Presence", "bin", None),
          ],
        },
        { "path": "VAgreement.csv",
          "features": [
              ("VAgrA.Presence", "bin", None),
              ("VAgrP.Presence", "bin", None),
              ("VAgrPOSS.Presence", "bin", None),
          ],
        },
        { "path": "VInfl_categories.csv",
          "features": [
              ("VInflAktionsart.Presence", "bin", None),
              ("VInflApplicative.Presence", "bin", None),
              ("VInflAspect.Presence", "bin", None),
              ("VInflCausative.Presence", "bin", None),
              ("VInflClassifier.Presence", "bin", None),
              ("VInflConnective.Presence", "bin", None),
              ("VInflConstruct.Presence", "bin", None),
              ("VInflControl.Presence", "bin", None),
              ("VInflDeixis.Presence", "bin", None),
              ("VInflVoice.Presence", "bin", None),
              ("VInflFocus.Presence", "bin", None),
              ("VInflPosterior.Presence", "bin", None),
              ("VInflGender.Presence", "bin", None),
              ("VInflIllocution.Presence", "bin", None),
              ("VInflMA.Presence", "bin", None),
              ("VInflMN.Presence", "bin", None),
              ("VInflModality.Presence", "bin", None),
              ("VInflMood.Presence", "bin", None),
              ("VInflMotion.Presence", "bin", None),
              ("VInflMSE.Presence", "bin", None),
              ("VInflNominalizer.Presence", None, None),
              ("VInflNumber.Presence", "bin", None),
              ("VInflPerson.Presence", "bin", None),
              ("VInflPolarity.Presence", "bin", None),
              ("VInflPotentialis.Presence", "bin", None),
              ("VInflQuantificational.Presence", "bin", None),
              ("VInflReflexive_and_Reciprocal.Presence", "bin", None),
              ("VInflRECIP.Presence", "bin", None),
              ("VInflReferential.Presence", "bin", None),
              ("VInflReflexive.Presence", "bin", None),
              ("VInflRepetition.Presence", "bin", None),
              ("VInflReversative.Presence", None, None),
              ("VInflSemistem.Presence", None, None),
              ("VInflSpatial.Presence", "bin", None),
              ("VInflStatus.Presence", "bin", None),
              ("VInflSwitch_reference.Presence", "bin", None),
              ("VInflTense_and_Evidential.Presence", "bin", None),
              ("VInflTA.Presence", "bin", None),
              ("VInflTAM.Presence", "bin", None),
              ("VInflTense.Presence", "bin", None),
              ("VInflTM.Presence", "bin", None),
              ("VInflValence.Presence", "bin", None),
              ("VInflVersion.Presence", "bin", None),
          ],
        },
        { "path": "VInfl_counts_per_position.csv",
          "features": [
              ("VInflCatAndAgrIn.n", "count", None),
              ("VInflCatAndAgrPost.n", "count", None),
              ("VInflCatAndAgrPrae.n", "count", None),
              ("VInflCatAndAgrSimul.n", "count", None),
              ("VInflCatAndAgrSplit.n", "count", None),
          ],
        },
        { "path": "VInfl_macrocategories.csv",
          "features": [
              ("macro.VInflClassification.Presence", "bin", None), # 2 types
              ("macro.VInflEvidential.Presence", "bin", None),
              ("macro.VInflInter_clausal.Presence", "bin", None),
              ("macro.VInflNP_related.Presence", "bin", None), # 2 types
              ("macro.VInflNumber.Presence", "bin", None),
              ("macro.VInflOperators.Presence", "bin", None),
              ("macro.VInflPerson.Presence", "bin", None), # 2 types
              ("macro.VInflPragmatic.Presence", "bin", None),
              ("macro.VInflPragmatics.Presence", "bin", None), # 2 types
              ("macro.VInflRole.Presence", "bin", None),
              ("macro.VInflTAM.Presence", "bin", None),
              ("macro.VInflValence.Presence", "bin", None),
          ],
        },
        { "path": "VInfl_cat_postposed.csv",
          "features": [
              ("VInflAktionsart.Post.Presence", "bin", None),
              ("VInflApplicative.Post.Presence", "bin", None), # 2
              ("VInflAspect.Post.Presence", "bin", None),
              ("VInflCausative.Post.Presence", "bin", None),
              ("VInflClassifier.Post.Presence", None, None), # 1 type
              ("VInflConnective.Post.Presence", "bin", None), # 2
              ("VInflConstruct.Post.Presence", "bin", None), # 2
              ("VInflControl.Post.Presence", None, None), # 1 type
              ("VInflDeixis.Post.Presence", "bin", None), # 2
              ("VInflVoice.Post.Presence", "bin", None),
              ("VInflFocus.Post.Presence", "bin", None),
              ("VInflPosterior.Post.Presence", "bin", None), # 2
              ("VInflGender.Post.Presence", "bin", None), # 2
              ("VInflIllocution.Post.Presence", "bin", None), # 2
              ("VInflMA.Post.Presence", "bin", None), # 2
              ("VInflMN.Post.Presence", "bin", None), # 2
              ("VInflModality.Post.Presence", "bin", None), # 2
              ("VInflMood.Post.Presence", "bin", None),
              ("VInflMotion.Post.Presence", None, None), # 1 type
              ("VInflMSE.Post.Presence", "bin", None), # 2
              ("VInflNominalizer.Post.Presence", None, None), # 1 type
              ("VInflNumber.Post.Presence", "bin", None),
              ("VInflPerson.Post.Presence", "bin", None), # 2
              ("VInflPolarity.Post.Presence", "bin", None),
              ("VInflPotentialis.Post.Presence", None, None), # 1 type
              ("VInflQuantificational.Post.Presence", "bin", None),
              ("VInflReflexive_and_Reciprocal.Post.Presence", "bin", None),
              ("VInflRECIP.Post.Presence", "bin", None),
              ("VInflReferential.Post.Presence", "bin", None), # 2
              ("VInflReflexive.Post.Presence", "bin", None),
              ("VInflRepetition.Post.Presence", "bin", None),
              ("VInflReversative.Post.Presence", None, None), # 1 type
              ("VInflSemistem.Post.Presence", None, None), # 1 type
              ("VInflSpatial.Post.Presence", "bin", None), # 2
              ("VInflStatus.Post.Presence", "bin", None),
              ("VInflSwitch_reference.Post.Presence", "bin", None), # 2
              ("VInflTense_and_Evidential.Post.Presence", "bin", None), # 2
              ("VInflTA.Post.Presence", "bin", None),
              ("VInflTAM.Post.Presence", "bin", None),
              ("VInflTense.Post.Presence", "bin", None),
              ("VInflTM.Post.Presence", "bin", None), # 2
              ("VInflValence.Post.Presence", "bin", None),
              ("VInflVersion.Post.Presence", None, None), # 1 type
          ],
        },
        { "path": "VInfl_cat_positions4.csv",
          "features": [
              ("VInflAktionsart.Position.binned4", "cat", None),
              ("VInflApplicative.Position.binned4", "cat", None), # 2
              ("VInflAspect.Position.binned4", "cat", None),
              ("VInflCausative.Position.binned4", "cat", None),
              ("VInflClassifier.Position.binned4", None, None), # 1 type
              ("VInflConnective.Position.binned4", "cat", None), # 2
              ("VInflConstruct.Position.binned4", "cat", None), # 2
              ("VInflControl.Position.binned4", None, None), # 1 type
              ("VInflDeixis.Position.binned4", "cat", None), # 2
              ("VInflVoice.Position.binned4", "cat", None),
              ("VInflFocus.Position.binned4", "cat", None),
              ("VInflPosterior.Position.binned4", "cat", None), # 2
              ("VInflGender.Position.binned4", "cat", None), # 2
              ("VInflIllocution.Position.binned4", "cat", None), # 2
              ("VInflMA.Position.binned4", "cat", None), # 2
              ("VInflMN.Position.binned4", "cat", None), # 2
              ("VInflModality.Position.binned4", "cat", None), # 2
              ("VInflMood.Position.binned4", "cat", None),
              ("VInflMotion.Position.binned4", None, None), # 1 type
              ("VInflMSE.Position.binned4", "cat", None), # 2
              ("VInflNominalizer.Position.binned4", None, None), # 1 type
              ("VInflNumber.Position.binned4", "cat", None),
              ("VInflPerson.Position.binned4", "cat", None), # 2
              ("VInflPolarity.Position.binned4", "cat", None),
              ("VInflPotentialis.Position.binned4", None, None), # 1 type
              ("VInflQuantificational.Position.binned4", "cat", None),
              ("VInflReflexive_and_Reciprocal.Position.binned4", "cat", None),
              ("VInflRECIP.Position.binned4", "cat", None),
              ("VInflReferential.Position.binned4", "cat", None), # 2
              ("VInflReflexive.Position.binned4", "cat", None),
              ("VInflRepetition.Position.binned4", "cat", None),
              ("VInflReversative.Position.binned4", None, None), # 1 type
              ("VInflSemistem.Position.binned4", None, None), # 1 type
              ("VInflSpatial.Position.binned4", "cat", None), # 2
              ("VInflStatus.Position.binned4", "cat", None),
              ("VInflSwitch_reference.Position.binned4", "cat", None), # 2
              ("VInflTense_and_Evidential.Position.binned4", "cat", None), # 2
              ("VInflTA.Position.binned4", "cat", None),
              ("VInflTAM.Position.binned4", "cat", None),
              ("VInflTense.Position.binned4", "cat", None),
              ("VInflTM.Position.binned4", "cat", None), # 2
              ("VInflValence.Position.binned4", "cat", None),
              ("VInflVersion.Position.binned4", None, None), # 1 type
          ],
        },
        { "path": "VInfl_cat_positions.csv",
          "features": [  # do not create "modbin" because they will be removed for low coverage
              ("VInflAktionsart.Position", "cat", None),
              ("VInflApplicative.Position", "cat", None), # 2
              ("VInflAspect.Position", "cat", None),
              ("VInflCausative.Position", "cat", None),
              ("VInflClassifier.Position", None, None), # 1 type
              ("VInflConnective.Position", "cat", None), # 2
              ("VInflConstruct.Position", "cat", None), # 2
              ("VInflControl.Position", None, None), # 1 type
              ("VInflDeixis.Position", "cat", None), # 2
              ("VInflVoice.Position", "cat", None),
              ("VInflFocus.Position", "cat", None),
              ("VInflPosterior.Position", "cat", None), # 2
              ("VInflGender.Position", "cat", None), # 2
              ("VInflIllocution.Position", "cat", None), # 2
              ("VInflMA.Position", "cat", None), # 2
              ("VInflMN.Position", "cat", None), # 2
              ("VInflModality.Position", "cat", None), # 2
              ("VInflMood.Position", "cat", None),
              ("VInflMotion.Position", None, None), # 1 type
              ("VInflMSE.Position", "cat", None), # 2
              ("VInflNominalizer.Position", None, None), # 1 type
              ("VInflNumber.Position", "cat", None),
              ("VInflPerson.Position", "cat", None), # 2
              ("VInflPolarity.Position", "cat", None),
              ("VInflPotentialis.Position", None, None), # 1 type
              ("VInflQuantificational.Position", "cat", None),
              ("VInflReflexive_and_Reciprocal.Position", "cat", None),
              ("VInflRECIP.Position", "cat", None),
              ("VInflReferential.Position", "cat", None), # 2
              ("VInflReflexive.Position", "cat", None),
              ("VInflRepetition.Position", "cat", None),
              ("VInflReversative.Position", None, None), # 1 type
              ("VInflSemistem.Position", None, None), # 1 type
              ("VInflSpatial.Position", "cat", None), # 2
              ("VInflStatus.Position", "cat", None),
              ("VInflSwitch_reference.Position", "cat", None), # 2
              ("VInflTense_and_Evidential.Position", "cat", None), # 2
              ("VInflTA.Position", "cat", None),
              ("VInflTAM.Position", "cat", None),
              ("VInflTense.Position", "cat", None),
              ("VInflTM.Position", "cat", None), # 2
              ("VInflValence.Position", "cat", None),
              ("VInflVersion.Position", None, None), # 1 type
          ],
        },
        { "path": "VInfl_cat_positions5.csv",
          "features": [
              ("VInflAktionsart.Position.binned5", "cat", None),
              ("VInflApplicative.Position.binned5", "cat", None), # 2
              ("VInflAspect.Position.binned5", "cat", None),
              ("VInflCausative.Position.binned5", "cat", None),
              ("VInflClassifier.Position.binned5", None, None), # 1 type
              ("VInflConnective.Position.binned5", "cat", None), # 2
              ("VInflConstruct.Position.binned5", "cat", None), # 2
              ("VInflControl.Position.binned5", None, None), # 1 type
              ("VInflDeixis.Position.binned5", "cat", None), # 2
              ("VInflVoice.Position.binned5", "cat", None),
              ("VInflFocus.Position.binned5", "cat", None),
              ("VInflPosterior.Position.binned5", "cat", None), # 2
              ("VInflGender.Position.binned5", "cat", None), # 2
              ("VInflIllocution.Position.binned5", "cat", None), # 2
              ("VInflMA.Position.binned5", "cat", None), # 2
              ("VInflMN.Position.binned5", "cat", None), # 2
              ("VInflModality.Position.binned5", "cat", None), # 2
              ("VInflMood.Position.binned5", "cat", None),
              ("VInflMotion.Position.binned5", None, None), # 1 type
              ("VInflMSE.Position.binned5", "cat", None), # 2
              ("VInflNominalizer.Position.binned5", None, None), # 1 type
              ("VInflNumber.Position.binned5", "cat", None),
              ("VInflPerson.Position.binned5", "cat", None), # 2
              ("VInflPolarity.Position.binned5", "cat", None),
              ("VInflPotentialis.Position.binned5", None, None), # 1 type
              ("VInflQuantificational.Position.binned5", "cat", None),
              ("VInflReflexive_and_Reciprocal.Position.binned5", "cat", None),
              ("VInflRECIP.Position.binned5", "cat", None),
              ("VInflReferential.Position.binned5", "cat", None), # 2
              ("VInflReflexive.Position.binned5", "cat", None),
              ("VInflRepetition.Position.binned5", "cat", None),
              ("VInflReversative.Position.binned5", None, None), # 1 type
              ("VInflSemistem.Position.binned5", None, None), # 1 type
              ("VInflSpatial.Position.binned5", "cat", None), # 2
              ("VInflStatus.Position.binned5", "cat", None),
              ("VInflSwitch_reference.Position.binned5", "cat", None), # 2
              ("VInflTense_and_Evidential.Position.binned5", "cat", None), # 2
              ("VInflTA.Position.binned5", "cat", None),
              ("VInflTAM.Position.binned5", "cat", None),
              ("VInflTense.Position.binned5", "cat", None),
              ("VInflTM.Position.binned5", "cat", None), # 2
              ("VInflValence.Position.binned5", "cat", None),
              ("VInflVersion.Position.binned5", None, None), # 1 type
          ],
        },
        { "path": "VInfl_cat_preposed.csv",
          "features": [
              ("VInflAktionsart.Pre.Presence", "bin", None),
              ("VInflApplicative.Pre.Presence", "bin", None), # 2
              ("VInflAspect.Pre.Presence", "bin", None),
              ("VInflCausative.Pre.Presence", "bin", None),
              ("VInflClassifier.Pre.Presence", None, None), # 1 type
              ("VInflConnective.Pre.Presence", "bin", None), # 2
              ("VInflConstruct.Pre.Presence", "bin", None), # 2
              ("VInflControl.Pre.Presence", None, None), # 1 type
              ("VInflDeixis.Pre.Presence", "bin", None), # 2
              ("VInflVoice.Pre.Presence", "bin", None),
              ("VInflFocus.Pre.Presence", "bin", None),
              ("VInflPosterior.Pre.Presence", "bin", None), # 2
              ("VInflGender.Pre.Presence", "bin", None), # 2
              ("VInflIllocution.Pre.Presence", "bin", None), # 2
              ("VInflMA.Pre.Presence", "bin", None), # 2
              ("VInflMN.Pre.Presence", "bin", None), # 2
              ("VInflModality.Pre.Presence", "bin", None), # 2
              ("VInflMood.Pre.Presence", "bin", None),
              ("VInflMotion.Pre.Presence", None, None), # 1 type
              ("VInflMSE.Pre.Presence", "bin", None), # 2
              ("VInflNominalizer.Pre.Presence", None, None), # 1 type
              ("VInflNumber.Pre.Presence", "bin", None),
              ("VInflPerson.Pre.Presence", "bin", None), # 2
              ("VInflPolarity.Pre.Presence", "bin", None),
              ("VInflPotentialis.Pre.Presence", "bin", None), # 1 type
              ("VInflQuantificational.Pre.Presence", "bin", None),
              ("VInflReflexive_and_Reciprocal.Pre.Presence", "bin", None),
              ("VInflRECIP.Pre.Presence", "bin", None),
              ("VInflReferential.Pre.Presence", "bin", None), # 2
              ("VInflReflexive.Pre.Presence", "bin", None),
              ("VInflRepetition.Pre.Presence", "bin", None),
              ("VInflReversative.Pre.Presence", None, None), # 1 type
              ("VInflSemistem.Pre.Presence", None, None), # 1 type
              ("VInflSpatial.Pre.Presence", "bin", None), # 2
              ("VInflStatus.Pre.Presence", "bin", None),
              ("VInflSwitch_reference.Pre.Presence", "bin", None), # 2
              ("VInflTense_and_Evidential.Pre.Presence", "bin", None), # 2
              ("VInflTA.Pre.Presence", "bin", None),
              ("VInflTAM.Pre.Presence", "bin", None),
              ("VInflTense.Pre.Presence", "bin", None),
              ("VInflTM.Pre.Presence", "bin", None), # 2
              ("VInflValence.Pre.Presence", "bin", None),
              ("VInflVersion.Pre.Presence", "bin", None), # 1 type
          ],
        },
        { "path": "VInfl_cat_multiexponence.csv",
          "features": [   # 1 type per feature
              ("multi.VInflAktionsart.Pre.Presence", "bin", None),
              ("multi.VInflApplicative.Pre.Presence", "bin", None),
              ("multi.VInflAspect.Pre.Presence", "bin", None),  # bin-min 1
              ("multi.VInflCausative.Pre.Presence", "bin", None),
              ("multi.VInflClassifier.Pre.Presence", "bin", None),
              ("multi.VInflConnective.Pre.Presence", "bin", None),
              ("multi.VInflConstruct.Pre.Presence", "bin", None),
              ("multi.VInflControl.Pre.Presence", "bin", None),
              ("multi.VInflDeixis.Pre.Presence", "bin", None),
              ("multi.VInflVoice.Pre.Presence", "bin", None),
              ("multi.VInflFocus.Pre.Presence", "bin", None),
              ("multi.VInflPosterior.Pre.Presence", "bin", None),
              ("multi.VInflGender.Pre.Presence", "bin", None),
              ("multi.VInflIllocution.Pre.Presence", "bin", None),
              ("multi.VInflMA.Pre.Presence", "bin", None),
              ("multi.VInflMN.Pre.Presence", "bin", None),
              ("multi.VInflModality.Pre.Presence", "bin", None),
              ("multi.VInflMood.Pre.Presence", "bin", None),
              ("multi.VInflMotion.Pre.Presence", "bin", None),
              ("multi.VInflMSE.Pre.Presence", "bin", None),
              ("multi.VInflNominalizer.Pre.Presence", "bin", None),
              ("multi.VInflNumber.Pre.Presence", "bin", None),
              ("multi.VInflPerson.Pre.Presence", "bin", None),
              ("multi.VInflPolarity.Pre.Presence", "bin", None),
              ("multi.VInflPotentialis.Pre.Presence", "bin", None),
              ("multi.VInflQuantificational.Pre.Presence", "bin", None),
              ("multi.VInflReflexive_and_Reciprocal.Pre.Presence", "bin", None),
              ("multi.VInflRECIP.Pre.Presence", "bin", None),
              ("multi.VInflReferential.Pre.Presence", "bin", None),
              ("multi.VInflReflexive.Pre.Presence", "bin", None),
              ("multi.VInflRepetition.Pre.Presence", "bin", None),
              ("multi.VInflReversative.Pre.Presence", "bin", None),
              ("multi.VInflSemistem.Pre.Presence", "bin", None),
              ("multi.VInflSpatial.Pre.Presence", "bin", None),
              ("multi.VInflStatus.Pre.Presence", "bin", None),
              ("multi.VInflSwitch_reference.Pre.Presence", "bin", None),
              ("multi.VInflTense_and_Evidential.Pre.Presence", "bin", None),
              ("multi.VInflTA.Pre.Presence", "bin", None),
              ("multi.VInflTAM.Pre.Presence", "bin", None),  # bin-min 2
              ("multi.VInflTense.Pre.Presence", "bin", None),
              ("multi.VInflTM.Pre.Presence", "bin", None),
              ("multi.VInflValence.Pre.Presence", "bin", None),
              ("multi.VInflVersion.Pre.Presence", "bin", None),
          ],
        },
        { "path": "VInfl_macrocat_postposed.csv",
          "features": [
              ("macro.VInflClassification.Post.Presence", "bin", None), # 2 types
              ("macro.VInflEvidential.Post.Presence", None, None),
              ("macro.VInflInter_clausal.Post.Presence", "bin", None),
              ("macro.VInflNP_related.Post.Presence", "bin", None), # 2 types
              ("macro.VInflNumber.Post.Presence", "bin", None),
              ("macro.VInflOperators.Post.Presence", "bin", None),
              ("macro.VInflPerson.Post.Presence", "bin", None), # 2 types
              ("macro.VInflPragmatic.Post.Presence", "bin", None),
              ("macro.VInflPragmatics.Post.Presence", "bin", None), # 2 types
              ("macro.VInflRole.Post.Presence", None, None),
              ("macro.VInflTAM.Post.Presence", "bin", None),
              ("macro.VInflValence.Post.Presence", "bin", None),
          ],
        },
        { "path": "VInfl_macrocat_position4.csv",
          "features": [
              ("macro.VInflClassification.Position.binned4", "cat", None), # 2 types
              ("macro.VInflEvidential.Position.binned4", None, None),
              ("macro.VInflInter_clausal.Position.binned4", "cat", None),
              ("macro.VInflNP_related.Position.binned4", "cat", None), # 2 types
              ("macro.VInflNumber.Position.binned4", "cat", None),
              ("macro.VInflOperators.Position.binned4", "cat", None),
              ("macro.VInflPerson.Position.binned4", "cat", None), # 2 types
              ("macro.VInflPragmatic.Position.binned4", "cat", None),
              ("macro.VInflPragmatics.Position.binned4", "cat", None), # 2 types
              ("macro.VInflRole.Position.binned4", None, None),
              ("macro.VInflTAM.Position.binned4", "cat", None),
              ("macro.VInflValence.Position.binned4", "cat", None),
          ],
        },
        { "path": "VInfl_macrocat_position.csv",
          "features": [
              ("macro.VInflClassification.Position", "cat", None), # 2 types
              ("macro.VInflEvidential.Position", None, None),
              ("macro.VInflInter_clausal.Position", "cat", None),
              ("macro.VInflNP_related.Position", "cat", None), # 2 types
              ("macro.VInflNumber.Position", "cat", None),
              ("macro.VInflOperators.Position", "cat", None),
              ("macro.VInflPerson.Position", "cat", None), # 2 types
              ("macro.VInflPragmatics.Position", "cat", None), # 2 types
              ("macro.VInflRole.Position", None, None),
              ("macro.VInflTAM.Position", "cat", None),
              ("macro.VInflValence.Position", "cat", None),
          ],
        },
        { "path": "VInfl_macrocat_position5.csv",
          "features": [
              ("macro.VInflClassification.Position.binned5", "cat", None), # 2 types
              ("macro.VInflEvidential.Position.binned5", None, None),
              ("macro.VInflInter_clausal.Position.binned5", "cat", None),
              ("macro.VInflNP_related.Position.binned5", "cat", None), # 2 types
              ("macro.VInflNumber.Position.binned5", "cat", None),
              ("macro.VInflOperators.Position.binned5", "cat", None),
              ("macro.VInflPerson.Position.binned5", "cat", None), # 2 types
              ("macro.VInflPragmatic.Position.binned5", "cat", None),
              ("macro.VInflPragmatics.Position.binned5", "cat", None), # 2 types
              ("macro.VInflRole.Position.binned5", None, None),
              ("macro.VInflTAM.Position.binned5", "cat", None),
              ("macro.VInflValence.Position.binned5", "cat", None),
          ],
        },
        { "path": "VInfl_macrocat_preposed.csv",
          "features": [
              ("macro.VInflClassification.Pre.Presence", "bin", None), # 2 types
              ("macro.VInflEvidential.Pre.Presence", None, None),
              ("macro.VInflInter_clausal.Pre.Presence", "bin", None),
              ("macro.VInflNP_related.Pre.Presence", "bin", None), # 2 types
              ("macro.VInflNumber.Pre.Presence", "bin", None),
              ("macro.VInflOperators.Pre.Presence", "bin", None),
              ("macro.VInflPerson.Pre.Presence", "bin", None), # 2 types
              ("macro.VInflPragmatic.Pre.Presence", "bin", None),
              ("macro.VInflPragmatics.Pre.Presence", "bin", None), # 2 types
              ("macro.VInflRole.Pre.Presence", None, None),
              ("macro.VInflTAM.Pre.Presence", "bin", None),
              ("macro.VInflValence.Pre.Presence", "bin", None),
          ],
        },
        { "path": "VInfl_macrocat_multiexponence.csv",
          "features": [ # 1 type per feature
              ("macromulti.VInflClassification.Pre.Presence", None, None),
              ("macromulti.VInflEvidential.Pre.Presence", None, None),
              ("macromulti.VInflInter_clausal.Pre.Presence", None, None),
              ("macromulti.VInflNP_related.Pre.Presence", None, None),
              ("macromulti.VInflNumber.Pre.Presence", None, None), # bin-min 1
              ("macromulti.VInflOperators.Pre.Presence", "bin", None),
              ("macromulti.VInflPerson.Pre.Presence", None, None),
              ("macromulti.VInflPragmatic.Pre.Presence", None, None),
              ("macromulti.VInflPragmatics.Pre.Presence", None, None),
              ("macromulti.VInflRole.Pre.Presence", None, None),
              ("macromulti.VInflTAM.Pre.Presence", None, None), # bin-min 3
              ("macromulti.VInflValence.Pre.Presence", None, None),
          ],
        },
        { "path": "VAgr_postposed.csv",
          "features": [
              ("VAgrA.Post.Presence", "bin", None),
              ("VAgrP.Post.Presence", "bin", None),
              ("VAgrPOSS.Post.Presence", None, None), # low-coverage
          ],
        },
        { "path": "VAgr_position4.csv",
          "features": [
              ("VAgrA.Position.binned4", "cat", None),
              ("VAgrP.Position.binned4", "cat", None),
              ("VAgrPOSS.Position.binned4", None, None), # low-coverage
          ],
        },
        { "path": "VAgr_position.csv",
          "features": [
              ("VAgrA.Position", "cat", None),
              ("VAgrP.Position", "cat", None),
              ("VAgrPOSS.Position", None, None), # low-coverage
              
          ],
        },
        { "path": "VAgr_position5.csv",
          "features": [
              ("VAgrA.Position.binned5", "cat", None),
              ("VAgrP.Position.binned5", "cat", None),
              ("VAgrPOSS.Position.binned5", None, None), # low-coverage
          ],
        },
        { "path": "VAgr_preposed.csv",
          "features": [
              ("VAgrA.Pre.Presence", "bin", None),
              ("VAgrP.Pre.Presence", "bin", None),
              ("VAgrPOSS.Pre.Presence", None, None), # low-coverage
          ],
        },
        { "path": "VAgr_multiexponence.csv",
          "features": [
              ("multi.VAgrA.Pre.Presence", "bin", None),
              ("multi.VAgrP.Pre.Presence", None, None),    # bin-min 2
              ("multi.VAgrPOSS.Pre.Presence", None, None), # all FALSE
          ],
        },
    ]
    flist = []
    fname2struct = {}
    for spec in specs:
        flist_each = {}
        for i, (k, ftype, vmap) in enumerate(spec["features"]):
            if ftype is None:
                continue
            fstruct = { "size": -1, "annotation": { "name": k, "source": spec["path"], "orig_idx": i, "total": 0, "count": defaultdict(int) }}
            if k in fname2struct:
                sys.stderr.write("non-unique feature name\t{}\n".format(k))
                exit(1)
            fname2struct[k] = fstruct
            flist.append(fstruct)
            flist_each[k] = fstruct
            if ftype == "bin":
                fstruct["type"] = "bin"
                fstruct["size"] = 1
            elif ftype == "modbin":
                fstruct["type"] = "bin"
                fstruct["size"] = 1
                assert(len(vmap) == 2)
                fstruct["annotation"]["vid2label"] = [None, None]
                fstruct["annotation"]["label2vid"] = {}
                for vk, vv in vmap.items():
                    if vv:
                        fstruct["annotation"]["vid2label"][1] = vk
                        fstruct["annotation"]["label2vid"][vk] = 1
                    else:
                        fstruct["annotation"]["vid2label"][0] = vk
                        fstruct["annotation"]["label2vid"][vk] = 0
            elif ftype in ("count", "pcount"):
                fstruct["type"] = "count"
                fstruct["size"] = 1                
            elif ftype == "cat":
                fstruct["type"] = "cat"
            else:
                sys.stderr.write("unsupported ftype\t{}\n".format(ftype))
                exit(1)

        with open(os.path.join(args.basedir, "data", spec["path"]), newline="", encoding='utf-8') as f:
            reader = csv.reader(f)
            next(reader) # skip legend
            dat = list(reader)
        for row in dat:
            lname = row.pop(0)
            if lname not in name2lang:
                continue
            lang = name2lang[lname]
            for (k, ftype, vmap), v in zip(spec["features"], row):
                if len(v) <= 0:
                    continue
                if ftype is None:
                    # ad-hoc amendments
                    if k == "Gender.Presence" and v == "FALSE":
                        v = 0
                        lang["annotation"]["features"]["Gender.n"] = v
                        # assume no double counting
                        fstruct2 = fname2struct["Gender.n"]
                        fstruct2["annotation"]["count"][v] += 1
                        fstruct2["annotation"]["total"] += 1

                    continue
                fstruct = flist_each[k]

                if ftype == "bin":
                    if v == "TRUE":
                        v = 1
                    elif v == "FALSE":
                        v = 0
                    else:
                        sys.stderr.write("not a binary value\t{}\n".format(v))
                        v = 0
                elif ftype in ("count", "pcount"):
                    v = int(v)
                    if ftype == "pcount":
                        v -= 1
                elif ftype == "modbin":
                    v = fstruct["annotation"]["label2vid"][v]
                elif ftype == "cat":
                    # temporarily assign raw value
                    pass

                if k in lang["annotation"]["features"]:
                    sys.stderr.write("multiple entries for the same datapoint\t{}\t{}\t{} == {}?\n".format(lang["annotation"]["name"], k, lang["annotation"]["features"][k], v))
                else:
                    lang["annotation"]["features"][k] = v
                    fstruct["annotation"]["count"][v] += 1
                    fstruct["annotation"]["total"] += 1


    # ordered counts
    for fstruct in flist:
        if not fstruct["type"] == "count":
            continue
        fstruct["annotation"]["vid2label"] = []
        fstruct["annotation"]["label2vid"] = {}
        for i, v in enumerate(sorted(fstruct["annotation"]["count"].keys())):
            fstruct["annotation"]["label2vid"][str(v)] = i
            fstruct["annotation"]["vid2label"].append(str(v))


    # assign vids to categorical features
    for fstruct in flist:
        if not fstruct["type"] == "cat":
            continue
        count = fstruct["annotation"]["count"]
        fstruct["size"] = len(count)
        fstruct["annotation"]["vid2label"] = [None] * len(count)
        fstruct["annotation"]["label2vid"] = {}
        for vid, k in enumerate(sorted(count.keys())):
            fstruct["annotation"]["label2vid"][k] = vid
            fstruct["annotation"]["vid2label"][vid] = k
        name = fstruct["annotation"]["name"]
        for lang in langs:
            if name in lang["annotation"]["features"]:
                k = lang["annotation"]["features"][name]
                # sys.stderr.write("{}\t{}\t{}\n".format(lang["annotation"]["name"], name, k))
                lang["annotation"]["features"][name] = fstruct["annotation"]["label2vid"][k]


    # remove features below frequency threshold
    sys.stderr.write("{} features before filtering\n".format(len(flist)))
    flist2 = []
    binsize = 0
    for fstruct in flist:
        is_good = True
        name = fstruct["annotation"]["name"]
        if fstruct["annotation"]["total"] < args.fthres:
            is_good = False
            sys.stderr.write("remove low-coverage feature\t{}\n".format(name))
        elif len(fstruct["annotation"]["count"]) < 2:
            is_good = False
            sys.stderr.write("remove single-valued feature\t{}\n".format(name))
        elif fstruct["type"] == "bin" and min(fstruct["annotation"]["count"].values()) < args.binminthres:
            is_good = False
            sys.stderr.write("remove binary feature whose minor value is below the frequency threshold\t{}\n".format(name))
        elif fstruct["type"] == "cat" and len(fstruct["annotation"]["count"]) <= 2:
            sys.stderr.write("WARNING: categorical feature should have more than two possible values\t{}\n".format(name))
        if is_good:
            fstruct["fid"] = len(flist2)
            flist2.append(fstruct)
            binsize += fstruct["size"]
        else:
            for lang in langs:
                if name in lang["annotation"]["features"]:
                    del lang["annotation"]["features"][name]
    flist = flist2
    sys.stderr.write("{} (binsize: {}) features after filtering\n".format(len(flist), binsize))


    # remove featureless language
    langs2 = []
    phylos = defaultdict(int)
    filled = 0
    for lang in langs:
        filled2 = len(lang["annotation"]["features"])
        filled += filled2
        if filled2 <= 0:
            pass
            # sys.stderr.write("remove featureless language\t{}\n".format(lang["annotation"]["name"]))
        else:
            langs2.append(lang)
            phylos[lang["phylo"]] += 1
    langs = langs2
    sys.stderr.write("{} languages remained\n".format(len(langs)))
    sys.stderr.write("{} phylogenetic groups\n".format(len(phylos)))

    total = len(langs) * len(flist)
    sys.stderr.write("coverage\t{} ({} / {})\n".format(filled / total, filled, total))

    # create catvect
    k2fid = {}
    for fstruct in flist:
        k2fid[fstruct["annotation"]["name"]] = fstruct["fid"]
    for lang in langs:
        catvect = [-1] * len(flist)
        for k, v in lang["annotation"]["features"].items():
            catvect[k2fid[k]] = v
        lang["catvect"] = catvect


    with open(args.flist, 'w') as f:
        f.write("%s\n" % json.dumps(flist, indent=4, sort_keys=True))

    with open(args.langs, 'w') as f:
        for lang in langs:
            f.write("%s\n" % json.dumps(lang))

if __name__ == "__main__":
    main()
