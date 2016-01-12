#!/usr/bin/python

from ROOT import *
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../L1AnalysisHelpers"))
from CreateHistograms import *

gROOT.Reset()
gROOT.SetBatch(kTRUE)

# Cut dicts
gmtCuts = {}
genCuts = {}
gmtCuts["diMu-pt1"] = ["((pT1 > 1) && (pT2 > 1))", "diMu-pt1"]
# TODO: Try different distance cut for each track finder/area.
gmtCuts["diMu-pt1_separated"] = ["((pT1 > 1) && (pT2 > 1) && (sqrt((eta1-eta2)**2+(phi1-phi2)**2) > 0.1))", "diMu-pt1_separated"]
# TODO: Have cuts to only look at certain track finders
genCuts["diMu-pt1"] = ["((pT1_gen > 1) && (pT2_gen > 1))", "diMu-pt1"]

efficiencyList = []
# TODO: Axis labels, think about more descriptive title.
# TODO: Can combine/infer stuff: Efficiency implies GMT/reco;
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
# Current emulators compared with each other (i.e. uGMT without cancel out system!)
efficiencyList.append([["deltaEta_gen", ""],
                       binningDict["distNarrow"],
                       "abs(eta1_gen-eta2_gen)",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["deltaPhi_gen", ""],
                       binningDict["distNarrow"],
                       "abs(phi1_gen-phi2_gen)",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["deltaR_gen", ""],
                       binningDict["distNarrow"],
                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2)",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu1_genEta", ""],
                       binningDict["etaFine"], "eta1_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu2_genEta", ""],
                       binningDict["etaFine"], "eta2_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu1_genPhi", ""],
                       binningDict["phiFine"], "phi1_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu2_genPhi", ""],
                       binningDict["phiFine"], "phi2_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu1_genPt", ""],
                       binningDict["ptFine"], "pT1_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu2_genPt", ""],
                       binningDict["ptFine"], "pT2_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["jPsi_genEta", ""],
                       binningDict["etaFine"], "eta_jpsi",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["jPsi_genPhi", ""],
                       binningDict["phiFine"], "phi_jpsi",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["jPsi_genPt", ""],
                       binningDict["ptFine"], "pT_jpsi",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1"],
                       genCuts["diMu-pt1"]])

# Basic coordinate-based cancel-out
efficiencyList.append([["deltaEta_gen", ""],
                       binningDict["distNarrow"],
                       "abs(eta1_gen-eta2_gen)",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["deltaPhi_gen", ""],
                       binningDict["distNarrow"],
                       "abs(phi1_gen-phi2_gen)",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["deltaR_gen", ""],
                       binningDict["distNarrow"],
                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2)",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu1_genEta", ""],
                       binningDict["etaFine"],
                       "eta1_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu2_genEta", ""],
                       binningDict["etaFine"],
                       "eta2_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu1_genPhi", ""],
                       binningDict["phiFine"],
                       "phi1_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu2_genPhi", ""],
                       binningDict["phiFine"],
                       "phi2_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu1_genPt", ""],
                       binningDict["ptFine"],
                       "pT1_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
efficiencyList.append([["mu2_genPt", ""],
                       binningDict["ptFine"],
                       "pT2_gen",
                       gmtCuts["diMu-pt1"], gmtCuts["diMu-pt1_separated"],
                       genCuts["diMu-pt1"]])
# TODO: Add invariant mass calculation.


# TODO: Get these via argparse
ntuple_files = []
ntuple_files.append("GMTDimuonNtuple.root")
ntuple_files.append("uGMTDimuonNtuple.root")
ntuple_names = []
ntuple_names.append("gmt_ntuple")
ntuple_names.append("ugmt_ntuple")
datasets = []
datasets.append("GMT")
datasets.append("uGMT")
distribution_labels = []
distribution_labels.append(["Gen muons", "GMT muons"])
distribution_labels.append(["Gen muons", "uGMT muons"])
line_colours = []
line_colours.append(38)
line_colours.append(46)

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, ntuple_files, datasets,
                                   ntuple_names, distribution_labels,
                                   line_colours)
