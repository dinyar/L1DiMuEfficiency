#!/usr/bin/python

import sys
import os
import argparse

desc = ''
parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--outDir', type=str, default='plots',
                    help='Folder to store plots in.')
opts = parser.parse_args()

sys.path.append(os.path.join(os.path.dirname(__file__),
                             "../../../L1AnalysisHelpers"))
from CreateHistograms import (binningDict, generateCombinedGhostPercHist,
                              generateCombinedEfficiencyHist,
                              generate2DRateHist)

from ROOT import gROOT, kTRUE, kGreen, kBlack

gROOT.Reset()
gROOT.SetBatch(kTRUE)

# Ntuples to use.
gmt_dimu_file = "legacy_gmt/GMTDimuonNtuple.root"
ugmt_dimu_file = "20161205_defaultTuning/uGMTDimuonNtuple.root"

# Cut dicts
genCuts = {}
genCuts["diMu-pt1"] = ["((pT1 > 1) && (pT2 > 1))", "diMu-ptGen1"]

genCuts["diMu-jpsiPtLow"] = ["((pT1 > 1) && (pT2 > 1)) && (pT_jpsi < 50)", "diMu-jpsiPtLow"]

gmtCuts = {}
gmtCuts["gmt_diMu-pt1_os"] = ["(pT1 > 1) && (pT2 > 1) && (ch1 != ch2)",
                                 "diMu-pt1_os"]
gmtCuts["ugmt_diMu-pt1_os"] = ["(pT1 > 1) && (pT2 > 1) && (ch1 != ch2)",
                                  "diMu-pt1_os"]


# twoDlist = []
# twoDlist.append(["dRvsPt_input", binningDict["pt140Fine"], binningDict["distWide"],
#                  "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2):pT_jpsi",
#                  gmtCuts["ugmt_diMu-pt1_q3"],
#                  ["#DeltaR(#mu^{-}#mu^{+})", "p_{T}(J/#Psi) [GeV/c]"]])
# twoDlist.append(["dRvsPt1_input", binningDict["pt140Fine"], binningDict["distWide"],
#                  "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2):pT1_gen",
#                  gmtCuts["ugmt_diMu-pt1_q3"],
#                  ["#DeltaR(#mu^{-}#mu^{+})", "p_{T}(leading #mu) [GeV/c]"]])
# twoDlist.append(["dRvsPt2_input", binningDict["pt140Fine"], binningDict["distWide"],
#                  "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2):pT2_gen",
#                  gmtCuts["ugmt_diMu-pt1_q3"],
#                  ["#DeltaR(#mu^{-}#mu^{+})", "p_{T}(trailing #mu) [GeV/c]"]])
# twoDlist.append(["Pt1vsPt2_input", binningDict["pt140Fine"], binningDict["pt140Fine"],
#                  "pT1_gen:pT2_gen",
#                  genCuts["diMu-pt1"],
#                  ["p_{T}(leading #mu) [GeV/c]", "p_{T}(trailing #mu) [GeV/c]"]])
# 
# for varList in twoDlist:
#     generate2DRateHist(varList, ugmt_dimu_file, "tf_ntuple")
# 
# twoDlist_ugmt = []
# twoDlist_ugmt.append(["dRvsPt_output", binningDict["pt140Fine"], binningDict["distWide"],
#                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2):pT_jpsi",
#                       gmtCuts["ugmt_diMu-pt1_q3"],
#                       ["#DeltaR(#mu^{-}#mu^{+})", "p_{T}(J/#Psi) [GeV/c]"]])
# twoDlist_ugmt.append(["dRvsPt1_output", binningDict["pt140Fine"], binningDict["distWide"],
#                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2):pT1_gen",
#                       gmtCuts["ugmt_diMu-pt1_q3"],
#                       ["#DeltaR(#mu^{-}#mu^{+})", "p_{T}(leading #mu) [GeV/c]"]])
# twoDlist_ugmt.append(["dRvsPt2_output", binningDict["pt140Fine"], binningDict["distWide"],
#                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2):pT2_gen",
#                       gmtCuts["ugmt_diMu-pt1_q3"],
#                       ["#DeltaR(#mu^{-}#mu^{+})", "p_{T}(trailing #mu) [GeV/c]"]])
# twoDlist_ugmt.append(["Pt1vsPt2_output", binningDict["pt140Fine"], binningDict["pt140Fine"],
#                       "pT1_gen:pT2_gen",
#                       genCuts["diMu-pt1"],
#                       ["p_{T}(leading #mu) [GeV/c]", "p_{T}(trailing #mu) [GeV/c]"]])
# 
# for varList in twoDlist_ugmt:
#     generate2DRateHist(varList, ugmt_dimu_file, "ugmt_ntuple")

efficiencyList = []
# Entries: Label for histogram (Will be used for filename and title) |
# binning | parameters used for project functions
efficiencyList.append([["mu1_genEta", "#eta(leading #mu)"],
                       binningDict["etaFineRestr"], "eta1_gen",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["mu2_genEta", "#eta(trailing #mu)"],
                       binningDict["etaFineRestr"], "eta2_gen",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["mu1_genPhi", "#phi(leading #mu)"],
                       binningDict["phiFineRestr"], "phi1_gen",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["mu2_genPhi", "#phi(trailing #mu)"],
                       binningDict["phiFineRestr"], "phi2_gen",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["mu1_genPt", "p_{T}(leading #mu) [GeV/c]"],
                       binningDict["pt140Fine"], "pT1_gen",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["mu2_genPt", "p_{T}(trailing #mu) [GeV/c]"],
                       binningDict["pt140Fine"], "pT2_gen",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["jPsi_genEta", "#eta(J/#Psi)"],
                       binningDict["etaFineRestr"], "eta_jpsi",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["jPsi_genPhi", "#phi(J/#Psi)"],
                       binningDict["phiFineRestr"], "phi_jpsi",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["jPsi_genPt", "p_{T}(J/#Psi) [GeV/c]"],
                       binningDict["pt140Fine"], "pT_jpsi",
                       genCuts["diMu-pt1"], [0, 1.4]])


# Plot efficiency for in- and ouput of uGMT w/o splitting into TF contributions

ugmt_inout_labels = []
ugmt_inout_labels.append(["RECO muons", "Legacy GMT, opposite sign", "GMT"])
ugmt_inout_labels.append(["RECO muons", "uGMT, opposite sign", "uGMT"])
jpsi_efficiency_ntuples = []
jpsi_efficiency_ntuples.append(gmt_dimu_file)
jpsi_efficiency_ntuples.append(ugmt_dimu_file)
ntuple_names = []
ntuple_names.append("gmt_ntuple")
ntuple_names.append("ugmt_ntuple")

line_colours = []
line_colours.append(kGreen+2)
line_colours.append(kBlack)
line_colours.append(30)
line_colours.append(38)
line_colours.append(8)
line_colours.append(28)
line_colours.append(7)
line_colours.append(9)
line_colours.append(32)
line_colours.append(35)
line_colours.append(42)
cuts = []
cuts.append(gmtCuts["gmt_diMu-pt1_os"])
cuts.append(gmtCuts["ugmt_diMu-pt1_os"])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, jpsi_efficiency_ntuples, ntuple_names,
                                   ugmt_inout_labels, line_colours, cuts,
                                   "jPsi_efficiency",
                                   rootFolder=opts.outDir, drawDistributions=False)

efficiencyList_lowPt = []
# Entries: Label for histogram (Will be used for filename and title) |
# binning | parameters used for project functions
efficiencyList_lowPt.append([["mu1_genEta", "#eta(leading #mu)"],
                       binningDict["etaFineRestr"], "eta1_gen",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["mu2_genEta", "#eta(trailing #mu)"],
                       binningDict["etaFineRestr"], "eta2_gen",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["mu1_genPhi", "#phi(leading #mu)"],
                       binningDict["phiFineRestr"], "phi1_gen",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["mu2_genPhi", "#phi(trailing #mu)"],
                       binningDict["phiFineRestr"], "phi2_gen",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["mu1_genPt", "p_{T}(leading #mu) [GeV/c]"],
                       binningDict["pt50Fine"], "pT1_gen",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["mu2_genPt", "p_{T}(trailing #mu) [GeV/c]"],
                       binningDict["pt50Fine"], "pT2_gen",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["jPsi_genEta", "#eta(J/#Psi)"],
                       binningDict["etaFineRestr"], "eta_jpsi",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["jPsi_genPhi", "#phi(J/#Psi)"],
                       binningDict["phiFineRestr"], "phi_jpsi",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])
efficiencyList_lowPt.append([["jPsi_genPt", "p_{T}(J/#Psi) [GeV/c]"],
                       binningDict["pt50Fine"], "pT_jpsi",
                       genCuts["diMu-jpsiPtLow"], [0, 1.4]])

for varList in efficiencyList_lowPt:
    generateCombinedEfficiencyHist(varList, jpsi_efficiency_ntuples, ntuple_names,
                                   ugmt_inout_labels, line_colours, cuts,
                                   "jPsi_efficiency_lowPt",
                                   rootFolder=opts.outDir, drawDistributions=False)

# # Plot efficiencies at uGMT ouput split by TF contributions.
# 
# tf_eff_labels = []
# tf_eff_labels.append(
#     ["RECO muons", "Legacy GMT", "GMT", "GMToutput"])
# tf_eff_labels.append(
#     ["RECO muons", "uGMT input, q>4", "TF", "uGMTinput"])
# tf_eff_labels.append(
#     ["RECO muons", "only BMTF muons, q>4", "uGMT", "diBMTF"])
# tf_eff_labels.append(
#     ["RECO muons", "only OMTF muons, q>4", "uGMT", "diOMTF"])
# tf_eff_labels.append(
#     ["RECO muons", "only EMTF muons, q>4", "uGMT", "diEMTF"])
# tf_eff_labels.append(["RECO muons", "BMTF+OMTF muons, q>4", "uGMT",
#                       "diBOMTF"])
# tf_eff_labels.append(
#     ["RECO muons", "OMTF+EMTF muons, q>4", "uGMT", "diOEMTF"])
# tf_eff_labels.append(
#     ["RECO muons", "BMTF+EMTF muons, q>4", "uGMT", "diBEMTF"])
# tf_eff_ntuples = []
# tf_eff_ntuples.append(gmt_dimu_file)
# tf_eff_ntuples.extend((len(tf_eff_labels) - 1) * [ugmt_dimu_file])
# tf_eff_ntuple_names = []
# tf_eff_ntuple_names.append("gmt_ntuple")
# tf_eff_ntuple_names.append("tf_ntuple")
# tf_eff_ntuple_names.extend((len(tf_eff_labels) - 2) * ["ugmt_ntuple"])
# tf_eff_line_colours = []
# tf_eff_line_colours.append(1)
# tf_eff_line_colours.append(46)
# tf_eff_line_colours.append(30)
# tf_eff_line_colours.append(38)
# tf_eff_line_colours.append(8)
# tf_eff_line_colours.append(28)
# tf_eff_line_colours.append(17)
# tf_eff_line_colours.append(7)
# tf_eff_cuts = []
# tf_eff_cuts.append(gmtCuts["gmt_diMu-pt1_q3"])
# tf_eff_cuts.append(gmtCuts["ugmt_diMu-pt1_q3"])
# tf_eff_cuts.append(gmtCuts["diBmtf_q4"])
# tf_eff_cuts.append(gmtCuts["diOmtf_q4"])
# tf_eff_cuts.append(gmtCuts["diEmtf_q4"])
# tf_eff_cuts.append(gmtCuts["diBOmtf_q4"])
# tf_eff_cuts.append(gmtCuts["diEOmtf_q4"])
# 
# for varList in efficiencyList:
#     generateCombinedEfficiencyHist(varList, tf_eff_ntuples,
#                                    tf_eff_ntuple_names, tf_eff_labels,
#                                    tf_eff_line_colours, tf_eff_cuts, "tf_eff",
#                                    drawGenMus=True, drawStackPlot=True,
#                                    rootFolder=opts.outDir)
# 
