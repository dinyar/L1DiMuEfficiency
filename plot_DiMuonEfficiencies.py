#!/usr/bin/python

import sys
import os
import argparse

desc = ''
parser = argparse.ArgumentParser(description=desc,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--qualityBasedCOU', default='False', action='store_true',
                    help='Flag to enable cancelling the muon with lower \
                    quality. Otherwise the muon with higher pT is cancelled.')
parser.add_argument('--outDir', type=str, default='plots',
                    help='Folder to store plots in.')
opts = parser.parse_args()

sys.path.append(os.path.join(os.path.dirname(__file__),
                             "../L1AnalysisHelpers"))
from CreateHistograms import (binningDict, generateCombinedGhostPercHist,
                              generateCombinedEfficiencyHist)

from ROOT import gROOT, kTRUE

gROOT.Reset()
gROOT.SetBatch(kTRUE)

# Cut dicts
genCuts = {}
genCuts["mu-pt1"] = ["(pT1_gen > 1)", "mu-ptGen1"]
genCuts["diMu-pt1"] = ["((pT1_gen > 1) && (pT2_gen > 1))", "diMu-ptGen1"]
gmtCuts = {}
gmtCuts["gmt_diMu-pt1"] = ["((pT1 > 1) && (pT2 > 1))",
                           "diMu-pt1"]

if opts.qualityBasedCOU is True:
    ghostSelector = "_q"
else:
    ghostSelector = "_pt"

gmtCuts["ugmt_diMu-pt1"] = ["((pT1" + ghostSelector + " > 1) && (pT2" +
                            ghostSelector + " > 1))",
                            "diMu-pt1"]

gmtCuts["bmtf"] = ["(tfType1" + ghostSelector + "==0)", "bmtf"]
gmtCuts["omtf"] = ["(tfType1" + ghostSelector + "==1)", "omtf"]
gmtCuts["emtf"] = ["(tfType1" + ghostSelector + "==2)", "emtf"]

gmtCuts["singleBmtf"] = ["((tfType1" + ghostSelector + "==0) && (tfType2" +
                         ghostSelector + "==-11))", "bmtf"]
gmtCuts["singleBmtfEtaFine"] = ["((tfType1" + ghostSelector +
                                "==0) && (tfType2" + ghostSelector +
                                "==-11) && (hf1" + ghostSelector + "==1))",
                                "bmtfEtaFine"]
gmtCuts["singleBmtfEtaCoarse"] = ["((tfType1" + ghostSelector +
                                  "==0) && (tfType2" + ghostSelector +
                                  "==-11) && (hf1" + ghostSelector + "==0))",
                                  "bmtfEtaCoarse"]
gmtCuts["singleOmtf"] = ["((tfType1" + ghostSelector + "==1) && (tfType2" +
                         ghostSelector + "==-11))", "omtf"]
gmtCuts["singleEmtf"] = ["((tfType1" + ghostSelector + "==2) && (tfType2" +
                         ghostSelector + "==-11))", "emtf"]

gmtCuts["diBmtf"] = ["((tfType1" + ghostSelector + "==0) && (tfType2" +
                     ghostSelector + "==0))", "diBmtf"]
gmtCuts["diOmtf"] = ["((tfType1" + ghostSelector + "==1) && (tfType2" +
                     ghostSelector + "==1))", "diOmtf"]
gmtCuts["diEmtf"] = ["((tfType1" + ghostSelector + "==2) && (tfType2" +
                     ghostSelector + "==2))", "diEmtf"]
gmtCuts["diBOmtf"] = ["(((tfType1" + ghostSelector + "==0) && (tfType2" +
                      ghostSelector + "==1)) || ((tfType1" + ghostSelector +
                      "==1) && (tfType2" + ghostSelector + "==0)))", "diBOmtf"]
gmtCuts["diBOmtfCoarse"] = ["(((tfType1" + ghostSelector + "==0) && (hf1" +
                            ghostSelector + "==0) && (tfType2" +
                            ghostSelector + "==1)) || ((tfType1" +
                            ghostSelector + "==1) && (tfType2" +
                            ghostSelector + "==0) && (hf2" + ghostSelector +
                            "==0)))", "diBOmtf"]
gmtCuts["diBOmtfFine"] = ["(((tfType1" + ghostSelector + "==0) && (hf1" +
                          ghostSelector + "==1) && (tfType2" + ghostSelector +
                          "==1)) || ((tfType1" + ghostSelector +
                          "==1) && (tfType2" + ghostSelector +
                          "==0) && (hf2" + ghostSelector + "==1)))",
                          "diBOmtf"]
gmtCuts["diEOmtf"] = ["(((tfType1" + ghostSelector + "==1) && (tfType2" +
                      ghostSelector + "==2)) || ((tfType1" + ghostSelector +
                      "==2) && (tfType2" + ghostSelector + "==1)))", "diEOmtf"]
gmtCuts["diBEmtf"] = ["(((tfType1" + ghostSelector + "==0) && (tfType2" +
                      ghostSelector + "==2)) || ((tfType1" + ghostSelector +
                      "==2) && (tfType2" + ghostSelector + "==0)))", "diBEmtf"]


# TODO: Do these plots also for maxQual!

efficiencyList = []
# TODO: For mu1/mu2 plot single mu efficiencies?
# Entries: Label for histogram (Will be used for filename and title) |
# binning | parameters used for project functions
efficiencyList.append([["deltaEta_gen", "#Delta#eta(#mu^{-}#mu^{+})"],
                       binningDict["distWide"],
                       "abs(eta1_gen-eta2_gen)",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["deltaPhi_gen", "#Delta#phi(#mu^{-}#mu^{+})"],
                       binningDict["distWide"],
                       "abs(phi1_gen-phi2_gen)",
                       genCuts["diMu-pt1"], [0, 1.4]])
efficiencyList.append([["deltaR_gen", "#DeltaR(#mu^{-}#mu^{+})"],
                       binningDict["distWide"],
                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2)",
                       genCuts["diMu-pt1"], [0, 1.4]])
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

gmt_singleMu_file = "GMTSingleMuNtuple.root"
gmt_dimu_file = "GMTDimuonNtuple.root"
ugmt_singleMu_file = "uGMTSingleMuNtuple.root"
ugmt_dimu_file = "uGMTDimuonNtuple.root"

jpsi_ntuples = []
jpsi_ntuples.append(gmt_dimu_file)
jpsi_ugmt_ntuples = []
# jpsi_ugmt_ntuples.append(ugmt_dimu_file)
# jpsi_ugmt_ntuples.append("uGMTDimuonNtuple-dR0_3.root")
# jpsi_ugmt_ntuples.append("uGMTDimuonNtuple-dR0_1.root")
# jpsi_ugmt_ntuples.append("uGMTDimuonNtuple-dR0_1-BMTFOMTFchargeMatch.root")
# jpsi_ugmt_ntuples.append("uGMTDimuonNtuple-dR0_05.root")
jpsi_ugmt_ntuples.append(
    "uGMTDimuonNtuple-dPhi0_05dEta0_1-BOMTF_dEtaFine0_1-dEtaCoarse0_3-EOMTF_dEta0_1-EMTF_dEta0_05.root")
jpsi_ugmt_ntuples.append(
    "uGMTDimuonNtuple-dPhi0_05dEta0_1-BMTF_OMTF_cM-BOMTF_dEtaFine0_1-dEtaCoarse0_3_cM-EOMTF_dEta0_1-EMTF_dEta0_05.root")

jpsi_ntuples.extend(jpsi_ugmt_ntuples)

ntuple_names = []
ntuple_names.append("gmt_ntuple")
ntuple_names.extend(len(jpsi_ugmt_ntuples) * ["ugmt_ntuple"])
labels = []
labels.append(["Gen muons", "GMT muons", "GMT"])
# labels.append(["Gen muons", "uGMT muons", "uGMT"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.3", "uGMT",
#                "dR0-3"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.1", "uGMT",
#                "dR0-1"])
# labels.append(["Gen muons",
#                "uGMT muons w/ cancel-out #DeltaR<0.1, \
# match charges in BMTF+OMTF", "uGMT",
#                "dR0-1-BMTFOMTFchargeMatch"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.05", "uGMT",
#                "dR0-05"])

# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01", "uGMT",
#                "dR0-01"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01,\
# match charges and #DeltaR<0.1 in OMTF", "uGMT",
#                "dR0-01_OMTF-dR0-1-chargeMatch"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01,\
# match charges and #DeltaR<0.3 in OMTF", "uGMT",
#                "dR0-01_OMTF-dR0-3-chargeMatch"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01,\
# #DeltaR<0.3 in OMTF", "uGMT",
#                "dR0-01_OMTF-dR0-3"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01,\
# #DeltaR<0.3 in BMTF/OMTF, #DeltaR<0.1 in EMTF/OMTF", "uGMT",
#                "dR0-01-BOMTF_dR0_3-EOMTF_dR0_1"])
# labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01,\
# match charges and #DeltaR<0.3 in BMTF/OMTF, #DeltaR<0.1 in EMTF/OMTF", "uGMT",
#                "dR0-01-BOMTF_dR0_3_chargeMatch-EOMTF_dR0_1"])
labels.append(["Gen muons", "uGMT muons w/ cancel-out #Delta#phi<0.05,\
#Delta#eta<0.1; #DeltaEta<0.05 in EMTF/OMTF, \
#DeltaEta<0.1 in  EMTF, #Delta#eta<0.3 in BMTF if no eta fine",
               "uGMT",
               "dPhi0_05dEta0_1-BOMTF_dEta0_3-EOMTF_dEta0_1-EMTF_dEta0_05"])
labels.append(["Gen muons", "uGMT muons w/ cancel-out #Delta#phi<0.05,\
#Delta#eta<0.1; charge match in BMTF+OMTF, charge match and #DeltaEta<0.05\
 in EMTF/OMTF, #DeltaEta<0.1 in EMTF, #Delta#eta<0.3 in BMTF if no eta \
fine",
               "uGMT",
               "dPhi0_05dEta0_1-BMTF_OMTF_cM-BOMTF_dEta0_3_cM-EOMTF_dEta0_1-EMTF_dEta0_05"])

line_colours = []
line_colours.append(1)
line_colours.append(46)
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
cuts.append(gmtCuts["gmt_diMu-pt1"])
cuts.extend(len(jpsi_ugmt_ntuples) * [gmtCuts["ugmt_diMu-pt1"]])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, jpsi_ntuples, ntuple_names,
                                   labels, line_colours, cuts, "jPsi",
                                   rootFolder=opts.outDir)


ccntuple = []
ccntuple.extend(3 * [ugmt_dimu_file])
ccntuple_name = []
ccntuple_name.extend(3 * ["ugmt_ntuple"])
ccdlabel = []
ccdlabel.append(["All muons", "BMTF muons", "uGMT"])
ccdlabel.append(["All muons", "OMTF muons", "uGMT"])
ccdlabel.append(["All muons", "EMTF muons", "uGMT"])
cclc = []
cclc.append(30)
cclc.append(36)
cclc.append(48)
cccuts = []
cccuts.append(gmtCuts["bmtf"])
cccuts.append(gmtCuts["omtf"])
cccuts.append(gmtCuts["emtf"])
chargeCheckList = []
chargeCheckList.append([["mu1_ch", "ch"],
                        binningDict["charge"], "ch1" + ghostSelector,
                        genCuts["diMu-pt1"], [0, 1.4]])
for varList in chargeCheckList:
    generateCombinedEfficiencyHist(varList, ccntuple, ccntuple_name,
                                   ccdlabel, cclc,
                                   cccuts, "charge_check", drawGenMus=True,
                                   drawStackPlot=True, rootFolder=opts.outDir)

ghostListWOgmt = []
# TODO: Add plot with deltaPt.
# NOTE: If no L1 muon at all deltaX will be 0! This will lead to inefficiencies
# TODO: Because of note above exclude events without an L1 muon.
# (we're not looking at efficiency, but ghost rate so we'd be underestimating
# ghost rate otherwise.)
# basically we want to know which percentage of singleMu events would be
# "poisoned" by ghosts
ghostListWOgmt.append([["deltaEta_L1", "#Delta#eta(#mu#mu_{Ghost})"],
                       binningDict["distVeryWide"],
                       "abs(eta1" + ghostSelector + "-eta2" + ghostSelector +
                       ")",
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["deltaPhi_L1", "#Delta#phi(#mu#mu_{Ghost})"],
                       binningDict["distVeryWide"],
                       "abs(phi1" + ghostSelector + "-phi2" + ghostSelector +
                       ")",
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["deltaR_L1", "#DeltaR(#mu#mu_{Ghost})"],
                       binningDict["distVeryWide"],
                       "sqrt((eta1" + ghostSelector + "-eta2" + ghostSelector +
                       ")**2+(phi1" + ghostSelector + "-phi2" + ghostSelector +
                       ")**2)",
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["deltaEta_L1-zoom", "#Delta#eta(#mu#mu_{Ghost})"],
                       binningDict["distNarrow"],
                       "abs(eta1" + ghostSelector + "-eta2" + ghostSelector +
                       ")",
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["deltaPhi_L1-zoom", "#Delta#phi(#mu#mu_{Ghost})"],
                       binningDict["distNarrow"],
                       "abs(phi1" + ghostSelector + "-phi2" + ghostSelector +
                       ")",
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["deltaR_L1-zoom", "#DeltaR(#mu#mu_{Ghost})"],
                       binningDict["distNarrow"],
                       "sqrt((eta1" + ghostSelector + "-eta2" +
                       ghostSelector + ")**2+(phi1" + ghostSelector +
                       "-phi2" + ghostSelector + ")**2)",
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["mu1_L1Eta", "#eta(leading #mu_{L1})"],
                       binningDict["etaFineRestr"], "eta1" + ghostSelector,
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["mu1_L1Phi", "#phi(leading #mu_{L1})"],
                       binningDict["phiFineRestr"], "phi1" + ghostSelector,
                       genCuts["mu-pt1"], [0, 0.3]])
ghostListWOgmt.append([["mu1_L1Pt", "p_{T}(leading #mu_{L1}) [GeV/c]"],
                       binningDict["pt140Fine"], "pT1" + ghostSelector,
                       genCuts["mu-pt1"], [0, 0.3]])

ghostListWgmt = []
ghostListWgmt.append([["mu1_genEta", "#eta(#mu)"],
                      binningDict["etaFineRestr"], "eta1_gen",
                      genCuts["mu-pt1"], [0, 0.3]])
ghostListWgmt.append([["mu1_genPhi", "#phi(#mu)"],
                      binningDict["phiFineRestr"], "phi1_gen",
                      genCuts["mu-pt1"], [0, 0.3]])
ghostListWgmt.append([["mu1_genPt", "p_{T}(#mu) [GeV/c]"],
                      binningDict["pt140Fine"], "pT1_gen",
                      genCuts["mu-pt1"], [0, 0.3]])

singleMu_ntuples = []
singleMu_ntuples.append(gmt_singleMu_file)
# singleMu_ntuples.append(ugmt_singleMu_file)
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_3.root")
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_1.root")
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_1-BMTFOMTFchargeMatch.root")
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_05.root")
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_01.root")

# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_01-OMTF_dR0_1_chargeMatch.root")
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_01-OMTF_dR0_3_chargeMatch.root")
# singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_01-OMTF_dR0_3.root")
# singleMu_ntuples.append(
#     "uGMTSingleMuNtuple-dR0_01-BOMTF_dR0_3-EOMTF_dR0_1.root")
# singleMu_ntuples.append(
#     "uGMTSingleMuNtuple-dR0_01-BOMTF_dR0_3_chargeMatch-EOMTF_dR0_1.root")
singleMu_ntuples.append(
    "uGMTSingleMuNtuple-dPhi0_05dEta0_1-BOMTF_dEtaFine0_1-dEtaCoarse0_3-EOMTF_dEta0_1-EMTF_dEta0_05.root")
singleMu_ntuples.append(
    "uGMTSingleMuNtuple-dPhi0_05dEta0_1-BMTF_OMTF_cM-BOMTF_dEtaFine0_1-dEtaCoarse0_3_cM-EOMTF_dEta0_1-EMTF_dEta0_05.root")

for varList in ghostListWgmt:
    generateCombinedGhostPercHist(varList, singleMu_ntuples, ntuple_names,
                                  labels, line_colours, cuts, "singleMu",
                                  rootFolder=opts.outDir)
for varList in ghostListWOgmt:
    generateCombinedGhostPercHist(varList, singleMu_ntuples[1:],
                                  ntuple_names[1:], labels[1:],
                                  line_colours[1:], cuts[1:], "singleMu",
                                  rootFolder=opts.outDir)


resolution_check_ntuple = []
resolution_check_ntuple.extend(5 * [ugmt_singleMu_file])
resolution_check_ntuple_name = []
resolution_check_ntuple_name.extend(
    len(resolution_check_ntuple) * ["ugmt_ntuple"])
resolution_check_dlabel = []
resolution_check_dlabel.append(["All muons", "BMTF muons", "uGMT"])
resolution_check_dlabel.append(["All muons", "BMTF muons, fine eta", "uGMT"])
resolution_check_dlabel.append(["All muons", "BMTF muons, coarse eta", "uGMT"])
resolution_check_dlabel.append(["All muons", "OMTF muons", "uGMT"])
resolution_check_dlabel.append(["All muons", "EMTF muons", "uGMT"])
resolution_check_line_colour = []
resolution_check_line_colour.append(28)
resolution_check_line_colour.append(30)
resolution_check_line_colour.append(32)
resolution_check_line_colour.append(36)
resolution_check_line_colour.append(48)
resolution_check_cuts = []
resolution_check_cuts.append(gmtCuts["singleBmtf"])
resolution_check_cuts.append(gmtCuts["singleBmtfEtaFine"])
resolution_check_cuts.append(gmtCuts["singleBmtfEtaCoarse"])
resolution_check_cuts.append(gmtCuts["singleOmtf"])
resolution_check_cuts.append(gmtCuts["singleEmtf"])
resolutionCheckList = []
resolutionCheckList.append([["phiResolution", "#Delta#phi(#mu_{L1}#mu_{Gen})"],
                            binningDict["distSym"], "phi1" + ghostSelector +
                            "-phi1_gen",
                            genCuts["mu-pt1"], [0, 1.4]])
resolutionCheckList.append([["etaResolution", "#Delta#eta(#mu_{L1}#mu_{Gen})"],
                            binningDict["distSym"], "eta1" + ghostSelector +
                            "-eta1_gen",
                            genCuts["mu-pt1"], [0, 1.4]])
resolutionCheckList.append([["mu1_l1Eta", "#eta(#mu_{L1})"],
                            binningDict["etaFineRestr"], "eta1" +
                            ghostSelector,
                            genCuts["mu-pt1"], [0, 1.4]])
resolutionCheckList.append([["mu1_genEta", "#eta(#mu_{Gen})"],
                            binningDict["etaFineRestr"], "eta1_gen",
                            genCuts["mu-pt1"], [0, 1.4]])
for varList in resolutionCheckList:
    generateCombinedEfficiencyHist(varList, resolution_check_ntuple,
                                   resolution_check_ntuple_name,
                                   resolution_check_dlabel,
                                   resolution_check_line_colour,
                                   resolution_check_cuts, "resolution_check",
                                   drawGenMus=False, drawStackPlot=True,
                                   rootFolder=opts.outDir)

ghost_distance_ntuple = []
ghost_distance_ntuple.extend(5 * [ugmt_singleMu_file])
ghost_distance_ntuple_name = []
ghost_distance_ntuple_name.extend(5 * ["ugmt_ntuple"])
ghost_distance_label = []
ghost_distance_label.append(["All muons", "BMTF muons", "uGMT"])
ghost_distance_label.append(["All muons", "OMTF muons", "uGMT"])
ghost_distance_label.append(["All muons", "EMTF muons", "uGMT"])
ghost_distance_label.append(["All muons", "BMTF+OMTF overlap muons", "uGMT"])
ghost_distance_label.append(["All muons", "EMTF+OMTF overlap muons", "uGMT"])
ghost_distance_line_colour = []
ghost_distance_line_colour.append(30)
ghost_distance_line_colour.append(36)
ghost_distance_line_colour.append(48)
ghost_distance_line_colour.append(8)
ghost_distance_line_colour.append(28)
ghost_distance_cuts = []
ghost_distance_cuts.append(gmtCuts["diBmtf"])
ghost_distance_cuts.append(gmtCuts["diOmtf"])
ghost_distance_cuts.append(gmtCuts["diEmtf"])
ghost_distance_cuts.append(gmtCuts["diBOmtf"])
ghost_distance_cuts.append(gmtCuts["diEOmtf"])
ghostDistanceList = []
ghostDistanceList.append([["phiResolution", "#Delta#phi(#mu_{L1}#mu_{Ghost})"],
                          binningDict["distSym"],
                          "phi1" + ghostSelector + "-phi2" + ghostSelector,
                          genCuts["mu-pt1"], [0, 1.4]])
ghostDistanceList.append([["etaResolution", "#Delta#eta(#mu_{L1}#mu_{Ghost})"],
                          binningDict["distSym"],
                          "eta1" + ghostSelector + "-eta2" + ghostSelector,
                          genCuts["mu-pt1"], [0, 1.4]])
ghostDistanceList.append([["ptResolution",
                           "#Delta p_{T}(#mu_{L1}#mu_{Ghost})"],
                          binningDict["pt140Fine"],
                          "pT1" + ghostSelector + "-pT2" + ghostSelector,
                          genCuts["mu-pt1"], [0, 1.4]])
for varList in ghostDistanceList:
    generateCombinedEfficiencyHist(varList, ghost_distance_ntuple,
                                   ghost_distance_ntuple_name,
                                   ghost_distance_label,
                                   ghost_distance_line_colour,
                                   ghost_distance_cuts, "ghost_distance",
                                   drawGenMus=False, drawStackPlot=True,
                                   rootFolder=opts.outDir)

tf_eff_ntuples = []
tf_eff_ntuples.append(gmt_dimu_file)
tf_eff_ntuples.extend(7 * [ugmt_dimu_file])
# Di BMTF/OMTF/EMTF, BMTF+OMTF, OMTF+EMTF, BMTF+EMTF (cross-check)
tf_eff_ntuple_names = []
tf_eff_ntuple_names.append("gmt_ntuple")
tf_eff_ntuple_names.extend(7 * ["ugmt_ntuple"])
tf_eff_labels = []
tf_eff_labels.append(["Gen muons", "GMT muons", "GMT"])
tf_eff_labels.append(["Gen muons", "only BMTF muons", "uGMT", "diBMTF"])
tf_eff_labels.append(["Gen muons", "only OMTF muons", "uGMT", "diOMTF"])
tf_eff_labels.append(["Gen muons", "only EMTF muons", "uGMT", "diEMTF"])
tf_eff_labels.append(["Gen muons", "BMTF coarse+OMTF muons", "uGMT",
                      "diBOMTFcoarse"])
tf_eff_labels.append(["Gen muons", "BMTF fine+OMTF muons", "uGMT",
                      "diBOMTFfine"])
tf_eff_labels.append(["Gen muons", "OMTF+EMTF muons", "uGMT", "diOEMTF"])
tf_eff_labels.append(["Gen muons", "BMTF+EMTF muons", "uGMT", "diBEMTF"])
tf_eff_line_colours = []
tf_eff_line_colours.append(1)
tf_eff_line_colours.append(46)
tf_eff_line_colours.append(30)
tf_eff_line_colours.append(38)
tf_eff_line_colours.append(8)
tf_eff_line_colours.append(28)
tf_eff_line_colours.append(17)
tf_eff_line_colours.append(7)
tf_eff_cuts = []
tf_eff_cuts.append(gmtCuts["gmt_diMu-pt1"])
tf_eff_cuts.append(gmtCuts["diBmtf"])
tf_eff_cuts.append(gmtCuts["diOmtf"])
tf_eff_cuts.append(gmtCuts["diEmtf"])
tf_eff_cuts.append(gmtCuts["diBOmtfCoarse"])
tf_eff_cuts.append(gmtCuts["diBOmtfFine"])
tf_eff_cuts.append(gmtCuts["diEOmtf"])
tf_eff_cuts.append(gmtCuts["diBEmtf"])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, tf_eff_ntuples,
                                   tf_eff_ntuple_names, tf_eff_labels,
                                   tf_eff_line_colours, tf_eff_cuts, "tf_eff",
                                   drawGenMus=True, drawStackPlot=True,
                                   rootFolder=opts.outDir)


tf_ghosts_ntuples = []
tf_ghosts_ntuples.append(gmt_singleMu_file)
tf_ghosts_ntuples.extend(7 * [ugmt_singleMu_file])

for varList in ghostListWgmt:
    generateCombinedGhostPercHist(varList, tf_ghosts_ntuples,
                                  tf_eff_ntuple_names, tf_eff_labels,
                                  tf_eff_line_colours, tf_eff_cuts,
                                  "tf_ghosts", drawGenMus=False,
                                  drawStackPlot=True, rootFolder=opts.outDir)
for varList in ghostListWOgmt:
    generateCombinedGhostPercHist(varList, tf_ghosts_ntuples[1:],
                                  tf_eff_ntuple_names[1:], tf_eff_labels[1:],
                                  tf_eff_line_colours[1:], tf_eff_cuts[1:],
                                  "tf_ghosts", drawGenMus=False,
                                  drawStackPlot=True, rootFolder=opts.outDir)


tf_eff_w_gb_ntuples = []
tf_eff_w_gb_ntuples.append(gmt_dimu_file)
tf_eff_w_gb_ntuples.extend(
    7 * ["uGMTDimuonNtuple-dPhi0_05dEta0_1-BOMTF_dEtaFine0_1-dEtaCoarse0_3-EOMTF_dEta0_1-EMTF_dEta0_05.root"])
tf_eff_w_gb_labels = []
tf_eff_w_gb_labels.append(["Gen muons", "GMT muons", "GMT"])
tf_eff_w_gb_labels.append(["Gen muons", "uGMT w/ cancel-out, only BMTF muons",
                           "uGMT", "diBMTF"])
tf_eff_w_gb_labels.append(["Gen muons", "uGMT w/ cancel-out, only OMTF muons",
                           "uGMT", "diOMTF"])
tf_eff_w_gb_labels.append(["Gen muons", "uGMT w/ cancel-out, only EMTF muons",
                           "uGMT", "diEMTF"])
tf_eff_w_gb_labels.append(["Gen muons",
                           "uGMT w/ cancel-out, BMTF coarse+OMTF muons",
                           "uGMT", "diBOMTFcoarse"])
tf_eff_w_gb_labels.append(["Gen muons",
                           "uGMT w/ cancel-out, BMTF fine+OMTF muons", "uGMT",
                           "diBOMTFfine"])
tf_eff_w_gb_labels.append(["Gen muons", "uGMT w/ cancel-out, OMTF+EMTF muons",
                           "uGMT", "diOEMTF"])
tf_eff_w_gb_labels.append(["Gen muons", "uGMT w/ cancel-out, BMTF+EMTF muons",
                           "uGMT", "diBEMTF"])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, tf_eff_w_gb_ntuples,
                                   tf_eff_ntuple_names, tf_eff_w_gb_labels,
                                   tf_eff_line_colours, tf_eff_cuts,
                                   "tf_eff_w_gb", drawGenMus=True,
                                   drawStackPlot=True, rootFolder=opts.outDir)


tf_ghosts_w_gb_ntuples = []
tf_ghosts_w_gb_ntuples.append(gmt_singleMu_file)
tf_ghosts_w_gb_ntuples.extend(
    7 * ["uGMTSingleMuNtuple-dPhi0_05dEta0_1-BOMTF_dEtaFine0_1-dEtaCoarse0_3-EOMTF_dEta0_1-EMTF_dEta0_05.root"])

for varList in ghostListWgmt:
    generateCombinedGhostPercHist(varList, tf_ghosts_w_gb_ntuples,
                                  tf_eff_ntuple_names, tf_eff_w_gb_labels,
                                  tf_eff_line_colours, tf_eff_cuts,
                                  "tf_ghosts_w_gb", drawGenMus=False,
                                  drawStackPlot=True, rootFolder=opts.outDir)
for varList in ghostListWOgmt:
    generateCombinedGhostPercHist(varList, tf_ghosts_w_gb_ntuples[1:],
                                  tf_eff_ntuple_names[1:],
                                  tf_eff_w_gb_labels[1:],
                                  tf_eff_line_colours[1:], tf_eff_cuts[1:],
                                  "tf_ghosts_w_gb", drawGenMus=False,
                                  drawStackPlot=True, rootFolder=opts.outDir)
