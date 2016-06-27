#!/usr/bin/python

# TODO: Merge with plot_... from 808.
# TODO: Use TF muons to judge ghost distance
# TODO: Also use TF muons for efficiency and ghost plots
# TODO: Add dimuon quality cuts

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
                             "../L1AnalysisHelpers"))
from CreateHistograms import (binningDict, generateCombinedGhostPercHist,
                              generateCombinedEfficiencyHist)

from ROOT import gROOT, kTRUE

gROOT.Reset()
gROOT.SetBatch(kTRUE)

# Ntuples to use.
ugmt_singleMu_file = "uGMTSingleMuNtuple.root"
ugmt_dimu_file = "uGMTDimuonNtuple.root"

# Cut dicts
genCuts = {}
genCuts["mu-pt1"] = ["(pT1_gen > 1)", "mu-ptGen1"]
genCuts["diMu-pt1"] = ["((pT1_gen > 1) && (pT2_gen > 1))", "diMu-ptGen1"]


gmtCuts = {}
gmtCuts["ugmt_diMu-pt1"] = ["((pT1 > 1) && (pT2 > 1))",
                            "diMu-pt1"]

gmtCuts["bmtf"] = ["(tfType1==0)", "bmtf"]
gmtCuts["omtf"] = ["(tfType1==1)", "omtf"]
gmtCuts["emtf"] = ["(tfType1==2)", "emtf"]

gmtCuts["singleBmtf"] = ["((tfType1==0) && (tfType2==-10))", "bmtf"]
gmtCuts["singleBmtfEtaFine"] = ["((tfType1==0) && (tfType2==-10) && (hf1==1))",
                                "bmtfEtaFine"]
gmtCuts["singleBmtfEtaCoarse"] = ["((tfType1==0) && (tfType2==-10) && (hf1==0))",
                                  "bmtfEtaCoarse"]
gmtCuts["singleOmtf"] = ["((tfType1==1) && (tfType2==-10))", "omtf"]
gmtCuts["singleEmtf"] = ["((tfType1==2) && (tfType2==-10))", "emtf"]

gmtCuts["diBmtf"] = ["((tfType1==0) && (tfType2==0))", "diBmtf"]
gmtCuts["diOmtf"] = ["((tfType1==1) && (tfType2==1))", "diOmtf"]
gmtCuts["diEmtf"] = ["((tfType1==2) && (tfType2==2))", "diEmtf"]
gmtCuts["diBOmtf"] = [
    "(((tfType1==0) && (tfType2==1)) || ((tfType1==1) && (tfType2==0)))", "diBOmtf"]
gmtCuts["diBOmtfCoarse"] = [
    "(((tfType1==0) && (hf1==0) && (tfType2==1)) || ((tfType1==1) && (tfType2==0) && (hf2==0)))", "diBOmtf"]
gmtCuts["diBOmtfFine"] = ["(((tfType1==0) && (hf1==1) && (tfType2==1)) || ((tfType1==1) && (tfType2==0) && (hf2==1)))",
                          "diBOmtf"]
gmtCuts["diEOmtf"] = [
    "(((tfType1==1) && (tfType2==2)) || ((tfType1==2) && (tfType2==1)))", "diEOmtf"]
gmtCuts["diBEmtf"] = [
    "(((tfType1==0) && (tfType2==2)) || ((tfType1==2) && (tfType2==0)))", "diBEmtf"]


# TODO: Do these plots also for maxQual!

efficiencyList = []
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


# Plot efficiency for in- and ouput of uGMT w/o splitting into TF contributions

jpsi_efficiency_ntuples = []
jpsi_efficiency_ntuples.append(ugmt_dimu_file)
jpsi_efficiency_ntuples.append(ugmt_dimu_file)
ntuple_names = []
ntuple_names.append("tf_ntuple")
ntuple_names.append("ugmt_ntuple")
ugmt_inout_labels = []
ugmt_inout_labels.append(["Gen muons", "uGMT input muons", "uGMT"])
ugmt_inout_labels.append(["Gen muons", "uGMT output muons", "uGMT"])

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
cuts.extend(len(jpsi_efficiency_ntuples) * [gmtCuts["ugmt_diMu-pt1"]])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, jpsi_efficiency_ntuples, ntuple_names,
                                   ugmt_inout_labels, line_colours, cuts,
                                   "jPsi_efficiency",
                                   rootFolder=opts.outDir)

# Plot ghosting probability for in- and output of uGMT w/o splitting into TFs

ghostList = []
# NOTE: If no L1 muon at all deltaX will be 0! This will lead to inefficiencies
# TODO: Because of note above exclude events without an L1 muon.
# (we're not looking at efficiency, but ghost rate so we'd be underestimating
# ghost rate otherwise.)
# basically we want to know which percentage of singleMu events would be
# "poisoned" by ghosts
ghostList.append([["deltaEta_L1", "#Delta#eta(#mu#mu_{Ghost})"],
                       binningDict["distVeryWide"],
                       "abs(eta1-eta2)",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["deltaPhi_L1", "#Delta#phi(#mu#mu_{Ghost})"],
                       binningDict["distVeryWide"],
                       "abs(phi1-phi2)",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["deltaR_L1", "#DeltaR(#mu#mu_{Ghost})"],
                       binningDict["distVeryWide"],
                       "sqrt((eta1-eta2)**2+(phi1-phi2)**2)",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["deltaEta_L1-zoom", "#Delta#eta(#mu#mu_{Ghost})"],
                       binningDict["distNarrow"],
                       "abs(eta1-eta2)",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["deltaPhi_L1-zoom", "#Delta#phi(#mu#mu_{Ghost})"],
                       binningDict["distNarrow"],
                       "abs(phi1-phi2)",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["deltaR_L1-zoom", "#DeltaR(#mu#mu_{Ghost})"],
                       binningDict["distNarrow"],
                       "sqrt((eta1-eta2)**2+(phi1-phi2)**2)",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["mu1_L1Eta", "#eta(leading #mu_{L1})"],
                       binningDict["etaFineRestr"], "eta1",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["mu1_L1Phi", "#phi(leading #mu_{L1})"],
                       binningDict["phiFineRestr"], "phi1",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["mu1_L1Pt", "p_{T}(leading #mu_{L1}) [GeV/c]"],
                       binningDict["pt140Fine"], "pT1",
                       genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["mu1_genEta", "#eta(#mu)"],
                      binningDict["etaFineRestr"], "eta1_gen",
                      genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["mu1_genPhi", "#phi(#mu)"],
                      binningDict["phiFineRestr"], "phi1_gen",
                      genCuts["mu-pt1"], [0, 0.5]])
ghostList.append([["mu1_genPt", "p_{T}(#mu) [GeV/c]"],
                      binningDict["pt140Fine"], "pT1_gen",
                      genCuts["mu-pt1"], [0, 0.5]])

singleMu_ghosting_ntuples = []
singleMu_ghosting_ntuples.append(ugmt_singleMu_file)
singleMu_ghosting_ntuples.append(ugmt_singleMu_file)

for varList in ghostList:
    generateCombinedGhostPercHist(varList, singleMu_ghosting_ntuples,
                                  ntuple_names, ugmt_inout_labels,
                                  line_colours, cuts,
                                  "singleMu_ghosts", rootFolder=opts.outDir)

# Plotting resolution of track finder inputs.

resolution_check_dlabel = []
resolution_check_dlabel.append(["All muons", "BMTF muons", "uGMT"])
resolution_check_dlabel.append(["All muons", "BMTF muons, fine eta", "uGMT"])
resolution_check_dlabel.append(["All muons", "BMTF muons, coarse eta", "uGMT"])
resolution_check_dlabel.append(["All muons", "OMTF muons", "uGMT"])
resolution_check_dlabel.append(["All muons", "EMTF muons", "uGMT"])
resolution_check_ntuple = []
resolution_check_ntuple.extend(len(resolution_check_dlabel) * [ugmt_singleMu_file])
resolution_check_ntuple_name = []
resolution_check_ntuple_name.extend(
    len(resolution_check_dlabel) * ["tf_ntuple"])
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
                            binningDict["distSym"], "phi1-phi1_gen",
                            genCuts["mu-pt1"], [0, 1.4]])
resolutionCheckList.append([["etaResolution", "#Delta#eta(#mu_{L1}#mu_{Gen})"],
                            binningDict["distSym"], "eta1-eta1_gen",
                            genCuts["mu-pt1"], [0, 1.4]])
resolutionCheckList.append([["mu1_l1Eta", "#eta(#mu_{L1})"],
                            binningDict["etaFineRestr"], "eta1",
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

# Plot distance of a ghost from the L1 muon at uGMT input

ghost_distance_label = []
ghost_distance_label.append(["All muons", "BMTF muons", "uGMT"])
ghost_distance_label.append(["All muons", "OMTF muons", "uGMT"])
ghost_distance_label.append(["All muons", "EMTF muons", "uGMT"])
ghost_distance_label.append(["All muons", "BMTF+OMTF overlap muons", "uGMT"])
ghost_distance_label.append(["All muons", "EMTF+OMTF overlap muons", "uGMT"])
ghost_distance_ntuple = []
ghost_distance_ntuple.extend(len(ghost_distance_label) * [ugmt_singleMu_file])
ghost_distance_ntuple_name = []
ghost_distance_ntuple_name.extend(len(ghost_distance_label) * ["tf_ntuple"])
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
                          "phi1-phi2",
                          genCuts["mu-pt1"], [0, 1.4]])
ghostDistanceList.append([["etaResolution", "#Delta#eta(#mu_{L1}#mu_{Ghost})"],
                          binningDict["distSym"],
                          "eta1-eta2",
                          genCuts["mu-pt1"], [0, 1.4]])
ghostDistanceList.append([["ptResolution",
                           "#Delta p_{T}(#mu_{L1}#mu_{Ghost})"],
                          binningDict["pt140Fine"],
                          "pT1-pT2",
                          genCuts["mu-pt1"], [0, 1.4]])
for varList in ghostDistanceList:
    generateCombinedEfficiencyHist(varList, ghost_distance_ntuple,
                                   ghost_distance_ntuple_name,
                                   ghost_distance_label,
                                   ghost_distance_line_colour,
                                   ghost_distance_cuts, "ghost_distance",
                                   drawGenMus=False, drawStackPlot=True,
                                   rootFolder=opts.outDir)

# Plot efficiencies at uGMT ouput split by TF contributions.

tf_eff_labels = []
tf_eff_labels.append(
["Gen muons", "uGMT input", "TF", "uGMTinput"])
tf_eff_labels.append(
["Gen muons", "only BMTF muons", "uGMT", "diBMTF"])
tf_eff_labels.append(
["Gen muons", "only OMTF muons", "uGMT", "diOMTF"])
tf_eff_labels.append(
["Gen muons", "only EMTF muons", "uGMT", "diEMTF"])
tf_eff_labels.append(["Gen muons", "BMTF+OMTF muons", "uGMT",
"diBOMTF"])
tf_eff_labels.append(
["Gen muons", "OMTF+EMTF muons", "uGMT", "diOEMTF"])
tf_eff_labels.append(
["Gen muons", "BMTF+EMTF muons", "uGMT", "diBEMTF"])
tf_eff_ntuples = []
tf_eff_ntuples.extend(len(tf_eff_labels) * [ugmt_dimu_file])
tf_eff_ntuple_names = []
tf_eff_ntuple_names.append("tf_ntuple")
tf_eff_ntuple_names.extend((len(tf_eff_labels)-1) * ["ugmt_ntuple"])
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
tf_eff_cuts.append(gmtCuts["ugmt_diMu-pt1"])
tf_eff_cuts.append(gmtCuts["diBmtf"])
tf_eff_cuts.append(gmtCuts["diOmtf"])
tf_eff_cuts.append(gmtCuts["diEmtf"])
tf_eff_cuts.append(gmtCuts["diBOmtf"])
tf_eff_cuts.append(gmtCuts["diEOmtf"])
tf_eff_cuts.append(gmtCuts["diBEmtf"])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, tf_eff_ntuples,
                                   tf_eff_ntuple_names, tf_eff_labels,
                                   tf_eff_line_colours, tf_eff_cuts, "tf_eff",
                                   drawGenMus=True, drawStackPlot=True,
                                   rootFolder=opts.outDir)


tf_ghosts_ntuples = []
tf_ghosts_ntuples.extend(len(tf_eff_labels) * [ugmt_singleMu_file])

for varList in ghostList:
    generateCombinedGhostPercHist(varList, tf_ghosts_ntuples,
                                  tf_eff_ntuple_names, tf_eff_labels,
                                  tf_eff_line_colours, tf_eff_cuts,
                                  "tf_ghosts", drawGenMus=False,
                                  drawStackPlot=True, rootFolder=opts.outDir)
