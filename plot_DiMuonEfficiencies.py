#!/usr/bin/python

from ROOT import *
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../L1AnalysisHelpers"))
from CreateHistograms import *

gROOT.Reset()
gROOT.SetBatch(kTRUE)

# Cut dicts
genCuts = {}
genCuts["mu-pt1"] = ["(pT1_gen > 1)", "mu-ptGen1"]
genCuts["L1GenMu-pt1"] = ["((pT1_gen > 1) && (pT1 > 1))", "mu-ptL1Gen1"]
genCuts["diMu-pt1"] = ["((pT1_gen > 1) && (pT2_gen > 1))", "diMu-ptGen1"]
gmtCuts = {}
gmtCuts["diMu-pt1"] = ["((pT1 > 1) && (pT2 > 1))", "diMu-pt1"]
gmtCuts["bmtf"] = ["(tfType1==0)", "bmtf"]
gmtCuts["omtf"] = ["(tfType1==1)", "omtf"]
gmtCuts["emtf"] = ["(tfType1==2)", "emtf"]

efficiencyList = []
# TODO: For mu1/mu2 plot single mu efficiencies?
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
efficiencyList.append([["deltaEta_gen", "#Delta#eta(#mu^{-}#mu^{+})"],
                       binningDict["distWide"],
                       "abs(eta1_gen-eta2_gen)",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["deltaPhi_gen", "#Delta#phi(#mu^{-}#mu^{+})"],
                       binningDict["distWide"],
                       "abs(phi1_gen-phi2_gen)",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["deltaR_gen", "#DeltaR(#mu^{-}#mu^{+})"],
                       binningDict["distWide"],
                       "sqrt((eta1_gen-eta2_gen)**2+(phi1_gen-phi2_gen)**2)",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["mu1_genEta", "#eta(leading #mu)"],
                       binningDict["etaFineRestr"], "eta1_gen",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["mu2_genEta", "#eta(trailing #mu)"],
                       binningDict["etaFineRestr"], "eta2_gen",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["mu1_genPhi", "#phi(leading #mu)"],
                       binningDict["phiFineRestr"], "phi1_gen",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["mu2_genPhi", "#phi(trailing #mu)"],
                       binningDict["phiFineRestr"], "phi2_gen",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["mu1_genPt", "p_{T}(leading #mu) [GeV/c]"],
                       binningDict["pt140Fine"], "pT1_gen",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["mu2_genPt", "p_{T}(trailing #mu) [GeV/c]"],
                       binningDict["pt140Fine"], "pT2_gen",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["jPsi_genEta", "#eta(J/#Psi)"],
                       binningDict["etaFineRestr"], "eta_jpsi",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["jPsi_genPhi", "#phi(J/#Psi)"],
                       binningDict["phiFineRestr"], "phi_jpsi",
                       genCuts["diMu-pt1"], [0, 1.2]])
efficiencyList.append([["jPsi_genPt", "p_{T}(J/#Psi) [GeV/c]"],
                       binningDict["pt140Fine"], "pT_jpsi",
                       genCuts["diMu-pt1"], [0, 1.2]])

jpsi_ntuples = []
jpsi_ntuples.append("GMTDimuonNtuple.root")
jpsi_ntuples.append("uGMTDimuonNtuple.root")
jpsi_ntuples.append("uGMTDimuonNtuple-dR0_3.root")
jpsi_ntuples.append("uGMTDimuonNtuple-dR0_1.root")
jpsi_ntuples.append("uGMTDimuonNtuple-dR0_05.root")
jpsi_ntuples.append("uGMTDimuonNtuple-dR0_01.root")
ntuple_names = []
ntuple_names.append("gmt_ntuple")
ntuple_names.append("ugmt_ntuple")
ntuple_names.append("ugmt_ntuple")
ntuple_names.append("ugmt_ntuple")
ntuple_names.append("ugmt_ntuple")
distribution_labels = []
distribution_labels.append(["Gen muons", "GMT muons", "GMT"])
distribution_labels.append(["Gen muons", "uGMT muons", "uGMT"])
distribution_labels.append(["Gen muons",
                            "uGMT muons w/ cancel-out #DeltaR<0.3", "uGMT"])
distribution_labels.append(["Gen muons",
                            "uGMT muons w/ cancel-out #DeltaR<0.1", "uGMT"])
distribution_labels.append(["Gen muons",
                            "uGMT muons w/ cancel-out #DeltaR<0.5", "uGMT"])
distribution_labels.append(["Gen muons",
                            "uGMT muons w/ cancel-out #DeltaR<0.01", "uGMT"])
line_colours = []
line_colours.append(38)
line_colours.append(46)
line_colours.append(30)
line_colours.append(1)
line_colours.append(8)
cuts = []
cuts.append(gmtCuts["diMu-pt1"])
cuts.append(gmtCuts["diMu-pt1"])
cuts.append(gmtCuts["diMu-pt1"])
cuts.append(gmtCuts["diMu-pt1"])
cuts.append(gmtCuts["diMu-pt1"])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, jpsi_ntuples, ntuple_names,
                                   distribution_labels, line_colours,
                                   cuts, "jPsi")


ccntuple = []
ccntuple.append("uGMTDimuonNtuple.root")
ccntuple.append("uGMTDimuonNtuple.root")
ccntuple.append("uGMTDimuonNtuple.root")
ccntuple_name = []
ccntuple_name.append("ugmt_ntuple")
ccntuple_name.append("ugmt_ntuple")
ccntuple_name.append("ugmt_ntuple")
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
                       binningDict["charge"], "ch1",
                       genCuts["diMu-pt1"], [0, 1.2]])
for varList in chargeCheckList:
    generateCombinedEfficiencyHist(varList, ccntuple, ccntuple_name,
                                   ccdlabel, cclc,
                                   cccuts, "charge_check")

ghostList = []
# TODO: If there is no L1 muon present deltaR etc are 0!
ghostList.append([["deltaEta_L1", "#Delta#eta(#mu#mu_{Ghost})"],
                  binningDict["distVeryWide"],
                  "abs(eta1-eta2)",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["deltaPhi_GMT", "#Delta#phi(#mu#mu_{Ghost})"],
                  binningDict["distVeryWide"],
                  "abs(phi1-phi2)",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["deltaR_L1", "#DeltaR(#mu#mu_{Ghost})"],
                  binningDict["distVeryWide"],
                  "sqrt((eta1-eta2)**2+(phi1-phi2)**2)",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["deltaEta_L1-zoom", "#Delta#eta(#mu#mu_{Ghost})"],
                  binningDict["distNarrow"],
                  "abs(eta1-eta2)",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["deltaPhi_L1-zoom", "#Delta#phi(#mu#mu_{Ghost})"],
                  binningDict["distNarrow"],
                  "abs(phi1-phi2)",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["deltaR_L1-zoom", "#DeltaR(#mu#mu_{Ghost})"],
                  binningDict["distNarrow"],
                  "sqrt((eta1-eta2)**2+(phi1-phi2)**2)",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["mu1_genEta", "#eta(#mu)"],
                  binningDict["etaFineRestr"], "eta1_gen",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["mu1_L1Eta", "#eta(leading #mu_{L1})"],
                  binningDict["etaFineRestr"], "eta1",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["mu1_genPhi", "#phi(#mu)"],
                  binningDict["phiFineRestr"], "phi1_gen",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["mu1_L1Phi", "#phi(leading #mu_{L1})"],
                  binningDict["phiFineRestr"], "phi1",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["mu1_genPt", "p_{T}(#mu) [GeV/c]"],
                  binningDict["pt140Fine"], "pT1_gen",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])
ghostList.append([["mu1_L1Pt", "p_{T}(leading #mu_{L1}) [GeV/c]"],
                  binningDict["pt140Fine"], "pT1",
                  genCuts["L1GenMu-pt1"], [0, 1.2]])

singleMu_ntuples = []
singleMu_ntuples.append("GMTSingleMuNtuple.root")
singleMu_ntuples.append("uGMTSingleMuNtuple.root")
singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_3.root")
singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_1.root")
singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_05.root")
singleMu_ntuples.append("uGMTSingleMuNtuple-dR0_01.root")
for varList in ghostList:
    generateCombinedEfficiencyHist(varList, singleMu_ntuples, ntuple_names,
                                   distribution_labels, line_colours,
                                   cuts, "singleMu")

resolution_check_ntuple = []
resolution_check_ntuple.append("uGMTSingleMuNtuple.root")
resolution_check_ntuple.append("uGMTSingleMuNtuple.root")
resolution_check_ntuple.append("uGMTSingleMuNtuple.root")
resolution_check_ntuple_name = []
resolution_check_ntuple_name.append("ugmt_ntuple")
resolution_check_ntuple_name.append("ugmt_ntuple")
resolution_check_ntuple_name.append("ugmt_ntuple")
resolution_check_dlabel = []
resolution_check_dlabel.append(["All muons", "BMTF muons", "uGMT"])
resolution_check_dlabel.append(["All muons", "OMTF muons", "uGMT"])
resolution_check_dlabel.append(["All muons", "EMTF muons", "uGMT"])
resolution_check_line_colour = []
resolution_check_line_colour.append(30)
resolution_check_line_colour.append(36)
resolution_check_line_colour.append(48)
resolution_check_cuts = []
resolution_check_cuts.append(gmtCuts["bmtf"])
resolution_check_cuts.append(gmtCuts["omtf"])
resolution_check_cuts.append(gmtCuts["emtf"])
resolutionCheckList = []
resolutionCheckList.append([["phiResolution", "#Delta#phi(#mu_{L1}#mu_{Gen})"],
                            binningDict["distSym"], "phi1-phi1_gen",
                            genCuts["L1GenMu-pt1"], [0, 1.2]])
resolutionCheckList.append([["etaResolution", "#Delta#eta(#mu_{L1}#mu_{Gen})"],
                            binningDict["distSym"], "eta1-eta1_gen",
                            genCuts["L1GenMu-pt1"], [0, 1.2]])
for varList in resolutionCheckList:
    generateCombinedEfficiencyHist(varList, resolution_check_ntuple,
                                   resolution_check_ntuple_name,
                                   resolution_check_dlabel,
                                   resolution_check_line_colour,
                                   resolution_check_cuts, "resolution_check")
