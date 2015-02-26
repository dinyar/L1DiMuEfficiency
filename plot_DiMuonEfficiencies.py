#!/usr/bin/python

from ROOT import *
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../L1AnalysisHelpers"))
from CreateHistograms import *

gROOT.Reset()
gROOT.SetBatch(kTRUE);

efficiencyList = []
# TODO: Axis labels, think about more descriptive title.
# TODO: Can combine/infer stuff: Efficiency implies GMT/reco;
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
efficiencyList.append(["deltaEta_reco", 50, 0, 0.4, "abs(Eta1_reco-Eta2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5"]])
efficiencyList.append(["deltaEta_reco", 50, 0, 0.4, "abs(Eta1_reco-Eta2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-brl"]])
efficiencyList.append(["deltaEta_reco", 50, 0, 0.4, "abs(Eta1_reco-Eta2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-ovl"]])
efficiencyList.append(["deltaEta_reco", 50, 0, 0.4, "abs(Eta1_reco-Eta2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-fwd"]])
efficiencyList.append(["deltaPhi_reco", 50, 0, 0.4, "abs(Phi1_reco-Phi2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5"]])
efficiencyList.append(["deltaPhi_reco", 50, 0, 0.4, "abs(Phi1_reco-Phi2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-brl"]])
efficiencyList.append(["deltaPhi_reco", 50, 0, 0.4, "abs(Phi1_reco-Phi2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-ovl"]])
efficiencyList.append(["deltaPhi_reco", 50, 0, 0.4, "abs(Phi1_reco-Phi2_reco)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-fwd"]])
efficiencyList.append(["deltaR_reco", 50, 0, 0.4, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5"]])
efficiencyList.append(["deltaR_reco", 50, 0, 0.4, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-brl"]])
efficiencyList.append(["deltaR_reco", 50, 0, 0.4, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-ovl"]])
efficiencyList.append(["deltaR_reco", 50, 0, 0.4, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["diMu-gmtPt5"], cutDict["diMu-recoPt5-fwd"]])
efficiencyList.append(["mu1_recoEta", 100, -2.8, 2.8, "Eta1_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"]])
efficiencyList.append(["mu2_recoEta", 100, -2.8, 2.8, "Eta2_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"]])
efficiencyList.append(["mu1_recoPhi", 100, -3.2, 3.2, "Phi1_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"]])
efficiencyList.append(["mu2_recoPhi", 100, -3.2, 3.2, "Phi2_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"]])
efficiencyList.append(["mu1_recoPt", 100, 0, 100, "pT1_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"]])
efficiencyList.append(["mu2_recoPt", 100, 0, 100, "pT2_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"]])
# TODO: Add invariant mass calculation.
# TODO: Add 2-D hist showing efficiency for pT of both muons.

effStackList = []
effStackList.append(["mu1_recoEta", 100, -2.8, 2.8, "Eta1_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
effStackList.append(["mu2_recoEta", 100, -2.8, 2.8, "Eta2_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu2"]])
effStackList.append(["mu1_recoPhi", 100, -3.2, 3.2, "Phi1_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
effStackList.append(["mu2_recoPhi", 100, -3.2, 3.2, "Phi2_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu2"]])
effStackList.append(["mu1_recoPt", 100, 0, 100, "pT1_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
effStackList.append(["mu2_recoPt", 100, 0, 100, "pT2_reco", cutDict["diMu-gmtPt1"], cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu2"]])


rateList = []
rateList.append(["deltaEta_reco", 50, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["diMu-recoPt5"]])
rateList.append(["deltaEta_reco", 50, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["diMu-gmtPt5"]])
rateList.append(["deltaPhi_reco", 50, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["diMu-recoPt5"]])
rateList.append(["deltaPhi_reco", 50, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["diMu-gmtPt5"]])
rateList.append(["deltaR_reco", 50, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["diMu-recoPt5"]])
rateList.append(["deltaR_reco", 50, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["diMu-gmtPt5"]])
rateList.append(["mu1_recoPt", 100, 0, 100, "pT1_reco", cutDict["diMu-recoPt1"]])
rateList.append(["mu2_recoPt", 100, 0, 100, "pT2_reco", cutDict["diMu-recoPt1"]])
rateList.append(["mu1_recoPt", 100, 0, 100, "pT1_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu2_recoPt", 100, 0, 100, "pT2_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu1_recoEta", 100, -2.8, 2.8, "Eta1_reco", cutDict["diMu-recoPt1"]])
rateList.append(["mu2_recoEta", 100, -2.8, 2.8, "Eta2_reco", cutDict["diMu-recoPt1"]])
rateList.append(["mu1_recoEta", 100, -2.8, 2.8, "Eta1_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu2_recoEta", 100, -2.8, 2.8, "Eta2_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu1_recoPhi", 100, -3.2, 3.2, "Phi1_reco", cutDict["diMu-recoPt1"]])
rateList.append(["mu2_recoPhi", 100, -3.2, 3.2, "Phi2_reco", cutDict["diMu-recoPt1"]])
rateList.append(["mu1_recoPhi", 100, -3.2, 3.2, "Phi1_reco", cutDict["diMu-gmtPt1"]])
rateList.append(["mu2_recoPhi", 100, -3.2, 3.2, "Phi2_reco", cutDict["diMu-gmtPt1"]])

rateStackList = []
rateStackList.append(["mu1_recoPt", 100, 0, 100, "pT1_reco", cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
rateStackList.append(["mu2_recoPt", 100, 0, 100, "pT2_reco", cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
rateStackList.append(["mu1_recoEta", 100, -2.8, 2.8, "Eta1_reco", cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
rateStackList.append(["mu2_recoEta", 100, -2.8, 2.8, "Eta2_reco", cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
rateStackList.append(["mu1_recoPhi", 100, -3.2, 3.2, "Phi1_reco", cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])
rateStackList.append(["mu2_recoPhi", 100, -3.2, 3.2, "Phi2_reco", cutDict["diMu-recoPt1"], stackCutDict["subsystems_mu1"]])


f = TFile.Open("DiMuNtuple.root")

ntuple = f.Get("ntuple")

for varList in efficiencyList:
    generateEfficiencyHist(varList)

for varList in effStackList:
    generateEfficiencyStack(varList)

for varList in rateList:
    generateRateHist(varList)

for varList in rateStackList:
    generateRateStack(varList)
