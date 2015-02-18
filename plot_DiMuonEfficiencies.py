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
efficiencyList.append(["deltaEta_reco", 63, 0, 0.4, "abs(Eta1_reco-Eta2_reco)", cutDict["gmtPt5"], cutDict["recoPt5"], []])
efficiencyList.append(["deltaPhi_reco", 63, 0, 0.4, "abs(Phi1_reco-Phi2_reco)", cutDict["gmtPt5"], cutDict["recoPt5"], []])
efficiencyList.append(["deltaR_reco", 63, 0, 0.4, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["gmtPt5"], cutDict["recoPt5"], []])
efficiencyList.append(["mu1_recoEta", 25, -2.5, 2.5, "Eta1_reco", cutDict["gmtPt1"], cutDict["recoPt1"], []])
efficiencyList.append(["mu1_recoEta", 25, -2.5, 2.5, "Eta1_reco", cutDict["gmtPt1"], cutDict["recoPt1"], stackCutDict["subsystems_mu1"]])
efficiencyList.append(["mu2_recoEta", 25, -2.5, 2.5, "Eta2_reco", cutDict["gmtPt1"], cutDict["recoPt1"], []])
efficiencyList.append(["mu2_recoEta", 25, -2.5, 2.5, "Eta2_reco", cutDict["gmtPt1"], cutDict["recoPt1"], stackCutDict["subsystems_mu2"]])
efficiencyList.append(["mu1_recoPhi", 25, -3.2, 3.2, "Phi1_reco", cutDict["gmtPt1"], cutDict["recoPt1"], []])
efficiencyList.append(["mu1_recoPhi", 25, -3.2, 3.2, "Phi1_reco", cutDict["gmtPt1"], cutDict["recoPt1"], stackCutDict["subsystems_mu1"]])
efficiencyList.append(["mu2_recoPhi", 25, -3.2, 3.2, "Phi2_reco", cutDict["gmtPt1"], cutDict["recoPt1"], []])
efficiencyList.append(["mu2_recoPhi", 25, -3.2, 3.2, "Phi2_reco", cutDict["gmtPt1"], cutDict["recoPt1"], stackCutDict["subsystems_mu2"]])
efficiencyList.append(["mu1_recoPt", 25, 0, 50, "pT1_reco", cutDict["gmtPt1"], cutDict["recoPt1"], []])
efficiencyList.append(["mu1_recoPt", 25, 0, 50, "pT1_reco", cutDict["gmtPt1"], cutDict["recoPt1"], stackCutDict["subsystems_mu1"]])
efficiencyList.append(["mu2_recoPt", 25, 0, 50, "pT2_reco", cutDict["gmtPt1"], cutDict["recoPt1"], []])
efficiencyList.append(["mu2_recoPt", 25, 0, 50, "pT2_reco", cutDict["gmtPt1"], cutDict["recoPt1"], stackCutDict["subsystems_mu2"]])
# TODO: Add invariant mass calculation.
# TODO: Add 2-D hist showing efficiency for pT of both muons.

rateList = []
rateList.append(["deltaEta_reco", 63, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["recoPt5"]])
rateList.append(["deltaEta_reco", 63, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["gmtPt5"]])
rateList.append(["deltaPhi_reco", 63, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["recoPt5"]])
rateList.append(["deltaPhi_reco", 63, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["gmtPt5"]])
rateList.append(["deltaR_reco", 63, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["recoPt5"]])
rateList.append(["deltaR_reco", 63, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["gmtPt5"]])
rateList.append(["mu1_recoPt", 25, 0, 50, "pT1_reco", cutDict["recoPt1"]])
rateList.append(["mu2_recoPt", 25, 0, 50, "pT2_reco", cutDict["recoPt1"]])
rateList.append(["mu1_recoPt", 25, 0, 50, "pT1_reco", cutDict["gmtPt1"]])
rateList.append(["mu2_recoPt", 25, 0, 50, "pT2_reco", cutDict["gmtPt1"]])
rateList.append(["mu1_recoEta", 25, -2.5, 2.5, "Eta1_reco", cutDict["recoPt1"]])
rateList.append(["mu2_recoEta", 25, -2.5, 2.5, "Eta2_reco", cutDict["recoPt1"]])
rateList.append(["mu1_recoEta", 25, -2.5, 2.5, "Eta1_reco", cutDict["gmtPt1"]])
rateList.append(["mu2_recoEta", 25, -2.5, 2.5, "Eta2_reco", cutDict["gmtPt1"]])
rateList.append(["mu1_recoPhi", 25, -3.2, 3.2, "Phi1_reco", cutDict["recoPt1"]])
rateList.append(["mu2_recoPhi", 25, -3.2, 3.2, "Phi2_reco", cutDict["recoPt1"]])
rateList.append(["mu1_recoPhi", 25, -3.2, 3.2, "Phi1_reco", cutDict["gmtPt1"]])
rateList.append(["mu2_recoPhi", 25, -3.2, 3.2, "Phi2_reco", cutDict["gmtPt1"]])

f = TFile.Open("DiMuNtuple.root")

ntuple = f.Get("ntuple")

for varList in efficiencyList:
    generateEfficiencyHist(varList)

for varList in rateList:
    generateRateHist(varList)
