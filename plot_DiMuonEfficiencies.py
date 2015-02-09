#!/usr/bin/python

from ROOT import *

def generateEfficiencyHist(varList):
    gStyle.SetOptStat(0)
    c1 = TCanvas('c1', "Efficiency vs. " + varList[0] + " - " + varList[5][1] + ", " + varList[6][1], 200, 10, 700, 500)
    tmpHist = TH1D("tmpHist", "", varList[1], varList[2], varList[3])
    efficiencyHist = TH1D("effHist", "Efficiency vs. " + varList[0] + " - " + varList[5][1] + ", " + varList[6][1], varList[1], varList[2], varList[3])
    ntuple.Project("tmpHist", varList[4], "1*" + varList[6][0])
    ntuple.Project("effHist", varList[4], "1*" + varList[5][0])
    efficiencyHist.Divide(tmpHist)
    efficiencyHist.DrawCopy()
    c1.Update()
    c1.Print("plots/hist_eff_"+varList[0]+"_"+varList[5][1] + "_" + varList[6][1]+".pdf", "pdf")

def generateRateHist(varList):
    c1 = TCanvas('c1', "Rate of " + varList[0] + " - " + varList[5][1], 200, 10, 700, 500)
    rateHist = TH1D("rateHist", "Rate of " + varList[0] + " - " + varList[5][1], varList[1], varList[2], varList[3])
    ntuple.Project("rateHist", varList[4], "1*" + varList[5][0])
    rateHist.DrawCopy()
    c1.Update()
    c1.Print("plots/hist_rate_"+varList[0]+"_"+varList[5][1]+".pdf", "pdf")

gROOT.Reset()

# etaScalePos= [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4];
# etaScale= [None] * 63
# for i in range(31):
#     etaScale[i]    = -etaScalePos[31-i]
#     etaScale[32+i] = etaScalePos[i+1]
# etaScale[31]=0;

cutDict = {}
cutDict["recoPt1"] = ["((pT1_reco>1) && (pT2_reco>1))", "DiRecoMu1"]
cutDict["gmtPt1"] = ["((pT1_GMT>1) && (pT2_GMT>1))", "DiGMTMu1"]
cutDict["recoPt5"] = ["((pT1_reco>5) && (pT2_reco>5))", "DiRecoMu5"]
cutDict["gmtPt5"] = ["((pT1_GMT>5) && (pT2_GMT>5))", "DiGMTMu5"]

efficiencyList = []
# TODO: Axis labels, think about more descriptive title.
# TODO: Can combine/infer stuff: Efficiency implies GMT/reco;
# Entries: Label for histogram (Will be used for filename and title) | binning | parameters used for project functions
efficiencyList.append(["deltaEta_reco", 63, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["gmtPt5"], cutDict["recoPt5"]])
efficiencyList.append(["deltaPhi_reco", 63, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["gmtPt5"], cutDict["recoPt5"]])
efficiencyList.append(["deltaR_reco", 63, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["gmtPt5"], cutDict["recoPt5"]])
efficiencyList.append(["pT1_reco", 25, 0, 50, "pT1_reco", cutDict["gmtPt1"], cutDict["recoPt1"]])
efficiencyList.append(["pT2_reco", 25, 0, 50, "pT2_reco", cutDict["gmtPt1"], cutDict["recoPt1"]])
# TODO: Add invariant calculation.
# TODO: Add 2-D hist showing efficiency for pT of both muons.

rateList = []
rateList.append(["deltaEta_reco", 63, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["recoPt5"]])
rateList.append(["deltaEta_reco", 63, 0, 2, "abs(Eta1_reco-Eta2_reco)", cutDict["gmtPt5"]])
rateList.append(["deltaPhi_reco", 63, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["recoPt5"]])
rateList.append(["deltaPhi_reco", 63, 0, 2, "abs(Phi1_reco-Phi2_reco)", cutDict["gmtPt5"]])
rateList.append(["deltaR_reco", 63, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["recoPt5"]])
rateList.append(["deltaR_reco", 63, 0, 2, "sqrt((Eta1_reco-Eta2_reco)**2+(Phi1_reco-Phi2_reco)**2)", cutDict["gmtPt5"]])
rateList.append(["pt1_reco", 25, 0, 50, "pT1_reco", cutDict["recoPt1"]])
rateList.append(["pt2_reco", 25, 0, 50, "pT2_reco", cutDict["recoPt1"]])
rateList.append(["pt1_reco", 25, 0, 50, "pT1_reco", cutDict["gmtPt1"]])
rateList.append(["pt2_reco", 25, 0, 50, "pT2_reco", cutDict["gmtPt1"]])

f = TFile.Open("DiMuNtuple.root")

ntuple = f.Get("ntuple")

for varList in efficiencyList:
    generateEfficiencyHist(varList)

for varList in rateList:
    generateRateHist(varList)
