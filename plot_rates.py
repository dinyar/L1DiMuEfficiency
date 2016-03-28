#!/usr/bin/python

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),
                             "../L1AnalysisHelpers"))
from CreateHistograms import (binningDict, simplePlotter)

from ROOT import gROOT, kTRUE

gROOT.Reset()
gROOT.SetBatch(kTRUE)

cuts = {}
cuts["mu-ptGen1"] = ["(pT1_gen > 1)", "mu-ptGen1"]
cuts["diMu-ptGen1"] = ["((pT1_gen > 1) && (pT2_gen > 1))", "diMu-ptGen1"]
cuts["mu-ptGMT1"] = ["(pT1 > 1)", "diMu-ptGMT1"]
cuts["diMu-ptGMT1"] = ["((pT1 > 1) && (pT2 > 1))", "diMu-ptGMT1"]

rateList = []
rateList.append([["mu1_genPt", "p_{T}(leading mu) [GeV/c]"],
                 binningDict["pt140Fine"], "pT1_gen"])
rateList.append([["mu2_genPt", "p_{T}(trailing mu) [GeV/c]"],
                 binningDict["pt140Fine"], "pT2_gen"])

gmt_singleNu_file = "GMTmuonRate.root"
ugmt_singleNu_file = "uGMTmuonRate.root"

ntuples = []
ntuples.extend(4 * [gmt_singleNu_file])  # Gen single/di-mus, GMT single/di-mus
ntuples.extend(2 * [ugmt_singleNu_file])  # For uGMT single muons, di-muons

ntuple_names = []
ntuple_names.extend(4 * ["gmt_ntuple"])
ntuple_names.extend(2 * ["ugmt_ntuple"])

labels = []
labels.append("Gen single muon rate")
labels.append("Gen di-muon rate")
labels.append("GMT single muon rate")
labels.append("GMT di-muon rate")
labels.append("uGMT single muon rate")
labels.append("uGMT di-muon rate")

line_colours = []
line_colours.append(1)
line_colours.append(46)
line_colours.append(30)
line_colours.append(38)
line_colours.append(8)
line_colours.append(28)

cutStrings = []
cutStrings.append(cuts["mu-ptGen1"])
cutStrings.append(cuts["diMu-ptGen1"])
cutStrings.append(cuts["mu-ptGMT1"])
cutStrings.append(cuts["diMu-ptGMT1"])
cutStrings.append(cuts["mu-ptGMT1"])
cutStrings.append(cuts["diMu-ptGMT1"])

for varList in rateList:
    simplePlotter(varList, ntuples, ntuple_names,
                  labels, cutStrings, line_colours, "ratePlots",
                  drawLogY=True)
