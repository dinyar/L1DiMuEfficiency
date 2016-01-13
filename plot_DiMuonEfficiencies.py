#!/usr/bin/python

from ROOT import *
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "../L1AnalysisHelpers"))
from CreateHistograms import *

gROOT.Reset()
gROOT.SetBatch(kTRUE)

# Cancel-out requirement for sector/wedge configuration
cancel_reqs = []
# Same track finder. Possibiliies:
# 1. Tivially neighbours
# 2. At wrap-around edge
# 2.a. Barrel
# 2.b. Overlap/endcap
cancel_reqs.append("((tfType1 == tfType2) &&\
                     ((abs(tfProcessor1-tfProcessor2) == 1) ||\
                      ((tfType1 == 0) &&\
                       (((tfProcessor1 == 0) && (tfProcessor2 == 11)) ||\
                        ((tfProcessor1 == 11) && (tfProcessor2 == 0)))) ||\
                      ((tfType1 > 0) &&\
                       (((tfProcessor1 == 0) && (tfProcessor2 == 11)) ||\
                        ((tfProcessor1 == 11) && (tfProcessor2 == 0))))\
                     ))")
# Different track finders, endcap/overlap. Possibilities:
# 1. Trivially neighbours
# 2. At wrap around edge
cancel_reqs.append("((((tfType1 == 1) && (tfType2 == 2)) ||\
                      ((tfType1 == 2) && (tfType2 == 1))) &&\
                     ((abs(tfProcessor1-tfProcessor2) < 2) ||\
                      (((tfProcessor1 == 0) && (tfProcessor2 == 5)) ||\
                       ((tfProcessor1 == 5) && (tfProcessor2 == 0)))\
                     ))")
# Different track finders, barrel/overlap. Possibilities:
# 1. Trivially beighbours
# 2. At wrap around edge
cancel_reqs.append("((((tfType1 == 0) && (tfType2 == 1)) &&\
                        (((2*tfProcessor2-1) <= tfProcessor1) &&\
                         ((2*tfProcessor2+2) >= tfProcessor1))) ||\
                      (((tfType1 == 1) && (tfType2 == 0)) &&\
                       (((2*tfProcessor1-1) <= tfProcessor2) &&\
                        ((2*tfProcessor1+2) >= tfProcessor2)))\
                    )")

cancel_requirement = ' || '.join(cancel_reqs)
cancel_requirement = '(' + cancel_requirement + ')'

# Cut dicts
genCuts = {}
genCuts["mu-pt1"] = ["(pT1_gen > 1)", "mu-ptGen1"]
genCuts["diMu-pt1"] = ["((pT1_gen > 1) && (pT2_gen > 1))", "diMu-ptGen1"]
gmtCuts = {}
gmtCuts["diMu-pt1"] = ["((pT1 > 1) && (pT2 > 1))", "diMu-pt1"]
# TODO: Try different distance cut for each track finder/area.
# TODO: Try different distance cuts for phi/eta
# TODO: Try using charge for matching.
gmtCuts["diMu-pt1_separatedFar"] = ["((pT1 > 1) && (pT2 > 1) && (sqrt((eta1-eta2)**2+(phi1-phi2)**2) > 0.1))", "diMu-pt1_separatedFar"]
gmtCuts["diMu-pt1_separatedFarFix"] = ["((pT1 > 1) && (pT2 > 1) && !((sqrt((eta1-eta2)**2+(phi1-phi2)**2) < 0.1) && " + cancel_requirement + "))", "diMu-pt1_separatedFarFix"]
gmtCuts["diMu-pt1_separatedNear"] = ["((pT1 > 1) && (pT2 > 1) && (sqrt((eta1-eta2)**2+(phi1-phi2)**2) > 0.01))", "diMu-pt1_separatedNear"]
gmtCuts["bmtf"] = ["(tfType1==0)", "bmtf"]
gmtCuts["omtf"] = ["(tfType1==1)", "omtf"]
gmtCuts["emtf"] = ["(tfType1==2)", "emtf"]
# TODO: Have cuts to only look at certain track finders

# TODO: Plot charge separately per TF.
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

# TODO: Get these via argparse
ntuple_files = []
ntuple_files.append("GMTDimuonNtuple.root")
ntuple_files.append("uGMTDimuonNtuple.root")
ntuple_files.append("uGMTDimuonNtuple.root")
ntuple_files.append("uGMTDimuonNtuple.root")
ntuple_files.append("uGMTDimuonNtuple.root")
ntuple_names = []
ntuple_names.append("gmt_ntuple")
ntuple_names.append("ugmt_ntuple")
ntuple_names.append("ugmt_ntuple")
ntuple_names.append("ugmt_ntuple")
ntuple_names.append("ugmt_ntuple")
distribution_labels = []
distribution_labels.append(["Gen muons", "GMT muons"])
distribution_labels.append(["Gen muons", "uGMT muons"])
distribution_labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.1"])
distribution_labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.1, w/ wedge comp"])
distribution_labels.append(["Gen muons", "uGMT muons w/ cancel-out #DeltaR<0.01"])
line_colours = []
line_colours.append(38)
line_colours.append(46)
line_colours.append(30)
line_colours.append(1)
line_colours.append(8)
gmt_cuts = []
gmt_cuts.append(gmtCuts["diMu-pt1"])
gmt_cuts.append(gmtCuts["diMu-pt1"])
gmt_cuts.append(gmtCuts["diMu-pt1_separatedFar"])
gmt_cuts.append(gmtCuts["diMu-pt1_separatedFarFix"])
gmt_cuts.append(gmtCuts["diMu-pt1_separatedNear"])

for varList in efficiencyList:
    generateCombinedEfficiencyHist(varList, ntuple_files, ntuple_names,
                                   distribution_labels, line_colours,
                                   gmt_cuts, "jPsi")


ccntuple = ["uGMTDimuonNtuple.root"]
ccntuple_name = ["ugmt_ntuple"]
ccdlabel = ["Gen muons", "uGMT muons"]
cclc = [1]
cccuts = []
cccuts.append(gmt_cuts["bmtf"])
cccuts.append(gmt_cuts["omtf"])
cccuts.append(gmt_cuts["emtf"])
chargeCheckList = []
chargeCheckList.append([["mu1_ch", "ch"],
                       binningDict["charge"], "ch1",
                       genCuts["diMu-pt1"], [0, 1.2]])
for varList in chargeCheckList:
    generateCombinedEfficiencyHist(varList, ccntuple, ccntuple_name,
                                   ccdlabel, cclc,
                                   cccuts, "charge_check")
