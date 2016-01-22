#!/usr/bin/env python
import sys
import os
from math import sqrt
from array import array
sys.path.append(os.path.join(os.path.dirname(__file__),
                "../L1TriggerDPG/L1Ntuples/macros/python"))
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

import ROOT as root
from L1Analysis import L1Ana, L1Ntuple

parser.add_argument("--NgenMu", dest="NgenMu", type=int, help="Number of generated muons to expect in input file.")
opts = parser.parse_args()

NgenMu = opts.NgenMu
if NgenMu == 1:
    ntupleName = "SingleMuNtuple"
elif NgenMu == 2:
    ntupleName = "DimuonNtuple"
cutList = []
cutList.append(["", None, None, None, None])
cutList.append(["-dR0_3", 5*[0.3], 5*[1], 5*[1], 5*[False]])
cutList.append(["-dR0_1", 5*[0.1], 5*[1], 5*[1], 5*[False]])
cutList.append(["-dR0_05", 5*[0.05], 5*[1], 5*[1], 5*[False]])
cutList.append(["-dR0_01", 5*[0.01], 5*[1], 5*[1], 5*[False]])
cutList.append(["-dR0_01-OMTF_dR0_1_chargeMatch",
                [0.01, 0.1, 0.01, 0.01, 0.01],
                5*[1], 5*[1], [False, True, False, False, False]])
cutList.append(["-dR0_01-OMTF_dR0_3_chargeMatch",
                [0.01, 0.3, 0.01, 0.01, 0.01],
                5*[1], 5*[1], [False, True, False, False, False]])


def checkMatchQuality(evt, mu1, mu2, dRcut, wEta, wPhi,
                      useChargeMatching=False, debug=False):
    if debug is True:
        print "evt: ", evt
        print "mu1: ", mu1
        print "mu1_eta: ", evt.ugmt.eta[mu1]
        print "mu2: ", mu2
        print "mu2_eta: ", evt.ugmt.eta[mu2]
        print "dRcut: ", dRcut
        print "wEta: ", wEta
        print "wPhi: ", wPhi
        print "useChargeMatching: ", useChargeMatching

    if useChargeMatching is True:
        if evt.ugmt.ch[mu1] != evt.ugmt.ch[mu2]:
            return False

    dEta = wEta * abs(evt.ugmt.eta[mu1]-evt.ugmt.eta[mu2])
    dPhi = wPhi * abs(evt.ugmt.phi[mu1]-evt.ugmt.phi[mu2])

    dR = sqrt(wEta * dEta**2 + wPhi * dPhi**2)

    if dR > dRcut:
        return False
    else:
        return True


def findCancelMus(evt, mu1, mu2):
    highPtMu = -1
    lowQualMu = -1

    if evt.ugmt.qual[mu1] <= evt.ugmt.qual[mu2]:
        lowQualMu = mu1
    else:
        lowQualMu = mu2

    if evt.ugmt.pt[mu1] > evt.ugmt.pt[mu2]:
        highPtMu = mu1
    else:
        highPtMu = mu2

    return lowQualMu, highPtMu


# TODO: At some point allow me to apply weights to dEta, dPhi depending on TF!
# TODO: Try using charge for matching.
# Order in TF-specific lists: BMTF, OMTF, EMTF, BMTF/OMTF, OMTF/EMTF
def doCancelOut(evt, dRcut, wEta=None, wPhi=None, useChargeMatching=None,
                debug=False):
    if debug is True:
        print "evt: ", evt
        print "dRcut: ", dRcut
        print "wEta: ", wEta
        print "wPhi: ", wPhi
        print "useChargeMatching: ", useChargeMatching

    if wEta is None:
        wEta = [1, 1, 1, 1, 1]
    if wPhi is None:
        wPhi = [1, 1, 1, 1, 1]
    if useChargeMatching is None:
        useChargeMatching = [False, False, False, False, False]

    ptCancelledMuons = set()
    qualCancelledMuons = set()

    if dRcut is None:
        return qualCancelledMuons, ptCancelledMuons

    for i in range(0, evt.ugmt.n):
        if evt.ugmt.bx[i] != 0:
            continue

        tfType1 = evt.ugmt.tfLink[i].tf
        tfIdx1 = evt.ugmt.tfLink[i].idx
        processor1 = evt.ugmt.tfInfo[tfType1].processor[tfIdx1]

        for j in range(i, evt.ugmt.n):
            if evt.ugmt.bx[i] != 0:
                continue

            tfType2 = evt.ugmt.tfLink[j].tf
            tfIdx2 = evt.ugmt.tfLink[j].idx
            processor2 = evt.ugmt.tfInfo[tfType2].processor[tfIdx2]

            match = -1
            # Same track finder. Possibilies for each one:
            # 1. Trivially neighbours
            # 2. At wrap-around edge
            if ((tfType1 == 0) and (tfType2 == 0)) and \
                ((processor1 == processor2+1) or
                 (processor2 == processor1+1) or
                 ((processor1 == 11) and (processor2 == 0)) or
                 ((processor2 == 11) and (processor1 == 0))):
                match = checkMatchQuality(evt, i, j, dRcut[0],
                                          wEta[0], wPhi[0],
                                          useChargeMatching[0], debug)
            elif ((tfType1 == 1) and (tfType2 == 1)) and \
                 ((processor1 == processor2+1) or
                  (processor2 == processor1+1) or
                  ((processor1 == 5) and (processor2 == 0)) or
                  ((processor2 == 5) and (processor1 == 0))):
                match = checkMatchQuality(evt, i, j, dRcut[1],
                                          wEta[1], wPhi[1],
                                          useChargeMatching[1], debug)
            elif ((tfType1 == 2) and (tfType2 == 2)) and \
                 ((processor1 == processor2+1) or
                  (processor2 == processor1+1) or
                  ((processor1 == 5) and (processor2 == 0)) or
                  ((processor2 == 5) and (processor1 == 0))):
                match = checkMatchQuality(evt, i, j, dRcut[2],
                                          wEta[2], wPhi[2],
                                          useChargeMatching[2], debug)
            # Different track finders, barrel/overlap. Possibilities:
            # 1. Trivially beighbours
            # 2. At wrap around edge
            elif ((tfType1 == 0) and (tfType2 == 1)) and \
                 ((processor1 == 2*processor2) or
                  (processor1 == 2*processor2+1) or
                  (processor1 == 2*processor2+2) or
                  (processor1 == 2*processor2+3) or
                  ((processor1 == 0) and (processor2 == 5)) or
                  ((processor1 == 1) and (processor2 == 5))):
                match = checkMatchQuality(evt, i, j, dRcut[3],
                                          wEta[3], wPhi[3],
                                          useChargeMatching[3], debug)
            elif((tfType1 == 1) and (tfType2 == 0)) and \
                ((processor2 == 2*processor1) or
                 (processor2 == 2*processor1+1) or
                 (processor2 == 2*processor1+2) or
                 (processor2 == 2*processor1+3) or
                 ((processor2 == 0) and (processor1 == 5)) or
                 ((processor2 == 1) and (processor1 == 5))):
                match = checkMatchQuality(evt, i, j, dRcut[3],
                                          wEta[3], wPhi[3],
                                          useChargeMatching[3], debug)
            # Different track finders, endcap/overlap. Possibilities:
            # 1. Trivially neighbours
            # 2. At wrap around edge
            elif ((tfType1 == 1) and (tfType2 == 2)) and \
                 ((processor1 == processor2+1) or
                  (processor1 == processor2-1) or
                  (processor1 == processor2) or
                  ((processor1 == 5) and (processor2 == 0)) or
                  ((processor2 == 5) and (processor1 == 0))):
                match = checkMatchQuality(evt, i, j, dRcut[4],
                                          wEta[4], wPhi[4],
                                          useChargeMatching[4], debug)
            elif ((tfType1 == 2) and (tfType2 == 1)) and \
                 ((processor1 == processor2+1) or
                  (processor1 == processor2-1) or
                  (processor1 == processor2) or
                  ((processor1 == 5) and (processor2 == 0)) or
                  ((processor2 == 5) and (processor1 == 0))):
                match = checkMatchQuality(evt, i, j, dRcut[4],
                                          wEta[4], wPhi[4],
                                          useChargeMatching[4], debug)

            if match is True:
                cancelledWqual, cancelledWpt = findCancelMus(evt, i, j)
                ptCancelledMuons.add(cancelledWpt)
                qualCancelledMuons.add(cancelledWqual)

    return qualCancelledMuons, ptCancelledMuons


def getGenMuons(evt):
    if evt.gen.pt[0] > evt.gen.pt[1]:
        return 0, 1
    else:
        return 1, 0


def getGmtMuons(evt):
    leadingPt = 0
    trailingPt = 0
    leadingGmtMu = -1
    trailingGmtMu = -1
    for i in range(0, evt.gmt.N):
        if evt.gmt.CandBx[i] != 0:
            continue

        if evt.gmt.Pt[i] > leadingPt:
            trailingPt = leadingPt
            trailingGmtMu = leadingGmtMu
            leadingPt = evt.gmt.Pt[i]
            leadingGmtMu = i
        elif evt.gmt.Pt[i] > trailingPt:
            trailingPt = evt.gmt.Pt[i]
            trailingGmtMu = i

    return leadingGmtMu, trailingGmtMu


# TODO: At some point possibly get TF muons here..
# (would also need to convert phi then)
def getUGmtMuons(evt, cancelledMuons):
    leadingPt = 0
    trailingPt = 0
    leadingUGmtMu = -1
    trailingUGmtMu = -1
    numMus = 0
    for i in range(0, evt.ugmt.n):
        if evt.ugmt.bx[i] != 0:
            continue
        elif i in cancelledMuons:
            continue

        numMus += 1
        if evt.ugmt.pt[i] > leadingPt:
            trailingPt = leadingPt
            trailingUGmtMu = leadingUGmtMu
            leadingPt = evt.ugmt.pt[i]
            leadingUGmtMu = i
        elif evt.ugmt.pt[i] > trailingPt:
            trailingPt = evt.ugmt.pt[i]
            trailingUGmtMu = i

    return leadingUGmtMu, trailingUGmtMu, numMus


def analyse(evt, gmt_content_list, ugmt_content_list,
            dRcut, wEta, wPhi, useChargeMatching, debug=False):
    count = 0
    for pdgId in evt.gen.id:
        if abs(pdgId) == 13:
            count += 1
    if count != NgenMu:
        print "Found {n} generated muons in event, skipping.".format(n=count)
        return [], []

    # Find muons with highest pT
    if NgenMu > 1:
        leadingMu, trailingMu = getGenMuons(evt)
    else:
        leadingMu = 0
        trailingMu = -1
    leadingGmtMu, trailingGmtMu = getGmtMuons(evt)
    highQualMuons, minPtMuons = doCancelOut(evt, dRcut, wEta, wPhi,
                                            useChargeMatching, debug)
    leadingUGmtMu_pt, trailingUGmtMu_pt, nL1Mus_pt = getUGmtMuons(evt,
                                                                  minPtMuons)
    leadingUGmtMu_q, trailingUGmtMu_q, nL1Mus_q = getUGmtMuons(evt,
                                                               highQualMuons)

    # Compute properties of J/Psi particle
    if NgenMu > 1:
        mu1 = root.TLorentzVector()
        mu2 = root.TLorentzVector()
        mu1.SetPxPyPzE(evt.gen.px[leadingMu], evt.gen.py[leadingMu],
                       evt.gen.pz[leadingMu], evt.gen.e[leadingMu])
        mu2.SetPxPyPzE(evt.gen.px[trailingMu], evt.gen.py[trailingMu],
                       evt.gen.pz[trailingMu], evt.gen.e[trailingMu])
        jPsi = root.TLorentzVector()
        jPsi = mu1 + mu2
    else:
        jPsi = root.TLorentzVector()

    gmt_content = array('f')
    for gmtVar in gmt_content_list:
        if gmtVar == "N":
            gmt_content.append(evt.gmt.N)
        elif gmtVar == "pT1" and (evt.gmt.N > 0):
            gmt_content.append(evt.gmt.Pt[leadingGmtMu])
        elif gmtVar == "eta1" and (evt.gmt.N > 0):
            gmt_content.append(evt.gmt.Eta[leadingGmtMu])
        elif gmtVar == "phi1" and (evt.gmt.N > 0):
            gmt_content.append(evt.gmt.Phi[leadingGmtMu])
        elif gmtVar == "qual1" and (evt.gmt.N > 0):
            gmt_content.append(evt.gmt.Qual[leadingGmtMu])
        elif gmtVar == "ch1" and (evt.gmt.N > 0):
            gmt_content.append(evt.gmt.Cha[leadingGmtMu])
        elif (gmtVar == "pT2") and (evt.gmt.N > 1):
            gmt_content.append(evt.gmt.Pt[trailingGmtMu])
        elif (gmtVar == "eta2") and (evt.gmt.N > 1):
            gmt_content.append(evt.gmt.Eta[trailingGmtMu])
        elif (gmtVar == "phi2") and (evt.gmt.N > 1):
            gmt_content.append(evt.gmt.Phi[trailingGmtMu])
        elif (gmtVar == "qual2") and (evt.gmt.N > 1):
            gmt_content.append(evt.gmt.Qual[trailingGmtMu])
        elif (gmtVar == "ch2") and (evt.gmt.N > 1):
            gmt_content.append(evt.gmt.Cha[trailingGmtMu])
        elif gmtVar == "pT1_gen":
            gmt_content.append(evt.gen.pt[leadingMu])
        elif gmtVar == "pT2_gen":
            gmt_content.append(evt.gen.pt[trailingMu])
        elif gmtVar == "eta1_gen":
            gmt_content.append(evt.gen.eta[leadingMu])
        elif gmtVar == "eta2_gen":
            gmt_content.append(evt.gen.eta[trailingMu])
        elif gmtVar == "phi1_gen":
            gmt_content.append(evt.gen.phi[leadingMu])
        elif gmtVar == "phi2_gen":
            gmt_content.append(evt.gen.phi[trailingMu])
        elif gmtVar == "pT_jpsi":
            gmt_content.append(jPsi.Pt())
        elif gmtVar == "eta_jpsi":
            gmt_content.append(jPsi.Eta())
        elif gmtVar == "phi_jpsi":
            gmt_content.append(jPsi.Phi())
        else:
            gmt_content.append(-11)

    ugmt_content = array('f')
    for ugmtVar in ugmt_content_list:
        # CAVEAT: This works only as long as uGMT doesn't perform cancel out
        if ugmtVar == "N":
            ugmt_content.append(evt.ugmt.n)
        elif ugmtVar == "pT1_pt" and (nL1Mus_pt > 0):
            ugmt_content.append(evt.ugmt.pt[leadingUGmtMu_pt])
        elif ugmtVar == "eta1_pt" and (nL1Mus_pt > 0):
            ugmt_content.append(evt.ugmt.eta[leadingUGmtMu_pt])
        elif ugmtVar == "phi1_pt" and (nL1Mus_pt > 0):
            ugmt_content.append(evt.ugmt.phi[leadingUGmtMu_pt])
        elif ugmtVar == "qual1_pt" and (nL1Mus_pt > 0):
            ugmt_content.append(evt.ugmt.qual[leadingUGmtMu_pt])
        elif ugmtVar == "ch1_pt" and (nL1Mus_pt > 0):
            ugmt_content.append(evt.ugmt.ch[leadingUGmtMu_pt])
        elif ugmtVar == "trkAddr1_pt" and (nL1Mus_pt > 0):
            tfType = evt.ugmt.tfLink[leadingUGmtMu_pt].tf
            tfIdx = evt.ugmt.tfLink[leadingUGmtMu_pt].idx
            trkAddr = evt.ugmt.tfInfo[tfType].trAddress[tfIdx]
            ugmt_content.append(trkAddr)
        elif ugmtVar == "tfType1_pt" and (nL1Mus_pt > 0):
            ugmt_content.append(evt.ugmt.tfLink[leadingUGmtMu_pt].tf)
        elif (ugmtVar == "tfProcessor1_pt") and (nL1Mus_pt > 0):
            tfType = evt.ugmt.tfLink[leadingUGmtMu_pt].tf
            tfIdx = evt.ugmt.tfLink[leadingUGmtMu_pt].idx
            processor = evt.ugmt.tfInfo[tfType].processor[tfIdx]
            ugmt_content.append(processor)
        elif (ugmtVar == "pT2_pt") and (nL1Mus_pt > 1):
            ugmt_content.append(evt.ugmt.pt[trailingUGmtMu_pt])
        elif (ugmtVar == "eta2_pt") and (nL1Mus_pt > 1):
            ugmt_content.append(evt.ugmt.eta[trailingUGmtMu_pt])
        elif (ugmtVar == "phi2_pt") and (nL1Mus_pt > 1):
            ugmt_content.append(evt.ugmt.phi[trailingUGmtMu_pt])
        elif (ugmtVar == "qual2_pt") and (nL1Mus_pt > 1):
            ugmt_content.append(evt.ugmt.qual[trailingUGmtMu_pt])
        elif (ugmtVar == "ch2_pt") and (nL1Mus_pt > 1):
            ugmt_content.append(evt.ugmt.ch[trailingUGmtMu_pt])
        elif (ugmtVar == "trkAddr2_pt") and (nL1Mus_pt > 1):
            tfType = evt.ugmt.tfLink[trailingUGmtMu_pt].tf
            tfIdx = evt.ugmt.tfLink[trailingUGmtMu_pt].idx
            trkAddr = evt.ugmt.tfInfo[tfType].trAddress[tfIdx]
            ugmt_content.append(trkAddr)
        elif (ugmtVar == "tfType2_pt") and (nL1Mus_pt > 1):
            ugmt_content.append(evt.ugmt.tfLink[trailingUGmtMu_pt].tf)
        elif (ugmtVar == "tfProcessor2_pt") and (nL1Mus_pt > 1):
            tfType = evt.ugmt.tfLink[trailingUGmtMu_pt].tf
            tfIdx = evt.ugmt.tfLink[trailingUGmtMu_pt].idx
            processor = evt.ugmt.tfInfo[tfType].processor[tfIdx]
            ugmt_content.append(processor)
        elif ugmtVar == "pT1_q" and (nL1Mus_q > 0):
            ugmt_content.append(evt.ugmt.pt[leadingUGmtMu_q])
        elif ugmtVar == "eta1_q" and (nL1Mus_q > 0):
            ugmt_content.append(evt.ugmt.eta[leadingUGmtMu_q])
        elif ugmtVar == "phi1_q" and (nL1Mus_q > 0):
            ugmt_content.append(evt.ugmt.phi[leadingUGmtMu_q])
        elif ugmtVar == "qual1_q" and (nL1Mus_q > 0):
            ugmt_content.append(evt.ugmt.qual[leadingUGmtMu_q])
        elif ugmtVar == "ch1_q" and (nL1Mus_q > 0):
            ugmt_content.append(evt.ugmt.ch[leadingUGmtMu_q])
        elif ugmtVar == "trkAddr1_q" and (nL1Mus_q > 0):
            tfType = evt.ugmt.tfLink[leadingUGmtMu_q].tf
            tfIdx = evt.ugmt.tfLink[leadingUGmtMu_q].idx
            trkAddr = evt.ugmt.tfInfo[tfType].trAddress[tfIdx]
            ugmt_content.append(trkAddr)
        elif ugmtVar == "tfType1_q" and (nL1Mus_q > 0):
            ugmt_content.append(evt.ugmt.tfLink[leadingUGmtMu_q].tf)
        elif (ugmtVar == "tfProcessor1_q") and (nL1Mus_q > 0):
            tfType = evt.ugmt.tfLink[leadingUGmtMu_q].tf
            tfIdx = evt.ugmt.tfLink[leadingUGmtMu_q].idx
            processor = evt.ugmt.tfInfo[tfType].processor[tfIdx]
            ugmt_content.append(processor)
        elif (ugmtVar == "pT2_q") and (nL1Mus_q > 1):
            ugmt_content.append(evt.ugmt.pt[trailingUGmtMu_q])
        elif (ugmtVar == "eta2_q") and (nL1Mus_q > 1):
            ugmt_content.append(evt.ugmt.eta[trailingUGmtMu_q])
        elif (ugmtVar == "phi2_q") and (nL1Mus_q > 1):
            ugmt_content.append(evt.ugmt.phi[trailingUGmtMu_q])
        elif (ugmtVar == "qual2_q") and (nL1Mus_q > 1):
            ugmt_content.append(evt.ugmt.qual[trailingUGmtMu_q])
        elif (ugmtVar == "ch2_q") and (nL1Mus_q > 1):
            ugmt_content.append(evt.ugmt.ch[trailingUGmtMu_q])
        elif (ugmtVar == "trkAddr2_q") and (nL1Mus_q > 1):
            tfType = evt.ugmt.tfLink[trailingUGmtMu_q].tf
            tfIdx = evt.ugmt.tfLink[trailingUGmtMu_q].idx
            trkAddr = evt.ugmt.tfInfo[tfType].trAddress[tfIdx]
            ugmt_content.append(trkAddr)
        elif (ugmtVar == "tfType2_q") and (nL1Mus_q > 1):
            ugmt_content.append(evt.ugmt.tfLink[trailingUGmtMu_q].tf)
        elif (ugmtVar == "tfProcessor2_q") and (nL1Mus_q > 1):
            tfType = evt.ugmt.tfLink[trailingUGmtMu_q].tf
            tfIdx = evt.ugmt.tfLink[trailingUGmtMu_q].idx
            processor = evt.ugmt.tfInfo[tfType].processor[tfIdx]
            ugmt_content.append(processor)
        elif ugmtVar == "pT1_gen":
            ugmt_content.append(evt.gen.pt[leadingMu])
        elif ugmtVar == "pT2_gen":
            ugmt_content.append(evt.gen.pt[trailingMu])
        elif ugmtVar == "eta1_gen":
            ugmt_content.append(evt.gen.eta[leadingMu])
        elif ugmtVar == "eta2_gen":
            ugmt_content.append(evt.gen.eta[trailingMu])
        elif ugmtVar == "phi1_gen":
            ugmt_content.append(evt.gen.phi[leadingMu])
        elif ugmtVar == "phi2_gen":
            ugmt_content.append(evt.gen.phi[trailingMu])
        elif ugmtVar == "ch1_gen":
            if evt.gen.id[leadingMu] > 0:
                # Muon
                ugmt_content.append(-1)
            else:
                # Anti muon
                ugmt_content.append(1)
        elif ugmtVar == "ch2_gen":
            if evt.gen.id[trailingMu] > 0:
                # Muon
                ugmt_content.append(-1)
            else:
                # Anti muon
                ugmt_content.append(1)
        elif ugmtVar == "pT_jpsi":
            ugmt_content.append(jPsi.Pt())
        elif ugmtVar == "eta_jpsi":
            ugmt_content.append(jPsi.Eta())
        elif ugmtVar == "phi_jpsi":
            ugmt_content.append(jPsi.Phi())
        else:
            ugmt_content.append(-11)

    return gmt_content, ugmt_content


def generate_content_lists():
    gmt = []
    gmt.append("N")
    gmt.append("pT1")
    gmt.append("pT2")
    gmt.append("eta1")
    gmt.append("eta2")
    gmt.append("phi1")
    gmt.append("phi2")
    gmt.append("qual1")
    gmt.append("qual2")
    gmt.append("ch1")
    gmt.append("ch2")
    gmt.append("pT1_gen")
    gmt.append("eta1_gen")
    gmt.append("phi1_gen")
    if NgenMu > 1:
        gmt.append("pT2_gen")
        gmt.append("eta2_gen")
        gmt.append("phi2_gen")
        gmt.append("pT_jpsi")
        gmt.append("eta_jpsi")
        gmt.append("phi_jpsi")
    ugmt = []
    ugmt.append("N")
    ugmt.append("pT1_pt")
    ugmt.append("pT2_pt")
    ugmt.append("eta1_pt")
    ugmt.append("eta2_pt")
    ugmt.append("phi1_pt")
    ugmt.append("phi2_pt")
    ugmt.append("qual1_pt")
    ugmt.append("qual2_pt")
    ugmt.append("ch1_pt")
    ugmt.append("ch2_pt")
    ugmt.append("trkAddr1_pt")
    ugmt.append("trkAddr2_pt")
    ugmt.append("tfType1_pt")
    ugmt.append("tfType2_pt")
    ugmt.append("tfProcessor1_pt")
    ugmt.append("tfProcessor2_pt")
    ugmt.append("pT1_q")
    ugmt.append("pT2_q")
    ugmt.append("eta1_q")
    ugmt.append("eta2_q")
    ugmt.append("phi1_q")
    ugmt.append("phi2_q")
    ugmt.append("qual1_q")
    ugmt.append("qual2_q")
    ugmt.append("ch1_q")
    ugmt.append("ch2_q")
    ugmt.append("trkAddr1_q")
    ugmt.append("trkAddr2_q")
    ugmt.append("tfType1_q")
    ugmt.append("tfType2_q")
    ugmt.append("tfProcessor1_q")
    ugmt.append("tfProcessor2_q")
    ugmt.append("pT1_gen")
    ugmt.append("eta1_gen")
    ugmt.append("phi1_gen")
    ugmt.append("ch1_gen")
    if NgenMu > 1:
        ugmt.append("pT2_gen")
        ugmt.append("eta2_gen")
        ugmt.append("phi2_gen")
        ugmt.append("ch2_gen")
        ugmt.append("pT_jpsi")
        ugmt.append("eta_jpsi")
        ugmt.append("phi_jpsi")

    return gmt, ugmt


def main():
    L1Ana.init_l1_analysis()
    print ""

    ntuple = L1Ntuple(opts.nevents)

    if opts.flist:
        ntuple.open_with_file_list(opts.flist)
    if opts.fname:
        ntuple.open_with_file(opts.fname)

    gmt_content_list, ugmt_content_list = generate_content_lists()
    gmt_content_string = ':'.join(gmt_content_list)
    ugmt_content_string = ':'.join(ugmt_content_list)

    gmt_ntuple_fname = "GMT" + ntupleName + ".root"
    gmt_f = root.TFile(gmt_ntuple_fname, 'recreate')
    gmt_f.cd()
    flat_gmt_tuple = root.TNtuple("gmt_ntuple", "ntupledump",
                                  gmt_content_string)

    ugmt_files = []
    ugmt_ntuples = []
    for cuts in cutList:
        ugmt_ntuple_fname = "uGMT" + ntupleName + cuts[0] + ".root"
        ugmt_f = root.TFile(ugmt_ntuple_fname, 'recreate')
        ugmt_f.cd()
        flat_ugmt_tuple = root.TNtuple("ugmt_ntuple", "ntupledump",
                                       ugmt_content_string)
        ugmt_files.append(ugmt_f)
        ugmt_ntuples.append(flat_ugmt_tuple)

    start_evt = opts.start_event
    end_evt = opts.start_event+ntuple.nevents
    gmt_ntuple_values = []
    ugmt_ntuple_values = []
    for i in range(start_evt, end_evt):
        event = ntuple[i]
        if (i+1) % 1000 == 0:
            L1Ana.log.info("Processing event: {n}".format(n=i))
            debug = True
        else:
            debug = False

        for cuts, ugmt_file, ugmt_tuple in zip(cutList, ugmt_files,
                                               ugmt_ntuples):
            gmt_ntuple_values, ugmt_ntuple_values = analyse(event,
                                                            gmt_content_list,
                                                            ugmt_content_list,
                                                            cuts[1], cuts[2],
                                                            cuts[3], cuts[4], debug)

            ugmt_file.cd()
            ugmt_tuple.Fill(ugmt_ntuple_values)
        gmt_f.cd()
        flat_gmt_tuple.Fill(gmt_ntuple_values)

    gmt_f.cd()
    gmt_f.Write()

    for ugmt_file in ugmt_files:
        ugmt_file.cd()
        ugmt_file.Write()

if __name__ == "__main__":
    main()
