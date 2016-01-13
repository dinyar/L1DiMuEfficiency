#!/usr/bin/env python
import sys
import os
from array import array
sys.path.append(os.path.join(os.path.dirname(__file__),
                "../L1TriggerDPG/L1Ntuples/macros/python"))
from ToolBox import parse_options_and_init_log
# have to do this first or ROOT masks the -h messages
opts, parser = parse_options_and_init_log()

import ROOT as root
from L1Analysis import L1Ana, L1Ntuple

# TODO: This should be passed in as command line option!
NgenMu = 1


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
def getUGmtMuons(evt):
    leadingPt = 0
    trailingPt = 0
    leadingUGmtMu = -1
    trailingUGmtMu = -1
    for i in range(0, evt.ugmt.n):
        if evt.ugmt.bx[i] != 0:
            continue

        if evt.ugmt.pt[i] > leadingPt:
            trailingPt = leadingPt
            trailingUGmtMu = leadingUGmtMu
            leadingPt = evt.ugmt.pt[i]
            leadingUGmtMu = i
        elif evt.ugmt.pt[i] > trailingPt:
            trailingPt = evt.ugmt.pt[i]
            trailingUGmtMu = i

    return leadingUGmtMu, trailingUGmtMu


def analyse(evt, gmt_content_list, ugmt_content_list):
    count = 0
    for pdgId in evt.gen.id:
        if abs(pdgId) == 13:
            count += 1
    if count != NgenMu:
        print "Found {n} generated muons in event, not processing event.".format(n=count)
        return [], []

    # Find muons with highest pT
    if NgenMu > 1:
        leadingMu, trailingMu = getGenMuons(evt)
    else:
        leadingMu = 0
        trailingMu = -1
    leadingGmtMu, trailingGmtMu = getGmtMuons(evt)
    leadingUGmtMu, trailingUGmtMu = getUGmtMuons(evt)

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
        jPsi = -1

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
        elif (gmtVar == "pT1_gen"):
            gmt_content.append(evt.gen.pt[leadingMu])
        elif (gmtVar == "pT2_gen"):
            gmt_content.append(evt.gen.pt[trailingMu])
        elif (gmtVar == "eta1_gen"):
            gmt_content.append(evt.gen.eta[leadingMu])
        elif (gmtVar == "eta2_gen"):
            gmt_content.append(evt.gen.eta[trailingMu])
        elif (gmtVar == "phi1_gen"):
            gmt_content.append(evt.gen.phi[leadingMu])
        elif (gmtVar == "phi2_gen"):
            gmt_content.append(evt.gen.phi[trailingMu])
        elif (gmtVar == "pT_jpsi"):
            gmt_content.append(jPsi.Pt())
        elif (gmtVar == "eta_jpsi"):
            gmt_content.append(jPsi.Eta())
        elif (gmtVar == "phi_jpsi"):
            gmt_content.append(jPsi.Phi())
        else:
            gmt_content.append(-11)

    ugmt_content = array('f')
    for ugmtVar in ugmt_content_list:
        # CAVEAT: This works only as long as uGMT doesn't perform cancel out
        if ugmtVar == "N":
            ugmt_content.append(evt.ugmt.n)
        elif ugmtVar == "pT1" and (evt.ugmt.n > 0):
            ugmt_content.append(evt.ugmt.pt[leadingUGmtMu])
        elif ugmtVar == "eta1" and (evt.ugmt.n > 0):
            ugmt_content.append(evt.ugmt.eta[leadingUGmtMu])
        elif ugmtVar == "phi1" and (evt.ugmt.n > 0):
            ugmt_content.append(evt.ugmt.phi[leadingUGmtMu])
        elif ugmtVar == "qual1" and (evt.ugmt.n > 0):
            ugmt_content.append(evt.ugmt.qual[leadingUGmtMu])
        elif ugmtVar == "ch1" and (evt.ugmt.n > 0):
            ugmt_content.append(evt.ugmt.ch[leadingUGmtMu])
        elif ugmtVar == "trkAddr1" and (evt.ugmt.n > 0):
            tfType = evt.ugmt.tfLink[leadingUGmtMu].tf
            tfIdx = evt.ugmt.tfLink[leadingUGmtMu].idx
            trkAddr = evt.ugmt.tfInfo[tfType].trAddress[tfIdx]
            ugmt_content.append(trkAddr)
        elif ugmtVar == "tfType1" and (evt.ugmt.n > 0):
            ugmt_content.append(evt.ugmt.tfLink[leadingUGmtMu].tf)
        elif (ugmtVar == "tfProcessor1") and (evt.ugmt.n > 0):
            tfType = evt.ugmt.tfLink[leadingUGmtMu].tf
            tfIdx = evt.ugmt.tfLink[leadingUGmtMu].idx
            processor = evt.ugmt.tfInfo[tfType].processor[tfIdx]
            ugmt_content.append(processor)
        elif (ugmtVar == "pT2") and (evt.ugmt.n > 1):
            ugmt_content.append(evt.ugmt.pt[trailingUGmtMu])
        elif (ugmtVar == "eta2") and (evt.ugmt.n > 1):
            ugmt_content.append(evt.ugmt.eta[trailingUGmtMu])
        elif (ugmtVar == "phi2") and (evt.ugmt.n > 1):
            ugmt_content.append(evt.ugmt.phi[trailingUGmtMu])
        elif (ugmtVar == "qual2") and (evt.ugmt.n > 1):
            ugmt_content.append(evt.ugmt.qual[trailingUGmtMu])
        elif (ugmtVar == "ch2") and (evt.ugmt.n > 1):
            ugmt_content.append(evt.ugmt.ch[trailingUGmtMu])
        elif (ugmtVar == "trkAddr2") and (evt.ugmt.n > 1):
            tfType = evt.ugmt.tfLink[trailingUGmtMu].tf
            tfIdx = evt.ugmt.tfLink[trailingUGmtMu].idx
            trkAddr = evt.ugmt.tfInfo[tfType].trAddress[tfIdx]
            ugmt_content.append(trkAddr)
        elif (ugmtVar == "tfType2") and (evt.ugmt.n > 1):
            ugmt_content.append(evt.ugmt.tfLink[trailingUGmtMu].tf)
        elif (ugmtVar == "tfProcessor2") and (evt.ugmt.n > 1):
            tfType = evt.ugmt.tfLink[trailingUGmtMu].tf
            tfIdx = evt.ugmt.tfLink[trailingUGmtMu].idx
            processor = evt.ugmt.tfInfo[tfType].processor[tfIdx]
            ugmt_content.append(processor)
        elif (ugmtVar == "pT1_gen"):
            ugmt_content.append(evt.gen.pt[leadingMu])
        elif (ugmtVar == "pT2_gen"):
            ugmt_content.append(evt.gen.pt[trailingMu])
        elif (ugmtVar == "eta1_gen"):
            ugmt_content.append(evt.gen.eta[leadingMu])
        elif (ugmtVar == "eta2_gen"):
            ugmt_content.append(evt.gen.eta[trailingMu])
        elif (ugmtVar == "phi1_gen"):
            ugmt_content.append(evt.gen.phi[leadingMu])
        elif (ugmtVar == "phi2_gen"):
            ugmt_content.append(evt.gen.phi[trailingMu])
        elif (ugmtVar == "ch1_gen"):
            if evt.gen.id[leadingMu] > 0:
                # Muon
                ugmt_content.append(-1)
            else:
                # Anti muon
                ugmt_content.append(1)
        elif (ugmtVar == "ch2_gen"):
            if evt.gen.id[trailingMu] > 0:
                # Muon
                ugmt_content.append(-1)
            else:
                # Anti muon
                ugmt_content.append(1)
        elif (ugmtVar == "pT_jpsi"):
            ugmt_content.append(jPsi.Pt())
        elif (ugmtVar == "eta_jpsi"):
            ugmt_content.append(jPsi.Eta())
        elif (ugmtVar == "phi_jpsi"):
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
    ugmt.append("pT1")
    ugmt.append("pT2")
    ugmt.append("eta1")
    ugmt.append("eta2")
    ugmt.append("phi1")
    ugmt.append("phi2")
    ugmt.append("qual1")
    ugmt.append("qual2")
    ugmt.append("ch1")
    ugmt.append("ch2")
    ugmt.append("trkAddr1")
    ugmt.append("trkAddr2")
    ugmt.append("tfType1")
    ugmt.append("tfType2")
    ugmt.append("tfProcessor1")
    ugmt.append("tfProcessor2")
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

    gmt_ntuple_fname = "GMTDimuonNtuple.root"
    gmt_f = root.TFile(gmt_ntuple_fname, 'recreate')
    gmt_f.cd()
    flat_gmt_tuple = root.TNtuple("gmt_ntuple", "ntupledump",
                                  gmt_content_string)

    ugmt_ntuple_fname = "uGMTDimuonNtuple.root"
    ugmt_f = root.TFile(ugmt_ntuple_fname, 'recreate')
    ugmt_f.cd()
    flat_ugmt_tuple = root.TNtuple("ugmt_ntuple", "ntupledump",
                                   ugmt_content_string)

    start_evt = opts.start_event
    end_evt = opts.start_event+ntuple.nevents
    gmt_ntuple_values = []
    ugmt_ntuple_values = []
    for i in range(start_evt, end_evt):
        event = ntuple[i]
        if (i+1) % 1000 == 0:
            L1Ana.log.info("Processing event: {n}".format(n=i))
        gmt_ntuple_values, ugmt_ntuple_values = analyse(event,
                                                        gmt_content_list,
                                                        ugmt_content_list)

        gmt_f.cd()
        flat_gmt_tuple.Fill(gmt_ntuple_values)
        ugmt_f.cd()
        flat_ugmt_tuple.Fill(ugmt_ntuple_values)

    gmt_f.cd()
    gmt_f.Write()
    ugmt_f.cd()
    ugmt_f.Write()

if __name__ == "__main__":
    main()
