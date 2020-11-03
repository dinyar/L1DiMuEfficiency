#!/bin/bash
cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_8_0_24/src/L1DiMuEfficiency/efficiencies/mc/edmFiles
eval `scram runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/d/dinyar/x509up_u57165
cmsRun l1NtupleMC_RAW2DIGI_GEN_Single.py
