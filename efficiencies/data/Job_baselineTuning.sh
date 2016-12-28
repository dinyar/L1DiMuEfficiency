#!/bin/bash
cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_8_0_9/src/L1DiMuEfficiency/efficiencies/data/
eval `scram runtime -sh`
root -l 'diMuEfficiency.C("fileList_singleMuData_Run2016G", "file_list_singleMu", "file_list_jPsi", "161227_defaultTuning", false, false, false, 11, 4)'
