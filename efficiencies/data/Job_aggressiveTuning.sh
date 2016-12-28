#!/bin/bash
cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_8_0_9/src/L1DiMuEfficiency/efficiencies/data/
eval `scram runtime -sh`
root -l 'diMuEfficiency.C("fileList_singleMuData_Run2016G", "file_list_singleMu-tuning_v5", "file_list_jPsi-tuning_v5", "161227_aggressiveTuning", false, false, false, 11, 4)'
