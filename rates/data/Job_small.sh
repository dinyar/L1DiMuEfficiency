#!/bin/bash
cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_8_0_24/src/L1DiMuEfficiency/rates/data/
eval `scram runtime -sh`
root -l 'compareDiMuRates.C("file_list_280017-280385_baselineTuning_small", "file_list_280017-280385_conservativeTuning_small", "file_list_280017-280385_aggressiveTuning_small", "161227_small", "280017-280385", 11, 4, 2208)'
