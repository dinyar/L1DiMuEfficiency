#!/bin/bash
cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_8_0_24/src/L1DiMuEfficiency/rates/data/
eval `scram runtime -sh`
root -l 'compareDiMuRates.C("file_list_280017-280385_baselineTuning", "file_list_280017-280385_conservativeTuning", "file_list_280017-280385_aggressiveTuning", "161228_doubleMu-4-4", "280017-280385", 4, 4, 2208)'
