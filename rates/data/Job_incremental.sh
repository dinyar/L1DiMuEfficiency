#!/bin/bash
cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_8_0_24/src/L1DiMuEfficiency/rates/data/
eval `scram runtime -sh`
root -l 'compareDiMuRates.C("fileList_fragments/file_list_279888-279991_noCOU_1, fileList_fragments/file_list_279588-280385_baselineTuning_7", "fileList_fragments/file_list_279588-280385_conservativeTuning_7", "fileList_fragments/file_list_279588-280385_aggressiveTuning_7", "161231_doubleMu-18-5", "279588-280385", 18, 5, 2208, false, true, true, false, 0, 15966299, 15035853, 14657510)'
