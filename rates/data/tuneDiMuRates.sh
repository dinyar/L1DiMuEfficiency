#!/bin/bash
#cd /afs/cern.ch/work/d/dinyar/cancel_out/mc/CMSSW_11_0_2/src/L1DiMuEfficiency/rates/data/
#cmsenv
root -l 'tuneDiMuRates.C("file_list_320838-321012_baselineTuning_short", "file_list_320838-321012_updatedTuning_short", "20201110_short", "320838-321012", 4, 4, 2544, false, false)'
