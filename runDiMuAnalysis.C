{
    gROOT->ProcessLine(".x ../L1TriggerDPG/L1Ntuples/macros/initL1Analysis.C");
    gROOT->ProcessLine(".L GMTDiMuEfficiencyNtupleizer.C+");
    gROOT->ProcessLine("GMTDiMuEfficiencyNtupleizer macro = GMTDiMuEfficiencyNtupleizer(\"file_list\");");
    gROOT->ProcessLine("macro.run(-1);");
}
