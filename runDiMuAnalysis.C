{
    gROOT->ProcessLine(".x ../L1TriggerDPG/L1Ntuples/macros/initL1Analysis.C");
    gROOT->ProcessLine(".L GMTDiMuEfficiency.C+");
    gROOT->ProcessLine("GMTDiMuEfficiency macro = GMTDiMuEfficiency(\"file_list\");");
    gROOT->ProcessLine("macro.run(-1);");
}
