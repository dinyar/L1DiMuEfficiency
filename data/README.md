== Scripts to run on data

= ghostDistances.C

Plots the distance between a L1 muon and its ghost for each track finder as well as for the boundaries between track finders. Distances plotted in dEta, dPhi, and dR.

Run with e.g. `.x ghostDistances.C("file_list", "[run number]")` in a Root interactive session.

Expects l1MuonRecoTree as well as either l1UpgradeTfMuonTree or l1UpgradeTfMuonEmuTree to be present.

= duMuRates.C

Creates plots vs. eta with unpacked and re-emulated data taken from the L1Upgrade tree. Plots created are for SingleMuOpen (but with a q > 4 cut), DoubleMuon0, and DoubleMuon with pT cuts that can be passed to the script (pT1 > 11, and pT2 > 4 by default). For the double muon plots both the leading and trailing muons are plotted.

Run with `.x diMuRates.C("file_list_275125", [mu1cut], [mu2cut])` in a Root interactive session. *Note:* The script expects mu1cut >= mu2cut.

Expects both L1UpgradeTree and L1UpgradeEmuTree to be present.

