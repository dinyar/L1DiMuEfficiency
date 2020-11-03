#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TTree.h"

#include "CMS_lumi.C"
#include "tdrstyle.C"

#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

bool readFList(std::string fname, std::vector<std::string>& listNtuples,
               std::string l1Treepath, std::string recoTreepath);
int setupTChain(const std::vector<std::string> listNtuples, TChain* l1Chain,
                TChain* truthChain);
void getTfMuonRates(int nCollBunches, int nevents, TChain* tfChain,
                    TChain* recoChain, const int pT1cut, const int pT2cut,
                    TH1D& doubleMuGhostRateHist,
                    TH1D& doubleMuGhostRateVsPtHist,
                    TH1D& doubleMuGhostRateOpenHist,
                    TH1D& doubleMuGhostRateOpenVsPtHist, TH1D& doubleMuRateHist,
                    TH1D& doubleMuRateVsPtHist, TH1D& doubleMuRateOpenHist,
                    TH1D& doubleMuRateOpenVsPtHist, bool retrieve_hists,
                    bool update_hists, TString folder, TString identifier);
void getMuonRates(int nCollBunches, int nevents, TChain* l1Chain,
                  TChain* recoChain, const int pT1cut, const int pT2cut,
                  TH1D& doubleMuGhostRateHist, TH1D& doubleMuGhostRateVsPtHist,
                  TH1D& doubleMuGhostRateOpenHist,
                  TH1D& doubleMuGhostRateOpenVsPtHist, TH1D& doubleMuRateHist,
                  TH1D& doubleMuRateVsPtHist, TH1D& doubleMuRateOpenHist,
                  TH1D& doubleMuRateOpenVsPtHist, bool retrieve_hists,
                  bool update_hists, TString folder, TString identifier);
void calcRates(int nCollBunches, int nevents, TH1D& doubleMuGhostRateHist,
               TH1D& doubleMuGhostRateVsPtHist, TH1D& doubleMuGhostRateOpenHist,
               TH1D& doubleMuGhostRateOpenVsPtHist, TH1D& doubleMuRateHist,
               TH1D& doubleMuRateVsPtHist, TH1D& doubleMuRateOpenHist,
               TH1D& doubleMuRateOpenVsPtHist);
void drawHistograms(TH1D& noCouHist, TH1D& baselineHist, TH1D& conservativeHist,
                    TH1D& aggressiveHist, TString filename, TString xAxisLabel,
                    TString descString, TString plotFolder, TString run,
                    bool logy = false);

void compareDiMuRates(
    const char* file_list_no_cou, const char* file_list_baseline,
    const char* file_list_conservative, const char* file_list_aggressive,
    TString folder = "tmp", TString run = "XXX", int mu1cut = 11,
    int mu2cut = 4, int nCollBunches = 2028, bool retrieve_noCOU_hists = false,
    bool update_noCOU_hists = false, bool retrieve_hists = false,
    bool update_hists = false, int noCouEntries = 0, int baselineEntries = 0,
    int conservativeEntries = 0, int aggressiveEntries = 0) {
  TString plotFolder = "plots/" + run + "/" + folder + "/";
  mkdir("plots/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("plots/" + run, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir(plotFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  TString histFolder = "hists/" + run + "/" + folder + "/";
  mkdir("hists/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("hists/" + run, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir(histFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  TH1::AddDirectory(kFALSE);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  setTDRStyle();
  writeExtraText = true;  // Add "Preliminary"
  // lumi_13TeV  = "4.9 fb^{-1}";
  lumi_sqrtS = "13 TeV";  // used with iPeriod = 0, e.g. for simulation-only
                          // plots (default is an empty string)

  int iPeriod = 0;  // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses
                    // lumi_sqrtS)

  std::string l1Treepath("l1UpgradeEmuTree/L1UpgradeTree");
  std::string tfTreepath("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree");
  std::string recoTreepath("l1MuonRecoTree/Muon2RecoTree");

  TChain* chainL1NoCOU = new TChain(tfTreepath.c_str());
  TChain* chainRecoNoCOU = new TChain(recoTreepath.c_str());
  TChain* chainL1Baseline = new TChain(l1Treepath.c_str());
  TChain* chainRecoBaseline = new TChain(recoTreepath.c_str());
  TChain* chainL1Conservative = new TChain(l1Treepath.c_str());
  TChain* chainRecoConservative = new TChain(recoTreepath.c_str());
  TChain* chainL1Aggressive = new TChain(l1Treepath.c_str());
  TChain* chainRecoAggressive = new TChain(recoTreepath.c_str());

  if (!retrieve_hists) {
    std::vector<std::string> listNtuplesBaseline;
    std::vector<std::string> listNtuplesConservative;
    std::vector<std::string> listNtuplesAggressive;

    bool success = readFList(file_list_baseline, listNtuplesBaseline,
                             l1Treepath, recoTreepath);
    success &= readFList(file_list_conservative, listNtuplesConservative,
                         l1Treepath, recoTreepath);
    success &= readFList(file_list_aggressive, listNtuplesAggressive,
                         l1Treepath, recoTreepath);

    if (!success) {
      std::cout << "Couldn't read NTuple file list. Exiting.. " << std::endl;
      return;
    }

    baselineEntries =
        setupTChain(listNtuplesBaseline, chainL1Baseline, chainRecoBaseline);
    conservativeEntries = setupTChain(
        listNtuplesConservative, chainL1Conservative, chainRecoConservative);
    aggressiveEntries = setupTChain(listNtuplesAggressive, chainL1Aggressive,
                                    chainRecoAggressive);
  }

  if (!retrieve_noCOU_hists) {
    std::vector<std::string> listNtuplesNoCOU;

    bool success =
        readFList(file_list_no_cou, listNtuplesNoCOU, tfTreepath, recoTreepath);

    if (!success) {
      std::cout << "Couldn't read TF NTuple file list. Exiting.. " << std::endl;
      return;
    }

    noCouEntries = setupTChain(listNtuplesNoCOU, chainL1NoCOU, chainRecoNoCOU);
  }

  // mu bins
  int nMuBins = 25;
  float muLo = -2.5;
  float muHi = 2.5;
  float muBinWidth = (muHi - muLo) / nMuBins;
  int nMuBinsPt = 30;
  float ptBins[nMuBinsPt + 1];
  for (int i = 0; i < 21; ++i) {
    ptBins[i] = i;
  }
  for (int i = 0; i < 5; ++i) {
    ptBins[21 + i] = 21 + 2 * i;
  }
  ptBins[26] = 32;
  ptBins[27] = 36;
  ptBins[28] = 40;
  ptBins[29] = 45;
  ptBins[30] = 50;

  // make histos
  TH1D doubleMuGhostRatesNoCOU("doubleMuGhostRatesNoCOU", "", nMuBins,
                               muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesVsPtNoCOU("doubleMuGhostRatesVsPtNoCOU", "", nMuBinsPt,
                                   ptBins);
  TH1D doubleMuGhostRatesOpenNoCOU("doubleMuGhostRatesOpenNoCOU", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenVsPtNoCOU("doubleMuGhostRatesOpenVsPtNoCOU", "",
                                       nMuBinsPt, ptBins);
  TH1D doubleMuRatesNoCOU("doubleMuRatesNoCOU", "", nMuBins, muLo - 0.1,
                          muHi + 0.1);
  TH1D doubleMuRatesVsPtNoCOU("doubleMuRatesVsPtNoCOU", "", nMuBinsPt, ptBins);
  TH1D doubleMuRatesOpenNoCOU("doubleMuRatesOpenNoCOU", "", nMuBins, muLo - 0.1,
                              muHi + 0.1);
  TH1D doubleMuRatesOpenVsPtNoCOU("doubleMuRatesOpenVsPtNoCOU", "", nMuBinsPt,
                                  ptBins);

  TH1D doubleMuGhostRatesBaseline("doubleMuGhostRatesBaseline", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesVsPtBaseline("doubleMuGhostRatesVsPtBaseline", "",
                                      nMuBinsPt, ptBins);
  TH1D doubleMuGhostRatesOpenBaseline("doubleMuGhostRatesOpenBaseline", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenVsPtBaseline("doubleMuGhostRatesOpenVsPtBaseline",
                                          "", nMuBinsPt, ptBins);
  TH1D doubleMuRatesBaseline("doubleMuRatesBaseline", "", nMuBins, muLo - 0.1,
                             muHi + 0.1);
  TH1D doubleMuRatesVsPtBaseline("doubleMuRatesVsPtBaseline", "", nMuBinsPt,
                                 ptBins);
  TH1D doubleMuRatesOpenBaseline("doubleMuRatesOpenBaseline", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenVsPtBaseline("doubleMuRatesOpenVsPtBaseline", "",
                                     nMuBinsPt, ptBins);

  TH1D doubleMuGhostRatesConservative("doubleMuGhostRatesConservative", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesVsPtConservative("doubleMuGhostRatesVsPtConservative",
                                          "", nMuBinsPt, ptBins);
  TH1D doubleMuGhostRatesOpenConservative("doubleMuGhostRatesOpenConservative",
                                          "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenVsPtConservative(
      "doubleMuGhostRatesOpenVsPtConservative", "", nMuBinsPt, ptBins);
  TH1D doubleMuRatesConservative("doubleMuRatesConservative", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesVsPtConservative("doubleMuRatesVsPtConservative", "",
                                     nMuBinsPt, ptBins);
  TH1D doubleMuRatesOpenConservative("doubleMuRatesOpenConservative", "",
                                     nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenVsPtConservative("doubleMuRatesOpenVsPtConservative",
                                         "", nMuBinsPt, ptBins);

  TH1D doubleMuGhostRatesAggressive("doubleMuGhostRatesAggressive", "", nMuBins,
                                    muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesVsPtAggressive("doubleMuGhostRatesVsPtAggressive", "",
                                        nMuBinsPt, ptBins);
  TH1D doubleMuGhostRatesOpenAggressive("doubleMuGhostRatesOpenAggressive", "",
                                        nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenVsPtAggressive(
      "doubleMuGhostRatesOpenVsPtAggressive", "", nMuBinsPt, ptBins);
  TH1D doubleMuRatesAggressive("doubleMuRatesAggressive", "", nMuBins,
                               muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesVsPtAggressive("doubleMuRatesVsPtAggressive", "", nMuBinsPt,
                                   ptBins);
  TH1D doubleMuRatesOpenAggressive("doubleMuRatesOpenAggressive", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenVsPtAggressive("doubleMuRatesOpenVsPtAggressive", "",
                                       nMuBinsPt, ptBins);

  getTfMuonRates(nCollBunches, noCouEntries, chainL1NoCOU, chainRecoNoCOU,
                 mu1cut, mu2cut, doubleMuGhostRatesNoCOU,
                 doubleMuGhostRatesVsPtNoCOU, doubleMuGhostRatesOpenNoCOU,
                 doubleMuGhostRatesOpenVsPtNoCOU, doubleMuRatesNoCOU,
                 doubleMuRatesVsPtNoCOU, doubleMuRatesOpenNoCOU,
                 doubleMuRatesOpenVsPtNoCOU, retrieve_noCOU_hists,
                 update_noCOU_hists, histFolder, "NoCOU");

  getMuonRates(nCollBunches, baselineEntries, chainL1Baseline,
               chainRecoBaseline, mu1cut, mu2cut, doubleMuGhostRatesBaseline,
               doubleMuGhostRatesVsPtBaseline, doubleMuGhostRatesOpenBaseline,
               doubleMuGhostRatesOpenVsPtBaseline, doubleMuRatesBaseline,
               doubleMuRatesVsPtBaseline, doubleMuRatesOpenBaseline,
               doubleMuRatesOpenVsPtBaseline, retrieve_hists, update_hists,
               histFolder, "Baseline");

  getMuonRates(
      nCollBunches, conservativeEntries, chainL1Conservative,
      chainRecoConservative, mu1cut, mu2cut, doubleMuGhostRatesConservative,
      doubleMuGhostRatesVsPtConservative, doubleMuGhostRatesOpenConservative,
      doubleMuGhostRatesOpenVsPtConservative, doubleMuRatesConservative,
      doubleMuRatesVsPtConservative, doubleMuRatesOpenConservative,
      doubleMuRatesOpenVsPtConservative, retrieve_hists, update_hists,
      histFolder, "Conservative");

  getMuonRates(nCollBunches, aggressiveEntries, chainL1Aggressive,
               chainRecoAggressive, mu1cut, mu2cut,
               doubleMuGhostRatesAggressive, doubleMuGhostRatesVsPtAggressive,
               doubleMuGhostRatesOpenAggressive,
               doubleMuGhostRatesOpenVsPtAggressive, doubleMuRatesAggressive,
               doubleMuRatesVsPtAggressive, doubleMuRatesOpenAggressive,
               doubleMuRatesOpenVsPtAggressive, retrieve_hists, update_hists,
               histFolder, "Aggressive");

  calcRates(nCollBunches, noCouEntries, doubleMuGhostRatesNoCOU,
            doubleMuGhostRatesVsPtNoCOU, doubleMuGhostRatesOpenNoCOU,
            doubleMuGhostRatesOpenVsPtNoCOU, doubleMuRatesNoCOU,
            doubleMuRatesVsPtNoCOU, doubleMuRatesOpenNoCOU,
            doubleMuRatesOpenVsPtNoCOU);

  calcRates(nCollBunches, baselineEntries, doubleMuGhostRatesBaseline,
            doubleMuGhostRatesVsPtBaseline, doubleMuGhostRatesOpenBaseline,
            doubleMuGhostRatesOpenVsPtBaseline, doubleMuRatesBaseline,
            doubleMuRatesVsPtBaseline, doubleMuRatesOpenBaseline,
            doubleMuRatesOpenVsPtBaseline);

  calcRates(nCollBunches, conservativeEntries, doubleMuGhostRatesConservative,
            doubleMuGhostRatesVsPtConservative,
            doubleMuGhostRatesOpenConservative,
            doubleMuGhostRatesOpenVsPtConservative, doubleMuRatesConservative,
            doubleMuRatesVsPtConservative, doubleMuRatesOpenConservative,
            doubleMuRatesOpenVsPtConservative);

  calcRates(nCollBunches, aggressiveEntries, doubleMuGhostRatesAggressive,
            doubleMuGhostRatesVsPtAggressive, doubleMuGhostRatesOpenAggressive,
            doubleMuGhostRatesOpenVsPtAggressive, doubleMuRatesAggressive,
            doubleMuRatesVsPtAggressive, doubleMuRatesOpenAggressive,
            doubleMuRatesOpenVsPtAggressive);

  std::ostringstream oss;
  oss << "ZeroBias, L1_DoubleMu_" << mu1cut << "_" << mu2cut;

  // doubleMuGhostRatesNoCOU.SetMaximum(0.043);
  // doubleMuGhostRatesBaseline.SetMaximum(0.043);
  // doubleMuGhostRatesConservative.SetMaximum(0.043);
  // doubleMuGhostRatesAggressive.SetMaximum(0.043);

  drawHistograms(doubleMuRatesNoCOU, doubleMuRatesBaseline,
                 doubleMuRatesConservative, doubleMuRatesAggressive,
                 "doubleMuRatesDoubleMuonLeading", "#eta (leading #mu)",
                 oss.str().c_str(), plotFolder, run);
  drawHistograms(doubleMuRatesVsPtNoCOU, doubleMuRatesVsPtBaseline,
                 doubleMuRatesVsPtConservative, doubleMuRatesVsPtAggressive,
                 "doubleMuRatesDoubleMuonVsPt", "p_{T} (leading #mu)",
                 oss.str().c_str(), plotFolder, run, true);
  drawHistograms(doubleMuRatesOpenNoCOU, doubleMuRatesOpenBaseline,
                 doubleMuRatesOpenConservative, doubleMuRatesOpenAggressive,
                 "doubleMuOpenRatesDoubleMuonLeading", "#eta (leading #mu)",
                 "Zero Bias, L1_DoubleMu0", plotFolder, run);
  drawHistograms(doubleMuRatesOpenVsPtNoCOU, doubleMuRatesOpenVsPtBaseline,
                 doubleMuRatesOpenVsPtConservative,
                 doubleMuRatesOpenVsPtAggressive,
                 "doubleMuOpenRatesDoubleMuonVsPt", "p_{T} (leading #mu)",
                 "Zero Bias, L1_DoubleMu0", plotFolder, run, true);

  oss << ", ghost rate";

  drawHistograms(doubleMuGhostRatesNoCOU, doubleMuGhostRatesBaseline,
                 doubleMuGhostRatesConservative, doubleMuGhostRatesAggressive,
                 "ghostRatesDoubleMuonLeading", "#eta (leading #mu)",
                 oss.str().c_str(), plotFolder, run);
  drawHistograms(doubleMuGhostRatesVsPtNoCOU, doubleMuGhostRatesVsPtBaseline,
                 doubleMuGhostRatesVsPtConservative,
                 doubleMuGhostRatesVsPtAggressive, "ghostRatesDoubleMuonVsPt",
                 "p_{T} (leading #mu)", oss.str().c_str(), plotFolder, run,
                 true);
  drawHistograms(doubleMuGhostRatesOpenNoCOU, doubleMuGhostRatesOpenBaseline,
                 doubleMuGhostRatesOpenConservative,
                 doubleMuGhostRatesOpenAggressive,
                 "ghostRatesDoubleMuonOpenLeading", "#eta (leading #mu)",
                 "Zero Bias, L1_DoubleMu0, ghost rate", plotFolder, run);
  drawHistograms(doubleMuGhostRatesOpenVsPtNoCOU,
                 doubleMuGhostRatesOpenVsPtBaseline,
                 doubleMuGhostRatesOpenVsPtConservative,
                 doubleMuGhostRatesOpenVsPtAggressive,
                 "ghostRatesDoubleMuonOpenVsPt", "p_{T} (leading #mu)",
                 "Zero Bias, L1_DoubleMu0, ghost rate", plotFolder, run, true);
}

bool readFList(std::string fname, std::vector<std::string>& listNtuples,
               std::string l1Treepath, std::string recoTreepath) {
  // OpenNtupleList
  std::ifstream flist(fname);
  if (!flist) {
    std::cout << "File " << fname << " is not found!" << std::endl;
    return false;
  }

  while (!flist.eof()) {
    std::string str;
    getline(flist, str);
    if (!flist.fail()) {
      if (str != "") listNtuples.push_back(str);
    }
  }

  // CheckFirstFile
  if (listNtuples.size() == 0) {
    return false;
  }

  TFile* rf = TFile::Open(listNtuples.at(0).c_str());

  if (rf == 0) {
    return false;
  }
  if (rf->IsOpen() == 0) {
    return false;
  }

  TTree* treeL1 = (TTree*)rf->Get(l1Treepath.c_str());
  TTree* treeReco = (TTree*)rf->Get(recoTreepath.c_str());

  if (!treeL1) {
    std::cout << l1Treepath << " tree not found.. " << std::endl;
    return false;
  }
  if (!treeReco) {
    std::cout << "No reco tree found.. " << std::endl;
    return false;
  }

  return true;
}

int setupTChain(const std::vector<std::string> listNtuples, TChain* l1Chain,
                TChain* truthChain) {
  for (unsigned int i = 0; i < listNtuples.size(); i++) {
    std::cout << " -- Adding " << listNtuples[i] << std::endl;
    l1Chain->Add(listNtuples[i].c_str());
    truthChain->Add(listNtuples[i].c_str());
  }

  // Init
  std::cout << "Estimate the number of entries... ";
  int nentries = l1Chain->GetEntries();
  std::cout << nentries << std::endl;

  return nentries;
}

bool matchL1toReco(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                   double mu1Eta, double mu2Eta, double mu1Phi, double mu2Phi) {
  double matchingWindowEta{0.5};
  double matchingWindowPhi{3.14159};
  bool matchedMu1{false};
  bool matchedMu2{false};
  for (int i = 0; i < reco_->nMuons; ++i) {
    if (reco_->isMediumMuon[i]) {
      if (!matchedMu1 && std::abs(reco_->phi[i] - mu1Phi) < matchingWindowPhi &&
          std::abs(reco_->eta[i] - mu1Eta) < matchingWindowEta) {
        matchedMu1 = true;
      } else if (!matchedMu2 &&
                 std::abs(reco_->phi[i] - mu2Phi) < matchingWindowPhi &&
                 std::abs(reco_->eta[i] - mu2Eta) < matchingWindowEta) {
        matchedMu2 = true;
      }
    }
  }
  if (matchedMu1 && matchedMu2) {
    return true;
  } else {
    return false;
  }
}

bool bestMatchL1toReco(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                       double mu1Eta, double mu2Eta, double mu1Phi,
                       double mu2Phi) {
  if (matchL1toReco(reco_, mu1Eta, mu2Eta, mu1Phi, mu2Phi) ||
      matchL1toReco(reco_, mu2Eta, mu1Eta, mu2Phi, mu1Phi)) {
    return true;
  } else {
    return false;
  }
}

void getTfMuonRates(int nCollBunches, int nevents, TChain* tfChain,
                    TChain* recoChain, const int mu1cut, const int mu2cut,
                    TH1D& doubleMuGhostRateHist,
                    TH1D& doubleMuGhostRateVsPtHist,
                    TH1D& doubleMuGhostRateOpenHist,
                    TH1D& doubleMuGhostRateOpenVsPtHist, TH1D& doubleMuRateHist,
                    TH1D& doubleMuRateVsPtHist, TH1D& doubleMuRateOpenHist,
                    TH1D& doubleMuRateOpenVsPtHist, bool retrieve_hists,
                    bool update_hists, TString folder, TString identifier) {
  TFile* f;
  if (retrieve_hists || update_hists) {
    std::cout << "Retrieving histograms.. " << std::endl;
    f = TFile::Open(folder + identifier + "Hists.root", "read");
    std::cout << "Opening file: " << folder + identifier + "Hists.root"
              << std::endl;
    doubleMuGhostRateHist =
        *(static_cast<TH1D*>(f->Get("doubleMuGhostRates" + identifier)));
    // doubleMuGhostRateHistTest->SetDirectory(0);
    doubleMuGhostRateVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuGhostRatesVsPt" + identifier)));
    // doubleMuGhostRateVsPtHist->SetDirectory(0);
    doubleMuGhostRateOpenHist =
        *(static_cast<TH1D*>(f->Get("doubleMuGhostRatesOpen" + identifier)));
    //->SetDirectory(0);
    doubleMuGhostRateOpenVsPtHist = *(
        static_cast<TH1D*>(f->Get("doubleMuGhostRatesOpenVsPt" + identifier)));
    //->SetDirectory(0);
    doubleMuRateHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRates" + identifier)));
    //->SetDirectory(0);
    doubleMuRateVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesVsPt" + identifier)));
    //->SetDirectory(0);
    doubleMuRateOpenHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesOpen" + identifier)));
    //->SetDirectory(0);
    doubleMuRateOpenVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesOpenVsPt" + identifier)));
    //->SetDirectory(0);
    f->Close();
    if (retrieve_hists) {
      return;
    }
  }

  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* bmtf_ =
      new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* omtf_ =
      new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* emtf_ =
      new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_ =
      new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  tfChain->SetBranchAddress("L1UpgradeBmtfMuon", &bmtf_);
  tfChain->SetBranchAddress("L1UpgradeBmtfMuon", &omtf_);
  tfChain->SetBranchAddress("L1UpgradeEmtfMuon", &emtf_);
  recoChain->SetBranchAddress("Muon", &reco_);

  for (Long64_t jentry = 0; jentry < nevents; ++jentry) {
    if ((jentry % 1000) == 0)
      std::cout << "Done " << jentry << " TF events.. \n";

    // std::cout << "Getting TF entry.." << std::endl;
    tfChain->GetEntry(jentry);
    int nTfMuons{bmtf_->nTfMuons};
    nTfMuons += omtf_->nTfMuons;
    nTfMuons += emtf_->nTfMuons;
    if (nTfMuons < 2) {
      continue;
    }

    // std::cout << "Getting reco entry.." << std::endl;
    recoChain->GetEntry(jentry);

    // get Mu rates
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1_ =
        new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2_ =
        new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
    int mu1{-1};
    int mu2{-1};
    double mu1Pt{0};
    double mu2Pt{0};
    // std::cout << "Looping over bmtf muons." << std::endl;
    for (int it = 0; it < bmtf_->nTfMuons; ++it) {
      if (bmtf_->tfMuonBx[it] != 0) {
        continue;
      }
      if (bmtf_->tfMuonHwQual[it] < 8) {
        continue;
      }
      double muPt{(bmtf_->tfMuonHwPt[it] - 1) * 0.5};
      if (muPt > mu1Pt) {
        // std::cout << "Using mu1 from bmtf" << std::endl;
        mu2 = mu1;
        mu2Pt = mu1Pt;
        tf2_ = tf1_;
        mu1 = it;
        mu1Pt = muPt;
        tf1_ = bmtf_;
      } else if (muPt > mu2Pt) {
        // std::cout << "Using mu2 from bmtf" << std::endl;
        mu2 = it;
        mu2Pt = muPt;
        tf2_ = bmtf_;
      }
    }
    // std::cout << "Looping over omtf muons." << std::endl;
    for (int it = 0; it < omtf_->nTfMuons; ++it) {
      if (omtf_->tfMuonBx[it] != 0) {
        continue;
      }
      if (omtf_->tfMuonHwQual[it] < 8) {
        continue;
      }
      double muPt{(omtf_->tfMuonHwPt[it] - 1) * 0.5};
      if (muPt > mu1Pt) {
        mu2 = mu1;
        mu2Pt = mu1Pt;
        tf2_ = tf1_;
        mu1 = it;
        mu1Pt = muPt;
        tf1_ = omtf_;
      } else if (muPt > mu2Pt) {
        mu2 = it;
        mu2Pt = muPt;
        tf2_ = omtf_;
      }
    }
    // std::cout << "Looping over emtf muons." << std::endl;
    for (int it = 0; it < emtf_->nTfMuons; ++it) {
      if (emtf_->tfMuonBx[it] != 0) {
        continue;
      }
      if (emtf_->tfMuonHwQual[it] < 8) {
        continue;
      }
      double muPt{(emtf_->tfMuonHwPt[it] - 1) * 0.5};
      if (muPt > mu1Pt) {
        mu2 = mu1;
        mu2Pt = mu1Pt;
        tf2_ = tf1_;
        mu1 = it;
        mu1Pt = muPt;
        tf1_ = emtf_;
      } else if (muPt > mu2Pt) {
        mu2 = it;
        mu2Pt = muPt;
        tf2_ = emtf_;
      }
    }

    // std::cout << "Filling rates." << std::endl;
    // Filling di muon rates
    double mu1Eta{-100};
    double mu2Eta{-100};
    if (mu1 != -1 && mu2 != -1) {
      mu1Eta = tf1_->tfMuonHwEta[mu1] * 0.010875;
      mu2Eta = tf2_->tfMuonHwEta[mu2] * 0.010875;
      if (mu1Pt >= mu1cut && mu2Pt >= mu2cut) {
        // std::cout << "mu1: " << mu1 << " mu2: " << mu2 << std::endl;
        doubleMuRateHist.Fill(mu1Eta);
        doubleMuRateVsPtHist.Fill(mu1Pt);
      }
      //   std::cout << "mu1: " << mu1 << " mu2: " << mu2 << std::endl;
      doubleMuRateOpenHist.Fill(mu1Eta);
      doubleMuRateOpenVsPtHist.Fill(mu1Pt);
    } else {
      continue;
    }

    // Computing ghost rates
    double mu1Phi{tf1_->tfMuonGlobalPhi[mu1] * 0.010908};
    double mu2Phi{tf2_->tfMuonGlobalPhi[mu2] * 0.010908};
    if (mu1 != -1 && mu2 != -1 &&
        (!bestMatchL1toReco(reco_, mu1Eta, mu2Eta, mu1Phi, mu2Phi))) {
      if (mu1Pt >= mu1cut && mu2Pt >= mu2cut) {
        doubleMuGhostRateHist.Fill(mu1Eta);
        doubleMuGhostRateVsPtHist.Fill(mu1Pt);
      }
      doubleMuGhostRateOpenHist.Fill(mu1Eta);
      doubleMuGhostRateOpenVsPtHist.Fill(mu1Pt);
    }
  }

  std::cout << "Writing histograms to " << folder + identifier + "Hists.root"
            << std::endl;
  f->Open(folder + identifier + "Hists.root", "recreate");
  doubleMuGhostRateHist.Write();
  doubleMuGhostRateVsPtHist.Write();
  doubleMuGhostRateOpenHist.Write();
  doubleMuGhostRateOpenVsPtHist.Write();
  doubleMuRateHist.Write();
  doubleMuRateVsPtHist.Write();
  doubleMuRateOpenHist.Write();
  doubleMuRateOpenVsPtHist.Write();
  std::cout << "Closing file." << std::endl;
  // f->Close();
  std::cout << "Done writing." << std::endl;
}

void getMuonRates(int nCollBunches, int nevents, TChain* l1Chain,
                  TChain* recoChain, const int mu1cut, const int mu2cut,
                  TH1D& doubleMuGhostRateHist, TH1D& doubleMuGhostRateVsPtHist,
                  TH1D& doubleMuGhostRateOpenHist,
                  TH1D& doubleMuGhostRateOpenVsPtHist, TH1D& doubleMuRateHist,
                  TH1D& doubleMuRateVsPtHist, TH1D& doubleMuRateOpenHist,
                  TH1D& doubleMuRateOpenVsPtHist, bool retrieve_hists,
                  bool update_hists, TString folder, TString identifier) {
  TFile* f;
  if (retrieve_hists || update_hists) {
    std::cout << "Retrieving histograms.. " << std::endl;
    f = TFile::Open(folder + identifier + "Hists.root", "read");
    std::cout << "Opening file: " << folder + identifier + "Hists.root"
              << std::endl;
    doubleMuGhostRateHist =
        *(static_cast<TH1D*>(f->Get("doubleMuGhostRates" + identifier)));
    // doubleMuGhostRateHistTest->SetDirectory(0);
    doubleMuGhostRateVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuGhostRatesVsPt" + identifier)));
    // doubleMuGhostRateVsPtHist->SetDirectory(0);
    doubleMuGhostRateOpenHist =
        *(static_cast<TH1D*>(f->Get("doubleMuGhostRatesOpen" + identifier)));
    //->SetDirectory(0);
    doubleMuGhostRateOpenVsPtHist = *(
        static_cast<TH1D*>(f->Get("doubleMuGhostRatesOpenVsPt" + identifier)));
    //->SetDirectory(0);
    doubleMuRateHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRates" + identifier)));
    //->SetDirectory(0);
    doubleMuRateVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesVsPt" + identifier)));
    //->SetDirectory(0);
    doubleMuRateOpenHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesOpen" + identifier)));
    //->SetDirectory(0);
    doubleMuRateOpenVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesOpenVsPt" + identifier)));
    //->SetDirectory(0);
    f->Close();
    if (retrieve_hists) {
      return;
    }
  }

  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* l1_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_ =
      new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  l1Chain->SetBranchAddress("L1Upgrade", &l1_);
  recoChain->SetBranchAddress("Muon", &reco_);

  for (Long64_t jentry = 0; jentry < nevents; ++jentry) {
    if ((jentry % 1000) == 0)
      std::cout << "Done " << jentry << " events..." << std::endl;

    l1Chain->GetEntry(jentry);
    if (l1_->nMuons < 2) {
      continue;
    }
    recoChain->GetEntry(jentry);

    // get Mu rates
    int mu1(-1);
    int mu2(-1);
    double mu1Pt(0);
    double mu2Pt(0);
    for (uint it = 0; it < l1_->nMuons; ++it) {
      if (l1_->muonBx[it] != 0) {
        continue;
      }
      if (l1_->muonQual[it] < 8) {
        continue;
      }
      if (l1_->muonEt[it] > mu1Pt) {
        mu2 = mu1;
        mu2Pt = mu1Pt;
        mu1 = it;
        mu1Pt = l1_->muonEt[it];
      } else if (l1_->muonEt[it] > mu2Pt) {
        mu2 = it;
        mu2Pt = l1_->muonEt[it];
      }
    }

    // Filling di muon rates
    if (mu1 != -1 && mu2 != -1) {
      if (mu1Pt >= mu1cut && mu2Pt >= mu2cut) {
        doubleMuRateHist.Fill(l1_->muonEta[mu1]);
        doubleMuRateVsPtHist.Fill(l1_->muonEt[mu1]);
      }
      doubleMuRateOpenHist.Fill(l1_->muonEta[mu1]);
      doubleMuRateOpenVsPtHist.Fill(l1_->muonEt[mu1]);
    }

    // Computing ghost rates
    if (mu1 != -1 && mu2 != -1 &&
        (!bestMatchL1toReco(reco_, l1_->muonEta[mu1], l1_->muonEta[mu2],
                            l1_->muonPhi[mu1], l1_->muonPhi[mu2]))) {
      if (mu1Pt >= mu1cut && mu2Pt >= mu2cut) {
        doubleMuGhostRateHist.Fill(l1_->muonEta[mu1]);
        doubleMuGhostRateVsPtHist.Fill(l1_->muonEt[mu1]);
      }
      doubleMuGhostRateOpenHist.Fill(l1_->muonEta[mu1]);
      doubleMuGhostRateOpenVsPtHist.Fill(l1_->muonEt[mu1]);
    }
  }

  std::cout << "Writing histograms to " << folder + identifier + "Hists.root"
            << std::endl;
  f->Open(folder + identifier + "Hists.root", "recreate");
  doubleMuGhostRateHist.Write();
  doubleMuGhostRateVsPtHist.Write();
  doubleMuGhostRateOpenHist.Write();
  doubleMuGhostRateOpenVsPtHist.Write();
  doubleMuRateHist.Write();
  doubleMuRateVsPtHist.Write();
  doubleMuRateOpenHist.Write();
  doubleMuRateOpenVsPtHist.Write();
  std::cout << "Closing file." << std::endl;
  // f->Close();
  std::cout << "Done writing." << std::endl;
}

void calcRates(int nCollBunches, int nevents, TH1D& doubleMuGhostRateHist,
               TH1D& doubleMuGhostRateVsPtHist, TH1D& doubleMuGhostRateOpenHist,
               TH1D& doubleMuGhostRateOpenVsPtHist, TH1D& doubleMuRateHist,
               TH1D& doubleMuRateVsPtHist, TH1D& doubleMuRateOpenHist,
               TH1D& doubleMuRateOpenVsPtHist) {
  // normalisation factor
  double norm =
      (11. * nCollBunches) / nevents;  // zb rate = n_colliding * 11 kHz
  std::cout << "Number of events examined: " << nevents << "\n";
  std::cout << "norm = " << norm << "\n";

  std::cout << "###########################################"
            << "\n";
  std::cout << "** Computed rates: **"
            << "\n";
  std::cout << "Rate with given pT cuts: "
            << doubleMuRateHist.GetEntries() * norm << "+/-"
            << TMath::Sqrt(doubleMuRateHist.GetEntries()) * norm << "\n";
  std::cout << "DoubleMuOpen rate: " << doubleMuRateOpenHist.GetEntries() * norm
            << "+/-" << TMath::Sqrt(doubleMuRateOpenHist.GetEntries()) * norm
            << "\n";
  std::cout << "Ghost rate with given pT cuts: "
            << doubleMuGhostRateHist.GetEntries() * norm << "+/-"
            << TMath::Sqrt(doubleMuGhostRateHist.GetEntries()) * norm << "\n";
  std::cout << "DoubleMuOpen ghost rate: "
            << doubleMuGhostRateOpenHist.GetEntries() * norm << "+/-"
            << TMath::Sqrt(doubleMuGhostRateOpenHist.GetEntries()) * norm
            << "\n";
  std::cout << "###########################################"
            << "\n";

  doubleMuRateHist.Sumw2();
  doubleMuRateVsPtHist.Sumw2();
  doubleMuRateOpenHist.Sumw2();
  doubleMuRateOpenVsPtHist.Sumw2();
  doubleMuGhostRateHist.Sumw2();
  doubleMuGhostRateVsPtHist.Sumw2();
  doubleMuGhostRateOpenHist.Sumw2();
  doubleMuGhostRateOpenVsPtHist.Sumw2();
  TF1* constant = new TF1("constant", "1", -5, 5);
  doubleMuRateHist.Multiply(constant, norm);
  doubleMuRateVsPtHist.Multiply(constant, norm);
  doubleMuRateOpenHist.Multiply(constant, norm);
  doubleMuRateOpenVsPtHist.Multiply(constant, norm);
  doubleMuGhostRateHist.Multiply(constant, norm);
  doubleMuGhostRateVsPtHist.Multiply(constant, norm);
  doubleMuGhostRateOpenHist.Multiply(constant, norm);
  doubleMuGhostRateOpenVsPtHist.Multiply(constant, norm);
}

void drawHistograms(TH1D& noCouHist, TH1D& baselineHist, TH1D& conservativeHist,
                    TH1D& aggressiveHist, TString filename, TString xAxisLabel,
                    TString descString, TString plotFolder, TString run,
                    bool logy = false) {
  TLatex n1;
  n1.SetNDC();
  n1.SetTextFont(52);
  n1.SetTextSize(0.035);

  TLatex n2;
  n2.SetNDC();
  n2.SetTextFont(52);
  n2.SetTextSize(0.035);

  TCanvas c1("c1", "", 525, 500);
  if (logy) {
    c1.SetLogy();
  }
  // c1.SetLeftMargin(0.18);

  //  muRatesOpenUnpack->SetLineWidth(2);

  noCouHist.SetLineColor(kRed + 1);
  noCouHist.SetMarkerStyle(25);
  noCouHist.SetMarkerColor(kRed + 1);
  noCouHist.Draw("E1HIST");
  noCouHist.GetXaxis()->SetTitle(xAxisLabel);
  noCouHist.GetYaxis()->SetTitle("Rate [kHz/bin]");
  // noCouHist.GetYaxis()->SetTitleOffset(1.5);

  aggressiveHist.SetLineColor(kGreen + 2);
  aggressiveHist.SetMarkerStyle(32);
  aggressiveHist.SetMarkerColor(kGreen + 2);
  aggressiveHist.Draw("same,E1HIST");
  aggressiveHist.GetXaxis()->SetTitle(xAxisLabel);
  aggressiveHist.GetYaxis()->SetTitle("Rate [kHz/bin]");
  // aggressiveHist.GetYaxis()->SetTitleOffset(1.5);

  baselineHist.SetLineColor(kCyan + 1);
  baselineHist.SetMarkerStyle(26);
  baselineHist.SetMarkerColor(kCyan + 1);
  baselineHist.Draw("same,E1HIST");
  baselineHist.GetXaxis()->SetTitle(xAxisLabel);
  baselineHist.GetYaxis()->SetTitle("Rate [kHz/bin]");
  // baselineHist.GetYaxis()->SetTitleOffset(1.5);

  conservativeHist.SetLineColor(kOrange + 7);
  conservativeHist.SetMarkerStyle(27);
  conservativeHist.SetMarkerColor(kOrange + 7);
  conservativeHist.Draw("same,E1HIST");
  conservativeHist.GetXaxis()->SetTitle(xAxisLabel);
  conservativeHist.GetYaxis()->SetTitle("Rate [kHz/bin]");
  // conservativeHist.GetYaxis()->SetTitleOffset(1.5);

  gPad->Modified();

  TLegend leg1(0.3, 0.72, 0.7, 0.92);
  leg1.SetFillColor(0);
  leg1.AddEntry(&noCouHist, "No uGMT cancel-out", "lp");
  leg1.AddEntry(&baselineHist, "Baseline tuning", "lp");
  leg1.AddEntry(&conservativeHist, "Conservative tuning", "lp");
  leg1.AddEntry(&aggressiveHist, "Aggressive tuning", "lp");
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.Draw();
  n1.DrawLatex(0.3, 0.68, "Run " + run + ", #sqrt{s} = 13 TeV");
  n2.DrawLatex(0.3, 0.63, descString);

  c1.SaveAs(plotFolder + filename + ".pdf");
  c1.SaveAs(plotFolder + filename + ".png");
}
