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
void getMuonRates(int nCollBunches, int nevents, TChain* l1Chain,
                  TChain* recoChain, const int pT1cut, const int pT2cut,
                  TH1D& doubleMuGhostRateHist,
                  TH1D& doubleMuGhostRateTrailingHist,
                  TH1D& doubleMuGhostRateOpenHist,
                  TH1D& doubleMuGhostRateOpenTrailingHist,
                  TH1D& doubleMuRateHist, TH1D& doubleMuRateTrailingHist,
                  TH1D& doubleMuRateOpenHist,
                  TH1D& doubleMuRateOpenTrailingHist);
void drawHistograms(TH1D& baselineHist, TH1D& conservativeHist,
                    TH1D& aggressiveHist, TString filename, TString xAxisLabel,
                    TString descString, TString plotFolder, TString run);

void compareDiMuRates(const char* file_list_baseline,
                      const char* file_list_conservative,
                      const char* file_list_aggressive, TString folder = "tmp",
                      TString run = "XXX", int mu1cut = 11, int mu2cut = 4,
                      int nCollBunches = 2028) {
  TString plotFolder = "plots/" + run + "/" + folder + "/";
  mkdir("plots/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("plots/run/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir(plotFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  setTDRStyle();
  writeExtraText = true;  // Add "Preliminary"
  // lumi_13TeV  = "4.9 fb^{-1}";
  lumi_sqrtS = "13 TeV";  // used with iPeriod = 0, e.g. for simulation-only
                          // plots (default is an empty string)

  int iPeriod = 0;  // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses
                    // lumi_sqrtS)

  std::vector<std::string> listNtuplesBaseline;
  std::vector<std::string> listNtuplesConservative;
  std::vector<std::string> listNtuplesAggressive;

  std::string l1Treepath("l1UpgradeEmuTree/L1UpgradeTree");
  std::string recoTreepath("l1MuonRecoTree/Muon2RecoTree");

  bool success = readFList(file_list_baseline, listNtuplesBaseline, l1Treepath,
                           recoTreepath);
  success &= readFList(file_list_conservative, listNtuplesConservative,
                       l1Treepath, recoTreepath);
  success &= readFList(file_list_aggressive, listNtuplesAggressive, l1Treepath,
                       recoTreepath);

  if (!success) {
    std::cout << "Couldn't read NTuple file list. Exiting.. " << std::endl;
    return;
  }

  TChain* chainL1Baseline = new TChain(l1Treepath.c_str());
  TChain* chainRecoBaseline = new TChain(recoTreepath.c_str());
  TChain* chainL1Conservative = new TChain(l1Treepath.c_str());
  TChain* chainRecoConservative = new TChain(recoTreepath.c_str());
  TChain* chainL1Aggressive = new TChain(l1Treepath.c_str());
  TChain* chainRecoAggressive = new TChain(recoTreepath.c_str());

  int baselineEntries{0};
  int conservativeEntries{0};
  int aggressiveEntries{0};

  baselineEntries =
      setupTChain(listNtuplesBaseline, chainL1Baseline, chainRecoBaseline);
  conservativeEntries = setupTChain(listNtuplesConservative,
                                    chainL1Conservative, chainRecoConservative);
  aggressiveEntries = setupTChain(listNtuplesAggressive, chainL1Aggressive,
                                  chainRecoAggressive);

  // mu bins
  int nMuBins = 25;
  float muLo = -2.5;
  float muHi = 2.5;
  float muBinWidth = (muHi - muLo) / nMuBins;
  int nTFbin = 110;
  int tflow = 0;
  int tfhigh = 109;
  int tfBinWidth = 1;

  // make histos
  TH1D doubleMuGhostRatesBaseline("doubleMuGhostRatesBaseline", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesTrailingBaseline("doubleMuGhostRatesTrailingBaseline",
                                          "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenBaseline("doubleMuGhostRatesOpenBaseline", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenTrailingBaseline(
      "doubleMuGhostRatesOpenTrailingBaseline", "", nMuBins, muLo - 0.1,
      muHi + 0.1);
  TH1D doubleMuRatesBaseline("doubleMuRatesBaseline", "", nMuBins, muLo - 0.1,
                             muHi + 0.1);
  TH1D doubleMuRatesTrailingBaseline("doubleMuRatesTrailingBaseline", "",
                                     nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenBaseline("doubleMuRatesOpenBaseline", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenTrailingBaseline("doubleMuRatesOpenTrailingBaseline",
                                         "", nMuBins, muLo - 0.1, muHi + 0.1);

  TH1D doubleMuGhostRatesConservative("doubleMuGhostRatesConservative", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesTrailingConservative(
      "doubleMuGhostRatesTrailingConservative", "", nMuBins, muLo - 0.1,
      muHi + 0.1);
  TH1D doubleMuGhostRatesOpenConservative("doubleMuGhostRatesOpenConservative",
                                          "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenTrailingConservative(
      "doubleMuGhostRatesOpenTrailingConservative", "", nMuBins, muLo - 0.1,
      muHi + 0.1);
  TH1D doubleMuRatesConservative("doubleMuRatesConservative", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesTrailingConservative("doubleMuRatesTrailingConservative",
                                         "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenConservative("doubleMuRatesOpenConservative", "",
                                     nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenTrailingConservative(
      "doubleMuRatesOpenTrailingConservative", "", nMuBins, muLo - 0.1,
      muHi + 0.1);

  TH1D doubleMuGhostRatesAggressive("doubleMuGhostRatesAggressive", "", nMuBins,
                                    muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesTrailingAggressive(
      "doubleMuGhostRatesTrailingAggressive", "", nMuBins, muLo - 0.1,
      muHi + 0.1);
  TH1D doubleMuGhostRatesOpenAggressive("doubleMuGhostRatesOpenAggressive", "",
                                        nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuGhostRatesOpenTrailingAggressive(
      "doubleMuGhostRatesOpenTrailingAggressive", "", nMuBins, muLo - 0.1,
      muHi + 0.1);
  TH1D doubleMuRatesAggressive("doubleMuRatesAggressive", "", nMuBins,
                               muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesTrailingAggressive("doubleMuRatesTrailingAggressive", "",
                                       nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenAggressive("doubleMuRatesOpenAggressive", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D doubleMuRatesOpenTrailingAggressive(
      "doubleMuRatesOpenTrailingAggressive", "", nMuBins, muLo - 0.1,
      muHi + 0.1);

  getMuonRates(
      nCollBunches, baselineEntries, chainL1Baseline, chainRecoBaseline, mu1cut,
      mu2cut, doubleMuGhostRatesBaseline, doubleMuGhostRatesTrailingBaseline,
      doubleMuGhostRatesOpenBaseline, doubleMuGhostRatesOpenTrailingBaseline,
      doubleMuRatesBaseline, doubleMuRatesTrailingBaseline,
      doubleMuRatesOpenBaseline, doubleMuRatesOpenTrailingBaseline);

  getMuonRates(
      nCollBunches, conservativeEntries, chainL1Conservative,
      chainRecoConservative, mu1cut, mu2cut, doubleMuGhostRatesConservative,
      doubleMuGhostRatesTrailingConservative,
      doubleMuGhostRatesOpenConservative,
      doubleMuGhostRatesOpenTrailingConservative, doubleMuRatesConservative,
      doubleMuRatesTrailingConservative, doubleMuRatesOpenConservative,
      doubleMuRatesOpenTrailingConservative);

  getMuonRates(
      nCollBunches, aggressiveEntries, chainL1Aggressive, chainRecoAggressive,
      mu1cut, mu2cut, doubleMuGhostRatesAggressive,
      doubleMuGhostRatesTrailingAggressive, doubleMuGhostRatesOpenAggressive,
      doubleMuGhostRatesOpenTrailingAggressive, doubleMuRatesAggressive,
      doubleMuRatesTrailingAggressive, doubleMuRatesOpenAggressive,
      doubleMuRatesOpenTrailingAggressive);

  std::ostringstream oss;
  oss << "ZeroBias, L1_DoubleMu_" << mu1cut << "_" << mu2cut;
  drawHistograms(doubleMuRatesBaseline, doubleMuRatesConservative,
                 doubleMuRatesAggressive, "doubleMuRatesDoubleMuonLeading",
                 "#eta (leading #mu)", oss.str().c_str(), plotFolder, run);
  drawHistograms(
      doubleMuRatesTrailingBaseline, doubleMuRatesTrailingConservative,
      doubleMuRatesTrailingAggressive, "doubleMuRatesDoubleMuonTrailing",
      "#eta (trailing #mu)", oss.str().c_str(), plotFolder, run);
  drawHistograms(doubleMuRatesOpenBaseline, doubleMuRatesOpenConservative,
                 doubleMuRatesOpenAggressive,
                 "doubleMuOpenRatesDoubleMuonLeading", "#eta (leading #mu)",
                 "Zero Bias, L1_SingleMuOpen", plotFolder, run);
  drawHistograms(doubleMuRatesOpenTrailingBaseline,
                 doubleMuRatesOpenTrailingConservative,
                 doubleMuRatesOpenTrailingAggressive,
                 "doubleMuOpenRatesDoubleMuonTrailing", "#eta (trailing #mu)",
                 "Zero Bias, L1_SingleMuOpen", plotFolder, run);

  oss << ", ghost rate";
  drawHistograms(doubleMuGhostRatesBaseline, doubleMuGhostRatesConservative,
                 doubleMuGhostRatesAggressive, "ghostRatesDoubleMuonLeading",
                 "#eta (leading #mu)", oss.str().c_str(), plotFolder, run);
  drawHistograms(doubleMuGhostRatesTrailingBaseline,
                 doubleMuGhostRatesTrailingConservative,
                 doubleMuGhostRatesTrailingAggressive,
                 "ghostRatesDiMuonTrailing", "#eta (trailing #mu)",
                 oss.str().c_str(), plotFolder, run);
  drawHistograms(doubleMuGhostRatesOpenBaseline,
                 doubleMuGhostRatesOpenConservative,
                 doubleMuGhostRatesOpenAggressive,
                 "ghostRatesDoubleMuonOpenLeading", "#eta (leading #mu)",
                 "Zero Bias, L1_DoubleMu0, ghost rate", plotFolder, run);
  drawHistograms(doubleMuGhostRatesOpenTrailingBaseline,
                 doubleMuGhostRatesOpenTrailingConservative,
                 doubleMuGhostRatesOpenTrailingAggressive,
                 "ghostRatesDoubleMuonOpenTrailing", "#eta (trailing #mu)",
                 "Zero Bias, L1_DoubleMu0, ghost rate", plotFolder, run);
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
    std::cout << "L1Upgrade tree not found.. " << std::endl;
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

void getMuonRates(int nCollBunches, int nevents, TChain* l1Chain,
                  TChain* recoChain, const int mu1cut, const int mu2cut,
                  TH1D& doubleMuGhostRateHist,
                  TH1D& doubleMuGhostRateTrailingHist,
                  TH1D& doubleMuGhostRateOpenHist,
                  TH1D& doubleMuGhostRateOpenTrailingHist,
                  TH1D& doubleMuRateHist, TH1D& doubleMuRateTrailingHist,
                  TH1D& doubleMuRateOpenHist,
                  TH1D& doubleMuRateOpenTrailingHist) {
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
        doubleMuRateTrailingHist.Fill(l1_->muonEta[mu2]);
      }
      doubleMuRateOpenHist.Fill(l1_->muonEta[mu1]);
      doubleMuRateOpenTrailingHist.Fill(l1_->muonEta[mu2]);
    }

    // Computing ghost rates
    int nRecoMus = 0;
    for (int i = 0; i < reco_->nMuons; ++i) {
      if (reco_->isTightMuon[i]) {
        ++nRecoMus;
      }
    }

    if (mu1 != -1 && mu2 != -1 && nRecoMus == 1) {
      if (mu1Pt >= mu1cut && mu2Pt >= mu2cut) {
        doubleMuGhostRateHist.Fill(l1_->muonEta[mu1]);
        doubleMuGhostRateTrailingHist.Fill(l1_->muonEta[mu2]);
      }
      doubleMuGhostRateOpenHist.Fill(l1_->muonEta[mu1]);
      doubleMuGhostRateOpenTrailingHist.Fill(l1_->muonEta[mu2]);
    }
  }

  // normalisation factor
  double norm =
      (11. * nCollBunches) / nevents;  // zb rate = n_colliding * 11 kHz
  std::cout << "norm = " << norm << std::endl;

  std::cout << "###########################################" << std::endl;
  std::cout << "** Computed rates: **" << std::endl;
  std::cout << "Rate with given pT cuts: "
            << doubleMuRateHist.GetEntries() * norm << std::endl;
  std::cout << "DoubleMuOpen rate: " << doubleMuRateOpenHist.GetEntries() * norm
            << std::endl;
  std::cout << "Ghost rate with given pT cuts: "
            << doubleMuGhostRateHist.GetEntries() * norm << std::endl;
  std::cout << "DoubleMuOpen ghost rate: "
            << doubleMuGhostRateOpenHist.GetEntries() * norm << std::endl;
  std::cout << "###########################################" << std::endl;

  doubleMuRateHist.Sumw2();
  doubleMuRateTrailingHist.Sumw2();
  doubleMuRateOpenHist.Sumw2();
  doubleMuRateOpenTrailingHist.Sumw2();
  doubleMuGhostRateHist.Sumw2();
  doubleMuGhostRateTrailingHist.Sumw2();
  doubleMuGhostRateOpenHist.Sumw2();
  doubleMuGhostRateOpenTrailingHist.Sumw2();
  TF1* constant = new TF1("constant", "1", -5, 5);
  doubleMuRateHist.Multiply(constant, norm);
  doubleMuRateTrailingHist.Multiply(constant, norm);
  doubleMuRateOpenHist.Multiply(constant, norm);
  doubleMuRateOpenTrailingHist.Multiply(constant, norm);
  doubleMuGhostRateHist.Multiply(constant, norm);
  doubleMuGhostRateTrailingHist.Multiply(constant, norm);
  doubleMuGhostRateOpenHist.Multiply(constant, norm);
  doubleMuGhostRateOpenTrailingHist.Multiply(constant, norm);
}

void drawHistograms(TH1D& baselineHist, TH1D& conservativeHist,
                    TH1D& aggressiveHist, TString filename, TString xAxisLabel,
                    TString descString, TString plotFolder, TString run) {
  TLatex n1;
  n1.SetNDC();
  n1.SetTextFont(52);
  n1.SetTextSize(0.035);

  TLatex n2;
  n2.SetNDC();
  n2.SetTextFont(52);
  n2.SetTextSize(0.035);

  TCanvas c1;

  //  muRatesOpenUnpack->SetLineWidth(2);

  baselineHist.SetLineColor(kOrange);
  baselineHist.GetXaxis()->SetTitle(xAxisLabel);
  baselineHist.GetYaxis()->SetTitle("Rate");
  baselineHist.SetMarkerStyle(23);
  baselineHist.SetMarkerColor(kOrange);
  baselineHist.Draw("E1HIST");
  baselineHist.GetXaxis()->SetTitle(xAxisLabel);
  baselineHist.GetYaxis()->SetTitle("Rate [kHz]");

  conservativeHist.SetLineColor(kBlue + 2);
  conservativeHist.GetXaxis()->SetTitle(xAxisLabel);
  conservativeHist.GetYaxis()->SetTitle("Rate");
  conservativeHist.SetMarkerStyle(20);
  conservativeHist.SetMarkerColor(kBlue + 2);
  conservativeHist.Draw("same,E1HIST");
  conservativeHist.GetXaxis()->SetTitle(xAxisLabel);
  conservativeHist.GetYaxis()->SetTitle("Rate [kHz]");

  aggressiveHist.SetLineColor(kGreen + 2);
  aggressiveHist.GetXaxis()->SetTitle(xAxisLabel);
  aggressiveHist.GetYaxis()->SetTitle("Rate");
  aggressiveHist.SetMarkerStyle(21);
  aggressiveHist.SetMarkerColor(kGreen + 2);
  aggressiveHist.Draw("same,E1HIST");
  aggressiveHist.GetXaxis()->SetTitle(xAxisLabel);
  aggressiveHist.GetYaxis()->SetTitle("Rate [kHz]");

  gPad->Modified();

  TLegend leg1(0.3, 0.72, 0.7, 0.92);
  leg1.SetFillColor(0);
  leg1.AddEntry(&baselineHist, "Baseline tuning", "lp");
  leg1.AddEntry(&conservativeHist, "Conservative tuning", "lp");
  leg1.AddEntry(&aggressiveHist, "Aggressive tuning", "lp");
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.Draw();
  n1.DrawLatex(0.3, 0.68, "Run " + run + " #sqrt{s} = 13 TeV");
  n2.DrawLatex(0.3, 0.63, descString);

  c1.SaveAs(plotFolder + filename + ".pdf");
  c1.SaveAs(plotFolder + filename + ".png");
}
