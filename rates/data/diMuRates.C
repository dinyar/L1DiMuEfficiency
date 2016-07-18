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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

// TODO: Make rate plots for uGMT inputs!

void diMuRates(const char* fname = "file_list_275125", TString folder = "tmp",
               TString run = "XXX", int mu1cut = 11, int mu2cut = 4,
               bool drawReEmu = false) {
  TString plotFolder = "plots/" + run + "/" + folder + "/";

  gStyle->SetOptStat(0);

  std::vector<std::string> listNtuples;

  // OpenNtupleList
  std::ifstream flist(fname);
  if (!flist) {
    std::cout << "File " << fname << " is not found!" << std::endl;
    return;
  }

  while (!flist.eof()) {
    std::string str;
    getline(flist, str);
    if (!flist.fail()) {
      if (str != "") listNtuples.push_back(str);
    }
  }

  // CheckFirstFile
  if (listNtuples.size() == 0) return;

  TFile* rf = TFile::Open(listNtuples[0].c_str());

  if (rf == 0) return;
  if (rf->IsOpen() == 0) return;

  std::string unpackTreepath("l1UpgradeTree/L1UpgradeTree");
  std::string reEmuTreepath("l1UpgradeEmuTree/L1UpgradeTree");
  std::string recoTreepath("l1MuonRecoTree/Muon2RecoTree");
  TTree* treeL1Unpack = (TTree*)rf->Get(unpackTreepath.c_str());
  TTree* treeL1reEmu = (TTree*)rf->Get(reEmuTreepath.c_str());
  TTree* treeReco = (TTree*)rf->Get(recoTreepath.c_str());

  if (!treeL1reEmu || !treeL1Unpack) {
    std::cout << "L1Upgrade trees not found.. " << std::endl;
    return;
  }
  if (!treeReco) {
    std::cout << "Reco tree not found.. " << std::endl;
    return;
  }

  // OpenWithoutInit
  TChain* l1UnpackChain = new TChain(unpackTreepath.c_str());
  TChain* l1reEmuChain = new TChain(reEmuTreepath.c_str());
  TChain* chainReco = new TChain(recoTreepath.c_str());
  for (unsigned int i = 0; i < listNtuples.size(); i++) {
    std::cout << " -- Adding " << listNtuples[i] << std::endl;
    l1UnpackChain->Add(listNtuples[i].c_str());
    l1reEmuChain->Add(listNtuples[i].c_str());
    chainReco->Add(listNtuples[i].c_str());
  }

  // Init
  std::cout << "Estimate the number of entries ..." << std::endl;
  int nentries = l1UnpackChain->GetEntries();
  std::cout << nentries << std::endl;
  int nevents = nentries;

  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* unpack_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisL1UpgradeDataFormat* reEmu_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_ =
      new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  l1UnpackChain->SetBranchAddress("L1Upgrade", &unpack_);
  l1reEmuChain->SetBranchAddress("L1Upgrade", &reEmu_);
  chainReco->SetBranchAddress("Muon", &reco_);

  // mu bins
  int nMuBins = 25;
  float muLo = -2.5;
  float muHi = 2.5;
  float muBinWidth = (muHi - muLo) / nMuBins;
  int nTFbin = 12;
  int tflow = 0;
  int tfhigh = 12;
  int tfBinWidth = 1;

  // make histos
  TH1D* singleMuRatesOpenUnpack =
      new TH1D("singleMuRatesOpenUnpack", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* singleMuRatesOpenReEmu =
      new TH1D("singleMuRatesOpenReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muGhostRatesUnpack =
      new TH1D("muGhostRatesUnpack", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muGhostRatesReEmu =
      new TH1D("muGhostRatesReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muGhostRatesOpenUnpack =
      new TH1D("muGhostRatesOpenUnpack", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muGhostRatesOpenReEmu =
      new TH1D("muGhostRatesOpenReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muGhostRatesOpenTrailingUnpack = new TH1D(
      "muGhostRatesOpenTrailingUnpack", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muGhostRatesOpenTrailingReEmu = new TH1D(
      "muGhostRatesOpenTrailingReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muRatesUnpack =
      new TH1D("muRatesUnpack", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muRatesReEmu =
      new TH1D("muRatesReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muRatesOpenUnpack =
      new TH1D("muRatesOpenUnpack", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muRatesOpenReEmu =
      new TH1D("muRatesOpenReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muRatesOpenTrailingUnpack = new TH1D("muRatesOpenTrailingUnpack", "",
                                             nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D* muRatesOpenTrailingReEmu =
      new TH1D("muRatesOpenTrailingReEmu", "", nMuBins, muLo - 0.1, muHi + 0.1);

  // TODO: eta vs. pt and eta vs tfMuIdx plots?
  TH1D* tfMuIdxUnpack = new TH1D("tfMuIdxUnpack", "", 110, 0, 109);
  TH1D* tfMuIdxReEmu = new TH1D("tfMuIdxReEmu", "", 110, 0, 109);
  TH1D* tfMuIdxAllUnpack = new TH1D("tfMuIdxAllUnpack", "", 110, 0, 109);
  TH1D* tfMuIdxAllReEmu = new TH1D("tfMuIdxAllReEmu", "", 110, 0, 109);

  int unpackCounts(0);
  int reEmuCounts(0);

  for (Long64_t jentry = 0; jentry < nevents; jentry++) {
    if ((jentry % 1000) == 0)
      std::cout << "Done " << jentry << " events..." << std::endl;

    l1UnpackChain->GetEntry(jentry);
    l1reEmuChain->GetEntry(jentry);
    chainReco->GetEntry(jentry);

    // get Mu rates
    int mu1Unpack(-1);
    int mu2Unpack(-1);
    double mu1PtUnpack(0);
    double mu2PtUnpack(0);
    int mu1ReEmu(-1);
    int mu2ReEmu(-1);
    double mu1PtReEmu(0);
    double mu2PtReEmu(0);
    for (uint it = 0; it < unpack_->nMuons; ++it) {
      if (unpack_->muonBx[it] != 0) {
        continue;
      }
      tfMuIdxAllUnpack->Fill(unpack_->muonTfMuonIdx[it]);

      if (unpack_->muonQual[it] <= 4) {
        continue;
      }
      if (unpack_->muonEt[it] > mu1PtUnpack) {
        mu2Unpack = mu1Unpack;
        mu2PtUnpack = mu1PtUnpack;
        mu1Unpack = it;
        mu1PtUnpack = unpack_->muonEt[it];
      } else if (unpack_->muonEt[it] > mu2PtUnpack) {
        mu2Unpack = it;
        mu2PtUnpack = unpack_->muonEt[it];
      }
    }
    if (drawReEmu) {
      for (uint it = 0; it < reEmu_->nMuons; ++it) {
        if (reEmu_->muonBx[it] != 0) {
          continue;
        }
        tfMuIdxAllReEmu->Fill(reEmu_->muonTfMuonIdx[it]);
        if (reEmu_->muonQual[it] <= 4) {
          continue;
        }
        if (reEmu_->muonEt[it] > mu1PtReEmu) {
          mu2ReEmu = mu1ReEmu;
          mu2PtReEmu = mu1PtReEmu;
          mu1ReEmu = it;
          mu1PtReEmu = reEmu_->muonEt[it];
        } else if (reEmu_->muonEt[it] > mu2PtReEmu) {
          mu2ReEmu = it;
          mu2PtReEmu = reEmu_->muonEt[it];
        }
      }
    }
    // Computing single muon rates
    if (mu1Unpack != -1) {
      singleMuRatesOpenUnpack->Fill(unpack_->muonEta[mu1Unpack]);
    }
    if (drawReEmu) {
      if (mu1ReEmu != -1) {
        singleMuRatesOpenReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
      }
    }

    // Computing di muon rates
    if (mu1Unpack != -1 && mu2Unpack != -1) {
      if (mu1PtUnpack >= mu1cut && mu2PtUnpack >= mu2cut) {
        muRatesUnpack->Fill(unpack_->muonEta[mu1Unpack]);
      }
      muRatesOpenUnpack->Fill(unpack_->muonEta[mu1Unpack]);
      muRatesOpenTrailingUnpack->Fill(unpack_->muonEta[mu2Unpack]);
      tfMuIdxUnpack->Fill(unpack_->muonTfMuonIdx[mu1Unpack]);
      ++unpackCounts;
    }
    if (drawReEmu) {
      if (mu1ReEmu != -1 && mu2ReEmu != -1) {
        if (mu1PtReEmu >= mu1cut && mu2PtReEmu >= mu2cut) {
          muRatesReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
        }
        muRatesOpenReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
        muRatesOpenTrailingReEmu->Fill(reEmu_->muonEta[mu2ReEmu]);
        tfMuIdxReEmu->Fill(reEmu_->muonTfMuonIdx[mu1ReEmu]);
        ++reEmuCounts;
      }
    }

    // Computing ghost rates
    int nRecoMus = 0;
    for (int i = 0; i < reco_->nMuons; ++i) {
      if (reco_->isTightMuon[i]) {
        ++nRecoMus;
      }
    }

    if (mu1Unpack != -1 && mu2Unpack != -1 && nRecoMus == 1) {
      if (mu1PtUnpack >= mu1cut && mu2PtUnpack >= mu2cut) {
        muGhostRatesUnpack->Fill(unpack_->muonEta[mu1Unpack]);
      }
      muGhostRatesOpenUnpack->Fill(unpack_->muonEta[mu1Unpack]);
      muGhostRatesOpenTrailingUnpack->Fill(unpack_->muonEta[mu2Unpack]);
      ++unpackCounts;
    }
    if (drawReEmu) {
      if (mu1ReEmu != -1 && mu2ReEmu != -1 && nRecoMus == 1) {
        if (mu1PtReEmu >= mu1cut && mu2PtReEmu >= mu2cut) {
          muGhostRatesReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
        }
        muGhostRatesOpenReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
        muGhostRatesOpenTrailingReEmu->Fill(reEmu_->muonEta[mu2ReEmu]);
        ++reEmuCounts;
      }
    }
  }

  // normalisation factor
  double norm = (11. * 2028.) / nevents;  // zb rate = n_colliding * 11 kHz
  std::cout << "norm = " << norm << std::endl;

  std::cout << "###########################################" << std::endl;
  std::cout << "** Computed rates: **" << std::endl;
  std::cout << "Unpacked rate: " << muRatesUnpack->GetEntries() * norm
            << std::endl;
  std::cout << "Unpacked open rate: " << muRatesOpenUnpack->GetEntries() * norm
            << std::endl;
  std::cout << "ReEmulated rate: " << muRatesReEmu->GetEntries() * norm
            << std::endl;
  std::cout << "ReEmulated open rate: " << muRatesOpenReEmu->GetEntries() * norm
            << std::endl;
  std::cout << "Unpacked SingleMuOpen rate: "
            << singleMuRatesOpenUnpack->GetEntries() * norm << std::endl;
  std::cout << "ReEmulated SingleMuOpen rate: "
            << singleMuRatesOpenReEmu->GetEntries() * norm << std::endl;
  std::cout << "###########################################" << std::endl;

  muRatesUnpack->Sumw2();
  muRatesReEmu->Sumw2();
  muRatesOpenUnpack->Sumw2();
  muRatesOpenReEmu->Sumw2();
  muRatesOpenTrailingUnpack->Sumw2();
  muRatesOpenTrailingReEmu->Sumw2();
  singleMuRatesOpenUnpack->Sumw2();
  singleMuRatesOpenReEmu->Sumw2();
  muGhostRatesUnpack->Sumw2();
  muGhostRatesReEmu->Sumw2();
  muGhostRatesOpenUnpack->Sumw2();
  muGhostRatesOpenReEmu->Sumw2();
  muGhostRatesOpenTrailingUnpack->Sumw2();
  muGhostRatesOpenTrailingReEmu->Sumw2();
  TF1* constant = new TF1("constant", "1", -5, 5);
  muRatesUnpack->Multiply(constant, norm);
  muRatesReEmu->Multiply(constant, norm);
  muRatesOpenUnpack->Multiply(constant, norm);
  muRatesOpenReEmu->Multiply(constant, norm);
  muRatesOpenTrailingUnpack->Multiply(constant, norm);
  muRatesOpenTrailingReEmu->Multiply(constant, norm);
  singleMuRatesOpenUnpack->Multiply(constant, norm);
  singleMuRatesOpenReEmu->Multiply(constant, norm);
  muGhostRatesUnpack->Multiply(constant, norm);
  muGhostRatesReEmu->Multiply(constant, norm);
  muGhostRatesOpenUnpack->Multiply(constant, norm);
  muGhostRatesOpenReEmu->Multiply(constant, norm);
  muGhostRatesOpenTrailingUnpack->Multiply(constant, norm);
  muGhostRatesOpenTrailingReEmu->Multiply(constant, norm);

  TLatex n3;
  n3.SetNDC();
  n3.SetTextFont(52);
  n3.SetTextSize(0.04);

  TLatex n4;
  n4.SetNDC();
  n4.SetTextFont(52);
  n4.SetTextSize(0.04);

  TCanvas* c1 = new TCanvas;

  //  muRatesOpenUnpack->SetLineWidth(2);
  muRatesOpenUnpack->SetLineColor(kOrange);
  muRatesOpenUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesOpenUnpack->GetYaxis()->SetTitle("Rate");
  muRatesOpenUnpack->SetMarkerStyle(23);
  muRatesOpenUnpack->SetMarkerColor(kOrange);
  muRatesOpenUnpack->Draw("E1HIST");
  muRatesOpenUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesOpenUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    muRatesOpenReEmu->SetTitle("");
    //  muRatesOpenReEmu->SetLineWidth(2);
    muRatesOpenReEmu->SetLineColor(kBlue);
    muRatesOpenReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muRatesOpenReEmu->GetYaxis()->SetTitle("Rate");
    muRatesOpenReEmu->SetMarkerStyle(20);
    muRatesOpenReEmu->SetMarkerColor(kBlue);
    muRatesOpenReEmu->Draw("same,E1HIST");
    muRatesOpenReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muRatesOpenReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }

  TLegend* leg1 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg1->SetFillColor(0);
  leg1->AddEntry(muRatesOpenUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg1->AddEntry(muRatesOpenReEmu, "ReEmu", "lp");
  }
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->Draw();
  leg1->Draw();
  n3.DrawLatex(0.4, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.4, 0.55, "Zero Bias, L1_DoubleMu0");

  c1->SaveAs(plotFolder + "ratesDiMuonOpenLeading.pdf");
  c1->SaveAs(plotFolder + "ratesDiMuonOpenLeading.png");

  TCanvas* c2 = new TCanvas;
  // c2->SetLogy();

  //  muRatesUnpack->SetLineWidth(2);
  muRatesUnpack->SetLineColor(kOrange);
  muRatesUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesUnpack->GetYaxis()->SetTitle("Rate");
  muRatesUnpack->SetMarkerStyle(23);
  muRatesUnpack->SetMarkerColor(kOrange);
  // muRatesUnpack->GetYaxis()->SetRangeUser(0, 1000);
  muRatesUnpack->Draw("E1HIST");
  muRatesUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    muRatesReEmu->SetTitle("");
    //  muRatesReEmu->SetLineWidth(2);
    muRatesReEmu->SetLineColor(kBlue);
    muRatesReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muRatesReEmu->GetYaxis()->SetTitle("Rate");
    muRatesReEmu->SetMarkerStyle(20);
    muRatesReEmu->SetMarkerColor(kBlue);
    // muRatesReEmu->GetYaxis()->SetRangeUser(0, 1000);
    muRatesReEmu->Draw("same,E1HIST");
    muRatesReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muRatesReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }

  TLegend* leg2 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg2->SetFillColor(0);
  leg2->AddEntry(muRatesUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg2->AddEntry(muRatesReEmu, "ReEmu", "lp");
  }
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw();
  leg2->Draw();
  n3.DrawLatex(0.4, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  std::ostringstream oss;
  oss << "ZeroBias, L1_DoubleMu_" << mu1cut << "_" << mu2cut;
  n4.DrawLatex(0.4, 0.55, oss.str().c_str());

  c2->SaveAs(plotFolder + "ratesDiMuonLeading.pdf");
  c2->SaveAs(plotFolder + "ratesDiMuonLeading.pdf");

  TCanvas* c3 = new TCanvas;

  //  singleMuRatesOpenUnpack->SetLineWidth(2);
  singleMuRatesOpenUnpack->SetLineColor(kOrange);
  singleMuRatesOpenUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  singleMuRatesOpenUnpack->GetYaxis()->SetTitle("Rate");
  singleMuRatesOpenUnpack->SetMarkerStyle(23);
  singleMuRatesOpenUnpack->SetMarkerColor(kOrange);
  singleMuRatesOpenUnpack->Draw("E1HIST");
  singleMuRatesOpenUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  singleMuRatesOpenUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    singleMuRatesOpenReEmu->SetTitle("");
    //  singleMuRatesOpenReEmu->SetLineWidth(2);
    singleMuRatesOpenReEmu->SetLineColor(kBlue);
    singleMuRatesOpenReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    singleMuRatesOpenReEmu->GetYaxis()->SetTitle("Rate");
    singleMuRatesOpenReEmu->SetMarkerStyle(20);
    singleMuRatesOpenReEmu->SetMarkerColor(kBlue);
    singleMuRatesOpenReEmu->Draw("same,E1HIST");
    singleMuRatesOpenReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    singleMuRatesOpenReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }

  TLegend* leg3 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg3->SetFillColor(0);
  leg3->AddEntry(singleMuRatesOpenUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg3->AddEntry(singleMuRatesOpenReEmu, "ReEmu", "lp");
  }
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->Draw();
  leg3->Draw();
  n3.DrawLatex(0.4, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.4, 0.55, "Zero Bias, L1_SingleMuOpen, q>4");

  c3->SaveAs(plotFolder + "ratesSingleMuonOpen.pdf");
  c3->SaveAs(plotFolder + "ratesSingleMuonOpen.png");

  TCanvas* c4 = new TCanvas;

  //  muRatesOpenTrailingUnpack->SetLineWidth(2);
  muRatesOpenTrailingUnpack->SetLineColor(kOrange);
  muRatesOpenTrailingUnpack->GetXaxis()->SetTitle("#eta (trailing #mu)");
  muRatesOpenTrailingUnpack->GetYaxis()->SetTitle("Rate");
  muRatesOpenTrailingUnpack->SetMarkerStyle(23);
  muRatesOpenTrailingUnpack->SetMarkerColor(kOrange);
  muRatesOpenTrailingUnpack->Draw("E1HIST");
  muRatesOpenTrailingUnpack->GetXaxis()->SetTitle("#eta (trailing #mu)");
  muRatesOpenTrailingUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    muRatesOpenTrailingReEmu->SetTitle("");
    //  muRatesOpenTrailingReEmu->SetLineWidth(2);
    muRatesOpenTrailingReEmu->SetLineColor(kBlue);
    muRatesOpenTrailingReEmu->GetXaxis()->SetTitle("#eta (trailing #mu)");
    muRatesOpenTrailingReEmu->GetYaxis()->SetTitle("Rate");
    muRatesOpenTrailingReEmu->SetMarkerStyle(20);
    muRatesOpenTrailingReEmu->SetMarkerColor(kBlue);
    muRatesOpenTrailingReEmu->Draw("same,E1HIST");
    muRatesOpenTrailingReEmu->GetXaxis()->SetTitle("#eta (trailing #mu)");
    muRatesOpenTrailingReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }

  TLegend* leg4 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg4->SetFillColor(0);
  leg4->AddEntry(muRatesOpenTrailingUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg4->AddEntry(muRatesOpenTrailingReEmu, "ReEmu", "lp");
  }
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->Draw();
  leg4->Draw();
  n3.DrawLatex(0.4, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.4, 0.55, "Zero Bias, L1_DoubleMu0");

  c4->SaveAs(plotFolder + "ratesDiMuonOpenTrailing.pdf");
  c4->SaveAs(plotFolder + "ratesDiMuonOpenTrailing.png");

  TCanvas* c5 = new TCanvas;

  muGhostRatesOpenUnpack->SetLineColor(kOrange);
  muGhostRatesOpenUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muGhostRatesOpenUnpack->GetYaxis()->SetTitle("Rate");
  muGhostRatesOpenUnpack->SetMarkerStyle(23);
  muGhostRatesOpenUnpack->SetMarkerColor(kOrange);
  muGhostRatesOpenUnpack->Draw("E1HIST");
  muGhostRatesOpenUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muGhostRatesOpenUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    muGhostRatesOpenReEmu->SetTitle("");
    muGhostRatesOpenReEmu->SetLineColor(kBlue);
    muGhostRatesOpenReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muGhostRatesOpenReEmu->GetYaxis()->SetTitle("Rate");
    muGhostRatesOpenReEmu->SetMarkerStyle(20);
    muGhostRatesOpenReEmu->SetMarkerColor(kBlue);
    muGhostRatesOpenReEmu->Draw("same,E1HIST");
    muGhostRatesOpenReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muGhostRatesOpenReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }

  TLegend* leg5 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg5->SetFillColor(0);
  leg5->AddEntry(muGhostRatesOpenUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg5->AddEntry(muGhostRatesOpenReEmu, "ReEmu", "lp");
  }
  leg5->SetBorderSize(0);
  leg5->SetFillStyle(0);
  leg5->Draw();
  leg5->Draw();
  n3.DrawLatex(0.3, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.3, 0.55, "Zero Bias, L1_DoubleMu0, ghost rate");

  c5->SaveAs(plotFolder + "ghostRatesDiMuonOpenLeading.pdf");
  c5->SaveAs(plotFolder + "ghostRatesDiMuonOpenLeading.png");

  TCanvas* c6 = new TCanvas;

  muGhostRatesUnpack->SetLineColor(kOrange);
  muGhostRatesUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muGhostRatesUnpack->GetYaxis()->SetTitle("Rate");
  muGhostRatesUnpack->SetMarkerStyle(23);
  muGhostRatesUnpack->SetMarkerColor(kOrange);
  muGhostRatesUnpack->Draw("E1HIST");
  muGhostRatesUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muGhostRatesUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    muGhostRatesReEmu->SetTitle("");
    muGhostRatesReEmu->SetLineColor(kBlue);
    muGhostRatesReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muGhostRatesReEmu->GetYaxis()->SetTitle("Rate");
    muGhostRatesReEmu->SetMarkerStyle(20);
    muGhostRatesReEmu->SetMarkerColor(kBlue);
    muGhostRatesReEmu->Draw("same,E1HIST");
    muGhostRatesReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
    muGhostRatesReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }

  TLegend* leg6 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg6->SetFillColor(0);
  leg6->AddEntry(muGhostRatesUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg6->AddEntry(muGhostRatesReEmu, "ReEmu", "lp");
  }
  leg6->SetBorderSize(0);
  leg6->SetFillStyle(0);
  leg6->Draw();
  leg6->Draw();
  n3.DrawLatex(0.4, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  oss << ", ghost rate";
  n4.DrawLatex(0.4, 0.55, oss.str().c_str());

  c6->SaveAs(plotFolder + "ghostRatesDiMuonLeading.pdf");
  c6->SaveAs(plotFolder + "ghostRatesDiMuonLeading.pdf");

  TCanvas* c7 = new TCanvas;

  muGhostRatesOpenTrailingUnpack->SetLineColor(kOrange);
  muGhostRatesOpenTrailingUnpack->GetXaxis()->SetTitle("#eta (trailing #mu)");
  muGhostRatesOpenTrailingUnpack->GetYaxis()->SetTitle("Rate");
  muGhostRatesOpenTrailingUnpack->SetMarkerStyle(23);
  muGhostRatesOpenTrailingUnpack->SetMarkerColor(kOrange);
  muGhostRatesOpenTrailingUnpack->Draw("E1HIST");
  muGhostRatesOpenTrailingUnpack->GetXaxis()->SetTitle("#eta (trailing #mu)");
  muGhostRatesOpenTrailingUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  if (drawReEmu) {
    muGhostRatesOpenTrailingReEmu->SetTitle("");
    muGhostRatesOpenTrailingReEmu->SetLineColor(kBlue);
    muGhostRatesOpenTrailingReEmu->GetXaxis()->SetTitle("#eta (trailing #mu)");
    muGhostRatesOpenTrailingReEmu->GetYaxis()->SetTitle("Rate");
    muGhostRatesOpenTrailingReEmu->SetMarkerStyle(20);
    muGhostRatesOpenTrailingReEmu->SetMarkerColor(kBlue);
    muGhostRatesOpenTrailingReEmu->Draw("same,E1HIST");
    muGhostRatesOpenTrailingReEmu->GetXaxis()->SetTitle("#eta (trailing #mu)");
    muGhostRatesOpenTrailingReEmu->GetYaxis()->SetTitle("Rate [kHz]");
    gPad->Modified();
  }
  TLegend* leg7 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg7->SetFillColor(0);
  leg7->AddEntry(muGhostRatesOpenTrailingUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg7->AddEntry(muGhostRatesOpenTrailingReEmu, "ReEmu", "lp");
  }

  leg7->SetBorderSize(0);
  leg7->SetFillStyle(0);
  leg7->Draw();
  leg7->Draw();
  n3.DrawLatex(0.3, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.3, 0.55, "Zero Bias, L1_DoubleMu0, ghost rate");

  c7->SaveAs(plotFolder + "ghostRatesDiMuonOpenTrailing.pdf");
  c7->SaveAs(plotFolder + "ghostRatesDiMuonOpenTrailing.png");

  TCanvas* c8 = new TCanvas;

  tfMuIdxUnpack->SetLineColor(kOrange);
  tfMuIdxUnpack->SetMarkerStyle(23);
  tfMuIdxUnpack->SetMarkerColor(kOrange);
  tfMuIdxUnpack->Draw("E1HIST");
  tfMuIdxUnpack->GetXaxis()->SetTitle("TF idx (leading #mu)");
  tfMuIdxUnpack->GetYaxis()->SetTitle("Counts");

  if (drawReEmu) {
    tfMuIdxReEmu->SetTitle("");
    tfMuIdxReEmu->SetLineColor(kBlue);
    tfMuIdxReEmu->SetMarkerStyle(20);
    tfMuIdxReEmu->SetMarkerColor(kBlue);
    tfMuIdxReEmu->Draw("same,E1HIST");
    tfMuIdxReEmu->GetXaxis()->SetTitle("TF idx (leading #mu)");
    tfMuIdxReEmu->GetYaxis()->SetTitle("Counts");
    gPad->Modified();
  }

  TLegend* leg8 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg8->SetFillColor(0);
  leg8->AddEntry(tfMuIdxUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg8->AddEntry(tfMuIdxReEmu, "ReEmu", "lp");
  }
  leg8->SetBorderSize(0);
  leg8->SetFillStyle(0);
  leg8->Draw();
  leg8->Draw();
  n3.DrawLatex(0.3, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.3, 0.55, "Zero Bias, L1_DoubleMu0");

  c8->SaveAs(plotFolder + "tfMuIdxLeading.pdf");
  c8->SaveAs(plotFolder + "tfMuIdxLeading.png");

  TCanvas* c9 = new TCanvas;

  tfMuIdxAllUnpack->SetLineColor(kOrange);
  tfMuIdxAllUnpack->SetMarkerStyle(23);
  tfMuIdxAllUnpack->SetMarkerColor(kOrange);
  tfMuIdxAllUnpack->Draw("E1HIST");
  tfMuIdxAllUnpack->GetXaxis()->SetTitle("TF idx (leading #mu)");
  tfMuIdxAllUnpack->GetYaxis()->SetTitle("Counts");

  if (drawReEmu) {
    tfMuIdxAllReEmu->SetTitle("");
    tfMuIdxAllReEmu->SetLineColor(kBlue);
    tfMuIdxAllReEmu->SetMarkerStyle(20);
    tfMuIdxAllReEmu->SetMarkerColor(kBlue);
    tfMuIdxAllReEmu->Draw("same,E1HIST");
    tfMuIdxAllReEmu->GetXaxis()->SetTitle("TF idx (leading #mu)");
    tfMuIdxAllReEmu->GetYaxis()->SetTitle("Counts");
    gPad->Modified();
  }
  TLegend* leg9 = new TLegend(0.4, 0.73, 0.6, 0.88);
  leg9->SetFillColor(0);
  leg9->AddEntry(tfMuIdxAllUnpack, "Unpack", "lp");
  if (drawReEmu) {
    leg9->AddEntry(tfMuIdxAllReEmu, "ReEmu", "lp");
  }
  leg9->SetBorderSize(0);
  leg9->SetFillStyle(0);
  leg9->Draw();
  leg9->Draw();
  n3.DrawLatex(0.3, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.3, 0.55, "Zero Bias");

  c9->SaveAs(plotFolder + "tfMuIdxAllLeading.pdf");
  c9->SaveAs(plotFolder + "tfMuIdxAllLeading.png");
}
