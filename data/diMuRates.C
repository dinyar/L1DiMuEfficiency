#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TChain.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"

// TODO: Make rate plots for uGMT inputs!

void diMuRates(const char * fname="file_list_275125", int mu1cut = 11, int mu2cut = 4){
  gStyle->SetOptStat(0);

  std::vector<std::string> listNtuples;

  // OpenNtupleList
  std::ifstream flist(fname);
  if (!flist) {
    std::cout << "File "<<fname<<" is not found!"<<std::endl;
    return;
  }

  while(!flist.eof())
  {
    std::string str;
    getline(flist,str);
    if (!flist.fail())
    {
      if (str!="") listNtuples.push_back(str);
    }
  }

  // CheckFirstFile
  if (listNtuples.size()==0) return;

  TFile* rf = TFile::Open(listNtuples[0].c_str());

  if (rf==0) return;
  if (rf->IsOpen()==0) return;

  std::string unpackTreepath("l1UpgradeTree/L1UpgradeTree");
  std::string reEmuTreepath("l1UpgradeEmuTree/L1UpgradeTree");
  TTree * treeL1Unpack  = (TTree*) rf->Get(unpackTreepath.c_str());
  TTree * treeL1reEmu  = (TTree*) rf->Get(reEmuTreepath.c_str());

  if (!treeL1reEmu || !treeL1Unpack) {
    std::cout<<"L1UpgradeTree not found.. "<<std::endl;
    return;
  }

  // OpenWithoutInit
  TChain* l1UnpackChain = new TChain(unpackTreepath.c_str());
  TChain* l1reEmuChain = new TChain(reEmuTreepath.c_str());
  for (unsigned int i=0;i<listNtuples.size();i++)
  {
    std::cout << " -- Adding " << listNtuples[i] << std::endl;
    l1UnpackChain->Add(listNtuples[i].c_str());
    l1reEmuChain->Add(listNtuples[i].c_str());
  }

  // Init
  std::cout << "Estimate the number of entries ..."<<std::endl;
  int nentries = l1UnpackChain->GetEntries();
  std::cout << nentries << std::endl;
  int nevents = nentries;

  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat *unpack_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisL1UpgradeDataFormat *reEmu_ = new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  l1UnpackChain->SetBranchAddress("L1Upgrade", &unpack_);
  l1reEmuChain->SetBranchAddress("L1Upgrade", &reEmu_);

  // mu bins
  int nMuBins = 25;
  float muLo = -2.5;
  float muHi = 2.5;
  float muBinWidth = (muHi-muLo)/nMuBins;

  //make histos
  TH1D* singleMuRatesOpenUnpack = new TH1D("singleMuRatesOpenUnpack", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* singleMuRatesOpenReEmu = new TH1D("singleMuRatesOpenReEmu", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* muRatesUnpack = new TH1D("muRatesUnpack", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* muRatesReEmu = new TH1D("muRatesReEmu", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* muRatesOpenUnpack = new TH1D("muRatesOpenUnpack", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* muRatesOpenReEmu = new TH1D("muRatesOpenReEmu", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* muRatesOpenTrailingUnpack = new TH1D("muRatesOpenTrailingUnpack", "", nMuBins, muLo-0.1, muHi+0.1);
  TH1D* muRatesOpenTrailingReEmu = new TH1D("muRatesOpenTrailingReEmu", "", nMuBins, muLo-0.1, muHi+0.1);
  int unpackCounts(0);
  int reEmuCounts(0);

  for (Long64_t jentry=0; jentry<nevents;jentry++){

    if((jentry%1000)==0) std::cout << "Done " << jentry  << " events..." << std::endl;

    l1UnpackChain->GetEntry(jentry);
    l1reEmuChain->GetEntry(jentry);

    // get Mu rates
    int mu1Unpack(-1);
    int mu2Unpack(-1);
    double mu1PtUnpack(0);
    double mu2PtUnpack(0);
    int mu1ReEmu(-1);
    int mu2ReEmu(-1);
    double mu1PtReEmu(0);
    double mu2PtReEmu(0);
    for(uint it=0; it<unpack_->nMuons; ++it){
      if (unpack_->muonQual[it]<=4 || unpack_->muonBx[it] != 0) { 
        continue;
      }
      if (unpack_->muonEt[it] > mu1PtUnpack) {
        mu2Unpack = mu1Unpack;
        mu2PtUnpack = mu1PtUnpack;
        mu1Unpack = it;
        mu1PtUnpack = unpack_->muonEt[it];
      } else if(unpack_->muonEt[it] > mu2PtUnpack) {
        mu2Unpack = it;
        mu2PtUnpack = unpack_->muonEt[it];
      }
    }
    for(uint it=0; it<reEmu_->nMuons; ++it){
      if (reEmu_->muonQual[it]<=4 || reEmu_->muonBx[it] != 0) { 
        continue;
      }
      if (reEmu_->muonEt[it] > mu1PtReEmu) {
        mu2ReEmu = mu1ReEmu;
        mu2PtReEmu = mu1PtReEmu;
        mu1ReEmu = it;
        mu1PtReEmu = reEmu_->muonEt[it];
      } else if(reEmu_->muonEt[it] > mu2PtReEmu) {
        mu2ReEmu = it;
        mu2PtReEmu = reEmu_->muonEt[it];
      }
    }
    if(mu1Unpack != -1) {
      singleMuRatesOpenUnpack->Fill(unpack_->muonEta[mu1Unpack]);
    }
    if(mu1ReEmu != -1) {
      singleMuRatesOpenReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
    }

    if(mu1Unpack != -1 && mu2Unpack != -1) {
      if(mu1PtUnpack >= mu1cut && mu2PtUnpack >= mu2cut) {
        muRatesUnpack->Fill(unpack_->muonEta[mu1Unpack]);
      }
      muRatesOpenUnpack->Fill(unpack_->muonEta[mu1Unpack]);
      muRatesOpenTrailingUnpack->Fill(unpack_->muonEta[mu2Unpack]);
      ++unpackCounts;
    }
    if(mu1ReEmu != -1 && mu2ReEmu != -1) {
      if(mu1PtReEmu >= mu1cut && mu2PtReEmu >= mu2cut) {
        muRatesReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
      }
      muRatesOpenReEmu->Fill(reEmu_->muonEta[mu1ReEmu]);
      muRatesOpenTrailingReEmu->Fill(reEmu_->muonEta[mu2ReEmu]);
      ++reEmuCounts;
    }
  }


  //normalisation factor
  double norm = (11.*2028.)/nevents; // zb rate = n_colliding * 11 kHz 
  std::cout << "norm = " << norm << std::endl;

  std::cout << "###########################################" << std::endl;
  std::cout << "** Computed rates: **" << std::endl;
  std::cout << "Unpacked rate: " << muRatesUnpack->GetEntries()*norm << std::endl;
  std::cout << "Unpacked open rate: " << muRatesOpenUnpack->GetEntries()*norm << std::endl;
  std::cout << "ReEmulated rate: " << muRatesReEmu->GetEntries()*norm << std::endl;
  std::cout << "ReEmulated open rate: " << muRatesOpenReEmu->GetEntries()*norm << std::endl;
  std::cout << "Unpacked SingleMuOpen rate: " << singleMuRatesOpenUnpack->GetEntries()*norm << std::endl;
  std::cout << "ReEmulated SingleMuOpen rate: " << singleMuRatesOpenReEmu->GetEntries()*norm << std::endl;
  std::cout << "###########################################" << std::endl;

  muRatesUnpack->Sumw2();
  muRatesReEmu->Sumw2();
  muRatesOpenUnpack->Sumw2();
  muRatesOpenReEmu->Sumw2();
  TF1 *constant = new TF1("constant", "1", -5, 5);
  muRatesUnpack->Multiply(constant, norm);
  muRatesReEmu->Multiply(constant, norm);
  muRatesOpenUnpack->Multiply(constant, norm);
  muRatesOpenReEmu->Multiply(constant, norm);

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


  TLegend* leg1 = new TLegend(0.4,0.73,0.6,0.88);
  leg1->SetFillColor(0);
  leg1->AddEntry(muRatesOpenUnpack,"Unpack","lp");
  leg1->AddEntry(muRatesOpenReEmu,"ReEmu","lp");
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->Draw();
  leg1->Draw();
  n3.DrawLatex(0.4, 0.6, "Run 275125 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.4, 0.55, "Zero Bias, L1_DoubleMu0");

  c1->SaveAs("ratesDiMuonOpenLeading.pdf");
  c1->SaveAs("ratesDiMuonOpenLeading.png");

  TCanvas* c2 = new TCanvas;
  //c2->SetLogy();

//  muRatesUnpack->SetLineWidth(2);
  muRatesUnpack->SetLineColor(kOrange);
  muRatesUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesUnpack->GetYaxis()->SetTitle("Rate");
  muRatesUnpack->SetMarkerStyle(23);
  muRatesUnpack->SetMarkerColor(kOrange);
  //muRatesUnpack->GetYaxis()->SetRangeUser(0, 1000);
  muRatesUnpack->Draw("E1HIST");
  muRatesUnpack->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesUnpack->GetYaxis()->SetTitle("Rate [kHz]");

  muRatesReEmu->SetTitle("");
//  muRatesReEmu->SetLineWidth(2);
  muRatesReEmu->SetLineColor(kBlue);
  muRatesReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesReEmu->GetYaxis()->SetTitle("Rate");
  muRatesReEmu->SetMarkerStyle(20);
  muRatesReEmu->SetMarkerColor(kBlue);
  //muRatesReEmu->GetYaxis()->SetRangeUser(0, 1000);
  muRatesReEmu->Draw("same,E1HIST");
  muRatesReEmu->GetXaxis()->SetTitle("#eta (leading #mu)");
  muRatesReEmu->GetYaxis()->SetTitle("Rate [kHz]");
  gPad->Modified();


  TLegend* leg2 = new TLegend(0.4,0.73,0.6,0.88);
  leg2->SetFillColor(0);
  leg2->AddEntry(muRatesUnpack,"Unpack","lp");
  leg2->AddEntry(muRatesReEmu,"ReEmu","lp");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw();
  leg2->Draw();
  n3.DrawLatex(0.4, 0.6, "Run 275125 #sqrt{s} = 13 TeV");
  std::ostringstream oss;
  oss << "ZeroBias, L1_DoubleMu_" << mu1cut << "_" << mu2cut;
  n4.DrawLatex(0.4, 0.55, oss.str().c_str());

  c2->SaveAs("ratesDiMuonLeading.pdf");
  c2->SaveAs("ratesDiMuonLeading.pdf");

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


  TLegend* leg3 = new TLegend(0.4,0.73,0.6,0.88);
  leg3->SetFillColor(0);
  leg3->AddEntry(singleMuRatesOpenUnpack,"Unpack","lp");
  leg3->AddEntry(singleMuRatesOpenReEmu,"ReEmu","lp");
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->Draw();
  leg3->Draw();
  n3.DrawLatex(0.4, 0.6, "Run 275125 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.4, 0.55, "Zero Bias, L1_SingleMuOpen, q>4");

  c3->SaveAs("ratesSingleMuonOpen.pdf");
  c3->SaveAs("ratesSingleMuonOpen.png");

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


  TLegend* leg4 = new TLegend(0.4,0.73,0.6,0.88);
  leg4->SetFillColor(0);
  leg4->AddEntry(muRatesOpenTrailingUnpack,"Unpack","lp");
  leg4->AddEntry(muRatesOpenTrailingReEmu,"ReEmu","lp");
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->Draw();
  leg4->Draw();
  n3.DrawLatex(0.4, 0.6, "Run 275125 #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.4, 0.55, "Zero Bias, L1_DoubleMu0");

  c4->SaveAs("ratesDiMuonOpenTrailing.pdf");
  c4->SaveAs("ratesDiMuonOpenTrailing.png");
}
