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

#include <iostream>
#include "TString.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeTfMuonDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

// TODO: Only check muons that are in neighbouring sectors
// TODO: Only check muons that are in neighbouring sector/wedges

double calcGlobalPhi(int locPhi, int tfType, int proc) {
  int globPhi = 0;
  if (tfType == 0) {
    // each BMTF processor corresponds to a 30 degree wedge = 48 in int-scale
    globPhi = (proc) * 48 + locPhi;
    // first processor starts at CMS phi = -15 degrees...
    globPhi += 576-24;
    // handle wrap-around (since we add the 576-24, the value will never be negative!)
    globPhi = globPhi%576;
  } else {
    // all others correspond to 60 degree sectors = 96 in int-scale
    globPhi = (proc) * 96 + locPhi;
    // first processor starts at CMS phi = 15 degrees (24 in int)... Handle wrap-around with %. Add 576 to make sure the number is positive
    globPhi = (globPhi + 600) % 576;
  }
  return globPhi;
}

void ghostDistances(const char * fname="L1Ntuple_list", TString run="XXXX"){

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // make trees and set branch addresses
  const char * ugmtUnpackTree="l1UpgradeTree/L1UpgradeTree";
  const char * ugmtMcTree="l1UpgradeEmuTree/L1UpgradeTree";
  const char * tfUnpackTree="l1UpgradeTfMuonTree/L1UpgradeTfMuonTree";
  const char * tfMcTree="l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree";
  const char * recoTree="l1MuonRecoTree/Muon2RecoTree";

  std::vector<std::string> listNtuples;

  // OpenNtupleList
  std::ifstream flist(fname);
  if (!flist) {
    std::cout << "File "<<fname<<" is not found!"<<std::endl;
    return false;
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
  if (listNtuples.size()==0) return false;

  TFile* file = TFile::Open(listNtuples[0].c_str());

  if (file==0) return false;
  if (file->IsOpen()==0) return false;

  bool foundTfUnpackTree = false;
  bool foundTfMcTree = false;
  bool foundReco = false;
  TTree * treeTfUnpack  = (TTree*) file->Get(tfUnpackTree);
  TTree * treeTfReEmu  = (TTree*) file->Get(tfMcTree);
  TTree * treeReco  = (TTree*) file->Get(recoTree);

  if (!(treeTfUnpack || treeTfReEmu) || !treeReco) {
    std::cout<<"Didn't find TF upgrade or reco tree, exiting.. "<<std::endl;
    return;
  }

  // OpenWithoutInit
  TChain* chainTfUnpack = new TChain(tfUnpackTree);
  TChain* chainTfReEmu = new TChain(tfMcTree);
  TChain* chainReco = new TChain(recoTree);
  for (unsigned int i=0;i<listNtuples.size();i++)
  {
          std::cout << " -- Adding " << listNtuples[i] << std::endl;
          chainTfUnpack->Add(listNtuples[i].c_str());
          chainTfReEmu->Add(listNtuples[i].c_str());
          chainReco->Add(listNtuples[i].c_str());
  }

  // Init
  std::cout << "Estimate the number of entries..."<<std::endl;
  int nevents = chainReco->GetEntries();
  std::cout << nevents << std::endl;

  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat *bmtfUnpack_ = new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat *bmtfReEmu_ = new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat *omtfUnpack_ = new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat *omtfReEmu_ = new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat *emtfUnpack_ = new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat *emtfReEmu_ = new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat *reco_ = new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  if (treeTfUnpack){
    std::cout << "Found TF upgrade tree from data, adding branches.." << std::endl;
    foundTfUnpackTree = true;
    chainTfUnpack->SetBranchAddress("L1UpgradeBmtfMuon", &bmtfUnpack_);
    chainTfUnpack->SetBranchAddress("L1UpgradeOmtfMuon", &omtfUnpack_);
    chainTfUnpack->SetBranchAddress("L1UpgradeEmtfMuon", &emtfUnpack_);
  }
  if (treeTfReEmu){
    std::cout << "Found TF upgrade tree from ReEmu, adding branches.." << std::endl;
    foundTfMcTree = true;
    chainTfReEmu->SetBranchAddress("L1UpgradeBmtfMuon", &bmtfReEmu_);
    chainTfReEmu->SetBranchAddress("L1UpgradeOmtfMuon", &omtfReEmu_);
    chainTfReEmu->SetBranchAddress("L1UpgradeEmtfMuon", &emtfReEmu_);
  }
  if (treeReco){
    std::cout << "Found reco tree, adding branches.." << std::endl;
    foundReco = true;
    chainReco->SetBranchAddress("Muon", &reco_);
  }

  // mu bins
  int nMuPtBins = 50;
  float muPtLo = 10.;
  float muPtHi = 50.;
  float muPtBinWidth = (muPtHi-muPtLo)/nMuPtBins;
  int nMuDetaBins = 32;
  float muDetaLo = -0.4;
  float muDetaHi = 0.4;
  float muDetaBinWidth = (muDetaHi-muDetaLo)/nMuDetaBins;
  int nMuDphiBins = 32;
  float muDphiLo = -0.4;
  float muDphiHi = 0.4;
  float muDphiBinWidth = (muDphiHi-muDphiLo)/nMuDphiBins;
  int nMuDrBins = 32;
  float muDrLo = -0.4;
  float muDrHi = 0.4;
  float muDrBinWidth = (muDrHi-muDrLo)/nMuDrBins;

  //make histos
  TH1D* resolutionUnpackEtaBMTF = new TH1D("resolutionUnpackEtaBMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* resolutionUnpackEtaOMTF = new TH1D("resolutionUnpackEtaOMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* resolutionUnpackEtaEMTF = new TH1D("resolutionUnpackEtaEMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* resolutionReEmuEtaBMTF = new TH1D("resolutionReEmuEtaBMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* resolutionReEmuEtaOMTF = new TH1D("resolutionReEmuEtaOMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* resolutionReEmuEtaEMTF = new TH1D("resolutionReEmuEtaEMTF", "", nMuDetaBins, muDetaLo, muDetaHi);

  TH1D* ghostDiffUnpackEtaBOMTFFine = new TH1D("ghostDiffUnpackEtaBOMTFFine", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffUnpackEtaBOMTFCoarse = new TH1D("ghostDiffUnpackEtaBOMTFCoarse", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffUnpackEtaEOMTF = new TH1D("ghostDiffUnpackEtaEOMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffUnpackEtaOMTF = new TH1D("ghostDiffUnpackEtaOMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffUnpackEtaEMTF = new TH1D("ghostDiffUnpackEtaEMTF", "", nMuDetaBins, muDetaLo, muDetaHi);

  TH1D* ghostDiffUnpackPhiBOMTFFine = new TH1D("ghostDiffUnpackPhiBOMTFFine", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffUnpackPhiBOMTFCoarse = new TH1D("ghostDiffUnpackPhiBOMTFCoarse", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffUnpackPhiEOMTF = new TH1D("ghostDiffUnpackPhiEOMTF", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffUnpackPhiOMTF = new TH1D("ghostDiffUnpackPhiOMTF", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffUnpackPhiEMTF = new TH1D("ghostDiffUnpackPhiEMTF", "", nMuDphiBins, muDphiLo, muDphiHi);

  TH1D* ghostDiffUnpackRBOMTFFine = new TH1D("ghostDiffUnpackRBOMTFFine", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffUnpackRBOMTFCoarse = new TH1D("ghostDiffUnpackRBOMTFCoarse", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffUnpackREOMTF = new TH1D("ghostDiffUnpackREOMTF", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffUnpackROMTF = new TH1D("ghostDiffUnpackROMTF", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffUnpackREMTF = new TH1D("ghostDiffUnpackREMTF", "", nMuDrBins, muDrLo, muDrHi);

  TH1D* ghostDiffReEmuEtaBOMTFFine = new TH1D("ghostDiffReEmuEtaBOMTFFine", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffReEmuEtaBOMTFCoarse = new TH1D("ghostDiffReEmuEtaBOMTFCoarse", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffReEmuEtaEOMTF = new TH1D("ghostDiffReEmuEtaEOMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffReEmuEtaOMTF = new TH1D("ghostDiffReEmuEtaOMTF", "", nMuDetaBins, muDetaLo, muDetaHi);
  TH1D* ghostDiffReEmuEtaEMTF = new TH1D("ghostDiffReEmuEtaEMTF", "", nMuDetaBins, muDetaLo, muDetaHi);

  TH1D* ghostDiffReEmuPhiBOMTFFine = new TH1D("ghostDiffReEmuPhiBOMTFFine", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffReEmuPhiBOMTFCoarse = new TH1D("ghostDiffReEmuPhiBOMTFCoarse", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffReEmuPhiEOMTF = new TH1D("ghostDiffReEmuPhiEOMTF", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffReEmuPhiOMTF = new TH1D("ghostDiffReEmuPhiOMTF", "", nMuDphiBins, muDphiLo, muDphiHi);
  TH1D* ghostDiffReEmuPhiEMTF = new TH1D("ghostDiffReEmuPhiEMTF", "", nMuDphiBins, muDphiLo, muDphiHi);

  TH1D* ghostDiffReEmuRBOMTFFine = new TH1D("ghostDiffReEmuRBOMTFFine", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffReEmuRBOMTFCoarse = new TH1D("ghostDiffReEmuRBOMTFCoarse", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffReEmuREOMTF = new TH1D("ghostDiffReEmuREOMTF", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffReEmuROMTF = new TH1D("ghostDiffReEmuROMTF", "", nMuDrBins, muDrLo, muDrHi);
  TH1D* ghostDiffReEmuREMTF = new TH1D("ghostDiffReEmuREMTF", "", nMuDrBins, muDrLo, muDrHi);

  unsigned nEMTFunpack = 0;
  unsigned nOMTFunpack = 0;
  unsigned nBOMTFFineUnpack = 0;
  unsigned nBOMTFCoarseUnpack = 0;
  unsigned nEOMTFunpack = 0;
  unsigned nEMTFreemu = 0;
  unsigned nOMTFreemu = 0;
  unsigned nBOMTFFineReemu = 0;
  unsigned nBOMTFCoarseReemu = 0;
  unsigned nEOMTFreemu = 0;
  
  // Plot distance from ghosts
  std::cout << "Running over " << nevents << std::endl;

  for (Long64_t jentry=0; jentry<nevents;jentry++){

    if((jentry%1000)==0) std::cout << "Done " << jentry  << " events..." << std::endl;

    chainTfUnpack->GetEntry(jentry);
    chainReco->GetEntry(jentry);
    
    // Check for only one reco muon
    if(reco_->nMuons != 1) {
      continue;
    }
    if(omtfUnpack_->nTfMuons == 1) {
      double eta = 0.010875 * omtfUnpack_->tfMuonHwEta[0];
      double dEta = eta - reco_->eta[0];
      resolutionUnpackEtaOMTF->Fill(dEta);
    }
    if(bmtfUnpack_->nTfMuons == 1) {
      double eta = 0.010875 * bmtfUnpack_->tfMuonHwEta[0];
      double dEta = eta - reco_->eta[0];
      resolutionUnpackEtaBMTF->Fill(dEta);
    }
    if(emtfUnpack_->nTfMuons == 1) {
      double eta = 0.010875 * emtfUnpack_->tfMuonHwEta[0];
      double dEta = eta - reco_->eta[0];
      resolutionUnpackEtaEMTF->Fill(dEta);
    }
    if(omtfUnpack_->nTfMuons == 2) {
      if(!(omtfUnpack_->tfMuonHwQual[0] > 4 && omtfUnpack_->tfMuonHwQual[1] > 4)) {
        continue;
      }
      ++nOMTFunpack;

      double phi1 = 0.010908 * calcGlobalPhi(omtfUnpack_->tfMuonHwPhi[0], 1, omtfUnpack_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(omtfUnpack_->tfMuonHwPhi[1], 1, omtfUnpack_->tfMuonProcessor[1]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * omtfUnpack_->tfMuonHwEta[0];
      double eta2 = 0.010875 * omtfUnpack_->tfMuonHwEta[1];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diOMTF event." << std::endl;
        std::cout << "eta1: " << omtfUnpack_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta2: " << omtfUnpack_->tfMuonHwEta[1] << std::endl;
      }

      ghostDiffUnpackEtaOMTF->Fill(dEta);
      ghostDiffUnpackPhiOMTF->Fill(dPhi);
      ghostDiffUnpackROMTF->Fill(dR);
    }
    if(emtfUnpack_->nTfMuons == 2) {
      if(!(emtfUnpack_->tfMuonHwQual[0] > 4 && emtfUnpack_->tfMuonHwQual[1] > 4)) {
        continue;
      }
      ++nEMTFunpack;

      double phi1 = 0.010908 * calcGlobalPhi(emtfUnpack_->tfMuonHwPhi[0], 2, emtfUnpack_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(emtfUnpack_->tfMuonHwPhi[1], 2, emtfUnpack_->tfMuonProcessor[1]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * emtfUnpack_->tfMuonHwEta[0];
      double eta2 = 0.010875 * emtfUnpack_->tfMuonHwEta[1];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diEMTF event." << std::endl;
        std::cout << "eta1: " << emtfUnpack_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta2: " << emtfUnpack_->tfMuonHwEta[1] << std::endl;
      }

      ghostDiffUnpackEtaEMTF->Fill(dEta);
      ghostDiffUnpackPhiEMTF->Fill(dPhi);
      ghostDiffUnpackREMTF->Fill(dR);
    }
    if(bmtfUnpack_->nTfMuons == 1 && omtfUnpack_->nTfMuons == 1) {
      if(!(bmtfUnpack_->tfMuonHwQual[0] > 4 && omtfUnpack_->tfMuonHwQual[0] > 4)) {
        continue;
      }

      double phi1 = 0.010908 * calcGlobalPhi(bmtfUnpack_->tfMuonHwPhi[0], 0, bmtfUnpack_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(omtfUnpack_->tfMuonHwPhi[0], 1, omtfUnpack_->tfMuonProcessor[0]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * bmtfUnpack_->tfMuonHwEta[0];
      double eta2 = 0.010875 * omtfUnpack_->tfMuonHwEta[0];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diBOMTF event." << std::endl;
        std::cout << "eta bmtf: " << bmtfUnpack_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta omtf: " << omtfUnpack_->tfMuonHwEta[0] << std::endl;
      }

      if(bmtfUnpack_->tfMuonHwHF[0] == 1) {
        ++nBOMTFFineUnpack;

        ghostDiffUnpackEtaBOMTFFine->Fill(dEta);
        ghostDiffUnpackPhiBOMTFFine->Fill(dPhi);
        ghostDiffUnpackRBOMTFFine->Fill(dR);
      } else {
        ++nBOMTFCoarseUnpack;

        ghostDiffUnpackEtaBOMTFCoarse->Fill(dEta);
        ghostDiffUnpackPhiBOMTFCoarse->Fill(dPhi);
        ghostDiffUnpackRBOMTFCoarse->Fill(dR);
      }
    }
    if(emtfUnpack_->nTfMuons == 1 && omtfUnpack_->nTfMuons == 1) {
      if(!(emtfUnpack_->tfMuonHwQual[0] > 4 && omtfUnpack_->tfMuonHwQual[0] > 4)) {
        continue;
      }
      ++nEOMTFunpack;

      double phi1 = 0.010908 * calcGlobalPhi(emtfUnpack_->tfMuonHwPhi[0], 2, emtfUnpack_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(omtfUnpack_->tfMuonHwPhi[0], 1, omtfUnpack_->tfMuonProcessor[0]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * emtfUnpack_->tfMuonHwEta[0];
      double eta2 = 0.010875 * omtfUnpack_->tfMuonHwEta[0];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diEOMTF event." << std::endl;
        std::cout << "eta emtf: " << emtfUnpack_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta omtf: " << omtfUnpack_->tfMuonHwEta[0] << std::endl;
      }

      ghostDiffUnpackEtaEOMTF->Fill(dEta);
      ghostDiffUnpackPhiEOMTF->Fill(dPhi);
      ghostDiffUnpackREOMTF->Fill(dR);
    }

    chainTfReEmu->GetEntry(jentry);
    
    if(omtfReEmu_->nTfMuons == 1) {
      double eta = 0.010875 * omtfReEmu_->tfMuonHwEta[0];
      double dEta = eta - reco_->eta[0];
      resolutionReEmuEtaOMTF->Fill(dEta);
    }
    if(bmtfReEmu_->nTfMuons == 1) {
      double eta = 0.010875 * bmtfReEmu_->tfMuonHwEta[0];
      double dEta = eta - reco_->eta[0];
      resolutionReEmuEtaBMTF->Fill(dEta);
    }
    if(emtfReEmu_->nTfMuons == 1) {
      double eta = 0.010875 * emtfReEmu_->tfMuonHwEta[0];
      double dEta = eta - reco_->eta[0];
      resolutionReEmuEtaEMTF->Fill(dEta);
    }
    if(omtfReEmu_->nTfMuons == 2) {
      if(!(omtfReEmu_->tfMuonHwQual[0] > 4 && omtfReEmu_->tfMuonHwQual[1] > 4)) {
        continue;
      }
      ++nOMTFreemu;

      double phi1 = 0.010908 * calcGlobalPhi(omtfReEmu_->tfMuonHwPhi[0], 1, omtfReEmu_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(omtfReEmu_->tfMuonHwPhi[1], 1, omtfReEmu_->tfMuonProcessor[1]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * omtfReEmu_->tfMuonHwEta[0];
      double eta2 = 0.010875 * omtfReEmu_->tfMuonHwEta[1];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diOMTF reemu event." << std::endl;
        std::cout << "eta1: " << omtfReEmu_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta2: " << omtfReEmu_->tfMuonHwEta[1] << std::endl;
      }

      ghostDiffReEmuEtaOMTF->Fill(dEta);
      ghostDiffReEmuPhiOMTF->Fill(dPhi);
      ghostDiffReEmuROMTF->Fill(dR);
    }
    if(emtfReEmu_->nTfMuons == 2) {
      if(!(emtfReEmu_->tfMuonHwQual[0] > 4 && emtfReEmu_->tfMuonHwQual[1] > 4)) {
        continue;
      }
      ++nEMTFreemu;

      double phi1 = 0.010908 * calcGlobalPhi(emtfReEmu_->tfMuonHwPhi[0], 2, emtfReEmu_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(emtfReEmu_->tfMuonHwPhi[1], 2, emtfReEmu_->tfMuonProcessor[1]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * emtfReEmu_->tfMuonHwEta[0];
      double eta2 = 0.010875 * emtfReEmu_->tfMuonHwEta[1];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diEMTF reemu event." << std::endl;
        std::cout << "eta1: " << emtfReEmu_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta2: " << emtfReEmu_->tfMuonHwEta[1] << std::endl;
      }

      ghostDiffReEmuEtaEMTF->Fill(dEta);
      ghostDiffReEmuPhiEMTF->Fill(dPhi);
      ghostDiffReEmuREMTF->Fill(dR);
    }
    if(bmtfReEmu_->nTfMuons == 1 && omtfReEmu_->nTfMuons == 1) {
      if(!(bmtfReEmu_->tfMuonHwQual[0] > 4 && omtfReEmu_->tfMuonHwQual[0] > 4)) {
        continue;
      }

      double phi1 = 0.010908 * calcGlobalPhi(bmtfReEmu_->tfMuonHwPhi[0], 0, bmtfReEmu_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(omtfReEmu_->tfMuonHwPhi[0], 1, omtfReEmu_->tfMuonProcessor[0]);
      double dPhi = phi1-phi2;
      double eta1 = 0.010875 * bmtfReEmu_->tfMuonHwEta[0];
      double eta2 = 0.010875 * omtfReEmu_->tfMuonHwEta[0];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diBOMTF reemu event." << std::endl;
        std::cout << "eta1 bmtf: " << bmtfReEmu_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta2 omtf: " << omtfReEmu_->tfMuonHwEta[0] << std::endl;
      }

      if(bmtfReEmu_->tfMuonHwHF[0] == 1) {
        ++nBOMTFFineReemu;

        ghostDiffReEmuEtaBOMTFFine->Fill(dEta);
        ghostDiffReEmuPhiBOMTFFine->Fill(dPhi);
        ghostDiffReEmuRBOMTFFine->Fill(dR);
      } else {
        ++nBOMTFCoarseReemu;

        ghostDiffReEmuEtaBOMTFCoarse->Fill(dEta);
        ghostDiffReEmuPhiBOMTFCoarse->Fill(dPhi);
        ghostDiffReEmuRBOMTFCoarse->Fill(dR);
      }
    }
    if(emtfReEmu_->nTfMuons == 1 && omtfReEmu_->nTfMuons == 1) {
      if(!(emtfReEmu_->tfMuonHwQual[0] > 4 && omtfReEmu_->tfMuonHwQual[0] > 4)) {
        continue;
      }
      ++nEOMTFreemu;

      double phi1 = 0.010908 * calcGlobalPhi(emtfReEmu_->tfMuonHwPhi[0], 2, emtfReEmu_->tfMuonProcessor[0]);
      double phi2 = 0.010908 * calcGlobalPhi(omtfReEmu_->tfMuonHwPhi[0], 1, omtfReEmu_->tfMuonProcessor[0]);
      double dPhi = phi1-phi2;
      double eta1 =  0.010875 * emtfReEmu_->tfMuonHwEta[0];
      double eta2 =  0.010875 * omtfReEmu_->tfMuonHwEta[0];
      double dEta = eta1-eta2;
      double dR = std::sqrt(TMath::Power(dEta, 2) + TMath::Power(dPhi, 2));

      if (dEta > 0.5) {
        std::cout << "dEta > 0.5 in diEOMTF reemu event." << std::endl;
        std::cout << "eta1 emtf: " << emtfReEmu_->tfMuonHwEta[0] << std::endl;
        std::cout << "eta2 omtf: " << omtfReEmu_->tfMuonHwEta[0] << std::endl;
      }

      ghostDiffReEmuEtaEOMTF->Fill(dEta);
      ghostDiffReEmuPhiEOMTF->Fill(dPhi);
      ghostDiffReEmuREOMTF->Fill(dR);
    }
  }



  // Report counts
  std::cout << "###############################################" << std::endl;
  std::cout << "Number of dimuons found:" << std::endl;
  std::cout << "###############################################" << std::endl;
  std::cout << std::endl;
  std::cout << "OMTF:" << std::endl;
  std::cout << "  - Unpack: " << nOMTFunpack << std::endl;
  std::cout << "  - ReEmu:  " << nOMTFreemu << std::endl;
  std::cout << "EMTF:" << std::endl;
  std::cout << "  - Unpack: " << nEMTFunpack << std::endl;
  std::cout << "  - ReEmu:  " << nEMTFreemu << std::endl;
  std::cout << "BMTF fine/OMTF:" << std::endl;
  std::cout << "  - Unpack: " << nBOMTFFineUnpack << std::endl;
  std::cout << "  - ReEmu:  " << nBOMTFFineReemu << std::endl;
  std::cout << "BMTF coarse/OMTF:" << std::endl;
  std::cout << "  - Unpack: " << nBOMTFCoarseUnpack << std::endl;
  std::cout << "  - ReEmu:  " << nBOMTFCoarseReemu << std::endl;
  std::cout << "OMTF/EMTF:" << std::endl;
  std::cout << "  - Unpack: " << nEOMTFunpack << std::endl;
  std::cout << "  - ReEmu:  " << nEOMTFreemu << std::endl;

  TLatex n1;
  n1.SetNDC();
  n1.SetTextFont(42);
  n1.SetTextSize(0.04);
   
  TLatex n3;
  n3.SetNDC();
  n3.SetTextFont(52);
  n3.SetTextSize(0.04);

  TLatex n4;
  n4.SetNDC();
  n4.SetTextFont(52);
  n4.SetTextSize(0.04);

  TString plotFolder = "plots/" + run + "/";

  // dEta
  
  TCanvas* c1 = new TCanvas;

  ghostDiffReEmuEtaEOMTF->SetTitle("");
  ghostDiffReEmuEtaEOMTF->SetLineWidth(2);
  ghostDiffReEmuEtaEOMTF->SetLineColor(kOrange);
  ghostDiffReEmuEtaEOMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuEtaEOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuEtaEOMTF->SetMarkerStyle(23);
  ghostDiffReEmuEtaEOMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuEtaEOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackEtaEOMTF->SetLineWidth(2);
  ghostDiffUnpackEtaEOMTF->SetLineColor(kBlue);
  ghostDiffUnpackEtaEOMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackEtaEOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackEtaEOMTF->SetMarkerStyle(20);
  ghostDiffUnpackEtaEOMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuEtaEOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackEtaEOMTF->Draw("E1HIST");
  ghostDiffReEmuEtaEOMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg1 = new TLegend(0.15,0.73,0.7,0.88);
  leg1->SetFillColor(0);
  leg1->AddEntry(ghostDiffReEmuEtaEOMTF, "ReEmu","lp");
  leg1->AddEntry(ghostDiffUnpackEtaEOMTF, "Unpack","lp");
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->Draw();
  leg1->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF/EMTF");

  c1->SaveAs(plotFolder + "dEtaEOMTF_" + run + ".pdf");
  c1->SaveAs(plotFolder + "dEtaEOMTF_" + run + ".png");


  TCanvas* c2 = new TCanvas;

  ghostDiffReEmuEtaBOMTFFine->SetTitle("");
  ghostDiffReEmuEtaBOMTFFine->SetLineWidth(2);
  ghostDiffReEmuEtaBOMTFFine->SetLineColor(kOrange);
  ghostDiffReEmuEtaBOMTFFine->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuEtaBOMTFFine->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuEtaBOMTFFine->SetMarkerStyle(23);
  ghostDiffReEmuEtaBOMTFFine->SetMarkerColor(kOrange);

  ghostDiffUnpackEtaBOMTFFine->SetLineWidth(2);
  ghostDiffUnpackEtaBOMTFFine->SetLineColor(kBlue);
  ghostDiffUnpackEtaBOMTFFine->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackEtaBOMTFFine->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackEtaBOMTFFine->SetMarkerStyle(20);
  ghostDiffUnpackEtaBOMTFFine->SetMarkerColor(kBlue);

  ghostDiffUnpackEtaBOMTFFine->Draw("E1HIST");
  ghostDiffReEmuEtaBOMTFFine->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg2 = new TLegend(0.15,0.73,0.7,0.88);
  leg2->SetFillColor(0);
  leg2->AddEntry(ghostDiffReEmuEtaBOMTFFine, "ReEmu","lp");
  leg2->AddEntry(ghostDiffUnpackEtaBOMTFFine, "Unpack","lp");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->Draw();
  leg2->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF fine/OMTF");

  c2->SaveAs(plotFolder + "dEtaBOMTFFine_" + run + ".pdf");
  c2->SaveAs(plotFolder + "dEtaBOMTFFine_" + run + ".png");

  TCanvas* c21= new TCanvas;

  ghostDiffReEmuEtaBOMTFCoarse->SetTitle("");
  ghostDiffReEmuEtaBOMTFCoarse->SetLineWidth(2);
  ghostDiffReEmuEtaBOMTFCoarse->SetLineColor(kOrange);
  ghostDiffReEmuEtaBOMTFCoarse->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuEtaBOMTFCoarse->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuEtaBOMTFCoarse->SetMarkerStyle(23);
  ghostDiffReEmuEtaBOMTFCoarse->SetMarkerColor(kOrange);

  ghostDiffUnpackEtaBOMTFCoarse->SetLineWidth(2);
  ghostDiffUnpackEtaBOMTFCoarse->SetLineColor(kBlue);
  ghostDiffUnpackEtaBOMTFCoarse->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackEtaBOMTFCoarse->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackEtaBOMTFCoarse->SetMarkerStyle(20);
  ghostDiffUnpackEtaBOMTFCoarse->SetMarkerColor(kBlue);

  ghostDiffUnpackEtaBOMTFCoarse->Draw("E1HIST");
  ghostDiffReEmuEtaBOMTFCoarse->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg21 = new TLegend(0.15,0.73,0.7,0.88);
  leg21->SetFillColor(0);
  leg21->AddEntry(ghostDiffReEmuEtaBOMTFCoarse, "ReEmu","lp");
  leg21->AddEntry(ghostDiffUnpackEtaBOMTFCoarse, "Unpack","lp");
  leg21->SetBorderSize(0);
  leg21->SetFillStyle(0);
  leg21->Draw();
  leg21->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF coarse/OMTF");

  c21->SaveAs(plotFolder + "dEtaBOMTFCoarse_" + run + ".pdf");
  c21->SaveAs(plotFolder + "dEtaBOMTFCoarse_" + run + ".png");


  TCanvas* c3 = new TCanvas;

  ghostDiffReEmuEtaOMTF->SetTitle("");
  ghostDiffReEmuEtaOMTF->SetLineWidth(2);
  ghostDiffReEmuEtaOMTF->SetLineColor(kOrange);
  ghostDiffReEmuEtaOMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuEtaOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuEtaOMTF->SetMarkerStyle(23);
  ghostDiffReEmuEtaOMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuEtaOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackEtaOMTF->SetLineWidth(2);
  ghostDiffUnpackEtaOMTF->SetLineColor(kBlue);
  ghostDiffUnpackEtaOMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackEtaOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackEtaOMTF->SetMarkerStyle(20);
  ghostDiffUnpackEtaOMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuEtaOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackEtaOMTF->Draw("E1HIST");
  ghostDiffReEmuEtaOMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg3 = new TLegend(0.15,0.73,0.7,0.88);
  leg3->SetFillColor(0);
  leg3->AddEntry(ghostDiffReEmuEtaOMTF, "ReEmu","lp");
  leg3->AddEntry(ghostDiffUnpackEtaOMTF, "Unpack","lp");
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->Draw();
  leg3->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF");

  c3->SaveAs(plotFolder + "dEtaOMTF_" + run + ".pdf");
  c3->SaveAs(plotFolder + "dEtaOMTF_" + run + ".png");


  TCanvas* c4 = new TCanvas;

  ghostDiffReEmuEtaEMTF->SetTitle("");
  ghostDiffReEmuEtaEMTF->SetLineWidth(2);
  ghostDiffReEmuEtaEMTF->SetLineColor(kOrange);
  ghostDiffReEmuEtaEMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuEtaEMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuEtaEMTF->SetMarkerStyle(23);
  ghostDiffReEmuEtaEMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuEtaEMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackEtaEMTF->SetLineWidth(2);
  ghostDiffUnpackEtaEMTF->SetLineColor(kBlue);
  ghostDiffUnpackEtaEMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackEtaEMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackEtaEMTF->SetMarkerStyle(20);
  ghostDiffUnpackEtaEMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuEtaEMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackEtaEMTF->Draw("E1HIST");
  ghostDiffReEmuEtaEMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg4 = new TLegend(0.15,0.73,0.7,0.88);
  leg4->SetFillColor(0);
  leg4->AddEntry(ghostDiffReEmuEtaEMTF, "ReEmu","lp");
  leg4->AddEntry(ghostDiffUnpackEtaEMTF, "Unpack","lp");
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->Draw();
  leg4->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, EMTF");

  c4->SaveAs(plotFolder + "dEtaEMTF_" + run + ".pdf");
  c4->SaveAs(plotFolder + "dEtaEMTF_" + run + ".png");


  // dPhi

  TCanvas* c5 = new TCanvas;

  ghostDiffReEmuPhiEOMTF->SetTitle("");
  ghostDiffReEmuPhiEOMTF->SetLineWidth(2);
  ghostDiffReEmuPhiEOMTF->SetLineColor(kOrange);
  ghostDiffReEmuPhiEOMTF->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuPhiEOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuPhiEOMTF->SetMarkerStyle(23);
  ghostDiffReEmuPhiEOMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuPhiEOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackPhiEOMTF->SetLineWidth(2);
  ghostDiffUnpackPhiEOMTF->SetLineColor(kBlue);
  ghostDiffUnpackPhiEOMTF->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackPhiEOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackPhiEOMTF->SetMarkerStyle(20);
  ghostDiffUnpackPhiEOMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuPhiEOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackPhiEOMTF->Draw("E1HIST");
  ghostDiffReEmuPhiEOMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg5 = new TLegend(0.15,0.73,0.7,0.88);
  leg5->SetFillColor(0);
  leg5->AddEntry(ghostDiffReEmuPhiEOMTF, "ReEmu","lp");
  leg5->AddEntry(ghostDiffUnpackPhiEOMTF, "Unpack","lp");
  leg5->SetBorderSize(0);
  leg5->SetFillStyle(0);
  leg5->Draw();
  leg5->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF/EMTF");

  c5->SaveAs(plotFolder + "dPhiEOMTF_" + run + ".pdf");
  c5->SaveAs(plotFolder + "dPhiEOMTF_" + run + ".png");

  TCanvas* c6 = new TCanvas;

  ghostDiffReEmuPhiBOMTFFine->SetTitle("");
  ghostDiffReEmuPhiBOMTFFine->SetLineWidth(2);
  ghostDiffReEmuPhiBOMTFFine->SetLineColor(kOrange);
  ghostDiffReEmuPhiBOMTFFine->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuPhiBOMTFFine->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuPhiBOMTFFine->SetMarkerStyle(23);
  ghostDiffReEmuPhiBOMTFFine->SetMarkerColor(kOrange);

  ghostDiffUnpackPhiBOMTFFine->SetLineWidth(2);
  ghostDiffUnpackPhiBOMTFFine->SetLineColor(kBlue);
  ghostDiffUnpackPhiBOMTFFine->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackPhiBOMTFFine->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackPhiBOMTFFine->SetMarkerStyle(20);
  ghostDiffUnpackPhiBOMTFFine->SetMarkerColor(kBlue);

  ghostDiffUnpackPhiBOMTFFine->Draw("E1HIST");
  ghostDiffReEmuPhiBOMTFFine->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg61 = new TLegend(0.15,0.73,0.7,0.88);
  leg61->SetFillColor(0);
  leg61->AddEntry(ghostDiffReEmuPhiBOMTFFine, "ReEmu","lp");
  leg61->AddEntry(ghostDiffUnpackPhiBOMTFFine, "Unpack","lp");
  leg61->SetBorderSize(0);
  leg61->SetFillStyle(0);
  leg61->Draw();
  leg61->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF fine/OMTF");

  c6->SaveAs(plotFolder + "dPhiBOMTFFine_" + run + ".pdf");
  c6->SaveAs(plotFolder + "dPhiBOMTFFine_" + run + ".png");

  TCanvas* c61 = new TCanvas;

  ghostDiffReEmuPhiBOMTFCoarse->SetTitle("");
  ghostDiffReEmuPhiBOMTFCoarse->SetLineWidth(2);
  ghostDiffReEmuPhiBOMTFCoarse->SetLineColor(kOrange);
  ghostDiffReEmuPhiBOMTFCoarse->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuPhiBOMTFCoarse->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuPhiBOMTFCoarse->SetMarkerStyle(23);
  ghostDiffReEmuPhiBOMTFCoarse->SetMarkerColor(kOrange);

  ghostDiffUnpackPhiBOMTFCoarse->SetLineWidth(2);
  ghostDiffUnpackPhiBOMTFCoarse->SetLineColor(kBlue);
  ghostDiffUnpackPhiBOMTFCoarse->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackPhiBOMTFCoarse->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackPhiBOMTFCoarse->SetMarkerStyle(20);
  ghostDiffUnpackPhiBOMTFCoarse->SetMarkerColor(kBlue);

  ghostDiffUnpackPhiBOMTFCoarse->Draw("E1HIST");
  ghostDiffReEmuPhiBOMTFCoarse->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg6 = new TLegend(0.15,0.73,0.7,0.88);
  leg6->SetFillColor(0);
  leg6->AddEntry(ghostDiffReEmuPhiBOMTFCoarse, "ReEmu","lp");
  leg6->AddEntry(ghostDiffUnpackPhiBOMTFCoarse, "Unpack","lp");
  leg6->SetBorderSize(0);
  leg6->SetFillStyle(0);
  leg6->Draw();
  leg6->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF coarse/OMTF");

  c61->SaveAs(plotFolder + "dPhiBOMTFCoarse_" + run + ".pdf");
  c61->SaveAs(plotFolder + "dPhiBOMTFCoarse_" + run + ".png");

  TCanvas* c7 = new TCanvas;

  ghostDiffReEmuPhiOMTF->SetTitle("");
  ghostDiffReEmuPhiOMTF->SetLineWidth(2);
  ghostDiffReEmuPhiOMTF->SetLineColor(kOrange);
  ghostDiffReEmuPhiOMTF->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuPhiOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuPhiOMTF->SetMarkerStyle(23);
  ghostDiffReEmuPhiOMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuPhiOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackPhiOMTF->SetLineWidth(2);
  ghostDiffUnpackPhiOMTF->SetLineColor(kBlue);
  ghostDiffUnpackPhiOMTF->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackPhiOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackPhiOMTF->SetMarkerStyle(20);
  ghostDiffUnpackPhiOMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuPhiOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackPhiOMTF->Draw("E1HIST");
  ghostDiffReEmuPhiOMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg7 = new TLegend(0.15,0.73,0.7,0.88);
  leg7->SetFillColor(0);
  leg7->AddEntry(ghostDiffReEmuPhiOMTF, "ReEmu","lp");
  leg7->AddEntry(ghostDiffUnpackPhiOMTF, "Unpack","lp");
  leg7->SetBorderSize(0);
  leg7->SetFillStyle(0);
  leg7->Draw();
  leg7->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF");

  c7->SaveAs(plotFolder + "dPhiOMTF_" + run + ".pdf");
  c7->SaveAs(plotFolder + "dPhiOMTF_" + run + ".png");


  TCanvas* c8 = new TCanvas;

  ghostDiffReEmuPhiEMTF->SetTitle("");
  ghostDiffReEmuPhiEMTF->SetLineWidth(2);
  ghostDiffReEmuPhiEMTF->SetLineColor(kOrange);
  ghostDiffReEmuPhiEMTF->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuPhiEMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuPhiEMTF->SetMarkerStyle(23);
  ghostDiffReEmuPhiEMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuPhiEMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackPhiEMTF->SetLineWidth(2);
  ghostDiffUnpackPhiEMTF->SetLineColor(kBlue);
  ghostDiffUnpackPhiEMTF->GetXaxis()->SetTitle("#Delta#phi(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackPhiEMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackPhiEMTF->SetMarkerStyle(20);
  ghostDiffUnpackPhiEMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuPhiEMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackPhiEMTF->Draw("E1HIST");
  ghostDiffReEmuPhiEMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg8 = new TLegend(0.15,0.73,0.7,0.88);
  leg8->SetFillColor(0);
  leg8->AddEntry(ghostDiffReEmuPhiEMTF, "ReEmu","lp");
  leg8->AddEntry(ghostDiffUnpackPhiEMTF, "Unpack","lp");
  leg8->SetBorderSize(0);
  leg8->SetFillStyle(0);
  leg8->Draw();
  leg8->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, EMTF");

  c8->SaveAs(plotFolder + "dPhiEMTF_" + run + ".pdf");
  c8->SaveAs(plotFolder + "dPhiEMTF_" + run + ".png");


  // dR

  TCanvas* c9 = new TCanvas;

  ghostDiffReEmuREOMTF->SetTitle("");
  ghostDiffReEmuREOMTF->SetLineWidth(2);
  ghostDiffReEmuREOMTF->SetLineColor(kOrange);
  ghostDiffReEmuREOMTF->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuREOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuREOMTF->SetMarkerStyle(23);
  ghostDiffReEmuREOMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuREOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackREOMTF->SetLineWidth(2);
  ghostDiffUnpackREOMTF->SetLineColor(kBlue);
  ghostDiffUnpackREOMTF->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackREOMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackREOMTF->SetMarkerStyle(20);
  ghostDiffUnpackREOMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuREOMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackREOMTF->Draw("E1HIST");
  ghostDiffReEmuREOMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg9 = new TLegend(0.15,0.73,0.7,0.88);
  leg9->SetFillColor(0);
  leg9->AddEntry(ghostDiffReEmuREOMTF, "ReEmu","lp");
  leg9->AddEntry(ghostDiffUnpackREOMTF, "Unpack","lp");
  leg9->SetBorderSize(0);
  leg9->SetFillStyle(0);
  leg9->Draw();
  leg9->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF/EMTF");

  c9->SaveAs(plotFolder + "dREOMTF_" + run + ".pdf");
  c9->SaveAs(plotFolder + "dREOMTF_" + run + ".png");

  TCanvas* c10 = new TCanvas;

  ghostDiffReEmuRBOMTFFine->SetTitle("");
  ghostDiffReEmuRBOMTFFine->SetLineWidth(2);
  ghostDiffReEmuRBOMTFFine->SetLineColor(kOrange);
  ghostDiffReEmuRBOMTFFine->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuRBOMTFFine->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuRBOMTFFine->SetMarkerStyle(23);
  ghostDiffReEmuRBOMTFFine->SetMarkerColor(kOrange);

  ghostDiffUnpackRBOMTFFine->SetLineWidth(2);
  ghostDiffUnpackRBOMTFFine->SetLineColor(kBlue);
  ghostDiffUnpackRBOMTFFine->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackRBOMTFFine->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackRBOMTFFine->SetMarkerStyle(20);
  ghostDiffUnpackRBOMTFFine->SetMarkerColor(kBlue);

  ghostDiffUnpackRBOMTFFine->Draw("E1HIST");
  ghostDiffReEmuRBOMTFFine->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg10 = new TLegend(0.15,0.73,0.7,0.88);
  leg10->SetFillColor(0);
  leg10->AddEntry(ghostDiffReEmuRBOMTFFine, "ReEmu","lp");
  leg10->AddEntry(ghostDiffUnpackRBOMTFFine, "Unpack","lp");
  leg10->SetBorderSize(0);
  leg10->SetFillStyle(0);
  leg10->Draw();
  leg10->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF fine/OMTF");

  c10->SaveAs(plotFolder + "dRBOMTFFine_" + run + ".pdf");
  c10->SaveAs(plotFolder + "dRBOMTFFine_" + run + ".png");

  TCanvas* c101 = new TCanvas;

  ghostDiffReEmuRBOMTFCoarse->SetTitle("");
  ghostDiffReEmuRBOMTFCoarse->SetLineWidth(2);
  ghostDiffReEmuRBOMTFCoarse->SetLineColor(kOrange);
  ghostDiffReEmuRBOMTFCoarse->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuRBOMTFCoarse->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuRBOMTFCoarse->SetMarkerStyle(23);
  ghostDiffReEmuRBOMTFCoarse->SetMarkerColor(kOrange);

  ghostDiffUnpackRBOMTFCoarse->SetLineWidth(2);
  ghostDiffUnpackRBOMTFCoarse->SetLineColor(kBlue);
  ghostDiffUnpackRBOMTFCoarse->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackRBOMTFCoarse->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackRBOMTFCoarse->SetMarkerStyle(20);
  ghostDiffUnpackRBOMTFCoarse->SetMarkerColor(kBlue);

  ghostDiffUnpackRBOMTFCoarse->Draw("E1HIST");
  ghostDiffReEmuRBOMTFCoarse->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg101 = new TLegend(0.15,0.73,0.7,0.88);
  leg101->SetFillColor(0);
  leg101->AddEntry(ghostDiffReEmuRBOMTFCoarse, "ReEmu","lp");
  leg101->AddEntry(ghostDiffUnpackRBOMTFCoarse, "Unpack","lp");
  leg101->SetBorderSize(0);
  leg101->SetFillStyle(0);
  leg101->Draw();
  leg101->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF coarse/OMTF");

  c101->SaveAs(plotFolder + "dRBOMTFCoarse_" + run + ".pdf");
  c101->SaveAs(plotFolder + "dRBOMTFCoarse_" + run + ".png");

  TCanvas* c11 = new TCanvas;

  ghostDiffReEmuROMTF->SetTitle("");
  ghostDiffReEmuROMTF->SetLineWidth(2);
  ghostDiffReEmuROMTF->SetLineColor(kOrange);
  ghostDiffReEmuROMTF->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuROMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuROMTF->SetMarkerStyle(23);
  ghostDiffReEmuROMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuROMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackROMTF->SetLineWidth(2);
  ghostDiffUnpackROMTF->SetLineColor(kBlue);
  ghostDiffUnpackROMTF->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackROMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackROMTF->SetMarkerStyle(20);
  ghostDiffUnpackROMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuROMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackROMTF->Draw("E1HIST");
  ghostDiffReEmuROMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg11 = new TLegend(0.15,0.73,0.7,0.88);
  leg11->SetFillColor(0);
  leg11->AddEntry(ghostDiffReEmuROMTF, "ReEmu","lp");
  leg11->AddEntry(ghostDiffUnpackROMTF, "Unpack","lp");
  leg11->SetBorderSize(0);
  leg11->SetFillStyle(0);
  leg11->Draw();
  leg11->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF");

  c11->SaveAs(plotFolder + "dROMTF_" + run + ".pdf");
  c11->SaveAs(plotFolder + "dROMTF_" + run + ".png");


  TCanvas* c12 = new TCanvas;

  ghostDiffReEmuREMTF->SetTitle("");
  ghostDiffReEmuREMTF->SetLineWidth(2);
  ghostDiffReEmuREMTF->SetLineColor(kOrange);
  ghostDiffReEmuREMTF->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffReEmuREMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffReEmuREMTF->SetMarkerStyle(23);
  ghostDiffReEmuREMTF->SetMarkerColor(kOrange);
  // ghostDiffReEmuREMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackREMTF->SetLineWidth(2);
  ghostDiffUnpackREMTF->SetLineColor(kBlue);
  ghostDiffUnpackREMTF->GetXaxis()->SetTitle("#DeltaR(#mu_{L1}, #mu_{Ghost})");
  ghostDiffUnpackREMTF->GetYaxis()->SetTitle("Counts");
  ghostDiffUnpackREMTF->SetMarkerStyle(20);
  ghostDiffUnpackREMTF->SetMarkerColor(kBlue);
  // ghostDiffReEmuREMTF->GetYaxis()->SetRangeUser(1, 1e3);

  ghostDiffUnpackREMTF->Draw("E1HIST");
  ghostDiffReEmuREMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg12 = new TLegend(0.15,0.73,0.7,0.88);
  leg12->SetFillColor(0);
  leg12->AddEntry(ghostDiffReEmuREMTF, "ReEmu","lp");
  leg12->AddEntry(ghostDiffUnpackREMTF, "Unpack","lp");
  leg12->SetBorderSize(0);
  leg12->SetFillStyle(0);
  leg12->Draw();
  leg12->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, EMTF");

  c12->SaveAs(plotFolder + "dREMTF_" + run + ".pdf");
  c12->SaveAs(plotFolder + "dREMTF_" + run + ".png");

  TCanvas* c20 = new TCanvas;

  resolutionReEmuEtaOMTF->SetTitle("");
  resolutionReEmuEtaOMTF->SetLineWidth(2);
  resolutionReEmuEtaOMTF->SetLineColor(kOrange);
  resolutionReEmuEtaOMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{reco})");
  resolutionReEmuEtaOMTF->GetYaxis()->SetTitle("Counts");
  resolutionReEmuEtaOMTF->SetMarkerStyle(23);
  resolutionReEmuEtaOMTF->SetMarkerColor(kOrange);

  resolutionUnpackEtaOMTF->SetLineWidth(2);
  resolutionUnpackEtaOMTF->SetLineColor(kBlue);
  resolutionUnpackEtaOMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{reco})");
  resolutionUnpackEtaOMTF->GetYaxis()->SetTitle("Counts");
  resolutionUnpackEtaOMTF->SetMarkerStyle(20);
  resolutionUnpackEtaOMTF->SetMarkerColor(kBlue);

  resolutionUnpackEtaOMTF->Draw("E1HIST");
  resolutionReEmuEtaOMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg20 = new TLegend(0.15,0.73,0.7,0.88);
  leg20->SetFillColor(0);
  leg20->AddEntry(resolutionReEmuEtaOMTF, "ReEmu","lp");
  leg20->AddEntry(resolutionUnpackEtaOMTF, "Unpack","lp");
  leg20->SetBorderSize(0);
  leg20->SetFillStyle(0);
  leg20->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, OMTF");

  c20->SaveAs(plotFolder + "dEtaRecoOMTF_" + run + ".pdf");
  c20->SaveAs(plotFolder + "dEtaRecoOMTF_" + run + ".png");

  TCanvas* c23 = new TCanvas;

  resolutionReEmuEtaBMTF->SetTitle("");
  resolutionReEmuEtaBMTF->SetLineWidth(2);
  resolutionReEmuEtaBMTF->SetLineColor(kOrange);
  resolutionReEmuEtaBMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{reco})");
  resolutionReEmuEtaBMTF->GetYaxis()->SetTitle("Counts");
  resolutionReEmuEtaBMTF->SetMarkerStyle(23);
  resolutionReEmuEtaBMTF->SetMarkerColor(kOrange);

  resolutionUnpackEtaBMTF->SetLineWidth(2);
  resolutionUnpackEtaBMTF->SetLineColor(kBlue);
  resolutionUnpackEtaBMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{reco})");
  resolutionUnpackEtaBMTF->GetYaxis()->SetTitle("Counts");
  resolutionUnpackEtaBMTF->SetMarkerStyle(20);
  resolutionUnpackEtaBMTF->SetMarkerColor(kBlue);

  resolutionUnpackEtaBMTF->Draw("E1HIST");
  resolutionReEmuEtaBMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg23 = new TLegend(0.15,0.73,0.7,0.88);
  leg23->SetFillColor(0);
  leg23->AddEntry(resolutionReEmuEtaBMTF, "ReEmu","lp");
  leg23->AddEntry(resolutionUnpackEtaBMTF, "Unpack","lp");
  leg23->SetBorderSize(0);
  leg23->SetFillStyle(0);
  leg23->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, BMTF");

  c23->SaveAs(plotFolder + "dEtaRecoBMTF_" + run + ".pdf");
  c23->SaveAs(plotFolder + "dEtaRecoBMTF_" + run + ".png");

  TCanvas* c22 = new TCanvas;

  resolutionReEmuEtaEMTF->SetTitle("");
  resolutionReEmuEtaEMTF->SetLineWidth(2);
  resolutionReEmuEtaEMTF->SetLineColor(kOrange);
  resolutionReEmuEtaEMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{reco})");
  resolutionReEmuEtaEMTF->GetYaxis()->SetTitle("Counts");
  resolutionReEmuEtaEMTF->SetMarkerStyle(23);
  resolutionReEmuEtaEMTF->SetMarkerColor(kOrange);

  resolutionUnpackEtaEMTF->SetLineWidth(2);
  resolutionUnpackEtaEMTF->SetLineColor(kBlue);
  resolutionUnpackEtaEMTF->GetXaxis()->SetTitle("#Delta#eta(#mu_{L1}, #mu_{reco})");
  resolutionUnpackEtaEMTF->GetYaxis()->SetTitle("Counts");
  resolutionUnpackEtaEMTF->SetMarkerStyle(20);
  resolutionUnpackEtaEMTF->SetMarkerColor(kBlue);

  resolutionUnpackEtaEMTF->Draw("E1HIST");
  resolutionReEmuEtaEMTF->Draw("same,E1HIST");
  gPad->Modified();

  TLegend* leg22 = new TLegend(0.15,0.73,0.7,0.88);
  leg22->SetFillColor(0);
  leg22->AddEntry(resolutionReEmuEtaEMTF, "ReEmu","lp");
  leg22->AddEntry(resolutionUnpackEtaEMTF, "Unpack","lp");
  leg22->SetBorderSize(0);
  leg22->SetFillStyle(0);
  leg22->Draw();
  n3.DrawLatex(0.15, 0.6, "Run " + run + " #sqrt{s} = 13 TeV");
  n4.DrawLatex(0.15, 0.55, "SingleMu, EMTF");

  c22->SaveAs(plotFolder + "dEtaRecoEMTF_" + run + ".pdf");
  c22->SaveAs(plotFolder + "dEtaRecoEMTF_" + run + ".png");

}
