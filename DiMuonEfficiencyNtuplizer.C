#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"

#include <iostream>
#include <string>
#include <vector>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeTfMuonDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

enum tftype { bmtf, omtf, emtf, none };

std::vector<std::string> getNtupleList(std::string fname);
double calcTFphi(int locPhi, tftype tfType, int proc);
double calcTFeta(int eta);
double calcTFpt(int pt);
void findUgmtMuons(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int& mu1,
                   int& mu2);
void findTFMuons(L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf_, int& mu1,
                 int& mu2, int& pt1, int& pt2,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2,
                 tftype tfType);
std::vector<std::string> generateGenPhysicsQuantities();
std::vector<std::string> generateUgmtPhysicsQuantities();
std::vector<std::string> generateTfPhysicsQuantities();
std::vector<std::string> createContentList(
    int nGenMu, std::vector<std::string> genPhysicsQuantities,
    std::vector<std::string> l1PhysicsQuantities);
void fillNtuple(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                std::vector<std::string> contentList, float ntupleValues[]);
void fillNtuple(L1Analysis::L1AnalysisGeneratorDataFormat* gen_, int genMu1,
                int genMu2, std::vector<std::string> contentList,
                float ntupleValues[]);
void fillNtuple(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int ugmtMu1,
                int ugmtMu2, std::vector<std::string> contentList,
                float ugmtNtupleValues[]);
void fillNtuple(int tfMu1,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1_,
                int tfMu2,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2_,
                std::vector<std::string> contentList, float tfNtupleValues[]);

void DiMuonEfficiencyNtuplizer(std::string fname = "L1Ntuple_list",
                               int nGenMu = 1, std::string outDir = "tmp/") {
  // make trees and set branch addresses
  const char* ugmtMcTree = "l1UpgradeEmuTree/L1UpgradeTree";
  const char* tfMcTree = "l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree";
  const char* genTree = "l1GeneratorTree/L1GenTree";
  const char* recoTree = "l1MuonRecoTree/Muon2RecoTree";

  std::vector<std::string> listNtuples = getNtupleList(fname);

  // CheckFirstFile
  if (listNtuples.size() == 0) return false;

  TFile* file = TFile::Open(listNtuples[0].c_str());

  if (file == 0) return false;
  if (file->IsOpen() == 0) return false;

  bool foundUgmtMcTree = false;
  bool foundTfMcTree = false;
  bool foundGen = false;
  bool foundReco = false;
  TTree* treeUgmtMC = (TTree*)file->Get(ugmtMcTree);
  TTree* treeTfMC = (TTree*)file->Get(tfMcTree);
  TTree* treeGen = (TTree*)file->Get(genTree);
  TTree* treeReco = (TTree*)file->Get(recoTree);

  if (!treeTfMC && !treeUgmtMC && !(treeReco || treeGen)) {
    std::cout << "Require both TF upgrade and upgrade trees as well as one of "
                 "reco or gen tree to run. Exiting.. "
              << std::endl;
    return;
  }

  // OpenWithoutInit
  TChain* chainUgmtMC = new TChain(ugmtMcTree);
  TChain* chainTfMC = new TChain(tfMcTree);
  TChain* chainGen = new TChain(genTree);
  TChain* chainReco = new TChain(recoTree);
  for (std::vector<string>::iterator it = listNtuples.begin();
       it != listNtuples.end(); ++it) {
    std::cout << " -- Adding " << *it << std::endl;
    chainUgmtMC->Add(it->c_str());
    chainTfMC->Add(it->c_str());
    chainGen->Add(it->c_str());
    chainReco->Add(it->c_str());
  }

  // Init
  std::cout << "Estimate the number of entries..." << std::endl;
  int nevents = chainReco->GetEntries();
  std::cout << nevents << std::endl;

  L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmtMC_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* bmtfMC_ =
      new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* omtfMC_ =
      new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* emtfMC_ =
      new L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_ =
      new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  L1Analysis::L1AnalysisGeneratorDataFormat* gen_ =
      new L1Analysis::L1AnalysisGeneratorDataFormat();

  if (treeUgmtMC) {
    std::cout << "Found upgrade tree from MC, adding branches.." << std::endl;
    foundUgmtTree = true;
    chainUgmtMC->SetBranchAddress("L1Upgrade", &ugmtMC_);
  }
  if (treeTfMC) {
    std::cout << "Found TF upgrade tree from MC, adding branches.."
              << std::endl;
    foundTfMcTree = true;
    chainTfMC->SetBranchAddress("L1UpgradeBmtfMuon", &bmtfMC_);
    chainTfMC->SetBranchAddress("L1UpgradeOmtfMuon", &omtfMC_);
    chainTfMC->SetBranchAddress("L1UpgradeEmtfMuon", &emtfMC_);
  }
  if (treeGen) {
    std::cout << "Found gen tree, adding branches.." << std::endl;
    foundGen = true;
    chainGen->SetBranchAddress("Generator", &gen_);
  }
  if (treeReco) {
    std::cout << "Found reco tree, adding branches.." << std::endl;
    foundReco = true;
    chainReco->SetBranchAddress("Muon", &reco_);
  }

  std::string ntupleName;
  if (nGenMu == 1) {
    ntupleName = "SingleMuNtuple";
  } else if (nGenMu == 2) {
    ntupleName = "DimuonNtuple";
  }

  std::string gmtNtupleFname(outDir + "/uGMT" + ntupleName + ".root");
  std::string tfNtupleFname(outDir + "/tf" + ntupleName + ".root");

  TFile* ugmtFile = new TFile(gmtNtupleFname.c_str(), "RECREATE");
  TFile* tfFile = new TFile(tfNtupleFname.c_str(), "RECREATE");

  std::vector<std::string> ugmtContentList = createContentList(
      nGenMu, generateGenPhysicsQuantities(), generateUgmtPhysicsQuantities());
  std::vector<std::string> tfContentList = createContentList(
      nGenMu, generateGenPhysicsQuantities(), generateTfPhysicsQuantities());

  std::ostringstream ugmtContentStream;
  for (std::vector<std::string>::iterator it = ugmtContentList.begin();
       it != ugmtContentList.end(); ++it) {
    ugmtContentStream << ":" << *it;
  }
  std::string ugmtContentStr(ugmtContentStream.str());

  std::ostringstream tfContentStream;
  for (std::vector<std::string>::iterator it = tfContentList.begin();
       it != tfContentList.end(); ++it) {
    tfContentStream << ":" << *it;
  }
  std::string tfContentStr(tfContentStream.str());

  TNtuple* ugmtNtuple =
      new TNtuple("ugmt_ntuple", "ntupledump", ugmtContentStr.c_str());
  TNtuple* tfNtuple =
      new TNtuple("tf_ntuple", "ntupledump", tfContentStr.c_str());

  Float_t ugmtNtupleValues[ugmtContentList.size()];
  Float_t tfNtupleValues[tfContentList.size()];

  // Fill ntuples
  std::cout << "Running over " << nevents << std::endl;

  for (Long64_t jentry = 0; jentry < nevents; ++jentry) {
    if ((jentry % 1000) == 0) {
      std::cout << "Done " << jentry << " events..." << std::endl;
    }

    chainUgmtMC->GetEntry(jentry);
    chainTfMC->GetEntry(jentry);

    if (foundReco) {
      chainReco->GetEntry(jentry);
      // Check for correct number of reco muons
      if (reco_->nMuons != nGenMu) {
        continue;
      }

      fillNtuple(reco_, ugmtContentList, ugmtNtupleValues);
      fillNtuple(reco_, tfContentList, tfNtupleValues);
    } else {
      chainGen->GetEntry(jentry);
      int genMu1 = -1;
      int genMu2 = -1;
      // Check for correct number of gen muons
      int nGenMus = 0;
      for (int i = 0; i < gen_->nPart; ++i) {
        if (abs(gen_->partId[i]) == 13) {
          ++nGenMus;
          genMu2 = genMu1;
          genMu1 = i;
        }
        if (nGenMus != nGenMu) {
          continue;
        }
      }
      fillNtuple(gen_, genMu1, genMu2, ugmtContentList, ugmtNtupleValues);
      fillNtuple(gen_, genMu1, genMu2, tfContentList, tfNtupleValues);
    }

    // Find two highest pT uGMT muons
    int ugmtMu1 = -1;
    int ugmtMu2 = -1;
    findUgmtMuons(ugmtMC_, ugmtMu1, ugmtMu2);

    // Find two highest pT TF muons
    int pt1 = 0;
    int pt2 = 0;
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1;
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2;
    int tfMu1 = -1;
    int tfMu2 = -1;
    findTFMuons(bmtfMC_, tfMu1, tfMu2, pt1, pt2, tf1, tf2, tftype::bmtf);
    findTFMuons(omtfMC_, tfMu1, tfMu2, pt1, pt2, tf1, tf2, tftype::omtf);
    findTFMuons(emtfMC_, tfMu1, tfMu2, pt1, pt2, tf1, tf2, tftype::emtf);

    // Fill ugmt ntuple
    ugmtFile->cd();
    fillNtuple(ugmtMC_, ugmtMu1, ugmtMu2, ugmtContentList, ugmtNtupleValues);
    ugmtNtuple->Fill(ugmtNtupleValues);
    // Fill tf ntuple
    tfFile->cd();
    fillNtuple(tfMu1, tf1, tfMu2, tf2, tfContentList, tfNtupleValues);
    tfNtuple->Fill(tfNtupleValues);
  }

  // Write ntuple file to disk.
  ugmtFile->cd();
  ugmtFile->Write();
  tfFile->cd();
  tfFile->Write();
}

std::vector<std::string> getNtupleList(std::string fname) {
  std::vector<std::string> listNtuples;
  std::ifstream flist(fname);
  if (!flist) {
    std::cout << "File " << fname << " is not found!" << std::endl;
    // TODO: Throw exception..
  }

  while (!flist.eof()) {
    std::string str;
    getline(flist, str);
    if (!flist.fail()) {
      if (str != "") listNtuples.push_back(str);
    }
  }

  return listNtuples;
}

double calcTFphi(int locPhi, tftype tfType, int proc) {
  int globPhi = 0;
  if (tfType == tftype::bmtf) {
    // each BMTF processor corresponds to a 30 degree wedge = 48 in int-scale
    globPhi = (proc)*48 + locPhi;
    // first processor starts at CMS phi = -15 degrees...
    globPhi += 576 - 24;
    // handle wrap-around (since we add the 576-24, the value will never be
    // negative!)
    globPhi = globPhi % 576;
  } else {
    // all others correspond to 60 degree sectors = 96 in int-scale
    globPhi = (proc)*96 + locPhi;
    // first processor starts at CMS phi = 15 degrees (24 in int)... Handle
    // wrap-around with %. Add 576 to make sure the number is positive
    globPhi = (globPhi + 600) % 576;
  }
  return 0.010908 * globPhi;
}

double calcTFeta(int eta) { return 0.010875 * eta; }

double calcTFpt(int pt) { return 0.5 * (pt - 1); }

void findUgmtMuons(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int& mu1,
                   int& mu2) {
  float pt1 = 0;
  float pt2 = 0;
  for (int i = 0; i < ugmt_->nMuons; ++i) {
    if (ugmt_->muonBx[i] != 0) {
      continue;
    }
    if (ugmt_->muonEt[i] > pt1) {
      pt2 = pt1;
      mu2 = mu1;
      pt1 = ugmt_->muonEt[i];
      mu1 = i;
    } else if (ugmt_->muonEt[i] > pt2) {
      pt2 = ugmt_->muonEt[i];
      mu2 = i;
    }
  }
}

void findTFMuons(L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf_, int& mu1,
                 int& mu2, int& pt1, int& pt2,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2,
                 tftype tfType) {
  for (int i = 0; i < tf_->nTfMuons; ++i) {
    if (tf_->tfMuonBx[i] != 0) {
      continue;
    }
    float pT = (tf_->tfMuonHwPt[i] - 1) * 0.5;
    if (pT > pt1) {
      pt2 = pt1;
      mu2 = mu1;
      tf2 = tf1;
      pt1 = pT;
      mu1 = i;
      tf1 = tf_;
    } else if (pT > pt2) {
      pt2 = pT;
      mu2 = i;
      tf2 = tf_;
    }
  }
}

std::vector<std::string> generateGenPhysicsQuantities() {
  std::vector<std::string> physicsQuantities;
  physicsQuantities.push_back("pT");
  physicsQuantities.push_back("eta");
  physicsQuantities.push_back("phi");
  physicsQuantities.push_back("ch");

  return physicsQuantities;
}

std::vector<std::string> generateUgmtPhysicsQuantities() {
  std::vector<std::string> physicsQuantities;
  physicsQuantities.push_back("pT");
  physicsQuantities.push_back("eta");
  physicsQuantities.push_back("hf");
  physicsQuantities.push_back("phi");
  physicsQuantities.push_back("qual");
  physicsQuantities.push_back("ch");
  physicsQuantities.push_back("trkAddr");
  physicsQuantities.push_back("tfType");
  physicsQuantities.push_back("tfProcessor");

  return physicsQuantities;
}

std::vector<std::string> generateTfPhysicsQuantities() {
  std::vector<std::string> physicsQuantities;
  physicsQuantities.push_back("pT");
  physicsQuantities.push_back("eta");
  physicsQuantities.push_back("hf");
  physicsQuantities.push_back("phi");
  physicsQuantities.push_back("qual");
  physicsQuantities.push_back("ch");
  physicsQuantities.push_back("tfType");
  physicsQuantities.push_back("tfProcessor");

  return physicsQuantities;
}

std::vector<std::string> createContentList(
    int nGenMu, std::vector<std::string> genPhysicsQuantities,
    std::vector<std::string> l1PhysicsQuantities) {
  std::vector<std::string> contentList;
  for (std::vector<std::string>::iterator it = l1PhysicsQuantities.begin();
       it != l1PhysicsQuantities.end(); ++it) {
    contentList.push_back(*it + "1");
    contentList.push_back(*it + "2");
  }

  for (std::vector<std::string>::iterator it = genPhysicsQuantities.begin();
       it != genPhysicsQuantities.end(); ++it) {
    contentList.push_back(*it + "1_gen");
    if (nGenMu > 1) {
      contentList.push_back(*it + "2_gen");
    }
  }

  return contentList;
}

void fillNtuple(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                std::vector<std::string> contentList, float ntupleValues[]) {
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1_gen") {
      ntupleValues[i] = reco_->pt[0];
    } else if (contentList.at(i) == "pT1_gen") {
      ntupleValues[i] = reco_->pt[1];
    } else if (contentList.at(i) == "eta1_gen") {
      ntupleValues[i] = reco_->eta[0];
    } else if (contentList.at(i) == "eta1_gen") {
      ntupleValues[i] = reco_->eta[1];
    } else if (contentList.at(i) == "phi1_gen") {
      ntupleValues[i] = reco_->phi[0];
    } else if (contentList.at(i) == "phi1_gen") {
      ntupleValues[i] = reco_->phi[1];
    } else if (contentList.at(i) == "ch1_gen") {
      ntupleValues[i] = reco_->charge[0];
    } else if (contentList.at(i) == "ch1_gen") {
      ntupleValues[i] = reco_->charge[1];
    } else {
      ntupleValues[i] = -10;
    }
  }
}

void fillNtuple(L1Analysis::L1AnalysisGeneratorDataFormat* gen_, int genMu1,
                int genMu2, std::vector<std::string> contentList,
                float ntupleValues[]) {
  // TODO: Missing charge!
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1_gen") {
      ntupleValues[i] = gen_->partPt[genMu1];
    } else if (contentList.at(i) == "pT1_gen") {
      ntupleValues[i] = gen_->partPt[genMu2];
    } else if (contentList.at(i) == "eta1_gen") {
      ntupleValues[i] = gen_->partEta[genMu1];
    } else if (contentList.at(i) == "eta1_gen") {
      ntupleValues[i] = gen_->partEta[genMu2];
    } else if (contentList.at(i) == "phi1_gen") {
      ntupleValues[i] = gen_->partPhi[genMu1];
    } else if (contentList.at(i) == "phi1_gen") {
      ntupleValues[i] = gen_->partPhi[genMu2];
    } else {
      ntupleValues[i] = -10;
    }
  }
}

void fillNtuple(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int ugmtMu1,
                int ugmtMu2, std::vector<std::string> contentList,
                float ugmtNtupleValues[]) {
  // TODO: Missing track addresses, HF bit, processor
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1") {
      ugmtNtupleValues[i] = ugmt_->muonEt[ugmtMu1];
    } else if (contentList.at(i) == "pT2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonEt[ugmtMu2];
    } else if (contentList.at(i) == "eta1") {
      ugmtNtupleValues[i] = ugmt_->muonEta[ugmtMu1];
    } else if (contentList.at(i) == "eta2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonEta[ugmtMu2];
    } else if (contentList.at(i) == "phi1") {
      ugmtNtupleValues[i] = ugmt_->muonPhi[ugmtMu1];
    } else if (contentList.at(i) == "phi2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonPhi[ugmtMu2];
    } else if (contentList.at(i) == "qual1") {
      ugmtNtupleValues[i] = ugmt_->muonQual[ugmtMu1];
    } else if (contentList.at(i) == "qual2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonQual[ugmtMu2];
    } else if (contentList.at(i) == "ch1") {
      ugmtNtupleValues[i] = ugmt_->muonChg[ugmtMu1];
    } else if (contentList.at(i) == "ch2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonChg[ugmtMu2];
    } else if (contentList.at(i) == "tfType1") {
      int muIdx = ugmt_->muonTfMuonIdx[ugmtMu1];
      if (muIdx < 18 || muIdx >= 90) {
        ugmtNtupleValues[i] = 3;  // EMTF
      } else if (muIdx < 36 || muIdx >= 72) {
        ugmtNtupleValues[i] = 2;  // OMTF
      } else if (muIdx >= 36 && muIdx < 72) {
        ugmtNtupleValues[i] = 1;  // BMTF
      }
    } else if (contentList.at(i) == "tfType2" && ugmtMu2 != -1) {
      int muIdx = ugmt_->muonTfMuonIdx[ugmtMu2];
      if (muIdx < 18 || muIdx >= 90) {
        ugmtNtupleValues[i] = 3;  // EMTF
      } else if (muIdx < 36 || muIdx >= 72) {
        ugmtNtupleValues[i] = 2;  // OMTF
      } else if (muIdx >= 36 && muIdx < 72) {
        ugmtNtupleValues[i] = 1;  // BMTF
      }
    } else {
      ugmtNtupleValues[i] = -10;
    }
  }
}

void fillNtuple(int tfMu1,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1_,
                int tfMu2,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2_,
                std::vector<std::string> contentList, float tfNtupleValues[]) {
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1") {
      tfNtupleValues[i] = calcTFpt(tf1_->tfMuonHwPt[tfMu1]);
    } else if (contentList.at(i) == "pT2" && tfMu2 != -1) {
      tfNtupleValues[i] = calcTFpt(tf2_->tfMuonHwPt[tfMu2]);
    } else if (contentList.at(i) == "eta1") {
      tfNtupleValues[i] = calcTFeta(tf1_->tfMuonHwEta[tfMu1]);
    } else if (contentList.at(i) == "eta2" && tfMu2 != -1) {
      tfNtupleValues[i] = calcTFeta(tf2_->tfMuonHwEta[tfMu2]);
    } else if (contentList.at(i) == "hf1") {
      tfNtupleValues[i] = tf1_->tfMuonHwHF[tfMu1];
    } else if (contentList.at(i) == "hf2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonHwHF[tfMu2];
    } else if (contentList.at(i) == "phi1") {
      tfNtupleValues[i] = calcTFphi(tf1_->tfMuonHwPhi[tfMu1],
                                    tf1_->tfMuonTrackFinderType[tfMu1],
                                    tf1_->tfMuonProcessor[tfMu1]);
    } else if (contentList.at(i) == "phi2" && tfMu2 != -1) {
      tfNtupleValues[i] = calcTFphi(tf2_->tfMuonHwPhi[tfMu2],
                                    tf2_->tfMuonTrackFinderType[tfMu2],
                                    tf2_->tfMuonProcessor[tfMu2]);
    } else if (contentList.at(i) == "qual1") {
      tfNtupleValues[i] = tf1_->tfMuonHwQual[tfMu1];
    } else if (contentList.at(i) == "qual2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonHwQual[tfMu2];
    } else if (contentList.at(i) == "ch1") {
      tfNtupleValues[i] = std::pow(-1, tf1_->tfMuonHwSign[tfMu1]);
    } else if (contentList.at(i) == "ch2" && tfMu2 != -1) {
      tfNtupleValues[i] = std::pow(-1, tf2_->tfMuonHwSign[tfMu2]);
    } else if (contentList.at(i) == "tfType1") {
      tfNtupleValues[i] = tf1_->tfMuonTrackFinderType[tfMu1];
    } else if (contentList.at(i) == "tfType2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonTrackFinderType[tfMu2];
    } else if (contentList.at(i) == "tfProcessor1") {
      tfNtupleValues[i] = tf1_->tfMuonProcessor[tfMu1];
    } else if (contentList.at(i) == "tfProcessor2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonProcessor[tfMu2];
    } else {
      tfNtupleValues[i] = -10;
    }
  }
}
