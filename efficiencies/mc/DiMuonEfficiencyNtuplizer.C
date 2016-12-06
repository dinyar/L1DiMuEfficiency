#include "TFile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TTree.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <vector>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeTfMuonDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

enum tftype { bmtf, omtf, emtf, none };

struct coordinateTruth {
  float eta;
  float phi;
};

std::vector<std::string> getNtupleList(std::string fname);
double calcTFphi(int locPhi, tftype tfType, int proc);
double calcTFeta(int eta);
double calcTFpt(int pt);
void findUgmtMuons(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int& mu1,
                   int& mu2, const coordinateTruth& truthCoords1,
                   const coordinateTruth& truthCoords2);
void findTFMuons(L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* bmtf_,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* omtf_,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* emtf_,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat*& tf1,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat*& tf2,
                 int& mu1, int& mu2, const coordinateTruth& truthCoords1,
                 const coordinateTruth& truthCoords2);
std::vector<std::string> generateGenPhysicsQuantities();
std::vector<std::string> generateUgmtPhysicsQuantities();
std::vector<std::string> generateTfPhysicsQuantities();
std::vector<std::string> createContentList(
    int nGenMu, std::vector<std::string> genPhysicsQuantities,
    std::vector<std::string> l1PhysicsQuantities);
void fillNtuple(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_, int recoMu1, int recoMu2,
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
  int nevents;
  if (treeReco) {
    nevents = chainReco->GetEntries();
  } else {
    nevents = chainGen->GetEntries();
  }
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

  const int dir_err =
      mkdir(outDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (-1 == dir_err) {
    std::cout << "Error creating directory or directory exists already."
              << std::endl;
  }

  std::string gmtNtupleFname(outDir + "/uGMT" + ntupleName + ".root");
  std::string tfNtupleFname(outDir + "/tf" + ntupleName + ".root");

  TFile* ugmtFile = new TFile(gmtNtupleFname.c_str(), "RECREATE");

  std::vector<std::string> ugmtContentList = createContentList(
      nGenMu, generateGenPhysicsQuantities(), generateUgmtPhysicsQuantities());
  std::vector<std::string> tfContentList = createContentList(
      nGenMu, generateGenPhysicsQuantities(), generateTfPhysicsQuantities());

  std::ostringstream ugmtContentStream;
  std::copy(ugmtContentList.begin(), ugmtContentList.end() - 1,
            std::ostream_iterator<std::string>(ugmtContentStream, ":"));
  ugmtContentStream << *(ugmtContentList.rbegin());
  std::string ugmtContentStr(ugmtContentStream.str());

  std::ostringstream tfContentStream;
  std::copy(tfContentList.begin(), tfContentList.end() - 1,
            std::ostream_iterator<std::string>(tfContentStream, ":"));
  tfContentStream << *(tfContentList.rbegin());
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

    coordinateTruth truthCoords1;
    coordinateTruth truthCoords2;

    if (foundReco) {
      chainReco->GetEntry(jentry);
      // Check for correct number of reco muons
      if (reco_->nMuons != nGenMu) {
        continue;
      }

      int recoMu1 {0};
      int recoMu2 {-1};

      if((nGenMu > 1) && (reco_->muonEt[0] < reco_->muonEt[1])) {
        recoMu1 = 1;
        recoMu2 = 0;
      } else if(nGenMu > 1) {
        recoMu2 = 1;
      } 

      truthCoords1.eta = reco_->eta[recoMu1];
      truthCoords1.phi = reco_->phi[recoMu1];

      if (nGenMu > 1) {
        truthCoords2.eta = reco_->eta[recoMu2];
        truthCoords2.phi = reco_->phi[recoMu2];
      } else {
        truthCoords2.eta = -99999;
        truthCoords2.phi = -99999;
      }

      fillNtuple(reco_, recoMu1, recoMu2, ugmtContentList, ugmtNtupleValues);
      fillNtuple(reco_, recoMu1, recoMu2, tfContentList, tfNtupleValues);
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
      }

      if (nGenMus != nGenMu) {
        continue;
      }

      if((genMu1 > -1) && (genMu2 > -1) && (gen_->partPt[genMu1] < gen_->partPt[genMu2])) {
        int tmp {genMu1};
        genMu1 = genMu2;
        genMu2 = tmp;
      }

      truthCoords1.eta = gen_->partEta[genMu1];
      truthCoords1.phi = gen_->partPhi[genMu1];

      if (nGenMu > 1) {
        truthCoords2.eta = gen_->partEta[genMu2];
        truthCoords2.phi = gen_->partPhi[genMu2];
      } else {
        truthCoords2.eta = -99999;
        truthCoords2.phi = -99999;
      }

      fillNtuple(gen_, genMu1, genMu2, ugmtContentList, ugmtNtupleValues);
      fillNtuple(gen_, genMu1, genMu2, tfContentList, tfNtupleValues);
    }

    // Find two highest pT uGMT muons
    int ugmtMu1 = -1;
    int ugmtMu2 = -1;
    findUgmtMuons(ugmtMC_, ugmtMu1, ugmtMu2, truthCoords1, truthCoords2);

    // Find two highest pT TF muons
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1;
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2;
    int tfMu1 = -1;
    int tfMu2 = -1;
    findTFMuons(bmtfMC_, omtfMC_, emtfMC_, tf1, tf2, tfMu1, tfMu2, truthCoords1,
                truthCoords2);

    // Fill ugmt ntuple
    ugmtFile->cd();
    fillNtuple(ugmtMC_, ugmtMu1, ugmtMu2, ugmtContentList, ugmtNtupleValues);
    ugmtNtuple->Fill(ugmtNtupleValues);
    // Fill tf ntuple
    fillNtuple(tfMu1, tf1, tfMu2, tf2, tfContentList, tfNtupleValues);
    tfNtuple->Fill(tfNtupleValues);
  }

  // Write ntuple file to disk.
  ugmtFile->cd();
  ugmtFile->Write();
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

void findUgmtMuon(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int& mu,
                  const coordinateTruth& truthCoord, double& dEta,
                  const int veto) {
  float pT = 0;
  dEta = 999;
  for (int i = 0; i < ugmt_->nMuons; ++i) {
    if (i == veto) {
      continue;
    }
    if (ugmt_->muonBx[i] != 0) {
      continue;
    }
    if (ugmt_->muonQual[i] < 8) {
      continue;
    }
    double dEtaTmp = std::abs(ugmt_->muonEta[i] - truthCoord.eta);
    // Match L1 mu to truth in eta only if the truth value is sane (i.e. that
    // muon exists)
    if ((dEtaTmp > 0.3) && (std::abs(truthCoord.eta) < 3)) {
      continue;
    }

    if (ugmt_->muonEt[i] > pT) {
      pT = ugmt_->muonEt[i];
      mu = i;
      dEta = dEtaTmp;
    }
  }
}

void findUgmtMuons(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int& mu1,
                   int& mu2, const coordinateTruth& truthCoords1,
                   const coordinateTruth& truthCoords2) {
  double dEta1_tmp = 999;
  double dEta2_tmp = 999;

  int mu11 = -1;
  int mu12 = -1;

  findUgmtMuon(ugmt_, mu11, truthCoords1, dEta1_tmp, -1);
  findUgmtMuon(ugmt_, mu12, truthCoords2, dEta2_tmp, mu11);
  double dEta1 = dEta1_tmp + dEta2_tmp;

  dEta1_tmp = 999;
  dEta2_tmp = 999;

  int mu21 = -1;
  int mu22 = -1;

  findUgmtMuon(ugmt_, mu22, truthCoords2, dEta2_tmp, -1);
  findUgmtMuon(ugmt_, mu21, truthCoords1, dEta1_tmp, mu22);
  double dEta2 = dEta1_tmp + dEta2_tmp;

  if (dEta1 < dEta2) {
    mu1 = mu11;
    mu2 = mu12;
  } else {
    mu1 = mu21;
    mu2 = mu22;
  }

  if ((mu1 > -1) && (mu2 > -1) && (ugmt_->muonEt[mu1] < ugmt_->muonEt[mu2])) {
    int tmp{mu1};
    mu1 = mu2;
    mu2 = tmp;
  }
}

void findTfMuon(L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf_,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat*& tf, int& mu,
                float& pt, const coordinateTruth& truthCoord, double& dEta,
                const int veto) {
  for (int i = 0; i < tf_->nTfMuons; ++i) {
    if (i == veto) {
      continue;
    }
    if (tf_->tfMuonBx[i] != 0) {
      continue;
    }
    if (tf_->tfMuonHwQual[i] < 8) {
      continue;
    }
    double dEtaTmp = std::abs(calcTFeta(tf_->tfMuonHwEta[i]) - truthCoord.eta);
    // Match L1 mu to truth in eta only if the truth value is sane (i.e. that
    // muon exists)
    if ((dEtaTmp > 0.3) && (std::abs(truthCoord.eta) < 3)) {
      continue;
    }

    // Select highest pT muon.
    float pT_tmp = (tf_->tfMuonHwPt[i] - 1) * 0.5;
    if (pT_tmp > pt) {
      pt = pT_tmp;
      mu = i;
      tf = tf_;
      dEta = dEtaTmp;
    }
  }
}

void findTFMuons(L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* bmtf_,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* omtf_,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* emtf_,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat*& tf1,
                 L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat*& tf2,
                 int& mu1, int& mu2, const coordinateTruth& truthCoords1,
                 const coordinateTruth& truthCoords2) {
  double dEta1_tmp = 999;
  double dEta2_tmp = 999;
  float pt1 = 0;
  float pt2 = 0;
  int bmtf_veto = -1;
  int omtf_veto = -1;
  int emtf_veto = -1;

  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf11;
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf12;
  int mu11 = -1;
  int mu12 = -1;

  findTfMuon(bmtf_, tf11, mu11, pt1, truthCoords1, dEta1_tmp, -1);
  findTfMuon(omtf_, tf11, mu11, pt1, truthCoords1, dEta1_tmp, -1);
  findTfMuon(emtf_, tf11, mu11, pt1, truthCoords1, dEta1_tmp, -1);
  if (tf11 == bmtf_) {
    bmtf_veto = mu11;
  } else if (tf11 == omtf_) {
    omtf_veto = mu11;
  } else if (tf11 == emtf_) {
    emtf_veto = mu11;
  }
  findTfMuon(bmtf_, tf12, mu12, pt2, truthCoords2, dEta2_tmp, bmtf_veto);
  findTfMuon(omtf_, tf12, mu12, pt2, truthCoords2, dEta2_tmp, omtf_veto);
  findTfMuon(emtf_, tf12, mu12, pt2, truthCoords2, dEta2_tmp, emtf_veto);
  double dEta1 = dEta1_tmp + dEta2_tmp;

  dEta1_tmp = 999;
  dEta2_tmp = 999;
  pt1 = 0;
  pt2 = 0;
  bmtf_veto = -1;
  omtf_veto = -1;
  emtf_veto = -1;

  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf21;
  L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf22;
  int mu21 = -1;
  int mu22 = -1;

  findTfMuon(bmtf_, tf22, mu22, pt2, truthCoords2, dEta2_tmp, -1);
  findTfMuon(omtf_, tf22, mu22, pt2, truthCoords2, dEta2_tmp, -1);
  findTfMuon(emtf_, tf22, mu22, pt2, truthCoords2, dEta2_tmp, -1);
  if (tf22 == bmtf_) {
    bmtf_veto = mu22;
  } else if (tf22 == omtf_) {
    omtf_veto = mu22;
  } else if (tf22 == emtf_) {
    emtf_veto = mu22;
  }
  findTfMuon(bmtf_, tf21, mu21, pt1, truthCoords1, dEta1_tmp, bmtf_veto);
  findTfMuon(omtf_, tf21, mu21, pt1, truthCoords1, dEta1_tmp, omtf_veto);
  findTfMuon(emtf_, tf21, mu21, pt1, truthCoords1, dEta1_tmp, emtf_veto);
  double dEta2 = dEta1_tmp + dEta2_tmp;

  if (dEta1 < dEta2) {
    mu1 = mu11;
    mu2 = mu12;
    tf1 = tf11;
    tf2 = tf12;
  } else {
    mu1 = mu21;
    mu2 = mu22;
    tf1 = tf21;
    tf2 = tf22;
  }

  if ((mu1 > -1) && (mu2 > -1) &&
      (tf1->tfMuonHwPt[mu1] < tf2->tfMuonHwPt[mu2])) {
    int tmpMu{mu1};
    L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tmpTf{tf1};
    mu1 = mu2;
    tf1 = tf2;
    mu2 = tmpMu;
    tf2 = tmpTf;
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
  physicsQuantities.push_back("invMass");

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
    if (nGenMu > 1) {
      contentList.push_back(*it + "_jpsi");
    }
  }

  for (std::vector<std::string>::iterator it = genPhysicsQuantities.begin();
       it != genPhysicsQuantities.end(); ++it) {
    contentList.push_back(*it + "1_gen");
    if (nGenMu > 1) {
      contentList.push_back(*it + "2_gen");
      contentList.push_back(*it + "_jpsi");
    }
  }

  return contentList;
}

void fillNtuple(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_, int recoMu1, int recoMu2,
                std::vector<std::string> contentList, float ntupleValues[]) {
  TLorentzVector jPsi;
  if (reco_->nMuons > 1) {
    TLorentzVector mu1;
    TLorentzVector mu2;
    mu1.SetPtEtaPhiE(reco_->pt[recoMu1], reco_->eta[recoMu1], reco_->phi[recoMu1], reco_->e[recoMu1]);
    mu2.SetPtEtaPhiE(reco_->pt[recoMu2], reco_->eta[recoMu2], reco_->phi[recoMu2], reco_->e[recoMu2]);
    jPsi = mu1 + mu2;
  }
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1_gen") {
      ntupleValues[i] = reco_->pt[recoMu1];
    } else if (contentList.at(i) == "pT2_gen") {
      ntupleValues[i] = reco_->pt[recoMu2];
    } else if (contentList.at(i) == "eta1_gen") {
      ntupleValues[i] = reco_->eta[recoMu1];
    } else if (contentList.at(i) == "eta2_gen") {
      ntupleValues[i] = reco_->eta[recoMu2];
    } else if (contentList.at(i) == "phi1_gen") {
      ntupleValues[i] = reco_->phi[recoMu1];
    } else if (contentList.at(i) == "phi2_gen") {
      ntupleValues[i] = reco_->phi[recoMu2];
    } else if (contentList.at(i) == "ch1_gen") {
      ntupleValues[i] = reco_->charge[recoMu1];
    } else if (contentList.at(i) == "ch2_gen") {
      ntupleValues[i] = reco_->charge[recoMu2];
    } else if (contentList.at(i) == "pT_jpsi") {
      ntupleValues[i] = jPsi.Pt();
    } else if (contentList.at(i) == "eta_jpsi") {
      ntupleValues[i] = jPsi.Eta();
    } else if (contentList.at(i) == "phi_jpsi") {
      ntupleValues[i] = jPsi.Phi();
    } else {
      ntupleValues[i] = -10;
    }
  }
}

void fillNtuple(L1Analysis::L1AnalysisGeneratorDataFormat* gen_, int genMu1,
                int genMu2, std::vector<std::string> contentList,
                float ntupleValues[]) {
  TLorentzVector jPsi;
  if (genMu2 > -1) {
    TLorentzVector mu1;
    TLorentzVector mu2;
    mu1.SetPtEtaPhiE(gen_->partPt[genMu1], gen_->partEta[genMu1],
                     gen_->partPhi[genMu1], gen_->partE[genMu1]);
    mu2.SetPtEtaPhiE(gen_->partPt[genMu2], gen_->partEta[genMu2],
                     gen_->partPhi[genMu2], gen_->partE[genMu2]);
    jPsi = mu1 + mu2;
  }
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1_gen") {
      ntupleValues[i] = gen_->partPt[genMu1];
    } else if (contentList.at(i) == "pT2_gen") {
      ntupleValues[i] = gen_->partPt[genMu2];
    } else if (contentList.at(i) == "eta1_gen") {
      ntupleValues[i] = gen_->partEta[genMu1];
    } else if (contentList.at(i) == "eta2_gen") {
      ntupleValues[i] = gen_->partEta[genMu2];
    } else if (contentList.at(i) == "phi1_gen") {
      ntupleValues[i] = gen_->partPhi[genMu1];
    } else if (contentList.at(i) == "phi2_gen") {
      ntupleValues[i] = gen_->partPhi[genMu2];
    } else if (contentList.at(i) == "ch1_gen") {
      ntupleValues[i] = gen_->partCh[0];
    } else if (contentList.at(i) == "ch2_gen") {
      ntupleValues[i] = gen_->partCh[1];
    } else if (contentList.at(i) == "pT_jpsi") {
      ntupleValues[i] = jPsi.Pt();
    } else if (contentList.at(i) == "eta_jpsi") {
      ntupleValues[i] = jPsi.Eta();
    } else if (contentList.at(i) == "phi_jpsi") {
      ntupleValues[i] = jPsi.Phi();
    } else {
      ntupleValues[i] = -10;
    }
  }
}

void fillNtuple(L1Analysis::L1AnalysisL1UpgradeDataFormat* ugmt_, int ugmtMu1,
                int ugmtMu2, std::vector<std::string> contentList,
                float ugmtNtupleValues[]) {
  TLorentzVector jPsi;
  if (ugmtMu2 > -1) {
    TLorentzVector mu1;
    TLorentzVector mu2;
    mu1.SetPtEtaPhiM(ugmt_->muonEt[ugmtMu1], ugmt_->muonEta[ugmtMu1],
                     ugmt_->muonPhi[ugmtMu1], 0.105);
    mu2.SetPtEtaPhiM(ugmt_->muonEt[ugmtMu2], ugmt_->muonEta[ugmtMu2],
                     ugmt_->muonPhi[ugmtMu2], 0.105);
    jPsi = mu1 + mu2;
  }

  // TODO: Missing track addresses, HF bit, processor
  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1" && ugmtMu1 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonEt[ugmtMu1];
    } else if (contentList.at(i) == "pT2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonEt[ugmtMu2];
    } else if (contentList.at(i) == "eta1" && ugmtMu1 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonEta[ugmtMu1];
    } else if (contentList.at(i) == "eta2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonEta[ugmtMu2];
    } else if (contentList.at(i) == "phi1" && ugmtMu1 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonPhi[ugmtMu1];
    } else if (contentList.at(i) == "phi2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonPhi[ugmtMu2];
    } else if (contentList.at(i) == "qual1" && ugmtMu1 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonQual[ugmtMu1];
    } else if (contentList.at(i) == "qual2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonQual[ugmtMu2];
    } else if (contentList.at(i) == "ch1" && ugmtMu1 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonChg[ugmtMu1];
    } else if (contentList.at(i) == "ch2" && ugmtMu2 != -1) {
      ugmtNtupleValues[i] = ugmt_->muonChg[ugmtMu2];
    } else if (contentList.at(i) == "tfType1" && ugmtMu1 != -1) {
      int muIdx = ugmt_->muonTfMuonIdx[ugmtMu1];
      if (muIdx < 18 || muIdx >= 90) {
        ugmtNtupleValues[i] = 2;  // EMTF
      } else if (muIdx < 36 || muIdx >= 72) {
        ugmtNtupleValues[i] = 1;  // OMTF
      } else if (muIdx >= 36 && muIdx < 72) {
        ugmtNtupleValues[i] = 0;  // BMTF
      }
    } else if (contentList.at(i) == "tfType2" && ugmtMu2 != -1) {
      int muIdx = ugmt_->muonTfMuonIdx[ugmtMu2];
      if (muIdx < 18 || muIdx >= 90) {
        ugmtNtupleValues[i] = 2;  // EMTF
      } else if (muIdx < 36 || muIdx >= 72) {
        ugmtNtupleValues[i] = 1;  // OMTF
      } else if (muIdx >= 36 && muIdx < 72) {
        ugmtNtupleValues[i] = 0;  // BMTF
      }
    }  // Don't need a catch-all case as we're calling the fill function for
       // gen/reco where all "non -filled" fields are initialized to -10
       // already.
  }
}

void fillNtuple(int tfMu1,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf1_,
                int tfMu2,
                L1Analysis::L1AnalysisL1UpgradeTfMuonDataFormat* tf2_,
                std::vector<std::string> contentList, float tfNtupleValues[]) {
  tftype tfType1 = tftype::none;
  if (tfMu1 != -1) {
    if (tf1_->tfMuonTrackFinderType[tfMu1] == 0) {
      tfType1 = tftype::bmtf;
    } else if (tf1_->tfMuonTrackFinderType[tfMu1] == 1 ||
               tf1_->tfMuonTrackFinderType[tfMu1] == 2) {
      tfType1 = tftype::omtf;
    } else if (tf1_->tfMuonTrackFinderType[tfMu1] == 3 ||
               tf1_->tfMuonTrackFinderType[tfMu1] == 4) {
      tfType1 = tftype::emtf;
    }
  }

  tftype tfType2 = tftype::none;
  if (tfMu2 != -1) {
    if (tf2_->tfMuonTrackFinderType[tfMu2] == 0) {
      tfType2 = tftype::bmtf;
    } else if (tf2_->tfMuonTrackFinderType[tfMu2] == 1 ||
               tf2_->tfMuonTrackFinderType[tfMu2] == 2) {
      tfType2 = tftype::omtf;
    } else if (tf2_->tfMuonTrackFinderType[tfMu2] == 3 ||
               tf2_->tfMuonTrackFinderType[tfMu2] == 4) {
      tfType2 = tftype::emtf;
    }
  }

  for (int i = 0; i < contentList.size(); ++i) {
    if (contentList.at(i) == "pT1" && tfMu1 != -1) {
      tfNtupleValues[i] = calcTFpt(tf1_->tfMuonHwPt[tfMu1]);
    } else if (contentList.at(i) == "pT2" && tfMu2 != -1) {
      tfNtupleValues[i] = calcTFpt(tf2_->tfMuonHwPt[tfMu2]);
    } else if (contentList.at(i) == "eta1" && tfMu1 != -1) {
      tfNtupleValues[i] = calcTFeta(tf1_->tfMuonHwEta[tfMu1]);
    } else if (contentList.at(i) == "eta2" && tfMu2 != -1) {
      tfNtupleValues[i] = calcTFeta(tf2_->tfMuonHwEta[tfMu2]);
    } else if (contentList.at(i) == "hf1" && tfMu1 != -1) {
      tfNtupleValues[i] = tf1_->tfMuonHwHF[tfMu1];
    } else if (contentList.at(i) == "hf2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonHwHF[tfMu2];
    } else if (contentList.at(i) == "phi1" && tfMu1 != -1) {
      tfNtupleValues[i] = calcTFphi(tf1_->tfMuonHwPhi[tfMu1], tfType1,
                                    tf1_->tfMuonProcessor[tfMu1]);
    } else if (contentList.at(i) == "phi2" && tfMu2 != -1) {
      tfNtupleValues[i] = calcTFphi(tf2_->tfMuonHwPhi[tfMu2], tfType2,
                                    tf2_->tfMuonProcessor[tfMu2]);
    } else if (contentList.at(i) == "qual1" && tfMu1 != -1) {
      tfNtupleValues[i] = tf1_->tfMuonHwQual[tfMu1];
    } else if (contentList.at(i) == "qual2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonHwQual[tfMu2];
    } else if (contentList.at(i) == "ch1" && tfMu1 != -1) {
      tfNtupleValues[i] = std::pow(-1, tf1_->tfMuonHwSign[tfMu1]);
    } else if (contentList.at(i) == "ch2" && tfMu2 != -1) {
      tfNtupleValues[i] = std::pow(-1, tf2_->tfMuonHwSign[tfMu2]);
    } else if (contentList.at(i) == "tfType1" && tfMu1 != -1) {
      tfNtupleValues[i] = tfType1;
    } else if (contentList.at(i) == "tfType2" && tfMu2 != -1) {
      tfNtupleValues[i] = tfType2;
    } else if (contentList.at(i) == "tfProcessor1" && tfMu1 != -1) {
      tfNtupleValues[i] = tf1_->tfMuonProcessor[tfMu1];
    } else if (contentList.at(i) == "tfProcessor2" && tfMu2 != -1) {
      tfNtupleValues[i] = tf2_->tfMuonProcessor[tfMu2];
    }  // Don't need a catch-all case as we're calling the fill function for
       // gen/reco where all "non-filled" fields are initialized to -10 already.
  }
}
