#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "L1Trigger/L1TNtuples/interface/L1AnalysisGeneratorDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"

const std::string unpackTreepath("l1UpgradeTree/L1UpgradeTree");
const std::string recoTreepath("l1MuonRecoTree/Muon2RecoTree");
const std::string genTreepath("l1GeneratorTree/L1GenTree");
// mu bins
const int nMuBins = 100;
const float muLo = 0;
const float muHi = 100;
// Tag and probe constants
const int tagPt = 27;

bool readFList(std::string fname, std::vector<std::string>& listNtuples);
int setupTChain(const std::vector<std::string> listNtuples, TChain* unpackChain,
                TChain* recoChain);
void getSingleMuDataEfficiency(int nentries, TChain* l1Chain, TChain* recoChain,
                               const int pTcut, const double etaLow,
                               const double etaHigh, TH1D& effHist,
                               TGraphAsymmErrors& effErrors);
void getSingleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pTcut, const double etaLow,
                             const double etaHigh, TH1D& effHist,
                             TGraphAsymmErrors& effErrors);
void getDoubleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pT1cut, const int pT2cut,
                             const double etaLow, const double etaHigh,
                             TH1D& effHist, TGraphAsymmErrors& effErrors);
bool findGenMuon(L1Analysis::L1AnalysisGeneratorDataFormat* gen_, int& mu1);
bool findGenMuon(L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                 const int nMus, int& mu1, int& mu2);
double recoDist(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_, const int mu1,
                const int mu2);
double matchL1toGen(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                    L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                    const int l1Mu, const int genMu);
bool findBestGenMatch(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                      L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                      int& l1mu, const int pTcut, const int genMu,
                      const double dRcut);
bool findBestGenMatches(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                        L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                        int& l1mu1, int& l1mu2, const int pT1cut,
                        const int pT2cut, const int genMu1, const int genMu2,
                        const double dRcut);
double matchL1toReco(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                     L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                     const int l1Mu, const int recoMu);
bool findBestRecoMatch(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                       L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                       int& l1mu, const int pTcut, const int recoMu,
                       const double dRcut);
// TODO:
// void DrawHistogram(const TH1D& hist, const TGraphAsymmErrors& err,
//                    const std::string name);
void DrawHistograms(std::vector<TH1D>& hists, const std::vector<int> colours,
                    const std::vector<int> markers,
                    const std::vector<std::string>& histnames,
                    const std::string& type, const std::string& name);
void DrawHistograms(std::vector<TH1D>& hists, const std::vector<int> colours,
                    const std::vector<int> markers,
                    const std::vector<std::string>& histnames,
                    const std::string& type,
                    std::vector<TGraphAsymmErrors>& errs,
                    const std::string& name);

void diMuEfficiency(std::string singleMuDataFile, std::string singleMuMcFile,
                    std::string diMuMcFile, std::string folder, int mu1cut = 2,
                    int mu2cut = 2) {
  std::string plotFolder = "plots/" + folder + "/";
  std::cout << "Creating directory: " << plotFolder << std::endl;
  const int dir_err =
      mkdir(plotFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (dir_err == -1) {
    std::cout << "Error creating directory or directory exists already."
              << std::endl;
    return;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  std::vector<std::string> listSingleDataNtuples;
  std::vector<std::string> listSingleMcNtuples;
  std::vector<std::string> listDoubleMcNtuples;

  bool success = readFList(singleMuDataFile, listSingleDataNtuples);
  success &= readFList(singleMuMcFile, listSingleMcNtuples);
  success &= readFList(diMuMcFile, listDoubleMcNtuples);

  if (!success) {
    std::cout << "Exiting.. " << std::endl;
    return;
  }

  // OpenWithoutInit
  TChain* l1SingleDataChain = new TChain(unpackTreepath.c_str());
  TChain* recoSingleDataChain = new TChain(recoTreepath.c_str());
  TChain* l1SingleMcChain = new TChain(unpackTreepath.c_str());
  TChain* genSingleMcChain = new TChain(genTreepath.c_str());
  TChain* l1DoubleMcChain = new TChain(unpackTreepath.c_str());
  TChain* genDoubleMcChain = new TChain(genTreepath.c_str());
  int singleDataEntries = setupTChain(listSingleDataNtuples, l1SingleDataChain,
                                      recoSingleDataChain);
  int singleMcEntries =
      setupTChain(listSingleMcNtuples, l1SingleMcChain, genSingleMcChain);
  int doubleMcEntries =
      setupTChain(listDoubleMcNtuples, l1DoubleMcChain, genDoubleMcChain);

  // make histos
  // BMTF
  const double bmtfLow = 0;
  const double bmtfHigh = 0.7;
  TH1D bmtfSingleMuDataEfficiency("bmtfSingleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D bmtfSingleMuMcEfficiency("bmtfSingleMuMcEfficiency", "", nMuBins,
                                muLo - 0.1, muHi + 0.1);
  TH1D bmtfDoubleMuDataEfficiency("bmtfDoubleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D bmtfDoubleMuMcEfficiency("bmtfDoubleMuMcEfficiency", "", nMuBins,
                                muLo - 0.1, muHi + 0.1);
  // BOMTF
  const double bomtfLow = 0.7;
  const double bomtfHigh = 0.9;
  TH1D bomtfSingleMuDataEfficiency("bomtfSingleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D bomtfSingleMuMcEfficiency("bomtfSingleMuMcEfficiency", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  TH1D bomtfDoubleMuDataEfficiency("bomtfDoubleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D bomtfDoubleMuMcEfficiency("bomtfDoubleMuMcEfficiency", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  // OMTF
  const double omtfLow = 0.9;
  const double omtfHigh = 1.15;
  TH1D omtfSingleMuDataEfficiency("omtfSingleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D omtfSingleMuMcEfficiency("omtfSingleMuMcEfficiency", "", nMuBins,
                                muLo - 0.1, muHi + 0.1);
  TH1D omtfDoubleMuDataEfficiency("omtfDoubleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D omtfDoubleMuMcEfficiency("omtfDoubleMuMcEfficiency", "", nMuBins,
                                muLo - 0.1, muHi + 0.1);
  // EOMTF
  const double eomtfLow = 1.15;
  const double eomtfHigh = 1.35;
  TH1D eomtfSingleMuDataEfficiency("eomtfSingleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D eomtfSingleMuMcEfficiency("eomtfSingleMuMcEfficiency", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  TH1D eomtfDoubleMuDataEfficiency("eomtfDoubleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D eomtfDoubleMuMcEfficiency("eomtfDoubleMuMcEfficiency", "", nMuBins,
                                 muLo - 0.1, muHi + 0.1);
  // EMTF
  const double emtfLow = 1.35;
  const double emtfHigh = 2.5;
  TH1D emtfSingleMuDataEfficiency("emtfSingleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D emtfSingleMuMcEfficiency("emtfSingleMuMcEfficiency", "", nMuBins,
                                muLo - 0.1, muHi + 0.1);
  TH1D emtfDoubleMuDataEfficiency("emtfDoubleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D emtfDoubleMuMcEfficiency("emtfDoubleMuMcEfficiency", "", nMuBins,
                                muLo - 0.1, muHi + 0.1);
  // For correct error bars
  TGraphAsymmErrors bmtfSingleMuDataErrors = TGraphAsymmErrors();
  TGraphAsymmErrors bmtfSingleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors bmtfDoubleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors bomtfSingleMuDataErrors = TGraphAsymmErrors();
  TGraphAsymmErrors bomtfSingleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors bomtfDoubleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors omtfSingleMuDataErrors = TGraphAsymmErrors();
  TGraphAsymmErrors omtfSingleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors omtfDoubleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors eomtfSingleMuDataErrors = TGraphAsymmErrors();
  TGraphAsymmErrors eomtfSingleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors eomtfDoubleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors emtfSingleMuDataErrors = TGraphAsymmErrors();
  TGraphAsymmErrors emtfSingleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors emtfDoubleMuMcErrors = TGraphAsymmErrors();

  getSingleMuDataEfficiency(singleDataEntries, l1SingleDataChain,
                            recoSingleDataChain, mu1cut, bmtfLow, bmtfHigh,
                            bmtfSingleMuDataEfficiency, bmtfSingleMuDataErrors);
  getSingleMuDataEfficiency(singleDataEntries, l1SingleDataChain,
                            recoSingleDataChain, mu1cut, bomtfLow, bomtfHigh,
                            bomtfSingleMuDataEfficiency,
                            bomtfSingleMuDataErrors);
  getSingleMuDataEfficiency(singleDataEntries, l1SingleDataChain,
                            recoSingleDataChain, mu1cut, omtfLow, omtfHigh,
                            omtfSingleMuDataEfficiency, omtfSingleMuDataErrors);
  getSingleMuDataEfficiency(singleDataEntries, l1SingleDataChain,
                            recoSingleDataChain, mu1cut, eomtfLow, eomtfHigh,
                            eomtfSingleMuDataEfficiency,
                            eomtfSingleMuDataErrors);
  getSingleMuDataEfficiency(singleDataEntries, l1SingleDataChain,
                            recoSingleDataChain, mu1cut, emtfLow, emtfHigh,
                            emtfSingleMuDataEfficiency, emtfSingleMuDataErrors);

  getSingleMuMcEfficiency(singleMcEntries, l1SingleMcChain, genSingleMcChain,
                          mu1cut, bmtfLow, bmtfHigh, bmtfSingleMuMcEfficiency,
                          bmtfSingleMuMcErrors);
  getSingleMuMcEfficiency(singleMcEntries, l1SingleMcChain, genSingleMcChain,
                          mu1cut, bomtfLow, bomtfHigh,
                          bomtfSingleMuMcEfficiency, bomtfSingleMuMcErrors);
  getSingleMuMcEfficiency(singleMcEntries, l1SingleMcChain, genSingleMcChain,
                          mu1cut, omtfLow, omtfHigh, omtfSingleMuMcEfficiency,
                          omtfSingleMuMcErrors);
  getSingleMuMcEfficiency(singleMcEntries, l1SingleMcChain, genSingleMcChain,
                          mu1cut, eomtfLow, eomtfHigh,
                          eomtfSingleMuMcEfficiency, eomtfSingleMuMcErrors);
  getSingleMuMcEfficiency(singleMcEntries, l1SingleMcChain, genSingleMcChain,
                          mu1cut, emtfLow, emtfHigh, emtfSingleMuMcEfficiency,
                          emtfSingleMuMcErrors);

  getDoubleMuMcEfficiency(doubleMcEntries, l1DoubleMcChain, genDoubleMcChain,
                          mu1cut, mu2cut, bmtfLow, bmtfHigh,
                          bmtfDoubleMuMcEfficiency, bmtfDoubleMuMcErrors);
  getDoubleMuMcEfficiency(doubleMcEntries, l1DoubleMcChain, genDoubleMcChain,
                          mu1cut, mu2cut, bomtfLow, bomtfHigh,
                          bomtfDoubleMuMcEfficiency, bomtfDoubleMuMcErrors);
  getDoubleMuMcEfficiency(doubleMcEntries, l1DoubleMcChain, genDoubleMcChain,
                          mu1cut, mu2cut, omtfLow, omtfHigh,
                          omtfDoubleMuMcEfficiency, omtfDoubleMuMcErrors);
  getDoubleMuMcEfficiency(doubleMcEntries, l1DoubleMcChain, genDoubleMcChain,
                          mu1cut, mu2cut, eomtfLow, eomtfHigh,
                          eomtfDoubleMuMcEfficiency, eomtfDoubleMuMcErrors);
  getDoubleMuMcEfficiency(doubleMcEntries, l1DoubleMcChain, genDoubleMcChain,
                          mu1cut, mu2cut, emtfLow, emtfHigh,
                          emtfDoubleMuMcEfficiency, emtfDoubleMuMcErrors);

  std::vector<std::string> regionNames;
  std::ostringstream oss1, oss2, oss3, oss4, oss5;
  oss1 << bmtfLow << " #leq |#eta| < " << bmtfHigh;
  oss2 << bomtfLow << " #leq |#eta| < " << bomtfHigh;
  oss3 << omtfLow << " #leq |#eta| < " << omtfHigh;
  oss4 << eomtfLow << " #leq |#eta| < " << eomtfHigh;
  oss5 << emtfLow << " #leq |#eta| < " << emtfHigh;
  regionNames.push_back(oss1.str());
  regionNames.push_back(oss2.str());
  regionNames.push_back(oss3.str());
  regionNames.push_back(oss4.str());
  regionNames.push_back(oss5.str());
  std::vector<int> colours;
  colours.push_back(41);
  colours.push_back(9);
  colours.push_back(36);
  colours.push_back(3);
  colours.push_back(48);
  std::vector<int> markers;
  markers.push_back(24);
  markers.push_back(25);
  markers.push_back(26);
  markers.push_back(27);
  markers.push_back(32);

  // Drawing the single mu efficiencies now as we're squaring the histograms
  // later.
  std::vector<TH1D> singleMuMcEffs;
  singleMuMcEffs.push_back(bmtfSingleMuMcEfficiency);
  singleMuMcEffs.push_back(eomtfSingleMuMcEfficiency);
  singleMuMcEffs.push_back(omtfSingleMuMcEfficiency);
  singleMuMcEffs.push_back(eomtfSingleMuMcEfficiency);
  singleMuMcEffs.push_back(emtfSingleMuMcEfficiency);
  std::vector<TGraphAsymmErrors> singleMuMcErrs;
  singleMuMcErrs.push_back(bmtfSingleMuMcErrors);
  singleMuMcErrs.push_back(bomtfSingleMuMcErrors);
  singleMuMcErrs.push_back(omtfSingleMuMcErrors);
  singleMuMcErrs.push_back(eomtfSingleMuMcErrors);
  singleMuMcErrs.push_back(emtfSingleMuMcErrors);

  DrawHistograms(singleMuMcEffs, colours, markers, regionNames,
                 "L1T Efficiency", singleMuMcErrs,
                 plotFolder + "singleMuonEfficiencies_MC");

  std::vector<TH1D> singleMuDataEffs;
  singleMuDataEffs.push_back(bmtfSingleMuDataEfficiency);
  singleMuDataEffs.push_back(eomtfSingleMuDataEfficiency);
  singleMuDataEffs.push_back(omtfSingleMuDataEfficiency);
  singleMuDataEffs.push_back(eomtfSingleMuDataEfficiency);
  singleMuDataEffs.push_back(emtfSingleMuDataEfficiency);
  std::vector<TGraphAsymmErrors> singleMuDataErrs;
  singleMuDataErrs.push_back(bmtfSingleMuDataErrors);
  singleMuDataErrs.push_back(bomtfSingleMuDataErrors);
  singleMuDataErrs.push_back(omtfSingleMuDataErrors);
  singleMuDataErrs.push_back(eomtfSingleMuDataErrors);
  singleMuDataErrs.push_back(emtfSingleMuDataErrors);

  DrawHistograms(singleMuDataEffs, colours, markers, regionNames,
                 "L1T Efficiency", singleMuDataErrs,
                 plotFolder + "singleMuonEfficiencies_Data");

  std::vector<TH1D> doubleMuMcEffs;
  doubleMuMcEffs.push_back(bmtfDoubleMuMcEfficiency);
  doubleMuMcEffs.push_back(eomtfDoubleMuMcEfficiency);
  doubleMuMcEffs.push_back(omtfDoubleMuMcEfficiency);
  doubleMuMcEffs.push_back(eomtfDoubleMuMcEfficiency);
  doubleMuMcEffs.push_back(emtfDoubleMuMcEfficiency);
  std::vector<TGraphAsymmErrors> doubleMuMcErrs;
  doubleMuMcErrs.push_back(bmtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(bomtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(omtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(eomtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(emtfDoubleMuMcErrors);

  DrawHistograms(doubleMuMcEffs, colours, markers, regionNames,
                 "L1T Efficiency", doubleMuMcErrs,
                 plotFolder + "doubleMuonEfficiencies_MC");

  TH1D bmtfRhoFactor("bmtfRhoFactor", "", nMuBins, muLo - 0.1, muHi + 0.1);
  bmtfRhoFactor.Sumw2();
  TH1D bomtfRhoFactor("bomtfRhoFactor", "", nMuBins, muLo - 0.1, muHi + 0.1);
  bomtfRhoFactor.Sumw2();
  TH1D omtfRhoFactor("omtfRhoFactor", "", nMuBins, muLo - 0.1, muHi + 0.1);
  omtfRhoFactor.Sumw2();
  TH1D eomtfRhoFactor("eomtfRhoFactor", "", nMuBins, muLo - 0.1, muHi + 0.1);
  eomtfRhoFactor.Sumw2();
  TH1D emtfRhoFactor("emtfRhoFactor", "", nMuBins, muLo - 0.1, muHi + 0.1);
  emtfRhoFactor.Sumw2();

  // Squaring single mu efficiencies to get "naive" double mu efficiencies.
  bmtfSingleMuMcEfficiency.Multiply(&bmtfSingleMuMcEfficiency);
  bomtfSingleMuMcEfficiency.Multiply(&bomtfSingleMuMcEfficiency);
  omtfSingleMuMcEfficiency.Multiply(&omtfSingleMuMcEfficiency);
  eomtfSingleMuMcEfficiency.Multiply(&eomtfSingleMuMcEfficiency);
  emtfSingleMuMcEfficiency.Multiply(&emtfSingleMuMcEfficiency);
  bmtfSingleMuDataEfficiency.Multiply(&bmtfSingleMuDataEfficiency);
  bomtfSingleMuDataEfficiency.Multiply(&bomtfSingleMuDataEfficiency);
  omtfSingleMuDataEfficiency.Multiply(&omtfSingleMuDataEfficiency);
  eomtfSingleMuDataEfficiency.Multiply(&eomtfSingleMuDataEfficiency);
  emtfSingleMuDataEfficiency.Multiply(&emtfSingleMuDataEfficiency);

  bmtfRhoFactor.Divide(&bmtfDoubleMuMcEfficiency, &bmtfSingleMuMcEfficiency, 1,
                       1,
                       "");  // Two different datasets, no binomial errors.
  bomtfRhoFactor.Divide(&bomtfDoubleMuMcEfficiency, &bomtfSingleMuMcEfficiency,
                        1, 1,
                        "");  // Two different datasets, no binomial errors.
  omtfRhoFactor.Divide(&omtfDoubleMuMcEfficiency, &omtfSingleMuMcEfficiency, 1,
                       1,
                       "");  // Two different datasets, no binomial errors.
  eomtfRhoFactor.Divide(&eomtfDoubleMuMcEfficiency, &eomtfSingleMuMcEfficiency,
                        1, 1,
                        "");  // Two different datasets, no binomial errors.
  emtfRhoFactor.Divide(&emtfDoubleMuMcEfficiency, &emtfSingleMuMcEfficiency, 1,
                       1,
                       "");  // Two different datasets, no binomial errors.

  bmtfDoubleMuDataEfficiency.Multiply(&bmtfSingleMuDataEfficiency,
                                      &bmtfRhoFactor);
  bomtfDoubleMuDataEfficiency.Multiply(&bomtfSingleMuDataEfficiency,
                                       &bomtfRhoFactor);
  omtfDoubleMuDataEfficiency.Multiply(&omtfSingleMuDataEfficiency,
                                      &omtfRhoFactor);
  eomtfDoubleMuDataEfficiency.Multiply(&eomtfSingleMuDataEfficiency,
                                       &eomtfRhoFactor);
  emtfDoubleMuDataEfficiency.Multiply(&emtfSingleMuDataEfficiency,
                                      &emtfRhoFactor);

  std::vector<TH1D> naiveDoubleMuMcEffs;
  naiveDoubleMuMcEffs.push_back(bmtfSingleMuMcEfficiency);
  naiveDoubleMuMcEffs.push_back(eomtfSingleMuMcEfficiency);
  naiveDoubleMuMcEffs.push_back(omtfSingleMuMcEfficiency);
  naiveDoubleMuMcEffs.push_back(eomtfSingleMuMcEfficiency);
  naiveDoubleMuMcEffs.push_back(emtfSingleMuMcEfficiency);

  DrawHistograms(naiveDoubleMuMcEffs, colours, markers, regionNames,
                 "L1T Efficiency",
                 plotFolder + "naiveDoubleMuonEfficiencies_MC");

  std::vector<TH1D> naiveDoubleMuDataEffs;
  naiveDoubleMuDataEffs.push_back(bmtfSingleMuDataEfficiency);
  naiveDoubleMuDataEffs.push_back(eomtfSingleMuDataEfficiency);
  naiveDoubleMuDataEffs.push_back(omtfSingleMuDataEfficiency);
  naiveDoubleMuDataEffs.push_back(eomtfSingleMuDataEfficiency);
  naiveDoubleMuDataEffs.push_back(emtfSingleMuDataEfficiency);

  DrawHistograms(naiveDoubleMuDataEffs, colours, markers, regionNames,
                 "L1T Efficiency",
                 plotFolder + "naiveDoubleMuonEfficiencies_Data");

  std::vector<TH1D> rhoFactors;
  rhoFactors.push_back(bmtfRhoFactor);
  rhoFactors.push_back(bomtfRhoFactor);
  rhoFactors.push_back(omtfRhoFactor);
  rhoFactors.push_back(eomtfRhoFactor);
  rhoFactors.push_back(emtfRhoFactor);

  DrawHistograms(rhoFactors, colours, markers, regionNames, "#rho",
                 plotFolder + "rhoFactors");

  std::vector<TH1D> doubleMuDataEffs;
  doubleMuDataEffs.push_back(bmtfDoubleMuDataEfficiency);
  doubleMuDataEffs.push_back(eomtfDoubleMuDataEfficiency);
  doubleMuDataEffs.push_back(omtfDoubleMuDataEfficiency);
  doubleMuDataEffs.push_back(eomtfDoubleMuDataEfficiency);
  doubleMuDataEffs.push_back(emtfDoubleMuDataEfficiency);

  DrawHistograms(doubleMuDataEffs, colours, markers, regionNames,
                 "L1T Efficiency", plotFolder + "doubleMuonEfficiencies_Data");
}

bool readFList(std::string fname, std::vector<std::string>& listNtuples) {
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

  TTree* treeL1Unpack = (TTree*)rf->Get(unpackTreepath.c_str());
  TTree* treeReco = (TTree*)rf->Get(recoTreepath.c_str());
  TTree* treeGen = (TTree*)rf->Get(genTreepath.c_str());

  if (!treeL1Unpack) {
    std::cout << "L1Upgrade trees not found.. " << std::endl;
    return false;
  }
  if (!treeReco && !treeGen) {
    std::cout << "No truth tree found.. " << std::endl;
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

void getSingleMuDataEfficiency(int nentries, TChain* l1Chain, TChain* recoChain,
                               const int pTcut, const double etaLow,
                               const double etaHigh, TH1D& effHist,
                               TGraphAsymmErrors& effErrors) {
  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* l1_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_ =
      new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  l1Chain->SetBranchAddress("L1Upgrade", &l1_);
  recoChain->SetBranchAddress("Muon", &reco_);

  // Method:
  // Event with at least two reco muons.
  // Find one L1 muon (matched to one of the reco muons) which is also a tight
  // muon
  // Now find second muon that passes your selection criteria
  // NOTE: Tag can be probe too.

  TH1D allEventsHist =
      TH1D("allEventsHist", "", nMuBins, muLo - 0.1, muHi + 0.1);

  effHist.Sumw2();
  allEventsHist.Sumw2();

  std::cout << "Running over " << nentries << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if ((jentry % 1000) == 0) {
      std::cout << "Done " << jentry << " events..." << std::endl;
    }

    l1Chain->GetEntry(jentry);
    recoChain->GetEntry(jentry);

    int nRecoMus(0);
    for (int i = 0; i < reco_->nMuons; ++i) {
      if (reco_->isTightMuon[i]) {
        ++nRecoMus;
      }
    }
    if (nRecoMus < 2) {
      continue;
    }

    std::vector<int> probeMus;

    for (int i = 0; i < reco_->nMuons; ++i) {
      if (!reco_->isTightMuon[i] || (reco_->pt[i] < tagPt)) {
        continue;
      }
      int tagMu(i);
      int l1mu1(-1);
      if (!findBestRecoMatch(l1_, reco_, l1mu1, pTcut, i, 0.5)) {
        continue;
      }
      // TODO: Possibly require that outside of eta region?
      for (int j = 0; j < reco_->nMuons; ++j) {
        if (!reco_->isTightMuon[j] || (tagMu == j) ||
            (std::find(probeMus.begin(), probeMus.end(), j) !=
             probeMus.end())) {
          continue;
        }
        // Check that probe and tag are more than dR=0.5 separated.
        if (recoDist(reco_, i, j) <= 0.5) {
          continue;
        }
        // TODO: Possibly check invariant mass.
        probeMus.push_back(j);
        if ((abs(reco_->eta[j]) < etaLow) || (abs(reco_->eta[j]) >= etaHigh)) {
          continue;
        }
        allEventsHist.Fill(reco_->pt[j]);
        int l1mu2(-1);
        if (!findBestRecoMatch(l1_, reco_, l1mu2, pTcut, j, 0.5)) {
          continue;
        }
        effHist.Fill(reco_->pt[j]);
      }
    }
  }
  effErrors.Divide(&effHist, &allEventsHist);
  effHist.Divide(&allEventsHist);
}

void getSingleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pTcut, const double etaLow,
                             const double etaHigh, TH1D& effHist,
                             TGraphAsymmErrors& effErrors) {
  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* l1_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisGeneratorDataFormat* gen_ =
      new L1Analysis::L1AnalysisGeneratorDataFormat();
  l1Chain->SetBranchAddress("L1Upgrade", &l1_);
  genChain->SetBranchAddress("Generator", &gen_);

  TH1D allEventsHist =
      TH1D("allEventsHist", "", nMuBins, muLo - 0.1, muHi + 0.1);

  effHist.Sumw2();
  allEventsHist.Sumw2();

  std::cout << "Running over " << nentries << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if ((jentry % 1000) == 0) {
      std::cout << "Done " << jentry << " events..." << std::endl;
    }

    l1Chain->GetEntry(jentry);
    genChain->GetEntry(jentry);

    int genMu1 = -1;
    int genMu2 = -1;
    // Check for correct number of gen muons
    if (!findGenMuon(gen_, genMu1)) {
      continue;
    }

    // Check if in eta region.
    if ((abs(gen_->partEta[genMu1]) < etaLow) ||
        (abs(gen_->partEta[genMu1]) >= etaHigh)) {
      continue;
    }

    allEventsHist.Fill(gen_->partPt[genMu1]);

    int mu = -1;
    if (!findBestGenMatch(l1_, gen_, mu, pTcut, genMu1, 0.5)) {
      continue;
    }

    effHist.Fill(gen_->partPt[genMu1]);
  }
  effErrors.Divide(&effHist, &allEventsHist);
  effHist.Divide(&allEventsHist);
}

void getDoubleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pT1cut, const int pT2cut,
                             const double etaLow, const double etaHigh,
                             TH1D& effHist, TGraphAsymmErrors& effErrors) {
  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* l1_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisGeneratorDataFormat* gen_ =
      new L1Analysis::L1AnalysisGeneratorDataFormat();
  l1Chain->SetBranchAddress("L1Upgrade", &l1_);
  genChain->SetBranchAddress("Generator", &gen_);

  TH1D allEventsHist =
      TH1D("allEventsHist", "", nMuBins, muLo - 0.1, muHi + 0.1);

  effHist.Sumw2();
  allEventsHist.Sumw2();

  std::cout << "Running over " << nentries << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if ((jentry % 1000) == 0) {
      std::cout << "Done " << jentry << " events..." << std::endl;
    }

    l1Chain->GetEntry(jentry);
    genChain->GetEntry(jentry);

    int genMu1 = -1;
    int genMu2 = -1;
    // Check for correct number of gen muons
    if (!findGenMuon(gen_, 2, genMu1, genMu2)) {
      continue;
    }

    // Check if in eta region.
    if ((abs(gen_->partEta[genMu1]) < etaLow) ||
        (abs(gen_->partEta[genMu1]) >= etaHigh) ||
        (abs(gen_->partEta[genMu2]) < etaLow) ||
        (abs(gen_->partEta[genMu2]) >= etaHigh)) {
      continue;
    }

    allEventsHist.Fill(gen_->partPt[genMu1]);
    allEventsHist.Fill(gen_->partPt[genMu2]);

    // Match L1 Muon to gen muon.
    int mu1 = -1;
    int mu2 = -1;
    if (!findBestGenMatches(l1_, gen_, mu1, mu2, pT1cut, pT2cut, genMu1, genMu2,
                            0.5)) {
      continue;
    }

    effHist.Fill(gen_->partPt[genMu1]);
    effHist.Fill(gen_->partPt[genMu2]);
  }
  effErrors.Divide(&effHist, &allEventsHist);
  effHist.Divide(&allEventsHist);
}

bool findGenMuon(L1Analysis::L1AnalysisGeneratorDataFormat* gen_, int& mu1) {
  int dummy;
  return findGenMuon(gen_, 1, mu1, dummy);
}

bool findGenMuon(L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                 const int nMus, int& mu1, int& mu2) {
  int nGenMus = 0;
  for (int i = 0; i < gen_->nPart; ++i) {
    if (abs(gen_->partId[i]) == 13) {
      ++nGenMus;
      mu2 = mu1;
      mu1 = i;
    }
  }
  if (nGenMus != nMus) {
    std::cout << "[ERROR] Found " << nGenMus << " generated muons instead of "
              << nMus << ". Are you sure you're running on the correct dataset?"
              << std::endl;

    return false;
  }
  return true;
}

double recoDist(L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_, const int mu1,
                const int mu2) {
  return sqrt(pow(reco_->eta[mu1] - reco_->eta[mu2], 2) +
              pow(reco_->phi[mu1] - reco_->phi[mu2], 2));
}

double matchL1toGen(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                    L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                    const int l1Mu, const int genMu) {
  return sqrt(pow(upgrade_->muonEta[l1Mu] - gen_->partEta[genMu], 2) +
              pow(upgrade_->muonPhi[l1Mu] - gen_->partPhi[genMu], 2));
}

bool findBestGenMatch(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                      L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                      int& l1mu, const int pTcut, const int genMu,
                      const double dRcut) {
  std::map<int, double> muCands;
  for (int i = 0; i < upgrade_->nMuons; ++i) {
    if (upgrade_->muonBx[i] != 0) {
      continue;
    }
    double dR = matchL1toGen(upgrade_, gen_, i, genMu);
    if (upgrade_->muonEt[i] > pTcut && dR <= dRcut) {
      muCands[i] = dR;
    }
  }
  int bestMu(-1);
  double bestDr(999);
  for (std::map<int, double>::iterator it = muCands.begin();
       it != muCands.end(); ++it) {
    if (it->second < bestDr) {
      bestDr = it->second;
      bestMu = it->first;
    }
  }
  l1mu = bestMu;
  if (bestDr > dRcut) {
    return false;
  } else {
    return true;
  }
}

bool findBestGenMatches(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                        L1Analysis::L1AnalysisGeneratorDataFormat* gen_,
                        int& l1mu1, int& l1mu2, const int pT1cut,
                        const int pT2cut, const int genMu1, const int genMu2,
                        const double dRcut) {
  std::map<int, double> mu1cands, mu2cands;

  // Find all L1 muons that can be matched to the two gen muons and pass the pT
  // cut.
  for (int i = 0; i < upgrade_->nMuons; ++i) {
    if (upgrade_->muonBx[i] != 0) {
      continue;
    }
    double dR1 = matchL1toGen(upgrade_, gen_, i, genMu1);
    double dR2 = matchL1toGen(upgrade_, gen_, i, genMu2);
    if (upgrade_->muonEt[i] > pT1cut && dR1 <= dRcut) {
      mu1cands[i] = dR1;
    }
    if (upgrade_->muonEt[i] > pT2cut && dR2 <= dRcut) {
      mu2cands[i] = dR2;
    }
  }

  // If none could be found we return with false.
  if (mu1cands.size() == 0 || mu2cands.size() == 0) {
    return false;
  }

  // Find the closest L1 candidate to the first gen muon.
  int bestMu1(-1);
  double bestDr1(999);
  for (std::map<int, double>::iterator it = mu1cands.begin();
       it != mu1cands.end(); ++it) {
    if (it->second < bestDr1) {
      bestDr1 = it->second;
      bestMu1 = it->first;
    }
  }

  // Find the closest L1 candidate (that wasn't matched to the first gen muon)
  // to the second gen muon.
  int bestMu2(-1);
  double bestDr2(999);
  for (std::map<int, double>::iterator it = mu2cands.begin();
       it != mu2cands.end(); ++it) {
    if (it->second < bestDr2 && it->second != bestMu1) {
      bestDr2 = it->second;
      bestMu2 = it->first;
    }
  }
  double dRglobal(bestDr1 + bestDr2);

  // Find the closest L1 candidate to the second gen muon (including the ones
  // matched before)
  int bestMu22(-1);
  double bestDr22(999);
  for (std::map<int, double>::iterator it = mu2cands.begin();
       it != mu2cands.end(); ++it) {
    if (it->second < bestDr22) {
      bestDr22 = it->second;
      bestMu22 = it->first;
    }
  }

  // Find the closest L1 candidate to the first gen muon (_excluding_ the one
  // matched to the second gen muon)
  int bestMu21(-1);
  double bestDr21(999);
  for (std::map<int, double>::iterator it = mu1cands.begin();
       it != mu1cands.end(); ++it) {
    if (it->second < bestDr21 && it->second != bestMu22) {
      bestDr21 = it->second;
      bestMu21 = it->first;
    }
  }
  double dRglobal2(bestDr21 + bestDr22);

  // Compare the sum of dRs for the two scenarios and match accordingly.
  if ((bestMu1 != bestMu2) && (dRglobal < dRglobal2)) {
    l1mu1 = bestMu1;
    l1mu2 = bestMu2;
    return true;
  } else if (bestMu21 != bestMu22) {
    l1mu1 = bestMu21;
    l1mu2 = bestMu22;
    return true;
  } else if (bestMu1 != bestMu2) {
    l1mu1 = bestMu1;
    l1mu2 = bestMu2;
    return true;
  }
  return false;
}

double matchL1toReco(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                     L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                     const int l1Mu, const int recoMu) {
  return sqrt(pow(upgrade_->muonEta[l1Mu] - reco_->eta[recoMu], 2) +
              pow(upgrade_->muonPhi[l1Mu] - reco_->phi[recoMu], 2));
}
bool findBestRecoMatch(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade_,
                       L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_,
                       int& l1mu, const int pTcut, const int recoMu,
                       const double dRcut) {
  std::map<int, double> muCands;
  for (int i = 0; i < upgrade_->nMuons; ++i) {
    if (upgrade_->muonBx[i] != 0) {
      continue;
    }
    double dR = matchL1toReco(upgrade_, reco_, i, recoMu);
    if (upgrade_->muonEt[i] > pTcut && dR <= dRcut) {
      muCands[i] = dR;
    }
  }
  int bestMu(-1);
  double bestDr(999);
  for (std::map<int, double>::iterator it = muCands.begin();
       it != muCands.end(); ++it) {
    if (it->second < bestDr) {
      bestDr = it->second;
      bestMu = it->first;
    }
  }
  l1mu = bestMu;
  if (bestDr > dRcut) {
    return false;
  } else {
    return true;
  }
}

void prepareHistograms(TLegend& l, std::vector<TH1D>& hists,
                       const std::vector<int> colours,
                       const std::vector<int> markers,
                       const std::vector<std::string>& histnames,
                       const std::string& type, const std::string& name) {
  std::vector<int>::const_iterator colour = colours.begin();
  std::vector<int>::const_iterator marker = markers.begin();
  std::vector<std::string>::const_iterator histname = histnames.begin();
  for (std::vector<TH1D>::iterator hist = hists.begin(); hist != hists.end();
       ++hist, ++colour, ++marker, ++histname) {
    // TODO: Set max and min for histograms.
    // TODO: Draw horizontal line at '1'?
    hist->SetMinimum(0);
    hist->SetMaximum(1.4);
    hist->SetLineColor(*colour);
    hist->GetXaxis()->SetTitle("p_{T}^{reco} (GeV/c)");
    // TODO: This needs to be configurable (e.g. for the rho factor)
    hist->GetYaxis()->SetTitle(type.c_str());
    hist->SetMarkerStyle(*marker);
    hist->SetMarkerColor(*colour);
    l.AddEntry(&(*hist), histname->c_str(), "lp");
  }

  l.SetFillColor(0);
  l.SetFillStyle(0);
  l.SetBorderSize(0);
}

void DrawHistograms(std::vector<TH1D>& hists, const std::vector<int> colours,
                    const std::vector<int> markers,
                    const std::vector<std::string>& histnames,
                    const std::string& type, const std::string& name) {
  // TODO: cmstdr style!!

  TCanvas c;
  TLegend l(0.4, 0.23, 0.6, 0.38);

  prepareHistograms(l, hists, colours, markers, histnames, type, name);

  for (std::vector<TH1D>::iterator hist = hists.begin(); hist != hists.end();
       ++hist) {
    hist->Draw("same,E1HIST");
  }

  l.Draw("same");

  std::ostringstream oss1, oss2, oss3;
  oss1 << name << ".pdf";
  oss2 << name << ".png";
  oss3 << name << ".root";
  c.SaveAs(oss1.str().c_str());
  c.SaveAs(oss2.str().c_str());
  c.SaveAs(oss3.str().c_str());
}

void DrawHistograms(std::vector<TH1D>& hists, const std::vector<int> colours,
                    const std::vector<int> markers,
                    const std::vector<std::string>& histnames,
                    const std::string& type,
                    std::vector<TGraphAsymmErrors>& errs,
                    const std::string& name) {
  // TODO: cmstdr style!!

  TCanvas c;
  TLegend l(0.4, 0.23, 0.6, 0.38);

  prepareHistograms(l, hists, colours, markers, histnames, type, name);

  std::vector<TGraphAsymmErrors>::iterator err = errs.begin();
  std::vector<int>::const_iterator colour = colours.begin();
  for (std::vector<TH1D>::iterator hist = hists.begin(); hist != hists.end();
       ++hist, ++err, ++colour) {
    hist->Draw("same,HIST");
    err->SetMarkerColor(*colour);
    err->SetLineColor(*colour);
    err->Draw("p,same");
  }

  l.Draw("same");

  std::ostringstream oss1, oss2, oss3;
  oss1 << name << ".pdf";
  oss2 << name << ".png";
  oss3 << name << ".root";
  c.SaveAs(oss1.str().c_str());
  c.SaveAs(oss2.str().c_str());
  c.SaveAs(oss3.str().c_str());
}
