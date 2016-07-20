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

#include "CMS_lumi.C"
#include "tdrstyle.C"

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
// pT bins
const int nMuBins = 35;
const float muLo = 0;
const float muHi = 70;
// eta bins
const int nEtaBins = 50;
const float etaLo = -2.5;
const float etaHi = 2.5;
// Tag and probe constants
const int tagPt = 27;

bool readFList(std::string fname, std::vector<std::string>& listNtuples);
int setupTChain(const std::vector<std::string> listNtuples, TChain* unpackChain,
                TChain* recoChain);
void getSingleMuDataEfficiency(int nentries, TChain* l1Chain, TChain* recoChain,
                               const int pTcut,
                               const std::vector<double>& etaLows,
                               const std::vector<double>& etaHighs,
                               std::vector<TH1D>& effHists,
                               std::vector<TGraphAsymmErrors>& effErrors,
                               bool retrieve_hists, std::string histFolder,
                               std::vector<std::string> histNames);
void getSingleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pTcut,
                             const std::vector<double>& etaLows,
                             const std::vector<double>& etaHighs,
                             std::vector<TH1D>& effHists,
                             std::vector<TGraphAsymmErrors>& effErrors,
                             bool retrieve_hists, std::string histFolder,
                             std::vector<std::string> histNames);
void getDoubleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pT1cut, const int pT2cut,
                             const std::vector<double>& etaLows,
                             const std::vector<double>& etaHighs,
                             std::vector<TH1D>& effHists,
                             std::vector<TGraphAsymmErrors>& effErrors,
                             bool retrieve_hists, std::string histFolder,
                             std::vector<std::string> histNames);
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
                    std::string diMuMcFile, std::string folder,
                    bool retrieve_hists = false, int mu1cut = 2,
                    int mu2cut = 2) {
  setTDRStyle();
  writeExtraText = true;  // Add "Preliminary"
  // lumi_13TeV  = "4.9 fb^{-1}";
  lumi_sqrtS = "13 TeV";  // used with iPeriod = 0, e.g. for simulation-only
                          // plots (default is an empty string)

  int iPeriod = 0;  // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses
                    // lumi_sqrtS)

  std::string plotFolder = "plots/" + folder + "/";
  std::cout << "Creating directory: " << plotFolder << std::endl;
  mkdir("plots/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  const int dir_err =
      mkdir(plotFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (dir_err == -1) {
    std::cout << "Error creating directory or directory exists already."
              << std::endl;
    return;
  }

  std::string histFolder = "hists/" + folder + "/";
  mkdir("hists/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir(histFolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // OpenWithoutInit
  TChain* l1SingleDataChain = new TChain(unpackTreepath.c_str());
  TChain* recoSingleDataChain = new TChain(recoTreepath.c_str());
  TChain* l1SingleMcChain = new TChain(unpackTreepath.c_str());
  TChain* genSingleMcChain = new TChain(genTreepath.c_str());
  TChain* l1DoubleMcChain = new TChain(unpackTreepath.c_str());
  TChain* genDoubleMcChain = new TChain(genTreepath.c_str());

  int singleDataEntries(0);
  int singleMcEntries(0);
  int doubleMcEntries(0);

  if (!retrieve_hists) {
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

    singleDataEntries = setupTChain(listSingleDataNtuples, l1SingleDataChain,
                                    recoSingleDataChain);
    singleMcEntries =
        setupTChain(listSingleMcNtuples, l1SingleMcChain, genSingleMcChain);
    doubleMcEntries =
        setupTChain(listDoubleMcNtuples, l1DoubleMcChain, genDoubleMcChain);
  }

  // make histos
  std::vector<std::string> singleMuMcEffNames;
  singleMuMcEffNames.push_back("totalSingleMuMcEfficiency");
  singleMuMcEffNames.push_back("bmtfSingleMuMcEfficiency");
  singleMuMcEffNames.push_back("bomtfSingleMuMcEfficiency");
  singleMuMcEffNames.push_back("omtfSingleMuMcEfficiency");
  singleMuMcEffNames.push_back("oemtfSingleMuMcEfficiency");
  singleMuMcEffNames.push_back("emtfSingleMuMcEfficiency");
  std::vector<std::string> singleMuDataEffNames;
  singleMuDataEffNames.push_back("totalSingleMuDataEfficiency");
  singleMuDataEffNames.push_back("bmtfSingleMuDataEfficiency");
  singleMuDataEffNames.push_back("bomtfSingleMuDataEfficiency");
  singleMuDataEffNames.push_back("omtfSingleMuDataEfficiency");
  singleMuDataEffNames.push_back("eomtfSingleMuDataEfficiency");
  singleMuDataEffNames.push_back("emtfSingleMuDataEfficiency");
  std::vector<std::string> doubleMuMcEffNames;
  doubleMuMcEffNames.push_back("totalDoubleMuMcEfficiency");
  doubleMuMcEffNames.push_back("bmtfDoubleMuMcEfficiency");
  doubleMuMcEffNames.push_back("bomtfDoubleMuMcEfficiency");
  doubleMuMcEffNames.push_back("omtfDoubleMuMcEfficiency");
  doubleMuMcEffNames.push_back("eomtfDoubleMuMcEfficiency");
  doubleMuMcEffNames.push_back("emtfDoubleMuMcEfficiency");

  std::vector<TH1D> singleMuMcEffs;
  std::vector<TH1D> singleMuDataEffs;
  std::vector<TH1D> doubleMuMcEffs;
  for (int i = 0; i < singleMuMcEffNames.size(); ++i) {
    TH1D singleMuMcEffHist(singleMuMcEffNames.at(i).c_str(), "", nMuBins,
                           muLo - 0.1, muHi + 0.1);
    singleMuMcEffs.push_back(singleMuMcEffHist);
    TH1D singleMuDataEffHist(singleMuDataEffNames.at(i).c_str(), "", nMuBins,
                             muLo - 0.1, muHi + 0.1);
    singleMuDataEffs.push_back(singleMuDataEffHist);
    TH1D doubleMuMcEffHist(doubleMuMcEffNames.at(i).c_str(), "", nMuBins,
                           muLo - 0.1, muHi + 0.1);
    doubleMuMcEffs.push_back(doubleMuMcEffHist);
  }

  // Full coverage
  const double totalLow = 0;
  const double totalHigh = 2.5;
  TH1D totalDoubleMuDataEfficiency("totalDoubleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D totalNaiveDoubleMuDataEfficiency("totalNaiveDoubleMuDataEfficiency", "",
                                        nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D totalNaiveDoubleMuMcEfficiency("totalNaiveDoubleMuMcEfficiency", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  // BMTF
  const double bmtfLow = 0;
  const double bmtfHigh = 0.7;
  TH1D bmtfDoubleMuDataEfficiency("bmtfDoubleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D bmtfNaiveDoubleMuDataEfficiency("bmtfNaiveDoubleMuDataEfficiency", "",
                                       nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D bmtfNaiveDoubleMuMcEfficiency("bmtfNaiveDoubleMuMcEfficiency", "",
                                     nMuBins, muLo - 0.1, muHi + 0.1);
  // BOMTF
  const double bomtfLow = 0.7;
  const double bomtfHigh = 0.9;
  TH1D bomtfDoubleMuDataEfficiency("bomtfDoubleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D bomtfNaiveDoubleMuDataEfficiency("bomtfNaiveDoubleMuDataEfficiency", "",
                                        nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D bomtfNaiveDoubleMuMcEfficiency("bomtfNaiveDoubleMuMcEfficiency", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  // OMTF
  const double omtfLow = 0.9;
  const double omtfHigh = 1.15;
  TH1D omtfDoubleMuDataEfficiency("omtfDoubleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D omtfNaiveDoubleMuDataEfficiency("omtfNaiveDoubleMuDataEfficiency", "",
                                       nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D omtfNaiveDoubleMuMcEfficiency("omtfNaiveDoubleMuMcEfficiency", "",
                                     nMuBins, muLo - 0.1, muHi + 0.1);
  // EOMTF
  const double eomtfLow = 1.15;
  const double eomtfHigh = 1.35;
  TH1D eomtfDoubleMuDataEfficiency("eomtfDoubleMuDataEfficiency", "", nMuBins,
                                   muLo - 0.1, muHi + 0.1);
  TH1D eomtfNaiveDoubleMuDataEfficiency("eomtfNaiveDoubleMuDataEfficiency", "",
                                        nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D eomtfNaiveDoubleMuMcEfficiency("eomtfNaiveDoubleMuMcEfficiency", "",
                                      nMuBins, muLo - 0.1, muHi + 0.1);
  // EMTF
  const double emtfLow = 1.35;
  const double emtfHigh = 2.5;
  TH1D emtfDoubleMuDataEfficiency("emtfDoubleMuDataEfficiency", "", nMuBins,
                                  muLo - 0.1, muHi + 0.1);
  TH1D emtfNaiveDoubleMuDataEfficiency("emtfNaiveDoubleMuDataEfficiency", "",
                                       nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D emtfNaiveDoubleMuMcEfficiency("emtfNaiveDoubleMuMcEfficiency", "",
                                     nMuBins, muLo - 0.1, muHi + 0.1);

  // For correct error bars
  TGraphAsymmErrors totalSingleMuDataErrors = TGraphAsymmErrors();
  TGraphAsymmErrors totalSingleMuMcErrors = TGraphAsymmErrors();
  TGraphAsymmErrors totalDoubleMuMcErrors = TGraphAsymmErrors();
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

  std::vector<double> etaLows;
  etaLows.push_back(bmtfLow);
  etaLows.push_back(bomtfLow);
  etaLows.push_back(omtfLow);
  etaLows.push_back(eomtfLow);
  etaLows.push_back(emtfLow);
  etaLows.push_back(totalLow);
  std::vector<double> etaHighs;
  etaHighs.push_back(bmtfHigh);
  etaHighs.push_back(bomtfHigh);
  etaHighs.push_back(omtfHigh);
  etaHighs.push_back(eomtfHigh);
  etaHighs.push_back(emtfHigh);
  etaHighs.push_back(totalHigh);

  std::vector<TGraphAsymmErrors> singleMuMcErrs;
  singleMuMcErrs.push_back(bmtfSingleMuMcErrors);
  singleMuMcErrs.push_back(bomtfSingleMuMcErrors);
  singleMuMcErrs.push_back(omtfSingleMuMcErrors);
  singleMuMcErrs.push_back(eomtfSingleMuMcErrors);
  singleMuMcErrs.push_back(emtfSingleMuMcErrors);
  singleMuMcErrs.push_back(totalSingleMuMcErrors);
  std::vector<TGraphAsymmErrors> singleMuDataErrs;
  singleMuDataErrs.push_back(bmtfSingleMuDataErrors);
  singleMuDataErrs.push_back(bomtfSingleMuDataErrors);
  singleMuDataErrs.push_back(omtfSingleMuDataErrors);
  singleMuDataErrs.push_back(eomtfSingleMuDataErrors);
  singleMuDataErrs.push_back(emtfSingleMuDataErrors);
  singleMuDataErrs.push_back(totalSingleMuDataErrors);
  std::vector<TGraphAsymmErrors> doubleMuMcErrs;
  doubleMuMcErrs.push_back(bmtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(bomtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(omtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(eomtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(emtfDoubleMuMcErrors);
  doubleMuMcErrs.push_back(totalDoubleMuMcErrors);

  getSingleMuDataEfficiency(singleDataEntries, l1SingleDataChain,
                            recoSingleDataChain, mu1cut, etaLows, etaHighs,
                            singleMuDataEffs, singleMuDataErrs, retrieve_hists,
                            histFolder, singleMuDataEffNames);

  getSingleMuMcEfficiency(singleMcEntries, l1SingleMcChain, genSingleMcChain,
                          mu1cut, etaLows, etaHighs, singleMuMcEffs,
                          singleMuMcErrs, retrieve_hists, histFolder,
                          singleMuMcEffNames);

  getDoubleMuMcEfficiency(doubleMcEntries, l1DoubleMcChain, genDoubleMcChain,
                          mu1cut, mu2cut, etaLows, etaHighs, doubleMuMcEffs,
                          doubleMuMcErrs, retrieve_hists, histFolder,
                          doubleMuMcEffNames);

  std::vector<TH1D> naiveDoubleMuMcEffs;
  std::vector<TH1D> naiveDoubleMuDataEffs;
  std::vector<TH1D> rhoFactors;
  std::vector<TH1D> doubleMuDataEffs;

  TH1D naiveDoubleMuDataEff("naiveDoubleMuDataEff", "", nMuBins, muLo - 0.1,
                            muHi + 0.1);
  TH1D naiveDoubleMuMcEff("naiveDoubleMuMcEff", "", nMuBins, muLo - 0.1,
                          muHi + 0.1);
  TH1D rhoFactor("rhoFactor", "", nMuBins, muLo - 0.1, muHi + 0.1);
  TH1D doubleMuDataEff("doubleMuDataEff", "", nMuBins, muLo - 0.1, muHi + 0.1);
  for (int i = 0; i < singleMuMcEffs.size(); ++i) {
    // Squaring single mu efficiencies to get "naive" double mu efficiencies.
    naiveDoubleMuDataEff.Multiply(&(singleMuDataEffs.at(i)),
                                  &(singleMuDataEffs.at(i)), 1, 1, "B");
    naiveDoubleMuDataEffs.push_back(naiveDoubleMuDataEff);
    naiveDoubleMuMcEff.Multiply(&(singleMuMcEffs.at(i)),
                                &(singleMuMcEffs.at(i)), 1, 1, "B");
    naiveDoubleMuMcEffs.push_back(naiveDoubleMuMcEff);

    // Calculating rho factor
    rhoFactor.Divide(
        &(doubleMuMcEffs.at(i)),
        &naiveDoubleMuMcEff);  // Two different datasets, no binomial errors.
    rhoFactors.push_back(rhoFactor);

    // Calculating double mu eff from data
    doubleMuDataEff.Multiply(&naiveDoubleMuDataEff, &rhoFactor);
    doubleMuDataEffs.push_back(doubleMuDataEff);
  }

  std::vector<std::string> regionNames;
  std::ostringstream oss1, oss2, oss3, oss4, oss5, ossTotal;
  oss1 << bmtfLow << " #leq |#eta| < " << bmtfHigh;
  oss2 << bomtfLow << " #leq |#eta| < " << bomtfHigh;
  oss3 << omtfLow << " #leq |#eta| < " << omtfHigh;
  oss4 << eomtfLow << " #leq |#eta| < " << eomtfHigh;
  oss5 << emtfLow << " #leq |#eta| < " << emtfHigh;
  ossTotal << totalLow << " #leq |#eta| < "
           << totalHigh;  // Will be plotted separately.
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

  std::vector<TH1D> totalSingleMuMcEff(singleMuMcEffs.end() - 1,
                                       singleMuMcEffs.end());
  std::vector<TGraphAsymmErrors> totalSingleMuMcErr(singleMuMcErrs.end() - 1,
                                                    singleMuMcErrs.end());
  std::vector<TH1D> totalSingleMuDataEff(singleMuDataEffs.end() - 1,
                                         singleMuDataEffs.end());
  std::vector<TGraphAsymmErrors> totalSingleMuDataErr(
      singleMuDataErrs.end() - 1, singleMuDataErrs.end());
  std::vector<TH1D> totalDoubleMuMcEff(doubleMuMcEffs.end() - 1,
                                       doubleMuMcEffs.end());
  std::vector<TGraphAsymmErrors> totalDoubleMuMcErr(doubleMuMcErrs.end() - 1,
                                                    doubleMuMcErrs.end());
  std::vector<TH1D> totalNaiveDoubleMuDataEff(naiveDoubleMuDataEffs.end() - 1,
                                              naiveDoubleMuDataEffs.end());
  std::vector<TH1D> totalNaiveDoubleMuMcEff(naiveDoubleMuMcEffs.end() - 1,
                                            naiveDoubleMuMcEffs.end());
  std::vector<TH1D> totalRhoFactor(rhoFactors.end() - 1, rhoFactors.end());
  std::vector<TH1D> totalDoubleMuDataEff(doubleMuDataEffs.end() - 1,
                                         doubleMuDataEffs.end());

  // Removing plot for entire coverage.
  singleMuMcEffs.pop_back();
  singleMuMcErrs.pop_back();
  singleMuDataEffs.pop_back();
  singleMuDataErrs.pop_back();
  doubleMuMcEffs.pop_back();
  doubleMuMcErrs.pop_back();
  naiveDoubleMuDataEffs.pop_back();
  naiveDoubleMuMcEffs.pop_back();
  rhoFactors.pop_back();
  doubleMuDataEffs.pop_back();

  DrawHistograms(singleMuMcEffs, colours, markers, regionNames,
                 "L1T Efficiency", singleMuMcErrs,
                 plotFolder + "singleMuonEfficiencies_MC");

  DrawHistograms(singleMuDataEffs, colours, markers, regionNames,
                 "L1T Efficiency", singleMuDataErrs,
                 plotFolder + "singleMuonEfficiencies_Data");

  DrawHistograms(doubleMuMcEffs, colours, markers, regionNames,
                 "L1T Efficiency", doubleMuMcErrs,
                 plotFolder + "doubleMuonEfficiencies_MC");

  DrawHistograms(naiveDoubleMuMcEffs, colours, markers, regionNames,
                 "L1T Efficiency",
                 plotFolder + "naiveDoubleMuonEfficiencies_MC");

  DrawHistograms(naiveDoubleMuDataEffs, colours, markers, regionNames,
                 "L1T Efficiency",
                 plotFolder + "naiveDoubleMuonEfficiencies_Data");

  DrawHistograms(rhoFactors, colours, markers, regionNames, "#rho",
                 plotFolder + "rhoFactors");

  DrawHistograms(doubleMuDataEffs, colours, markers, regionNames,
                 "L1T Efficiency", plotFolder + "doubleMuonEfficiencies_Data");

  // Draw individual plots for histograms by being lazy..

  std::vector<TH1D> bmtfDoubleMuDataEff;
  bmtfDoubleMuDataEff.push_back(doubleMuDataEffs.at(0));
  std::vector<TH1D> bomtfDoubleMuDataEff;
  bomtfDoubleMuDataEff.push_back(doubleMuDataEffs.at(1));
  std::vector<TH1D> omtfDoubleMuDataEff;
  omtfDoubleMuDataEff.push_back(doubleMuDataEffs.at(2));
  std::vector<TH1D> eomtfDoubleMuDataEff;
  eomtfDoubleMuDataEff.push_back(doubleMuDataEffs.at(3));
  std::vector<TH1D> emtfDoubleMuDataEff;
  emtfDoubleMuDataEff.push_back(doubleMuDataEffs.at(4));

  std::vector<std::string> totalName;
  std::vector<std::string> bmtfName;
  std::vector<std::string> bomtfName;
  std::vector<std::string> omtfName;
  std::vector<std::string> eomtfName;
  std::vector<std::string> emtfName;
  totalName.push_back(ossTotal.str());
  bmtfName.push_back(oss1.str());
  bomtfName.push_back(oss2.str());
  omtfName.push_back(oss3.str());
  eomtfName.push_back(oss4.str());
  emtfName.push_back(oss5.str());

  std::vector<int> totalColour;
  totalColour.push_back(kRed);
  std::vector<int> bmtfColour;
  bmtfColour.push_back(41);
  std::vector<int> bomtfColour;
  bomtfColour.push_back(9);
  std::vector<int> omtfColour;
  omtfColour.push_back(36);
  std::vector<int> eomtfColour;
  eomtfColour.push_back(3);
  std::vector<int> emtfColour;
  emtfColour.push_back(48);

  std::vector<int> totalMarker;
  totalMarker.push_back(2);
  std::vector<int> bmtfMarker;
  bmtfMarker.push_back(24);
  std::vector<int> bomtfMarker;
  bomtfMarker.push_back(25);
  std::vector<int> omtfMarker;
  omtfMarker.push_back(26);
  std::vector<int> eomtfMarker;
  eomtfMarker.push_back(27);
  std::vector<int> emtfMarker;
  emtfMarker.push_back(32);

  DrawHistograms(totalSingleMuMcEff, totalColour, markers, totalName,
                 "L1T Efficiency", totalSingleMuMcErr,
                 plotFolder + "totalSingleMuonEfficiencies_MC");
  DrawHistograms(totalSingleMuDataEff, totalColour, markers, totalName,
                 "L1T Efficiency", totalSingleMuDataErr,
                 plotFolder + "totalSingleMuonEfficiencies_Data");
  DrawHistograms(totalDoubleMuMcEff, totalColour, markers, totalName,
                 "L1T Efficiency", totalDoubleMuMcErr,
                 plotFolder + "totalDoubleMuonEfficiencies_MC");
  DrawHistograms(totalNaiveDoubleMuDataEff, totalColour, markers, totalName,
                 "L1T Efficiency",
                 plotFolder + "totalNaiveDoubleMuonEfficiencies_Data");
  DrawHistograms(totalNaiveDoubleMuMcEff, totalColour, markers, totalName,
                 "L1T Efficiency",
                 plotFolder + "totalNaiveDoubleMuonEfficiencies_MC");
  DrawHistograms(totalRhoFactor, totalColour, markers, totalName,
                 "L1T Efficiency", plotFolder + "totalRhoFactor");
  DrawHistograms(totalDoubleMuDataEff, totalColour, markers, totalName,
                 "L1T Efficiency",
                 plotFolder + "totalDoubleMuonEfficiencies_Data");

  DrawHistograms(bmtfDoubleMuDataEff, bmtfColour, markers, bmtfName,
                 "L1T Efficiency",
                 plotFolder + "bmtfDoubleMuonEfficiencies_Data");
  DrawHistograms(bomtfDoubleMuDataEff, bomtfColour, markers, bomtfName,
                 "L1T Efficiency",
                 plotFolder + "bomtfDoubleMuonEfficiencies_Data");
  DrawHistograms(omtfDoubleMuDataEff, omtfColour, markers, omtfName,
                 "L1T Efficiency",
                 plotFolder + "omtfDoubleMuonEfficiencies_Data");
  DrawHistograms(eomtfDoubleMuDataEff, eomtfColour, markers, eomtfName,
                 "L1T Efficiency",
                 plotFolder + "eomtfDoubleMuonEfficiencies_Data");
  DrawHistograms(emtfDoubleMuDataEff, emtfColour, markers, emtfName,
                 "L1T Efficiency",
                 plotFolder + "emtfDoubleMuonEfficiencies_Data");
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
                               const int pTcut,
                               const std::vector<double>& etaLows,
                               const std::vector<double>& etaHighs,
                               std::vector<TH1D>& effHists,
                               std::vector<TGraphAsymmErrors>& effErrors,
                               bool retrieve_hists, std::string histFolder,
                               std::vector<std::string> histNames) {
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

  std::vector<std::string> allEventsHistNames;
  std::vector<TH1D> allEventsHists;
  for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
    std::ostringstream oss;
    oss << "singleMuDataAllEventsHist" << nRegion;
    allEventsHistNames.push_back(oss.str());
    TH1D allEventsHist =
        TH1D(oss.str().c_str(), "", nMuBins, muLo - 0.1, muHi + 0.1);
    allEventsHist.Sumw2();
    allEventsHists.push_back(allEventsHist);
    effHists.at(nRegion).Sumw2();
  }

  std::cout << "Running over " << nentries << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    // If we're getting the data from stored histograms we exit.
    if (retrieve_hists) {
      break;
    }

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
      if (reco_->hlt_isomu[i] != 1) {
        continue;
      }
      if (reco_->iso[i] >= 0.1) {
        continue;
      }
      if (reco_->hlt_isoDeltaR[i] >= 0.3) {
        continue;
      }
      int tagMu(i);

      // TODO: Possibly require that outside of eta region?
      for (int j = 0; j < reco_->nMuons; ++j) {
        if (!reco_->isTightMuon[j] || (tagMu == j) ||
            (std::find(probeMus.begin(), probeMus.end(), j) !=
             probeMus.end())) {
          continue;
        }
        // Check that probe and tag are more than dR=0.5 separated.
        if (recoDist(reco_, tagMu, j) <= 0.5) {
          continue;
        }
        // TODO: Possibly check invariant mass.
        probeMus.push_back(j);
        for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
          if ((abs(reco_->eta[j]) >= etaLows.at(nRegion)) &&
              (abs(reco_->eta[j]) < etaHighs.at(nRegion))) {
            allEventsHists.at(nRegion).Fill(reco_->pt[j]);
            int l1mu(-1);
            if (!findBestRecoMatch(l1_, reco_, l1mu, pTcut, j, 0.5)) {
              continue;
            }

            if (l1_->muonQual[l1mu] < 8) {
              continue;
            }

            effHists.at(nRegion).Fill(reco_->pt[j]);
          }
        }
      }
    }
  }

  TFile f;
  for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
    std::ostringstream oss;
    oss << histFolder << "singleMuDataEff_" << nRegion << ".root";
    if (retrieve_hists) {
      std::cout << "Opening file " << oss.str() << std::endl;
      f.Open(oss.str().c_str(), "read");

      std::cout << "Retrieving histogram:" << histNames.at(nRegion) << " and "
                << allEventsHistNames.at(nRegion) << std::endl;
      effHists.at(nRegion) =
          *(static_cast<TH1D*>(f.Get(histNames.at(nRegion).c_str())));
      std::cout << "..." << std::endl;
      allEventsHists.at(nRegion) =
          *(static_cast<TH1D*>(f.Get(allEventsHistNames.at(nRegion).c_str())));
      std::cout << "Done" << std::endl;
    } else {
      f.Open(oss.str().c_str(), "new");

      effHists.at(nRegion).Write();
      allEventsHists.at(nRegion).Write();
    }
    f.Close();

    effErrors.at(nRegion).Divide(&(effHists.at(nRegion)),
                                 &(allEventsHists.at(nRegion)));
    effHists.at(nRegion).Divide(&(allEventsHists.at(nRegion)));
  }
}

void getSingleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pTcut,
                             const std::vector<double>& etaLows,
                             const std::vector<double>& etaHighs,
                             std::vector<TH1D>& effHists,
                             std::vector<TGraphAsymmErrors>& effErrors,
                             bool retrieve_hists, std::string histFolder,
                             std::vector<std::string> histNames) {
  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* l1_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisGeneratorDataFormat* gen_ =
      new L1Analysis::L1AnalysisGeneratorDataFormat();
  l1Chain->SetBranchAddress("L1Upgrade", &l1_);
  genChain->SetBranchAddress("Generator", &gen_);

  std::vector<std::string> allEventsHistNames;
  std::vector<TH1D> allEventsHists;
  for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
    std::ostringstream oss;
    oss << "singleMuMcAllEventsHist" << nRegion;
    allEventsHistNames.push_back(oss.str());
    TH1D allEventsHist =
        TH1D(oss.str().c_str(), "", nMuBins, muLo - 0.1, muHi + 0.1);
    allEventsHist.Sumw2();
    allEventsHists.push_back(allEventsHist);
    effHists.at(nRegion).Sumw2();
  }

  std::cout << "Running over " << nentries << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    // If we're getting the data from stored histograms we exit.
    if (retrieve_hists) {
      break;
    }

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
    for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
      if ((abs(gen_->partEta[genMu1]) >= etaLows.at(nRegion)) &&
          (abs(gen_->partEta[genMu1]) < etaHighs.at(nRegion))) {
        allEventsHists.at(nRegion).Fill(gen_->partPt[genMu1]);

        // TODO: DEBUG
        if (l1_->nMuons < 1) {
          continue;
        }

        int qual(0);
        for (int i = 0; i < l1_->nMuons; ++i) {
          if (l1_->muonQual[i] > qual) {
            qual = l1_->muonQual[i];
          }
        }

        if (qual < 8) {
          continue;
        }

        effHists.at(nRegion).Fill(gen_->partPt[genMu1]);
      }
    }
  }

  TFile f;
  for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
    std::ostringstream oss;
    oss << histFolder << "singleMuMcEff_" << nRegion << ".root";
    if (retrieve_hists) {
      f.Open(oss.str().c_str(), "read");

      effHists.at(nRegion) =
          *(static_cast<TH1D*>(f.Get(histNames.at(nRegion).c_str())));
      allEventsHists.at(nRegion) =
          *(static_cast<TH1D*>(f.Get(allEventsHistNames.at(nRegion).c_str())));
    } else {
      f.Open(oss.str().c_str(), "new");

      effHists.at(nRegion).Write();
      allEventsHists.at(nRegion).Write();
    }
    f.Close();

    effErrors.at(nRegion).Divide(&(effHists.at(nRegion)),
                                 &(allEventsHists.at(nRegion)));
    effHists.at(nRegion).Divide(&(allEventsHists.at(nRegion)));
  }
}

void getDoubleMuMcEfficiency(int nentries, TChain* l1Chain, TChain* genChain,
                             const int pT1cut, const int pT2cut,
                             const std::vector<double>& etaLows,
                             const std::vector<double>& etaHighs,
                             std::vector<TH1D>& effHists,
                             std::vector<TGraphAsymmErrors>& effErrors,
                             bool retrieve_hists, std::string histFolder,
                             std::vector<std::string> histNames) {
  // set branch addresses
  L1Analysis::L1AnalysisL1UpgradeDataFormat* l1_ =
      new L1Analysis::L1AnalysisL1UpgradeDataFormat();
  L1Analysis::L1AnalysisGeneratorDataFormat* gen_ =
      new L1Analysis::L1AnalysisGeneratorDataFormat();
  l1Chain->SetBranchAddress("L1Upgrade", &l1_);
  genChain->SetBranchAddress("Generator", &gen_);

  std::vector<std::string> allEventsHistNames;
  std::vector<TH1D> allEventsHists;
  for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
    std::ostringstream oss;
    oss << "doubleMuAllEventsHist" << nRegion;
    allEventsHistNames.push_back(oss.str());
    TH1D allEventsHist =
        TH1D(oss.str().c_str(), "", nMuBins, muLo - 0.1, muHi + 0.1);
    allEventsHist.Sumw2();
    allEventsHists.push_back(allEventsHist);
    effHists.at(nRegion).Sumw2();
  }

  std::cout << "Running over " << nentries << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    // If we're getting the data from stored histograms we exit.
    if (retrieve_hists) {
      break;
    }

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
    for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
      if ((abs(gen_->partEta[genMu1]) >= etaLows.at(nRegion)) &&
          (abs(gen_->partEta[genMu1]) < etaHighs.at(nRegion)) &&
          (abs(gen_->partEta[genMu2]) >= etaLows.at(nRegion)) &&
          (abs(gen_->partEta[genMu2]) < etaHighs.at(nRegion))) {
        allEventsHists.at(nRegion).Fill(gen_->partPt[genMu1]);
        allEventsHists.at(nRegion).Fill(gen_->partPt[genMu2]);

        // TODO: DEBUG
        if (l1_->nMuons < 2) {
          continue;
        }

        int qual1(0);
        int qual2(0);
        for (int i = 0; i < l1_->nMuons; ++i) {
          if (l1_->muonQual[i] > qual1) {
            qual2 = qual1;
            qual1 = l1_->muonQual[i];
          } else if (l1_->muonQual[i] > qual2) {
            qual2 = l1_->muonQual[i];
          }
        }
        if ((qual1 < 8) && (qual2 < 8)) {
          continue;
        }

        effHists.at(nRegion).Fill(gen_->partPt[genMu1]);
        effHists.at(nRegion).Fill(gen_->partPt[genMu2]);
      }
    }
  }

  TFile f;
  for (int nRegion = 0; nRegion < effHists.size(); ++nRegion) {
    std::ostringstream oss;
    oss << histFolder << "doubleMuMcEff_" << nRegion << ".root";
    if (retrieve_hists) {
      f.Open(oss.str().c_str(), "read");

      effHists.at(nRegion) =
          *(static_cast<TH1D*>(f.Get(histNames.at(nRegion).c_str())));
      allEventsHists.at(nRegion) =
          *(static_cast<TH1D*>(f.Get(allEventsHistNames.at(nRegion).c_str())));
    } else {
      f.Open(oss.str().c_str(), "new");

      effHists.at(nRegion).Write();
      allEventsHists.at(nRegion).Write();
    }
    f.Close();

    effErrors.at(nRegion).Divide(&(effHists.at(nRegion)),
                                 &(allEventsHists.at(nRegion)));
    effHists.at(nRegion).Divide(&(allEventsHists.at(nRegion)));
  }
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

  // Find all L1 muons that can be matched to the two gen muons and pass the
  // pT
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
    hist->SetMarkerSize(0.5);
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
  TLegend l(0.25, 0.2, 0.45, 0.4);

  prepareHistograms(l, hists, colours, markers, histnames, type, name);

  for (std::vector<TH1D>::iterator hist = hists.begin(); hist != hists.end();
       ++hist) {
    hist->Draw("same,E1HIST");
  }

  l.Draw("same");

  CMS_lumi(&c, 0, 11);

  c.Update();
  c.RedrawAxis();
  c.GetFrame()->Draw();

  std::ostringstream oss1, oss2;
  oss1 << name << ".pdf";
  oss2 << name << ".png";
  c.SaveAs(oss1.str().c_str());
  c.SaveAs(oss2.str().c_str());
}

void DrawHistograms(std::vector<TH1D>& hists, const std::vector<int> colours,
                    const std::vector<int> markers,
                    const std::vector<std::string>& histnames,
                    const std::string& type,
                    std::vector<TGraphAsymmErrors>& errs,
                    const std::string& name) {
  // TODO: cmstdr style!!

  TCanvas c;
  TLegend l(0.25, 0.2, 0.45, 0.4);

  prepareHistograms(l, hists, colours, markers, histnames, type, name);

  std::vector<TGraphAsymmErrors>::iterator err = errs.begin();
  std::vector<int>::const_iterator colour = colours.begin();
  std::vector<int>::const_iterator marker = markers.begin();
  for (std::vector<TH1D>::iterator hist = hists.begin(); hist != hists.end();
       ++hist, ++err, ++colour, ++marker) {
    hist->Draw("same,HIST");
    err->SetMarkerStyle(*marker);
    err->SetMarkerSize(0.5);
    err->SetMarkerColor(*colour);
    err->SetLineColor(*colour);
    err->Draw("p,same");
  }

  l.Draw("same");

  CMS_lumi(&c, 0, 11);

  c.Update();
  c.RedrawAxis();
  c.GetFrame()->Draw();

  std::ostringstream oss1, oss2;
  oss1 << name << ".pdf";
  oss2 << name << ".png";
  c.SaveAs(oss1.str().c_str());
  c.SaveAs(oss2.str().c_str());
}
