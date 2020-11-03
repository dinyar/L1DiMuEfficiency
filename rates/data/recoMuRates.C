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

int nMuBinsPt = 18;
float ptBins[19];

bool readFList(std::string fname, std::vector<std::string>& listNtuples,
               std::string recoTreepath);
int setupTChain(const std::vector<std::string> listNtuples, TChain* truthChain);
void getMuonRates(int nCollBunches, int nevents, TChain* recoChain,
                  TH1D& singleMuRateOpenVsPtHist,
                  TH1D& doubleMuRateOpenVsPtHist, bool retrieve_hists,
                  bool update_hists, TString folder, TString identifier);
void calcRates(int nCollBunches, int nevents, TH1D& singleMuRateOpenVsPtHist,
               TH1D& doubleMuRateOpenVsPtHist);
void drawHistograms(TH1D& singleMuHist, TH1D& doubleMuHist, TString filename,
                    TString xAxisLabel, TString descString, TString plotFolder,
                    TString run, TString singleMuDesc, TString doubleMuDesc);

void recoMuRates(const char* file_list, TString folder = "tmp",
                 TString run = "XXX", int nCollBunches = 2028,
                 bool retrieve_hists = false, bool update_hists = false,
                 int nEntries = 0) {
  TString plotFolder = "ratePlots/" + run + "/" + folder + "/";
  mkdir("ratePlots/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("ratePlots/" + run, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir(plotFolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  TString histFolder = "rateHists/" + run + "/" + folder + "/";
  mkdir("rateHists/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  mkdir("rateHists/" + run, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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

  std::string recoTreepath("l1MuonRecoTree/Muon2RecoTree");

  TChain* chainReco = new TChain(recoTreepath.c_str());

  if (!retrieve_hists) {
    std::vector<std::string> listNtuples;

    bool success = readFList(file_list, listNtuples, recoTreepath);

    if (!success) {
      std::cout << "Couldn't read NTuple file list. Exiting.. " << std::endl;
      return;
    }

    nEntries = setupTChain(listNtuples, chainReco);
  }

  for (int i = 0; i < 14; ++i) {
    ptBins[i] = i;
  }
  for (int i = 0; i < 4; ++i) {
    ptBins[14 + i] = 14 + 3 * i;
  }
  ptBins[18] = 30;

  // make histos
  TH1D singleMuRatesOpenVsPt("singleMuRatesOpenVsPt", "", nMuBinsPt, ptBins);
  TH1D doubleMuRatesOpenVsPt("doubleMuRatesOpenVsPt", "", nMuBinsPt, ptBins);
  singleMuRatesOpenVsPt.Sumw2();
  doubleMuRatesOpenVsPt.Sumw2();

  getMuonRates(nCollBunches, nEntries, chainReco, singleMuRatesOpenVsPt,
               doubleMuRatesOpenVsPt, retrieve_hists, update_hists, histFolder,
               "Reco");

  calcRates(nCollBunches, nEntries, singleMuRatesOpenVsPt,
            doubleMuRatesOpenVsPt);

  std::ostringstream singleMuNoCutDesc;
  singleMuNoCutDesc << "Single Muon rate";
  std::ostringstream doubleMuNoCutDesc;
  doubleMuNoCutDesc << "Double Muon rate";
  drawHistograms(singleMuRatesOpenVsPt, doubleMuRatesOpenVsPt,
                 "openRatesVsPtMu", "p_{T} threshold [GeV/c]", "ZeroBias",
                 plotFolder, run, singleMuNoCutDesc.str().c_str(),
                 doubleMuNoCutDesc.str().c_str());
}

bool readFList(std::string fname, std::vector<std::string>& listNtuples,
               std::string recoTreepath) {
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

  TTree* treeReco = (TTree*)rf->Get(recoTreepath.c_str());

  if (!treeReco) {
    std::cout << "No reco tree found.. " << std::endl;
    return false;
  }

  return true;
}

int setupTChain(const std::vector<std::string> listNtuples,
                TChain* truthChain) {
  for (unsigned int i = 0; i < listNtuples.size(); i++) {
    std::cout << " -- Adding " << listNtuples[i] << std::endl;
    truthChain->Add(listNtuples[i].c_str());
  }

  // Init
  std::cout << "Estimate the number of entries... ";
  int nentries = truthChain->GetEntries();
  std::cout << nentries << std::endl;

  return nentries;
}

void getMuonRates(int nCollBunches, int nevents, TChain* recoChain,
                  TH1D& singleMuRateOpenVsPtHist,
                  TH1D& doubleMuRateOpenVsPtHist, bool retrieve_hists,
                  bool update_hists, TString folder, TString identifier) {
  TFile* f;
  if (retrieve_hists || update_hists) {
    std::cout << "Retrieving histograms.. " << std::endl;
    f = TFile::Open(folder + identifier + "Hists.root", "read");
    std::cout << "Opening file: " << folder + identifier + "Hists.root"
              << std::endl;
    singleMuRateOpenVsPtHist =
        *(static_cast<TH1D*>(f->Get("singleMuRatesOpenVsPt")));
    doubleMuRateOpenVsPtHist =
        *(static_cast<TH1D*>(f->Get("doubleMuRatesOpenVsPt")));
    f->Close();
    if (retrieve_hists) {
      return;
    }
  }

  // set branch addresses
  L1Analysis::L1AnalysisRecoMuon2DataFormat* reco_ =
      new L1Analysis::L1AnalysisRecoMuon2DataFormat();
  recoChain->SetBranchAddress("Muon", &reco_);

  for (Long64_t jentry = 0; jentry < nevents; ++jentry) {
    if ((jentry % 100000) == 0)
      std::cout << "Done " << jentry << " events..." << std::endl;

    recoChain->GetEntry(jentry);

    // get Mu rates
    int mu1{-1};
    int mu2{-1};
    double mu1Pt{0};
    double mu2Pt{0};
    for (uint it = 0; it < reco_->nMuons; ++it) {
      if ((!reco_->isLooseMuon[it]) || reco_->phiSt1[it] == -9999) {
        continue;
      }

      if (reco_->pt[it] > mu1Pt) {
        mu2 = mu1;
        mu2Pt = mu1Pt;
        mu1 = it;
        mu1Pt = reco_->pt[it];
      } else if (reco_->pt[it] > mu2Pt) {
        mu2 = it;
        mu2Pt = reco_->pt[it];
      }
    }
    if (mu1 != -1) {
      for (int i = 0; i < nMuBinsPt + 1; ++i) {
        if (mu1Pt >= ptBins[i]) {
          singleMuRateOpenVsPtHist.Fill(ptBins[i]);
        }
      }
    }
    if (mu1 != -1 && mu2 != -1) {
      for (int i = 0; i < nMuBinsPt + 1; ++i) {
        if (mu1Pt >= ptBins[i] && mu2Pt >= ptBins[i]) {
          doubleMuRateOpenVsPtHist.Fill(ptBins[i]);
        }
      }
    }
  }

  std::cout << "Writing histograms to " << folder + identifier + "Hists.root"
            << std::endl;
  f->Open(folder + identifier + "Hists.root", "recreate");
  singleMuRateOpenVsPtHist.Write();
  doubleMuRateOpenVsPtHist.Write();
  std::cout << "Closing file." << std::endl;
  // f->Close();
  std::cout << "Done writing." << std::endl;
}

void calcRates(int nCollBunches, int nevents, TH1D& singleMuRateOpenVsPtHist,
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
  std::cout << "DoubleMuOpen rate: "
            << doubleMuRateOpenVsPtHist.GetBinContent(1) * norm << "+/-"
            << TMath::Sqrt(doubleMuRateOpenVsPtHist.GetBinContent(1)) * norm
            << "\n";
  std::cout << "SingleMuOpen rate: "
            << singleMuRateOpenVsPtHist.GetBinContent(1) * norm << "+/-"
            << TMath::Sqrt(singleMuRateOpenVsPtHist.GetBinContent(1)) * norm
            << "\n";
  std::cout << "###########################################"
            << "\n";

  doubleMuRateOpenVsPtHist.Scale(norm);
  singleMuRateOpenVsPtHist.Scale(norm);
}

void drawHistograms(TH1D& singleMuHist, TH1D& doubleMuHist, TString filename,
                    TString xAxisLabel, TString descString, TString plotFolder,
                    TString run, TString singleMuDesc, TString doubleMuDesc) {
  TLatex n1;
  n1.SetNDC();
  n1.SetTextFont(52);
  n1.SetTextSize(0.035);

  TLatex n2;
  n2.SetNDC();
  n2.SetTextFont(52);
  n2.SetTextSize(0.035);

  TCanvas c1("c1", "", 800, 800);
  TPad pad1("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1.SetLogy();
  pad1.SetBottomMargin(0);
  pad1.Draw();
  pad1.cd();
  // c1.SetLeftMargin(0.18);

  //  muRatesOpenUnpack->SetLineWidth(2);
  singleMuHist.SetMinimum(0.002);

  singleMuHist.SetLineColor(kGreen + 2);
  singleMuHist.SetMarkerStyle(32);
  singleMuHist.SetMarkerColor(kGreen + 2);
  singleMuHist.Draw("EP");
  // singleMuHist.GetXaxis()->SetTitle(xAxisLabel);
  singleMuHist.GetXaxis()->SetTitleFont(43);
  singleMuHist.GetXaxis()->SetTitleSize(23);
  singleMuHist.GetXaxis()->SetTitleOffset(3.5);
  singleMuHist.GetYaxis()->SetTitle("Rate [kHz]");

  // singleMuHist.GetYaxis()->SetTitleOffset(1.55);

  doubleMuHist.SetLineColor(kOrange + 7);
  doubleMuHist.SetMarkerStyle(27);
  doubleMuHist.SetMarkerColor(kOrange + 7);
  doubleMuHist.Draw("same,EP");
  // doubleMuHist.GetXaxis()->SetTitle(xAxisLabel);
  doubleMuHist.GetYaxis()->SetTitle("Rate [kHz]");
  // doubleMuHist.GetYaxis()->SetTitleOffset(1.5);

  gPad->Modified();

  TLegend leg1(0.5, 0.72, 0.9, 0.92);
  leg1.SetFillColor(0);
  leg1.AddEntry(&singleMuHist, singleMuDesc, "lp");
  leg1.AddEntry(&doubleMuHist, doubleMuDesc, "lp");
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.Draw();
  n1.DrawLatex(0.5, 0.68, "Runs " + run + ", #sqrt{s} = 13 TeV");
  n2.DrawLatex(0.5, 0.63, descString);

  // Do not draw the Y axis label on the upper plot and redraw a small
  // axis instead, in order to avoid the first label (0) to be clipped.
  // singleMuHist.GetYaxis()->SetLabelSize(0.);
  // TGaxis* axis = new TGaxis(-5, 20, -5, 220, 20, 220, 510, "");
  // axis->SetLabelFont(43);  // Absolute font size in pixel (precision
  // 3) axis->SetLabelSize(15); axis->Draw();

  // lower plot will be in pad
  c1.cd();  // Go back to the main canvas before defining pad2
  TPad pad2("pad2", "pad2", 0, 0.01, 1, 0.3);
  pad2.SetTopMargin(0);
  pad2.SetBottomMargin(0.25);
  pad2.Draw();
  pad2.cd();  // pad2 becomes the current pad

  TH1D* ratioplot = (TH1D*)singleMuHist.Clone("ratioplot");
  ratioplot->Sumw2();
  ratioplot->Divide(&doubleMuHist);
  ratioplot->Draw("ep");
  ratioplot->SetTitle("");
  ratioplot->SetMinimum(0);

  ratioplot->GetYaxis()->SetTitle("#frac{Single Muons}{Double Muons}");
  ratioplot->SetLineColor(kBlack);
  ratioplot->SetMarkerStyle(21);
  ratioplot->SetMarkerColor(kBlack);
  ratioplot->GetYaxis()->SetNdivisions(505);
  // ratioplot->GetYaxis()->SetTitleSize(20);
  // ratioplot->GetYaxis()->SetTitleFont(43);
  // ratioplot->GetYaxis()->SetTitleOffset(1.55);
  ratioplot->GetYaxis()->SetLabelFont(
      43);  // Absolute font size in pixel (precision 3)
  ratioplot->GetYaxis()->SetLabelSize(15);
  // X axis ratio plot settings
  ratioplot->GetXaxis()->SetTitle(xAxisLabel);
  ratioplot->GetXaxis()->SetTitleSize(20);
  ratioplot->GetXaxis()->SetTitleFont(43);
  ratioplot->GetXaxis()->SetTitleOffset(4.);
  ratioplot->GetXaxis()->SetLabelFont(
      43);  // Absolute font size in pixel (precision 3)
  ratioplot->GetXaxis()->SetLabelSize(15);

  c1.SaveAs(plotFolder + filename + ".pdf");
  c1.SaveAs(plotFolder + filename + ".png");
}
