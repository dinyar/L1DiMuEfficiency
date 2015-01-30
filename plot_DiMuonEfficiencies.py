from ROOT import *

gROOT.Reset()

etaScalePos= [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.75,1.8,1.85,1.9,1.95,2.,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4];
etaScale= [];
for i in range(31):
    etaScale[i]    = -etaScalePos[31-i]
    etaScale[32+i] = etaScalePos[i+1]
etaScale[31]=0;


hist_effVsDr       = TH1D("hist_effVsDr", "Efficiency vs. deltaR", 200, 0, 10)
hist_effVsDphi     = TH1D("hist_effVsDphi", "Efficiency vs. deltaPhi", 200, 0, 6)
hist_effVsDeta     = TH1D("hist_effVsDeta", "Efficiency vs. deltaEta", 200, 0, 5)

hist_effVsPt_DT    = TH1D("hist_effVsPt_DT", "Efficiency vs. pT, DT", 25, 0, 50)
hist_effVsPt_CSC   = TH1D("hist_effVsPt_CSC", "Efficiency vs. pT, CSC", 25, 0, 50)
hist_effVsPt_RPCb  = TH1D("hist_effVsPt_RPCb", "Efficiency vs. pT, RPCb", 25, 0, 50)
hist_effVsPt_RPCf  = TH1D("hist_effVsPt_RPCf", "Efficiency vs. pT, RPCf", 25, 0, 50)

hist_effVsPhi_DT   = TH1D("hist_effVsPhi_DT", "Efficiency vs. Phi, DT", 25, 0, 7)
hist_effVsPhi_CSC  = TH1D("hist_effVsPhi_CSC", "Efficiency vs. Phi, CSC", 25, 0, 7)
hist_effVsPhi_RPCb = TH1D("hist_effVsPhi_RPCb", "Efficiency vs. Phi, RPCb", 25, 0, 7)
hist_effVsPhi_RPCf = TH1D("hist_effVsPhi_RPCf", "Efficiency vs. Phi, RPCf", 25, 0, 7)

hist_effVsEta_DT   = TH1D("hist_effVsEta_DT", "Efficiency vs. Eta, DT", 62, etaScale)
hist_effVsEta_CSC  = TH1D("hist_effVsEta_CSC", "Efficiency vs. Eta, CSC", 62, etaScale)
hist_effVsEta_RPCb = TH1D("hist_effVsEta_RPCb", "Efficiency vs. Eta, RPCb", 62, etaScale)
hist_effVsEta_RPCf = TH1D("hist_effVsEta_RPCf", "Efficiency vs. Eta, RPCf", 62, etaScale)
