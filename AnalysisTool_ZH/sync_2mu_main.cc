// own
#include "AnalysisTool/Analyse.h"
// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
// C & C++
#include <iostream>
#include <vector>
#include <sstream>


using namespace std;
//You should derive your own class from Analyse.
class MyAnalysis : public Analyse {
private:
    // pile-up
    std::vector < double > dataPileUp;
    std::vector < double > dataPileUp_Up;
    std::vector < double > dataPileUp_Down;

    UInt_t currun;
    UInt_t curlumi;
    Int_t  hltmu9;

    // vertex
    int    cVtxNdf;
    double cVtxZ;
    //_________muon cuts_________
    double cPtMu;
    double cEtaMu;
    double cPtMuMax;
    double cEtaMuMax;
    double cChi2NdfMu;
    double cNHitsMuonMu;
    double cNMatchedStations;
    double cDxyMu;
    double cDzMu;
    int    cNHitsPixelMu;
    int    cNTrackerLayersMu;
    double cIsoMuPU;
    double cIsoMuTk;
    int    muidchoice;

    double cInvMassB_low;
    double cInvMassB_high;
    double c2MuPtCut;

    double cInvMass;
    double massZ;

    //_________jet cuts_________
    double cPtJet;
    double cEtaJet;
    double cDR;

    //_________sync cuts_________
    double cSyncWindowLow;
    double cSyncWindowHigh;
    double c_Pre_PtJetLeading;
    double c_Pre_PtJetSubleading;
    double c_Pre_METMax;
    double c_VBFTight_InvMassJetMin;
    double c_VBFTight_DiJetdEtaMin;
    double c_ggFTight_InvMassJetMin;
    double c_ggFTightDiMuPtMin;
    double c_01jet_DiMuPtMin;

    //_________MET cuts_________
    double cMET;


    // define outputfile
    TFile *histfile;
    // define the ROOT trees
    TTree *ftreeCat1;
    // define variables to be stored in ROOT tree
    int    trunNr;
    double tevNr;
    float  tevWt;
    float  tgenWt;
    float  tInvMass;

    TH1D *hVtxN;
    TH1D *hPUweight;

    // muons
    TH1D *hMuPt;
    TH1D *hLeadingMuPt;
    TH1D *hSubleadingMuPt;
    TH1D *hMuEta;
    TH1D *hLeadingMuEta;
    TH1D *hSubleadingMuEta;
    TH1D *hInvMass2Mu; // all combos

    TH1D *h2MuEta;
    TH1D *h2MuPt;

    TH1D *hInvMass2Mu_sync_0; // events in in window 100-110 gev
    TH1D *hInvMass2Mu_sync_1; // events in in window 100-110 gev
    TH1D *hInvMass2Mu_sync_2; // events in in window 100-110 gev
    TH1D *hInvMass2Mu_sync_3; // events in in window 100-110 gev
    TH1D *hInvMass2Mu_sync_4; // events in in window 100-110 gev
    TH1D *hInvMass2Mu_sync_5; // events in in window 100-110 gev

    TH1D *hNJets_sync_0;
    TH1D *h2JetM_sync_0;
    TH1D *h2JetDEta_sync_0;
    TH1D *h2MuPt_sync_0;
    TH1D *hNJets_sync_1;
    TH1D *h2JetM_sync_1;
    TH1D *h2JetDEta_sync_1;
    TH1D *h2MuPt_sync_1;
    TH1D *hNJets_sync_2;
    TH1D *h2JetM_sync_2;
    TH1D *h2JetDEta_sync_2;
    TH1D *h2MuPt_sync_2;
    TH1D *hNJets_sync_3;
    TH1D *h2JetM_sync_3;
    TH1D *h2JetDEta_sync_3;
    TH1D *h2MuPt_sync_3;
    TH1D *hNJets_sync_4;
    TH1D *h2JetM_sync_4;
    TH1D *h2JetDEta_sync_4;
    TH1D *h2MuPt_sync_4;
    TH1D *hNJets_sync_5;
    TH1D *h2JetM_sync_5;
    TH1D *h2JetDEta_sync_5;
    TH1D *h2MuPt_sync_5;

    // jets
    TH1D *hJetPt;
    TH1D *hJetEta;
    TH1D *hNJets;
    TH1D *hLeadingJetPt;
    TH1D *hSubleadingJetPt;
    TH1D *hLeadingJetEta;
    TH1D *hSubleadingJetEta;
    TH1D *h2JetM;
    TH1D *h2JetPt;
    TH1D *h2JetDEta;

    // met
    TH1D *hPFMETType1;
    TH1D *hPFMETPuppiType1;

    // efficiencies
    TH1D *hEfficiencies;
    TH1D *hEfficiencies_w;
    vector<TString> eff_names;
    vector<double *> eff_counters;
    vector<double *> eff_counters_w;


    // counters
    double nEv_Skim,      nEv_Skim_w;
    double nEv_TriggerMu, nEv_TriggerMu_w;
    double nEv_PV,        nEv_PV_w;
    double nEv_PVNdf,     nEv_PVNdf_w;
    double nEv_PVZ,       nEv_PVZ_w;
    double nEv_GAndTr,    nEv_GAndTr_w;
    double nEv_Pt,        nEv_Pt_w;
    double nEv_Eta,       nEv_Eta_w;
    double nEv_MuID,      nEv_MuID_w;
    double nEv_IsoPU,     nEv_IsoPU_w;
    double nEv_PtEtaMax,  nEv_PtEtaMax_w;
    double nEv_2SGMu,     nEv_2SGMu_w;
    double nEv_MuID2,     nEv_MuID2_w;
    double nEv_V3rdJet,   nEv_V3rdJet_w;
    double nEv_InvMassMu, nEv_InvMassMu_w;

    // debug
    bool debug;
    bool isMC;


public:
    MyAnalysis();
    virtual ~MyAnalysis();
    //AnalyseEvent is a virtual function which is called for each event.
    virtual Int_t AnalyseEvent();
    double getAoverBError(double nA, double nB);
    double getIsoPU(Muon mu, double rho);
    double getIsoPUE(Electron mu, double rho);
    void   AddPUWeightFile(string filename);
    bool   isGoodBJet(Jet j);
    bool   isGoodJet(Jet j);
    bool   hasJetIDLoose(Jet j);
    bool   hasJetSamePV(Jet j, Muon mu1, Muon mu2, std::vector < Vertex > PVs);
    double GetHLTEffScale();
    double GetMuEffScale(double muPt, double muEta);
    double GetElEffScale(double elPt, double elEta);
    string asString(double f);
    double sumWeights;
};
//Constructor:
MyAnalysis::MyAnalysis() : Analyse(), currun(0), curlumi(0) {
    // don't touch, these are changed with the output file name // // //
    isMC=false; sumWeights = 0.0;
    // // // //

    // vertex
    cVtxNdf = 4;
    cVtxZ   = 24.;

    cPtMu        = 10.0; // GeV;
    cPtMuMax     = 20.0;
    cEtaMu       = 2.4;  // bbZ: 2.1
    cEtaMuMax    = 2.4;  //

    cChi2NdfMu   = 10.0;
    cNHitsMuonMu = 0;
    cNMatchedStations = 1;
    cDxyMu            = 0.02; // cm
    cDzMu             = 0.14; // cm
    cNHitsPixelMu     = 0;
    cNTrackerLayersMu = 5;
    cIsoMuPU          = 0.25;
    cIsoMuTk          = 0.10;
    //muidchoice = 3; // medium
    muidchoice = 4; // tight

    cInvMassB_low  = 120.;
    cInvMassB_high = 130.;
    cInvMass       = 60.; //GeV
    c2MuPtCut	   = 38.; // GeV

    massZ  = 91.1876; //GeV

    // jet cuts
    cPtJet = 30.;  // GeV;
    cEtaJet = 2.4;
    cDR     = 0.5;

    // sync cuts
    cSyncWindowLow  = 100.;
    cSyncWindowHigh = 110.;
    c_Pre_PtJetLeading    = 40.;
    c_Pre_PtJetSubleading = 30.;
    c_Pre_METMax          = 40.;
    c_VBFTight_InvMassJetMin = 650.;
    c_VBFTight_DiJetdEtaMin  = 3.5;
    c_ggFTight_InvMassJetMin = 250.;
    c_ggFTightDiMuPtMin      = 50;
    c_01jet_DiMuPtMin = 10;

    // load needed informations
    LoadTrigger();
    LoadBeamSpot();
    LoadPrimVertices();
    LoadMuons();
    //LoadElectrons();
    //LoadTracks();
    LoadAK4PFCHSJets();
    LoadMET();
    //LoadGenParticles();
    //LoadGenInfo();
    //LoadAllGenParticles();
    UsePileUpInfo();

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IS IT DATA OR MONTE CARLO? /////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //isMC=true;
    isMC=false;

    // output file
    if (!isMC)    histfile = new TFile("2mu_data_ana_out.root", "RECREATE");
    else if (isMC) histfile = new TFile("2mu_MC_ana_out.root", "RECREATE");
    else histfile = new TFile("2mu_noinfo_ana_out.root", "RECREATE");
    histfile->cd();

    TVector3 zDir(0,0,1);
    // root tree:
    ftreeCat1 = new TTree("Category1", "Category1", 1);
    // variables to be stored into ROOT tree
    ftreeCat1->Branch("tevNr",        &tevNr,        "tevNr/I");
    ftreeCat1->Branch("trunNr",       &trunNr,       "trunNr/I");
    ftreeCat1->Branch("tInvMass",     &tInvMass,     "tInvMass/F");

    ftreeCat1->Branch("tevWt",        &tevWt,        "tevWt/F");
    ftreeCat1->Branch("tgenWt",       &tgenWt,       "tgenWt/F");

    // vertex
    hVtxN       = new TH1D("hVtxN", "N Vtx", 100, 0., 100.);
    hVtxN->GetXaxis()->SetTitle("N_{PV}");
    hVtxN->GetYaxis()->SetTitle("Candidates");

    hPUweight   = new TH1D("hPUweight", "PU weight", 100, 0., 10.);
    hPUweight->GetXaxis()->SetTitle("PU weight");
    hPUweight->GetYaxis()->SetTitle("Candidates");

    // muons
    hMuPt      = new TH1D("hMuPt", "mu Pt",    160, 0., 800.);
    hMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hMuEta     = new TH1D("hMuEta", "mu Eta",  44, -2.2, 2.2);
    hMuEta->GetXaxis()->SetTitle("#eta_{#mu}");
    hMuEta->GetYaxis()->SetTitle("Candidates/0.1");

    hLeadingMuPt  = new TH1D("hLeadingMuPt", "leading #mu Pt",    160, 0., 800.);
    hLeadingMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hLeadingMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hLeadingMuEta = new TH1D("hLeadingMuEta", "leading #mu Eta",  44, -2.2, 2.2);
    hLeadingMuEta->GetXaxis()->SetTitle("#eta_{#mu}");
    hLeadingMuEta->GetYaxis()->SetTitle("Candidates/0.1");

    hSubleadingMuPt  = new TH1D("hSubleadingMuPt", "subleading #mu Pt",    160, 0., 800.);
    hSubleadingMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hSubleadingMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hSubleadingMuEta = new TH1D("hSubleadingMuEta", "subleading #mu Eta",  44, -2.2, 2.2);
    hSubleadingMuEta->GetXaxis()->SetTitle("#eta_{#mu}");
    hSubleadingMuEta->GetYaxis()->SetTitle("Candidates/0.1");

    hInvMass2Mu    = new TH1D("hInvMass2Mu", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    h2MuEta        = new TH1D("h2MuEta", "2mu Eta",  132, -6.6, 6.6);
    h2MuEta->GetXaxis()->SetTitle("#eta_{#mu^{+}#mu^{-}}");
    h2MuEta->GetYaxis()->SetTitle("Candidates/0.1");

    h2MuPt         = new TH1D("h2MuPt", "2mu pT",  2000, 0, 1000);
    h2MuPt->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt->GetYaxis()->SetTitle("Candidates/0.5 GeV");

    // sync cat plots
    hInvMass2Mu_sync_0 = new TH1D("hInvMass2Mu_sync_0", "M mumu cat_0", 60, 90., 120.);
    hInvMass2Mu_sync_0->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_sync_0->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_sync_1 = new TH1D("hInvMass2Mu_sync_1", "M mumu cat_1", 60, 90., 120.);
    hInvMass2Mu_sync_1->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_sync_1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_sync_2 = new TH1D("hInvMass2Mu_sync_2", "M mumu cat_2", 60, 90., 120.);
    hInvMass2Mu_sync_2->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_sync_2->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_sync_3 = new TH1D("hInvMass2Mu_sync_3", "M mumu cat_3", 60, 90., 120.);
    hInvMass2Mu_sync_3->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_sync_3->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_sync_4 = new TH1D("hInvMass2Mu_sync_4", "M mumu cat_4", 60, 90., 120.);
    hInvMass2Mu_sync_4->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_sync_4->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_sync_5 = new TH1D("hInvMass2Mu_sync_5", "M mumu cat_5", 60, 90., 120.);
    hInvMass2Mu_sync_5->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_sync_5->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hNJets_sync_0 = new TH1D("hNJets_sync_0", "jet multiplicity cat_0", 40, 0., 40.);
    hNJets_sync_0->GetXaxis()->SetTitle("N_{jets}");
    hNJets_sync_0->GetYaxis()->SetTitle("Events");
    hNJets_sync_1 = new TH1D("hNJets_sync_1", "jet multiplicity cat_1", 40, 0., 40.);
    hNJets_sync_1->GetXaxis()->SetTitle("N_{jets}");
    hNJets_sync_1->GetYaxis()->SetTitle("Events");
    hNJets_sync_2 = new TH1D("hNJets_sync_2", "jet multiplicity cat_2", 40, 0., 40.);
    hNJets_sync_2->GetXaxis()->SetTitle("N_{jets}");
    hNJets_sync_2->GetYaxis()->SetTitle("Events");
    hNJets_sync_3 = new TH1D("hNJets_sync_3", "jet multiplicity cat_3", 40, 0., 40.);
    hNJets_sync_3->GetXaxis()->SetTitle("N_{jets}");
    hNJets_sync_3->GetYaxis()->SetTitle("Events");
    hNJets_sync_4 = new TH1D("hNJets_sync_4", "jet multiplicity cat_4", 40, 0., 40.);
    hNJets_sync_4->GetXaxis()->SetTitle("N_{jets}");
    hNJets_sync_4->GetYaxis()->SetTitle("Events");
    hNJets_sync_5 = new TH1D("hNJets_sync_5", "jet multiplicity cat_5", 40, 0., 40.);
    hNJets_sync_5->GetXaxis()->SetTitle("N_{jets}");
    hNJets_sync_5->GetYaxis()->SetTitle("Events");

    h2JetM_sync_0 = new TH1D("h2JetM_sync_0", "M jj cat_0", 4000, 0., 2000.);
    h2JetM_sync_0->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM_sync_0->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    h2JetM_sync_1 = new TH1D("h2JetM_sync_1", "M jj cat_1", 4000, 0., 2000.);
    h2JetM_sync_1->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM_sync_1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    h2JetM_sync_2 = new TH1D("h2JetM_sync_2", "M jj cat_2", 4000, 0., 2000.);
    h2JetM_sync_2->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM_sync_2->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    h2JetM_sync_3 = new TH1D("h2JetM_sync_3", "M jj cat_3", 4000, 0., 2000.);
    h2JetM_sync_3->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM_sync_3->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    h2JetM_sync_4 = new TH1D("h2JetM_sync_4", "M jj cat_4", 4000, 0., 2000.);
    h2JetM_sync_4->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM_sync_4->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    h2JetM_sync_5 = new TH1D("h2JetM_sync_5", "M jj cat_5", 4000, 0., 2000.);
    h2JetM_sync_5->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM_sync_5->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    h2JetDEta_sync_0 = new TH1D("h2JetDEta_sync_0", "2jet dEta cat_0",  100, -5, 5);
    h2JetDEta_sync_0->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta_sync_0->GetYaxis()->SetTitle("Candidates/0.1");
    h2JetDEta_sync_1 = new TH1D("h2JetDEta_sync_1", "2jet dEta cat_1",  100, -5, 5);
    h2JetDEta_sync_1->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta_sync_1->GetYaxis()->SetTitle("Candidates/0.1");
    h2JetDEta_sync_2 = new TH1D("h2JetDEta_sync_2", "2jet dEta cat_2",  100, -5, 5);
    h2JetDEta_sync_2->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta_sync_2->GetYaxis()->SetTitle("Candidates/0.1");
    h2JetDEta_sync_3 = new TH1D("h2JetDEta_sync_3", "2jet dEta cat_3",  100, -5, 5);
    h2JetDEta_sync_3->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta_sync_3->GetYaxis()->SetTitle("Candidates/0.1");
    h2JetDEta_sync_4 = new TH1D("h2JetDEta_sync_4", "2jet dEta cat_4",  100, -5, 5);
    h2JetDEta_sync_4->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta_sync_4->GetYaxis()->SetTitle("Candidates/0.1");
    h2JetDEta_sync_5 = new TH1D("h2JetDEta_sync_5", "2jet dEta cat_5",  100, -5, 5);
    h2JetDEta_sync_5->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta_sync_5->GetYaxis()->SetTitle("Candidates/0.1");

    h2MuPt_sync_0 = new TH1D("h2MuPt_sync_0", "2mu pT cat_0",  2000, 0, 1000);
    h2MuPt_sync_0->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt_sync_0->GetYaxis()->SetTitle("Candidates/0.5 GeV");
    h2MuPt_sync_1 = new TH1D("h2MuPt_sync_1", "2mu pT cat_1",  2000, 0, 1000);
    h2MuPt_sync_1->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt_sync_1->GetYaxis()->SetTitle("Candidates/0.5 GeV");
    h2MuPt_sync_2 = new TH1D("h2MuPt_sync_2", "2mu pT cat_2",  2000, 0, 1000);
    h2MuPt_sync_2->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt_sync_2->GetYaxis()->SetTitle("Candidates/0.5 GeV");
    h2MuPt_sync_3 = new TH1D("h2MuPt_sync_3", "2mu pT cat_3",  2000, 0, 1000);
    h2MuPt_sync_3->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt_sync_3->GetYaxis()->SetTitle("Candidates/0.5 GeV");
    h2MuPt_sync_4 = new TH1D("h2MuPt_sync_4", "2mu pT cat_4",  2000, 0, 1000);
    h2MuPt_sync_4->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt_sync_4->GetYaxis()->SetTitle("Candidates/0.5 GeV");
    h2MuPt_sync_5 = new TH1D("h2MuPt_sync_5", "2mu pT cat_5",  2000, 0, 1000);
    h2MuPt_sync_5->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuPt_sync_5->GetYaxis()->SetTitle("Candidates/0.5 GeV");

    // jets
    hJetPt         = new TH1D("hJetPt", "jet Pt",    160, 0., 800.);
    hJetPt->GetXaxis()->SetTitle("p_{T HJet}[GeV/c]");
    hJetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hJetEta        = new TH1D("hJetEta", "jet Eta",  52, -2.6, 2.6);
    hJetEta->GetXaxis()->SetTitle("#eta_{HJet}");
    hJetEta->GetYaxis()->SetTitle("Candidates/0.1");

    hNJets         = new TH1D("hNJets", "jet multiplicity", 40, 0., 40.);
    hNJets->GetXaxis()->SetTitle("N_{jets}");
    hNJets->GetYaxis()->SetTitle("Events");

    hLeadingJetPt         = new TH1D("hLeadingJetPt", "leading jet Pt",    160, 0., 800.);
    hLeadingJetPt->GetXaxis()->SetTitle("p_{T leading Jet}[GeV/c]");
    hLeadingJetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hLeadingJetEta        = new TH1D("hLeadingJetEta", "leading jet Eta",  52, -2.6, 2.6);
    hLeadingJetEta->GetXaxis()->SetTitle("#eta_{HLeadingJet}");
    hLeadingJetEta->GetYaxis()->SetTitle("Candidates/0.1");

    hSubleadingJetPt         = new TH1D("hSubleadingJetPt", "subleading jet Pt",    160, 0., 800.);
    hSubleadingJetPt->GetXaxis()->SetTitle("p_{T subleading Jet}[GeV/c]");
    hSubleadingJetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hSubleadingJetEta        = new TH1D("hSubleadingJetEta", "jet Eta",  52, -2.6, 2.6);
    hSubleadingJetEta->GetXaxis()->SetTitle("#eta_{HSubleadingJet}");
    hSubleadingJetEta->GetYaxis()->SetTitle("Candidates/0.1");

    h2JetM         = new TH1D("h2JetM", "M jj", 4000, 0., 2000.);
    h2JetM->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    h2JetPt        = new TH1D("h2JetPt", "2jet Pt",    160, 0., 800.);
    h2JetPt->GetXaxis()->SetTitle("p_{T jj}[GeV/c]");
    h2JetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    h2JetDEta     = new TH1D("h2JetDEta", "2jet dEta",  100, -5, 5);
    h2JetDEta->GetXaxis()->SetTitle("#Delta #eta_{leading jet - subl. jet}");
    h2JetDEta->GetYaxis()->SetTitle("Candidates/0.1");

    hPFMETType1         = new TH1D("hPFMETType1", "PFMETType1", 4000, 0., 2000.);
    hPFMETType1->GetXaxis()->SetTitle("MET [GeV/c^{2}]");
    hPFMETType1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hPFMETPuppiType1         = new TH1D("hPFMETPuppiType1", "PFMETType1", 4000, 0., 2000.);
    hPFMETPuppiType1->GetXaxis()->SetTitle("PuppiPuppi  MET [GeV/c^{2}]");
    hPFMETPuppiType1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EFFICIENCIES ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    string muidchoice_ = "none";
    if (muidchoice==4) muidchoice_ = "tight";
    else if (muidchoice==3) muidchoice_ = "medium";
    // make sure that these are listed in the order that they are implemented
    eff_names.push_back("Skim");                          eff_counters.push_back(&nEv_Skim); eff_counters_w.push_back(&nEv_Skim_w);
    eff_names.push_back("Trigger");                       eff_counters.push_back(&nEv_TriggerMu); eff_counters_w.push_back(&nEv_TriggerMu_w);
    eff_names.push_back("nPV > 0");                       eff_counters.push_back(&nEv_PV); eff_counters_w.push_back(&nEv_PV_w);
    eff_names.push_back("PVNdf > "+asString(cVtxNdf));    eff_counters.push_back(&nEv_PVNdf); eff_counters_w.push_back(&nEv_PVNdf_w);
    eff_names.push_back("PVZ < "+asString(cVtxZ));        eff_counters.push_back(&nEv_PVZ); eff_counters_w.push_back(&nEv_PVZ_w);
    eff_names.push_back("Global+Tracker muon");           eff_counters.push_back(&nEv_GAndTr); eff_counters_w.push_back(&nEv_GAndTr_w);
    eff_names.push_back("mu pT > "+asString(cPtMu));      eff_counters.push_back(&nEv_Pt); eff_counters_w.push_back(&nEv_Pt_w);
    eff_names.push_back("mu eta < "+asString(cEtaMu));    eff_counters.push_back(&nEv_Eta); eff_counters_w.push_back(&nEv_Eta_w);
    eff_names.push_back("single mu ID");                  eff_counters.push_back(&nEv_MuID); eff_counters_w.push_back(&nEv_MuID_w);
    eff_names.push_back("mu iso PU < "+asString(cIsoMuPU)); eff_counters.push_back(&nEv_IsoPU); eff_counters_w.push_back(&nEv_IsoPU_w);
    eff_names.push_back("at least 1 mu with pT > "+asString(cPtMuMax)+", eta < "+asString(cEtaMuMax)); eff_counters.push_back(&nEv_PtEtaMax); eff_counters_w.push_back(&nEv_PtEtaMax_w);
    eff_names.push_back("2 tight or loose mus");          eff_counters.push_back(&nEv_2SGMu); eff_counters_w.push_back(&nEv_2SGMu_w);
    eff_names.push_back("2 "+muidchoice_+" mus");          eff_counters.push_back(&nEv_MuID2); eff_counters_w.push_back(&nEv_MuID2_w);
    eff_names.push_back("nJets =< 2");                    eff_counters.push_back(&nEv_V3rdJet); eff_counters_w.push_back(&nEv_V3rdJet_w);
    eff_names.push_back("dimuon pair");                   eff_counters.push_back(&nEv_InvMassMu); eff_counters_w.push_back(&nEv_InvMassMu_w);

    // set up efficiencies histogram
    hEfficiencies    = new TH1D("hEfficiencies", "Efficiencies", eff_counters.size(), 0., (double)eff_counters.size());
    hEfficiencies_w  = new TH1D("hEfficiencies_w", "Efficiencies_w", eff_counters_w.size(), 0., (double)eff_counters_w.size());
    hEfficiencies   ->GetYaxis()->SetTitle("Events");
    hEfficiencies_w ->GetYaxis()->SetTitle("Events");
    for(unsigned int i = 0; i < eff_counters.size(); i++) {
        *eff_counters[i] = 0;
        *eff_counters_w[i] = 0;
        hEfficiencies   ->GetXaxis()->SetBinLabel((i+1) , eff_names[i]);
        hEfficiencies_w ->GetXaxis()->SetBinLabel((i+1) , eff_names[i]);
    }

    //debug = true;
    debug = false;

    if (debug) { 
        cerr<<"*\n*\n*\n*\n*\n*\n*\n*\n*\n*\n*\n*"<<endl;
        cerr<<"Warning: debug is true. outfile will be giant"<<endl;
        cerr<<"*\n*\n*\n*\n*\n*\n*\n*\n*\n*\n*\n*\n"<<endl;
    }
}
// Destructor:
MyAnalysis::~MyAnalysis() {
    // set up efficiencies histogram
    for(unsigned int i = 0; i < eff_counters.size(); i++) {
        hEfficiencies->SetBinContent((i+1), *eff_counters[i]);
        hEfficiencies_w->SetBinContent((i+1), *eff_counters_w[i]);
    }

    histfile->Write();
    histfile->Close();
    //print out
    cout<<"====================================================================================="<<endl;
    cout<<setw(40)<<"Selection"<<setw(15)<<setprecision(8)<<"Events"<<setw(15)<<"Eff.(Skim) [%]"<<setw(15)<<"Rel.Eff. [%]"<<endl;
    cout<<"-------------------------------------------------------------------------------------"<<endl;
    cout<<setw(40)<<eff_names[0]<<setw(15)<<(*eff_counters[0])<<setw(15)<<" --- "<<setw(15)<<" --- "<<endl;
    for(unsigned int i = 1; i < eff_names.size(); i++) {
        cout<<setw(40)<<eff_names[i]<<setw(15)<<setprecision(8)<<(*eff_counters[i])<<setw(15)<<setprecision(4)<<100.*(*eff_counters[i])/(*eff_counters[0])<<setw(15)<<100.*(*eff_counters[i])/(*eff_counters[i-1])<<endl;
    }
    cout<<"====================================================================================="<<endl;
    cerr<<setprecision(18)<<"******************** sumWeights = "<<sumWeights<<endl;
    cout<<setprecision(18)<<"******************** sumWeights = "<<sumWeights<<endl;
    cout<<setprecision(18)<<"******************** nEvents    = "<<*eff_counters[0]<<endl;
}

// Analysis
Int_t MyAnalysis::AnalyseEvent() {

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PILEUP WEIGHTING ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    double pileupweight = 1.;
    if (isMC) {
        pileupweight =  GenWeight() * GetPileUpWeight(dataPileUp);
        sumWeights += GenWeight();
    }
    ++nEv_Skim; nEv_Skim_w += pileupweight;

    hPUweight->Fill(pileupweight);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TRIGGER ////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    vector<string> triggernames;
    if (triggernames.size()) triggernames.clear();

    // Turn on triggers
    bool useHLTIsoMu20 = true;
    //bool useHLTIsoMu24 = false;
    bool useHLTIsoTkMu20 = true;

    if (useHLTIsoMu20) {
        triggernames.push_back("HLT_IsoMu20");
        triggernames.push_back("HLT_IsoMu20_v1");
        triggernames.push_back("HLT_IsoMu20_v2");
        triggernames.push_back("HLT_IsoMu20_v3");
        triggernames.push_back("HLT_IsoMu20_v4");
    }
    if (useHLTIsoTkMu20) {
        triggernames.push_back("HLT_IsoTkMu20");
        triggernames.push_back("HLT_IsoTkMu20_v1");
        triggernames.push_back("HLT_IsoTkMu20_v2");
        triggernames.push_back("HLT_IsoTkMu20_v3");
        triggernames.push_back("HLT_IsoTkMu20_v4");
    }

    Int_t trigger = GetHLTrigger(triggernames);
        if (debug) {
            for (int i = 0; i < GetNumHLTriggers(); ++i) {
                cout<<"GetHLTriggerName(i):"<<GetHLTriggerName(i)<<endl;
                cout<<"GetHLTriggerIndex(GetHLTriggerName(i)): "<<GetHLTriggerIndex(GetHLTriggerName(i))<<endl;
                cout<<endl;
            }
        }
    // check trigger
    if (trigger == -1) {
        cout << "Trigger Problem" << endl;
        if (debug) {
            for(int i = 0; i < GetNumHLTriggers(); ++i) {
                cout<<"GetHLTriggerName(i):"<<GetHLTriggerName(i)<<endl;
            }
        }
        return(1);
    } else if (trigger == 0) {
        if (debug) cout<<"trigger did not fire!"<<endl;
        return(1);
    }
    if (debug) cout<<"trigger fired!!! whooooooooooooooooooooooooooooooooooooooooooo"<<endl;
    ++nEv_TriggerMu; nEv_TriggerMu_w += pileupweight;


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PRIMARY VERTEX /////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // select event with a valid PV
    // loop over prim vertices
    std::vector < Vertex > PVs;
    bool isVtxNdfOK = false;
    bool isVtxZOK   = false;
    for(unsigned int v = 0; v < NumPrimVertices(); ++v) {
        if(!isVtxNdfOK) isVtxNdfOK = PrimVertices(v).Ndof() > cVtxNdf;
        if(!isVtxZOK) isVtxZOK = TMath::Abs(PrimVertices(v).Z()) < cVtxZ;
        if(!(isVtxNdfOK && isVtxZOK)) continue;
        PVs.push_back(PrimVertices(v));
    }
    if(!(PVs.size() > 0)) return (1);
    hVtxN->Fill(PVs.size(), pileupweight);

    if(isVtxNdfOK) ++nEv_PVNdf; nEv_PVNdf_w += pileupweight;
    if(isVtxZOK)   ++nEv_PVZ;   nEv_PVZ_w += pileupweight;
    if(PVs.size() < 1) return (1);
    ++nEv_PV; nEv_PV_w += pileupweight;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MUONS //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // select event with 2 good muons
    // get the muons
    std::vector< Muon > muonsVec;
    bool isGAndTr        = false;
    bool isPtCutOK       = false;
    bool isEtaCutOK      = false;
    int  nMuPtEtaMax     = 0;
    bool isChi2NdfOK     = false;
    bool isNMuonHitsOK   = false;
    bool isNMatchedStationsOK = false;
    bool isDxyOK              = false;
    bool isDzOK               = false;
    bool isNPixelHitsOK       = false;
    bool isNTrackerLayersOK   = false;
    bool isMuIDOk             = false;
    bool isIsoPUOK            = false;
    bool isIsoTkOK            = false;
    int nLooseMus = 0;
    int nTightMus = 0;
    int nMedMus   = 0;
    int nVetoMus  = 0;

    // loop over the muons
    for(unsigned int i = 0; i < NumMuons(); ++i) {
        // all muons pass these cuts:
        if(!(Muons(i).IsGlobal() && Muons(i).IsTracker())) continue;
        ////if( !Muons(i).isPFMuon())continue;
        isGAndTr = true;
        if(!(Muons(i).Pt() > cPtMu)) continue;
        isPtCutOK = true;
        if(!(TMath::Abs(Muons(i).Eta()) <= cEtaMu)) continue;
        isEtaCutOK = true;
        // collect muon ID
        if (Muons(i).MuID() <= 1) {
            nVetoMus++;
            continue;
        }

        if (Muons(i).MuID() > 2) { isMuIDOk = true; }

        if (Muons(i).MuID() == 4) nTightMus++;
        if (Muons(i).MuID() >= 3) nMedMus++;
        if (Muons(i).MuID() >= 2) nLooseMus++;
        if (!(this->getIsoPU(Muons(i), AK4PFRho()) < cIsoMuPU)) continue; isIsoPUOK = true;
        //if (!(Muons(i).IsoR3TrackRel() < cIsoMuTk)) continue; isIsoTkOK = true;
        // at least one muon in the event must ALSO pass this cut:
        if((Muons(i).Pt() > cPtMuMax && TMath::Abs(Muons(i).Eta()) <= cEtaMuMax)) ++nMuPtEtaMax;
        // save the good mus
        muonsVec.push_back(Muons(i));
    }

    if(isGAndTr)   { ++nEv_GAndTr;        nEv_GAndTr_w += pileupweight; }
    if(isPtCutOK)  { ++nEv_Pt;            nEv_Pt_w += pileupweight; }
    if(isEtaCutOK) { ++nEv_Eta;           nEv_Eta_w += pileupweight; }
    if(isMuIDOk)   { ++nEv_MuID;          nEv_MuID_w += pileupweight; }
    if(isIsoPUOK)  { ++nEv_IsoPU; nEv_IsoPU_w += pileupweight; }
    //if(isIsoTkOK)  { ++nEv_IsoTk;         nEv_IsoTk_w += pileupweight; }
    if(nMuPtEtaMax > 0) { ++nEv_PtEtaMax; nEv_PtEtaMax_w += pileupweight; }
    if(nMuPtEtaMax < 1) return(1);

    // select events with at least 2 good mus
    if(muonsVec.size() < 2) return(1);
    ++nEv_2SGMu; nEv_2SGMu_w += pileupweight;

    //if (nTightMus >= 1) { ++nEv_MuID1; nEv_MuID1_w += pileupweight; }


    // ask for them both to have same ID
    if (muidchoice == 4) {
        if (nTightMus >= 2) { ++nEv_MuID2; nEv_MuID2_w += pileupweight; }
        else return (1);
    } else if (muidchoice == 3) {
        if (nMedMus >= 2) { ++nEv_MuID2; nEv_MuID2_w += pileupweight; }
        else return (1);
    } else {
        cerr<<"no muidchoice"<<endl;
        return (1);
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // JETS ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // for each jet, check against all of the e and mu to make sure it is far enough away
    std::vector< Jet > jetsVec;
    for (unsigned int jet = 0; jet < NumAK4PFCHSJets(); ++jet) {
        bool isdROK = true;
        for (unsigned int i = 0; i < muonsVec.size(); ++i) {
            double dR = AK4PFCHSJets(jet).DeltaR(muonsVec[i]);
            if(dR < cDR) isdROK = false;
        }
        //if (!isdROK) continue;
        if (!(this->hasJetIDLoose(AK4PFCHSJets(jet)))) continue;
        if (!(AK4PFCHSJets(jet).Pt() > cPtJet)) continue;
        if (!(TMath::Abs(AK4PFCHSJets(jet).Eta()) <= cEtaJet)) continue;
        // store the jets
        jetsVec.push_back(AK4PFCHSJets(jet));
    }

    // veto on the 3rd jet
    //if(jetsVec.size() > 2)     return(1);
    //++nEv_V3rdJet;


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MET ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    hPFMETType1->Fill(PFMETTYPE1().Pt(), pileupweight);
    hPFMETPuppiType1->Fill(PFMETPUPPITYPE1().Pt(), pileupweight);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PAIRS //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // DIMUON PAIRS ///////////////////////////////////////////////////////////////////////////////
    bool isChargeMuCutOK  = false;
    bool isSamePVMuCutOK  = false;
    bool isInvMassMuCutOK = false;
    std::vector< TLorentzVector > diMuonVec;
    std::vector< int > Muon_idx1;
    std::vector< int > Muon_idx2;

    // 1st loop over muons
    for (unsigned int i = 0; i < muonsVec.size(); ++i) {
        // 2nd loop over muons
        for (unsigned int j = i+1; j < muonsVec.size(); ++j) {
            // check muon charge
            if (muonsVec[i].Charge()*muonsVec[j].Charge() > 0) continue;
            isChargeMuCutOK = true;     
            // dz cut
            double dzMu = muonsVec[i].Dz() - muonsVec[j].Dz();
            if (!(TMath::Abs(dzMu) < cDzMu)) continue;
            isSamePVMuCutOK = true;
            // check the invariant mass
            TLorentzVector diMuon = muonsVec[i] + muonsVec[j];
            if (!(diMuon.M() > cInvMass)) continue;
            isInvMassMuCutOK = true;

            diMuonVec.push_back(diMuon);
            // make two arrays: if i and j make a pair, put each index in a separate array
            // so Muon_idx1[0] and Muon_idx2[0] make a pair, and Muon_idx1[1] and Muon_idx2[1] make another, and so on
            Muon_idx1.push_back(i);
            Muon_idx2.push_back(j);
        } // end 2nd loop over muons
    } //end 1st loop over muons

    if (!isSamePVMuCutOK) return(1);
    //++nEv_SamePVMu; nEv_SamePVMu_w += pileupweight;
    if (!isChargeMuCutOK) return(1);
    //++nEv_ChargeMu; nEv_ChargeMu_w += pileupweight;
    if (!isInvMassMuCutOK) return(1);
    //++nEv_InvMassMu; nEv_InvMassMu_w += pileupweight;

    // only look at 100-110 GeV for now
    if (!(diMuonVec[0].M() <= cSyncWindowHigh && diMuonVec[0].M() >= cSyncWindowLow)) { return (1); }
    ++nEv_InvMassMu; nEv_InvMassMu_w += pileupweight;

    // DIJET PAIRS ////////////////////////////////////////////////////////////////////////////////
    std::vector< TLorentzVector > diJetVec;
    for(unsigned int i = 0; i < jetsVec.size(); ++i) {
        for(unsigned int j = i+1; j < jetsVec.size(); ++j) {
            TLorentzVector diJet = jetsVec[i] + jetsVec[j];
    //        if(!(diJet.M() > cInvMass)) continue;
            h2JetM->Fill(diJet.M(), pileupweight);
            h2JetPt->Fill(diJet.Pt(), pileupweight);
            if(jetsVec[i].Eta()*jetsVec[j].Eta() > 0) continue;
            diJetVec.push_back(diJet);
        }
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Selection is done - FILL PLOTS /////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // fill muon plots
    for (unsigned int i = 0; i < muonsVec.size(); ++i) {
        hMuPt       ->Fill(muonsVec[i].Pt(),  pileupweight);
        hMuEta      ->Fill(muonsVec[i].Eta(), pileupweight);
    }

    if (muonsVec.size() > 0) {
        hLeadingMuPt->Fill(muonsVec[0].Pt(),  pileupweight);
        hLeadingMuEta->Fill(muonsVec[0].Eta(),  pileupweight);
    }
    if (muonsVec.size() > 1) {
        hSubleadingMuPt->Fill(muonsVec[1].Pt(),  pileupweight);
        hSubleadingMuEta->Fill(muonsVec[1].Eta(),  pileupweight);
    }

    hInvMass2Mu  ->Fill(diMuonVec[0].M(), pileupweight);
    h2MuEta  ->Fill(diMuonVec[0].Eta(), pileupweight);
    h2MuPt  ->Fill(diMuonVec[0].Pt(), pileupweight);

    // fill jet plots
    for(unsigned int jet = 0; jet < jetsVec.size(); ++jet) {
        hJetPt ->Fill(jetsVec[jet].Pt() , pileupweight);
        hJetEta->Fill(jetsVec[jet].Eta(), pileupweight);
    }
    hNJets->Fill(jetsVec.size(), pileupweight);
    if (diJetVec.size()) {
        h2JetM  ->Fill(diJetVec[0].M(), pileupweight);
        h2JetPt ->Fill(diJetVec[0].Pt(), pileupweight);
    }

    if (jetsVec.size() > 0) {
        hLeadingJetPt->Fill(jetsVec[0].Pt(), pileupweight);
        hLeadingJetEta->Fill(jetsVec[0].Eta(), pileupweight);
    }
    if (jetsVec.size() > 1) {
        hSubleadingJetPt->Fill(jetsVec[1].Pt(), pileupweight);
        hSubleadingJetEta->Fill(jetsVec[1].Eta(), pileupweight);
        h2JetDEta->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
    }


    // fill sync plots
    bool passesPreselection = false;
    bool passesVBFTight = false;
    bool passesggFTight = false;
    bool passes01jetTight = false;


    // which selections does it pass?
    if (jetsVec.size() > 1 && jetsVec[0].Pt() > c_Pre_PtJetLeading && jetsVec[1].Pt() > c_Pre_PtJetSubleading && PFMETTYPE1().Pt() < c_Pre_METMax) passesPreselection = true;
    if (diJetVec.size() && diJetVec[0].M() > c_VBFTight_InvMassJetMin && TMath::Abs(jetsVec[0].Eta() - jetsVec[1].Eta()) > c_VBFTight_DiJetdEtaMin) passesVBFTight = true;
    if (diJetVec.size() && diJetVec[0].M() > c_ggFTight_InvMassJetMin && diMuonVec[0].Pt() > c_ggFTightDiMuPtMin) passesggFTight = true;
    if (diMuonVec[0].Pt() > c_01jet_DiMuPtMin) passes01jetTight = true;



    // all events in the window that get sent to preselection
    if (diMuonVec[0].M() < cSyncWindowHigh && diMuonVec[0].M() > cSyncWindowLow) {
        hNJets_sync_0->Fill(jetsVec.size(), pileupweight);
        if (diJetVec.size()) {
            h2JetM_sync_0->Fill(diJetVec[0].M(), pileupweight);
            h2JetDEta_sync_0->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
        }
        h2MuPt_sync_0->Fill(diMuonVec[0].Pt(), pileupweight);
        hInvMass2Mu_sync_0->Fill(diMuonVec[0].M(), pileupweight);
    }

    int MyEventCategory = 0;

    if (diMuonVec[0].M() < cSyncWindowHigh && diMuonVec[0].M() > cSyncWindowLow) { // only look at 100-110 GeV for now
        if (passesPreselection) {
            if (passesVBFTight) {
                MyEventCategory = 1;
            } else if (passesggFTight) { // fails VBFTight
                MyEventCategory = 2;
            } else { // fails ggFTight
                MyEventCategory = 3;
            }
        } else { // fails preselection
            if (passes01jetTight) {
                MyEventCategory = 4;
            } else { // fails 01jetTight
                MyEventCategory = 5;
            }
        }
    }


    if (diMuonVec[0].M() < cSyncWindowHigh && diMuonVec[0].M() > cSyncWindowLow) {
        if (MyEventCategory==1) {
            cout<<setprecision(15)<<"SYNC_CAT1: Run number "<<Run()<<", event number "<<Number()<<endl;
            hNJets_sync_1->Fill(jetsVec.size(), pileupweight);
            if (diJetVec.size()) {
                h2JetM_sync_1->Fill(diJetVec[0].M(), pileupweight);
                h2JetDEta_sync_1->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
            }
            h2MuPt_sync_1->Fill(diMuonVec[0].Pt(), pileupweight);
            hInvMass2Mu_sync_1->Fill(diMuonVec[0].M(), pileupweight);
        } else if (MyEventCategory==2) {
            cout<<setprecision(15)<<"SYNC_CAT2: Run number "<<Run()<<", event number "<<Number()<<endl;
            hNJets_sync_2->Fill(jetsVec.size(), pileupweight);
            if (diJetVec.size()) {
                h2JetM_sync_2->Fill(diJetVec[0].M(), pileupweight);
                h2JetDEta_sync_2->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
            }
            h2MuPt_sync_2->Fill(diMuonVec[0].Pt(), pileupweight);
            hInvMass2Mu_sync_2->Fill(diMuonVec[0].M(), pileupweight);
        } else if (MyEventCategory==3) {
            cout<<setprecision(15)<<"SYNC_CAT3: Run number "<<Run()<<", event number "<<Number()<<endl;
            hNJets_sync_3->Fill(jetsVec.size(), pileupweight);
            if (diJetVec.size()) {
                h2JetM_sync_3->Fill(diJetVec[0].M(), pileupweight);
                h2JetDEta_sync_3->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
            }
            h2MuPt_sync_3->Fill(diMuonVec[0].Pt(), pileupweight);
            hInvMass2Mu_sync_3->Fill(diMuonVec[0].M(), pileupweight);
        } else if (MyEventCategory==4) {
            cout<<setprecision(15)<<"SYNC_CAT4: Run number "<<Run()<<", event number "<<Number()<<endl;
            hNJets_sync_4->Fill(jetsVec.size(), pileupweight);
            if (diJetVec.size()) {
                h2JetM_sync_4->Fill(diJetVec[0].M(), pileupweight);
                h2JetDEta_sync_4->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
            }
            h2MuPt_sync_4->Fill(diMuonVec[0].Pt(), pileupweight);
            hInvMass2Mu_sync_4->Fill(diMuonVec[0].M(), pileupweight);
        } else if (MyEventCategory==5) {
            cout<<setprecision(15)<<"SYNC_CAT5: Run number "<<Run()<<", event number "<<Number()<<endl;
            hNJets_sync_5->Fill(jetsVec.size(), pileupweight);
            if (diJetVec.size()) {
                h2JetM_sync_5->Fill(diJetVec[0].M(), pileupweight);
                h2JetDEta_sync_5->Fill(jetsVec[0].Eta() - jetsVec[1].Eta(), pileupweight);
            }
            h2MuPt_sync_5->Fill(diMuonVec[0].Pt(), pileupweight);
            hInvMass2Mu_sync_5->Fill(diMuonVec[0].M(), pileupweight);
        } else {
            cerr<<"sync_cat logic error B"<<endl;
        }
    }

    // fill the tree for limits calculations
    trunNr   = Run();
    tevNr    = Number();
    tgenWt   = GenWeight();
    tevWt    = pileupweight;
    tInvMass = diMuonVec[0].M();
    ftreeCat1->Fill();

    return(1);
}

int main() {
    //Create an instance of your analysis class.
    MyAnalysis ana;
    string namebuf;
    string filename;
    cout<<"Enter one filename and press enter. To stop this enter END exactly."
        << endl;
    while(true) {
        cin >> namebuf;
        if(namebuf == "END") {
            break;
        } else {
            UInt_t slashpos = namebuf.find_last_of("/");
            if(slashpos == namebuf.size()) {
                filename = namebuf;
            } else {
                filename = namebuf.substr(slashpos+1);
            }
            cout << filename << endl;
// NOTE
            if(filename.find("LUMI_") == 0) {
                ana.AddLumiFile(namebuf.c_str());
                //}else if(filename.find("ReWeight1D") == 0){
                //ana.AddPUWeightFile(namebuf.c_str());
            } else {
                ana.AddFile(namebuf.c_str());
            }
        }
    }
        //ana.AddWeightFile("/afs/cern.ch/work/e/ekennedy/work/fsanalysis/CMSSW_7_2_5/src/AnalysisTool/AnalysisTool_ZH/ReWeight1DWhole.root");
        ana.AddPUWeightFile("/afs/cern.ch/work/e/ekennedy/work/fsanalysis/CMSSW_7_2_5/src/AnalysisTool/AnalysisTool_ZH/MyDataPileupHistogram_v2.root");
    // Loop will start to run the analysis on the specified range or on
    // all events if no range is given.
    ana.SetPrintInfo(10000);
    //ana.EnableDuplicateCheck();
    //ana.Loop(0,5000);
    ana.Loop();
    ana.PrintLumiOfRuns();
    cout << "Lumi: " << ana.GetLumi() << endl;
}
//________________________________
bool MyAnalysis::isGoodJet(Jet j) {
    bool isOk = false;
    bool isPtOk  = j.Pt() >= cPtJet;
    bool isEtaOk = TMath::Abs(j.Eta()) <= cEtaJet;

    isOk = isPtOk && isEtaOk && this->hasJetIDLoose(j);
    return isOk;
}
//________________________________
bool MyAnalysis::hasJetIDLoose(Jet j) {
    bool chargedHadFrac = j.ChargedHadEnergyFraction() > 0.;
    bool neutralHadFrac = ((j.HadEnergyFraction() - j.ChargedHadEnergyFraction()) < 0.99);
    bool chargedEMFrac  = j.ChargedEMEnergyFraction() < 0.99;
    bool neutralEMFrac  = ((j.EMEnergyFraction() - j.ChargedEMEnergyFraction()) < 0.99);
    bool nConstituents  = ((j.ChargedMulti() + j.NeutralMulti()) > 1);
    bool nCharged       = j.ChargedMulti() > 0;

    bool jetID = (chargedHadFrac && neutralHadFrac && chargedEMFrac && neutralEMFrac  && nConstituents && nCharged);
    return jetID;
    //return j.PU_Jet_full_loose();
}
//________________________________
double MyAnalysis::getAoverBError(double nA, double nB) {
    double e_nA = TMath::Sqrt(nA);
    double e_nB = TMath::Sqrt(nB);
    double nB4  = TMath::Power(nB,4);
    double nB2e_nA2_p_nA2e_nB2 = nB*nB*e_nA*e_nA + nA*nA*e_nB*e_nB;
    if(nB4 == 0) {
        cout<<"error nB:"<<nB<<endl;
        return 0.;
    }
    if(nA == nB) {
        return 0.;
    }
    double e_AoverB = TMath::Sqrt(nB2e_nA2_p_nA2e_nB2/nB4);
    if(nA/nB + e_AoverB > 1) {
        return (1-nA/nB);
    }

    return e_AoverB;
}
//________________________________
double MyAnalysis::getIsoPU(Muon mu, double rho) {
    double corrPU = mu.IsoR3ECal() + mu.IsoR3HCal() - rho*Pi()*0.09;
    if(corrPU < 0) {
        corrPU = 0.;
    }
    return (mu.IsoR3Track() + corrPU)/mu.Pt();
}

//________________________________
double MyAnalysis::getIsoPUE(Electron el, double rho) {
    double corrPU = el.IsoR3ECal() + el.IsoR3HCal() - rho*Pi()*0.09;
    if(corrPU < 0) corrPU = 0.;
    return (el.IsoR3Track() + corrPU)/el.Pt();
}
//________________________________
void MyAnalysis::AddPUWeightFile(string filename) {
cout<<"AddPUWeightFile"<<endl;
    TFile weightfile(filename.c_str(), "READ");
    if(weightfile.IsZombie()) {
        cerr << "ERROR AddLumiFile: " << filename << " is not a valid file."<<endl;
    }
    TH1D *Vertex      = (TH1D *)weightfile.Get("pileup");
    //TH1D *Vertex_Up   = (TH1D *)weightfile.Get("pileup72450");
    //TH1D *Vertex_Down = (TH1D *)weightfile.Get("pileup65550");

    Double_t tot      = Vertex->Integral();
    //Double_t tot_up   = Vertex_Up->Integral();
    //Double_t tot_down = Vertex_Down->Integral();

    for(int i=0; i<Vertex->GetNbinsX(); i++) {
        dataPileUp.push_back(Vertex->GetBinContent(i+1)/tot);
    }

    //for(int i=0; i<Vertex_Up->GetNbinsX(); i++) {
    //    dataPileUp_Up.push_back(Vertex_Up->GetBinContent(i+1)/tot_up);
    //}

    //for(int i=0; i<Vertex_Down->GetNbinsX(); i++) {
    //    dataPileUp_Down.push_back(Vertex_Down->GetBinContent(i+1)/tot_down);
    //}
    weightfile.Close();
}
//________________________________
// work in progress
double MyAnalysis::GetMuEffScale(double muPt, double muEta) {
    double my_ratio = 1.;
    double eff_data;
    double eff_MC;
    return my_ratio;
}
//________________________________
// work in progress
double MyAnalysis::GetElEffScale(double elPt, double elEta) {
    double my_ratio = 1.;
    double eff_data;
    double eff_MC;
    return my_ratio;
}
/*
//________________________________
double MyAnalysis::GetMCWeight(int nPVMC, string nameMC, string nameData) {
//cerr<<"GetMCWeight"<<endl;
    double my_weight = 1.;
    if(nameMC == "DYJetsToLL" && nameData == "Run2015") return w2015_DY[nPVMC];
    if(nameMC == "TTJets" && nameData  == "Run2015") return w2015_TT[nPVMC];
//    if(nameMC == "QCD" && nameData == "Run2015")        return w2015_QCD[nPVMC];
//    if(nameMC == "Tbar_s" && nameData  == "Run2015")    return w2015_Tbar_s[nPVMC];
//    if(nameMC == "Tbar_t" && nameData  == "Run2015")    return w2015_Tbar_t[nPVMC];
//    if(nameMC == "Tbar_tW" && nameData == "Run2015")    return w2015_Tbar_tW[nPVMC];
//    if(nameMC == "T_s" && nameData  == "Run2015")       return w2015_T_s[nPVMC];
//    if(nameMC == "T_t" && nameData  == "Run2015")       return w2015_T_t[nPVMC];
//    if(nameMC == "T_tW" && nameData == "Run2015")       return w2015_T_tW[nPVMC];
//    if(nameMC == "WJetsToLNu" && nameData == "Run2015") return w2015_WJets[nPVMC];
//    if(nameMC == "WW" && nameData == "Run2015")         return w2015_WW[nPVMC];
//    if(nameMC == "WZ" && nameData == "Run2015")         return w2015_WZ[nPVMC];
//    if(nameMC == "ZZ" && nameData == "Run2015")         return w2015_ZZ[nPVMC];
//    if(nameMC == "H" && nameData == "Run2015")          return w2015_H[nPVMC];
cerr<<"pileupweight = "<<my_weight<<endl;
    return my_weight;
}
//________________________________
void MyAnalysis::AddWeightFile(string filename) {
cerr<<"AddWeightFile"<<endl;

    TFile *weightfile = new TFile(filename.c_str());
    if(weightfile->IsZombie()) {
        cerr << "ERROR AddLumiFile: " << filename << " is not a valid file."<<endl;
    }

    TH1D *nPV_Run2015     = (TH1D *) weightfile->Get("pileup");

    TH1D *nPV_DYJetsToLL  = (TH1D *) weightfile->Get("nPV_DYJetsToLL");
    TH1D *nPV_TTJets      = (TH1D *) weightfile->Get("nPV_TTJets");

//    TH1D *nPV_QCD         = (TH1D *) weightfile->Get("nPV_QCD");
//    TH1D *nPV_Tbar_s      = (TH1D *) weightfile->Get("nPV_Tbar_s");
//    TH1D *nPV_Tbar_t      = (TH1D *) weightfile->Get("nPV_Tbar_t");
//    TH1D *nPV_Tbar_tW     = (TH1D *) weightfile->Get("nPV_Tbar_tW");
//    TH1D *nPV_T_s         = (TH1D *) weightfile->Get("nPV_T_s");
//    TH1D *nPV_T_t         = (TH1D *) weightfile->Get("nPV_T_t");
//    TH1D *nPV_T_tW        = (TH1D *) weightfile->Get("nPV_T_tW");
//    TH1D *nPV_WJets       = (TH1D *) weightfile->Get("nPV_WJetsToLNu");
//    TH1D *nPV_WW          = (TH1D *) weightfile->Get("nPV_WW");
//    TH1D *nPV_WZ          = (TH1D *) weightfile->Get("nPV_WZ");
//    TH1D *nPV_ZZ          = (TH1D *) weightfile->Get("nPV_ZZ");
//    TH1D *nPV_H           = (TH1D *) weightfile->Get("nPV_Higgs");


    nPV_Run2015   ->Scale(1./(nPV_Run2015   ->Integral()));

    nPV_DYJetsToLL->Scale(1./(nPV_DYJetsToLL->Integral()));
    nPV_TTJets    ->Scale(1./(nPV_TTJets    ->Integral()));

//    nPV_QCD       ->Scale(1./(nPV_QCD       ->Integral()));
//    nPV_Tbar_s    ->Scale(1./(nPV_Tbar_s    ->Integral()));
//    nPV_Tbar_t    ->Scale(1./(nPV_Tbar_t    ->Integral()));
//    nPV_Tbar_tW   ->Scale(1./(nPV_Tbar_tW   ->Integral()));
//    nPV_T_s       ->Scale(1./(nPV_T_s       ->Integral()));
//    nPV_T_t       ->Scale(1./(nPV_T_t       ->Integral()));
//    nPV_T_tW      ->Scale(1./(nPV_T_tW      ->Integral()));
//    nPV_WJets     ->Scale(1./(nPV_WJets     ->Integral()));
//    nPV_WW        ->Scale(1./(nPV_WW        ->Integral()));
//    nPV_WZ        ->Scale(1./(nPV_WZ        ->Integral()));
//    nPV_ZZ        ->Scale(1./(nPV_ZZ        ->Integral()));
//    nPV_H         ->Scale(1./(nPV_H         ->Integral()));

    int nBins = nPV_Run2015->GetXaxis()->GetNbins();
    double sum2015  = 0.;
//TH1D *nPV_DYJetsToLL_ratio = (TH1D *)weightfile->Get("nPV_DYJetsToLL_ratio");
//TH1D *nPV_TTJets_ratio = (TH1D *)weightfile->Get("nPV_TTJets_ratio");

// these are plots of the scale factors (mc/data)
TH1D *nPV_Run2015_ratio =    (TH1D*)(nPV_Run2015)->Clone("nPV_Run2015_ratio");
TH1D *nPV_DYJetsToLL_ratio = (TH1D*)(nPV_DYJetsToLL)->Clone("nPV_DYJetsToLL_ratio");
TH1D *nPV_TTJets_ratio =     (TH1D*)(nPV_TTJets)->Clone("nPV_TTJets_ratio");

    // loop over number of PV
    for(int ibin = 0; ibin < nBins; ++ibin) {
        double nPV2015       = nPV_Run2015   ->GetBinContent(nPV_Run2015   ->GetXaxis()->FindBin(ibin));
        double nPVDYJetsToLL = nPV_DYJetsToLL->GetBinContent(nPV_DYJetsToLL->GetXaxis()->FindBin(ibin));
        double nPVTTJets     = nPV_TTJets    ->GetBinContent(nPV_TTJets    ->GetXaxis()->FindBin(ibin));

//        double nPVQCD        = nPV_QCD       ->GetBinContent(nPV_QCD       ->GetXaxis()->FindBin(ibin));
//        double nPVTbar_s     = nPV_Tbar_s    ->GetBinContent(nPV_Tbar_s    ->GetXaxis()->FindBin(ibin));
//        double nPVTbar_t     = nPV_Tbar_t    ->GetBinContent(nPV_Tbar_t    ->GetXaxis()->FindBin(ibin));
//        double nPVTbar_tW    = nPV_Tbar_tW   ->GetBinContent(nPV_Tbar_tW   ->GetXaxis()->FindBin(ibin));
//        double nPVT_s        = nPV_Tbar_s    ->GetBinContent(nPV_T_s       ->GetXaxis()->FindBin(ibin));
//        double nPVT_t        = nPV_Tbar_t    ->GetBinContent(nPV_T_t       ->GetXaxis()->FindBin(ibin));
//        double nPVT_tW       = nPV_Tbar_tW   ->GetBinContent(nPV_T_tW      ->GetXaxis()->FindBin(ibin));
//        double nPVWJets      = nPV_WJets     ->GetBinContent(nPV_WJets     ->GetXaxis()->FindBin(ibin));
//        double nPVWW         = nPV_WW        ->GetBinContent(nPV_WW        ->GetXaxis()->FindBin(ibin));
//        double nPVWZ         = nPV_WZ        ->GetBinContent(nPV_WZ        ->GetXaxis()->FindBin(ibin));
//        double nPVZZ         = nPV_ZZ        ->GetBinContent(nPV_ZZ        ->GetXaxis()->FindBin(ibin));
//        double nPVH          = nPV_H         ->GetBinContent(nPV_H         ->GetXaxis()->FindBin(ibin));

        if(nPV2015 != 0) {
            nPV_Run2015_ratio->SetBinContent(ibin,nPV2015/nPV2015);
        } else {
            nPV_Run2015_ratio->SetBinContent(ibin, nPV2015);
        }

        if(nPVDYJetsToLL != 0) {
            w2015_DY.push_back(nPV2015/nPVDYJetsToLL);
            //w2015_DY.push_back(nPV2015/nPVDYJetsToLL);
            //nPV_DYJetsToLL_ratio->SetBinContent(ibin, nPV2015/nPVDYJetsToLL);
        } else {
            w2015_DY.push_back(nPVDYJetsToLL);
            //nPV_DYJetsToLL_ratio->SetBinContent(ibin, nPVDYJetsToLL);
        }

        if(nPVTTJets !=0) {
            w2015_TT.push_back(nPV2015/nPVTTJets);
            //nPV_TTJets_ratio->SetBinContent(ibin, nPV2015/nPVTTJets);
        } else {
            w2015_TT.push_back(nPVTTJets);
            //nPV_TTJets_ratio->SetBinContent(ibin, nPVTTJets);
        }


//        if(nPVQCD !=0) {
//            w2015_QCD.push_back(nPV2015/nPVQCD);
//        } else {
//            w2015_QCD.push_back(nPVQCD);
//        }
//        if(nPVTbar_s !=0) {
//            w2015_Tbar_s.push_back(nPV2015/nPVTbar_s);
//        } else {
//            w2015_Tbar_s.push_back(nPVTbar_s);
//        }
//        if(nPVTbar_t !=0) {
//            w2015_Tbar_t.push_back(nPV2015/nPVTbar_t);
//        } else {
//            w2015_Tbar_t.push_back(nPVTbar_t);
//        }
//        if(nPVTbar_tW !=0) {
//            w2015_Tbar_tW.push_back(nPV2015/nPVTbar_tW);
//        } else {
//            w2015_Tbar_tW.push_back(nPVTbar_tW);
//        }
//        if(nPVT_s !=0) {
//            w2015_T_s.push_back(nPV2015/nPVT_s);
//        } else {
//            w2015_T_s.push_back(nPVT_s);
//        }
//        if(nPVT_t !=0) {
//            w2015_T_t.push_back(nPV2015/nPVT_t);
//        } else {
//            w2015_T_t.push_back(nPVT_t);
//        }
//        if(nPVT_tW !=0) {
//            w2015_T_tW.push_back(nPV2015/nPVT_tW);
//        } else {
//            w2015_T_tW.push_back(nPVT_tW);
//        }
//        if(nPVWJets !=0) {
//            w2015_WJets.push_back(nPV2015/nPVWJets);
//        } else {
//            w2015_WJets.push_back(nPVWJets);
//        }
//        if(nPVWW !=0) {
//            w2015_WW.push_back(nPV2015/nPVWW);
//        } else {
//            w2015_WW.push_back(nPVWW);
//        }
//        if(nPVWZ !=0) {
//            w2015_WZ.push_back(nPV2015/nPVWZ);
//        } else {
//            w2015_WZ.push_back(nPVWZ);
//        }
//        if(nPVZZ !=0) {
//            w2015_ZZ.push_back(nPV2015/nPVZZ);
//        } else {
//            w2015_ZZ.push_back(nPVZZ);
//        }
//        if(nPVH !=0) {
//            w2015_H.push_back(nPV2015/nPVH);
//        } else {
//            w2015_H.push_back(nPVH);
//        }

        sum2015  += nPV2015;
    //w2015_DY.push_back(nPV_DYJetsToLL_ratio->GetBinContent(ibin));
    //w2015_TT.push_back(nPV_TTJets_ratio->GetBinContent(ibin));

    } // end loop over number of PV

    // loop over number of PV
//sum2015 = nPV_Run2015->Integral();
cerr<<"sum2015 = "<<sum2015<<endl;
    for(int ibin = 0; ibin < nBins; ++ibin) {
        w2015_DY[ibin] /= sum2015;
        w2015_TT[ibin] /= sum2015;
//        w2015_Tbar_s[ibin] /= sum2015;
//        w2015_Tbar_t[ibin] /= sum2015;
//        w2015_Tbar_tW[ibin] /= sum2015;
//        w2015_T_s[ibin] /= sum2015;
//        w2015_T_t[ibin] /= sum2015;
//        w2015_T_tW[ibin] /= sum2015;
//        w2015_QCD[ibin] /= sum2015;
//        w2015_WJets[ibin] /= sum2015;
//        w2015_WW[ibin] /= sum2015;
//        w2015_WZ[ibin] /= sum2015;
//        w2015_ZZ[ibin] /= sum2015;
//        w2015_H[ibin] /= sum2015;
    } // end of loop over the number of PV
    weightfile->Close();
}
*/
// _________________________________________________
string MyAnalysis::asString(double f) {
    std::ostringstream oss;
    oss<<f;
    string s = oss.str();
    return s;
}
