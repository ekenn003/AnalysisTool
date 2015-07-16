// own
//#include "AnalysisTool/Analyse.h"
#include "AnalysisTool/Analyse.h"
// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
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

    // pile-up weight
//
    vector< double > w2012_DY;
    vector< double > w2012_TT;
    vector< double > w2012_QCD;
    vector< double > w2012_Tbar_s;
    vector< double > w2012_Tbar_t;
    vector< double > w2012_Tbar_tW;
    vector< double > w2012_T_s;
    vector< double > w2012_T_t;
    vector< double > w2012_T_tW;
    vector< double > w2012_WJets;
    vector< double > w2012_WW;
    vector< double > w2012_WZ;
    vector< double > w2012_ZZ;
    vector< double > w2012_H;

    UInt_t currun;
    UInt_t curlumi;
    Int_t  hltmu9;
    /*
        // cut variables
        double etaB;
        double etaO;
        double etaE;
        double etaBarrel;
    */
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

    double cInvMassB_low;
    double cInvMassB_high;
    double c2MuPtCut;

    //_________electron cuts_________
    double cPtE;
    double cEtaE1;
    double cEtaE2;
    double cEtaE3;

    double cInvMass;
    double massZ;

    //_________jet cuts_________
    double cPtJet;
    double cEtaJet;
    double cDR;

    //_________MET cuts_________
    double cMET;


    // define outputfile
    TFile *histfile;
    // define the ROOT trees
    TTree *ftreeCat1;
    // define variables to be stored in ROOT tree
    int    trunNr;
    double tevNr;
    float  tInvMass;

    TH1D *hVtxN;
    TH1D *hVtxN_u;
    TH1D *hVtxNdf;
    TH1D *hVtxR;
    TH1D *hVtxZ;
    TH1D *hPUweight;

    // muons
    TH1D *hMuPt;
    TH1D *hMuEta;
    TH1D *hMuPhi;
    TH1D *hMuIso;
    TH1D *hMuIsoPU;
    TH1D *hRho;
    TH1D *hMuChi2Ndf;
    TH1D *hMuChi2NdfIT;
    TH1D *hMuITHits;
    TH1D *hMuOTHits;
    TH1D *hMuDxy;
    TH1D *hMuDz;
    TH1D *hMuPHits;
    TH1D *hMuDPtZ;
    TH1D *hMuDEtaZ;
    TH1D *hMuDPhiZ;
    TH1D *h2MuPtZ;
    TH1D *h2MuEtaZ;
    TH1D *h2MuPhiZ;
    TH1D *hMuDPtH;
    TH1D *hMuDEtaH;
    TH1D *hMuDPhiH;
    TH1D *h2MuPtH;
    TH1D *h2MuEtaH;
    TH1D *h2MuPhiH;
    TH1D *hDPhi2MuZ2MuH;
    TH1D *hInvMass2MuZ;
    TH1D *hInvMass2MuH;
    // Baseline histos
    TH1D *hInvMass2Mu; // all combos
    TH1D *hInvMass2Mu_wo; // pick two

    TH1D *hInvMass2Mu_A_cat1;
    TH1D *hInvMass2Mu_A_cat2;
    TH1D *hInvMass2Mu_A_cat3; 

    TH1D *hInvMass2Mu_B_cat1;
    TH1D *hInvMass2Mu_B_cat2;
    TH1D *hInvMass2Mu_B_cat3; 

    TH1D *hInvMass2Mu_C_cat1;
    TH1D *hInvMass2Mu_C_cat2;
    TH1D *hInvMass2Mu_C_cat3; 

    TH1D *h2MuEta;
    TH1D *hInvMass4Mu;

    // electrons
//    TH1D *hEPt;
//    TH1D *hEEta;
//    TH1D *hEPhi;
//    TH1D *hEIsoPU;
//    TH1D *hEChi2Ndf;
//    TH1D *hEITHits;
//    TH1D *hEDxy;
//    TH1D *hEDz;
//    TH1D *hEPHits;
//    TH1D *hEDPt;
//    TH1D *hEDEta;
//    TH1D *hEDPhi;
//    TH1D *h2EPt;
//    TH1D *h2EEta;
//    TH1D *h2EPhi;
//    TH1D *hDPhi2E2Mu;

    // jets
    TH1D *hJetPt;
    TH1D *hJetEta;
    TH1D *hJetPhi;
    TH1D *hNJets;
    TH1D *h2JetM;
    TH1D *h2JetPt;

    TH1D *hEfficiencies;
    TH1D *hAngleProd;
    TH1D *hCosAngleProd;

    // efficiencies
    vector<TString> eff_names;
    vector<double *> eff_counters;

    // pileup
    double weight;

    // counters
    double nEv_Skim;
    double nEv_GenZ2Mu;
    double nEv_TriggerSMu;
    double nEv_PV;
    double nEv_PVNdf;
    double nEv_PVZ;
    double nEv_MuID;
    double nEv_MET;
    double nEv_GAndTr;
    double nEv_Pt;
    double nEv_Eta;
    double nEv_PtEtaMax;
    double nEv_Chi2Ndf;
    double nEv_MuonHits;
    double nEv_MatchedStations;
    double nEv_Dxy;
    double nEv_Dz;
    double nEv_PixelHits;
    double nEv_TrackerLayers;
    double nEv_IsoPU;
//    double nEv_GSF;
//    double nEv_PtE;
//    double nEv_EtaE;
//    double nEv_Chi2NdfE;
//    double nEv_DxyE;
//    double nEv_DzE;
//    double nEv_PixelHitsE;
//    double nEv_TrackerLayersE;
//    double nEv_IsoPUE;
    double nEv_2SGMu;
    double nEv_4SGMu;
    double nEv_1TMu;
    double nEv_2TMu;
    double nEv_3TMu;
    double nEv_4TMu;
//    double nEv_2SGE;
//    double nEv_3rdMuVeto;
//    double nEv_3rdEVeto;
//    double nEv_SamePVE;
//    double nEv_ChargeE;
//    double nEv_InvMassE;
    double nEv_SamePVMu;
    double nEv_ChargeMu;
    double nEv_InvMassMu;
    double nEv_etaGapJetVeto;

TString trigNameString;

    // debug
    bool debug;
    bool isData;
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
    void   AddWeightFile(string filename);
    bool   isGoodBJet(Jet j);
    bool   isGoodJet(Jet j);
    bool   hasJetIDLoose(Jet j);
    bool   hasJetSamePV(Jet j, Muon mu1, Muon mu2, std::vector < Vertex > PVs);
    double GetMCWeight(int nPVMC, string nameMC, string nameData);

};
//Constructor:
MyAnalysis::MyAnalysis() : Analyse(), currun(0), curlumi(0) {
    // don;t touch, these are changed with the output file name
    isData=false;
    isMC=false;

    // vertex
    cVtxNdf = 4;
    cVtxZ   = 24.;

    cPtMu        = 10.0; // GeV;
    cPtMuMax     = 20.0; // Run2012B pT>30GeV; Run2012 pT>24GeV;
    cEtaMu       = 2.4;  // bbZ: 2.1
    cEtaMuMax    = 2.4;  //

    cChi2NdfMu   = 10.0;
    cNHitsMuonMu = 0;
    cNMatchedStations = 1;
    cDxyMu            = 0.02; // cm
    cDzMu             = 0.1; // cm
    cNHitsPixelMu     = 0;
    cNTrackerLayersMu = 5;
    cIsoMuPU          = 0.12;//0.10; // corrected for pile-up

    cInvMassB_low     = 120.;
    cInvMassB_high    = 130.;
    cInvMass          = 60.; //GeV
    c2MuPtCut	    = 38.; // GeV

    massZ             = 91.1876; //GeV
    cPtE   = 20.; //GeV
    cEtaE1 = 1.4442;
    cEtaE2 = 1.566;
    cEtaE3 = 2.5;

    // jet cuts
    cPtJet       = 30.;  // GeV;
    cEtaJet      = 2.4;  //
    cDR          = 0.5;

    // load needed informations
    LoadTrigger();
    LoadBeamSpot();
    LoadPrimVertices();
    LoadMuons();
    LoadElectrons();
    LoadTracks();
    LoadAK4PFJets();
    LoadMET();
    LoadGenParticles();
    LoadAllGenParticles();
    UsePileUpInfo();
    //UsePrimVertexInfo();

    // output file
    histfile = new TFile("4mu_"+trigNameString+"_test_mini.root", "RECREATE");
    histfile->cd();
    //isData=true;
    isMC=true;

    TVector3 zDir(0,0,1);
    // root tree:
    ftreeCat1 = new TTree("Category1", "Category1", 1);
    // variables to be stored into ROOT tree
    ftreeCat1->Branch("tevNr",        &tevNr,        "tevNr/I");
    ftreeCat1->Branch("trunNr",       &trunNr,       "trunNr/I");
    ftreeCat1->Branch("tInvMass",     &tInvMass,     "tInvMass/F");


    // vertex
    hVtxN       = new TH1D("hVtxN", "N Vtx", 100, 0., 100.);
    hVtxN->GetXaxis()->SetTitle("N_{PV}");
    hVtxN->GetYaxis()->SetTitle("Candidates");
    hVtxN_u     = new TH1D("hVtxN_u", "N Vtx", 100, 0., 100.);
    hVtxN_u->GetXaxis()->SetTitle("N_{PV}");
    hVtxN_u->GetYaxis()->SetTitle("Candidates");
    hVtxNdf     = new TH1D("hVtxNdf", "Vtx Ndf", 200, 0., 200.);
    hVtxNdf->GetXaxis()->SetTitle("Ndf_{PV}");
    hVtxNdf->GetYaxis()->SetTitle("Candidates");
    hVtxZ       = new TH1D("hVtxZ", "Vtx Z", 60, -30., 30.);
    hVtxZ->GetXaxis()->SetTitle("z_{PV} [cm]");
    hVtxZ->GetYaxis()->SetTitle("Candidates/1.0 cm");
    hPUweight   = new TH1D("hPUweight", "PU weight", 100, 0., 10.);
    hPUweight->GetXaxis()->SetTitle("PU weight");
    hPUweight->GetYaxis()->SetTitle("Candidates");

//    // electron
//    hEPt        = new TH1D("hEPt", "e Pt",    160, 0., 800.);
//    hEPt->GetXaxis()->SetTitle("p_{T e}[GeV/c]");
//    hEPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
//    hEEta       = new TH1D("hEEta", "e Eta",  44, -2.2, 2.2);
//    hEEta->GetXaxis()->SetTitle("#eta_{e}");
//    hEEta->GetYaxis()->SetTitle("Candidates/0.1");
//    hEPhi       = new TH1D("hEPhi", "e Phi", 34, -3.4, 3.4);
//    hEPhi->GetXaxis()->SetTitle("#varphi_{e} [rad]");
//    hEPhi->GetYaxis()->SetTitle("Candidates/0.2[rad]");
//
//    hEIsoPU    = new TH1D("hEIsoPU", "E IsoPUCombR", 300, -0.15, 0.15);
//    hEIsoPU->GetXaxis()->SetTitle("IsoPU^{PU}_{e}");
//    hEIsoPU->GetYaxis()->SetTitle("Candidates/0.001");
//    hEChi2Ndf  = new TH1D("hEChi2Ndf", "E Chi2Ndf", 100, 0., 10.);
//    hEChi2Ndf->GetXaxis()->SetTitle("#chi^{2}/Ndf_{e}");
//    hEChi2Ndf->GetYaxis()->SetTitle("Candidates/0.01");
//    hEITHits   = new TH1D("hEITHits", "E ITHits", 40, 0., 40.);
//    hEITHits->GetXaxis()->SetTitle("Hits_{e IT}");
//    hEITHits->GetYaxis()->SetTitle("Candidates");
//    hEPHits    = new TH1D("hEPHits", "E PHits", 12, 0., 12.);
//    hEPHits->GetXaxis()->SetTitle("Hits_{e Pixel}");
//    hEPHits->GetYaxis()->SetTitle("Candidates");
//    hEDxy      = new TH1D("hEDxy", "E Dxy", 100, -0.02, 0.02);
//    hEDxy->GetXaxis()->SetTitle("d_{xy e} [cm]");
//    hEDxy->GetYaxis()->SetTitle("Candidates/0.0004 [cm]");
//    hEDz       = new TH1D("hEDz", "E Dz", 2000, -1., 1.);
//    hEDz->GetXaxis()->SetTitle("d_{z e} [cm]");
//    hEDz->GetYaxis()->SetTitle("Candidates/0.001 [cm]");


    // mu
    hMuPt      = new TH1D("hMuPt", "mu Pt",    160, 0., 800.);
    hMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hMuEta     = new TH1D("hMuEta", "mu Eta",  44, -2.2, 2.2);
    hMuEta->GetXaxis()->SetTitle("#eta_{#mu}");
    hMuEta->GetYaxis()->SetTitle("Candidates/0.1");
    hMuPhi     = new TH1D("hMuPhi", "mu Phi", 34, -3.4, 3.4);
    hMuPhi->GetXaxis()->SetTitle("#varphi_{#mu} [rad]");
    hMuPhi->GetYaxis()->SetTitle("Candidates/0.2[rad]");

    hMuIso       = new TH1D("hMuIso", "Mu IsoCombR", 150, 0., 0.15);
    hMuIso->GetXaxis()->SetTitle("Iso^{PU}_{#mu}");
    hMuIso->GetYaxis()->SetTitle("Candidates/0.001");
    hMuIsoPU     = new TH1D("hMuIsoPU", "Mu IsoPUCombR", 300, -0.15, 0.15);
    hMuIsoPU->GetXaxis()->SetTitle("IsoPU^{PU}_{#mu}");
    hMuIsoPU->GetYaxis()->SetTitle("Candidates/0.001");
    hRho         = new TH1D("hRho", "Rho", 100, 0., 100.);
    hRho->GetXaxis()->SetTitle("#rho [GeV/c]");
    hRho->GetYaxis()->SetTitle("Candidates/1.0");
    hMuChi2Ndf   = new TH1D("hMuChi2Ndf", "Mu Chi2Ndf", 100, 0., 10.);
    hMuChi2Ndf->GetXaxis()->SetTitle("#chi^{2}/Ndf_{#mu}");
    hMuChi2Ndf->GetYaxis()->SetTitle("Candidates/0.01");
    hMuChi2NdfIT = new TH1D("hMuChi2NdfIT", "Mu Chi2Ndf", 100, 0., 10.);
    hMuChi2NdfIT->GetXaxis()->SetTitle("#chi^{2}/Ndf_{#mu}");
    hMuChi2NdfIT->GetYaxis()->SetTitle("Candidates/0.01");
    hMuITHits    = new TH1D("hMuITHits", "Mu ITHits", 40, 0., 40.);
    hMuITHits->GetXaxis()->SetTitle("Hits_{#mu IT}");
    hMuITHits->GetYaxis()->SetTitle("Candidates");
    hMuOTHits    = new TH1D("hMuOTHits", "Mu OTHits", 60, 0., 60.);
    hMuOTHits->GetXaxis()->SetTitle("Hits_{#mu OT}");
    hMuOTHits->GetYaxis()->SetTitle("Candidates");
    hMuPHits     = new TH1D("hMuPHits", "Mu OTHits", 12, 0., 12.);
    hMuPHits->GetXaxis()->SetTitle("Hits_{#mu Pixel}");
    hMuPHits->GetYaxis()->SetTitle("Candidates");
    hMuDxy       = new TH1D("hMuDxy", "Mu Dxy", 100, -0.02, 0.02);
    hMuDxy->GetXaxis()->SetTitle("d_{xy #mu} [cm]");
    hMuDxy->GetYaxis()->SetTitle("Candidates/0.0004 [cm]");
    hMuDz        = new TH1D("hMuDz", "Mu Dz", 2000, -1., 1.);
    hMuDz->GetXaxis()->SetTitle("d_{z #mu} [cm]");
    hMuDz->GetYaxis()->SetTitle("Candidates/0.001 [cm]");

    hMuDPtZ    = new TH1D("hMuDPtZ", "mu Pt",    320, -800., 800.);
    hMuDPtZ->GetXaxis()->SetTitle("#Delta p_{T #mu^{+} - #mu^{-}}[GeV/c]");
    hMuDPtZ->GetYaxis()->SetTitle("Candidates/5.0GeV");
    hMuDEtaZ   = new TH1D("hMuDEtaZ", "mu Eta",  88, -4.4, 4.4);
    hMuDEtaZ->GetXaxis()->SetTitle("#Delta #eta_{#mu^{+} - #mu^{-}}");
    hMuDEtaZ->GetYaxis()->SetTitle("Candidates/0.1");
    hMuDPhiZ   = new TH1D("hMuDPhiZ", "mu Phi", 34, -3.4, 3.4);
    hMuDPhiZ->GetXaxis()->SetTitle("#Delta #varphi_{#mu^{+} - #mu^{-}} [rad]");
    hMuDPhiZ->GetYaxis()->SetTitle("Candidates/0.2[rad]");
    h2MuPtZ    = new TH1D("h2MuPtZ", "2mu Pt",    160, 0., 800.);
    h2MuPtZ->GetXaxis()->SetTitle("p_{T #mu^{+}#mu^{-}}[GeV/c]");
    h2MuPtZ->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    h2MuEtaZ   = new TH1D("h2MuEtaZ", "2mu Eta",  132, -6.6, 6.6);
    h2MuEtaZ->GetXaxis()->SetTitle("#eta_{#mu^{+}#mu^{-}}");
    h2MuEtaZ->GetYaxis()->SetTitle("Candidates/0.1");
    h2MuPhiZ   = new TH1D("h2MuPhiZ", "2mu Phi", 34, -3.4, 3.4);
    h2MuPhiZ->GetXaxis()->SetTitle("#varphi_{#mu^{+}#mu^{-}} [rad]");
    h2MuPhiZ->GetYaxis()->SetTitle("Candidates/0.2[rad]");

    hMuDPtH    = new TH1D("hMuDPtH", "mu Pt",    320, -800., 800.);
    hMuDPtH->GetXaxis()->SetTitle("#Delta p_{T #mu^{+} - #mu^{-}}[GeV/c]");
    hMuDPtH->GetYaxis()->SetTitle("Candidates/5.0GeV");
    hMuDEtaH   = new TH1D("hMuDEtaH", "mu Eta",  88, -4.4, 4.4);
    hMuDEtaH->GetXaxis()->SetTitle("#Delta #eta_{#mu^{+} - #mu^{-}}");
    hMuDEtaH->GetYaxis()->SetTitle("Candidates/0.1");
    hMuDPhiH   = new TH1D("hMuDPhiH", "mu Phi", 34, -3.4, 3.4);
    hMuDPhiH->GetXaxis()->SetTitle("#Delta #varphi_{#mu^{+} - #mu^{-}} [rad]");
    hMuDPhiH->GetYaxis()->SetTitle("Candidates/0.2[rad]");
    h2MuPtH    = new TH1D("h2MuPtH", "2mu Pt",    160, 0., 800.);
    h2MuPtH->GetXaxis()->SetTitle("p_{T #mu^{+}#mu^{-}}[GeV/c]");
    h2MuPtH->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    h2MuEtaH   = new TH1D("h2MuEtaH", "2mu Eta",  132, -6.6, 6.6);
    h2MuEtaH->GetXaxis()->SetTitle("#eta_{#mu^{+}#mu^{-}}");
    h2MuEtaH->GetYaxis()->SetTitle("Candidates/0.1");
    h2MuPhiH   = new TH1D("h2MuPhiH", "2mu Phi", 34, -3.4, 3.4);
    h2MuPhiH->GetXaxis()->SetTitle("#varphi_{#mu^{+}#mu^{-}} [rad]");
    h2MuPhiH->GetYaxis()->SetTitle("Candidates/0.2[rad]");

    hDPhi2MuZ2MuH  = new TH1D("hDPhi2MuZ2MuH", " ", 34, -3.4, 3.4);
    hDPhi2MuZ2MuH->GetXaxis()->SetTitle("#Delta#varphi_{{#mu^{+}#mu^{-} - {#mu^{+}#mu^{-}} [rad]");
    hDPhi2MuZ2MuH->GetYaxis()->SetTitle("Candidates/0.2[rad]");

    // pairs
    hInvMass2MuZ   = new TH1D("hInvMass2MuZ", "M mumu (Z)", 4000, 0., 2000.);
    hInvMass2MuZ->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2MuZ->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2MuH   = new TH1D("hInvMass2MuH", "M mumu (H)", 4000, 0., 2000.);
    hInvMass2MuH->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2MuH->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hInvMass2Mu    = new TH1D("hInvMass2Mu", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_wo = new TH1D("hInvMass2Mu_wo", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu_wo->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_wo->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    h2MuEta        = new TH1D("h2MuEta", "2mu Eta",  132, -6.6, 6.6);
    h2MuEta->GetXaxis()->SetTitle("#eta_{#mu^{+}#mu^{-}}");
    h2MuEta->GetYaxis()->SetTitle("Candidates/0.1");
    hInvMass4Mu    = new TH1D("hInvMass4Mu", "M mumu", 4000, 0., 2000.);
    hInvMass4Mu->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass4Mu->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    // jets
    hJetPt         = new TH1D("hJetPt", "jet Pt",    160, 0., 800.);
    hJetPt->GetXaxis()->SetTitle("p_{T HJet}[GeV/c]");
    hJetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hJetEta        = new TH1D("hJetEta", "jet Eta",  52, -2.6, 2.6);
    hJetEta->GetXaxis()->SetTitle("#eta_{HJet}");
    hJetEta->GetYaxis()->SetTitle("Candidates/0.1");
    hJetPhi        = new TH1D("hJetPhi", "jet Phi", 34, -3.4, 3.4);
    hJetPhi->GetXaxis()->SetTitle("#varphi_{HJet} [rad]");
    hJetPhi->GetYaxis()->SetTitle("Candidates/0.2[rad]");
    hNJets         = new TH1D("hNJets", "jet multiplicity", 40, 0., 40.);
    hNJets->GetXaxis()->SetTitle("N_{jets}");
    hNJets->GetYaxis()->SetTitle("Events");
    h2JetPt        = new TH1D("h2JetPt", "2jet Pt",    160, 0., 800.);
    h2JetPt->GetXaxis()->SetTitle("p_{T jj}[GeV/c]");
    h2JetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    h2JetM         = new TH1D("h2JetM", "M jj", 4000, 0., 2000.);
    h2JetM->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
    h2JetM->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hAngleProd = new TH1D("hAngleProd", "M ", 32, 0., 3.2);
    hAngleProd ->GetXaxis()->SetTitle("#alpha_{#tau^{+}#tau^{-} - #muZ}[rad]");
    hAngleProd ->GetYaxis()->SetTitle("Candidates/0.1[rad]");
    hCosAngleProd = new TH1D("hCosAngleProd", "M ", 50, -1., 1.);
    hCosAngleProd ->GetXaxis()->SetTitle("cos#alpha_{#tau^{+}#tau^{-} - #muZ}[rad]");
    hCosAngleProd ->GetYaxis()->SetTitle("Candidates/0.04");



    // CATEGORIES
    hInvMass2Mu_A_cat1    = new TH1D("hInvMass2Mu_A_cat1", "M mumu A cat 1", 4000, 0., 2000.);
    hInvMass2Mu_A_cat1->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_A_cat1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_A_cat2    = new TH1D("hInvMass2Mu_A_cat2", "M mumu A cat 2", 4000, 0., 2000.);
    hInvMass2Mu_A_cat2->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_A_cat2->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_A_cat3    = new TH1D("hInvMass2Mu_A_cat3", "M mumu A cat 3", 4000, 0., 2000.);
    hInvMass2Mu_A_cat3->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_A_cat3->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    
    hInvMass2Mu_B_cat1    = new TH1D("hInvMass2Mu_B_cat1", "M mumu B cat 1", 4000, 0., 2000.);
    hInvMass2Mu_B_cat1->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GV/c^{2}]");
    hInvMass2Mu_B_cat1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_B_cat2    = new TH1D("hInvMass2Mu_B_cat2", "M mumu B cat 2", 4000, 0., 2000.);
    hInvMass2Mu_B_cat2->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_B_cat2->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_B_cat3    = new TH1D("hInvMass2Mu_B_cat3", "M mumu B cat 3", 4000, 0., 2000.);
    hInvMass2Mu_B_cat3->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_B_cat3->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    
    hInvMass2Mu_C_cat1    = new TH1D("hInvMass2Mu_C_cat1", "M mumu C cat 1", 4000, 0., 2000.);
    hInvMass2Mu_C_cat1->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_C_cat1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_C_cat2    = new TH1D("hInvMass2Mu_C_cat2", "M mumu C cat 2", 4000, 0., 2000.);
    hInvMass2Mu_C_cat2->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_C_cat2->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");
    hInvMass2Mu_C_cat3    = new TH1D("hInvMass2Mu_C_cat3", "M mumu C cat 3", 4000, 0., 2000.);
    hInvMass2Mu_C_cat3->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_C_cat3->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");


    // EFFICIENCIES
    eff_names.push_back("Skim");                                  eff_counters.push_back(&nEv_Skim);
    eff_names.push_back("Gen Z->2#Mu");                           eff_counters.push_back(&nEv_GenZ2Mu);
    //eff_names.push_back("HLT_IsoMu24/30*");                       eff_counters.push_back(&nEv_TriggerSMu);
    eff_names.push_back("Trigger");                               eff_counters.push_back(&nEv_TriggerSMu);
    eff_names.push_back("Ndf_{PV}>4");                            eff_counters.push_back(&nEv_PVNdf);
    eff_names.push_back("|z_{PV}|<24cm");                         eff_counters.push_back(&nEv_PVZ);
    eff_names.push_back("PV");                                    eff_counters.push_back(&nEv_PV);
    eff_names.push_back("G & Tr.");                               eff_counters.push_back(&nEv_GAndTr);
    eff_names.push_back("p_{T#mu}>10 GeV/c");                     eff_counters.push_back(&nEv_Pt);
    eff_names.push_back("|#eta_{#mu}|< 2.4");                     eff_counters.push_back(&nEv_Eta);
    eff_names.push_back("p_{T#mu}>24 GeV/c & |#eta_{#mu}|< 2.1"); eff_counters.push_back(&nEv_PtEtaMax);
    eff_names.push_back("Muon ID");                               eff_counters.push_back(&nEv_MuID);
//    eff_names.push_back("#chi^{2}/ndf < 10");                     eff_counters.push_back(&nEv_Chi2Ndf);
//    eff_names.push_back("N_{Muon Hits}>0");                       eff_counters.push_back(&nEv_MuonHits);
//    eff_names.push_back("N_{Matched Stations}>1");                eff_counters.push_back(&nEv_MatchedStations);
//    eff_names.push_back("|d_{xy}|<0.02cm");                       eff_counters.push_back(&nEv_Dxy);
//    eff_names.push_back("|d_{z}|<0.1cm");                         eff_counters.push_back(&nEv_Dz);
//    eff_names.push_back("N_{Pixel Hits}>0");                      eff_counters.push_back(&nEv_PixelHits);
//    eff_names.push_back("N_{Tracker Layers}>5");                  eff_counters.push_back(&nEv_TrackerLayers);
    eff_names.push_back("Iso^{PU}_{#mu} < 0.12");                 eff_counters.push_back(&nEv_IsoPU);
    eff_names.push_back("2#mu");                                  eff_counters.push_back(&nEv_2SGMu);
    eff_names.push_back("1 tight #mu");                           eff_counters.push_back(&nEv_1TMu);
    eff_names.push_back("2 tight #mu");                           eff_counters.push_back(&nEv_1TMu);
    eff_names.push_back("4#mu");                                  eff_counters.push_back(&nEv_4SGMu);
    eff_names.push_back("3 tight #mu");                           eff_counters.push_back(&nEv_3TMu);
    eff_names.push_back("4 tight #mu");                           eff_counters.push_back(&nEv_4TMu);
//    eff_names.push_back("3^{rd}#mu-veto");                        eff_counters.push_back(&nEv_3rdMuVeto);
//    eff_names.push_back("GSF 90%");                               eff_counters.push_back(&nEv_GSF);
//    eff_names.push_back("p_{Te}>20 GeV/c");                       eff_counters.push_back(&nEv_PtE);
//    eff_names.push_back("|#eta_{e}|< 2.5");                       eff_counters.push_back(&nEv_EtaE);
//    eff_names.push_back("#chi^{2}/ndf_{e} < 10");                 eff_counters.push_back(&nEv_Chi2NdfE);
//    eff_names.push_back("|d_{xy e}|<0.02cm");                     eff_counters.push_back(&nEv_DxyE);
//    eff_names.push_back("|d_{z e}|<0.1cm");                       eff_counters.push_back(&nEv_DzE);
//    eff_names.push_back("N_{Pixel Hits e}>0");                    eff_counters.push_back(&nEv_PixelHitsE);
//    eff_names.push_back("N_{Tracker Layers e}>5");                eff_counters.push_back(&nEv_TrackerLayersE);
//    eff_names.push_back("Iso^{PU}_{e} < 0.12");                   eff_counters.push_back(&nEv_IsoPUE);
//    eff_names.push_back("2e");                                    eff_counters.push_back(&nEv_2SGE);
//    eff_names.push_back("3^{rd}e-veto");                          eff_counters.push_back(&nEv_3rdEVeto);
//    eff_names.push_back("e: |dz| < 0.14cm");                      eff_counters.push_back(&nEv_SamePVE);
//    eff_names.push_back("e: Q_{#mu1}*Q_{#mu2}<0");                eff_counters.push_back(&nEv_ChargeE);
//    eff_names.push_back("e: M_{e^{+}e^{-}>60GeV/c^{s}}");         eff_counters.push_back(&nEv_InvMassE);
    eff_names.push_back("#mu: |dz| < 0.14cm");                    eff_counters.push_back(&nEv_SamePVMu);
    eff_names.push_back("#mu: Q_{#mu1}*Q_{#mu2}<0");              eff_counters.push_back(&nEv_ChargeMu);
    eff_names.push_back("#mu: M_{#mu^{+}#mu^{-}>60GeV/c^{s}}");   eff_counters.push_back(&nEv_InvMassMu);
    eff_names.push_back("#eta gap jet veto");                     eff_counters.push_back(&nEv_etaGapJetVeto);


    // set up efficiencies histogram
    hEfficiencies  = new TH1D("hEfficiencies", "Efficiencies", eff_counters.size(), 0., (double)eff_counters.size());
    hEfficiencies->GetYaxis()->SetTitle("Events");
    for(unsigned int i = 0; i < eff_counters.size(); i++) {
        *eff_counters[i] = 0;
        hEfficiencies->GetXaxis()->SetBinLabel((i+1) , eff_names[i]);
    }

    // debug = true;
    debug = false;
}
// Destructor:
MyAnalysis::~MyAnalysis() {
    // set up efficiencies histogram
    for(unsigned int i = 0; i < eff_counters.size(); i++) {
        hEfficiencies->SetBinContent((i+1), *eff_counters[i]);
    }
    histfile->Write();
    histfile->Close();
    //print out
    cout<<"====================================================================================="<<endl;
    cout<<setw(40)<<"Selection"<<setw(15)<<setprecision(5)<<"Events"<<setw(15)<<"Eff.(Z2Mu) [%]"<<setw(15)<<"Rel.Eff. [%]"<<endl;
    cout<<"-------------------------------------------------------------------------------------"<<endl;
    cout<<setw(40)<<eff_names[0]<<setw(15)<<(*eff_counters[0])<<setw(15)<<" --- "<<setw(15)<<" --- "<<endl;
    for(unsigned int i = 1; i < eff_names.size(); i++) {
        cout<<setw(40)<<eff_names[i]<<setw(15)<<(*eff_counters[i])<<setw(15)<<100.*(*eff_counters[i])/nEv_GenZ2Mu<<setw(15)<<100.*(*eff_counters[i])/(*eff_counters[i-1])<<endl;
    }
    cout<<"====================================================================================="<<endl;
}

// Analysis
Int_t MyAnalysis::AnalyseEvent() {

    //pileup
    double pileupweight      = GetPileUpWeight(dataPileUp); //GetPrimVertexWeight(dataPrimVert);
    double pileupweight_up   = GetPileUpWeight(dataPileUp_Up);
    double pileupweight_down = GetPileUpWeight(dataPileUp_Down);
    //if( Run() == 191057) cout<<"c0"<<" Run:"<<Run()<<" Event:"<<Number()<<endl;
    ++nEv_Skim;



// require gen Z that goes to 2 muons


    for(unsigned int i = 0; i < NumAllGenParticles(); ++i) {
        //if (AllGenParticles(i).PDGId() == 25 || AllGenParticles(i).PDGId() == 23 || AllGenParticles(i).PDGId() == 24) {
            cout<<AllGenParticles(i).PDGId()<<endl;
          //  cout<<AllGenParticles(i).NumDaughters()<<endl;
        //}
    }
    cout<<endl;        



//++nEv_GenZ2Mu





    // get the trigger
    //////////////////////////////////////////////////////////////////////
    vector<string> triggernames;
    bool useHLTIsoMu24 = false;
    bool useHLTMu17Mu8= false;
    bool useHLTMu17TkMu8 = false;

    // Turn on triggers

    // Single muon triggers
    //useHLTIsoMu24 = true;

    // Double muon triggers
    //useHLTMu17Mu8 = true;
    //useHLTMu17TkMu8 = true;
    
    //HLT_IsoMu20_eta2p1_v8
    //HLT_IsoMu30_v12
    //HLT_IsoMu40_eta2p1_v11

    triggernames.push_back("HLT_IsoMu24_IterTrk02_v1");

    //triggernames.push_back("HLT_DoubleMu5_IsoMu5");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v1");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v2");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v3");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v4");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v5");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v6");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v7");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v8");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v9");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v10");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v11");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v12");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v13");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v14");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v15");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v16");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v17");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v18");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v19");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v20");
    //triggernames.push_back("HLT_DoubleMu5_IsoMu5_v21");
    //trigNameString = "HLT_DoubleMu5_IsoMu5";

    //triggernames.push_back("HLT_IsoMu24");
    //triggernames.push_back("HLT_IsoMu24_v1");
    //triggernames.push_back("HLT_IsoMu24_v2");
    //triggernames.push_back("HLT_IsoMu24_v3");
    //triggernames.push_back("HLT_IsoMu24_v4");
    //triggernames.push_back("HLT_IsoMu24_v5");
    //triggernames.push_back("HLT_IsoMu24_v6");
    //triggernames.push_back("HLT_IsoMu24_v7");
    //triggernames.push_back("HLT_IsoMu24_v8");
    //triggernames.push_back("HLT_IsoMu24_v9");
    //triggernames.push_back("HLT_IsoMu24_v10");
    //triggernames.push_back("HLT_IsoMu24_v11");
    //triggernames.push_back("HLT_IsoMu24_v12");
    //triggernames.push_back("HLT_IsoMu24_v13");
    //triggernames.push_back("HLT_IsoMu24_v14");
    //triggernames.push_back("HLT_IsoMu24_v15");
    //triggernames.push_back("HLT_IsoMu24_v16");
    //triggernames.push_back("HLT_IsoMu24_v17");
    //triggernames.push_back("HLT_IsoMu24_v18");
    //trigNameString = "HLT_IsoMu24";

    if(useHLTIsoMu24) {
        // PromptReco-v1 Run2012
        triggernames.push_back("HLT_IsoMu24_eta2p1_v23");
    }

    if(useHLTMu17Mu8) {
        triggernames.push_back("HLT_Mu17_Mu8_v18");
    }

    if(useHLTMu17TkMu8) {
        triggernames.push_back("HLT_Mu17_TkMu8_v11");
    }
    Int_t trigger = GetHLTrigger(triggernames);
    // check trigger
    if(trigger == -1) {
        cout << "Trigger Problem" << endl;
        for(int i = 0; i < GetNumHLTriggers(); ++i) {
            if(debug) {
                cout<<"GetHLTriggerName(i):"<<GetHLTriggerName(i)<<endl;
            }
        }
//        return(1);
    } else if(trigger == 0) {
        cout<<"trigger did not fire!"<<endl;
        return(1);
    }
    ++nEv_TriggerSMu;




    //////////////////////////////////////////////////////////////////////
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
    hVtxN_u->Fill(PVs.size());
    if( isData) pileupweight = 1;
    if (!isData) pileupweight=this->GetMCWeight(PVs.size(), "myMC", "myData");

    hPUweight->Fill(pileupweight);

    for(unsigned int v = 0; v < PVs.size(); ++v) {
        hVtxNdf->Fill(PVs[v].Ndof(), pileupweight);
        hVtxZ  ->Fill(PVs[v].Z(), pileupweight);
    }
    hVtxN->Fill(PVs.size(), pileupweight);

    if(isVtxNdfOK) ++nEv_PVNdf;
    if(isVtxZOK) ++nEv_PVZ;
    if(PVs.size() < 1) return (1);
    ++nEv_PV;

//  pileupweight = 1.; //data
    pileupweight = this->GetMCWeight(PVs.size(), "myMC", "myData");  //MC
    //////////////////////////////////////////////////////////////////////
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
    bool isMuIDOk = false;
    bool isIsoPUOK            = false;
    int nLooseMus = 0;
    int nTightMus = 0;
    int nVetoMus = 0;

    for(unsigned int i = 0; i < NumMuons(); ++i) {
cout<<"muon ["<<i<<"] id = "<<Muons(i).MuID()<<endl; cout<<endl; cout<<endl;
}
    // loop over the muons
    for(unsigned int i = 0; i < NumMuons(); ++i) {
        if(!(Muons(i).IsGlobal() && Muons(i).IsTracker())) continue;
        ////if( !Muons(i).isPFMuon())continue;
        isGAndTr = true;
        if(!(Muons(i).Pt() > cPtMu)) continue;
        isPtCutOK = true;
        if(!(TMath::Abs(Muons(i).Eta()) <= cEtaMu)) continue;
        isEtaCutOK = true;
        if((Muons(i).Pt() > cPtMuMax && TMath::Abs(Muons(i).Eta()) <= cEtaMuMax)) ++nMuPtEtaMax;
        //if(!(Muons(i).Chi2OverNdof() < cChi2NdfMu)) continue;
        //isChi2NdfOK = true;
        //if(!(Muons(i).OuterTrack().NHits() > cNHitsMuonMu)) continue;
        //isNMuonHitsOK = true;
        //if(!(Muons(i).NumChambersWithSegments() > cNMatchedStations)) continue;
        //isNMatchedStationsOK= true;
        //if(!(TMath::Abs(Muons(i).InnerTrack().Dxy()) < cDxyMu)) continue;
        //isDxyOK = true;
        //if(!(TMath::Abs(Muons(i).InnerTrack().Dz()) < cDzMu)) continue;
        //isDzOK = true;
        //if(!(Muons(i).InnerTrack().NPixelHits() > cNHitsPixelMu)) continue;
        //isNPixelHitsOK = true;
        //if(!((Muons(i).InnerTrack().NPixelLayers() + Muons(i).InnerTrack().NStripLayers())> cNTrackerLayersMu)) continue;
        //isNTrackerLayersOK = true;
        if (!(Muons(i).MuID() > 1 && Muons(i).MuID() < 5)) {
            nVetoMus++;
            continue;
        }
        isMuIDOk = true;
        if (Muons(i).MuID() == 4) nTightMus++;
        if (Muons(i).MuID() == 2) nLooseMus++;
        if(!(this->getIsoPU(Muons(i), AK4PFRho()) < cIsoMuPU)) continue;
        isIsoPUOK = true;
        // save the good mus
        muonsVec.push_back(Muons(i));
    }

    if(isGAndTr) ++nEv_GAndTr;
    if(isPtCutOK) ++nEv_Pt;
    if(isEtaCutOK) ++nEv_Eta;
    if(nMuPtEtaMax > 0) ++nEv_PtEtaMax;
    if(nMuPtEtaMax < 1) return(1);
//    if(isChi2NdfOK) ++nEv_Chi2Ndf;
//    if(isNMuonHitsOK) ++nEv_MuonHits;
//    if(isNMatchedStationsOK) ++nEv_MatchedStations;
//    if(isDxyOK) ++nEv_Dxy;
//    if(isDzOK) ++nEv_Dz;
//    if(isNPixelHitsOK) ++nEv_PixelHits;
//    if(isNTrackerLayersOK) ++nEv_TrackerLayers;
    if(isMuIDOk) ++nEv_MuID;
    if(isIsoPUOK) ++nEv_IsoPU;

    // select events with at least 2 good mus
    if(muonsVec.size() < 2) return(1);
    ++nEv_2SGMu;
    if (nTightMus >= 1) ++nEv_1TMu; else return(1);
    if (nTightMus >= 2) ++nEv_2TMu; else return(1);
    // select events with at least 4 good mus
    if( muonsVec.size() < 4) return(1);
    ++nEv_4SGMu;
    if (nTightMus >= 3) ++nEv_3TMu; else return(1);
    if (nTightMus >= 4) ++nEv_4TMu; else return(1);

    //get the Z->mumu
    bool isChargeMuCutOK  = false;
    bool isSamePVMuCutOK  = false;
    bool isInvMassMuCutOK = false;
    std::vector< TLorentzVector > diMuonVec;
    std::vector< int > Muon_idx1;
    std::vector< int > Muon_idx2;

    // 1st loop over muons
    for(unsigned int i = 0; i < muonsVec.size(); ++i) {

        // 2nd loop over muons
        for(unsigned int j = i+1; j < muonsVec.size(); ++j) {
            // check muon charge
            if(muonsVec[i].Charge()*muonsVec[j].Charge() > 0) continue;
            isChargeMuCutOK = true;     
            // dz cut
            double dzMu = muonsVec[i].Dz() - muonsVec[j].Dz();
            if(!(TMath::Abs(dzMu) < 0.14)) continue;
            isSamePVMuCutOK = true;
            // check the invariant mass
            TLorentzVector diMuon = muonsVec[i] + muonsVec[j];
            if(!(diMuon.M() > cInvMass)) continue;
            isInvMassMuCutOK = true;

            // All possible opposite-sign combinations
            hInvMass2Mu  ->Fill(diMuon.M(), pileupweight);
            h2MuEta  ->Fill(diMuon.Eta(), pileupweight);

            diMuonVec.push_back(diMuon);
            // make two arrays: if i and j make a pair, put each index in a separate array
            // so Muon_idx1[0] and Muon_idx2[0] make a pair, and Muon_idx1[1] and Muon_idx2[1] make another, and so on
            Muon_idx1.push_back(i);
            Muon_idx2.push_back(j);
        } // end 2nd loop over muons

    } //end 1st loop over muons

    if(!isSamePVMuCutOK) return(1);
    ++nEv_SamePVMu;
    if(!isChargeMuCutOK) return(1);
    ++nEv_ChargeMu;
    if(!isInvMassMuCutOK) return(1);
    ++nEv_InvMassMu;

    // get the Z dimuon pair
    int Zmu1_idx    = -1;
    int Zmu2_idx    = -1;
    int diMuonZ_idx = -1;
    double dMassMax = 1000.; //GeV // i had to increase this from 100 to prevent segfaults when Z index = -1
    for(int d = 0; d < diMuonVec.size(); ++d) { //find the pair whose diMuon.M is closest to Mz
        if(TMath::Abs(diMuonVec[d].M() - massZ) < dMassMax) {   //find out if this pair is better than the last one
            dMassMax    = TMath::Abs(diMuonVec[d].M() - massZ);
            Zmu1_idx    = Muon_idx1[d];
            Zmu2_idx    = Muon_idx2[d];
            diMuonZ_idx = d;
        }
    }

    // set the other one to the H pair
    int diMuonH_idx = -1;
    int Hmu1_idx    = -1;
    int Hmu2_idx    = -1;

    for(int d = 0; d < diMuonVec.size(); ++d) {
        if(d == diMuonZ_idx) continue;
        if(Muon_idx1[d] == Zmu1_idx) continue;
        if(Muon_idx2[d] == Zmu2_idx) continue;
        diMuonH_idx = d;
        Hmu1_idx = Muon_idx1[d];
        Hmu2_idx = Muon_idx2[d];
    }

    if(diMuonH_idx == -1 || diMuonZ_idx == -1) return(1);

    hInvMass2MuZ  ->Fill(diMuonVec[diMuonZ_idx].M(), pileupweight);
    hInvMass2Mu_wo->Fill(diMuonVec[diMuonZ_idx].M(), pileupweight);
    hMuDPtZ ->Fill(muonsVec[Zmu1_idx].Pt() - muonsVec[Zmu2_idx].Pt(),  pileupweight);
    hMuDEtaZ->Fill(muonsVec[Zmu1_idx].Eta()- muonsVec[Zmu2_idx].Eta(), pileupweight);
    hMuDPhiZ->Fill(muonsVec[Zmu1_idx].DeltaPhi(muonsVec[Zmu2_idx]),  pileupweight);
    h2MuPtZ  ->Fill(diMuonVec[diMuonZ_idx].Pt(),  pileupweight);
    h2MuEtaZ ->Fill(diMuonVec[diMuonZ_idx].Eta(), pileupweight);
    h2MuPhiZ ->Fill(diMuonVec[diMuonZ_idx].Phi(), pileupweight);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // JETS ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // for each jet, check against all of the mu to make sure it is far enough away
    std::vector< Jet > jetsVec;
    for(unsigned int jet = 0; jet < NumAK4PFJets(); ++jet) {
        bool isdROK = true;
        for(unsigned int i = 0; i < muonsVec.size(); ++i) {
            double dR = AK4PFJets(jet).DeltaR(muonsVec[i]);
            if(dR < cDR) isdROK = false;
        }
        if(!isdROK) continue;
        if(!(this->hasJetIDLoose(AK4PFJets(jet)))) continue;
        if(!(AK4PFJets(jet).Pt() > cPtJet)) continue;
        if(!(TMath::Abs(AK4PFJets(jet).Eta()) <= cEtaJet)) continue;
        // store the jets
        jetsVec.push_back(AK4PFJets(jet));
    }

    // eta gap jet veto
    // if there are any jets inside the eta window, skip the event
    double dMinEtaGap = TMath::Min(diMuonVec[diMuonZ_idx].Eta(), diMuonVec[diMuonH_idx].Eta());
    double dMaxEtaGap = TMath::Max(diMuonVec[diMuonZ_idx].Eta(), diMuonVec[diMuonH_idx].Eta());
    int nJetsEtaGap   = 0;
    for(unsigned int jet = 0; jet < jetsVec.size(); ++jet) {
        hJetPt ->Fill(jetsVec[jet].Pt() , weight);
        hJetEta->Fill(jetsVec[jet].Eta(), weight);
        hJetPhi->Fill(jetsVec[jet].Phi(), weight);
        if(jetsVec[jet].Eta() > dMinEtaGap && jetsVec[jet].Eta() < dMaxEtaGap)
            ++nJetsEtaGap;
    }
    hNJets->Fill(jetsVec.size(), weight);

    if(nJetsEtaGap > 0) return(1);
    ++nEv_etaGapJetVeto;
    hRho->Fill(AK4PFRho());


    ////////////////////////////////////////////////////////////////////////////////
    // muon plots
    //cout<<"new ev."<<endl;
    for(unsigned int i = 0; i < muonsVec.size(); ++i) {
        hMuPt       ->Fill(muonsVec[i].Pt(),  pileupweight);
        hMuEta      ->Fill(muonsVec[i].Eta(), pileupweight);
        hMuPhi      ->Fill(muonsVec[i].Phi(), pileupweight);
        hMuIso      ->Fill(this->getIsoPU(muonsVec[i], AK4PFRho()),  pileupweight);
        hMuIsoPU    ->Fill(muonsVec[i].IsoR3ECal() + muonsVec[i].IsoR3HCal() - AK4PFRho()*Pi()*0.09);
        hMuChi2Ndf  ->Fill(muonsVec[i].Chi2OverNdof(),  pileupweight);
        hMuChi2NdfIT->Fill(muonsVec[i].InnerTrack().Chi2OverNdof(),  pileupweight);
        hMuITHits   ->Fill(muonsVec[i].InnerTrack().NHits(),  pileupweight);
        hMuOTHits   ->Fill(muonsVec[i].OuterTrack().NHits(),  pileupweight);
        hMuPHits    ->Fill(muonsVec[i].InnerTrack().NPixelHits(),  pileupweight);
        hMuDxy      ->Fill(muonsVec[i].InnerTrack().Dxy(),  pileupweight);
        hMuDz       ->Fill(muonsVec[i].InnerTrack().Dz(),  pileupweight);
    }


    if(diMuonH_idx == -1) return(1);


    // how many muons are in the barrel?
    int barrelcntr = 0;
    for (unsigned int muidx = 0; muidx < muonsVec.size(); muidx++) {
        if (muonsVec.size() != 4 ) continue;
        if (TMath::Abs(muonsVec[muidx].Eta()) <= 0.8 ) barrelcntr++;
    }

    // Fill categorisation A plots
    if (barrelcntr == 4) hInvMass2Mu_A_cat1->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else if (barrelcntr == 2 || barrelcntr == 3) hInvMass2Mu_A_cat2->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else if (barrelcntr == 0 || barrelcntr == 1) hInvMass2Mu_A_cat3->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else cout<<"something went wrong with cat A logic"<<endl;

    // Fill categorisation B plots
    if (barrelcntr == 4 || barrelcntr == 3) hInvMass2Mu_B_cat1->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else if (barrelcntr == 2 || barrelcntr == 1) hInvMass2Mu_B_cat2->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else if (barrelcntr == 0) hInvMass2Mu_B_cat3->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else cout<<"something went wrong with cat B logic"<<endl;

    // Fill categorisation C plots
    double etaH1 = TMath::Abs(muonsVec[Hmu1_idx].Eta());
    double etaH2 = TMath::Abs(muonsVec[Hmu2_idx].Eta());
    if (etaH1 <= 0.8 && etaH2 <= 0.8) hInvMass2Mu_C_cat1->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else if (etaH1 <= 0.8 || etaH2 <= 0.8) hInvMass2Mu_C_cat2->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else if (!(etaH1 <= 0.8 || etaH2 <= 0.8)) hInvMass2Mu_C_cat3->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    else cout<<"something went wrong with cat C logic"<<endl;


    if(isData && !(diMuonVec[diMuonH_idx].M() >= cInvMassB_low &&  diMuonVec[diMuonH_idx].M() <= cInvMassB_high)) {
        hInvMass2MuH  ->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
        hInvMass2Mu_wo->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    }

    if(!isData) {
        hInvMass2MuH  ->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
        hInvMass2Mu_wo->Fill(diMuonVec[diMuonH_idx].M(), pileupweight);
    }

    hMuDPtH  ->Fill(muonsVec[Hmu1_idx].Pt() - muonsVec[Hmu2_idx].Pt(),  pileupweight);
    hMuDEtaH ->Fill(muonsVec[Hmu1_idx].Eta()- muonsVec[Hmu2_idx].Eta(), pileupweight);
    hMuDPhiH ->Fill(muonsVec[Hmu1_idx].DeltaPhi(muonsVec[Hmu2_idx]),  pileupweight);
    h2MuPtH  ->Fill(diMuonVec[diMuonH_idx].Pt(),  pileupweight);
    h2MuEtaH ->Fill(diMuonVec[diMuonH_idx].Eta(), pileupweight);
    h2MuPhiH ->Fill(diMuonVec[diMuonH_idx].Phi(), pileupweight);

    TLorentzVector fourMuons = (muonsVec[Zmu1_idx] + muonsVec[Zmu2_idx] +
                                muonsVec[Hmu1_idx] + muonsVec[Hmu2_idx]);
    hInvMass4Mu  ->Fill(fourMuons.M(), pileupweight);
    hDPhi2MuZ2MuH->Fill(diMuonVec[diMuonZ_idx].DeltaPhi(diMuonVec[diMuonH_idx]), pileupweight);

    // fill the tree for limits calculations
    trunNr   = Run();
    tevNr    = Number();
    tInvMass = diMuonVec[0].M();
    ftreeCat1->Fill();

    cout<<endl;
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
            if(filename.find("LUMI_") == 0) {
                ana.AddLumiFile(namebuf.c_str());
                //}else if(filename.find("ReWeight1D") == 0){
            } else if(filename.find("PileUp69400_22jan_2012") == 0) {
                ana.AddPUWeightFile(namebuf.c_str());
            } else {
                ana.AddFile(namebuf.c_str());
            }
        }
    }
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
void MyAnalysis::AddPUWeightFile(string filename) {
    TFile weightfile(filename.c_str(), "READ");
    if(weightfile.IsZombie()) {
        cerr << "ERROR AddLumiFile: " << filename << " is not a valid file."<<endl;
    }

    TH1D *Vertex      = (TH1D *)weightfile.Get("pileup69400");
    TH1D *Vertex_Up   = (TH1D *)weightfile.Get("pileup72870");
    TH1D *Vertex_Down = (TH1D *)weightfile.Get("pileup65930");

    Double_t tot      = Vertex->Integral();
    Double_t tot_up   = Vertex_Up->Integral();
    Double_t tot_down = Vertex_Down->Integral();

    for(int i=0; i<Vertex->GetNbinsX(); i++) {
        dataPileUp.push_back(Vertex->GetBinContent(i+1)/tot);
        // cout<<dataPileUp[i]<<endl;
    }

    for(int i=0; i<Vertex_Up->GetNbinsX(); i++) {
        dataPileUp_Up.push_back(Vertex_Up->GetBinContent(i+1)/tot_up);
        // cout<<dataPileUp_Up[i]<<endl;
    }

    for(int i=0; i<Vertex_Down->GetNbinsX(); i++) {
        dataPileUp_Down.push_back(Vertex_Down->GetBinContent(i+1)/tot_down);
        // cout<<dataPileUp_Down[i]<<endl;
    }

    weightfile.Close();
}
//________________________________
double MyAnalysis::GetMCWeight(int nPVMC, string nameMC, string nameData) {
    double my_weight = 1.;
    if(nameMC == "DYJetsToLL" && nameData == "Run2012") return w2012_DY[nPVMC];
    if(nameMC == "TTJets" && nameData  == "Run2012") return w2012_TT[nPVMC];
    if(nameMC == "QCD" && nameData == "Run2012") return w2012_QCD[nPVMC];
    if(nameMC == "Tbar_s" && nameData  == "Run2012") return w2012_Tbar_s[nPVMC];
    if(nameMC == "Tbar_t" && nameData  == "Run2012") return w2012_Tbar_t[nPVMC];
    if(nameMC == "Tbar_tW" && nameData == "Run2012") return w2012_Tbar_tW[nPVMC];
    if(nameMC == "T_s" && nameData  == "Run2012") return w2012_T_s[nPVMC];
    if(nameMC == "T_t" && nameData  == "Run2012") return w2012_T_t[nPVMC];
    if(nameMC == "T_tW" && nameData == "Run2012") return w2012_T_tW[nPVMC];
    if(nameMC == "WJetsToLNu" && nameData == "Run2012") return w2012_WJets[nPVMC];
    if(nameMC == "WW" && nameData == "Run2012") return w2012_WW[nPVMC];
    if(nameMC == "WZ" && nameData == "Run2012") return w2012_WZ[nPVMC];
    if(nameMC == "ZZ" && nameData == "Run2012") return w2012_ZZ[nPVMC];
    if(nameMC == "H" && nameData == "Run2012") return w2012_H[nPVMC];
    return my_weight;
}
//________________________________
void MyAnalysis::AddWeightFile(string filename) {
    TFile *weightfile = new TFile(filename.c_str());
    if(weightfile->IsZombie()) {
        cerr << "ERROR AddLumiFile: " << filename << " is not a valid file."<<endl;
    }
    TH1D *nPV_Run2012     = (TH1D *) weightfile->Get("nPV_Run2012");
    TH1D *nPV_DYJetsToLL  = (TH1D *) weightfile->Get("nPV_DYJetsToLL");
    TH1D *nPV_TTJets      = (TH1D *) weightfile->Get("nPV_TTJets");
    TH1D *nPV_QCD         = (TH1D *) weightfile->Get("nPV_QCD");
    TH1D *nPV_Tbar_s      = (TH1D *) weightfile->Get("nPV_Tbar_s");
    TH1D *nPV_Tbar_t      = (TH1D *) weightfile->Get("nPV_Tbar_t");
    TH1D *nPV_Tbar_tW     = (TH1D *) weightfile->Get("nPV_Tbar_tW");
    TH1D *nPV_T_s         = (TH1D *) weightfile->Get("nPV_T_s");
    TH1D *nPV_T_t         = (TH1D *) weightfile->Get("nPV_T_t");
    TH1D *nPV_T_tW        = (TH1D *) weightfile->Get("nPV_T_tW");
    TH1D *nPV_WJets       = (TH1D *) weightfile->Get("nPV_WJetsToLNu");
    TH1D *nPV_WW          = (TH1D *) weightfile->Get("nPV_WW");
    TH1D *nPV_WZ          = (TH1D *) weightfile->Get("nPV_WZ");
    TH1D *nPV_ZZ          = (TH1D *) weightfile->Get("nPV_ZZ");
    TH1D *nPV_H           = (TH1D *) weightfile->Get("nPV_Higgs");

    nPV_Run2012   ->Scale(1./(nPV_Run2012   ->Integral()));
    nPV_DYJetsToLL->Scale(1./(nPV_DYJetsToLL->Integral()));
    nPV_TTJets    ->Scale(1./(nPV_TTJets    ->Integral()));
    nPV_QCD       ->Scale(1./(nPV_QCD       ->Integral()));
    nPV_Tbar_s    ->Scale(1./(nPV_Tbar_s    ->Integral()));
    nPV_Tbar_t    ->Scale(1./(nPV_Tbar_t    ->Integral()));
    nPV_Tbar_tW   ->Scale(1./(nPV_Tbar_tW   ->Integral()));
    nPV_T_s       ->Scale(1./(nPV_T_s       ->Integral()));
    nPV_T_t       ->Scale(1./(nPV_T_t       ->Integral()));
    nPV_T_tW      ->Scale(1./(nPV_T_tW      ->Integral()));
    nPV_WJets     ->Scale(1./(nPV_WJets     ->Integral()));
    nPV_WW        ->Scale(1./(nPV_WW        ->Integral()));
    nPV_WZ        ->Scale(1./(nPV_WZ        ->Integral()));
    nPV_ZZ        ->Scale(1./(nPV_ZZ        ->Integral()));
    nPV_H         ->Scale(1./(nPV_H         ->Integral()));


    int nBins = nPV_Run2012->GetXaxis()->GetNbins();
    double sum2012  = 0.;

    // loop over number of PV
    for(int ibin = 0; ibin < nBins; ++ibin) {
        double nPV2012       = nPV_Run2012   ->GetBinContent(nPV_Run2012   ->GetXaxis()->FindBin(ibin));
        double nPVDYJetsToLL = nPV_DYJetsToLL->GetBinContent(nPV_DYJetsToLL->GetXaxis()->FindBin(ibin));
        double nPVTTJets     = nPV_TTJets    ->GetBinContent(nPV_TTJets    ->GetXaxis()->FindBin(ibin));
        double nPVQCD        = nPV_QCD       ->GetBinContent(nPV_QCD       ->GetXaxis()->FindBin(ibin));
        double nPVTbar_s     = nPV_Tbar_s    ->GetBinContent(nPV_Tbar_s    ->GetXaxis()->FindBin(ibin));
        double nPVTbar_t     = nPV_Tbar_t    ->GetBinContent(nPV_Tbar_t    ->GetXaxis()->FindBin(ibin));
        double nPVTbar_tW    = nPV_Tbar_tW   ->GetBinContent(nPV_Tbar_tW   ->GetXaxis()->FindBin(ibin));
        double nPVT_s        = nPV_Tbar_s    ->GetBinContent(nPV_T_s       ->GetXaxis()->FindBin(ibin));
        double nPVT_t        = nPV_Tbar_t    ->GetBinContent(nPV_T_t       ->GetXaxis()->FindBin(ibin));
        double nPVT_tW       = nPV_Tbar_tW   ->GetBinContent(nPV_T_tW      ->GetXaxis()->FindBin(ibin));
        double nPVWJets      = nPV_WJets     ->GetBinContent(nPV_WJets     ->GetXaxis()->FindBin(ibin));
        double nPVWW         = nPV_WW        ->GetBinContent(nPV_WW        ->GetXaxis()->FindBin(ibin));
        double nPVWZ         = nPV_WZ        ->GetBinContent(nPV_WZ        ->GetXaxis()->FindBin(ibin));
        double nPVZZ         = nPV_ZZ        ->GetBinContent(nPV_ZZ        ->GetXaxis()->FindBin(ibin));
        double nPVH          = nPV_H         ->GetBinContent(nPV_H         ->GetXaxis()->FindBin(ibin));

        if(nPVDYJetsToLL != 0) {
            w2012_DY.push_back(nPV2012/nPVDYJetsToLL);
        } else {
            w2012_DY.push_back(nPVDYJetsToLL);
        }
        if(nPVTTJets !=0) {
            w2012_TT.push_back(nPV2012/nPVTTJets);
        } else {
            w2012_TT.push_back(nPVTTJets);
        }
        if(nPVQCD !=0) {
            w2012_QCD.push_back(nPV2012/nPVQCD);
        } else {
            w2012_QCD.push_back(nPVQCD);
        }
        if(nPVTbar_s !=0) {
            w2012_Tbar_s.push_back(nPV2012/nPVTbar_s);
        } else {
            w2012_Tbar_s.push_back(nPVTbar_s);
        }
        if(nPVTbar_t !=0) {
            w2012_Tbar_t.push_back(nPV2012/nPVTbar_t);
        } else {
            w2012_Tbar_t.push_back(nPVTbar_t);
        }
        if(nPVTbar_tW !=0) {
            w2012_Tbar_tW.push_back(nPV2012/nPVTbar_tW);
        } else {
            w2012_Tbar_tW.push_back(nPVTbar_tW);
        }
        if(nPVT_s !=0) {
            w2012_T_s.push_back(nPV2012/nPVT_s);
        } else {
            w2012_T_s.push_back(nPVT_s);
        }
        if(nPVT_t !=0) {
            w2012_T_t.push_back(nPV2012/nPVT_t);
        } else {
            w2012_T_t.push_back(nPVT_t);
        }
        if(nPVT_tW !=0) {
            w2012_T_tW.push_back(nPV2012/nPVT_tW);
        } else {
            w2012_T_tW.push_back(nPVT_tW);
        }
        if(nPVWJets !=0) {
            w2012_WJets.push_back(nPV2012/nPVWJets);
        } else {
            w2012_WJets.push_back(nPVWJets);
        }
        if(nPVWW !=0) {
            w2012_WW.push_back(nPV2012/nPVWW);
        } else {
            w2012_WW.push_back(nPVWW);
        }
        if(nPVWZ !=0) {
            w2012_WZ.push_back(nPV2012/nPVWZ);
        } else {
            w2012_WZ.push_back(nPVWZ);
        }
        if(nPVZZ !=0) {
            w2012_ZZ.push_back(nPV2012/nPVZZ);
        } else {
            w2012_ZZ.push_back(nPVZZ);
        }
        if(nPVH !=0) {
            w2012_H.push_back(nPV2012/nPVH);
        } else {
            w2012_H.push_back(nPVH);
        }
        sum2012  += nPV2012;
    } // end loop over number of PV
    // loop over number of PV
    for(int ibin = 0; ibin < nBins; ++ibin) {
        w2012_DY[ibin] /= sum2012;
        w2012_TT[ibin] /= sum2012;
        w2012_Tbar_s[ibin] /= sum2012;
        w2012_Tbar_t[ibin] /= sum2012;
        w2012_Tbar_tW[ibin] /= sum2012;
        w2012_T_s[ibin] /= sum2012;
        w2012_T_t[ibin] /= sum2012;
        w2012_T_tW[ibin] /= sum2012;
        w2012_QCD[ibin] /= sum2012;
        w2012_WJets[ibin] /= sum2012;
        w2012_WW[ibin] /= sum2012;
        w2012_WZ[ibin] /= sum2012;
        w2012_ZZ[ibin] /= sum2012;
        w2012_H[ibin] /= sum2012;
    } // end of loop over the number of PV
    weightfile->Close();
}

