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
    TH1D *hMuPt_noCuts;
    TH1D *hMuPt_basicCuts;
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
    // Baseline histos
    TH1D *hInvMass2Mu; // all combos
    TH1D *hInvMass2Mu_wo; // pick two

    TH1D *hInvMass2Mu_none;
    TH1D *hInvMass2Mu_one; 
    TH1D *hInvMass2Mu_both;

    TH1D *h2MuEta;
    TH1D *h2MuPt;
    TH1D *hInvMass4Mu;

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
    double nEv_TriggerMu;
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
    double nEv_2SGMu;
    double nEv_1TMu;
    double nEv_2TMu;
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
    string asString(double f);

};
//Constructor:
MyAnalysis::MyAnalysis() : Analyse(), currun(0), curlumi(0) {
    // don;t touch, these are changed with the output file name
    isData=false; isMC=false;

    // vertex
    cVtxNdf = 4;
    cVtxZ   = 24.;

    cPtMu        = 10.0; // GeV;
    cPtMuMax     = 24.0; // Run2012B pT>30GeV; Run2012 pT>24GeV;
    cEtaMu       = 2.4;  // bbZ: 2.1
    cEtaMuMax    = 2.1;  //

    cChi2NdfMu   = 10.0;
    cNHitsMuonMu = 0;
    cNMatchedStations = 1;
    cDxyMu            = 0.02; // cm
    cDzMu             = 0.14; // cm
    cNHitsPixelMu     = 0;
    cNTrackerLayersMu = 5;
    cIsoMuPU          = 0.15;//0.10; // corrected for pile-up

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

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IS IT DATA OR MONTE CARLO? /////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //isData=true;
    isMC=true;
    // output file
    if (isData)    histfile = new TFile("2mu_data_ana_out.root", "RECREATE");
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

    // mu
    hMuPt      = new TH1D("hMuPt", "mu Pt",    160, 0., 800.);
    hMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hMuPt_noCuts      = new TH1D("hMuPt_noCuts", "mu Pt trig only",    160, 0., 800.);
    hMuPt_noCuts->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt_noCuts->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
    hMuPt_basicCuts      = new TH1D("hMuPt_basicCuts", "mu Pt kin. cuts only",    160, 0., 800.);
    hMuPt_basicCuts->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt_basicCuts->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
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

    hInvMass2Mu    = new TH1D("hInvMass2Mu", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hInvMass2Mu_wo = new TH1D("hInvMass2Mu_wo", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu_wo->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_wo->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hInvMass2Mu_none    = new TH1D("hInvMass2Mu_none", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu_none->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_none->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hInvMass2Mu_one    = new TH1D("hInvMass2Mu_one", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu_one->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_one->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    hInvMass2Mu_both    = new TH1D("hInvMass2Mu_both", "M mumu", 4000, 0., 2000.);
    hInvMass2Mu_both->GetXaxis()->SetTitle("M_{#mu^{+}#mu^{-}} [GeV/c^{2}]");
    hInvMass2Mu_both->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

    h2MuEta        = new TH1D("h2MuEta", "2mu Eta",  132, -6.6, 6.6);
    h2MuEta->GetXaxis()->SetTitle("#eta_{#mu^{+}#mu^{-}}");
    h2MuEta->GetYaxis()->SetTitle("Candidates/0.1");

    h2MuPt         = new TH1D("h2MuPt", "2mu pT",  2000, 0, 1000);
    h2MuEta->GetXaxis()->SetTitle("p_{T}_{#mu^{+}#mu^{-}}");
    h2MuEta->GetYaxis()->SetTitle("Candidates/0.5 GeV");

   // // jets
   // hJetPt         = new TH1D("hJetPt", "jet Pt",    160, 0., 800.);
   // hJetPt->GetXaxis()->SetTitle("p_{T HJet}[GeV/c]");
   // hJetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
   // hJetEta        = new TH1D("hJetEta", "jet Eta",  52, -2.6, 2.6);
   // hJetEta->GetXaxis()->SetTitle("#eta_{HJet}");
   // hJetEta->GetYaxis()->SetTitle("Candidates/0.1");
   // hJetPhi        = new TH1D("hJetPhi", "jet Phi", 34, -3.4, 3.4);
   // hJetPhi->GetXaxis()->SetTitle("#varphi_{HJet} [rad]");
   // hJetPhi->GetYaxis()->SetTitle("Candidates/0.2[rad]");
    hNJets         = new TH1D("hNJets", "jet multiplicity", 40, 0., 40.);
    hNJets->GetXaxis()->SetTitle("N_{jets}");
    hNJets->GetYaxis()->SetTitle("Events");
   // h2JetPt        = new TH1D("h2JetPt", "2jet Pt",    160, 0., 800.);
   // h2JetPt->GetXaxis()->SetTitle("p_{T jj}[GeV/c]");
   // h2JetPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");
   // h2JetM         = new TH1D("h2JetM", "M jj", 4000, 0., 2000.);
   // h2JetM->GetXaxis()->SetTitle("M_{jj} [GeV/c^{2}]");
   // h2JetM->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

   // hAngleProd = new TH1D("hAngleProd", "M ", 32, 0., 3.2);
   // hAngleProd ->GetXaxis()->SetTitle("#alpha_{#tau^{+}#tau^{-} - #muZ}[rad]");
   // hAngleProd ->GetYaxis()->SetTitle("Candidates/0.1[rad]");
   // hCosAngleProd = new TH1D("hCosAngleProd", "M ", 50, -1., 1.);
   // hCosAngleProd ->GetXaxis()->SetTitle("cos#alpha_{#tau^{+}#tau^{-} - #muZ}[rad]");
   // hCosAngleProd ->GetYaxis()->SetTitle("Candidates/0.04");



    // EFFICIENCIES
    eff_names.push_back("Skim");                                  eff_counters.push_back(&nEv_Skim);
    eff_names.push_back("Gen Z->2#Mu");                           eff_counters.push_back(&nEv_GenZ2Mu);
    eff_names.push_back("Trigger");                               eff_counters.push_back(&nEv_TriggerMu);
    eff_names.push_back("Ndf_{PV} > "+asString(cVtxNdf));         eff_counters.push_back(&nEv_PVNdf);
    eff_names.push_back("|z_{PV}| < "+asString(cVtxZ)+"cm");      eff_counters.push_back(&nEv_PVZ);
    eff_names.push_back("nPV > 0");                               eff_counters.push_back(&nEv_PV);
    eff_names.push_back("G & Tr.");                               eff_counters.push_back(&nEv_GAndTr);
    eff_names.push_back("p_{T#mu}> "+asString(cPtMu)+" GeV");     eff_counters.push_back(&nEv_Pt);
    eff_names.push_back("|#eta_{#mu}|< "+asString(cEtaMu));       eff_counters.push_back(&nEv_Eta);
    eff_names.push_back("Muon ID");                               eff_counters.push_back(&nEv_MuID);
    eff_names.push_back("Iso #mu < "+asString(cIsoMuPU));         eff_counters.push_back(&nEv_IsoPU);
    eff_names.push_back("p_{T#mu} > "+asString(cPtMuMax) + " GeV & |#eta_{#mu}| < "+asString(cEtaMuMax)); eff_counters.push_back(&nEv_PtEtaMax);
    eff_names.push_back("2#mu");                                  eff_counters.push_back(&nEv_2SGMu);
    eff_names.push_back("1 tight #mu");                           eff_counters.push_back(&nEv_1TMu);
    eff_names.push_back("2 tight #mu");                           eff_counters.push_back(&nEv_2TMu);
    eff_names.push_back("#mu: |dz| < "+asString(cDzMu)+"cm");     eff_counters.push_back(&nEv_SamePVMu);
    eff_names.push_back("#mu: Q_{#mu1}*Q_{#mu2} < 0");            eff_counters.push_back(&nEv_ChargeMu);
    eff_names.push_back("#mu: M_{#mu^{+}#mu^{-} > "+asString(cInvMass)+"GeV");   eff_counters.push_back(&nEv_InvMassMu);
    eff_names.push_back("#eta gap jet veto");                     eff_counters.push_back(&nEv_etaGapJetVeto);


    // set up efficiencies histogram
    hEfficiencies  = new TH1D("hEfficiencies", "Efficiencies", eff_counters.size(), 0., (double)eff_counters.size());
    hEfficiencies->GetYaxis()->SetTitle("Events");
    for(unsigned int i = 0; i < eff_counters.size(); i++) {
        *eff_counters[i] = 0;
        hEfficiencies->GetXaxis()->SetBinLabel((i+1) , eff_names[i]);
    }

    //debug = true;
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
    cout<<setw(40)<<"Selection"<<setw(15)<<setprecision(8)<<"Events"<<setw(15)<<"Eff.(Z2Mu) [%]"<<setw(15)<<"Rel.Eff. [%]"<<endl;
    cout<<"-------------------------------------------------------------------------------------"<<endl;
    cout<<setw(40)<<eff_names[0]<<setw(15)<<(*eff_counters[0])<<setw(15)<<" --- "<<setw(15)<<" --- "<<endl;
    for(unsigned int i = 1; i < eff_names.size(); i++) {
        cout<<setw(40)<<eff_names[i]<<setw(15)<<setprecision(8)<<(*eff_counters[i])<<setw(15)<<setprecision(4)<<100.*(*eff_counters[i])/nEv_GenZ2Mu<<setw(15)<<100.*(*eff_counters[i])/(*eff_counters[i-1])<<endl;
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
    bool genZCheck = false;
    if (isMC) {
        for(unsigned int i = 0; i < NumAllGenParticles(); ++i) {
            if (TMath::Abs(AllGenParticles(i).PDGId()) == 23) {
                if (AllGenParticles(i).NumDaughters() == 2) {
                        for (unsigned int iii = 0; iii < AllGenParticles(i).NumDaughters(); iii++) {
                            if (TMath::Abs(AllGenParticles(i).GetDaughter(iii).PDGId()) == 13 ) genZCheck = true;
                        }
                 }
            }
        }
    } else { genZCheck = true; }

    if (isMC && !genZCheck) return 1;
    ++nEv_GenZ2Mu;

    // get the trigger
    //////////////////////////////////////////////////////////////////////
    vector<string> triggernames;

    // Turn on triggers

    // Single muon triggers
    bool useHLTIsoMu20 = true;
//    bool useHLTIsoMu24 = false;
//
//    // Double muon triggers
//    bool useHLTDoubleIsoMu17 = false;
//    bool useHLTMu27TkMu8 = true;
//    bool useHLTMu17Mu8 = false;

    if (useHLTIsoMu20) {
        triggernames.push_back("HLT_IsoMu20");
        triggernames.push_back("HLT_IsoMu20_v1");
        triggernames.push_back("HLT_IsoMu20_v2");
        triggernames.push_back("HLT_IsoMu20_v3");
        triggernames.push_back("HLT_IsoMu20_v4");
        triggernames.push_back("HLT_IsoMu20_v5");
    }
//    if (useHLTIsoMu24) {
///////        triggernames.push_back("HLT_IsoMu24_eta2p1");
///////        triggernames.push_back("HLT_IsoMu24_eta2p1_v1");
///////        triggernames.push_back("HLT_IsoMu24_eta2p1_v2");
///////        triggernames.push_back("HLT_IsoMu24_eta2p1_v3");
///////        triggernames.push_back("HLT_IsoMu24_eta2p1_v4");
///////        triggernames.push_back("HLT_IsoMu24_eta2p1_v5");
//    }
//    if (useHLTDoubleIsoMu17) {
//        triggernames.push_back("HLT_DoubleIsoMu17_eta2p1");
//        triggernames.push_back("HLT_DoubleIsoMu17_eta2p1_v1");
//        triggernames.push_back("HLT_DoubleIsoMu17_eta2p1_v2");
//        triggernames.push_back("HLT_DoubleIsoMu17_eta2p1_v3");
//        triggernames.push_back("HLT_DoubleIsoMu17_eta2p1_v4");
//    }
//    if (useHLTMu27TkMu8) {
        triggernames.push_back("HLT_Mu27_TkMu8");
        triggernames.push_back("HLT_Mu27_TkMu8_v1");
        triggernames.push_back("HLT_Mu27_TkMu8_v2");
        triggernames.push_back("HLT_Mu27_TkMu8_v3");
//    }
//    if (useHLTMu17Mu8) {
////////        triggernames.push_back("HLT_Mu17_Mu8");
////////        triggernames.push_back("HLT_Mu17_Mu8_v1");
////////        triggernames.push_back("HLT_Mu17_Mu8_v2");
////////        triggernames.push_back("HLT_Mu17_Mu8_v3");
//    }

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
    ++nEv_TriggerMu;




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
    if( isData) pileupweight = 1.;
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
    int nMedMus = 0;
    int nVetoMus = 0;

    for(unsigned int i = 0; i < NumMuons(); ++i) {
        hMuPt_noCuts->Fill(Muons(i).Pt(),  pileupweight);
    }
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
        hMuPt_basicCuts->Fill(Muons(i).Pt(),  pileupweight);
        // collect muon ID
        if (Muons(i).MuID() > 4) {
            cout<<"line something went very wrong with muID here"<<endl;
            nVetoMus++;
            continue;
        } else if (Muons(i).MuID() <= 1) {
            nVetoMus++;
            continue;
        }
        isMuIDOk = true;
        if (Muons(i).MuID() == 4) nTightMus++;
        if (Muons(i).MuID() >= 3) nMedMus++;
        if (Muons(i).MuID() >= 2) nLooseMus++;
        //if (!(this->getIsoPU(Muons(i), AK4PFRho()) < cIsoMuPU)) continue;
        if (!(Muons(i).IsoR3TrackRel() < 0.1)) continue;
        isIsoPUOK = true;
        // at least one muon in the event must ALSO pass this cut:
        if((Muons(i).Pt() > cPtMuMax && TMath::Abs(Muons(i).Eta()) <= cEtaMuMax)) ++nMuPtEtaMax;
        // save the good mus
        muonsVec.push_back(Muons(i));
    }

    if(isGAndTr) ++nEv_GAndTr;
    if(isPtCutOK) ++nEv_Pt;
    if(isEtaCutOK) ++nEv_Eta;
    if(isMuIDOk) ++nEv_MuID;
    if(isIsoPUOK) ++nEv_IsoPU;
    if(nMuPtEtaMax > 0) ++nEv_PtEtaMax;
    if(nMuPtEtaMax < 1) return(1);

    // select events with at least 2 good mus
    if(muonsVec.size() < 2) return(1);
    ++nEv_2SGMu;
    if (nTightMus >= 1) ++nEv_1TMu;
    else return (1);
    if (nTightMus >= 2) ++nEv_2TMu;
    else return (1);

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

            // All possible opposite-sign combinations
            hInvMass2Mu  ->Fill(diMuon.M(), pileupweight);
            h2MuEta  ->Fill(diMuon.Eta(), pileupweight);
            h2MuPt  ->Fill(diMuon.Pt(), pileupweight);

            diMuonVec.push_back(diMuon);
            // make two arrays: if i and j make a pair, put each index in a separate array
            // so Muon_idx1[0] and Muon_idx2[0] make a pair, and Muon_idx1[1] and Muon_idx2[1] make another, and so on
            Muon_idx1.push_back(i);
            Muon_idx2.push_back(j);
        } // end 2nd loop over muons
    } //end 1st loop over muons

    if (!isSamePVMuCutOK) return(1);
    ++nEv_SamePVMu;
    if (!isChargeMuCutOK) return(1);
    ++nEv_ChargeMu;
    if (!isInvMassMuCutOK) return(1);
    ++nEv_InvMassMu;


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // JETS ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // for each jet, check against all of the mu to make sure it is far enough away
    std::vector< Jet > jetsVec;
    for (unsigned int jet = 0; jet < NumAK4PFJets(); ++jet) {
        bool isdROK = true;
        for (unsigned int i = 0; i < muonsVec.size(); ++i) {
            double dR = AK4PFJets(jet).DeltaR(muonsVec[i]);
            if(dR < cDR) isdROK = false;
        }
        if (!isdROK) continue;
        if (!(this->hasJetIDLoose(AK4PFJets(jet)))) continue;
        if (!(AK4PFJets(jet).Pt() > cPtJet)) continue;
        if (!(TMath::Abs(AK4PFJets(jet).Eta()) <= cEtaJet)) continue;
        // store the jets
        jetsVec.push_back(AK4PFJets(jet));
    }

    //// eta gap jet veto
    //// if there are any jets inside the eta window, skip the event
    //double dMinEtaGap = TMath::Min(diMuonVec[diMuonZ_idx].Eta(), diMuonVec[1].Eta());
    //double dMaxEtaGap = TMath::Max(diMuonVec[diMuonZ_idx].Eta(), diMuonVec[diMuonH_idx].Eta());
    //int nJetsEtaGap   = 0;
    //for(unsigned int jet = 0; jet < jetsVec.size(); ++jet) {
    //    hJetPt ->Fill(jetsVec[jet].Pt() , weight);
    //    hJetEta->Fill(jetsVec[jet].Eta(), weight);
    //    hJetPhi->Fill(jetsVec[jet].Phi(), weight);
    //    if(jetsVec[jet].Eta() > dMinEtaGap && jetsVec[jet].Eta() < dMaxEtaGap)
    //        ++nJetsEtaGap;
    //}
    hNJets->Fill(jetsVec.size(), weight);

    //if(nJetsEtaGap > 0) return(1);
    ++nEv_etaGapJetVeto;
    hRho->Fill(AK4PFRho());

    ////////////////////////////////////////////////////////////////////////////////
    // muon plots
    for (unsigned int i = 0; i < muonsVec.size(); ++i) {
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



    // how many muons are in the barrel?
    int barrelcntr = 0;
    for (unsigned int muidx = 0; muidx < muonsVec.size(); muidx++) {
        if (TMath::Abs(muonsVec[muidx].Eta()) <= 0.8 ) barrelcntr++;
    }
    if (barrelcntr == 2) hInvMass2Mu_both->Fill(diMuonVec[0].M(), pileupweight);
    else if (barrelcntr == 1) hInvMass2Mu_one->Fill(diMuonVec[0].M(), pileupweight);
    else if (barrelcntr == 0) hInvMass2Mu_none->Fill(diMuonVec[0].M(), pileupweight);
    // fill the tree for limits calculations
    trunNr   = Run();
    tevNr    = Number();
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
cout<<"AddPUWeightFile"<<endl;
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
        cout<<dataPileUp[i]<<endl;
    }

    for(int i=0; i<Vertex_Up->GetNbinsX(); i++) {
        dataPileUp_Up.push_back(Vertex_Up->GetBinContent(i+1)/tot_up);
        cout<<dataPileUp_Up[i]<<endl;
    }

    for(int i=0; i<Vertex_Down->GetNbinsX(); i++) {
        dataPileUp_Down.push_back(Vertex_Down->GetBinContent(i+1)/tot_down);
        cout<<dataPileUp_Down[i]<<endl;
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
cout<<"AddWeightFile"<<endl;
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

// _________________________________________________
string MyAnalysis::asString(double f) {
    std::ostringstream oss;
    oss<<f;
    string s = oss.str();
    return s;
}
