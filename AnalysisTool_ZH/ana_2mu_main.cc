////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
// ana_2mu_main.cc                                                                    //
//                                                                                    //
// This gives examples of how to access things with class Analyse running on          //
// the new (76X) RootMaker.                                                           //
//                                                                                    //
// 03 May 2016                                                                        //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////

// own
#include "AnalysisTool/Analyse.h"
// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TROOT.h>
// C & C++
#include <iostream>
#include <vector>
#include <sstream>
using namespace std;

// you should derive your own class from Analyse.
//________________________________
class MyAnalysis : public Analyse
{
  private:
    // pile-up
    std::vector < double > dataPileUp;
    std::vector < double > dataPileUp_Up;
    std::vector < double > dataPileUp_Down;

    UInt_t currun;
    UInt_t curlumi;
    double sumWeights;

    // cuts
    int    cVtxNdf;
    double cVtxZ;
    //_________muon cuts_________
    double cPtMu;
    double cEtaMu;
    double cPtMuMax;
    double cEtaMuMax;
    double cDxyMu;
    double cDzMu;
    double cIsoMu;
    // dimuon pair cuts
    double cInvMass;
    double c2MuPtCut;
    //_________electron cuts_________
    double cPtE;
    double cEtaE;
    double cIsoE;
    //_________jet cuts_________
    double cPtJet;
    double cEtaJet;
    //_________MET cuts_________
    double cMET;
    //_________misc cuts_________
    double massZ;

    // define output file
    TFile *histfile;

    // event histograms
    TH1D *hVtxN;
    TH1D *hVtxN_u;
    TH1D *hVtxN_after;
    TH1D *hVtxN_after_u;
    TH1D *hPUweight;
    // muon histograms
    TH1D *hMuPt;
    TH1D *hMuEta;
    TH1D *hMuCorrPt;
    // dimuon histograms
    TH1D *h2MuPt;
    TH1D *h2MuEta;
    TH1D *h2MuInvM;
    // electron histograms
    TH1D *hElPt;
    TH1D *hElEta;
    // MET histograms
    TH1D *hPFMETType1;
    // jet histograms
    TH1D *hJetPt;
    TH1D *hJetEta;

    // counters
    double nEv_Skim;
    double nEv_PV;
    double nEv_HLT;
    double nEv_Pt;
    double nEv_Eta;
    double nEv_MuID;
    double nEv_MuIso;
    double nEv_PtEtaMax;
    double nEv_JetID;


    // debug
    bool debug;

  public:
    MyAnalysis();
    virtual ~MyAnalysis();
    virtual Int_t AnalyseEvent(); // AnalyseEvent is a virtual function which is called for each event.

    double getAoverBError(double nA, double nB);
    void   AddPUWeightFile(string filename);
    bool   hasJetSamePV(Jet j, Muon mu1, Muon mu2, std::vector < Vertex > PVs);
    double GetHLTEffScale();
    double GetEffScale(double lepPt, double lepEta);
    string asString(double f);
};

// constructor
//________________________________
MyAnalysis::MyAnalysis() :
    Analyse(),
    currun(0),
    curlumi(0)
{
    gROOT->ProcessLine("#include <vector>");
    sumWeights = 0.0;

    debug = true;
    //debug = false;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // CUTS ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////





    // cuts
    cVtxNdf = 4;
    cVtxZ   = 24.;
    //_________muon cuts_________
    cPtMu     = 10.0; // GeV;
    cEtaMu    = 2.4; // choice here should depend on HLT
    cPtMuMax  = 20.0; // cut on trigger-matched mu
    cEtaMuMax = 2.4; // cut on trigger-matched mu
    cDxyMu    = 0.02; // cm
    cDzMu     = 0.14; // cm
    // isolation
    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
    cIsoMu = 0.15; // PF combined w/dB correction Tight
    //cIsoMu = 0.25; // PF combined w/dB correction Loose
    //cIsoMu = 0.05; // Tracker-based Tight
    //cIsoMu = 0.10; // Tracker-based Loose
    // dimuon pair cuts
    cInvMass  = 60.; // GeV
    c2MuPtCut = 38.; // GeV
    //_________electron cuts_________
    cPtE = 10.;
    cEtaE = 2.4;
    //_________jet cuts_________
    cPtJet = 30.;  // GeV;
    cEtaJet = 4.7;
    //_________MET cuts_________
    cMET = 40.; // GeV
    //_________misc cuts_________
    massZ = 91.1876; // GeV


    // load needed informations
    LoadTrigger();
    LoadBeamSpot();
    LoadPrimVertices();
    LoadMuons();
    LoadElectrons();
    LoadAK4PFCHSJets();
    LoadMET();
    LoadGenParticles();
    LoadGenInfo();
    UsePileUpInfo();
    TVector3 zDir(0,0,1);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // OUTPUT FILE ////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    histfile = new TFile("2mu_ana_out.root", "RECREATE");
    histfile->cd();


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // BOOK HISTOGRAMS ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    hPUweight   = new TH1D("hPUweight", "PU weight", 100, 0., 10.);
    hPUweight->GetXaxis()->SetTitle("PU weight");
    hPUweight->GetYaxis()->SetTitle("Candidates");

    // vertex
    hVtxN        = new TH1D("hVtxN", "N Vtx", 100, 0., 100.);
    hVtxN->GetXaxis()->SetTitle("N_{PV}");
    hVtxN->GetYaxis()->SetTitle("Candidates");

    hVtxN_u      = new TH1D("hVtxN_u", "N Vtx unweighted", 100, 0., 100.);
    hVtxN_u->GetXaxis()->SetTitle("N_{PV}");
    hVtxN_u->GetYaxis()->SetTitle("Candidates");

    hVtxN_after  = new TH1D("hVtxN_after", "N Vtx after selection", 100, 0., 100.);
    hVtxN_after->GetXaxis()->SetTitle("N_{PV}");
    hVtxN_after->GetYaxis()->SetTitle("Candidates");

    hVtxN_after_u = new TH1D("hVtxN_after_u", "N Vtx unweighted after selection", 100, 0., 100.);
    hVtxN_after_u->GetXaxis()->SetTitle("N_{PV}");
    hVtxN_after_u->GetYaxis()->SetTitle("Candidates");

    // muons
    hMuPt      = new TH1D("hMuPt", "mu Pt", 160, 0., 800.);
    hMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    hMuEta     = new TH1D("hMuEta", "mu #eta", 100, -5., 5.);
    hMuEta->GetXaxis()->SetTitle("#eta_{#mu}");
    hMuEta->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    hMuCorrPt      = new TH1D("hMuCorrPt", "mu rochester Pt", 160, 0., 800.);
    hMuCorrPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuCorrPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    // dimuons
    h2MuInvM      = new TH1D("h2MuInvM", "dimuon inv mass", 400, 0., 800.);
    h2MuInvM->GetXaxis()->SetTitle("M_{#mu#mu}[GeV/c]");
    h2MuInvM->GetYaxis()->SetTitle("Candidates/2GeV/c");

    h2MuPt      = new TH1D("h2MuPt", "dimuon Pt", 160, 0., 800.);
    h2MuPt->GetXaxis()->SetTitle("p_{T #mu#mu}[GeV/c]");
    h2MuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    h2MuEta     = new TH1D("h2MuEta", "dimuon #eta", 100, -5., 5.);
    h2MuEta->GetXaxis()->SetTitle("#eta_{#mu#mu}");
    h2MuEta->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    // electrons
    hElPt      = new TH1D("hElPt", "e Pt", 160, 0., 800.);
    hElPt->GetXaxis()->SetTitle("p_{T e}[GeV/c]");
    hElPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    hElEta     = new TH1D("hElEta", "e #eta", 100, -5., 5.);
    hElEta->GetXaxis()->SetTitle("#eta_{e}");
    hElEta->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    // MET
    hPFMETType1         = new TH1D("hPFMETType1", "PFMETType1", 4000, 0., 2000.);
    hPFMETType1->GetXaxis()->SetTitle("MET [GeV/c^{2}]");
    hPFMETType1->GetYaxis()->SetTitle("Candidates/0.5[GeV/c^{2}]");

}

MyAnalysis::~MyAnalysis()
{
    histfile->Write();
    histfile->Close();
}

// Analysis
//________________________________
Int_t MyAnalysis::AnalyseEvent()
{
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                           //
    // EVENT AND OBJECT SELECTION + EVENT PLOTS                                                  //
    //                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PILEUP WEIGHTING ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    double pileupweight = 0.;
    if (IsData()) {
        pileupweight = 1.;
        sumWeights += 1.;
    } else {
        pileupweight =  GenWeight() * GetPileUpWeight(dataPileUp);
        sumWeights += GenWeight();
    }
    ++nEv_Skim;

    hPUweight->Fill(pileupweight);



    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TRIGGER ////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // vector of HLT paths to check - event will pass if any of these triggers fired
    std::vector<string> HLTs;

    HLTs.push_back("IsoMu20");
    HLTs.push_back("IsoTkMu20");

    if (not (EventPassesHLT(HLTs)) ) return (1);
    ++nEv_HLT;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PRIMARY VERTEX /////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // select event with a valid PV
    // loop over prim vertices
    std::vector < Vertex > PVs;
    for(size_t v = 0; v < NumPrimVertices(); ++v) {
        PVs.push_back(PrimVertices(v));
    }

    hVtxN_u->Fill(PVs.size());
    hVtxN->Fill(PVs.size(), pileupweight);

    if(PVs.size() < 1) return (1);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MUONS //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // select event with 2 good muons
    // vector to store good muons
    std::vector< Muon > muonsVec;
    // things we will check for each muon
    // isChi2NdfOK, isNMuonHitsOK, isNMatchedStationsOK, isNPixelHitsOK, isNTrackerLayersOK, all taken 
    // care of by Tight Muon ID isDxyOK and isDzOK are also considered in the Tight ID, but are not 
    // identical to the cuts used previously, we can consider removing them too
    bool isGAndTr   = false;
    bool isPtCutOK  = false;
    bool isEtaCutOK = false;
    bool isDxyOK    = false;
    bool isDzOK     = false;
    bool isMuIDOk   = false;
    bool isMuIsoOK  = false;
    // event counters
    int nMuPtEtaMax = 0;
    int nLooseMus   = 0;
    int nTightMus   = 0;
    int nMedMus     = 0;

    // loop over muons
    for(size_t i = 0; i < NumMuons(); ++i) {
        // kinematic cuts
        // require Global and Tracker muon (also available: IsStandAlone(), IsCalo())
        if( not(Muons(i).IsGlobal()) && Muons(i).IsTracker() ) continue;
        isGAndTr = true;
        // also available: rochester corrected pT, RochCorrPt()
        if( not(Muons(i).Pt() >= cPtMu) ) continue;
        isPtCutOK = true;
        if( not(TMath::Abs(Muons(i).Eta()) <= cEtaMu) ) continue;
        isEtaCutOK = true;

        // check muon ID
        // require tight muID (also available: IsMediumMuon() and IsLooseMuon())
        if( not(Muons(i).IsTightMuon()) ) continue;
        isMuIDOk = true;

        // muon isolation
        // below are examples for IsoPFR4dBCombRel() and IsoR4TrackRel()
        // also available: IsoR3CombinedRel()
        // relative PF combined w/dB correction, cone of r=0.4 (also available: IsoPFR3dBCombRel()) 
        if( not(Muons(i).IsoPFR4dBCombRel() < cIsoMu) ) continue;

        // relative tracker-based cone of r=0.4 (also available: IsoR3TrackRel())
        //if( not(Muons(i).IsoR4TrackRel() < cIsoMu) ) continue;

        isMuIsoOK = true;

        // at least one trigger-matched muon in the event must ALSO pass this cut:
        if( Muons(i).MatchesHLT(HLTs) && (Muons(i).Pt() > cPtMuMax && TMath::Abs(Muons(i).Eta()) <= cEtaMuMax)) ++nMuPtEtaMax;

        // save the good mus
        muonsVec.push_back(Muons(i));
    }

    if(isPtCutOK) ++nEv_Pt;
    if(isEtaCutOK) ++nEv_Eta;
    if(isMuIDOk) ++nEv_MuID;
    if(isMuIsoOK) ++nEv_MuIso;

    if(nMuPtEtaMax > 0) ++nEv_PtEtaMax;
    else return(1);

    // keep events with at least two good muons
    if(muonsVec.size() < 2) return(1);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ELECTRONS //////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // vector to store good electrons
    std::vector< Electron > electronsVec;

    // loop over electrons
    for(size_t e = 0; e < NumElectrons(); ++e) {
        // kinematic cuts
        if( not(Electrons(e).Pt() >= cPtE) ) continue;
        if( not(TMath::Abs(Electrons(e).Eta()) <= cEtaE) ) continue;

        // check electron ID
        // cut-based ID: IsTightElectron(), IsLooseElectron(), IsMediumElectron()
        // mva: WP80_v1(), WP90_v1()
        if( not(Electrons(e).IsMediumElectron()) ) continue;

        // electron isolation
        // available: IsoPFR3dBCombRel(), IsoPFR3RhoCombRel()
        // HOWEVER there is isolation requirements built into the electron ID 
        // so you would have to find your own value for cIsoE
        //if( not(Electrons(e).IsoPFR3RhoCombRel() < cIsoE) ) continue;

        // save the good es
        electronsVec.push_back(Electrons(e));
    }

    // keep events with at least two good electrons
    //if(electronsVec.size() < 2) return(1);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // JETS ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // jet cleaning is done in RootMaker against electrons, taus, and muons with dR of 0.3
    // tighter cleaning can be done here if desired
    std::vector< Jet > jetsVec;
    // cuts on jets
    bool isJetIDOK = false;

    for (size_t j = 0; j < NumAK4PFCHSJets(); ++j) {
        // kinematic cuts
        if( not(AK4PFCHSJets(j).Pt() > cPtJet)) continue;
        if( not(TMath::Abs(AK4PFCHSJets(j).Eta()) <= cEtaJet)) continue;

        // tighter cleaning against muons and electrons we selected:
	//bool isdROK = true;
	//for (size_t i = 0; i < muonsVec.size(); ++i) {
        //    double dR = AK4PFCHSJets(j).DeltaR(muonsVec[i]);
        //    if(dR < cDR) isdROK = false;
        //}
        //if (!isdROK) continue;
        //for (size_t i = 0; i < electronsVec.size(); ++i) {
        //    double dR = AK4PFCHSJets(j).DeltaR(electronsVec[i]);
        //    if(dR < cDR) isdROK = false;
        //}
        //if (!isdROK) continue;

        // jet ID
        // require loose jet ID (also available: IsTightJet(), IsTightLepVetoJet())
        if( not(AK4PFCHSJets(j).IsLooseJet()) ) continue;
        isJetIDOK = true;

        // store the jets
        jetsVec.push_back(AK4PFCHSJets(j));
    }

    if (isJetIDOK) ++nEv_JetID;

    // keep events with at least one good jet
    //if(jetsVec.size() < 1) return(1);


    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MET ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // the only MET available in miniAOD is type1
    // MET in the event is stored in a TLorentzVector called PFMETTYPE1()
    float evtMET = PFMETTYPE1().Pt();

    // cut on MET
    //if( not(evtMET < cMET) ) return(1);

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
    for (size_t i = 0; i < muonsVec.size(); ++i) {
        // 2nd loop over muons
        for (size_t j = i+1; j < muonsVec.size(); ++j) {
            // check muon charge
            if( not(muonsVec[i].Charge()*muonsVec[j].Charge() < 0) ) continue;
            isChargeMuCutOK = true;     
            // dz cut
            double dzMu = muonsVec[i].Dz() - muonsVec[j].Dz();
            //if( not(dzMu < cDzMu) ) continue;
            isSamePVMuCutOK = true;
            // check the invariant mass
            TLorentzVector diMuon = muonsVec[i] + muonsVec[j];
            if( not(diMuon.M() > cInvMass) ) continue;
            isInvMassMuCutOK = true;

            diMuonVec.push_back(diMuon);
            // make two arrays: if i and j make a pair, put each index in a separate array
            // so Muon_idx1[0] and Muon_idx2[0] make a pair, and Muon_idx1[1] and Muon_idx2[1] make another, and so on
            Muon_idx1.push_back(i);
            Muon_idx2.push_back(j);
        } // end 2nd loop over muons
    } //end 1st loop over muons

    if (!isSamePVMuCutOK) return(1);
    if (!isChargeMuCutOK) return(1);
    if (!isInvMassMuCutOK) return(1);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                           //
    // FILL PLOTS                                                                                //
    //                                                                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // fill nPV plots
    hVtxN_after_u->Fill(PVs.size());
    hVtxN_after->Fill(PVs.size(), pileupweight);

    // fill muon plots
    for (size_t i = 0; i < muonsVec.size(); ++i) {
        hMuPt->Fill(muonsVec[i].Pt(), pileupweight);
        hMuEta->Fill(muonsVec[i].Eta(), pileupweight);
        hMuCorrPt->Fill(muonsVec[i].RochCorrPt(), pileupweight);
    }
    // fill dimuon plots
    for (size_t i = 0; i < diMuonVec.size(); ++i) {
        h2MuInvM->Fill(diMuonVec[i].M(), pileupweight);
        h2MuPt->Fill(diMuonVec[i].Pt(), pileupweight);
        h2MuEta->Fill(diMuonVec[i].Eta(), pileupweight);
    }

    // fill electron plots
    for (size_t i = 0; i < electronsVec.size(); ++i) {
        hElPt->Fill(electronsVec[i].Pt(), pileupweight);
        hElEta->Fill(electronsVec[i].Eta(), pileupweight);
    }

    // fill jet plots
    for (size_t i = 0; i < jetsVec.size(); ++i) {
        hJetPt->Fill(jetsVec[i].Pt(), pileupweight);
        hJetEta->Fill(jetsVec[i].Eta(), pileupweight);
    }

    // fill MET plots
    hPFMETType1->Fill(evtMET, pileupweight);




    return(1);
}

//________________________________
int main()
{
    // create an instance of your analysis class:
    MyAnalysis ana;
    string namebuf;
    string filename;
    cout<<"Enter one filename and press enter. To stop this enter END exactly."<<endl;
    while(true) {
        cin >> namebuf;
        if(namebuf == "END") break;
        else {
            UInt_t slashpos = namebuf.find_last_of("/");
            if(slashpos == namebuf.size()) {
                filename = namebuf;
            } else {
                filename = namebuf.substr(slashpos+1);
            }
            cout<<filename<<endl;
// NOTE
            if(filename.find("LUMI_") == 0) ana.AddLumiFile(namebuf.c_str());
            else if(filename.find("MyDataPileupHistogram") == 0) ana.AddPUWeightFile(namebuf.c_str());
            else ana.AddFile(namebuf.c_str());
        }
    }
//    ana.AddPUWeightFile("MyDataPileupHistogram_v2.root");

    // Loop will start to run the analysis on the specified range or on all events if no range is given.
    ana.SetPrintInfo(10000);
    //ana.EnableDuplicateCheck();
    //ana.Loop(0,5000);
    ana.Loop();
    ana.PrintLumiOfRuns();
    cout<<"Lumi: "<<ana.GetLumi()<<endl;
}

//________________________________
double MyAnalysis::getAoverBError(double nA, double nB)
{
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
void MyAnalysis::AddPUWeightFile(string filename)
{
    TFile weightfile(filename.c_str(), "READ");
    if(weightfile.IsZombie()) {
        cerr << "ERROR AddLumiFile: " << filename << " is not a valid file."<<endl;
    }
    TH1D *Vertex      = (TH1D *)weightfile.Get("pileup");
    Double_t tot      = Vertex->Integral();

    for(int i=0; i<Vertex->GetNbinsX(); i++) {
//        dataPileUp.push_back(Vertex->GetBinContent(i+1)/tot);
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
double MyAnalysis::GetEffScale(double lepPt, double lepEta)
{
    double my_ratio = 1.;
    double eff_data;
    double eff_MC;
    return my_ratio;
}

// _________________________________________________
string MyAnalysis::asString(double f)
{
    std::ostringstream oss;
    oss<<f;
    string s = oss.str();
    return s;
}
