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
// C & C++
#include <iostream>
#include <vector>
#include <sstream>

#include <TROOT.h>

using namespace std;
//You should derive your own class from Analyse.
class MyAnalysis : public Analyse {
private:

    UInt_t currun;
    UInt_t curlumi;

    // define outputfile
    TFile *histfile;

    TH1D *hVtxN;
    TH1D *hVtxN_u;
    TH1D *hMuPt;
    TH1D *hNMuons;

    // debug
    bool debug;
    bool isMC;


public:
    MyAnalysis();
    virtual ~MyAnalysis();
    virtual Int_t AnalyseEvent(); //AnalyseEvent is a virtual function which is called for each event.
    double getAoverBError(double nA, double nB);
//    double getIsoPUE(Electron mu, double rho);
    void   AddPUWeightFile(string filename);
    bool   isGoodBJet(Jet j);
    bool   isGoodJet(Jet j);
    bool   hasJetIDLoose(Jet j);
    bool   hasJetSamePV(Jet j, Muon mu1, Muon mu2, std::vector < Vertex > PVs);
    double GetHLTEffScale();
    double GetMuEffScale(double muPt, double muEta);
    double GetElEffScale(double elPt, double elEta);
    string asString(double f);
};
//Constructor:
MyAnalysis::MyAnalysis() : Analyse(), currun(0), curlumi(0) {
    // don't touch, these are changed with the output file name // // //
gROOT->ProcessLine("#include <vector>");
    isMC=false;
    // // // //


    // load needed informations
    LoadTrigger();
    LoadBeamSpot();
    LoadPrimVertices();
    LoadMuons();
    LoadElectrons();
    LoadTracks();
    LoadAK4PFCHSJets();
    LoadMET();
    LoadGenParticles();
    LoadGenInfo();
    LoadAllGenParticles();
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
    // vertex
    hVtxN       = new TH1D("hVtxN", "N Vtx", 100, 0., 100.);
    hVtxN->GetXaxis()->SetTitle("N_{PV}");
    hVtxN->GetYaxis()->SetTitle("Candidates");

    hVtxN_u     = new TH1D("hVtxN_u", "N Vtx", 100, 0., 100.);
    hVtxN_u->GetXaxis()->SetTitle("N_{PV}");
    hVtxN_u->GetYaxis()->SetTitle("Candidates");

    hMuPt      = new TH1D("hMuPt", "mu Pt",    160, 0., 800.);
    hMuPt->GetXaxis()->SetTitle("p_{T #mu}[GeV/c]");
    hMuPt->GetYaxis()->SetTitle("Candidates/5.0GeV/c");

    hNMuons        = new TH1D("hNMuons", "muon multiplicity", 40, 0., 40.);
    hNMuons->GetXaxis()->SetTitle("N_{#mu}");
    hNMuons->GetYaxis()->SetTitle("Events");

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
    histfile->Write();
    histfile->Close();
}

// Analysis
Int_t MyAnalysis::AnalyseEvent() {
    hNMuons->Fill(NumMuons(), 1.);

    double pileupweight = 1.;

    // loop over prim vertices
    std::vector < Vertex > PVs;
    for(unsigned int v = 0; v < NumPrimVertices(); ++v) {
        PVs.push_back(PrimVertices(v));
    }
    hVtxN_u->Fill(PVs.size());
    hVtxN->Fill(PVs.size(), pileupweight);

    if(PVs.size() < 1) return (1);

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // MUONS //////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // select event with 2 good muons
    // get the muons
    std::vector< Muon > muonsVec;
    for(unsigned int i = 0; i < NumMuons(); ++i) {
        muonsVec.push_back(Muons(i));
    }

    if(muonsVec.size() < 2) return(1);


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
            isSamePVMuCutOK = true;
            // check the invariant mass
            TLorentzVector diMuon = muonsVec[i] + muonsVec[j];
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







    // fill muon plots
    for (unsigned int i = 0; i < muonsVec.size(); ++i) {
        hMuPt       ->Fill(muonsVec[i].Pt(),  pileupweight);
    }
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
    bool isPtOk  = j.Pt() >= 30;
    bool isEtaOk = TMath::Abs(j.Eta()) <= 4.7;

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

////________________________________
//double MyAnalysis::getIsoPUE(Electron el, double rho) {
//    double corrPU = el.IsoR3ECal() + el.IsoR3HCal() - rho*Pi()*0.09;
//    if(corrPU < 0) corrPU = 0.;
//    return (el.IsoR3Track() + corrPU)/el.Pt();
//}
//________________________________
void MyAnalysis::AddPUWeightFile(string filename) {
cout<<"AddPUWeightFile"<<endl;
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
