#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

int createReWeightPieces(TString sample, TString pileUpDataFile) {
    gROOT->SetBatch();
    TFile* outFile = new TFile("ReWeightPiece_" + sample + ".root","RECREATE");

    TFile *nPVtrees = new TFile("nPV_trees.root","READ");
    TFile *pileUpData = new TFile(pileUpDataFile,"READ");

    // Data pileup distribution
    TH1D *nPV_Run2012 = pileUpData->Get("pileup");
    nPV_Run2012->SetName("nPV_Run2012");
    int nbins = nPV_Run2012->GetNbinsX();
    double xmin = nPV_Run2012->GetXaxis()->GetXmin();
    double xmax = nPV_Run2012->GetXaxis()->GetXmax();
    outFile->cd();

    Float_t genweight = 0;
    Float_t weight_0 = 0;
    Float_t weight_1 = 0;
    TTree *nPVs = (TTree*)nPVtrees->Get(sample);
    nPVs->SetBranchAddress("genweight",&genweight);
    int nEntries = nPVs->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        if (i%10000==0) cout<<sample<<" event "<<i<<"... "<<i*100/nEntries<<"\% done"<<endl;
        nPVs->GetEntry(i);
        Float_t nW = genweight;
        weight_0 += nW;
        if (nW > 0.00) weight_1 += 1.0;
        else if (nW < 0.00) weight_1 += -1.0;
    }
    TH1D *nPV_sample = new TH1D("nPV_" + sample,"nPV_" + sample, nbins, xmin, xmax);
    nPVs->Draw("numtruepileupinteractions>>nPV_" + sample);
    nPV_sample->SetName("nPV_" + sample);
    TH1F *sumWeights = new TH1F("sumWeights_" + sample, "sumWeights_" + sample, 5,0.,5.);
    sumWeights->SetBinContent(2,weight_0);
    sumWeights->SetBinContent(4,weight_1);

    TH1D *nPV_Run2012_norm = (TH1D*)(nPV_Run2012->Clone("nPV_Run2012_norm"));
    TH1D *nPV_sample_norm = (TH1D*)(nPV_sample->Clone("nPV_" + sample + "_norm"));
    nPV_sample_norm->Scale(1./(nPV_sample_norm->Integral()));


    TH1D *nPV_sample_ratio = (TH1D*)(nPV_sample->Clone("nPV_" + sample + "_ratio"));
    for (int ibin = 0; ibin < nbins; ++ibin) {
        double nPV2012   = nPV_Run2012_norm->GetBinContent(nPV_Run2012_norm->GetXaxis()->FindBin(ibin));
        double nPVSample = nPV_sample_norm->GetBinContent( nPV_sample_norm->GetXaxis()->FindBin(ibin));

        if (nPVSample != 0) nPV_sample_ratio->SetBinContent(ibin, nPV2012/nPVSample);
        else nPV_sample_ratio->SetBinContent(ibin, nPVSample);
    }

    nPV_Run2012->Delete();
    nPV_Run2012_norm->Delete();
    nPV_sample->Write();
    nPV_sample_norm->Write();
    nPV_sample_ratio->Write();
    sumWeights->Write();
    outFile->Write();
    outFile->Close();
    return 0;
}
