#include "TFile.h"
#include "TH1.h"

int combinePileUpDists() {
    gROOT->SetBatch();

    TFile* outFile = new TFile("PileUpDistRun2015D.root","RECREATE");

    TFile *file_65550 = new TFile("PileUpData_65550.root","READ");
    TFile *file_69000 = new TFile("PileUpData_69000.root","READ");
    TFile *file_72450 = new TFile("PileUpData_72450.root","READ");

    TH1D *h_65550 = (TH1D *)file_65550->Get("pileup");
    TH1D *h_69000 = (TH1D *)file_69000->Get("pileup");
    TH1D *h_72450 = (TH1D *)file_72450->Get("pileup");

    h_65550->SetName("pileup65550");
    h_69000->SetName("pileup69000");
    h_72450->SetName("pileup72450");

    outFile->cd();
    h_65550->Write();
    h_69000->Write();
    h_72450->Write();
    outFile->Write();
    outFile->Close();

    return 0;
}
