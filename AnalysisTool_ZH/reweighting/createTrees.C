#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

int createTrees(TString indir, TString sample) {
    gROOT->SetBatch();
  
    cout<<"sample name = "<<sample<<endl;
    cout<<"indir = "<<indir<<endl;

    TFile* infoFile = new TFile("ReWeightInfo_" + sample + ".root","RECREATE");
    infoFile->cd();

    TChain* chain = new TChain("makeroottree/AC1B");
    chain->Add(indir);
    int nEntries = chain->GetEntries();
    chain->SetBranchStatus("*", 0);
    chain->SetBranchStatus("numtruepileupinteractions", 1);
    chain->SetBranchStatus("genweight", 1);
    TTree* chain_new = chain->CloneTree(0);
    for (int i = 0; i < nEntries; i++) {
        if (i%10000==0) cout<<sample<<" event "<<i<<"... "<<i*100/nEntries<<"\% done"<<endl;
        chain->GetEntry(i);
        chain_new->Fill(); }
    chain_new->SetName(sample);
    chain_new->Write();

    infoFile->Write();
    infoFile->Close();

    return 0;
}
