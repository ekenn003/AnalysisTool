#ifndef ACANALYSISBOX
#define ACANALYSISBOX
#include "Analyse.h"
#include "DataForm.h"

using namespace std;
using namespace TMath;

class AnalysisBox {
private:
    Analyse *an;

public:
    AnalysisBox(Analyse *analyse);
    Long64_t GetNumAddedEvents() const {
        return(an->GetNumAddedEvents());
    }
    string GetCurrentFileName() const {
        return(an->GetCurrentFileName());
    }
    Int_t Batch_MyId() const {
        return(an->Batch_MyId());
    }
    Int_t Batch_NumJobs() const {
        return(an->Batch_NumJobs());
    }

    Double_t Number() const {
        return(an->Number());
    }
    UInt_t Run() const {
        return(an->Run());
    }
    UInt_t LumiBlock() const {
        return(an->LumiBlock());
    }
    UInt_t TimeUnix() const {
        return(an->TimeUnix());
    }
    UInt_t TimeMicroSec() const {
        return(an->TimeMicroSec());
    }
    Double_t AK4PFRho() const {
        return(an->AK4PFRho());
    }
    Double_t AK4PFSigma() const {
        return(an->AK4PFSigma());
    }

    //RECO-level information
    BeamSpot GetBeamSpot() const {
        return(an->GetBeamSpot());
    }

    Muon Muons(UInt_t n, int correction = 0) const {
        return(an->Muons(n, correction));
    }
    UInt_t NumMuons() const {
        return(an->NumMuons());
    }

    Electron Electrons(UInt_t n, int correction = 0) const {
        return(an->Electrons(n, correction));
    }
    UInt_t NumElectrons() const {
        return(an->NumElectrons());
    }

    Tau Taus(UInt_t n) const {
        return(an->Taus(n));
    }
    UInt_t NumTaus() const {
        return(an->NumTaus());
    }

    Photon Photons(UInt_t n) const {
        return(an->Photons(n));
    }
    UInt_t NumPhotons() const {
        return(an->NumPhotons());
    }

    TLorentzVector PFMET() const {
        return(an->PFMET());
    }
    TLorentzVector PFMETTYPE1() const {
        return(an->PFMETTYPE1());
    }
    TLorentzVector PFMETPUPPITYPE1() const {
        return(an->PFMETPUPPITYPE1());
    }
    TLorentzVector PFMETTYPE0TYPE1() const {
        return(an->PFMETTYPE0TYPE1());
    }

    Jet AK4PFCHSJets(UInt_t n) const {
        return(an->AK4PFCHSJets(n));
    }
    UInt_t NumAK4PFCHSJets() const {
        return(an->NumAK4PFCHSJets());
    }

/*
    Track Tracks(UInt_t n) const {
        return(an->Tracks(n));
    }
    UInt_t NumTracks() const {
        return(an->NumTracks());
    }
*/
    Vertex PrimVertices(UInt_t n) const {
        return(an->PrimVertices(n));
    }
    UInt_t NumPrimVertices() const {
        return(an->NumPrimVertices());
    }
/*
    SuperCluster SuperClusters(UInt_t n, TVector3 refpoint) const {
        return(an->SuperClusters(n, refpoint));
    }
    SuperCluster SuperClusters(UInt_t n) const {
        return(an->SuperClusters(n));
    }
    UInt_t NumSuperClusters() const {
        return(an->NumSuperClusters());
    }
*/
    Int_t NumPileUpInteractionsMinus() const {
        return(an->NumPileUpInteractionsMinus());
    }
    Int_t NumPileUpInteractions() const {
        return(an->NumPileUpInteractions());
    }
    Int_t NumPileUpInteractionsPlus() const {
        return(an->NumPileUpInteractionsPlus());
    }
    Float_t NumTruePileUpInteractions() const {
        return(an->NumTruePileUpInteractions());
    }
    Double_t GetPileUpMaxWeight(vector<Double_t> &datadist) const {
        return(an->GetPileUpMaxWeight(datadist));
    }
    Double_t GetPileUpWeight(vector<Double_t> &datadist) const {
        return(an->GetPileUpWeight(datadist));
    }
    Double_t GetPrimVertexMaxWeight(vector<Double_t> &datadist) const {
        return(an->GetPrimVertexMaxWeight(datadist));
    }
    Double_t GetPrimVertexWeight(vector<Double_t> &datadist) const {
        return(an->GetPrimVertexWeight(datadist));
    }
    Int_t NumGoodPrimVertices() const {
        return(an->NumGoodPrimVertices());
    }

    //generator-level information
    Double_t GenWeight() const {
        return(an->GenWeight());
    }
    Double_t GenId1() const {
        return(an->GenId1());
    }
    Double_t Genx1() const {
        return(an->Genx1());
    }
    Double_t GenId2() const {
        return(an->GenId2());
    }
    Double_t Genx2() const {
        return(an->Genx2());
    }
    Double_t GenScale() const {
        return(an->GenScale());
    }

    GenParticle AllGenParticles(UInt_t n) const {
        return(an->AllGenParticles(n));
    }
    UInt_t NumAllGenParticles() const {
        return(an->NumAllGenParticles());
    }

    GenLightParticle GenParticles(UInt_t n) const {
        return(an->GenParticles(n));
    }
    UInt_t NumGenParticles() const {
        return(an->NumGenParticles());
    }

    GenJet GenAK4Jets(UInt_t n) const {
        return(an->GenAK4Jets(n));
    }
    UInt_t NumGenAK4Jets() const {
        return(an->NumGenAK4Jets());
    }

    TLorentzVector GenMETCalo() const {
        return(an->GenMETCalo());
    }
    TLorentzVector GenMETTrue() const {
        return(an->GenMETTrue());
    }

    Long64_t Processed() const {
        return(an->Processed());
    }

    //trigger information
    bool GetL1Trigger(UInt_t bit) const {
        return(an->GetL1Trigger(bit));
    }
    bool GetL1TriggerBits(UInt_t bit) const {
        return(an->GetL1TriggerBits(bit));
    }
    bool GetHLTrigger(UInt_t index) const {
        return(an->GetHLTrigger(index));
    }
    Int_t GetHLTrigger(vector<string> triggernames) const {
        return(an->GetHLTrigger(triggernames));
    }
    Int_t GetHLTriggerIndex(string triggername) {
        return(an->GetHLTriggerIndex(triggername));
    }
    string GetHLTriggerName(UInt_t index) {
        return(an->GetHLTriggerName(index));
    }
    Int_t GetHLTPrescale(UInt_t triggerindex) {
        return(an->GetHLTPrescale(triggerindex));
    }
    Int_t GetNumHLTriggers() {
        return(an->GetNumHLTriggers());
    }

    TriggerSelection *GetTriggerSelection(string id) {
        return(an->GetTriggerSelection(id));
    }

    Int_t IsLumiAvailable() const {
        return(an->IsLumiAvailable());
    }
    Double_t GetInstLumi() const {
        return(an->GetInstLumi());
    }
    Double_t GetAvgPU() const {
        return(an->GetAvgPU());
    }
    Double_t GetLumi(Int_t format = 0) {
        return(an->GetLumi(format));
    }
    //Double_t GetLumiBlockLumi() {return(an->GetLumiBlockLumi());}
    bool IsData() const {
        return(an->IsData());
    }
    bool IsMC() const {
        return(an->IsMC());
    }
};

#endif
