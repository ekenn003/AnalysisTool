#include "DataForm.h"
#include "Analyse.h"


//GenJet
void GenJet::ScaleVisible(bool visible) {
    if(myvisible == false && visible == true) {
        Double_t f = 1. - InvisibleEnergy()/E();
        this->SetPxPyPzE(f*Px(), f*Py(), f*Pz(), f*E());
        myvisible = true;
    }
    if(myvisible == true && visible == false) {
        Double_t f = 1. - InvisibleEnergy()/E();
        this->SetPxPyPzE(Px()/f, Py()/f, Pz()/f, E()/f);
        myvisible = false;
    }
}

GenBasicParticle::GenBasicParticle(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid) :
    TLorentzVector(Px, Py, Pz, E),
    vertex(X, Y, Z),
    status(Status),
    pdgid(Pdgid) {

}

GenLightParticle::GenLightParticle(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Info, Int_t Indirectmother) :
    GenBasicParticle(E, Px, Py, Pz, X, Y, Z, Status, Pdgid),
    info(Info),
    indirectmother(Indirectmother) {
}

GenParticle::GenParticle(const Analyse *Myanalyse, UInt_t Myindex, Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Motherfirst, UInt_t Mothernum, UInt_t Daughterfirst, UInt_t Daughternum) :
    GenBasicParticle(E, Px, Py, Pz, X, Y, Z, Status, Pdgid),
    myanalyse(Myanalyse),
    myindex(Myindex),
    motherfirst(Motherfirst),
    mothernum(Mothernum),
    daughterfirst(Daughterfirst),
    daughternum(Daughternum) {
}

UInt_t GenParticle::GetMotherIndex(UInt_t num) const {
    return(myanalyse->genallparticles_mothers[motherfirst+num]);
}

UInt_t GenParticle::GetDaughterIndex(UInt_t num) const {
    return(myanalyse->genallparticles_daughters[daughterfirst+num]);
}

GenParticle GenParticle::GetMother(UInt_t num) const {
    if(num >= mothernum) {
        cout << "GenParticle::GetMother: mother index out of range. return this." << endl;
        return(*this);
    }
    return(myanalyse->AllGenParticles(GetMotherIndex(num)));
}


GenParticle GenParticle::GetDaughter(UInt_t num) const {
    if(num >= daughternum) {
        cout << "GenParticle::GetDaughter: daughter index out of range. return this." << endl;
        return(*this);
    }
    return(myanalyse->AllGenParticles(GetDaughterIndex(num)));
}

bool GenParticle::HasAnyMotherPDGId(Int_t pdgid, bool antiparticle) {
    vector<UInt_t> visited;
    return(HasAnyMotherPDGId(visited, pdgid, antiparticle));
}

bool GenParticle::HasAnyMotherPDGId(vector<UInt_t> &visited, Int_t pdgid, bool antiparticle) {
    for(UInt_t i = 0 ; i < NumMothers() ; i++) {
        bool isvisited = false;
        for(UInt_t u = 0 ; u < visited.size() ; u++) {
            if(visited[u] == GetMotherIndex(i)) {
                isvisited = true;
                break;
            }
        }
        if(isvisited) continue;
        visited.push_back(GetMotherIndex(i));
        if(antiparticle && Abs(pdgid) == GetMother(i).PDGId()) return(true);
        if(!antiparticle && pdgid == GetMother(i).PDGId()) return(true);
        if(GetMother(i).HasAnyMotherPDGId(visited, pdgid, antiparticle)) return(true);
    }
    return(false);
}

TrackComposedParticle::TrackComposedParticle(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t *Cov) :
    TVector3(Vx,Vy,Vz),
    covmatrix(3),
    chi2(Chi2),
    ndof(Ndof) {
    covmatrix(0,0) = Cov[0];
    covmatrix(0,1) = Cov[1];
    covmatrix(0,2) = Cov[2];
    covmatrix(1,0) = Cov[1];
    covmatrix(1,1) = Cov[3];
    covmatrix(1,2) = Cov[4];
    covmatrix(2,0) = Cov[2];
    covmatrix(2,1) = Cov[4];
    covmatrix(2,2) = Cov[5];
}



void TrackComposedParticle::AddTrack(Track &newtrack) {
    tracks.push_back(newtrack);
}

TVector3 TrackComposedParticle::Momentum() const {
    TVector3 momentum(0., 0., 0.);
    for(int i = 0 ; i < NumTracks() ; ++i) {
        momentum += Tracks(i).Vect();
    }
    return(momentum);
}

TVector3 TrackComposedParticle::Distance(Vertex &vertex) const {
    TVector3 mom=Momentum();
    double l = (mom*(*this-vertex))/mom.Mag2();
    return(*this-l*mom);
}

Double_t TrackComposedParticle::VertexSig3D(Vertex &vertex) const {
    TVector3 dist = *this - vertex;
    return(dist.Mag2()/Sqrt(((vertex.CovMatrix()+CovMatrix())* dist)*dist));
}

Double_t TrackComposedParticle::VertexSig2D(Vertex &vertex) const {
    TVector3 dist = *this - vertex;
    dist.SetZ(0);
    return(dist.Mag2()/Sqrt(((vertex.CovMatrix()+CovMatrix())* dist)*dist));
}

TriggerObject::TriggerObject(const Analyse *ma, const vector<string> *Triggernames, UInt_t Trigger) :
    MA(ma),
    trigger(Trigger),
    triggernames(Triggernames) {

}

Int_t TriggerObject::Trigger(string triggername) const {
    Int_t index = -1;
    for(UInt_t i = 0 ; i < triggernames->size() ; i++) {
        if(triggername == triggernames->at(i)) {
            index = i;
            break;
        }
    }

    if(index == -1) return(0);
    bool result = (trigger & 1<<index) != 0;
    size_t pos = triggername.find(":");
    if(pos != string::npos) {
        triggername = triggername.substr(0, pos);
    }
    Int_t glindex = MA->GetHLTriggerIndex(triggername);
    if(MA->GetHLTPrescale(glindex) == 1) {
        if(result) return(1);
        else return(-1);
    }
    if(result) return(2);
    else return(-2);
}



Muon::Muon(const Analyse *ma, UInt_t n, int correction) :
    TLorentzVector(ma->muon_px->at(n), ma->muon_py->at(n), ma->muon_pz->at(n), sqrt(ma->muon_px->at(n)*ma->muon_px->at(n)+ma->muon_py->at(n)*ma->muon_py->at(n)+ma->muon_pz->at(n)*ma->muon_pz->at(n)+MuonMassQ)),
    TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTMuonNames(), ma->muon_trigger->at(n)),
    pterror(ma->muon_pterror->at(n)),
    chi2(ma->muon_chi2->at(n)),
    ndof(ma->muon_ndof->at(n)),
    is_global(ma->muon_is_global->at(n)),
    is_tracker(ma->muon_is_tracker->at(n)),
    is_standalone(ma->muon_is_standalone->at(n)),
    is_pf_muon(ma->muon_is_pf_muon->at(n)),
    is_tight_muon(ma->muon_is_tight_muon->at(n)),
    is_medium_muon(ma->muon_is_medium_muon->at(n)),
    is_loose_muon(ma->muon_is_loose_muon->at(n)),
    isolationr3track(ma->muon_isolationr3track->at(n)),
    isolationr3ntrack(ma->muon_isolationr3ntrack->at(n)),
    isolationr3ecal(ma->muon_isolationr3ecal->at(n)),
    isolationr3hcal(ma->muon_isolationr3hcal->at(n)),
    pfisolationr3_sumchargedhadronpt(ma->muon_pfisolationr3_sumchargedhadronpt->at(n)),
    pfisolationr3_sumchargedparticlept(ma->muon_pfisolationr3_sumchargedparticlept->at(n)),
    pfisolationr3_sumneutralhadronet(ma->muon_pfisolationr3_sumneutralhadronet->at(n)),
    pfisolationr3_sumphotonet(ma->muon_pfisolationr3_sumphotonet->at(n)),
    pfisolationr3_sumneutralhadronethighthreshold(ma->muon_pfisolationr3_sumneutralhadronethighthreshold->at(n)),
    pfisolationr3_sumphotonethighthreshold(ma->muon_pfisolationr3_sumphotonethighthreshold->at(n)),
    pfisolationr3_sumpupt(ma->muon_pfisolationr3_sumpupt->at(n)),
    pfisolationr4_sumchargedhadronpt(ma->muon_pfisolationr4_sumchargedhadronpt->at(n)),
    pfisolationr4_sumchargedparticlept(ma->muon_pfisolationr4_sumchargedparticlept->at(n)),
    pfisolationr4_sumneutralhadronet(ma->muon_pfisolationr4_sumneutralhadronet->at(n)),
    pfisolationr4_sumphotonet(ma->muon_pfisolationr4_sumphotonet->at(n)),
    pfisolationr4_sumneutralhadronethighthreshold(ma->muon_pfisolationr4_sumneutralhadronethighthreshold->at(n)),
    pfisolationr4_sumphotonethighthreshold(ma->muon_pfisolationr4_sumphotonethighthreshold->at(n)),
    pfisolationr4_sumpupt(ma->muon_pfisolationr4_sumpupt->at(n)),
    pfisolationr4_dBrel(ma->muon_pfisolationr4_dBrel->at(n)),
    ecalenergy(ma->muon_ecalenergy->at(n)),
    hcalenergy(ma->muon_hcalenergy->at(n)),
    charge(ma->muon_charge->at(n)),
    numchambers(ma->muon_numchambers->at(n)),
    numchamberswithsegments(ma->muon_numchamberswithsegments->at(n)),
    numvalidmuonhits(ma->muon_numvalidmuonhits->at(n)),
    nummatchedstations(ma->muon_nummatchedstations->at(n)),

    matches_IsoMu20(ma->muon_matches_IsoMu20->at(n)),
    matches_IsoTkMu20(ma->muon_matches_IsoTkMu20->at(n)),

    trackermuonquality(ma->muon_trackermuonquality->at(n)) {
}

bool Muon::MatchesHLT(std::vector<string> hltnames) const {
    bool matchedAny = false;
    // range-based 'for' loops are not allowed in C++98 mode :(
    //for (auto &path : hltnames) {
    for (size_t i = 0; i < hltnames.size(); i++) {
        if (hltnames[i] == "IsoMu20")        matchedAny = matchedAny || (bool(matches_IsoMu20));
        else if (hltnames[i] == "IsoTkMu20") matchedAny = matchedAny || (bool(matches_IsoTkMu20));

        else { cerr<<"HLT name "<<hltnames[i]<<" not recognized."<<endl; throw; }
    }

    return(matchedAny);
}

Int_t Muon::NumStations() const {
    Int_t result(0);
    if(trackermuonquality & 1<<24) result++;
    if(trackermuonquality & 1<<25) result++;
    if(trackermuonquality & 1<<26) result++;
    if(trackermuonquality & 1<<27) result++;
    if(trackermuonquality & 1<<28) result++;
    if(trackermuonquality & 1<<29) result++;
    if(trackermuonquality & 1<<30) result++;
    if(trackermuonquality & 1<<31) result++;
    return(result);
}

Electron::Electron(const Analyse *ma, UInt_t n, int correction) :
    TLorentzVector(ma->electron_px->at(n), ma->electron_py->at(n), ma->electron_pz->at(n), sqrt(ma->electron_px->at(n)*ma->electron_px->at(n)+ma->electron_py->at(n)*ma->electron_py->at(n)+ma->electron_pz->at(n)*ma->electron_pz->at(n)+ElectronMassQ)),
    TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTElectronNames(), ma->electron_trigger->at(n)),
    //outerpoint(ma->electron_outerx->at(n), ma->electron_outery->at(n), ma->electron_outerz->at(n)),
    //closestpoint(ma->electron_closestpointx->at(n), ma->electron_closestpointy->at(n), ma->electron_closestpointz->at(n)),
    dxy(ma->electron_dxy->at(n)),
    dxyerr(ma->electron_dxyerr->at(n)),
    dz(ma->electron_dz->at(n)),
    dzerr(ma->electron_dzerr->at(n)),
    correctedecalenergy(ma->electron_correctedecalenergy->at(n)),
    esuperclusterovertrack(ma->electron_esuperclusterovertrack->at(n)),
    eseedclusterovertrack(ma->electron_eseedclusterovertrack->at(n)),
    deltaetasuperclustertrack(ma->electron_deltaetasuperclustertrack->at(n)),
    deltaphisuperclustertrack(ma->electron_deltaphisuperclustertrack->at(n)),
    e1x5(ma->electron_e1x5->at(n)),
    e2x5(ma->electron_e2x5->at(n)),
    e5x5(ma->electron_e5x5->at(n)),
    r9(ma->electron_r9->at(n)),
    sigmaetaeta(ma->electron_sigmaetaeta->at(n)),
    sigmaietaieta(ma->electron_sigmaietaieta->at(n)),
    sigmaiphiiphi(ma->electron_sigmaiphiiphi->at(n)),
    ehcaloverecaldepth1(ma->electron_ehcaloverecaldepth1->at(n)),
    ehcaloverecaldepth2(ma->electron_ehcaloverecaldepth2->at(n)),
    ehcaltoweroverecaldepth1(ma->electron_ehcaltoweroverecaldepth1->at(n)),
    ehcaltoweroverecaldepth2(ma->electron_ehcaltoweroverecaldepth2->at(n)),
    isolationr3track(ma->electron_isolationr3track->at(n)),
    isolationr3ecal(ma->electron_isolationr3ecal->at(n)),
    isolationr3hcal(ma->electron_isolationr3hcal->at(n)),
    isolationr4track(ma->electron_isolationr4track->at(n)),
    isolationr4ecal(ma->electron_isolationr4ecal->at(n)),
    isolationr4hcal(ma->electron_isolationr4hcal->at(n)),
    isolationpfr3charged(ma->electron_isolationpfr3charged->at(n)),
    isolationpfr3photon(ma->electron_isolationpfr3photon->at(n)),
    isolationpfr3neutral(ma->electron_isolationpfr3neutral->at(n)),
    trackchi2(ma->electron_trackchi2->at(n)),
    trackndof(ma->electron_trackndof->at(n)),
    nhits(ma->electron_nhits->at(n)),
    nmissinghits(ma->electron_nmissinghits->at(n)),
    npixelhits(ma->electron_npixelhits->at(n)),
    npixellayers(ma->electron_npixellayers->at(n)),
    nstriplayers(ma->electron_nstriplayers->at(n)),
    nhitsexpected(ma->electron_nhitsexpected->at(n)),
    //convdist(ma->electron_convdist->at(n)),
    //convdcot(ma->electron_convdcot->at(n)),
    //convradius(ma->electron_convradius->at(n)),
    iselectron(ma->electron_iselectron->at(n)),
    passconversionveto(ma->electron_passconversionveto->at(n)),
    ecaldrivenseed(ma->electron_ecaldrivenseed->at(n)),
    trackerdrivenseed(ma->electron_trackerdrivenseed->at(n)),
    gapinfo(ma->electron_gapinfo->at(n)),
    fbrems(ma->electron_fbrems->at(n)),
    numbrems(ma->electron_numbrems->at(n)),
    charge(ma->electron_charge->at(n)) {
    //info(ma->electron_info->at(n)) {
    //vtx(ma->electron_vtx->at(n)) {
    //supercluster.push_back(SuperCluster(ma->electron_supercluster_e->at(n), ma->electron_supercluster_x->at(n), ma->electron_supercluster_y->at(n), ma->electron_supercluster_z->at(n), ma->electron_supercluster_rawe->at(n), ma->electron_supercluster_phiwidth->at(n), ma->electron_supercluster_etawidth->at(n)));

    if(correction == 2012) {
        int cat = -1;
        if(IsEB() && Abs(Eta()) < 1. && R9() < 0.94) cat = 0;
        if(IsEB() && Abs(Eta()) < 1. && R9() >= 0.94) cat = 1;
        if(IsEB() && Abs(Eta()) >= 1. && R9() < 0.94) cat = 2;
        if(IsEB() && Abs(Eta()) >= 1. && R9() >= 0.94) cat = 3;
        if(!IsEB() && Abs(Eta()) < 1. && R9() < 0.94) cat = 4;
        if(!IsEB() && Abs(Eta()) < 1. && R9() >= 0.94) cat = 5;
        if(!IsEB() && Abs(Eta()) >= 1. && R9() < 0.94) cat = 6;
        if(!IsEB() && Abs(Eta()) >= 1. && R9() >= 0.94) cat = 7;
        double momscale = 1.;
        if(ma->IsData()) {
            double scaleval[][10] = {
                {190645,190781,1.0020,0.9980,1.0032,0.9919,0.9945,0.9881,0.9965,0.9862},
                {190782,191042,1.0079,1.0039,1.0063,0.9951,0.9996,0.9932,1.0010,0.9907},
                {191043,193555,0.9989,0.9949,0.9998,0.9885,0.9968,0.9904,0.9987,0.9884},
                {193556,194150,0.9974,0.9934,0.9954,0.9841,0.9969,0.9905,0.9988,0.9885},
                {194151,194532,0.9980,0.9940,0.9965,0.9852,0.9986,0.9922,0.9994,0.9891},
                {194533,195113,0.9983,0.9943,0.9984,0.9872,1.0006,0.9943,0.9999,0.9896},
                {195114,195915,0.9984,0.9944,0.9977,0.9864,1.0010,0.9946,1.0004,0.9900},
                {195916,198115,0.9975,0.9936,0.9965,0.9852,1.0020,0.9956,0.9992,0.9889},
                {198116,199803,1.0010,0.9970,0.9999,0.9886,0.9963,0.9899,1.0044,0.9941},
                {199804,200048,1.0021,0.9982,1.0008,0.9895,0.9965,0.9901,1.0060,0.9957},
                {200049,200151,1.0035,0.9996,1.0017,0.9905,0.9992,0.9928,1.0101,0.9999},
                {200152,200490,1.0013,0.9973,1.0003,0.9890,0.9991,0.9927,1.0073,0.9970},
                {200491,200531,1.0035,0.9995,1.0017,0.9905,0.9995,0.9931,1.0106,1.0004},
                {200532,201656,1.0017,0.9978,0.9999,0.9887,0.9978,0.9914,1.0069,0.9967},
                {201657,202305,1.0026,0.9986,1.0003,0.9891,0.9987,0.9923,1.0121,1.0018},
                {202305,203002,1.0037,0.9998,1.0010,0.9897,1.0003,0.9940,1.0144,1.0042},
                {203003,203984,1.0061,1.0024,1.0021,0.9923,0.9994,0.9914,1.0122,1.0033},
                {203985,205085,1.0050,1.0012,1.0045,0.9947,1.0004,0.9924,1.0108,1.0018},
                {205086,205310,1.0062,1.0025,1.0045,0.9947,1.0060,0.9981,1.0208,1.0119},
                {205311,206207,1.0056,1.0018,1.0033,0.9935,1.0015,0.9936,1.0157,1.0068},
                {206208,206483,1.0060,1.0022,1.0036,0.9938,0.9993,0.9913,1.0181,1.0092},
                {206484,206597,1.0062,1.0025,1.0033,0.9935,1.0043,0.9964,1.0172,1.0082},
                {206598,206896,1.0060,1.0023,1.0021,0.9923,0.9980,0.9900,1.0169,1.0079},
                {206897,207220,1.0063,1.0025,1.0033,0.9935,0.9985,0.9905,1.0180,1.0091},
                {207221,208686,1.0064,1.0026,1.0036,0.9938,1.0027,0.9948,1.0190,1.0100}
            };

            for(int i = 0 ; i < 25 ; ++i) {
                if(scaleval[i][0] <= ma->Run() && scaleval[i][1] > ma->Run()) {
                    momscale = scaleval[i][cat+2];
                    break;
                }

            }
        } else {
            double range = ma->analysisrandom->Uniform();
            double dsigMC = 0.;
            if(range < 0.63) {
                if(cat == 0) dsigMC = 0.0103;
                else if(cat == 1) dsigMC = 0.0090;
                else if(cat == 2) dsigMC = 0.0190;
                else if(cat == 3) dsigMC = 0.0156;
                else if(cat == 4) dsigMC = 0.0269;
                else if(cat == 5) dsigMC = 0.0287;
                else if(cat == 6) dsigMC = 0.0364;
                else if(cat == 7) dsigMC = 0.0321;
            } else {
                if(cat == 0) dsigMC = 0.0109;
                if(cat == 1) dsigMC = 0.0099;
                if(cat == 2) dsigMC = 0.0182;
                if(cat == 3) dsigMC = 0.0200;
                if(cat == 4) dsigMC = 0.0282;
                if(cat == 5) dsigMC = 0.0309;
                if(cat == 6) dsigMC = 0.0386;
                if(cat == 7) dsigMC = 0.0359;
            }
            momscale = ma->analysisrandom->Gaus(1., dsigMC);
        }
        (*this) *= momscale;
    }
}

/*
bool Electron::WorkingPoint(Int_t Missinghits, Float_t Convdist, Float_t Convdcot, Float_t Sigmaietaieta, Float_t Deltaphisctrack, Float_t Deltaetasctrack, Float_t Ehcaloverecal) const {
    return(NMissingHits() <= Missinghits && !(Abs(ConversionDist()) < Convdist && Abs(ConversionDCot()) < Convdcot) && SigmaIEtaIEta() < Sigmaietaieta && DeltaEtaSuperClusterTrack() < Deltaetasctrack && DeltaPhiSuperClusterTrack() < Deltaphisctrack && EHcalOverECal() < Ehcaloverecal);
}

bool Electron::WP95_v1(Int_t combined) const {
    if(combined == 1) {
        if(IsEE()) {
            return(WorkingPoint(1, -1., -1., 0.03, 0.7, 0.01, 0.07) && IsoR3CombinedRel() < 0.1);
        } else if(IsEB()) {
            return(WorkingPoint(1, -1., -1., 0.01, 0.8, 0.007, 0.15) && IsoR3CombinedRel() < 0.15);
        }
    } else {
        if(IsEE()) {
            return(WorkingPoint(1, -1., -1., 0.03, 0.7, 0.01, 0.07) && IsoR3TrackRel() < 0.08 && IsoR3ECalRel() < 0.06 && IsoR3HCalRel() < 0.05);
        } else if(IsEB()) {
            return(WorkingPoint(1, -1., -1., 0.01, 0.8, 0.007, 0.15) && IsoR3TrackRel() < 0.15 && IsoR3ECalRel() < 2. && IsoR3HCalRel() < 0.12);
        }
    }
    return(false);
}

bool Electron::WP90_v1(Int_t combined) const {
    if(combined == 1) {
        if(IsEE()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.7, 0.009, 0.05) && IsoR3CombinedRel() < 0.07);
        } else if(IsEB()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.8, 0.007, 0.12) && IsoR3CombinedRel() < 0.1);
        }
    } else {
        if(IsEE()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.7, 0.009, 0.05) && IsoR3TrackRel() < 0.05 && IsoR3ECalRel() < 0.06 && IsoR3HCalRel() < 0.03);
        } else if(IsEB()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.8, 0.007, 0.12) && IsoR3TrackRel() < 0.12 && IsoR3ECalRel() < 0.09 && IsoR3HCalRel() < 0.1);
        }
    }
    return(false);
}

bool Electron::WP85_v1(Int_t combined) const {
    if(combined == 1) {
        if(IsEE()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.04, 0.007, 0.025) && IsoR3CombinedRel() < 0.06);
        } else if(IsEB()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.06, 0.006, 0.04) && IsoR3CombinedRel() < 0.09);
        }
    } else {
        if(IsEE()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.03, 0.04, 0.007, 0.025) && IsoR3TrackRel() < 0.05 && IsoR3ECalRel() < 0.05 && IsoR3HCalRel() < 0.025);
        } else if(IsEB()) {
            return(WorkingPoint(1, 0.02, 0.02, 0.01, 0.06, 0.006, 0.04) && IsoR3TrackRel() < 0.09 && IsoR3ECalRel() < 0.08 && IsoR3HCalRel() < 0.1);
        }
    }
    return(false);
}

bool Electron::WP80_v1(Int_t combined) const {
    if(combined == 1) {
        if(IsEE()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.03, 0.007, 0.025) && IsoR3CombinedRel() < 0.06);
        } else if(IsEB()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.06, 0.004, 0.04) && IsoR3CombinedRel() < 0.07);
        }
    } else {
        if(IsEE()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.03, 0.007, 0.025) && IsoR3TrackRel() < 0.04 && IsoR3ECalRel() < 0.05 && IsoR3HCalRel() < 0.025);
        } else if(IsEB()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.06, 0.004, 0.04) && IsoR3TrackRel() < 0.09 && IsoR3ECalRel() < 0.07 && IsoR3HCalRel() < 0.1);
        }
    }
    return(false);
}

bool Electron::WP70_v1(Int_t combined) const {
    if(combined == 1) {
        if(IsEE()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3CombinedRel() < 0.03);
        } else if(IsEB()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.03, 0.004, 0.025) && IsoR3CombinedRel() < 0.04);
        }
    } else {
        if(IsEE()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3TrackRel() < 0.025 && IsoR3ECalRel() < 0.025 && IsoR3HCalRel() < 0.02);
        } else if(IsEB()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.03, 0.004, 0.025) && IsoR3TrackRel() < 0.05 && IsoR3ECalRel() < 0.06 && IsoR3HCalRel() < 0.03);
        }
    }
    return(false);
}

bool Electron::WP60_v1(Int_t combined) const {
    if(combined == 1) {
        if(IsEE()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3CombinedRel() < 0.02);
        } else if(IsEB()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.025, 0.004, 0.025) && IsoR3CombinedRel() < 0.03);
        }
    } else {
        if(IsEE()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.03, 0.02, 0.005, 0.025) && IsoR3TrackRel() < 0.025 && IsoR3ECalRel() < 0.02 && IsoR3HCalRel() < 0.02);
        } else if(IsEB()) {
            return(WorkingPoint(0, 0.02, 0.02, 0.01, 0.025, 0.004, 0.025) && IsoR3TrackRel() < 0.04 && IsoR3ECalRel() < 0.04 && IsoR3HCalRel() < 0.03);
        }
    }
    return(false);
}
*/
//bool Electron::Loose2012v1()
//{
//	if(IsEE())
//	{
//if(DeltaEtaSuperClusterTrack() < 0.007 && DeltaPhiSuperClusterTrack() < 0.8 && SigmaIEtaIEta() < 0.01 && EHcalOverECal() < 0.15 && Dxy() < 0.04 && Dz() < 0.2 && Abs(1.- ESuperClusterOverTrack())/SCs(0).E() < 1000000. )
//	}
//	else if(IsEB())
//	{
//	}
//	return(false);
//}

//Conversion::Conversion(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t *Cov, Float_t mymvaout, TVector3 Ecalpoint1, TVector3 Ecalpoint2) : TrackComposedParticle(Vx, Vy, Vz, Chi2, Ndof, Cov),
//    mvaout(mymvaout) {
//    ecalpoints.push_back(Ecalpoint1);
//    ecalpoints.push_back(Ecalpoint2);
//}

Photon::Photon(const Analyse *ma, UInt_t n) :
    TLorentzVector(ma->photon_px->at(n), ma->photon_py->at(n), ma->photon_pz->at(n), sqrt(ma->photon_px->at(n)*ma->photon_px->at(n)+ma->photon_py->at(n)*ma->photon_py->at(n)+ma->photon_pz->at(n)*ma->photon_pz->at(n))),
    TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTPhotonNames(), 0),
    e1x5(ma->photon_e1x5->at(n)),
    e2x5(ma->photon_e2x5->at(n)),
    e3x3(ma->photon_e3x3->at(n)),
    e5x5(ma->photon_e5x5->at(n)),
    sigmaietaieta(ma->photon_sigmaietaieta->at(n)),
    sigmaietaiphi(ma->photon_sigmaietaiphi->at(n)),
    sigmaiphiiphi(ma->photon_sigmaiphiiphi->at(n)),
    ehcaloverecaldepth1(ma->photon_ehcaloverecaldepth1->at(n)),
    ehcaloverecaldepth2(ma->photon_ehcaloverecaldepth2->at(n)),
    ehcaltoweroverecaldepth1(ma->photon_ehcaltoweroverecaldepth1->at(n)),
    ehcaltoweroverecaldepth2(ma->photon_ehcaltoweroverecaldepth2->at(n)),
    maxenergyxtal(ma->photon_maxenergyxtal->at(n)),
    isolationr3track(ma->photon_isolationr3track->at(n)),
    isolationr3trackhollow(ma->photon_isolationr3trackhollow->at(n)),
    isolationr3ntrack(ma->photon_isolationr3ntrack->at(n)),
    isolationr3ntrackhollow(ma->photon_isolationr3ntrackhollow->at(n)),
    isolationr3ecal(ma->photon_isolationr3ecal->at(n)),
    isolationr3hcal(ma->photon_isolationr3hcal->at(n)),
    isolationr4track(ma->photon_isolationr4track->at(n)),
    isolationr4trackhollow(ma->photon_isolationr4trackhollow->at(n)),
    isolationr4ntrack(ma->photon_isolationr4ntrack->at(n)),
    isolationr4ntrackhollow(ma->photon_isolationr4ntrackhollow->at(n)),
    isolationr4ecal(ma->photon_isolationr4ecal->at(n)),
    isolationr4hcal(ma->photon_isolationr4hcal->at(n)),
    isolationpfr3charged(ma->photon_isolationpfr3charged->at(n)),
    isolationpfr3photon(ma->photon_isolationpfr3photon->at(n)),
    isolationpfr3neutral(ma->photon_isolationpfr3neutral->at(n)),
    //isolationpfr4charged(ma->photon_isolationpfr4charged->at(n)),
    //isolationpfr4photon(ma->photon_isolationpfr4photon->at(n)),
    //isolationpfr4neutral(ma->photon_isolationpfr4neutral->at(n)),
    //isolationpfr4noscfootprintcharged(ma->photon_isolationpfr4noscfootprintcharged->at(n)),
    //isolationpfr4noscfootprintphoton(ma->photon_isolationpfr4noscfootprintphoton->at(n)),
    //isolationpfr4noscfootprintneutral(ma->photon_isolationpfr4noscfootprintneutral->at(n)),
    //info(ma->photon_info->at(n)),
    isphoton(ma->photon_isphoton->at(n)),
    hasconversiontracks(ma->photon_hasconversiontracks->at(n)),
    haspixelseed(ma->photon_haspixelseed->at(n)),
    passelectronveto(ma->photon_passelectronveto->at(n)),
    ispfphoton(ma->photon_ispfphoton->at(n)),
    gapinfo(ma->photon_gapinfo->at(n)) {
    //supercluster.push_back(SuperCluster(ma->photon_supercluster_e->at(n), ma->photon_supercluster_x->at(n), ma->photon_supercluster_y->at(n), ma->photon_supercluster_z->at(n), ma->photon_supercluster_rawe->at(n), ma->photon_supercluster_phiwidth->at(n), ma->photon_supercluster_etawidth->at(n)));
}

Jet::Jet(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t Hadronicenergy, Float_t Chargedhadronicenergy, Float_t Emenergy, Float_t Chargedemenergy, Float_t Hfemenergy, Float_t Hfhadronicenergy, Float_t Electronenergy, Float_t Muonenergy, Int_t Chargedmulti, Int_t Neutralmulti, Int_t Hfemmulti, Int_t Hfhadronicmulti, Int_t Electronmulti, Int_t Muonmulti, Float_t Chargeda, Float_t Chargedb, Float_t Neutrala, Float_t Neutralb, Float_t Alla, Float_t Allb, Float_t Chargedfractionmv, Float_t Energycorr,Float_t Energycorrunc, Float_t Energycorrl7uds, Float_t Energycorrl7bottom, Int_t Btag, Int_t Mcflavour, Float_t Puidfull, Float_t Puidsimple, Float_t Puidcutbased) :
    TLorentzVector(Px, Py, Pz, E),
    hadronicenergy(Hadronicenergy),
    chargedhadronicenergy(Chargedhadronicenergy),
    emenergy(Emenergy),
    chargedemenergy(Chargedemenergy),
    hfemenergy(Hfemenergy),
    hfhadronicenergy(Hfhadronicenergy),
    electronenergy(Electronenergy),
    muonenergy(Muonenergy),
    chargedmulti(Chargedmulti),
    neutralmulti(Neutralmulti),
    hfemmulti(Hfemmulti),
    hfhadronicmulti(Hfhadronicmulti),
    electronmulti(Electronmulti),
    muonmulti(Muonmulti),
    chargeda(Chargeda),
    chargedb(Chargedb),
    neutrala(Neutrala),
    neutralb(Neutralb),
    alla(Alla),
    allb(Allb),
    chargedfractionmv(Chargedfractionmv),
    energycorr(Energycorr),
    energycorrunc(Energycorrunc),
    energycorrl7uds(Energycorrl7uds),
    energycorrl7bottom(Energycorrl7bottom),
    mcflavour(Mcflavour),
    btag(Btag),
    puidfull(Puidfull),
    puidsimple(Puidsimple),
    puidcutbased(Puidcutbased),
    mymvonly(false) {
}

void Jet::ScaleMV(bool mvonly) {
    if(mymvonly == false && mvonly == true) {
        Double_t f = ChargedMomentumFractionFromMV();
        this->SetPxPyPzE(f*Px(), f*Py(), f*Pz(), f*E());
        mymvonly = true;
    }
    if(mymvonly == true && mvonly == false) {
        Double_t f = ChargedMomentumFractionFromMV();
        this->SetPxPyPzE(Px()/f, Py()/f, Pz()/f, E()/f);
        mymvonly = false;
    }
}

bool Jet::BTag(string disc) const {
    if     (disc == "JPL")     return (btag & 1);
    else if(disc == "JPM")     return (btag & 2);
    else if(disc == "JPT")     return (btag & 3);
    else if(disc == "CSVv2L")  return (btag & 4);
    else if(disc == "CSVv2M")  return (btag & 5);
    else if(disc == "CSVv2T")  return (btag & 6);
    else if(disc == "CMVAv2L") return (btag & 7);
    else if(disc == "CMVAv2M") return (btag & 8);
    else if(disc == "CMVAv2T") return (btag & 9);
    return(false);
}


Track::Track(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Double_t Ox, Double_t Oy, Double_t Oz, Double_t Cx, Double_t Cy, Double_t Cz, Double_t Chi2, Double_t Ndof, Double_t Dxy, Double_t Dxyerr, Double_t Dz, Double_t Dzerr, Int_t Charge, Int_t Nhits, Int_t Nmissinghits, Int_t Npixelhits, Int_t Npixellayers, Int_t Nstriplayers, Int_t Vtx, Float_t Dedxharmonic2) :
    TLorentzVector(Px, Py, Pz, E),
    outerpoint(Ox, Oy, Oz),
    closestpoint(Cx, Cy, Cz),
    chi2(Chi2),
    ndof(Ndof),
    dxy(Dxy),
    dxyerr(Dxyerr),
    dz(Dz),
    dzerr(Dzerr),
    charge(Charge),
    nhits(Nhits),
    nmissinghits(Nmissinghits),
    npixelhits(Npixelhits),
    npixellayers(Npixellayers),
    nstriplayers(Nstriplayers),
    //vtx(Vtx),
    dedxharmonic2(Dedxharmonic2) {
}

EcalHit::EcalHit(Double_t Energy, Double_t X, Double_t Y, Double_t Z) : TVector3(X, Y, Z), energy(Energy) {
}

Cluster::Cluster(Double_t E, Double_t x, Double_t y, Double_t z, Int_t Size) :
    TLorentzVector(E *x/sqrt(x*x+y*y+z *z), E *y/sqrt(x*x+y*y+z *z), E *z/sqrt(x*x+y*y+z *z), E),
    position(x,y,z),
    size(Size) {
}

void Cluster::AddHit(EcalHit &hit) {
    hits.push_back(hit);
}
/*
SuperCluster::SuperCluster(Double_t E, Double_t x, Double_t y, Double_t z, Float_t Rawe, Float_t Phiwidth, Float_t Etawidth) :
    TLorentzVector(E *x/sqrt(x*x+y*y+z *z), E *y/sqrt(x*x+y*y+z *z), E *z/sqrt(x*x+y*y+z *z), E),
    position(x,y,z),
    rawe(Rawe),
    phiwidth(Phiwidth),
    etawidth(Etawidth) {
}

void SuperCluster::AddCluster(Cluster &newcluster) {
    clusters.push_back(newcluster);
}

void SuperCluster::AddESCluster(Cluster &newcluster) {
    esclusters.push_back(newcluster);
}
*/
std::vector<Float_t> dummyvec;
Tau::Tau(const Analyse *ma, UInt_t n) :
    TLorentzVector(ma->tau_px->at(n), ma->tau_py->at(n), ma->tau_pz->at(n), Sqrt(ma->tau_px->at(n)*ma->tau_px->at(n)+ma->tau_py->at(n)*ma->tau_py->at(n)+ma->tau_pz->at(n)*ma->tau_pz->at(n)+TauMassQ)),
    TriggerObject(ma, ma->runlist.find(ma->Run())->second.GetHLTTauNames(), ma->tau_trigger->at(n)),
    isolationneutralspt(ma->tau_isolationneutralspt->at(n)),
    isolationneutralsnum(ma->tau_isolationneutralsnum->at(n)),
    isolationchargedpt(ma->tau_isolationchargedpt->at(n)),
    isolationchargednum(ma->tau_isolationchargednum->at(n)),
    isolationgammapt(ma->tau_isolationgammapt->at(n)),
    isolationgammanum(ma->tau_isolationgammanum->at(n)),
    charge(ma->tau_charge->at(n)),
    dishps(ma->tau_dishps->at(n)),
    emfraction(ma->tau_emfraction->at(n)),
    hcaltotoverplead(ma->tau_hcaltotoverplead->at(n)),
    hcal3x3overplead(ma->tau_hcal3x3overplead->at(n)),
    ecalstripsumeoverplead(ma->tau_ecalstripsumeoverplead->at(n)),
    bremsrecoveryeoverplead(ma->tau_bremsrecoveryeoverplead->at(n)),
    calocomp(ma->tau_calocomp->at(n)),
    segcomp(ma->tau_segcomp->at(n)),
    jet(ma->tau_ak4pfjet_e->at(n), ma->tau_ak4pfjet_px->at(n), ma->tau_ak4pfjet_py->at(n), ma->tau_ak4pfjet_pz->at(n), ma->tau_ak4pfjet_hadronicenergy->at(n), ma->tau_ak4pfjet_chargedhadronicenergy->at(n), ma->tau_ak4pfjet_emenergy->at(n), ma->tau_ak4pfjet_chargedemenergy->at(n), -1.,-1.,-1.,-1., ma->tau_ak4pfjet_chargedmulti->at(n), ma->tau_ak4pfjet_neutralmulti->at(n),-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1., -1., -1., -1, 0, 0, 0., 0., 0.),
    //jet(ma->tau_ak4pfjet_e->at(n), ma->tau_ak4pfjet_px->at(n), ma->tau_ak4pfjet_py->at(n), ma->tau_ak4pfjet_pz->at(n), ma->tau_ak4pfjet_hadronicenergy->at(n), ma->tau_ak4pfjet_chargedhadronicenergy->at(n), ma->tau_ak4pfjet_emenergy->at(n), ma->tau_ak4pfjet_chargedemenergy->at(n), -1.,-1.,-1.,-1., ma->tau_ak4pfjet_chargedmulti->at(n), ma->tau_ak4pfjet_neutralmulti->at(n),-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1., -1., -1., -1, 0, 0., 0., 0.),
    taudiscriminators(ma->runlist.find(ma->Run())->second.GetTauDiscriminators()) {
    UInt_t begin = ma->tau_chargedbegin->at(n);
    UInt_t end = ma->tau_charged_count;
    if(n < ma->tau_count - 1) {
        end = ma->tau_chargedbegin->at(n+1);
    }

    for(UInt_t i = begin ; i < end ; i++) {
        tracks.push_back(Track(Sqrt(ma->tau_charged_px->at(i)*ma->tau_charged_px->at(i)+ma->tau_charged_py->at(i)*ma->tau_charged_py->at(i)+ma->tau_charged_pz->at(i)*ma->tau_charged_pz->at(i)), ma->tau_charged_px->at(i), ma->tau_charged_py->at(i), ma->tau_charged_pz->at(i), ma->tau_charged_outerx->at(i), ma->tau_charged_outery->at(i), ma->tau_charged_outerz->at(i), ma->tau_charged_closestpointx->at(i), ma->tau_charged_closestpointy->at(i), ma->tau_charged_closestpointz->at(i), ma->tau_charged_chi2->at(i), ma->tau_charged_ndof->at(i), ma->tau_charged_dxy->at(i), ma->tau_charged_dxyerr->at(i), ma->tau_charged_dz->at(i), ma->tau_charged_dzerr->at(i), ma->tau_charged_charge->at(i), ma->tau_charged_nhits->at(i), ma->tau_charged_nmissinghits->at(i), ma->tau_charged_npixelhits->at(i), ma->tau_charged_npixellayers->at(i), ma->tau_charged_nstriplayers->at(i), -1, ma->tau_charged_dedxharmonic2->at(i)));
    }

}

Track Tau::LeadingTrack() const {
    Float_t ptmax = 0.;
    Int_t index = -1;
    for(UInt_t n = 0 ; n < tracks.size() ; n++) {
        if(ptmax < tracks[n].Pt()) {
            ptmax = tracks[n].Pt();
            index = n;
        }
    }
    return(tracks[index]);
}

Int_t Tau::TauDiscriminator(string disname) const {
    Int_t pos = -1;
    for(UInt_t i = 0 ; i < taudiscriminators->size() ; i++) {
        if(disname == taudiscriminators->at(i)) {
            pos = i;
            break;
        }
    }
    if(pos == -1) return(-1);
    if((dishps & 1<<pos) != 0) return(1);
    return(0);
}

Vertex::Vertex(Double_t X, Double_t Y, Double_t Z, UInt_t Isvalid, UInt_t Isfake, Double_t Chi2, Double_t Ndof, Int_t Ntracks, Float_t Cov0, Float_t Cov1, Float_t Cov2, Float_t Cov3, Float_t Cov4, Float_t Cov5) :
    TVector3(X, Y, Z),
    chi2(Chi2),
    ndof(Ndof),
    ntracks(Ntracks),
    isvalid(Isvalid),
    isfake(Isfake),
    covmatrix(3) {
    covmatrix(0,0) = Cov0;
    covmatrix(0,1) = Cov1;
    covmatrix(0,2) = Cov2;
    covmatrix(1,0) = Cov1;
    covmatrix(1,1) = Cov3;
    covmatrix(1,2) = Cov4;
    covmatrix(2,0) = Cov2;
    covmatrix(2,1) = Cov4;
    covmatrix(2,2) = Cov5;
}

BeamSpot::BeamSpot(Double_t X, Double_t Y, Double_t Z, const Float_t *Cov, Double_t Xwidth, Double_t Ywidth, Double_t Zsigma) :
    TVector3(X, Y, Z),
    covmatrix(3),
    xwidth(Xwidth),
    ywidth(Ywidth),
    zsigma(Zsigma) {
    covmatrix(0,0) = Cov[0];
    covmatrix(0,1) = Cov[1];
    covmatrix(0,2) = Cov[2];
    covmatrix(1,0) = Cov[1];
    covmatrix(1,1) = Cov[3];
    covmatrix(1,2) = Cov[4];
    covmatrix(2,0) = Cov[2];
    covmatrix(2,1) = Cov[4];
    covmatrix(2,2) = Cov[5];
}

void splitstring(string input, vector<string> &output) {
    UInt_t posstart = 0;
    for(UInt_t i = 0 ; i < input.size() ; i++) {
        if(input[i] == ' ') {
            output.push_back(input.substr(posstart, i-posstart));
            posstart=i+1;
        }
    }
}

string combinestring(const vector<string> &input) {
    string output;
    for(UInt_t i = 0 ; i < input.size() ; i++) {
        output += input[i] + string(" ");
    }
    return(output);
}

//Luminosity
void Luminosity::AddFiles(const vector<string> &newfiles) {
    for(size_t i = 0 ; i < newfiles.size() ; i++) {
        if(find(filenames.begin(), filenames.end(), newfiles[i]) == filenames.end()) {
            filenames.push_back(newfiles[i]);
        }
    }
}

string Luminosity::GetFilesString() const {
    string files(filenames[0].substr(filenames[0].find_last_of("/")+1));
    for(size_t i = 1 ; i < filenames.size() ; i++) {
        files += " " + filenames[i].substr(filenames[i].find_last_of("/")+1);
    }
    return(files);
}

//RunInfo
RunInfo::RunInfo(UInt_t Number, UInt_t Hltcount, string Hltnames, string Hltmunames, string Hltelnames, string Hltphotonnames, string Hlttaunames, string Hltjetnames, string Taudiscriminators, UInt_t Hltprescaletablescount, UInt_t *Hltprescaletables, UInt_t L1algocount, UInt_t L1algoprescaletablescount, UInt_t *L1algoprescaletables, UInt_t L1techcount, UInt_t L1techprescaletablescount, UInt_t *L1techprescaletables) :
    runnumber(Number),
    hltnum(Hltcount),
    hlttablesnum(Hltprescaletablescount/Hltcount),
    l1algonum(L1algocount),
    l1algotablesnum(L1algoprescaletablescount/L1algocount),
    l1technum(L1techcount),
    l1techtablesnum(L1techprescaletablescount/L1techcount) {
    splitstring(Hltnames, hltnames);
    splitstring(Hltmunames, hltmunames);
    splitstring(Hltelnames, hltelnames);
    splitstring(Hlttaunames, hlttaunames);
    splitstring(Hltphotonnames, hltphotonnames);
    splitstring(Hltjetnames, hltjetnames);
    splitstring(Taudiscriminators, taudiscriminators);
    hltprescales.resize(NumHLT());
    for(UInt_t i = 0 ; i < NumHLT() ; i++) {
        for(UInt_t j = 0 ; j < NumHLTTables() ; j++) {
            hltprescales[i].push_back(Hltprescaletables[i+Hltcount*j]);
        }
    }

    l1algoprescales.resize(NumL1Algo());
    for(UInt_t i = 0 ; i < NumL1Algo() ; i++) {
        for(UInt_t j = 0 ; j < NumL1AlgoTables() ; j++) {
            l1algoprescales[i].push_back(L1algoprescaletables[i+L1algocount*j]);
        }
    }

    l1techprescales.resize(NumL1Tech());
    for(UInt_t i = 0 ; i < NumL1Tech() ; i++) {
        for(UInt_t j = 0 ; j < NumL1TechTables() ; j++) {
            l1techprescales[i].push_back(L1techprescaletables[i+L1techcount*j]);
        }
    }

}

vector<string> RunInfo::MatchTriggerNames(string name) {
    boost::cmatch what;
    vector<string> result;

    boost::regex trigregex(name.c_str());
    for(UInt_t i = 0 ; i < NumHLT() ; i++) {
        if(boost::regex_match(HLTName(i).c_str(), what, trigregex)) {
            result.push_back(HLTName(i));
        }
    }
    return(result);
}

Int_t RunInfo::HLTIndex(string name) const {
    for(Int_t i = 0 ; i < Int_t(hltnames.size()) ; i++) {
        if(hltnames[i] == name) return(i);
    }
    return(-1);
}

string RunInfo::HLTAllNames() const {
    return(combinestring(hltnames));
}

string RunInfo::HLTMuonAllNames() const {
    return(combinestring(hltmunames));
}

string RunInfo::HLTElectronAllNames() const {
    return(combinestring(hltelnames));
}

string RunInfo::HLTTauAllNames() const {
    return(combinestring(hlttaunames));
}

string RunInfo::HLTPhotonAllNames() const {
    return(combinestring(hltphotonnames));
}

string RunInfo::HLTJetAllNames() const {
    return(combinestring(hltjetnames));
}

string RunInfo::TauDiscriminatorsAllNames() const {
    return(combinestring(taudiscriminators));
}

void TriggerRun::SetBlock(UInt_t block, Int_t index, UInt_t prescale) {
    if(lumiinfo.find(block) == lumiinfo.end()) {
        lumiinfo[block] = TriggerLumi(index, prescale);
    }
}

TriggerLumi TriggerRun::GetBlock(UInt_t block) {
    //return(lumiinfo[block]);
    map<UInt_t, TriggerLumi>::const_iterator blockinfo = lumiinfo.find(block);
    if(blockinfo == lumiinfo.end()) {
        return(TriggerLumi(-1, 0));
    }
    return(blockinfo->second);
}

TriggerSelection::TriggerSelection(Analyse *an, vector<string> names, bool useprescaled) : AN(an) {
    for(map<UInt_t, RunInfo>::iterator a = AN->runlist.begin() ; a != AN->runlist.end(); ++a) {
        UInt_t runnumber = a->first;
        Float_t runlumi = 0.;
        vector<UInt_t> indices;
        for(UInt_t i = 0 ; i < names.size() ; i++) {
            vector<string> goodnames = a->second.MatchTriggerNames(names[i]);
            if(goodnames.size() == 1) {
                indices.push_back(a->second.HLTIndex(goodnames[0]));
            } else if(goodnames.size() > 1) {
                cout << "WARNING: TriggerSelection::TriggerSelection: " << names[i] << " is not a unique selection!" << endl;
            }
        }

        map<UInt_t, map<UInt_t, Luminosity> >::iterator blocklist = AN->lumilist.find(runnumber);
        for(map<UInt_t, Luminosity>::iterator b = blocklist->second.begin() ; b != blocklist->second.end() ; ++b) {
            Int_t minindex = -1;
            UInt_t minprescale = 10000000;

            for(UInt_t i = 0 ; i < indices.size() ; i++) {
                UInt_t prescale = a->second.HLTPrescale(indices[i], b->second.HLTTable());
                if(prescale <= 0) {
                    continue;
                }
                if(prescale < minprescale) {
                    minindex = indices[i];
                    minprescale = prescale;
                }
                if(prescale == 1) {
                    break;
                }
            }
            if(useprescaled == false && minprescale != 1) {
                runinfo[runnumber].SetBlock(b->first, -1, minprescale);
                continue;
            }
            if(b->second.LumiValue() != -1) {
                runlumi += b->second.LumiValue() * minprescale;
            }
            runinfo[runnumber].SetBlock(b->first, minindex, minprescale);
        }
        runinfo[runnumber].Lumi(runlumi);
    }
}

Int_t TriggerSelection::Result() {
    const TriggerLumi &triggerlumi = runinfo[AN->Run()].GetBlock(AN->LumiBlock());
    if(triggerlumi.Index() < 0) {
        return(0);
    }
    if(AN->GetHLTrigger(triggerlumi.Index())) {
        return(triggerlumi.Prescale());
    } else {
        return(triggerlumi.Prescale() * -1);
    }
}

Float_t TriggerSelection::LumiUsed(Int_t format) {
    Double_t lumi = 0.;
    Double_t alllumi = 0.;
    Double_t zerolumi = 0.;
    UInt_t nolumiinfo = 0;
    Double_t events = 0.;
    Double_t eventsprocessed = 0.;
    for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = AN->lumilist.begin() ; a != AN->lumilist.end() ; ++a) {
        Double_t runlumi = 0.;
        Double_t runalllumi = 0.;
        Double_t runzerolumi = 0.;
        Double_t runevents = 0.;
        Double_t runeventsprocessed = 0.;
        UInt_t runnolumiinfo = 0;
        UInt_t numblocks = 0;
        for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++) {
            const TriggerLumi &triggerlumi = runinfo[a->first].GetBlock(b->first);
            if(b->second.LumiValue() == -1 || triggerlumi.Index() < 0) {
                continue;
            }
            if(AN->IsInRange(a->first, b->first)) {
                numblocks++;
                //	cout << a->first << " " << b->first << " " << minLumi << " " << maxLumi << endl;
                if(b->second == true) {
                    if(b->second.NumEvents() > 0) {
                        runlumi += b->second.ProcessedLumi()*triggerlumi.Prescale();
                        runalllumi += b->second.LumiValue()*triggerlumi.Prescale();
                        runevents += b->second.NumEvents();
                        runeventsprocessed += b->second.NumEventsProcessed();
                        if(b->second.ProcessedFraction() > 1) {
                            cerr << "WARNING GetLumi: You ran on " << b->second.ProcessedFraction()*100 <<"% of available events in Run " << a->first << ", Lumiblock " << b->first <<". :-O" << endl;
                        }
                        if(format >= 3) cout << "    Block: " << b->first << ", fraction: " << b->second.ProcessedFraction() << ", lumi: " << b->second.ProcessedLumi() << " pb^-1" << endl;
                    } else {
                        runzerolumi += b->second.LumiValue()*triggerlumi.Prescale();
                        runalllumi += b->second.LumiValue()*triggerlumi.Prescale();
                    }
                } else if(b->second == false) {
                    runnolumiinfo += b->second.NumEventsProcessed();
                }
            }
        }
        if(a->first <= AN->maxRun && a->first >= AN->minRun) {
            if(runevents != 0) zerolumi *= runeventsprocessed/runevents;
            if(format >= 2)  cout << "Run: " << a->first << ", blocks: " << numblocks << ", fraction: " << runeventsprocessed/runevents << ", total lumi: " << runalllumi << " pb^-1, processed lumi: " << runlumi+zerolumi << " zero event lumi: " << zerolumi << " pb^-1, no lumi info for: " << runnolumiinfo << " event(s)."<< endl;
            lumi += runlumi;
            zerolumi += runzerolumi;
            alllumi += runalllumi;
            events += runevents;
            eventsprocessed += runeventsprocessed;
            nolumiinfo += runnolumiinfo;
        }
    }
    if(format >= 1) {
        cout << "From run " << AN->minRun << ", block " << AN->minLumi << " to run " << AN->maxRun << ", block " << AN->maxLumi << "." << endl;
        cout << "Total lumi: " << lumi+zerolumi << " pb^-1, zero event lumi: " << zerolumi << " pb^-1, (lumi in Range: " << alllumi << " pb^-1), no lumi-info for " << nolumiinfo << " event(s), events: " << eventsprocessed << "(" << events << ")"  << endl;
    }
    return(lumi+zerolumi);

}

Float_t TriggerSelection::LumiBeforeEvent() {
    Float_t lumi = 0.;
    for(map<UInt_t, TriggerRun>::iterator a = runinfo.begin() ; a != runinfo.end(); ++a) {
        if(a->first < AN->Run()) {
            lumi+=a->second.Lumi();
        } else {
            break;
        }
    }

    map<UInt_t, map<UInt_t, Luminosity> >::iterator blocklist = AN->lumilist.find(AN->Run());
    for(map<UInt_t, Luminosity>::iterator b = blocklist->second.begin() ; b != blocklist->second.end() ; ++b) {
        if(b->first <= AN->LumiBlock()) {
            const TriggerLumi &triggerlumi = runinfo[AN->Run()].GetBlock(b->first);
            if(triggerlumi.Index() >= 0 && b->second.LumiValue() != -1) {
                lumi+=b->second.LumiValue()*triggerlumi.Prescale();
            }
        } else {
            break;
        }
    }
    return(lumi);
}


string TriggerSelection::GetTriggerName(UInt_t run, UInt_t lumiblock) {
    if(run == 0 && lumiblock == 0) {
        run = AN->Run();
        lumiblock = AN->LumiBlock();
    }

    Int_t index = runinfo[run].GetBlock(lumiblock).Index();
    if(index != -1) {
        return(AN->runlist[run].HLTName(index));
    }
    cerr << "TriggerSelection::GetTriggerName: Invalid run and/or lumiblock." << endl;
    return("ERROR: TriggerSelection::GetTriggerName: Invalid run and/or lumiblock.");
}

void TriggerSelection::PrintInfo() {
    Float_t lumi = 0.;
    Float_t lumiuse = 0.;
    Float_t curprescale;
    string curname;
    bool beg = true;
    UInt_t begrun, begblock;
    Float_t beglumi = 0.;
    UInt_t run, block;

    for(map<UInt_t, TriggerRun>::iterator runit = runinfo.begin() ; runit != runinfo.end() ; ++runit) {
        for(map< UInt_t, TriggerLumi >::iterator lumiit = runit->second.Begin() ; lumiit != runit->second.End() ; ++lumiit) {
            run = runit->first;
            block = lumiit->first;
            Int_t index = lumiit->second.Index();
            Float_t prescale = 0;
            string name("NA");
            if(index != -1) {
                name = AN->runlist[run].HLTName(index);
                prescale = lumiit->second.Prescale();
                if(AN->lumilist[run][block].LumiValue() != -1) {
                    lumiuse += AN->lumilist[run][block].LumiValue();
                }
            }
            if(beg) {
                begrun = run;
                begblock = block;
                beglumi = lumi;
                curname = name;
                curprescale = prescale;
                beg = false;
            }
            if(AN->lumilist[run][block].LumiValue() != -1) {
                lumi += AN->lumilist[run][block].LumiValue();
            }
            if(curname != name || curprescale != prescale) {
                cout << begrun << ":" << begblock << " - " << run << ":" << block << ", " << curname << "(" << curprescale << "), Lumi: " << lumi << "(" << lumi - beglumi << ")/pb" << endl;
                begrun = run;
                begblock = block;
                beglumi = lumi;
                curname = name;
                curprescale = prescale;
            }
        }
    }
    cout << begrun << ":" << begblock << " - " << run << ":" << block << ", " << curname << "(" << curprescale << "), Lumi: " << lumi << "(" << lumi - beglumi << ")/pb" << endl;
    cout << "Triggered Lumi: " << lumiuse << "/pb" << endl;

}




