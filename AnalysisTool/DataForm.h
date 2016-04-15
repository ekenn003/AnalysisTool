#ifndef ACDATA
#define ACDATA

#include <TLorentzVector.h>
#include <TMatrixFSym.h>

#include <iostream>
#include <string>
#include <map>

class Analyse;
class Vertex;

using namespace std;
using namespace TMath;

const Int_t M_btagmax = 6;
const Double_t MuonMass = 0.105658;
const Double_t MuonMassQ = MuonMass*MuonMass;
const Double_t ElectronMass = 0.000511;
const Double_t ElectronMassQ = ElectronMass*ElectronMass;
const Double_t TauMass = 1.77699;
const Double_t TauMassQ = TauMass*TauMass;


//////////////////////////////////////////////////////////////////////
// CLASS: GENJET /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class GenJet : public TLorentzVector {
  private:
    Float_t einvisible;
    Int_t   flavour;
    UInt_t  info;
    bool    myvisible;
  public:
    // constructors
    GenJet(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t Einvisible, Int_t Flavour, UInt_t Info) : TLorentzVector(Px, Py, Pz, E), einvisible(Einvisible), flavour(Flavour), info(Info), myvisible(false) {}
    GenJet() : flavour(0) {}
    // methods
    Float_t InvisibleEnergy() const { return(einvisible); }
    Int_t   Flavour() const { return(flavour); }
    void    ScaleVisible(bool visible); //include or exlude invisible energy fraction
    enum    CONSTITUENTS { CONS_b=1<<0, CONS_bbar=1<<1, CONS_c=1<<2, CONS_cbar=1<<3, CONS_light=1<<4, CONS_gluon=1<<5 };
    bool    HasConstituent(UInt_t constituent) const { return(info & constituent); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: GENBASICPARTICLE ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class GenBasicParticle : public TLorentzVector {
  private:
    TVector3 vertex;
    Int_t    status;
    Int_t    pdgid;

  public:
    // constructors
    GenBasicParticle(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid);
    GenBasicParticle() : pdgid(0) {}
    // methods
    Int_t    Status() const { return(status); }
    Int_t    PDGId() const { return(pdgid); }
    TVector3 Vertex() const { return(vertex); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: GENLIGHTPARTICLE ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class GenLightParticle : public GenBasicParticle {
  private:
    UInt_t info;
    Int_t  indirectmother;

  public:
    // constructors
    GenLightParticle(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Info, Int_t Indirectmother);
    GenLightParticle() {}
    // methods
    UInt_t MoreInfo() const { return(info); }
    UInt_t MoreInfo(INFO test) const { return(info &test); }
    Int_t  IndirectMother() const { return(indirectmother); }
    enum   INFO {FromWp=1<<0, FromWm=1<<1, FromGamma=1<<2, FromZ=1<<3, Fromh=1<<4, FromH=1<<5, FromA=1<<6, FromHp=1<<7, FromHm=1<<8, Fromt=1<<9, Fromtbar=1<<10, Fromb=1<<11, Frombbar=1<<12, FromZprime=1<<13, Fromtprime=1<<14, Fromtprimebar=1<<15, Fromtaup=1<<16, Fromtaum=1<<17};
};

//////////////////////////////////////////////////////////////////////
// CLASS: GENPARTICLE ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class GenParticle : public GenBasicParticle {
  private:
    const  Analyse *myanalyse;
    UInt_t myindex;
    UInt_t motherfirst;
    UInt_t mothernum;
    UInt_t daughterfirst;
    UInt_t daughternum;
    bool   HasAnyMotherPDGId(vector<UInt_t> &visited, Int_t pdgid, bool antiparticle = true);
    UInt_t GetMotherIndex(UInt_t num) const;
    UInt_t GetDaughterIndex(UInt_t num) const;

  public:
    // constructors
    GenParticle(const Analyse *Myanalyse, UInt_t Myindex, Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t X, Float_t Y, Float_t Z, Int_t Status, Int_t Pdgid, UInt_t Motherfirst, UInt_t Mothernum, UInt_t Daughterfirst, UInt_t Daughternum);
    GenParticle() : myanalyse(0) {}
    // methods
    UInt_t      NumMothers() const { return(mothernum); }
    GenParticle GetMother(UInt_t num) const;
    UInt_t      NumDaughters() const { return(daughternum); }
    UInt_t      GetMyIndex() const { return(myindex); }
    GenParticle GetDaughter(UInt_t num) const;
    bool        HasAnyMotherPDGId(Int_t pdgid, bool antiparticle = true);
    bool operator ==(GenParticle &other) const { return(GetMyIndex() == other.GetMyIndex()); }
    bool operator !=(GenParticle &other) const { return(GetMyIndex() != other.GetMyIndex()); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: ECALHIT ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class EcalHit : public TVector3 {
  private:
    Double_t energy;

  public:
    // constructors
    EcalHit() {}
    EcalHit(Double_t Energy, Double_t X, Double_t Y, Double_t Z);
    // methods
    Double_t E() const { return(energy); }
};


//////////////////////////////////////////////////////////////////////
// CLASS: CLUSTER ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Cluster : public TLorentzVector {
    friend class Analyse;
  private:
    TVector3 position;
    Int_t    size;
    vector<EcalHit> hits;
    void     AddHit(EcalHit &hit);

  public:
    // constructors
    Cluster() {}
    Cluster(Double_t E, Double_t x, Double_t y, Double_t z, Int_t Size);
    // methods
    Int_t    Size() const { return(size); }
    TVector3 Position() const { return(position); }
    Int_t    NumHits() const { return(hits.size()); }
    EcalHit  Hits(UInt_t n) const { return(hits[n]); };
};

//////////////////////////////////////////////////////////////////////
// CLASS: TRACK //////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Track : public TLorentzVector {
  private:
    TVector3 outerpoint;
    TVector3 closestpoint;
    Double_t chi2;
    Double_t ndof;
    Double_t dxy;
    Double_t dxyerr;
    Double_t dz;
    Double_t dzerr;
    Int_t    charge;
    Int_t    nhits;
    Int_t    nmissinghits;
    Int_t    npixelhits;
    Int_t    npixellayers;
    Int_t    nstriplayers;
    Int_t    vtx;
    Float_t  dedxharmonic2;

  public:
    // constructors
    Track() {}
    Track(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Double_t Ox, Double_t Oy, Double_t Oz, Double_t Cx, Double_t Cy, Double_t Cz, Double_t Chi2, Double_t Ndof, Double_t Dxy, Double_t Dxyerr, Double_t Dz, Double_t Dzerr, Int_t Charge, Int_t Nhits, Int_t Nmissinghits, Int_t Npixelhits, Int_t Npixellayers, Int_t Nstriplayers, Int_t Vtx, Float_t Dedxharmonic2);
    // methods
    Double_t Chi2() const { return(chi2); }
    Double_t Ndof() const { return(ndof); }
    Double_t Chi2OverNdof() const { return(chi2/ndof); }
    Int_t    Charge() const { return(charge); }
    Int_t    NHits() const { return(nhits); }
    Int_t    NMissingHits() const { return(nmissinghits); }
    Int_t    NPixelHits() const { return(npixelhits); }
    Int_t    NPixelLayers() const { return(npixellayers); }
    Int_t    NStripLayers() const { return(nstriplayers); }
    Int_t    VertexNumber() const { return(vtx); }
    TVector3 ECalPoint() const { return(outerpoint); }
    TVector3 ClosestPoint() const { return(closestpoint); }
    TVector3 MomentumSpace() const { return(TVector3(Px(), Py(), Pz())); }
    Float_t  dEdxHarmonic2() const { return(dedxharmonic2); }
    Double_t Dxy() const { return(dxy); }
    Double_t DxyError() const { return(dxyerr); }
    Double_t DxySig() const { return(dxy/dxyerr); }
    Double_t Dz() const { return(dz); }
    Double_t DzError() const { return(dzerr); }
    Double_t DzSig() const { return(dz/dzerr); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: TRACKCOMPOSEDPARTICLE //////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class TrackComposedParticle : public TVector3 {
    friend class Analyse;
  private:
    TMatrixFSym   covmatrix;
    Double_t      chi2;
    Double_t      ndof;
    vector<Track> tracks;
    void          AddTrack(Track &newtrack);

  public:
    // constructors
    TrackComposedParticle() {}
    TrackComposedParticle(Double_t Vx, Double_t Vy, Double_t Vz, Double_t Chi2, Double_t Ndof, const Float_t *Cov);
    // methods
    Double_t    Chi2() const { return(chi2); }
    Double_t    Ndof() const { return(ndof); }
    Double_t    Chi2OverNdof() const { return(chi2/ndof); }
    Track       Tracks(Int_t n) const { return(tracks[n]); }
    Int_t       NumTracks() const { return(tracks.size()); }
    Double_t    XError() const { return(Sqrt(covmatrix(0,0))); }
    Double_t    YError() const { return(Sqrt(covmatrix(1,1))); }
    Double_t    ZError() const { return(Sqrt(covmatrix(2,2))); }
    TMatrixFSym CovMatrix() const { return(covmatrix); }
    Double_t    CovMatrix(Int_t i, Int_t j) const { return(covmatrix(i,j)); }
    Double_t    VertexSig3D(Vertex &vertex) const;
    Double_t    VertexSig2D(Vertex &vertex) const;
    TVector3    Momentum() const;
    TVector3    Distance(Vertex &vertex) const;
};

//////////////////////////////////////////////////////////////////////
// CLASS: TRIGGEROBJECT //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class TriggerObject {
  private:
    const Analyse *MA;
    UInt_t trigger;
    const vector<string> *triggernames;

  public:
    // constructors
    TriggerObject() {}
    TriggerObject(const Analyse *ma, const vector<string> *Triggernames, UInt_t Trigger);
    // methods
    Int_t Trigger(string name) const; //0: Trigger not available or matched. 1 not prescale fired, -1 not prescaled not fired, 2 prescaled fired, -2 prescalded not fired
    const vector<string> *TriggerNames() const { return(triggernames); }
    const Analyse *MyAn() const { return(MA); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: MUON ///////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Muon : public TLorentzVector, public TriggerObject {
    friend class Analyse;
  private:
    Float_t pterror;
    Float_t chi2;
    Float_t ndof;
    Track   innertrack;
    Track   outertrack;
    Float_t isolationr3track;
    Int_t   isolationr3ntrack;
    Float_t isolationr3ecal;
    Float_t isolationr3hcal;
    Float_t pfisolationr3_sumchargedhadronpt;
    Float_t pfisolationr3_sumchargedparticlept;
    Float_t pfisolationr3_sumneutralhadronet;
    Float_t pfisolationr3_sumphotonet;
    Float_t pfisolationr3_sumneutralhadronethighthreshold;
    Float_t pfisolationr3_sumphotonethighthreshold;
    Float_t pfisolationr3_sumpupt;
    Float_t pfisolationr4_sumchargedhadronpt;
    Float_t pfisolationr4_sumchargedparticlept;
    Float_t pfisolationr4_sumneutralhadronet;
    Float_t pfisolationr4_sumphotonet;
    Float_t pfisolationr4_sumneutralhadronethighthreshold;
    Float_t pfisolationr4_sumphotonethighthreshold;
    Float_t pfisolationr4_sumpupt;
    Float_t ecalenergy;
    Float_t hcalenergy;
    Int_t   charge;
    Int_t   muID;
    Int_t   numchambers;
    Int_t   numchamberswithsegments;
    Int_t   numvalidmuonhits;
    Int_t   nummatchedstations;
    Int_t   is_global;
    Int_t   is_tracker;
    Int_t   is_calo;
    Int_t   is_standalone;
    Int_t   is_pf_muon;
    Int_t   trackermuonquality;
    void    SetInnerTrack(Track &Innertrack) { innertrack = Innertrack; }
    void    SetOuterTrack(Track &Outertrack) { outertrack = Outertrack; }
  public:
    // constructors
    Muon() {}
    Muon(const Analyse *ma, UInt_t n, int correction = 0);
    // methods
    Float_t PtError() const { return(pterror); }
    Float_t IsoR3Track() const { return(isolationr3track); }
    Float_t IsoR3TrackRel() const { return(isolationr3track/Pt()); }
    Float_t IsoPFR3ChargedHadrons() const { return(pfisolationr3_sumchargedhadronpt); }
    Float_t IsoPFR3ChargedParticles() const { return(pfisolationr3_sumchargedparticlept); }
    Float_t IsoPFR3NeutralHadrons() const { return(pfisolationr3_sumneutralhadronet); }
    Float_t IsoPFR3SumPhotonEt() const { return(pfisolationr3_sumphotonet); }
    Float_t IsoPFR3SumPUPt() const { return(pfisolationr3_sumpupt); }
    Float_t IsoPFR4ChargedHadrons() const { return(pfisolationr4_sumchargedhadronpt); }
    Float_t IsoPFR4ChargedParticles() const { return(pfisolationr4_sumchargedparticlept); }
    Float_t IsoPFR4NeutralHadrons() const { return(pfisolationr4_sumneutralhadronet); }
    Float_t IsoPFR4Photons() const { return(pfisolationr4_sumphotonet); }
    Float_t IsoPFR4SumPUPt() const { return(pfisolationr4_sumpupt); }
    Int_t   IsoR3NTrack() const { return(isolationr3ntrack); }
    Float_t IsoR3ECal() const { return(isolationr3ecal); }
    Float_t IsoR3HCal() const { return(isolationr3hcal); }
    Float_t IsoR3Combined() const { return(IsoR3Track()+IsoR3ECal()+IsoR3HCal()); }
    Float_t IsoR3CombinedRel() const { return(IsoR3Combined()/Pt()); }
    Float_t ECalEnergy() const { return(ecalenergy); }
    Float_t HCalEnergy() const { return(hcalenergy); }
    Int_t   Charge() const { return(charge); }
    Int_t   NumChambers() const { return(numchambers); }
    Int_t   NumChambersWithSegments() const { return(numchamberswithsegments); }
    Int_t   NumValidMuonHits() const { return(numvalidmuonhits); }
    Int_t   NumMatchedStations() const { return(nummatchedstations); }
    bool    IsGlobal() const { return(is_global); }
    bool    IsTracker() const { return(is_tracker); }
    bool    IsStandAlone() const { return(is_standalone); }
    bool    IsCalo() const { return(is_calo); }
    bool    HasInnerTrack() const { return(hasinnertrack); }
    bool    HasOuterTrack() const { return(hasoutertrack); }
    bool    InOut() const { return((trackermuonquality & 1<<30) > 0); }
    bool    OutIn() const { return((trackermuonquality & 1<<31) > 0); }
    bool    IsPFMuon() const { return(is_pf_muon); }
    Float_t Dxy() const { return(innertrack.Dxy()); }
    Float_t DxyError() const { return(innertrack.DxyError()); }
    Float_t Dz() const { return(innertrack.Dz()); }
    Float_t DzError() const { return(innertrack.DzError()); }
    Float_t Chi2() const { return(chi2); }
    Float_t Ndof() const { return(ndof); }
    Float_t Chi2OverNdof() const { return(chi2/ndof); }
    Track   InnerTrack() const { return(innertrack); }
    Track   OuterTrack() const { return(outertrack); }
    Int_t   NumStations() const;
    bool    IsAll() const { return(trackermuonquality & 1 << 0); }
    bool    IsAllGlobalMuons() const { return(trackermuonquality & 1 << 1); }
    bool    IsAllStandAloneMuons() const { return(trackermuonquality & 1 << 2); }
    bool    IsAllTrackerMuons() const { return(trackermuonquality & 1 << 3); }
    bool    IsTrackerMuonArbitrated() const { return(trackermuonquality & 1 << 4); }
    bool    IsAllArbitrated() const { return(trackermuonquality & 1 << 5); }
    bool    IsGlobalMuonPromptTight() const { return(trackermuonquality & 1 << 6); }
    bool    IsTMLastStationLoose() const { return(trackermuonquality & 1 << 7); }
    bool    IsTMLastStationTight() const { return(trackermuonquality & 1 << 8); }
    bool    IsTM2DCompatibilityLoose() const { return(trackermuonquality & 1 << 9); }
    bool    IsTM2DCompatibilityTight() const { return(trackermuonquality & 1 << 10); }
    bool    IsTMOneStationLoose() const { return(trackermuonquality & 1 << 11); }
    bool    IsTMOneStationTight() const { return(trackermuonquality & 1 << 12); }
    bool    IsTMLastStationOptimizedLowPtLoose() const { return(trackermuonquality & 1 << 13); }
    bool    IsTMLastStationOptimizedLowPtTight() const { return(trackermuonquality & 1 << 14); }
    bool    IsGMTkChiCompatibility() const { return(trackermuonquality & 1 << 15); }
    bool    IsGMStaChiCompatibility() const { return(trackermuonquality & 1 << 16); }
    bool    IsGMTkKinkTight() const { return(trackermuonquality & 1 << 17); }
    bool    IsTMLastStationAngLoose() const { return(trackermuonquality & 1 << 18); }
    bool    IsTMLastStationAngTight() const { return(trackermuonquality & 1 << 19); }
    bool    IsTMOneStationAngLoose() const { return(trackermuonquality & 1 << 20); }
    bool    IsTMOneStationAngTight() const { return(trackermuonquality & 1 << 21); }
    bool    IsTMLastStationOptimizedBarrelLowPtLoose() const { return(trackermuonquality & 1 << 22); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: ELECTRON ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Electron : public TLorentzVector, public TriggerObject {
    friend class Analyse;
  private:
    TVector3 outerpoint;
    TVector3 closestpoint;
    Double_t dxy;
    Double_t dxyerr;
    Double_t dz;
    Double_t dzerr;
    Float_t  correctedecalenergy;
    Float_t  esuperclusterovertrack;
    Float_t  eseedclusterovertrack;
    Float_t  deltaetasuperclustertrack;
    Float_t  deltaphisuperclustertrack;
    Float_t  e1x5;
    Float_t  e2x5;
    Float_t  e5x5;
    Float_t  r9;
    Float_t  sigmaetaeta;
    Float_t  sigmaietaieta;
    Float_t  sigmaiphiiphi;
    Float_t  ehcaloverecaldepth1;
    Float_t  ehcaloverecaldepth2;
    Float_t  ehcaltoweroverecaldepth1;
    Float_t  ehcaltoweroverecaldepth2;
    Float_t  isolationr3track;
    Float_t  isolationr3ecal;
    Float_t  isolationr3hcal;
    Float_t  isolationr4track;
    Float_t  isolationr4ecal;
    Float_t  isolationr4hcal;
    Float_t  isolationpfr3charged;
    Float_t  isolationpfr3photon;
    Float_t  isolationpfr3neutral;
    Float_t  trackchi2;
    Float_t  trackndof;
    Int_t    nhits;
    Int_t    nmissinghits;
    Int_t    npixelhits;
    Int_t    npixellayers;
    Int_t    nstriplayers;
    Int_t    nhitsexpected;
    Float_t  convdist;
    Float_t  convdcot;
    Float_t  convradius;
    UInt_t   gapinfo;
    Float_t  fbrems;
    Int_t    numbrems;
    Int_t    charge;
    Byte_t   info;
    Byte_t   eID;
    Int_t    vtx;
    //vector<SuperCluster> supercluster;
    bool     WorkingPoint(Int_t Missinghits, Float_t Convdist, Float_t Convdcot, Float_t Sigmaietaieta, Float_t Deltaphisctrack, Float_t Deltaetasctrack, Float_t Ehcaloverecal) const;

  public:
    // constructors
    Electron() {}
    Electron(const Analyse *ma, UInt_t n, int correction = 0);
    // methods
    TVector3 EcalPoint() const { return(outerpoint); }
    TVector3 ClosestPoint() const { return(closestpoint); }
    Float_t  CorrectedEcalEnergy() const { return(correctedecalenergy); }
    Float_t  ESuperClusterOverTrack() const { return(esuperclusterovertrack); }
    Float_t  ESeedClusterOverTrack() const { return(eseedclusterovertrack); }
    Float_t  DeltaEtaSuperClusterTrack() const { return(deltaetasuperclustertrack); }
    Float_t  DeltaPhiSuperClusterTrack() const { return(deltaphisuperclustertrack); }
    Float_t  E1x5() const { return(e1x5); }
    Float_t  E2x5() const { return(e2x5); }
    Float_t  E5x5() const { return(e5x5); }
    Float_t  R9() const { return(r9); }
    Float_t  SigmaEtaEta() const { return(sigmaetaeta); }
    Float_t  SigmaIEtaIEta() const { return(sigmaietaieta); }
    Float_t  SigmaIPhiIPhi() const { return(sigmaiphiiphi); }
    Float_t  EHcalOverECal() const { return(ehcaloverecaldepth1+ehcaloverecaldepth2); }
    Float_t  EHcalOverECalDepth1() const { return(ehcaloverecaldepth1); }
    Float_t  EHcalOverECalDepth2() const { return(ehcaloverecaldepth2); }
    Float_t  EHcalTowerOverECal() const { return(ehcaltoweroverecaldepth1+ehcaltoweroverecaldepth2); }
    Float_t  EHcalTowerOverECalDepth1() const { return(ehcaltoweroverecaldepth1); }
    Float_t  EHcalTowerOverECalDepth2() const { return(ehcaltoweroverecaldepth2); }
    Float_t  IsoR3Track() const { return(isolationr3track); }
    Float_t  IsoR3TrackRel() const { return(isolationr3track/Pt()); }
    Float_t  IsoR3ECal() const { return(isolationr3ecal); }
    Float_t  IsoR3ECalRel() const { return(IsoR3ECal()/Pt()); }
    Float_t  IsoR3HCal() const { return(isolationr3hcal); }
    Float_t  IsoR3HCalRel() const { return(IsoR3HCal()/Pt()); }
    Float_t  IsoR3Combined() const { return(IsoR3Track() + IsoR3ECal() + IsoR3HCal()); }
    Float_t  IsoR3CombinedRel() const { return(IsoR3Combined()/Pt()); }
    Float_t  IsoR4Track() const { return(isolationr4track); }
    Float_t  IsoR4ECal() const { return(isolationr4ecal); }
    Float_t  IsoR4HCal() const { return(isolationr4hcal); }
    Float_t  IsoPFR3Charged() const { return(isolationpfr3charged); }
    Float_t  IsoPFR3Photon() const { return(isolationpfr3photon); }
    Float_t  IsoPFR3Neutral() const { return(isolationpfr3neutral); }
    Float_t  TrackChi2() const { return(trackchi2); }
    Float_t  TrackNdof() const { return(trackndof); }
    Float_t  TrackChi2OverNdof() const { return(trackchi2/trackndof); }
    Int_t    NHits() const { return(nhits); }
    Int_t    NMissingHits() const { return(nmissinghits); }
    Int_t    NPixelHits() const { return(npixelhits); }
    Int_t    NPixelLayers() const { return(npixellayers); }
    Int_t    NStripLayers() const { return(nstriplayers); }
    Int_t    NHitsExpected() const { return(nhitsexpected); }
    Float_t  ConversionDist() const { return(convdist); }
    Float_t  ConversionDCot() const { return(convdcot); }
    Float_t  ConversionRadius() const { return(convradius); }
    Float_t  FractionBrems() const { return(fbrems); }
    Int_t    NumBrems() const { return(numbrems); }
    Int_t    Charge() const { return(charge); }
    Int_t    VertexNumber() const { return(vtx); }
    bool     HasMatchedConversion() const { return(info & 1<<1); }
    Double_t Dxy() const { return(dxy); }
    Double_t DxyError() const { return(dxyerr); }
    Double_t Dz() const { return(dz); }
    Double_t DzError() const { return(dzerr); }
    bool     WP95_v1(Int_t combined = 0) const;
    bool     WP90_v1(Int_t combined = 0) const;
    bool     WP85_v1(Int_t combined = 0) const;
    bool     WP80_v1(Int_t combined = 0) const;
    bool     WP70_v1(Int_t combined = 0) const;
    bool     WP60_v1(Int_t combined = 0) const;
    bool     IsEB() const { return(gapinfo & 1<<0); }
    bool     IsEE() const { return(gapinfo & 1<<1); }
    bool     IsEBGap() const { return(gapinfo & 1<<2); }
    bool     IsEBEtaGap() const { return(gapinfo & 1<<3); }
    bool     IsEBPhiGap() const { return(gapinfo & 1<<4); }
    bool     IsEEGap() const { return(gapinfo & 1<<5); }
    bool     IsEERingGap() const { return(gapinfo & 1<<6); }
    bool     IsEEDeeGap() const { return(gapinfo & 1<<7); }
    bool     IsEBEEGap() const { return(gapinfo & 1<<8); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: VERTEX /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Vertex : public TVector3 {
  private:
    Double_t    chi2;
    Double_t    ndof;
    UInt_t      isfake;
    UInt_t      isvalid;
    Int_t       ntracks;
    TMatrixFSym covmatrix;

  public:
    // constructors
    Vertex(Double_t X, Double_t Y, Double_t Z, UInt_t Isvalid, UInt_t Isfake, Double_t Chi2, Double_t Ndof, Int_t Ntracks, Float_t Cov0, Float_t Cov1, Float_t Cov2, Float_t Cov3, Float_t Cov4, Float_t Cov5);
    Vertex() : covmatrix(3) {}
    // methods
    Double_t    Chi2() const { return(chi2); }
    Double_t    Ndof() const { return(ndof); }
    Double_t    Chi2OverNdof() const { return(chi2/ndof); }
    Int_t       NTracks() const { return(ntracks); }
    Double_t    XError() const { return(Sqrt(covmatrix(0,0))); }
    Double_t    YError() const { return(Sqrt(covmatrix(1,1))); }
    Double_t    ZError() const { return(Sqrt(covmatrix(2,2))); }
    bool        IsFake() const { return(isfake); }
    bool        IsValid() const { return(isvalid); }
    TMatrixFSym CovMatrix() const { return(covmatrix); }
    Double_t    CovMatrix(Int_t i, Int_t j) const { return(covmatrix(i,j)); }
    bool        IsGood() const { return(!isfake && ndof > 4 && Z() < 24. && Perp() < 2.); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: PHOTON /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Photon : public TLorentzVector, public TriggerObject {
    friend class Analyse;
  private:
    Float_t e1x5;
    Float_t e2x5;
    Float_t e3x3;
    Float_t e5x5;
    Float_t sigmaietaieta;
    Float_t sigmaietaiphi;
    Float_t sigmaiphiiphi;
    Float_t ehcaloverecaldepth1;
    Float_t ehcaloverecaldepth2;
    Float_t ehcaltoweroverecaldepth1;
    Float_t ehcaltoweroverecaldepth2;
    Float_t maxenergyxtal;
    Float_t isolationr3track;
    Float_t isolationr3trackhollow;
    UInt_t  isolationr3ntrack;
    UInt_t  isolationr3ntrackhollow;
    Float_t isolationr3ecal;
    Float_t isolationr3hcal;
    Float_t isolationr4track;
    Float_t isolationr4trackhollow;
    UInt_t  isolationr4ntrack;
    UInt_t  isolationr4ntrackhollow;
    Float_t isolationr4ecal;
    Float_t isolationr4hcal;
    Float_t isolationpfr3charged;
    Float_t isolationpfr3photon;
    Float_t isolationpfr3neutral;
    Float_t isolationpfr4charged;
    Float_t isolationpfr4photon;
    Float_t isolationpfr4neutral;
    Float_t isolationpfr4noscfootprintcharged;
    Float_t isolationpfr4noscfootprintphoton;
    Float_t isolationpfr4noscfootprintneutral;
    UChar_t info;
    UInt_t  gapinfo;
    //vector<SuperCluster> supercluster;

  public:
    // constructors
    Photon() {}
    Photon(const Analyse *ma, UInt_t n);
    // methods
    Float_t E1x5() const { return(e1x5); }
    Float_t E2x5() const { return(e2x5); }
    Float_t E3x3() const { return(e3x3); }
    Float_t E5x5() const { return(e5x5); }
    Float_t R9() const { return(0); }
    Float_t R25() const { return(0); }
    Float_t SigmaIEtaIEta() const { return(sigmaietaieta); }
    Float_t SigmaIEtaIPhi() const { return(sigmaietaiphi); }
    Float_t SigmaIPhiIPhi() const { return(sigmaiphiiphi); }
    Float_t MaxEnergyXtal() const { return(maxenergyxtal); }
    Float_t EHCalOverECal() const { return(ehcaloverecaldepth1+ehcaloverecaldepth2); }
    Float_t EHCalOverECalDepth1() const { return(ehcaloverecaldepth1); }
    Float_t EHCalOverECalDepth2() const { return(ehcaloverecaldepth2); }
    Float_t EHCalTowerOverECal() const { return(ehcaltoweroverecaldepth1+ehcaltoweroverecaldepth2); }
    Float_t EHCalTowerOverECalDepth1() const { return(ehcaltoweroverecaldepth1); }
    Float_t EHCalTowerOverECalDepth2() const { return(ehcaltoweroverecaldepth2); }
    Float_t IsoR3Track() const { return(isolationr3track); }
    Float_t IsoR3TrackHollow() const { return(isolationr3trackhollow); }
    UInt_t  IsoR3NTrack() const { return(isolationr3ntrack); }
    UInt_t  IsoR3NTrackHollow() const { return(isolationr3ntrackhollow); }
    Float_t IsoR3ECal() const { return(isolationr3ecal); }
    Float_t IsoR3HCal() const { return(isolationr3hcal); }
    Float_t IsoR4Track() const { return(isolationr4track); }
    Float_t IsoR4TrackHollow() const { return(isolationr4trackhollow); }
    UInt_t  IsoR4NTrackHollow() const { return(isolationr4ntrackhollow); }
    Float_t IsoR4ECal() const { return(isolationr4ecal); }
    Float_t IsoR4HCal() const { return(isolationr4hcal); }
    Float_t IsoPFR3Charged() const { return(isolationpfr3charged); }
    Float_t IsoPFR3Photon() const { return(isolationpfr3photon); }
    Float_t IsoPFR3Neutral() const { return(isolationpfr3neutral); }
    Float_t IsoPFR4Charged() const { return(isolationpfr4charged); }
    Float_t IsoPFR4Photon() const { return(isolationpfr4photon); }
    Float_t IsoPFR4Neutral() const { return(isolationpfr4neutral); }
    Float_t IsoPFR3NoFPCharged() const { return(isolationpfr4noscfootprintcharged); }
    Float_t IsoPFR3NoFPPhoton() const { return(isolationpfr4noscfootprintphoton); }
    Float_t IsoPFR3NoFPNeutral() const { return(isolationpfr4noscfootprintneutral); }
    bool    isPhoton() const { return((info & 1) > 0); }
    bool    HasConversionTracks() const { return((info & 1<<1) > 0); }
    bool    HasPixelSeed() const { return((info & 1<<2) > 0); }
    bool    HasPromptElectron() const { return((info & 1<<3) > 0); }
    bool    IsEB() const { return((gapinfo & 1<<0) > 0); }
    bool    IsEE() const { return((gapinfo & 1<<1) > 0); }
    bool    IsEBGap() const { return((gapinfo & 1<<2) > 0); }
    bool    IsEBEtaGap() const { return((gapinfo & 1<<3) > 0); }
    bool    IsEBPhiGap() const { return((gapinfo & 1<<4) > 0); }
    bool    IsEEGap() const { return((gapinfo & 1<<5) > 0); }
    bool    IsEERingGap() const { return((gapinfo & 1<<6) > 0); }
    bool    IsEEDeeGap() const { return((gapinfo & 1<<7) > 0); }
    bool    IsEBEEGap() const { return((gapinfo & 1<<8) > 0); }
};


//////////////////////////////////////////////////////////////////////
// CLASS: JET ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//Not all Variables are filled for all Jet types.
class Jet : public TLorentzVector {
  private:
    Float_t hadronicenergy;
    Float_t chargedhadronicenergy;
    Float_t emenergy;
    Float_t chargedemenergy;
    Float_t hfemenergy;
    Float_t hfhadronicenergy;
    Float_t electronenergy;
    Float_t muonenergy;
    Int_t   chargedmulti;
    Int_t   neutralmulti;
    Int_t   hfemmulti;
    Int_t   hfhadronicmulti;
    Int_t   electronmulti;
    Int_t   muonmulti;
    Float_t chargeda;
    Float_t chargedb;
    Float_t neutrala;
    Float_t neutralb;
    Float_t alla;
    Float_t allb;
    Float_t chargedfractionmv;
    Float_t energycorr;
    Float_t energycorrunc;
    Float_t energycorrl7uds;
    Float_t energycorrl7bottom;
    Int_t   mcflavour;
    Float_t puidfull;
    Float_t puidsimple;
    Float_t puidcutbased;
    bool    mymvonly;
    vector<Float_t> btag;

  public:
    // contsructors
    Jet() {}
    //Jet(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t Hadronicenergy, Float_t Chargedhadronicenergy, Float_t Emenergy, Float_t Chargedemenergy, Float_t Hfemenergy, Float_t Hfhadronicenergy, Float_t Electronenergy, Float_t Muonenergy, Int_t Chargedmulti, Int_t Neutralmulti, Int_t Hfemmulti, Int_t Hfhadronicmulti, Int_t Electronmulti, Int_t Muonmulti, Float_t Chargeda, Float_t Chargedb, Float_t Neutrala, Float_t Neutralb, Float_t Alla, Float_t Allb, Float_t Chargedfractionmv, Float_t Energycorr,Float_t Energycorrunc, Float_t Energycorrl7uds, Float_t Energycorrl7bottom, const Float_t *Btag, Int_t Mcflavour, Float_t Puidfull, Float_t Puidsimple, Float_t Puidcutbased);
    Jet(Double_t E, Double_t Px, Double_t Py, Double_t Pz, Float_t Hadronicenergy, Float_t Chargedhadronicenergy, Float_t Emenergy, Float_t Chargedemenergy, Float_t Hfemenergy, Float_t Hfhadronicenergy, Float_t Electronenergy, Float_t Muonenergy, Int_t Chargedmulti, Int_t Neutralmulti, Int_t Hfemmulti, Int_t Hfhadronicmulti, Int_t Electronmulti, Int_t Muonmulti, Float_t Chargeda, Float_t Chargedb, Float_t Neutrala, Float_t Neutralb, Float_t Alla, Float_t Allb, Float_t Chargedfractionmv, Float_t Energycorr,Float_t Energycorrunc, Float_t Energycorrl7uds, Float_t Energycorrl7bottom, const vector<Float_t> Btag, Int_t Mcflavour, Float_t Puidfull, Float_t Puidsimple, Float_t Puidcutbased);
    // methods
    Float_t HadEnergyFraction() const { return(hadronicenergy/E()/energycorr); }
    Float_t HadEnergy() const { return(hadronicenergy); }
    Float_t EMEnergyFraction() const { return(emenergy/E()/energycorr); }
    Float_t EMEnergy() const { return(emenergy); }
    Float_t ChargedHadEnergyFraction() const { return(chargedhadronicenergy/E()/energycorr); }
    Float_t ChargedHadEnergy() const { return(chargedhadronicenergy); }
    Float_t ChargedEMEnergyFraction() const { return(chargedemenergy/E()/energycorr); }
    Float_t ChargedEMEnergy() const { return(chargedemenergy); }
    Float_t HFEMEnergy() const { return(hfemenergy); }
    Float_t HFHadEnergy() const { return(hfhadronicenergy); }
    Float_t ElectronEnergy() const { return(electronenergy); }
    Float_t MuonEnergy() const { return(muonenergy); }
    Int_t   ChargedMulti() const { return(chargedmulti); }
    Int_t   NeutralMulti() const { return(neutralmulti); }
    Int_t   HFEMMulti() const { return(hfemmulti); }
    Int_t   HFHadMulti() const { return(hfhadronicmulti); }
    Int_t   ElectronMulti() const { return(electronmulti); }
    Int_t   MuonMulti() const { return(muonmulti); }
    Float_t ChargedPrincipalAxisA() const { return(chargeda); }
    Float_t ChargedPrincipalAxisB() const { return(chargedb); }
    Float_t NeutralPrincipalAxisA() const { return(neutrala); }
    Float_t NeutralPrincipalAxisB() const { return(neutralb); }
    Float_t PrincipalAxisA() const { return(alla); }
    Float_t PrincipalAxisB() const { return(allb); }
    Float_t ChargedMomentumFractionFromMV() const { return(chargedfractionmv); }
    void    ScaleMV(bool mvonly);
    Float_t BTag(Int_t n) const;
    Float_t JECRaw() const { return(energycorr); }
    Float_t JECUncertainty() const { return(energycorrunc); }
    Float_t JECL7UDSCorretion() const { return(energycorrl7uds); }
    Float_t JECL7BottomCorretion() const { return(energycorrl7bottom); }
    Float_t trackCountingHighPurBJetTags() const { return(btag.at(0)); }
    Float_t trackCountingHighEffBJetTags() const { return(btag.at(1)); }
    Float_t simpleSecondaryVertexHighPurBJetTags() const { return(btag.at(2)); }
    Float_t simpleSecondaryVertexHighEffBJetTags() const { return(btag.at(3)); }
    Float_t combinedSecondaryVertexBJetTags() const { return(btag.at(4)); }
    Float_t combinedSecondaryVertexMVABJetTags() const { return(btag.at(5)); }
    Float_t puIdFull() const { return(puidfull); }
    Float_t puIdSimple() const { return(puidsimple); }
    Float_t puIdCutBased() const { return(puidcutbased); }
    Int_t   MCFlavour() const { return(mcflavour); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: TAU ////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Tau : public TLorentzVector, public TriggerObject {
    friend class Analyse;
  private:
    Float_t isolationneutralspt;
    UInt_t  isolationneutralsnum;
    Float_t isolationchargedpt;
    UInt_t  isolationchargednum;
    Float_t isolationgammapt;
    UInt_t  isolationgammanum;
    Char_t  charge;
    UInt_t  dishps;
    Float_t emfraction;
    Float_t hcaltotoverplead;
    Float_t hcal3x3overplead;
    Float_t ecalstripsumeoverplead;
    Float_t bremsrecoveryeoverplead;
    Float_t calocomp;
    Float_t segcomp;
    Jet     jet;
    vector<Track> tracks;
    const vector<string> *taudiscriminators;
  public:
    // constructor
    Tau() {}
    Tau(const Analyse *ma, UInt_t n);
    // methods
    Float_t IsoNeutralsPt() const { return(isolationneutralspt); }
    UInt_t  IsoNeutralsNum() const { return(isolationneutralsnum); }
    Float_t IsoChargedPt() const { return(isolationchargedpt); }
    UInt_t  IsoChargedNum() const { return(isolationchargednum); }
    Float_t IsoGammaPt() const { return(isolationgammapt); }
    UInt_t  IsoGammaNum() const { return(isolationgammanum); }
    Int_t   Charge() const { return(charge); }
    Float_t EMFraction() const { return(emfraction); }
    Float_t HCalTotOverPLead() const { return(hcaltotoverplead); }
    Float_t HCal3x3OverPLead() const { return(hcal3x3overplead); }
    Float_t ECalStripSumEOverPLead() const { return(ecalstripsumeoverplead); }
    Float_t BremsRecoveryEOverPLead() const { return(bremsrecoveryeoverplead); }
    Float_t CaloComp() const { return(calocomp); }
    Float_t SegComp() const { return(segcomp); }
    Track   LeadingTrack() const;
    UInt_t  NumTracks() const { return(tracks.size()); }
    Track   GetTrack(UInt_t num) const { return(tracks[num]); }
    Int_t   TauDiscriminator(string disname) const; 
};

//////////////////////////////////////////////////////////////////////
// CLASS: BEAMSPOT ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class BeamSpot : public TVector3 {
  private:
    TMatrixFSym covmatrix;
    Double_t    xwidth;
    Double_t    ywidth;
    Double_t    zsigma;

  public:
    // constructors
    BeamSpot() {}
    BeamSpot(Double_t X, Double_t Y, Double_t Z, const Float_t *Cov, Double_t Xwidth, Double_t Ywidth, Double_t Zsigma);
    // methods
    Double_t    XError() const { return(Sqrt(covmatrix(0,0))); }
    Double_t    YError() const { return(Sqrt(covmatrix(1,1))); }
    Double_t    ZError() const { return(Sqrt(covmatrix(2,2))); }
    TMatrixFSym CovMatrix() const { return(covmatrix); }
    Double_t    CovMatrix(Int_t i, Int_t j) const { return(covmatrix(i,j)); }
    Double_t    XWidth() const { return(xwidth); }
    Double_t    YWidth() const { return(ywidth); }
    Double_t    ZSigma() const { return(zsigma); }
};

//////////////////////////////////////////////////////////////////////
// CLASS: LUMINOSITY /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class Luminosity {
  private:
    Float_t lumival;
    Float_t lumierr;
    Float_t livefraction;
    Float_t deadfraction;
    Float_t avgpu;
    UInt_t  quality;
    UInt_t  numevents;
    UInt_t  numeventsorig;
    UInt_t  hlttable;
    UInt_t  l1techtable;
    UInt_t  l1algotable;
    UInt_t  counter;
    vector<string> filenames;

  public:
    // constructors
    Luminosity(Float_t Lumival, Float_t Lumierr, Float_t Livefraction, Float_t Deadfraction, Float_t Avgpu, UInt_t Quality, UInt_t Numevents, UInt_t Numeventsprocessed, UInt_t Hlttable, UInt_t L1techtable, UInt_t L1algotable, vector<string> newfiles) : lumival(Lumival), lumierr(Lumierr), livefraction(Livefraction), deadfraction(Deadfraction), avgpu(Avgpu), quality(Quality), numevents(Numevents), numeventsorig(Numeventsprocessed), hlttable(Hlttable), l1techtable(L1techtable), l1algotable(L1algotable), counter(0) {
        AddFiles(newfiles);
    }
    Luminosity() : lumival(-1), numevents(0), numeventsorig(0), counter(0) {}
    // methods
    bool operator ==(const bool test) {
        bool val = true;
        if(lumival == -1) val = false;
        return(test == val);
    }
    void operator ++() {
        counter++;
    }
    Luminosity &operator +=(Luminosity &other) {
        numevents += other.NumEvents();
        numeventsorig += other.NumEventsOrig();
        AddFiles(other.GetFiles());
        return(*this);
    }
    void AddFiles(const vector<string> &newfiles);
    const vector<string> &GetFiles() const {
        return(filenames);
    }
    string GetFilesString() const;
    Float_t LumiValue() const {
        return(lumival);
    }
    Float_t AvgPU() const {
        return(avgpu);
    }
    Float_t LumiError() const {
        return(lumierr);
    }
    Float_t LiveFraction() const {
        return(livefraction);
    }
    Float_t DeadFraction() const {
        return(deadfraction);
    }
    UInt_t Quality() const {
        return(quality);
    }
    UInt_t NumEvents() const {
        return(numevents);
    }
    UInt_t NumEventsOrig() const {
        return(numeventsorig);
    }
    UInt_t NumEventsProcessed() const {
        return(counter);
    }
    Float_t ProcessedLumi() const {
        if(numevents != 0) {
            return(lumival*counter/numevents);
        } else {
            return(lumival);
        }
    }
    Float_t ProcessedFraction() const {
        if(numevents != 0) {
            return(Float_t(counter)/numevents);
        } else {
            return(1.);
        }
    }
    UInt_t HLTTable() const {
        return(hlttable);
    }
    UInt_t L1TechTable() const {
        return(l1techtable);
    }
    UInt_t L1AlgoTable() const {
        return(l1algotable);
    }
    void Value(Float_t myvalue) {
        lumival = myvalue;
    }
    void ValueErr(Float_t myvalueerr) {
        lumierr = myvalueerr;
    }
    void LiveFrac(Float_t mylivefrac) {
        livefraction = mylivefrac;
    }
    void DeadFrac(Float_t mydeadfrac) {
        deadfraction = mydeadfrac;
    }
    void AvgPU(Float_t myavgpu) {
        avgpu = myavgpu;
    }
    void Quality(UInt_t myquality) {
        quality = myquality;
    }
    void EventsFiltered(UInt_t myeventsfiltered) {
        numevents = myeventsfiltered;
    }
    void EventsOriginal(UInt_t myeventsorig) {
        numeventsorig = myeventsorig;
    }
    void EventsProcessed(UInt_t myeventsprocessed) {
        counter = myeventsprocessed;
    }
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class RunInfo {
private:
    UInt_t runnumber;
    vector<string> hltnames;
    Int_t hltnum;
    Int_t hlttablesnum;
    vector<vector<UInt_t> > hltprescales;
    Int_t l1algonum;
    Int_t l1algotablesnum;
    vector<vector<UInt_t> > l1algoprescales;
    Int_t l1technum;
    Int_t l1techtablesnum;
    vector<vector<UInt_t> > l1techprescales;
    vector<string> hltmunames;
    vector<string> hltelnames;
    vector<string> hltphotonnames;
    vector<string> hlttaunames;
    vector<string> hltjetnames;
    vector<string> taudiscriminators;
public:
    RunInfo(UInt_t Number, UInt_t Hltcount, string Hltnames, string Hltmunames, string Hltelnames, string Hltphotonnames, string Hlttaunames, string Hltjetnames, string Taudiscriminators, UInt_t Hltprescaletablescount, UInt_t *Hltprescaletables, UInt_t L1algocount, UInt_t L1algoprescaletablescount, UInt_t *L1algoprescaletables, UInt_t L1techcount, UInt_t L1techprescaletablescount, UInt_t *L1techprescaletables);
    RunInfo() : runnumber(0) {}
    UInt_t Run() const {
        return(runnumber);
    }
    UInt_t RunNumber() const {
        return(runnumber);
    }
    UInt_t NumHLTTables() const {
        return(hlttablesnum);
    }
    UInt_t NumHLT() const {
        return(hltnum);
    }
    string HLTAllNames() const;
    string HLTMuonAllNames() const;
    string HLTElectronAllNames() const;
    string HLTTauAllNames() const;
    string HLTPhotonAllNames() const;
    string HLTJetAllNames() const;
    string TauDiscriminatorsAllNames() const;
    UInt_t HLTPrescale(UInt_t trigger, UInt_t table) const {
        if(trigger < NumHLT() && table < NumHLTTables()) {
            return(hltprescales[trigger][table]);
        } else {
            return(1);
        }
    }
    string HLTName(UInt_t trigger) const {
        return(hltnames[trigger]);
    }
    vector<string> MatchTriggerNames(string name);
    Int_t HLTIndex(string name) const;
    UInt_t NumL1TechTables() const {
        return(l1techtablesnum);
    }
    UInt_t NumL1Tech() const {
        return(l1technum);
    }
    UInt_t L1TechPrescale(UInt_t trigger, UInt_t table) const {
        if(trigger < NumL1Tech() && table < NumL1TechTables()) {
            return(l1techprescales[trigger][table]);
        } else {
            return(1);
        }
    }
    UInt_t NumL1AlgoTables() const {
        return(l1algotablesnum);
    }
    UInt_t NumL1Algo() const {
        return(l1algonum);
    }
    UInt_t L1AlgoPrescale(UInt_t trigger, UInt_t table) const {
        if(trigger < NumL1Algo() && table < NumL1AlgoTables()) {
            return(l1algoprescales[trigger][table]);
        } else {
            return(1);
        }
    }
    const vector<string> *GetHLTMuonNames() const {
        return(&hltmunames);
    }
    const vector<string> *GetHLTElectronNames() const {
        return(&hltelnames);
    }
    const vector<string> *GetHLTTauNames() const {
        return(&hlttaunames);
    }
    const vector<string> *GetHLTPhotonNames() const {
        return(&hltphotonnames);
    }
    const vector<string> *GetHLTJetNames() const {
        return(&hltjetnames);
    }
    const vector<string> *GetTauDiscriminators() const {
        return(&taudiscriminators);
    }
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class TriggerLumi {
private:
    Int_t index;
    UInt_t prescale;
public:
    TriggerLumi() : index(-1), prescale(0) {}
    TriggerLumi(Int_t Index, UInt_t Prescale) : index(Index), prescale(Prescale) {}
    Int_t Index() const {
        return(index);
    }
    UInt_t Prescale() const {
        return(prescale);
    }
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class TriggerRun {
private:
    Float_t lumi;
    map< UInt_t, TriggerLumi> lumiinfo;
public:
    TriggerRun() : lumi(-1.) {}
    void Lumi(Float_t runlumi) {
        lumi = runlumi;
    }
    void SetBlock(UInt_t block, Int_t index, UInt_t prescale);
    TriggerLumi GetBlock(UInt_t block);
    Float_t Lumi() const {
        return(lumi);
    }
    map< UInt_t, TriggerLumi>::iterator Begin() {
        return(lumiinfo.begin());
    }
    map< UInt_t, TriggerLumi>::iterator End() {
        return(lumiinfo.end());
    }
};

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
class TriggerSelection {
private:
    Analyse *AN;
    map<UInt_t, TriggerRun> runinfo;
public:
    TriggerSelection(Analyse *an, vector<string> names, bool useprescaled = false);
    Int_t Result();
    Float_t LumiUsed(Int_t format = 0);
    Float_t LumiBeforeEvent();
    string GetTriggerName(UInt_t run = 0, UInt_t lumiblock = 0);
    void PrintInfo();
};

#endif
