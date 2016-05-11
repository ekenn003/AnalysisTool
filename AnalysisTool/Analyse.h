#ifndef ACANALYSE
#define ACANALYSE
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixFSym.h>
#include <TStopwatch.h>
#include <TRandom3.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>

#include "DataForm.h"
//#include "rochester/rochcor2012v2.h"

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include <sys/types.h>
#include <dirent.h>

using namespace std;
using namespace TMath;

const UInt_t M_trackmaxcount = 500;
const UInt_t M_superclustermaxcount = 1000;
const UInt_t M_superclustermembermaxcount = 1000;
const UInt_t M_superclusterhitmaxcount = 5000;
const UInt_t M_primvertexmaxcount = 500;
const UInt_t M_muonmaxcount = 500;
const UInt_t M_taumaxcount = 500;
const UInt_t M_electronmaxcount = 500;
const UInt_t M_photonmaxcount = 500;
const UInt_t M_conversionmaxcount = 500;
const UInt_t M_jetmaxcount = 500;
const UInt_t M_musecverticesmaxcount = 500;
const UInt_t M_secverticesmaxcount = 1000;
const UInt_t M_genallparticlesmaxcount = 10000;
const UInt_t M_genparticlesmaxcount = 500;
const UInt_t M_genjetmaxcount = 500;
const UInt_t M_genmotherdaughtermaxcount = 100000;

bool sortac1b(string a, string b);
bool sortroot(string a, string b);

class Analyse
{
    friend class GenParticle;
    friend class Tau;
    friend class Muon;
    friend class Electron;
    friend class Jet;
    friend class Photon;
    friend class TriggerSelection;

  private:
    vector<string>   filenames;
    vector<Long64_t> filelimits;
    Long64_t currentmin, currentmax;
    TFile    *currentfile;
    Int_t    currentloadtype;
    TTree    *tree;
    string   treename;
    Int_t    printinfo;
    Long64_t processed;
    Int_t    loadbeamspot;
    Int_t    loadmuons;
    Int_t    loadelectrons;
    Int_t    loadphotons;
    Int_t    loadtaus;
    Int_t    loadmet;
    Int_t    loadak4pfchsjets;
    Int_t    loadtracks;
    Int_t    loadprimvertices;
    Int_t    loadtrigger;
    Int_t    loadgeninfo;
    Int_t    loadgenparticles;
    Int_t    loadgenjets;
    Int_t    loadallgenparticles;

    void SetLoad();
    void Load();
    void GetEvent(Long64_t num, Int_t loadtype = 0);

    //dupcheck
    bool duplicatecheck;
    map< UInt_t, map< UInt_t, map< Double_t, UInt_t> > > eventlist;

    //Skimming
    TTree  *skimtree;
    map< UInt_t, map< UInt_t, UInt_t > > selected;
    string skimfilename;
    TFile  *skimfile;

    //Lumi calculation
    map< UInt_t, map< UInt_t, Luminosity > > lumilist;
    map< UInt_t, RunInfo> runlist;
    map< string, TriggerSelection *> triggerselections;
    bool   lumicalculation;
    bool   IsInRange(UInt_t Run, UInt_t LumiBlock);
    UInt_t minRun;
    UInt_t minLumi;
    UInt_t maxRun;
    UInt_t maxLumi;
    //JSON filter
    bool jsonfilter;
    map< UInt_t, map< UInt_t, bool > > jsonlist;

    //MPI
    Int_t batch_myid;
    Int_t batch_numjobs;
    bool  batch_emptyjob;
    map<  UInt_t, map< UInt_t, bool > > batchselection;
    bool  IsBatchSelected(UInt_t run, UInt_t lumiblock);

    //Pileup
    bool   usepileupinfo;
    double pileupscale;
    vector<Double_t> pileUpDistMinus;
    vector<Double_t> pileUpDist;
    vector<Double_t> pileUpDistPlus;
    bool useprimvertexinfo;
    vector<Double_t> primVertexDist;

    TRandom3 *analysisrandom;

    //Data
    Bool_t   isdata;
    UInt_t   errors;
    Double_t event_nr;
    UInt_t   event_luminosityblock;
    UInt_t   event_run;
    UInt_t   event_timeunix;
    UInt_t   event_timemicrosec;
    UChar_t  trigger_level1bits[8];
    UChar_t  trigger_level1[128];
    UChar_t  trigger_HLT[128];

    Int_t   IsoMu20Pass;
    Int_t   IsoTkMu20Pass;

    Float_t beamspot_x;
    Float_t beamspot_y;
    Float_t beamspot_z;
    Float_t beamspot_xwidth;
    Float_t beamspot_ywidth;
    Float_t beamspot_zsigma;
    Float_t beamspot_cov[6];

    UInt_t primvertex_count;
    vector<Float_t> *primvertex_x;
    vector<Float_t> *primvertex_y;
    vector<Float_t> *primvertex_z;
    vector<UInt_t>  *primvertex_isvalid;
    vector<UInt_t>  *primvertex_isfake;
    vector<Float_t> *primvertex_chi2;
    vector<Float_t> *primvertex_ndof;
    vector<Int_t>   *primvertex_ntracks;
    vector<Float_t> *primvertex_cov0;
    vector<Float_t> *primvertex_cov1;
    vector<Float_t> *primvertex_cov2;
    vector<Float_t> *primvertex_cov3;
    vector<Float_t> *primvertex_cov4;
    vector<Float_t> *primvertex_cov5;

    UInt_t muon_count;
    vector<Float_t> *muon_px;
    vector<Float_t> *muon_rochesterPx;
    vector<Float_t> *muon_py;
    vector<Float_t> *muon_rochesterPy;
    vector<Float_t> *muon_pz;
    vector<Float_t> *muon_rochesterPz;
    vector<Float_t> *muon_pt;
    vector<Float_t> *muon_rochesterPt;
    vector<Float_t> *muon_pterror;
    vector<Float_t> *muon_chi2;
    vector<Float_t> *muon_ndof;
    vector<Int_t>   *muon_is_global;
    vector<Int_t>   *muon_is_tracker;
    vector<Int_t>   *muon_is_standalone;
    vector<Int_t>   *muon_is_pf_muon;
    vector<Int_t>   *muon_is_tight_muon;
    vector<Int_t>   *muon_is_medium_muon;
    vector<Int_t>   *muon_is_loose_muon;
    //vector<Int_t>   *muon_innertrack_vtx;
    vector<Float_t> *muon_innertrack_px;
    vector<Float_t> *muon_innertrack_py;
    vector<Float_t> *muon_innertrack_pz;
    //vector<Float_t> *muon_innertrack_outerx;
    //vector<Float_t> *muon_innertrack_outery;
    //vector<Float_t> *muon_innertrack_outerz;
    //vector<Float_t> *muon_innertrack_closestpointx;
    //vector<Float_t> *muon_innertrack_closestpointy;
    //vector<Float_t> *muon_innertrack_closestpointz;
    vector<Float_t> *muon_innertrack_chi2;
    vector<Float_t> *muon_innertrack_ndof;
    vector<Float_t> *muon_innertrack_dxy;
    vector<Float_t> *muon_innertrack_dxyerr;
    vector<Float_t> *muon_innertrack_dz;
    vector<Float_t> *muon_innertrack_dzerr;
    //vector<Float_t> *muon_innertrack_dedxharmonic2;
    vector<Int_t>   *muon_innertrack_charge;
    vector<UChar_t> *muon_innertrack_nhits;
    vector<UChar_t> *muon_innertrack_nmissinghits;
    vector<UChar_t> *muon_innertrack_npixelhits;
    vector<UChar_t> *muon_innertrack_npixellayers;
    vector<UChar_t> *muon_innertrack_nstriplayers;
    vector<Float_t> *muon_outertrack_px;
    vector<Float_t> *muon_outertrack_py;
    vector<Float_t> *muon_outertrack_pz;
    vector<UChar_t> *muon_outertrack_hits;
    vector<UChar_t> *muon_outertrack_missinghits;
    vector<Float_t> *muon_outertrack_chi2;
    vector<Float_t> *muon_outertrack_ndof;
    vector<Float_t> *muon_isolationr3track;
    vector<Int_t>   *muon_isolationr3ntrack;
    vector<Float_t> *muon_isolationr3ecal;
    vector<Float_t> *muon_isolationr3hcal;
    vector<Float_t> *muon_ecalenergy;
    vector<Float_t> *muon_hcalenergy;
    vector<Float_t> *muon_pfisolationr3_sumchargedhadronpt;
    vector<Float_t> *muon_pfisolationr3_sumchargedparticlept;
    vector<Float_t> *muon_pfisolationr3_sumneutralhadronet;
    vector<Float_t> *muon_pfisolationr3_sumphotonet;
    vector<Float_t> *muon_pfisolationr3_sumneutralhadronethighthreshold;
    vector<Float_t> *muon_pfisolationr3_sumphotonethighthreshold;
    vector<Float_t> *muon_pfisolationr3_sumpupt;
    vector<Float_t> *muon_pfisolationr4_sumchargedhadronpt;
    vector<Float_t> *muon_pfisolationr4_sumchargedparticlept;
    vector<Float_t> *muon_pfisolationr4_sumneutralhadronet;
    vector<Float_t> *muon_pfisolationr4_sumphotonet;
    vector<Float_t> *muon_pfisolationr4_sumneutralhadronethighthreshold;
    vector<Float_t> *muon_pfisolationr4_sumphotonethighthreshold;
    vector<Float_t> *muon_pfisolationr4_sumpupt;
    //vector<Float_t> *muon_pfisolationr4_dBrel;
    vector<Int_t>   *muon_charge;
    vector<Int_t>   *muon_numchambers;
    vector<Int_t>   *muon_numchamberswithsegments;
    vector<Int_t>   *muon_numvalidmuonhits;
    vector<Int_t>   *muon_nummatchedstations;
    vector<UInt_t>  *muon_trigger;
    vector<UInt_t>  *muon_trackermuonquality;
    vector<Int_t>   *muon_matches_IsoMu20;
    vector<Int_t>   *muon_matches_IsoTkMu20;

    UInt_t ak4pfchsjet_count;
    vector<Float_t> *ak4pfchsjet_energy;
    vector<Float_t> *ak4pfchsjet_px;
    vector<Float_t> *ak4pfchsjet_py;
    vector<Float_t> *ak4pfchsjet_pz;
    //vector<Float_t> *ak4pfchsjet_pt;
    vector<Float_t> *ak4pfchsjet_area;
    vector<Float_t> *ak4pfchsjet_hadronicenergy;
    vector<Float_t> *ak4pfchsjet_chargedhadronicenergy;
    vector<Float_t> *ak4pfchsjet_emenergy;
    vector<Float_t> *ak4pfchsjet_chargedemenergy;
    vector<Float_t> *ak4pfchsjet_hfemenergy;
    vector<Float_t> *ak4pfchsjet_hfhadronicenergy;
    vector<Float_t> *ak4pfchsjet_electronenergy;
    vector<Float_t> *ak4pfchsjet_muonenergy;
    vector<UInt_t>  *ak4pfchsjet_chargedmulti;
    vector<UInt_t>  *ak4pfchsjet_neutralmulti;
    vector<UInt_t>  *ak4pfchsjet_hfhadronicmulti;
    vector<UInt_t>  *ak4pfchsjet_hfemmulti;
    vector<UInt_t>  *ak4pfchsjet_electronmulti;
    vector<UInt_t>  *ak4pfchsjet_muonmulti;
    vector<Float_t> *ak4pfchsjet_chargeda;
    vector<Float_t> *ak4pfchsjet_chargedb;
    vector<Float_t> *ak4pfchsjet_neutrala;
    vector<Float_t> *ak4pfchsjet_neutralb;
    vector<Float_t> *ak4pfchsjet_alla;
    vector<Float_t> *ak4pfchsjet_allb;
    vector<Float_t> *ak4pfchsjet_chargedfractionmv;
    vector<Float_t> *ak4pfchsjet_energycorr;
    vector<Float_t> *ak4pfchsjet_energycorrunc;
    //vector<Float_t> *ak4pfchsjet_energycorrl7uds;
    //vector<Float_t> *ak4pfchsjet_energycorrl7bottom;
    vector<Int_t>   *ak4pfchsjet_btag;
    vector<Int_t>   *ak4pfchsjet_idLoose;
    vector<Int_t>   *ak4pfchsjet_idTight;
    vector<Int_t>   *ak4pfchsjet_idTightLepVeto;
    vector<UInt_t>  *ak4pfchsjet_trigger;
    vector<Int_t>   *ak4pfchsjet_mcflavour;

    UInt_t electron_count;
    //vector<Int_t>   *electron_vtx;
    vector<Float_t> *electron_px;
    vector<Float_t> *electron_py;
    vector<Float_t> *electron_pz;
    vector<Float_t> *electron_correctedecalenergy;
    vector<Float_t> *electron_trackchi2;
    vector<Float_t> *electron_trackndof;
    //vector<Float_t> *electron_outerx;
    //vector<Float_t> *electron_outery;
    //vector<Float_t> *electron_outerz;
    //vector<Float_t> *electron_closestpointx;
    //vector<Float_t> *electron_closestpointy;
    //vector<Float_t> *electron_closestpointz;
    vector<Float_t> *electron_esuperclusterovertrack;
    vector<Float_t> *electron_eseedclusterovertrack;
    vector<Float_t> *electron_deltaetasuperclustertrack;
    vector<Float_t> *electron_deltaphisuperclustertrack;
    vector<Float_t> *electron_e1x5;
    vector<Float_t> *electron_e2x5;
    vector<Float_t> *electron_e5x5;
    vector<Float_t> *electron_r9;
    vector<Float_t> *electron_sigmaetaeta;
    vector<Float_t> *electron_sigmaietaieta;
    vector<Float_t> *electron_sigmaiphiiphi;
    vector<Float_t> *electron_ehcaloverecaldepth1;
    vector<Float_t> *electron_ehcaloverecaldepth2;
    vector<Float_t> *electron_ehcaltoweroverecaldepth1;
    vector<Float_t> *electron_ehcaltoweroverecaldepth2;
    vector<Float_t> *electron_isolationr3track;
    vector<Float_t> *electron_isolationr3ecal;
    vector<Float_t> *electron_isolationr3hcal;
    vector<Float_t> *electron_isolationr4track;
    vector<Float_t> *electron_isolationr4ecal;
    vector<Float_t> *electron_isolationr4hcal;
    vector<Float_t> *electron_isolationpfr3charged;
    vector<Float_t> *electron_isolationpfr3photon;
    vector<Float_t> *electron_isolationpfr3neutral;
    vector<Int_t>   *electron_charge;
    vector<Float_t> *electron_effectiveArea;

    vector<Int_t>   *electron_nhits;
    vector<Int_t>   *electron_nmissinghits;
    vector<Int_t>   *electron_npixelhits;
    vector<Int_t>   *electron_npixellayers;
    vector<Int_t>   *electron_nstriplayers;
    vector<Int_t>   *electron_nhitsexpected;
    vector<Float_t> *electron_dxy;
    vector<Float_t> *electron_dxyerr;
    vector<Float_t> *electron_dz;
    vector<Float_t> *electron_dzerr;
    //vector<Float_t> *electron_convdist;
    //vector<Float_t> *electron_convdcot;
    //vector<Float_t> *electron_convradius;
    vector<UInt_t>  *electron_gapinfo;
    vector<Float_t> *electron_fbrems;
    vector<Int_t>   *electron_numbrems;
    //vector<UChar_t> *electron_info;
    vector<Int_t>   *electron_cutBasedLoose;
    vector<Int_t>   *electron_cutBasedMedium;
    vector<Int_t>   *electron_cutBasedTight;
    vector<Int_t>   *electron_mvaNonTrigWP90;
    vector<Int_t>   *electron_mvaNonTrigWP80;

    vector<Int_t>   *electron_iselectron;
    vector<Int_t>   *electron_passconversionveto;
    vector<Int_t>   *electron_ecaldrivenseed;
    vector<Int_t>   *electron_trackerdrivenseed;
    vector<UInt_t>  *electron_trigger;
    vector<Float_t> *electron_supercluster_e;
    vector<Float_t> *electron_supercluster_x;
    vector<Float_t> *electron_supercluster_y;
    vector<Float_t> *electron_supercluster_z;
    vector<Float_t> *electron_supercluster_rawe;
    vector<Float_t> *electron_supercluster_phiwidth;
    vector<Float_t> *electron_supercluster_etawidth;
    vector<Int_t>   *electron_supercluster_nbasiccluster;

    UInt_t photon_count;
    vector<Float_t> *photon_px;
    vector<Float_t> *photon_py;
    vector<Float_t> *photon_pz;
    vector<Float_t> *photon_e1x5;
    vector<Float_t> *photon_e2x5;
    vector<Float_t> *photon_e3x3;
    vector<Float_t> *photon_e5x5;
    vector<Float_t> *photon_sigmaietaieta;
    vector<Float_t> *photon_sigmaietaiphi;
    vector<Float_t> *photon_sigmaiphiiphi;
    vector<Float_t> *photon_ehcaloverecaldepth1;
    vector<Float_t> *photon_ehcaloverecaldepth2;
    vector<Float_t> *photon_ehcaltoweroverecaldepth1;
    vector<Float_t> *photon_ehcaltoweroverecaldepth2;
    vector<Float_t> *photon_maxenergyxtal;
    vector<Float_t> *photon_isolationr3track;
    vector<Float_t> *photon_isolationr3trackhollow;
    vector<UInt_t>  *photon_isolationr3ntrack;
    vector<UInt_t>  *photon_isolationr3ntrackhollow;
    vector<Float_t> *photon_isolationr3ecal;
    vector<Float_t> *photon_isolationr3hcal;
    vector<Float_t> *photon_isolationr4track;
    vector<Float_t> *photon_isolationr4trackhollow;
    vector<UInt_t>  *photon_isolationr4ntrack;
    vector<UInt_t>  *photon_isolationr4ntrackhollow;
    vector<Float_t> *photon_isolationr4ecal;
    vector<Float_t> *photon_isolationr4hcal;
    vector<Float_t> *photon_isolationpfr3charged;
    vector<Float_t> *photon_isolationpfr3photon;
    vector<Float_t> *photon_isolationpfr3neutral;
    //vector<Float_t> *photon_isolationpfr4charged;
    //vector<Float_t> *photon_isolationpfr4photon;
    //vector<Float_t> *photon_isolationpfr4neutral;
    //vector<Float_t> *photon_isolationpfr4noscfootprintcharged;
    //vector<Float_t> *photon_isolationpfr4noscfootprintphoton;
    //vector<Float_t> *photon_isolationpfr4noscfootprintneutral;
    vector<Float_t> *photon_supercluster_e;
    vector<Float_t> *photon_supercluster_x;
    vector<Float_t> *photon_supercluster_y;
    vector<Float_t> *photon_supercluster_z;
    vector<Float_t> *photon_supercluster_rawe;
    vector<Float_t> *photon_supercluster_phiwidth;
    vector<Float_t> *photon_supercluster_etawidth;
    vector<Int_t>   *photon_supercluster_nbasiccluster;
    //vector<UChar_t> *photon_info;
    vector<Int_t>   *photon_isphoton;
    vector<Int_t>   *photon_hasconversiontracks;
    vector<Int_t>   *photon_haspixelseed;
    vector<Int_t>   *photon_passelectronveto;
    vector<Int_t>   *photon_ispfphoton;
    vector<UInt_t>  *photon_gapinfo;
    vector<UInt_t>  *photon_trigger;
    //vector<UInt_t>  *photon_conversionbegin;

    UInt_t tau_count;
    //UInt_t tau_charged_count;
    vector<Float_t> *tau_px;
    vector<Float_t> *tau_py;
    vector<Float_t> *tau_pz;
    vector<Float_t> *tau_isolationneutralspt;
    vector<UInt_t>  *tau_isolationneutralsnum;
    vector<Float_t> *tau_isolationchargedpt;
    vector<UInt_t>  *tau_isolationchargednum;
    vector<Float_t> *tau_isolationgammapt;
    vector<UInt_t>  *tau_isolationgammanum;
    vector<Int_t>   *tau_charge;
    vector<UInt_t>  *tau_disc;
    //vector<Float_t> *tau_emfraction;
    //vector<Float_t> *tau_hcaltotoverplead;
    //vector<Float_t> *tau_hcal3x3overplead;
    //vector<Float_t> *tau_ecalstripsumeoverplead;
    //vector<Float_t> *tau_bremsrecoveryeoverplead;
    //vector<Float_t> *tau_calocomp;
    //vector<Float_t> *tau_segcomp;
    vector<UInt_t>  *tau_trigger;
    //vector<Float_t> *tau_ak4pfjet_e;
    //vector<Float_t> *tau_ak4pfjet_px;
    //vector<Float_t> *tau_ak4pfjet_py;
    //vector<Float_t> *tau_ak4pfjet_pz;
    //vector<Float_t> *tau_ak4pfjet_hadronicenergy;
    //vector<Float_t> *tau_ak4pfjet_chargedhadronicenergy;
    //vector<Float_t> *tau_ak4pfjet_emenergy;
    //vector<Float_t> *tau_ak4pfjet_chargedemenergy;
    //vector<UInt_t>  *tau_ak4pfjet_chargedmulti;
    //vector<UInt_t>  *tau_ak4pfjet_neutralmulti;
    //vector<UInt_t>  *tau_ak4pfjet_trigger;
    //vector<UInt_t>  *tau_chargedbegin;
    //vector<Float_t> *tau_charged_px;
    //vector<Float_t> *tau_charged_py;
    //vector<Float_t> *tau_charged_pz;
    //vector<Float_t> *tau_charged_outerx;
    //vector<Float_t> *tau_charged_outery;
    //vector<Float_t> *tau_charged_outerz;
    //vector<Float_t> *tau_charged_closestpointx;
    //vector<Float_t> *tau_charged_closestpointy;
    //vector<Float_t> *tau_charged_closestpointz;
    //vector<Float_t> *tau_charged_chi2;
    //vector<Float_t> *tau_charged_ndof;
    //vector<Float_t> *tau_charged_dxy;
    //vector<Float_t> *tau_charged_dxyerr;
    //vector<Float_t> *tau_charged_dz;
    //vector<Float_t> *tau_charged_dzerr;
    //vector<Float_t> *tau_charged_dedxharmonic2;
    //vector<Int_t>   *tau_charged_charge;
    //vector<UChar_t> *tau_charged_nhits;
    //vector<UChar_t> *tau_charged_nmissinghits;
    //vector<UChar_t> *tau_charged_npixelhits;
    //vector<UChar_t> *tau_charged_npixellayers;
    //vector<UChar_t> *tau_charged_nstriplayers;

    Float_t event_rho;
    //Float_t ak4pfjet_sigma;

    //Float_t pfmet_ex;
    //Float_t pfmet_ey;
    //Float_t pfmet_et;
    vector<Float_t> *pfmettype1_ex;
    vector<Float_t> *pfmettype1_ey;
    vector<Float_t> *pfmettype1_et;
    //Float_t pfmetpuppitype1_ex;
    //Float_t pfmetpuppitype1_ey;
    //Float_t pfmettype0type1_ex;
    //Float_t pfmettype0type1_ey;

    //Generator Information
    Float_t genweight;
    Float_t genid1;
    Float_t genx1;
    Float_t genid2;
    Float_t genx2;
    Float_t genScale;

    Int_t numpileupinteractionsminus;
    Int_t numpileupinteractions;
    Int_t numpileupinteractionsplus;
    Float_t numtruepileupinteractions;

    //Float_t genmetcalo_ex;
    //Float_t genmetcalo_ey;
    //Float_t genmettrue_ex;
    //Float_t genmettrue_ey;

    vector<Float_t> *genmet_ex;
    vector<Float_t> *genmet_ey;
    //UInt_t genak4jet_count;
    //Float_t genak4jet_e[M_genjetmaxcount];
    //Float_t genak4jet_px[M_genjetmaxcount];
    //Float_t genak4jet_py[M_genjetmaxcount];
    //Float_t genak4jet_pz[M_genjetmaxcount];
    //Float_t genak4jet_einvisible[M_genjetmaxcount];
    //Int_t genak4jet_flavour[M_genjetmaxcount];
    //UInt_t genak4jet_info[M_genjetmaxcount];

    UInt_t genjet_count;
    vector<Float_t> *genjet_e;
    vector<Float_t> *genjet_px;
    vector<Float_t> *genjet_py;
    vector<Float_t> *genjet_pz;
    vector<Float_t> *genjet_einvisible;

    UInt_t genparticles_count;
    Float_t genparticles_e[M_genparticlesmaxcount];
    Float_t genparticles_px[M_genparticlesmaxcount];
    Float_t genparticles_py[M_genparticlesmaxcount];
    Float_t genparticles_pz[M_genparticlesmaxcount];
    Float_t genparticles_vx[M_genparticlesmaxcount];
    Float_t genparticles_vy[M_genparticlesmaxcount];
    Float_t genparticles_vz[M_genparticlesmaxcount];
    Int_t genparticles_pdgid[M_genparticlesmaxcount];
    Int_t genparticles_status[M_genparticlesmaxcount];
    Int_t genparticles_indirectmother[M_genparticlesmaxcount];
    UInt_t genparticles_info[M_genparticlesmaxcount];

    UInt_t genparticles_motherbeg[M_genallparticlesmaxcount];
    UInt_t genparticles_daughterbeg[M_genallparticlesmaxcount];
    UInt_t genparticlesmother_count;
    UInt_t genparticles_mothers[M_genmotherdaughtermaxcount];
    UInt_t genparticlesdaughter_count;
    UInt_t genparticles_daughters[M_genmotherdaughtermaxcount];

    //UInt_t genallparticles_count;
    //Float_t genallparticles_e[M_genallparticlesmaxcount];
    //Float_t genallparticles_px[M_genallparticlesmaxcount];
    //Float_t genallparticles_py[M_genallparticlesmaxcount];
    //Float_t genallparticles_pz[M_genallparticlesmaxcount];
    //Float_t genallparticles_vx[M_genallparticlesmaxcount];
    //Float_t genallparticles_vy[M_genallparticlesmaxcount];
    //Float_t genallparticles_vz[M_genallparticlesmaxcount];
    //Int_t genallparticles_pdgid[M_genallparticlesmaxcount];
    //Int_t genallparticles_status[M_genallparticlesmaxcount];
    //UInt_t genallparticles_motherbeg[M_genallparticlesmaxcount];
    //UInt_t genallparticles_daughterbeg[M_genallparticlesmaxcount];

    //UInt_t genallparticlesmother_count;
    //UInt_t genallparticles_mothers[M_genmotherdaughtermaxcount];

    //UInt_t genallparticlesdaughter_count;
    //UInt_t genallparticles_daughters[M_genmotherdaughtermaxcount];

public:
    Analyse(int argc = 0, char **argv = 0, bool batchmode = false);
    virtual ~Analyse();
    Int_t AddFile(string filename);
    Int_t AddDir(string path);
    void Batch_Prepare(bool simple = false); //will run on a file based splitting. dublicate checking will not work.
//do not use this unless you are running on 2011 MC where the lumi based splitting will not work.
    Int_t Batch_MyId() const {
        return(batch_myid);
    }
    Int_t Batch_NumJobs() const {
        return(batch_numjobs);
    }
    Long64_t GetNumAddedEvents() const {
        return(filelimits[filelimits.size()-1]);
    };
    string GetCurrentFileName() const {
        if(currentfile != 0) {
            return(currentfile->GetName());
        } else {
            return("NO FILE OPENED");
        }
    }

    // general event information
    Double_t Number() const { return(event_nr); }
    UInt_t   Run() const { return(event_run); }
    UInt_t   LumiBlock() const { return(event_luminosityblock); }
    UInt_t   TimeUnix() const { return(event_timeunix); }
    UInt_t   TimeMicroSec() const { return(event_timemicrosec); }
    Double_t Rho() const { return(event_rho); }
    bool     IsData() const { return(isdata); }
    bool     IsMC() const { return(!(isdata)); }
    //Double_t Sigma() const { return(ak4pfjet_sigma); }

    // RECO-level information
    void LoadBeamSpot(bool select = true);
    BeamSpot GetBeamSpot() const;

    void LoadMuons(bool select = true);
    Muon Muons(UInt_t n, int correction = 0) const;
    UInt_t NumMuons() const { return(muon_count); }

    void LoadElectrons(bool select = true);
    Electron Electrons(UInt_t n, int correction = 0) const;
    UInt_t NumElectrons() const { return(electron_count); }

    void LoadPhotons(bool select = true);
    Photon Photons(UInt_t n) const;
    UInt_t NumPhotons() const { return(photon_count); }

    void LoadTaus(bool select = true);
    Tau Taus(UInt_t n) const;
    UInt_t NumTaus() const { return(tau_count); }

    void LoadMET(bool select = true);
    //TLorentzVector PFMET() const;
    TLorentzVector PFMETTYPE1() const;
    //TLorentzVector PFMETPUPPITYPE1() const;
    //TLorentzVector PFMETTYPE0TYPE1() const;

    void LoadAK4PFCHSJets(bool select = true);
    Jet AK4PFCHSJets(UInt_t n) const;
    UInt_t NumAK4PFCHSJets() const {
        return(ak4pfchsjet_count);
    }

    void LoadPrimVertices(bool select = true);
    Vertex PrimVertices(UInt_t n) const;
    UInt_t NumPrimVertices() const {
        return(primvertex_count);
    }

    // generator-level information
    void LoadGenInfo(bool select = true);
    Double_t GenWeight() const {
        return(genweight);
    }
    Double_t GenId1() const {
        return(genid1);
    }
    Double_t Genx1() const {
        return(genx1);
    }
    Double_t GenId2() const {
        return(genid2);
    }
    Double_t Genx2() const {
        return(genx2);
    }
    Double_t GenScale() const {
        return(genScale);
    }

    Int_t NumPileUpInteractionsMinus() const {
        return(numpileupinteractionsminus);
    }
    Int_t NumPileUpInteractions() const {
        return(numpileupinteractions);
    }
    Int_t NumPileUpInteractionsPlus() const {
        return(numpileupinteractionsplus);
    }
    Float_t NumTruePileUpInteractions() const {
        return(numtruepileupinteractions);
    }
    void UsePileUpInfo();
    //Double_t GetPileUpWeight(Double_t meaninteractions) const;
    void SetPileUpScale(double puscale) {
        pileupscale = puscale;
    }
    Double_t GetPileUpMaxWeight(vector<Double_t> &datadist) const;
    Double_t GetPileUpWeight(vector<Double_t> &datadist) const;
    void UsePrimVertexInfo();
    //Double_t GetPrimVertexWeight(Double_t meaninteractions) const;
    Double_t GetPrimVertexMaxWeight(vector<Double_t> &datadist) const;
    Double_t GetPrimVertexWeight(vector<Double_t> &datadist) const;
    Int_t NumGoodPrimVertices() const;

    //void LoadAllGenParticles(bool select = true);
    //GenParticle AllGenParticles(UInt_t n) const;
    //UInt_t NumAllGenParticles() const {
    //    return(genallparticles_count);
    //}

    void LoadGenParticles(bool select = true);
    GenLightParticle GenParticles(UInt_t n) const;
    UInt_t NumGenParticles() const {
        return(genparticles_count);
    }

    void LoadGenJets(bool select = true);
    GenJet GenJets(UInt_t n) const;
    UInt_t NumGenJets() const {
        //return(genak4jet_count);
        return(genjet_count);
    }

    //TLorentzVector GenMETCalo() const;
    //TLorentzVector GenMETTrue() const;
    TLorentzVector GenMET() const;

    //virtual dummies:
    //executed for each event in the loop
    virtual Int_t AnalyseEvent() {
        return(1);
    };
    //executed once at the beginning of a the loop
    virtual void BeginLoop() {};
    //executed once at the end of a loop
    virtual void EndLoop() {};

    //use this to start your analysis
    Long64_t Loop(Long64_t start = 0, Long64_t end = -1);
    //Number of Processed Events
    Long64_t Processed() const {
        return(processed);
    }

    //setting output: -1 no output, the p-th events are printed
    void SetPrintInfo(Int_t p = -1) {
        printinfo = p;
    }

    //tool to check for duplicated events
    void EnableDuplicateCheck(bool switchval = true) {
        duplicatecheck = switchval;
    }
    UInt_t CheckDuplicate() {
        return(eventlist[Run()][LumiBlock()][Number()]);
    }

    //Skiming
    Int_t PrepareSkimming(string filename);
    Int_t SkimEvent();

    //trigger information
    void LoadTrigger(bool select = true);
    TriggerSelection *AddTriggerSelection(string id, vector<string> triggernames, bool useprescaled = false);
    TriggerSelection *GetTriggerSelection(string id);
    bool GetL1Trigger(UInt_t bit) const;
    bool GetL1TriggerBits(UInt_t bit) const;
    bool GetHLTrigger(UInt_t index) const;
    Int_t GetHLTrigger(vector<string> triggernames) const; //looks if one of the tiggers fired. Only unprescaled trigger are considered. returns -1 if non of the trigger is unprescaled, 0 non fired, 1 at least one fired.
    Int_t GetHLTriggerIndex(string triggername) const;
    string GetHLTriggerName(UInt_t index) const;
    Int_t GetHLTPrescale(UInt_t triggerindex) const;
    Int_t GetNumHLTriggers() const;

    //Lumi calculation
    void AddLumiFile(string filename, string dir = "");
    Int_t IsLumiAvailable() const;
    Double_t GetInstLumi() const;
    Double_t GetAvgPU() const;
    Double_t GetLumi(Int_t format = 0);
    //Double_t GetLumiBlockLumi(); //Lumi in dataset up to current event.
    void PrintPrescaleInfo(string triggername);
    void PrintPrescaleInfoB(string triggername);
    void PrintLumiOfRuns();
    void PrintLumiOfLumiSectionsInRun(UInt_t runnumber);
    //very special don't use!
    Int_t SetLumi(UInt_t run, UInt_t block, Float_t lumival, Float_t avgpu = -1);
    void ResetLumiValues();
    void WriteLumiFile(string filename);
    //Use JSON filter
    bool LoadJSON(string filename);


// v. hacky and should be removed asap
    Bool_t passesisomu20() const { return(IsoMu20Pass); }
    Bool_t passesisotkmu20() const { return(IsoTkMu20Pass); }

    bool EventPassesHLT(std::vector<string> hltnames) const;
};

extern Analyse *GLAN;
Long64_t mem_usage();
#endif
