#include "Analyse.h"
#include "DataForm.h"
#ifdef USE_MPI
#include "mpi.h"
MPI_Status mpistatus;
#endif
Analyse *GLAN = 0;


//____________________________________________________________________________________
Analyse::Analyse(int argc, char **argv, bool batchmode) :
    isdata(isdata),
    currentfile(0),
    currentloadtype(-1),
    printinfo(1),
    processed(0),
    loadbeamspot(0),
    loadmuons(0),
    loadelectrons(0),
    loadphotons(0),
    loadtaus(0),
    loadmet(0),
    loadak4pfchsjets(0),
    loadtracks(0),
    loadprimvertices(0),
    loadtrigger(0),
    loadgeninfo(0),
    loadgenparticles(0),
    loadgenak4jets(0),
    loadallgenparticles(0),
    duplicatecheck(false),
    skimtree(0),
    skimfilename("None"),
    skimfile(0),
    lumicalculation(false),
    minRun(1000000000),
    minLumi(100000000),
    maxRun(0),
    maxLumi(0),
    jsonfilter(false),
    batch_myid(-1),
    batch_numjobs(-1),
    batch_emptyjob(false),
    usepileupinfo(false),
    pileupscale(1.),
    pileUpDistMinus(200, 0.),
    pileUpDist(200, 0.),
    pileUpDistPlus(200, 0.),
    useprimvertexinfo(false),
    primVertexDist(200, 0.),
    primvertex_count(0),
    muon_count(0),
    electron_count(0),
    photon_count(0),
    tau_count(0),
    //tau_charged_count(0),
    genparticles_count(0)
    //genallparticles_count(0),
    //genallparticlesmother_count(0),
    //genallparticlesdaughter_count(0) 
{
    GLAN = this;
#ifdef USE_MPI
    if(batchmode) {
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &batch_numjobs);
        MPI_Comm_rank(MPI_COMM_WORLD, &batch_myid);
        cout << "INIT " << batch_numjobs << " " << batch_myid << endl;
    }
#endif
    analysisrandom = new TRandom3;
}

//____________________________________________________________________________________
Analyse::~Analyse() {
    if(skimfilename != string("None")) {
        TFile *file = skimtree->GetCurrentFile();
        file->Write();
        file->Close();
        string dirname;
        string outfilename;
        string lumifilename;

        UInt_t slashpos = skimfilename.find_last_of("/");
        if(slashpos != skimfilename.size()) {
            dirname = skimfilename.substr(0, slashpos);
            outfilename = skimfilename.substr(slashpos+1);
        } else {
            outfilename = skimfilename;
        }

        lumifilename = string("LUMI_INFO.root");
        cout << "Writing LUMI-file: " << lumifilename << endl;

        WriteLumiFile(dirname + string("/") + lumifilename);
    }

#ifdef USE_MPI
    if(Batch_MyId() != -1) {
        MPI_Finalize();
    }
#endif
}

//____________________________________________________________________________________
void Analyse::SetLoad() {
    tree->SetBranchStatus("*", false);
    //tree->SetBranchStatus("errors", true);
    tree->SetBranchStatus("event_*", true);
    tree->SetBranchStatus("numpileupinteractions*", true);
    tree->SetBranchStatus("numtruepileupinteractions", true);
    tree->SetBranchStatus("event_rho", true);
    tree->SetBranchStatus("isdata", true);
    tree->SetBranchStatus("*Pass", true);
    //tree->SetBranchStatus("ak4pfjet_sigma", true);

    if(loadbeamspot == 1) {
        tree->SetBranchStatus("beamspot_*", true);
    } else {
        tree->SetBranchStatus("beamspot_*", false);
    }

    if(loadmuons == 1) {
        tree->SetBranchStatus("muon_*", true);
    } else {
        tree->SetBranchStatus("muon_*", false);
    }

    if(loadelectrons == 1) {
        tree->SetBranchStatus("electron_*", true);
    } else {
        tree->SetBranchStatus("electron_*", false);
    }

    if(loadphotons == 1) {
        tree->SetBranchStatus("photon_*", true);
    } else {
        tree->SetBranchStatus("photon_*", false);
    }

    if(loadtaus == 1) {
        tree->SetBranchStatus("tau_*", true);
    } else {
        tree->SetBranchStatus("tau_*", false);
    }

    if(loadmet == 1) {
        tree->SetBranchStatus("pfmet*", true);
    } else {
        tree->SetBranchStatus("pfmet*", false);
    }

    if(loadak4pfchsjets == 1) {
        tree->SetBranchStatus("ak4pfchsjet_*", true);
    } else {
        tree->SetBranchStatus("ak4pfchsjet_*", false);
    }

/*
    if(loadtracks == 1) {
        tree->SetBranchStatus("track_*", true);
    } else {
        tree->SetBranchStatus("track_*", false);
    }
    if(loadsuperclusters > 0) {
        tree->SetBranchStatus("supercluster*", true);
        if(loadsuperclusters != 3 && loadsuperclusters != 4) {
            tree->SetBranchStatus("supercluster_basiccluster_hit*", false);
            tree->SetBranchStatus("supercluster_escluster_hit*", false);
        }
        if(loadsuperclusters != 2 && loadsuperclusters != 4) {
            tree->SetBranchStatus("supercluster_basiccluster*", false);
            tree->SetBranchStatus("supercluster_escluster*", false);
        }
    } else {
        tree->SetBranchStatus("supercluster*", false);
    }
*/
    if(loadprimvertices == 1) {
        tree->SetBranchStatus("primvertex_*", true);
    } else {
        tree->SetBranchStatus("primvertex_*", false);
    }

    //if(loadtrigger == 1) {
    //    tree->SetBranchStatus("trigger_*", true);
    //} else {
    //    tree->SetBranchStatus("trigger_*", false);
    //}

    if(loadgeninfo == 1) {
        tree->SetBranchStatus("genweight", true);
        tree->SetBranchStatus("genid1", true);
        tree->SetBranchStatus("genx1", true);
        tree->SetBranchStatus("genid2", true);
        tree->SetBranchStatus("genx2", true);
        tree->SetBranchStatus("genScale", true);
        tree->SetBranchStatus("genmet*", true);
    } else {
        tree->SetBranchStatus("genweight", false);
        tree->SetBranchStatus("genid1", false);
        tree->SetBranchStatus("genx1", false);
        tree->SetBranchStatus("genid2", false);
        tree->SetBranchStatus("genx2", false);
        tree->SetBranchStatus("genScale", false);
        tree->SetBranchStatus("genmet*", false);
    }

    //if(loadallgenparticles == 1) {
    //    tree->SetBranchStatus("genall*", true);
    //} else {
    //    tree->SetBranchStatus("genall*", false);
    //}

    if(loadgenparticles == 1) {
        tree->SetBranchStatus("genparticle*", true);
    } else {
        tree->SetBranchStatus("genparticle*", false);
    }

    if(loadgenak4jets == 1) {
        tree->SetBranchStatus("genak4jet*", true);
    } else {
        tree->SetBranchStatus("genak4jet*", false);
    }
}

//____________________________________________________________________________________
bool sortac1b(string a, string b) {
    vector<string> strs;
    boost::split(strs, a, boost::is_any_of("_"));
    int na = atoi(strs[1].c_str());
    boost::split(strs, b, boost::is_any_of("_"));
    int nb = atoi(strs[1].c_str());
    return(na < nb);
}

//____________________________________________________________________________________
bool sortroot(string a, string b) {
    vector<string> strs;
    boost::split(strs, a, boost::is_any_of("_."));
    int na = atoi(strs[strs.size()-2].c_str());
    boost::split(strs, b, boost::is_any_of("_."));
    int nb = atoi(strs[strs.size()-2].c_str());
    return(na < nb);
}

//____________________________________________________________________________________
bool sortfiles(string a, string b) {
    if(a.find("AC1B") != string::npos) {
        vector<string> strs;
        boost::split(strs, a, boost::is_any_of("_"));
        int na = atoi(strs[strs.size()-3].c_str());
        boost::split(strs, b, boost::is_any_of("_"));
        int nb = atoi(strs[strs.size()-3].c_str());
        return(na < nb);
    } else {
        vector<string> strs;
        boost::split(strs, a, boost::is_any_of("_."));
        int na = atoi(strs[strs.size()-2].c_str());
        boost::split(strs, b, boost::is_any_of("_."));
        int nb = atoi(strs[strs.size()-2].c_str());
        return(na < nb);
    }
}

//____________________________________________________________________________________
Int_t Analyse::AddDir(string path) {
    if(Batch_MyId() == 0 || Batch_MyId() == -1) {
        DIR *dir = opendir(path.c_str());
        if(dir == 0) {
            cerr << "ERROR AddDir: could not open directory!" << endl;
            exit(0);
        }
        dirent *currfile;
        vector<string> filenames;
        size_t orig = 0;
        while((currfile = readdir(dir)) != 0) {
            string currname(currfile->d_name);
            if(currname.find(".root") != string::npos) {
                if(currname.find("LUMI") != string::npos) {
                    AddLumiFile(path+string("/")+currname);
                    cout << "Added LUMI-file: " << currname << endl;
                    continue;
                }
                if(Batch_MyId() != 0 || filenames.size() == 0) { //In batchmode only one file is loaded for job 0
                    if(currname.find("AC1B") != string::npos) {
                        orig++;
                    }
                    filenames.push_back(currname);
                }
            }
        }

        if(filenames.size() == orig) {
            sort(filenames.begin(), filenames.end(), sortac1b);
        } else {
            sort(filenames.begin(), filenames.end(), sortroot);

        }

        for(size_t i = 0 ; i < filenames.size() ; i++) {
            cout << i+1 << ". Add: " << filenames[i] << endl;
            AddFile(path+string("/")+filenames[i]);
        }

        closedir(dir);
        return(filenames.size());
    }
    return(-1);
}

//____________________________________________________________________________________
Int_t Analyse::AddFile(string filename) {
    TTree *test = 0;
    TFile *file = new TFile(filename.c_str(), "READ");
    if(file->IsZombie()) {
        cerr << "ERROR AddFile: " << filename << " is not a valid file." << endl;
    }
    file->GetObject("makeroottree/AC1B", test);
    if(test != 0) {
        if(treename.size() == 0) {
            treename = string("makeroottree/AC1B");
        }
        filenames.push_back(filename);
        if(filelimits.size() == 0) {
            filelimits.push_back(test->GetEntries());
        } else {
            filelimits.push_back(test->GetEntries() + filelimits[filelimits.size()-1]);
        }
        file->Close();
        return(1);
    }

    file->GetObject("AC1B", test);
    if(test != 0) {
        if(treename.size() == 0) {
            treename = string("AC1B");
        }
        filenames.push_back(filename);
        if(filelimits.size() == 0) {
            filelimits.push_back(test->GetEntries());
        } else {
            filelimits.push_back(test->GetEntries() + filelimits[filelimits.size()-1]);
        }
        file->Close();
        return(1);
    }

    cerr <<"ERROR Addfile: " << filename << " no AC1B tree found." << endl;
    return(-1);

}

//____________________________________________________________________________________
void Analyse::GetEvent(Long64_t num, Int_t loadtype) {
    bool changed = false;
    if(!(num >= currentmin && num < currentmax)) {
        if(currentfile != 0) {
            delete currentfile;
        }
        for(UInt_t i = 0 ; i < filelimits.size() ; i++) {
            if(num < filelimits[i]) {
                currentfile = new TFile(filenames[i].c_str(), "READ");
                currentfile->GetObject(treename.c_str(), tree);
                currentmax = filelimits[i];
                if(i == 0) {
                    currentmin = 0;
                } else {
                    currentmin = filelimits[i-1];
                }
                Load();

                break;
            }
        }
        changed = true;
    }

    if(changed || currentloadtype != loadtype) {
        if(loadtype == 0) {
            SetLoad();
        } else if(loadtype == 3) {
            tree->SetBranchStatus("*", false);
            tree->SetBranchStatus("isdata", true);
            tree->SetBranchStatus("event_nr", true);
            tree->SetBranchStatus("event_luminosityblock", true);
            tree->SetBranchStatus("event_run", true);
        } else if(loadtype == 1) {
            tree->SetBranchStatus("*", false);
            tree->SetBranchStatus("numtruepileupinteractions", true);
            tree->SetBranchStatus("numpileupinteractions", true);
            tree->SetBranchStatus("isdata", true);
            tree->SetBranchStatus("event_nr", true);
            tree->SetBranchStatus("event_luminosityblock", true);
            tree->SetBranchStatus("event_run", true);
        } else if(loadtype == 2) {
            tree->SetBranchStatus("*", false);
            tree->SetBranchStatus("primvertex_*", true);
            tree->SetBranchStatus("isdata", true);
            tree->SetBranchStatus("event_nr", true);
            tree->SetBranchStatus("event_luminosityblock", true);
            tree->SetBranchStatus("event_run", true);
        } else {
            cerr << "ERROR GetEvent: This is an unknown loadtype." << endl;
        }
        currentloadtype = loadtype;
    }
    tree->GetEntry(num - currentmin);
}

//____________________________________________________________________________________
void Analyse::Load() {
    tree->SetBranchAddress("isdata", &isdata);
    //tree->SetBranchAddress("errors", &errors);
    tree->SetBranchAddress("event_nr", &event_nr);
    tree->SetBranchAddress("event_luminosityblock", &event_luminosityblock);
    tree->SetBranchAddress("event_run", &event_run);
    tree->SetBranchAddress("event_timeunix", &event_timeunix);
    tree->SetBranchAddress("event_timemicrosec", &event_timemicrosec);
    //tree->SetBranchAddress("trigger_level1bits", trigger_level1bits);
    //tree->SetBranchAddress("trigger_level1", trigger_level1);
    //tree->SetBranchAddress("trigger_HLT", trigger_HLT);

    tree->SetBranchAddress("IsoMu20Pass", &IsoMu20Pass);
    tree->SetBranchAddress("IsoTkMu20Pass", &IsoTkMu20Pass);

    tree->SetBranchAddress("beamspot_x", &beamspot_x);
    tree->SetBranchAddress("beamspot_y", &beamspot_y);
    tree->SetBranchAddress("beamspot_z", &beamspot_z);
    tree->SetBranchAddress("beamspot_xwidth", &beamspot_xwidth);
    tree->SetBranchAddress("beamspot_ywidth", &beamspot_ywidth);
    tree->SetBranchAddress("beamspot_zsigma", &beamspot_zsigma);
    tree->SetBranchAddress("beamspot_cov", &beamspot_cov);

    tree->SetBranchAddress("primvertex_count", &primvertex_count);
    tree->SetBranchAddress("primvertex_x", &primvertex_x);
    tree->SetBranchAddress("primvertex_y", &primvertex_y);
    tree->SetBranchAddress("primvertex_z", &primvertex_z);
    //tree->SetBranchAddress("primvertex_info", primvertex_info);
    tree->SetBranchAddress("primvertex_isvalid", &primvertex_isvalid);
    tree->SetBranchAddress("primvertex_isfake", &primvertex_isfake);
    tree->SetBranchAddress("primvertex_chi2", &primvertex_chi2);
    tree->SetBranchAddress("primvertex_ndof", &primvertex_ndof);
    //tree->SetBranchAddress("primvertex_ptq", primvertex_ptq);
    tree->SetBranchAddress("primvertex_ntracks", &primvertex_ntracks);
    //tree->SetBranchAddress("primvertex_cov", primvertex_cov);
    tree->SetBranchAddress("primvertex_cov0", &primvertex_cov0);
    tree->SetBranchAddress("primvertex_cov1", &primvertex_cov1);
    tree->SetBranchAddress("primvertex_cov2", &primvertex_cov2);
    tree->SetBranchAddress("primvertex_cov3", &primvertex_cov3);
    tree->SetBranchAddress("primvertex_cov4", &primvertex_cov4);
    tree->SetBranchAddress("primvertex_cov5", &primvertex_cov5);

    tree->SetBranchAddress("muon_count", &muon_count);
    tree->SetBranchAddress("muon_px", &muon_px);
    tree->SetBranchAddress("muon_py", &muon_py);
    tree->SetBranchAddress("muon_pz", &muon_pz);
    tree->SetBranchAddress("muon_pt", &muon_pt);
    tree->SetBranchAddress("muon_rochesterPx", &muon_rochesterPx);
    tree->SetBranchAddress("muon_rochesterPy", &muon_rochesterPy);
    tree->SetBranchAddress("muon_rochesterPz", &muon_rochesterPz);
    tree->SetBranchAddress("muon_rochesterPt", &muon_rochesterPt);
    tree->SetBranchAddress("muon_pterror", &muon_pterror);
    tree->SetBranchAddress("muon_chi2", &muon_chi2);
    tree->SetBranchAddress("muon_ndof", &muon_ndof);
    tree->SetBranchAddress("muon_is_global", &muon_is_global);
    tree->SetBranchAddress("muon_is_tracker", &muon_is_tracker);
    tree->SetBranchAddress("muon_is_standalone", &muon_is_standalone);
    tree->SetBranchAddress("muon_is_pf_muon", &muon_is_pf_muon);
    tree->SetBranchAddress("muon_is_tight_muon", &muon_is_tight_muon);
    tree->SetBranchAddress("muon_is_medium_muon", &muon_is_medium_muon);
    tree->SetBranchAddress("muon_is_loose_muon", &muon_is_loose_muon);
    //tree->SetBranchAddress("muon_innertrack_vtx", &muon_innertrack_vtx);
    tree->SetBranchAddress("muon_innertrack_px", &muon_innertrack_px);
    tree->SetBranchAddress("muon_innertrack_py", &muon_innertrack_py);
    tree->SetBranchAddress("muon_innertrack_pz", &muon_innertrack_pz);
    //tree->SetBranchAddress("muon_innertrack_outerx", &muon_innertrack_outerx);
    //tree->SetBranchAddress("muon_innertrack_outery", &muon_innertrack_outery);
    //tree->SetBranchAddress("muon_innertrack_outerz", &muon_innertrack_outerz);
    //tree->SetBranchAddress("muon_innertrack_closestpointx", &muon_innertrack_closestpointx);
    //tree->SetBranchAddress("muon_innertrack_closestpointy", &muon_innertrack_closestpointy);
    //tree->SetBranchAddress("muon_innertrack_closestpointz", &muon_innertrack_closestpointz);
    tree->SetBranchAddress("muon_innertrack_chi2", &muon_innertrack_chi2);
    tree->SetBranchAddress("muon_innertrack_ndof", &muon_innertrack_ndof);
    tree->SetBranchAddress("muon_innertrack_dxy", &muon_innertrack_dxy);
    tree->SetBranchAddress("muon_innertrack_dxyerr", &muon_innertrack_dxyerr);
    tree->SetBranchAddress("muon_innertrack_dz", &muon_innertrack_dz);
    tree->SetBranchAddress("muon_innertrack_dzerr", &muon_innertrack_dzerr);
    //tree->SetBranchAddress("muon_innertrack_dedxharmonic2", &muon_innertrack_dedxharmonic2);
    tree->SetBranchAddress("muon_innertrack_charge", &muon_innertrack_charge);
    tree->SetBranchAddress("muon_innertrack_nhits", &muon_innertrack_nhits);
    tree->SetBranchAddress("muon_innertrack_npixelhits", &muon_innertrack_npixelhits);
    tree->SetBranchAddress("muon_innertrack_nmissinghits", &muon_innertrack_nmissinghits);
    tree->SetBranchAddress("muon_innertrack_npixellayers", &muon_innertrack_npixellayers);
    tree->SetBranchAddress("muon_innertrack_nstriplayers", &muon_innertrack_nstriplayers);
    tree->SetBranchAddress("muon_outertrack_px", &muon_outertrack_px);
    tree->SetBranchAddress("muon_outertrack_py", &muon_outertrack_py);
    tree->SetBranchAddress("muon_outertrack_pz", &muon_outertrack_pz);
    tree->SetBranchAddress("muon_outertrack_hits", &muon_outertrack_hits);
    tree->SetBranchAddress("muon_outertrack_missinghits", &muon_outertrack_missinghits);
    tree->SetBranchAddress("muon_outertrack_chi2", &muon_outertrack_chi2);
    tree->SetBranchAddress("muon_outertrack_ndof", &muon_outertrack_ndof);
    tree->SetBranchAddress("muon_isolationr3track", &muon_isolationr3track);
    tree->SetBranchAddress("muon_isolationr3ntrack", &muon_isolationr3ntrack);
    tree->SetBranchAddress("muon_isolationr3ecal", &muon_isolationr3ecal);
    tree->SetBranchAddress("muon_isolationr3hcal", &muon_isolationr3hcal);
    tree->SetBranchAddress("muon_pfisolationr3_sumchargedhadronpt", &muon_pfisolationr3_sumchargedhadronpt);
    tree->SetBranchAddress("muon_pfisolationr3_sumchargedparticlept", &muon_pfisolationr3_sumchargedparticlept);
    tree->SetBranchAddress("muon_pfisolationr3_sumneutralhadronet", &muon_pfisolationr3_sumneutralhadronet);
    tree->SetBranchAddress("muon_pfisolationr3_sumphotonet", &muon_pfisolationr3_sumphotonet);
    tree->SetBranchAddress("muon_pfisolationr3_sumneutralhadronethighthreshold", &muon_pfisolationr3_sumneutralhadronethighthreshold);
    tree->SetBranchAddress("muon_pfisolationr3_sumphotonethighthreshold", &muon_pfisolationr3_sumphotonethighthreshold);
    tree->SetBranchAddress("muon_pfisolationr3_sumpupt", &muon_pfisolationr3_sumpupt);
    tree->SetBranchAddress("muon_pfisolationr4_sumneutralhadronethighthreshold", &muon_pfisolationr4_sumneutralhadronethighthreshold);
    tree->SetBranchAddress("muon_pfisolationr4_sumphotonethighthreshold", &muon_pfisolationr4_sumphotonethighthreshold);
    tree->SetBranchAddress("muon_pfisolationr4_sumpupt", &muon_pfisolationr4_sumpupt);
    tree->SetBranchAddress("muon_pfisolationr4_sumchargedhadronpt", &muon_pfisolationr4_sumchargedhadronpt);
    tree->SetBranchAddress("muon_pfisolationr4_sumchargedparticlept", &muon_pfisolationr4_sumchargedparticlept);
    tree->SetBranchAddress("muon_pfisolationr4_sumneutralhadronet", &muon_pfisolationr4_sumneutralhadronet);
    tree->SetBranchAddress("muon_pfisolationr4_sumphotonet", &muon_pfisolationr4_sumphotonet);
    //tree->SetBranchAddress("muon_pfisolationr4_dBrel", &muon_pfisolationr4_dBrel);
    tree->SetBranchAddress("muon_ecalenergy", &muon_ecalenergy);
    tree->SetBranchAddress("muon_hcalenergy", &muon_hcalenergy);
    tree->SetBranchAddress("muon_charge", &muon_charge);
    tree->SetBranchAddress("muon_numchambers", &muon_numchambers);
    tree->SetBranchAddress("muon_numchamberswithsegments", &muon_numchamberswithsegments);
    tree->SetBranchAddress("muon_numvalidmuonhits", &muon_numvalidmuonhits);
    tree->SetBranchAddress("muon_nummatchedstations", &muon_nummatchedstations);
    tree->SetBranchAddress("muon_trigger", &muon_trigger);
    tree->SetBranchAddress("muon_trackermuonquality", &muon_trackermuonquality);
    tree->SetBranchAddress("muon_matches_IsoMu20", &muon_matches_IsoMu20);
    tree->SetBranchAddress("muon_matches_IsoTkMu20", &muon_matches_IsoTkMu20);

    tree->SetBranchAddress("ak4pfchsjet_count", &ak4pfchsjet_count);
    tree->SetBranchAddress("ak4pfchsjet_energy", &ak4pfchsjet_energy);
    tree->SetBranchAddress("ak4pfchsjet_px", &ak4pfchsjet_px);
    tree->SetBranchAddress("ak4pfchsjet_py", &ak4pfchsjet_py);
    tree->SetBranchAddress("ak4pfchsjet_pz", &ak4pfchsjet_pz);
    tree->SetBranchAddress("ak4pfchsjet_area", &ak4pfchsjet_area);
    tree->SetBranchAddress("ak4pfchsjet_hadronicenergy", &ak4pfchsjet_hadronicenergy);
    tree->SetBranchAddress("ak4pfchsjet_chargedhadronicenergy", &ak4pfchsjet_chargedhadronicenergy);
    tree->SetBranchAddress("ak4pfchsjet_emenergy", &ak4pfchsjet_emenergy);
    tree->SetBranchAddress("ak4pfchsjet_chargedemenergy", &ak4pfchsjet_chargedemenergy);
    tree->SetBranchAddress("ak4pfchsjet_hfemenergy", &ak4pfchsjet_hfemenergy);
    tree->SetBranchAddress("ak4pfchsjet_hfhadronicenergy", &ak4pfchsjet_hfhadronicenergy);
    tree->SetBranchAddress("ak4pfchsjet_electronenergy", &ak4pfchsjet_electronenergy);
    tree->SetBranchAddress("ak4pfchsjet_muonenergy", &ak4pfchsjet_muonenergy);
    tree->SetBranchAddress("ak4pfchsjet_chargedmulti", &ak4pfchsjet_chargedmulti);
    tree->SetBranchAddress("ak4pfchsjet_neutralmulti", &ak4pfchsjet_neutralmulti);
    tree->SetBranchAddress("ak4pfchsjet_hfhadronicmulti", &ak4pfchsjet_hfhadronicmulti);
    tree->SetBranchAddress("ak4pfchsjet_hfemmulti", &ak4pfchsjet_hfemmulti);
    tree->SetBranchAddress("ak4pfchsjet_electronmulti", &ak4pfchsjet_electronmulti);
    tree->SetBranchAddress("ak4pfchsjet_muonmulti", &ak4pfchsjet_muonmulti);
    tree->SetBranchAddress("ak4pfchsjet_chargeda", &ak4pfchsjet_chargeda);
    tree->SetBranchAddress("ak4pfchsjet_chargedb", &ak4pfchsjet_chargedb);
    tree->SetBranchAddress("ak4pfchsjet_neutrala", &ak4pfchsjet_neutrala);
    tree->SetBranchAddress("ak4pfchsjet_neutralb", &ak4pfchsjet_neutralb);
    tree->SetBranchAddress("ak4pfchsjet_alla", &ak4pfchsjet_alla);
    tree->SetBranchAddress("ak4pfchsjet_allb", &ak4pfchsjet_allb);
    tree->SetBranchAddress("ak4pfchsjet_chargedfractionmv", &ak4pfchsjet_chargedfractionmv);
    tree->SetBranchAddress("ak4pfchsjet_energycorr", &ak4pfchsjet_energycorr);
    tree->SetBranchAddress("ak4pfchsjet_energycorrunc", &ak4pfchsjet_energycorrunc);
    //tree->SetBranchAddress("ak4pfchsjet_energycorrl7uds", &ak4pfchsjet_energycorrl7uds);
    //tree->SetBranchAddress("ak4pfchsjet_energycorrl7bottom", &ak4pfchsjet_energycorrl7bottom);
    tree->SetBranchAddress("ak4pfchsjet_btag", &ak4pfchsjet_btag);
    tree->SetBranchAddress("ak4pfchsjet_mcflavour", &ak4pfchsjet_mcflavour);

    tree->SetBranchAddress("electron_count", &electron_count);
    //tree->SetBranchAddress("electron_vtx", &electron_vtx);
    tree->SetBranchAddress("electron_px", &electron_px);
    tree->SetBranchAddress("electron_py", &electron_py);
    tree->SetBranchAddress("electron_pz", &electron_pz);
    tree->SetBranchAddress("electron_correctedecalenergy", &electron_correctedecalenergy);
    tree->SetBranchAddress("electron_trackchi2", &electron_trackchi2);
    tree->SetBranchAddress("electron_trackndof", &electron_trackndof);
    //tree->SetBranchAddress("electron_outerx", &electron_outerx);
    //tree->SetBranchAddress("electron_outery", &electron_outery);
    //tree->SetBranchAddress("electron_outerz", &electron_outerz);
    //tree->SetBranchAddress("electron_closestpointx", &electron_closestpointx);
    //tree->SetBranchAddress("electron_closestpointy", &electron_closestpointy);
    //tree->SetBranchAddress("electron_closestpointz", &electron_closestpointz);
    tree->SetBranchAddress("electron_esuperclusterovertrack", &electron_esuperclusterovertrack);
    tree->SetBranchAddress("electron_eseedclusterovertrack", &electron_eseedclusterovertrack);
    tree->SetBranchAddress("electron_deltaetasuperclustertrack", &electron_deltaetasuperclustertrack);
    tree->SetBranchAddress("electron_deltaphisuperclustertrack", &electron_deltaphisuperclustertrack);
    tree->SetBranchAddress("electron_e1x5", &electron_e1x5);
    tree->SetBranchAddress("electron_e2x5", &electron_e2x5);
    tree->SetBranchAddress("electron_e5x5", &electron_e5x5);
    tree->SetBranchAddress("electron_r9", &electron_r9);
    tree->SetBranchAddress("electron_sigmaetaeta", &electron_sigmaetaeta);
    tree->SetBranchAddress("electron_sigmaietaieta", &electron_sigmaietaieta);
    tree->SetBranchAddress("electron_sigmaiphiiphi", &electron_sigmaiphiiphi);
    tree->SetBranchAddress("electron_ehcaloverecaldepth1", &electron_ehcaloverecaldepth1);
    tree->SetBranchAddress("electron_ehcaltoweroverecaldepth2", &electron_ehcaltoweroverecaldepth2);
    tree->SetBranchAddress("electron_ehcaltoweroverecaldepth1", &electron_ehcaltoweroverecaldepth1);
    tree->SetBranchAddress("electron_ehcaloverecaldepth2", &electron_ehcaloverecaldepth2);
    tree->SetBranchAddress("electron_isolationr3track", &electron_isolationr3track);
    tree->SetBranchAddress("electron_isolationr3ecal", &electron_isolationr3ecal);
    tree->SetBranchAddress("electron_isolationr3hcal", &electron_isolationr3hcal);
    //tree->SetBranchAddress("electron_pfisolationr3_dBrel", &electron_pfisolationr3_dBrel);
    tree->SetBranchAddress("electron_isolationr4track", &electron_isolationr4track);
    tree->SetBranchAddress("electron_isolationr4ecal", &electron_isolationr4ecal);
    tree->SetBranchAddress("electron_isolationr4hcal", &electron_isolationr4hcal);
    tree->SetBranchAddress("electron_isolationpfr3charged", &electron_isolationpfr3charged);
    tree->SetBranchAddress("electron_isolationpfr3photon", &electron_isolationpfr3photon);
    tree->SetBranchAddress("electron_isolationpfr3neutral", &electron_isolationpfr3neutral);
    tree->SetBranchAddress("electron_nhits", &electron_nhits);
    tree->SetBranchAddress("electron_npixelhits", &electron_npixelhits);
    tree->SetBranchAddress("electron_nmissinghits", &electron_nmissinghits);
    tree->SetBranchAddress("electron_npixellayers", &electron_npixellayers);
    tree->SetBranchAddress("electron_nstriplayers", &electron_nstriplayers);
    tree->SetBranchAddress("electron_nhitsexpected", &electron_nhitsexpected);
    tree->SetBranchAddress("electron_dxy", &electron_dxy);
    tree->SetBranchAddress("electron_dxyerr", &electron_dxyerr);
    tree->SetBranchAddress("electron_dz", &electron_dz);
    tree->SetBranchAddress("electron_dzerr", &electron_dzerr);
    //tree->SetBranchAddress("electron_convdist", &electron_convdist);
    //tree->SetBranchAddress("electron_convdcot", &electron_convdcot);
    //tree->SetBranchAddress("electron_convradius", &electron_convradius);
    tree->SetBranchAddress("electron_gapinfo", &electron_gapinfo);
    tree->SetBranchAddress("electron_fbrems", &electron_fbrems);
    tree->SetBranchAddress("electron_numbrems", &electron_numbrems);
    tree->SetBranchAddress("electron_charge", &electron_charge);

    tree->SetBranchAddress("electron_effectiveArea", &electron_effectiveArea);

    //tree->SetBranchAddress("electron_info", &electron_info);
    tree->SetBranchAddress("electron_cutBasedLoose", &electron_cutBasedLoose);
    tree->SetBranchAddress("electron_cutBasedMedium", &electron_cutBasedMedium);
    tree->SetBranchAddress("electron_cutBasedTight", &electron_cutBasedTight);
    tree->SetBranchAddress("electron_mvaNonTrigWP90", &electron_mvaNonTrigWP90);
    tree->SetBranchAddress("electron_mvaNonTrigWP80", &electron_mvaNonTrigWP80);
    tree->SetBranchAddress("electron_iselectron", &electron_iselectron);
    tree->SetBranchAddress("electron_passconversionveto", &electron_passconversionveto);
    tree->SetBranchAddress("electron_ecaldrivenseed", &electron_ecaldrivenseed);
    tree->SetBranchAddress("electron_trackerdrivenseed", &electron_trackerdrivenseed);
    tree->SetBranchAddress("electron_trigger", &electron_trigger);
    tree->SetBranchAddress("electron_supercluster_e", &electron_supercluster_e);
    tree->SetBranchAddress("electron_supercluster_x", &electron_supercluster_x);
    tree->SetBranchAddress("electron_supercluster_y", &electron_supercluster_y);
    tree->SetBranchAddress("electron_supercluster_z", &electron_supercluster_z);
    tree->SetBranchAddress("electron_supercluster_rawe", &electron_supercluster_rawe);
    tree->SetBranchAddress("electron_supercluster_phiwidth", &electron_supercluster_phiwidth);
    tree->SetBranchAddress("electron_supercluster_etawidth", &electron_supercluster_etawidth);
    tree->SetBranchAddress("electron_supercluster_nbasiccluster", &electron_supercluster_nbasiccluster);

    tree->SetBranchAddress("photon_count", &photon_count);
    tree->SetBranchAddress("photon_px", &photon_px);
    tree->SetBranchAddress("photon_py", &photon_py);
    tree->SetBranchAddress("photon_pz", &photon_pz);
    tree->SetBranchAddress("photon_e1x5", &photon_e1x5);
    tree->SetBranchAddress("photon_e2x5", &photon_e2x5);
    tree->SetBranchAddress("photon_e3x3", &photon_e3x3);
    tree->SetBranchAddress("photon_e5x5", &photon_e5x5);
    tree->SetBranchAddress("photon_sigmaietaieta", &photon_sigmaietaieta);
    tree->SetBranchAddress("photon_sigmaietaiphi", &photon_sigmaietaiphi);
    tree->SetBranchAddress("photon_sigmaiphiiphi", &photon_sigmaiphiiphi);
    tree->SetBranchAddress("photon_ehcaloverecaldepth1", &photon_ehcaloverecaldepth1);
    tree->SetBranchAddress("photon_ehcaloverecaldepth2", &photon_ehcaloverecaldepth2);
    tree->SetBranchAddress("photon_ehcaltoweroverecaldepth1", &photon_ehcaltoweroverecaldepth1);
    tree->SetBranchAddress("photon_ehcaltoweroverecaldepth2", &photon_ehcaltoweroverecaldepth2);
    tree->SetBranchAddress("photon_maxenergyxtal", &photon_maxenergyxtal);
    tree->SetBranchAddress("photon_isolationr3track", &photon_isolationr3track);
    tree->SetBranchAddress("photon_isolationr3trackhollow", &photon_isolationr3trackhollow);
    tree->SetBranchAddress("photon_isolationr3ntrack", &photon_isolationr3ntrack);
    tree->SetBranchAddress("photon_isolationr3ntrackhollow", &photon_isolationr3ntrackhollow);
    tree->SetBranchAddress("photon_isolationr3ecal", &photon_isolationr3ecal);
    tree->SetBranchAddress("photon_isolationr3hcal", &photon_isolationr3hcal);
    tree->SetBranchAddress("photon_isolationr4track", &photon_isolationr4track);
    tree->SetBranchAddress("photon_isolationr4trackhollow", &photon_isolationr4trackhollow);
    tree->SetBranchAddress("photon_isolationr4ntrack", &photon_isolationr4ntrack);
    tree->SetBranchAddress("photon_isolationr4ntrackhollow", &photon_isolationr4ntrackhollow);
    tree->SetBranchAddress("photon_isolationr4ecal", &photon_isolationr4ecal);
    tree->SetBranchAddress("photon_isolationr4hcal", &photon_isolationr4hcal);
    tree->SetBranchAddress("photon_isolationpfr3charged", &photon_isolationpfr3charged);
    tree->SetBranchAddress("photon_isolationpfr3photon", &photon_isolationpfr3photon);
    tree->SetBranchAddress("photon_isolationpfr3neutral", &photon_isolationpfr3neutral);
    //tree->SetBranchAddress("photon_isolationpfr4charged", &photon_isolationpfr4charged);
    //tree->SetBranchAddress("photon_isolationpfr4photon", &photon_isolationpfr4photon);
    //tree->SetBranchAddress("photon_isolationpfr4neutral", &photon_isolationpfr4neutral);
    //tree->SetBranchAddress("photon_isolationpfr4noscfootprintcharged", &photon_isolationpfr4noscfootprintcharged);
    //tree->SetBranchAddress("photon_isolationpfr4noscfootprintphoton", &photon_isolationpfr4noscfootprintphoton);
    //tree->SetBranchAddress("photon_isolationpfr4noscfootprintneutral", &photon_isolationpfr4noscfootprintneutral);
    tree->SetBranchAddress("photon_supercluster_e", &photon_supercluster_e);
    tree->SetBranchAddress("photon_supercluster_x", &photon_supercluster_x);
    tree->SetBranchAddress("photon_supercluster_y", &photon_supercluster_y);
    tree->SetBranchAddress("photon_supercluster_z", &photon_supercluster_z);
    tree->SetBranchAddress("photon_supercluster_rawe", &photon_supercluster_rawe);
    tree->SetBranchAddress("photon_supercluster_phiwidth", &photon_supercluster_phiwidth);
    tree->SetBranchAddress("photon_supercluster_etawidth", &photon_supercluster_etawidth);
    tree->SetBranchAddress("photon_supercluster_nbasiccluster", &photon_supercluster_nbasiccluster);
    //tree->SetBranchAddress("photon_info", &photon_info);
    tree->SetBranchAddress("photon_isphoton", &photon_isphoton);
    tree->SetBranchAddress("photon_hasconversiontracks", &photon_hasconversiontracks);
    tree->SetBranchAddress("photon_haspixelseed", &photon_haspixelseed);
    tree->SetBranchAddress("photon_passelectronveto", &photon_passelectronveto);
    tree->SetBranchAddress("photon_ispfphoton", &photon_ispfphoton);
    tree->SetBranchAddress("photon_gapinfo", &photon_gapinfo);
    //tree->SetBranchAddress("photon_conversionbegin", &photon_conversionbegin);

    tree->SetBranchAddress("tau_count", &tau_count);
    //tree->SetBranchAddress("tau_charged_count", &tau_charged_count);
    tree->SetBranchAddress("tau_px", &tau_px);
    tree->SetBranchAddress("tau_py", &tau_py);
    tree->SetBranchAddress("tau_pz", &tau_pz);
    tree->SetBranchAddress("tau_isolationneutralspt", &tau_isolationneutralspt);
    tree->SetBranchAddress("tau_isolationneutralsnum", &tau_isolationneutralsnum);
    tree->SetBranchAddress("tau_isolationchargedpt", &tau_isolationchargedpt);
    tree->SetBranchAddress("tau_isolationchargednum", &tau_isolationchargednum);
    tree->SetBranchAddress("tau_isolationgammapt", &tau_isolationgammapt);
    tree->SetBranchAddress("tau_isolationgammanum", &tau_isolationgammanum);
    tree->SetBranchAddress("tau_charge", &tau_charge);
    tree->SetBranchAddress("tau_disc", &tau_disc);
    //tree->SetBranchAddress("tau_emfraction", &tau_emfraction);
    //tree->SetBranchAddress("tau_hcaltotoverplead", &tau_hcaltotoverplead);
    //tree->SetBranchAddress("tau_hcal3x3overplead", &tau_hcal3x3overplead);
    //tree->SetBranchAddress("tau_ecalstripsumeoverplead", &tau_ecalstripsumeoverplead);
    //tree->SetBranchAddress("tau_bremsrecoveryeoverplead", &tau_bremsrecoveryeoverplead);
    //tree->SetBranchAddress("tau_calocomp", &tau_calocomp);
    //tree->SetBranchAddress("tau_segcomp", &tau_segcomp);
    tree->SetBranchAddress("tau_trigger", &tau_trigger);
    //tree->SetBranchAddress("tau_ak4pfjet_e", &tau_ak4pfjet_e);
    //tree->SetBranchAddress("tau_ak4pfjet_px", &tau_ak4pfjet_px);
    //tree->SetBranchAddress("tau_ak4pfjet_py", &tau_ak4pfjet_py);
    //tree->SetBranchAddress("tau_ak4pfjet_pz", &tau_ak4pfjet_pz);
    //tree->SetBranchAddress("tau_ak4pfjet_hadronicenergy", &tau_ak4pfjet_hadronicenergy);
    //tree->SetBranchAddress("tau_ak4pfjet_chargedhadronicenergy", &tau_ak4pfjet_chargedhadronicenergy);
    //tree->SetBranchAddress("tau_ak4pfjet_emenergy", &tau_ak4pfjet_emenergy);
    //tree->SetBranchAddress("tau_ak4pfjet_chargedemenergy", &tau_ak4pfjet_chargedemenergy);
    //tree->SetBranchAddress("tau_ak4pfjet_chargedmulti", &tau_ak4pfjet_chargedmulti);
    //tree->SetBranchAddress("tau_ak4pfjet_neutralmulti", &tau_ak4pfjet_neutralmulti);
/*
    tree->SetBranchAddress("tau_chargedbegin", &tau_chargedbegin);
    tree->SetBranchAddress("tau_charged_px", &tau_charged_px);
    tree->SetBranchAddress("tau_charged_py", &tau_charged_py);
    tree->SetBranchAddress("tau_charged_pz", &tau_charged_pz);
    tree->SetBranchAddress("tau_charged_outerx", &tau_charged_outerx);
    tree->SetBranchAddress("tau_charged_outery", &tau_charged_outery);
    tree->SetBranchAddress("tau_charged_outerz", &tau_charged_outerz);
    tree->SetBranchAddress("tau_charged_closestpointx", &tau_charged_closestpointx);
    tree->SetBranchAddress("tau_charged_closestpointy", &tau_charged_closestpointy);
    tree->SetBranchAddress("tau_charged_closestpointz", &tau_charged_closestpointz);
    tree->SetBranchAddress("tau_charged_chi2", &tau_charged_chi2);
    tree->SetBranchAddress("tau_charged_ndof", &tau_charged_ndof);
    tree->SetBranchAddress("tau_charged_dxy", &tau_charged_dxy);
    tree->SetBranchAddress("tau_charged_dxyerr", &tau_charged_dxyerr);
    tree->SetBranchAddress("tau_charged_dz", &tau_charged_dz);
    tree->SetBranchAddress("tau_charged_dzerr", &tau_charged_dzerr);
    tree->SetBranchAddress("tau_charged_dedxharmonic2", &tau_charged_dedxharmonic2);
    tree->SetBranchAddress("tau_charged_charge", &tau_charged_charge);
    tree->SetBranchAddress("tau_charged_nhits", &tau_charged_nhits);
    tree->SetBranchAddress("tau_charged_nmissinghits", &tau_charged_nmissinghits);
    tree->SetBranchAddress("tau_charged_npixelhits", &tau_charged_npixelhits);
    tree->SetBranchAddress("tau_charged_npixellayers", &tau_charged_npixellayers);
    tree->SetBranchAddress("tau_charged_nstriplayers", &tau_charged_nstriplayers);
*/

    tree->SetBranchAddress("event_rho", &event_rho);
    //tree->SetBranchAddress("ak4pfjet_sigma", &ak4pfjet_sigma);

    //tree->SetBranchAddress("pfmet_ex", &pfmet_ex);
    //tree->SetBranchAddress("pfmet_ey", &pfmet_ey);
    tree->SetBranchAddress("pfmettype1_ex", &pfmettype1_ex);
    tree->SetBranchAddress("pfmettype1_ey", &pfmettype1_ey);
    tree->SetBranchAddress("pfmettype1_et", &pfmettype1_et);
    //tree->SetBranchAddress("pfmetpuppitype1_ex", &pfmetpuppitype1_ex);
    //tree->SetBranchAddress("pfmetpuppitype1_ey", &pfmetpuppitype1_ey);
    //tree->SetBranchAddress("pfmettype0type1_ex", &pfmettype0type1_ex);
    //tree->SetBranchAddress("pfmettype0type1_ey", &pfmettype0type1_ey);

    tree->SetBranchAddress("genweight", &genweight);
    tree->SetBranchAddress("genx1",     &genx1);
    tree->SetBranchAddress("genid1",    &genid1);
    tree->SetBranchAddress("genx2",     &genx2);
    tree->SetBranchAddress("genid2",    &genid2);
    tree->SetBranchAddress("genScale",  &genScale);

    tree->SetBranchAddress("numpileupinteractionsminus", &numpileupinteractionsminus);
    tree->SetBranchAddress("numpileupinteractions",      &numpileupinteractions);
    tree->SetBranchAddress("numpileupinteractionsplus",  &numpileupinteractionsplus);
    tree->SetBranchAddress("numtruepileupinteractions",  &numtruepileupinteractions);

    tree->SetBranchAddress("genparticles_count", &genparticles_count);
    tree->SetBranchAddress("genparticles_e",      &genparticles_e);
    tree->SetBranchAddress("genparticles_px",     &genparticles_px);
    tree->SetBranchAddress("genparticles_py",     &genparticles_py);
    tree->SetBranchAddress("genparticles_pz",     &genparticles_pz);
    tree->SetBranchAddress("genparticles_vx",     &genparticles_vx);
    tree->SetBranchAddress("genparticles_vy",     &genparticles_vy);
    tree->SetBranchAddress("genparticles_vz",     &genparticles_vz);
    tree->SetBranchAddress("genparticles_pdgid",  &genparticles_pdgid);
    tree->SetBranchAddress("genparticles_status", &genparticles_status);
    tree->SetBranchAddress("genparticles_indirectmother", &genparticles_indirectmother);
    tree->SetBranchAddress("genparticles_info",        &genparticles_info);

    //tree->SetBranchAddress("genparticlesmother_count",    &genparticlesmother_count);
    //tree->SetBranchAddress("genparticlesdaughter_count",  &genparticlesdaughter_count);
    //tree->SetBranchAddress("genparticles_motherbeg",   &genparticles_motherbeg);
    //tree->SetBranchAddress("genparticles_daughterbeg", &genparticles_daughterbeg);
    //tree->SetBranchAddress("genparticles_mothers",     &genparticles_mothers);
    //tree->SetBranchAddress("genparticles_daughters",   &genparticles_daughters);

    tree->SetBranchAddress("genallparticles_count",         &genallparticles_count);
    tree->SetBranchAddress("genallparticlesmother_count",   &genallparticlesmother_count);
    tree->SetBranchAddress("genallparticlesdaughter_count", &genallparticlesdaughter_count);
    tree->SetBranchAddress("genallparticles_e",           &genallparticles_e);
    tree->SetBranchAddress("genallparticles_px",          &genallparticles_px);
    tree->SetBranchAddress("genallparticles_py",          &genallparticles_py);
    tree->SetBranchAddress("genallparticles_pz",          &genallparticles_pz);
    tree->SetBranchAddress("genallparticles_vx",          &genallparticles_vx);
    tree->SetBranchAddress("genallparticles_vy",          &genallparticles_vy);
    tree->SetBranchAddress("genallparticles_vz",          &genallparticles_vz);
    tree->SetBranchAddress("genallparticles_pdgid",       &genallparticles_pdgid);
    tree->SetBranchAddress("genallparticles_status",      &genallparticles_status);
    tree->SetBranchAddress("genallparticles_motherbeg",   &genallparticles_motherbeg);
    tree->SetBranchAddress("genallparticles_daughterbeg", &genallparticles_daughterbeg);
    tree->SetBranchAddress("genallparticles_mothers",     &genallparticles_mothers);
    tree->SetBranchAddress("genallparticles_daughters",   &genallparticles_daughters);

    tree->SetBranchAddress("genmet_ex", &genmet_ex);
    tree->SetBranchAddress("genmet_ey", &genmet_ey);

    tree->SetBranchAddress("genak4jet_count", &genak4jet_count);
    tree->SetBranchAddress("genak4jet_e",          &genak4jet_e);
    tree->SetBranchAddress("genak4jet_px",         &genak4jet_px);
    tree->SetBranchAddress("genak4jet_py",         &genak4jet_py);
    tree->SetBranchAddress("genak4jet_pz",         &genak4jet_pz);
    tree->SetBranchAddress("genak4jet_einvisible", &genak4jet_einvisible);
    tree->SetBranchAddress("genak4jet_flavour",    &genak4jet_flavour);
    tree->SetBranchAddress("genak4jet_info",       &genak4jet_info);

}

//____________________________________________________________________________________
void Analyse::UsePileUpInfo() {
    usepileupinfo = true;
}

//____________________________________________________________________________________
void Analyse::UsePrimVertexInfo() {
    useprimvertexinfo = true;
}

//____________________________________________________________________________________
void Analyse::LoadBeamSpot(bool select) {
    if(select) {
        loadbeamspot = 1;
    } else {
        loadbeamspot = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadMuons(bool select) {
    if(select) {
        loadmuons = 1;
    } else {
        loadmuons = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadElectrons(bool select) {
    if(select) {
        loadelectrons = 1;
    } else {
        loadelectrons = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadPhotons(bool select) {
    if(select) {
        loadphotons = 1;
    } else {
        loadphotons = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadTaus(bool select) {
    if(select) {
        loadtaus = 1;
    } else {
        loadtaus = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadMET(bool select) {
    if(select) {
        loadmet = 1;
    } else {
        loadmet = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadAK4PFCHSJets(bool select) {
    if(select) loadak4pfchsjets = 1;
    else loadak4pfchsjets = 0;
}

//____________________________________________________________________________________
void Analyse::LoadPrimVertices(bool select) {
    if(select) {
        loadprimvertices = 1;
    } else {
        loadprimvertices = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadTrigger(bool select) {
    if(select) {
        loadtrigger = 1;
    } else {
        loadtrigger = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadGenInfo(bool select) {
    if(select) {
        loadgeninfo = 1;
    } else {
        loadgeninfo = 0;
    }
}

/*
//____________________________________________________________________________________
void Analyse::LoadAllGenParticles(bool select) {
    if(select) {
        loadallgenparticles = 1;
    } else {
        loadallgenparticles = 0;
    }
}
*/
//____________________________________________________________________________________
void Analyse::LoadGenParticles(bool select) {
    if(select) {
        loadgenparticles = 1;
    } else {
        loadgenparticles = 0;
    }
}

//____________________________________________________________________________________
void Analyse::LoadGenJets(bool select) {
    if(select) {
        loadgenak4jets = 1;
    } else {
        loadgenak4jets = 0;
    }
}
/*
//____________________________________________________________________________________
bool Analyse::GetL1Trigger(UInt_t bit) const {
    if((trigger_level1[bit/8] & 1<<(bit % 8)) > 0) {
        return(true);
    } else {
        return(false);
    }
}

//____________________________________________________________________________________
bool Analyse::GetL1TriggerBits(UInt_t bit) const {
    if((trigger_level1bits[bit/8] & 1<<(bit % 8)) > 0) {
        return(true);
    } else {
        return(false);
    }
}

//____________________________________________________________________________________
bool Analyse::GetHLTrigger(UInt_t index) const {
    if((trigger_HLT[index/8] & 1<<(index % 8)) > 0) {
        return(true);
    } else {
        return(false);
    }
}
*/
//____________________________________________________________________________________
BeamSpot Analyse::GetBeamSpot() const {
    return(BeamSpot(beamspot_x, beamspot_y, beamspot_z, beamspot_cov, beamspot_xwidth, beamspot_ywidth, beamspot_zsigma));
}

//____________________________________________________________________________________
Muon Analyse::Muons(UInt_t n, int correction) const {
    //vector<string> a = runlist.find(Run())->second.GetHLTMuonNames();
    Muon newmuon(this, n, correction);
    Track Outertrack(sqrt(pow(muon_outertrack_px->at(n),2) + pow(muon_outertrack_py->at(n),2) + pow(muon_outertrack_pz->at(n),2)), muon_outertrack_px->at(n), muon_outertrack_py->at(n), muon_outertrack_pz->at(n), 10000., 10000., 10000., 1000., 1000., 1000., muon_outertrack_chi2->at(n), muon_outertrack_ndof->at(n), 0., 0., 0., 0., muon_charge->at(n), muon_outertrack_hits->at(n), muon_outertrack_missinghits->at(n), 0, 0, 0, -1, 0.);
    newmuon.SetOuterTrack(Outertrack);

   // removed pending the fix of muon_innertrack_vtx
    //Track Innertrack(sqrt(muon_innertrack_px->at(n)*muon_innertrack_px->at(n)+muon_innertrack_py->at(n)*muon_innertrack_py->at(n)+muon_innertrack_pz->at(n)*muon_innertrack_pz->at(n)), muon_innertrack_px->at(n), muon_innertrack_py->at(n), muon_innertrack_pz->at(n), muon_innertrack_outerx->at(n), muon_innertrack_outery->at(n), muon_innertrack_outerz->at(n), muon_innertrack_closestpointx->at(n), muon_innertrack_closestpointy->at(n), muon_innertrack_closestpointz->at(n), muon_innertrack_chi2->at(n), muon_innertrack_ndof->at(n), muon_innertrack_dxy->at(n), muon_innertrack_dxyerr->at(n), muon_innertrack_dz->at(n), muon_innertrack_dzerr->at(n), muon_innertrack_charge->at(n), muon_innertrack_nhits->at(n), muon_innertrack_nmissinghits->at(n), muon_innertrack_npixelhits->at(n), muon_innertrack_npixellayers->at(n), muon_innertrack_nstriplayers->at(n), muon_innertrack_vtx->at(n), muon_innertrack_dedxharmonic2->at(n));
    //newmuon.SetInnerTrack(Innertrack);

    return(newmuon);
}

//____________________________________________________________________________________
Electron Analyse::Electrons(UInt_t n, int correction) const {
    return(Electron(this, n, correction));
}

//____________________________________________________________________________________
Photon Analyse::Photons(UInt_t n) const {
    Photon newphoton(this, n);
//
//    UInt_t begin = photon_conversionbegin->at(n);
//    UInt_t end = conversion_count;
//    if(n < photon_count - 1) {
//        end = photon_conversionbegin->at(n+1);
//    }
//
//    for(UInt_t i = begin ; i < end ; i++) {
//        Conversion newconversion(conversion_vx[i], conversion_vy[i], conversion_vz[i], conversion_chi2[i], conversion_ndof[i], conversion_cov[i], conversion_mvaout[i], TVector3(conversion_trackecalpointx[i][0], conversion_trackecalpointy[i][0], conversion_trackecalpointz[i][0]), TVector3(conversion_trackecalpointx[i][1], conversion_trackecalpointy[i][1], conversion_trackecalpointz[i][1]));
//        if(conversion_trackndof[i][0] != -1) {
//            Double_t E = sqrt(conversion_trackpx[i][0]*conversion_trackpx[i][0]+conversion_trackpy[i][0]*conversion_trackpy[i][0]+conversion_trackpz[i][0]*conversion_trackpz[i][0]);
//            Track newtrack(E, conversion_trackpx[i][0], conversion_trackpy[i][0], conversion_trackpz[i][0], conversion_trackecalpointx[i][0], conversion_trackecalpointy[i][0], conversion_trackecalpointz[i][0], conversion_trackclosestpointx[i][0], conversion_trackclosestpointy[i][0], conversion_trackclosestpointz[i][0], conversion_trackchi2[i][0], conversion_trackndof[i][0], conversion_trackdxy[i][0], conversion_trackdxyerr[i][0], conversion_trackdz[i][0], conversion_trackdzerr[i][0], conversion_trackcharge[i][0], conversion_tracknhits[i][0], conversion_tracknmissinghits[i][0], conversion_tracknpixelhits[i][0], conversion_tracknpixellayers[i][0], conversion_tracknstriplayers[i][0], -1, 0.);
//            newconversion.AddTrack(newtrack);
//        }
//        if(conversion_trackndof[i][1] != -1) {
//            Double_t E = sqrt(conversion_trackpx[i][1]*conversion_trackpx[i][1]+conversion_trackpy[i][1]*conversion_trackpy[i][1]+conversion_trackpz[i][1]*conversion_trackpz[i][1]);
//            Track newtrack(E, conversion_trackpx[i][1], conversion_trackpy[i][1], conversion_trackpz[i][1], conversion_trackecalpointx[i][1], conversion_trackecalpointy[i][1], conversion_trackecalpointz[i][1], conversion_trackclosestpointx[i][1], conversion_trackclosestpointy[i][1], conversion_trackclosestpointz[i][1], conversion_trackchi2[i][1], conversion_trackndof[i][1], conversion_trackdxy[i][1], conversion_trackdxyerr[i][1], conversion_trackdz[i][1], conversion_trackdzerr[i][1], conversion_trackcharge[i][1], conversion_tracknhits[i][1], conversion_tracknmissinghits[i][1], conversion_tracknpixelhits[i][1], conversion_tracknpixellayers[i][1], conversion_tracknstriplayers[i][1], -1, 0.);
//            newconversion.AddTrack(newtrack);
//        }
//        if(newconversion.Mag() > 0.01) {
//            newphoton.AddConversion(newconversion);
//        }
//
//    }
    return(newphoton);
}

//____________________________________________________________________________________
Tau Analyse::Taus(UInt_t n) const {
    return(Tau(this, n));
}


//TLorentzVector Analyse::PFMET() const {
//    return(TLorentzVector(pfmet_ex, pfmet_ey, 0., sqrt(pow(pfmet_ex, 2) + pow(pfmet_ey, 2))));
//}

//____________________________________________________________________________________
TLorentzVector Analyse::PFMETTYPE1() const {
    // because of how RootMaker works, met is always a vector<Float_t> with only one entry
    return(TLorentzVector(pfmettype1_ex->at(0), pfmettype1_ey->at(0), 0., sqrt(pow(pfmettype1_ex->at(0), 2) + pow(pfmettype1_ey->at(0), 2))));
}

//____________________________________________________________________________________
Jet Analyse::AK4PFCHSJets(UInt_t n) const {
    //return(Jet(ak4pfchsjet_energy->at(n), ak4pfchsjet_px->at(n), ak4pfchsjet_py->at(n), ak4pfchsjet_pz->at(n), ak4pfchsjet_hadronicenergy->at(n), ak4pfchsjet_chargedhadronicenergy->at(n), ak4pfchsjet_emenergy->at(n), ak4pfchsjet_chargedemenergy->at(n), ak4pfchsjet_hfemenergy->at(n), ak4pfchsjet_hfhadronicenergy->at(n), ak4pfchsjet_electronenergy->at(n), ak4pfchsjet_muonenergy->at(n), ak4pfchsjet_chargedmulti->at(n), ak4pfchsjet_neutralmulti->at(n), ak4pfchsjet_hfhadronicmulti->at(n), ak4pfchsjet_hfemmulti->at(n), ak4pfchsjet_electronmulti->at(n), ak4pfchsjet_muonmulti->at(n), ak4pfchsjet_chargeda->at(n), ak4pfchsjet_chargedb->at(n), ak4pfchsjet_neutrala->at(n), ak4pfchsjet_neutralb->at(n), ak4pfchsjet_alla->at(n), ak4pfchsjet_allb->at(n), ak4pfchsjet_chargedfractionmv->at(n), ak4pfchsjet_energycorr->at(n), ak4pfchsjet_energycorrunc->at(n), ak4pfchsjet_energycorrl7uds->at(n), ak4pfchsjet_energycorrl7bottom->at(n), ak4pfchsjet_btag->at(n), ak4pfchsjet_mcflavour->at(n), 0., 0., 0.));
    return(Jet(ak4pfchsjet_energy->at(n), ak4pfchsjet_px->at(n), ak4pfchsjet_py->at(n), ak4pfchsjet_pz->at(n), ak4pfchsjet_hadronicenergy->at(n), ak4pfchsjet_chargedhadronicenergy->at(n), ak4pfchsjet_emenergy->at(n), ak4pfchsjet_chargedemenergy->at(n), ak4pfchsjet_hfemenergy->at(n), ak4pfchsjet_hfhadronicenergy->at(n), ak4pfchsjet_electronenergy->at(n), ak4pfchsjet_muonenergy->at(n), ak4pfchsjet_chargedmulti->at(n), ak4pfchsjet_neutralmulti->at(n), ak4pfchsjet_hfhadronicmulti->at(n), ak4pfchsjet_hfemmulti->at(n), ak4pfchsjet_electronmulti->at(n), ak4pfchsjet_muonmulti->at(n), ak4pfchsjet_chargeda->at(n), ak4pfchsjet_chargedb->at(n), ak4pfchsjet_neutrala->at(n), ak4pfchsjet_neutralb->at(n), ak4pfchsjet_alla->at(n), ak4pfchsjet_allb->at(n), ak4pfchsjet_chargedfractionmv->at(n), ak4pfchsjet_energycorr->at(n), ak4pfchsjet_energycorrunc->at(n), ak4pfchsjet_btag->at(n), ak4pfchsjet_mcflavour->at(n), 0., 0., 0.));
}

//____________________________________________________________________________________
Vertex Analyse::PrimVertices(UInt_t n) const {
    return(Vertex(primvertex_x->at(n), primvertex_y->at(n), primvertex_z->at(n), primvertex_isfake->at(n), primvertex_isvalid->at(n), primvertex_chi2->at(n), primvertex_ndof->at(n), primvertex_ntracks->at(n), primvertex_cov0->at(n), primvertex_cov1->at(n), primvertex_cov2->at(n), primvertex_cov3->at(n), primvertex_cov4->at(n), primvertex_cov5->at(n)));
}

//____________________________________________________________________________________
GenLightParticle Analyse::GenParticles(UInt_t n) const {
    return(GenLightParticle(genparticles_e->at(n), genparticles_px->at(n), genparticles_py->at(n), genparticles_pz->at(n), genparticles_vx->at(n), genparticles_vy->at(n), genparticles_vz->at(n), genparticles_status->at(n), genparticles_pdgid->at(n), genparticles_info->at(n), genparticles_indirectmother->at(n)));
    //return(GenLightParticle(genparticles_e->at(n), genparticles_px->at(n), genparticles_py->at(n), genparticles_pz->at(n), genparticles_vx->at(n), genparticles_vy->at(n), genparticles_vz->at(n), genparticles_status->at(n), genparticles_pdgid->at(n)));
}

//____________________________________________________________________________________
GenJet Analyse::GenJets(UInt_t n) const {
    return(GenJet(genak4jet_e->at(n), genak4jet_px->at(n), genak4jet_py->at(n), genak4jet_pz->at(n), genak4jet_einvisible->at(n)));
}

//____________________________________________________________________________________
GenParticle Analyse::AllGenParticles(UInt_t n) const {
    UInt_t motherend(n == genallparticles_count-1 ? genallparticlesmother_count : genallparticles_motherbeg->at(n+1));
    UInt_t daughterend(n == genallparticles_count-1 ? genallparticlesdaughter_count : genallparticles_daughterbeg->at(n+1));

    GenParticle newpart(this, n, genallparticles_e->at(n), genallparticles_px->at(n), genallparticles_py->at(n), genallparticles_pz->at(n), genallparticles_vx->at(n), genallparticles_vy->at(n), genallparticles_vz->at(n), genallparticles_status->at(n), genallparticles_pdgid->at(n), genallparticles_motherbeg->at(n), motherend - genallparticles_motherbeg->at(n), genallparticles_daughterbeg->at(n), daughterend - genallparticles_daughterbeg->at(n));

    return(newpart);
}

////____________________________________________________________________________________
//TLorentzVector Analyse::GenMET() const {
//    return(TLorentzVector(genmet_ex, genmet_ey, 0., sqrt(pow(genmet_ex, 2) + pow(genmet_ey, 2))));
//}

//____________________________________________________________________________________
bool Analyse::IsBatchSelected(UInt_t run, UInt_t lumiblock) {

    if(Batch_MyId() <= 0) return(true);
    if(batchselection.find(run) != batchselection.end() && batchselection.find(run)->second.find(lumiblock) != batchselection.find(run)->second.end()) {
        return(true);
    }
    return(false);
}

//____________________________________________________________________________________
void Analyse::Batch_Prepare(bool simple) {
#ifdef USE_MPI
    cout << "START Batch_Prepare " << Batch_MyId() << endl;
    if(Batch_MyId() == 0) {
        if(simple) {
            map<string, vector<pair<UInt_t, UInt_t> > > filelumimap;
            UInt_t numlumisections = 0;
            for(map< UInt_t, map< UInt_t, Luminosity > >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a) {
                for(map< UInt_t, Luminosity >::iterator b = a->second.begin() ; b != a->second.end() ; ++b) {
                    if((jsonfilter && jsonlist[a->first][b->first]) || !jsonfilter) {
                        numlumisections++;
                        const vector<string> &filesofls = b->second.GetFiles();
                        for(int i = 0 ; i < filesofls.size() ; i++) {
                            filelumimap[filesofls[i]].push_back(pair<UInt_t, UInt_t>(a->first, b->first));
                        }
                    }
                }
            }
            cerr << "JOB " << Batch_MyId() << ": " << "Total number of lumisections: " << numlumisections << endl;

            UInt_t numfiles = filelumimap.size();
            vector<UInt_t> lumiout;
            vector<string> files;
            string filesout;
            UInt_t jobnum = 1;
            UInt_t stepval = numfiles/(Batch_NumJobs()-1);
            UInt_t stepmod = numfiles % (Batch_NumJobs()-1);
            UInt_t i = 0;
            for(map<string, vector<pair<UInt_t, UInt_t> > >::iterator a = filelumimap.begin() ; a != filelumimap.end() ; ++a) {
                i++;
                files.push_back(a->first);
                for(vector<pair<UInt_t, UInt_t> >::iterator b = a->second.begin() ; b != a->second.end() ; ++b) {
                    if((jsonfilter && jsonlist[b->first][b->second]) || !jsonfilter) {
                        lumiout.push_back(b->first);
                        lumiout.push_back(b->second);
                    }
                }
                if(i == jobnum*stepval + Min(jobnum, stepmod)) {
                    filesout += files[0];
                    for(size_t j = 1 ; j < files.size() ; j++) {
                        filesout += " " + files[j];
                    }
                    int size = lumiout.size();
                    MPI_Send(&size, 1, MPI_INTEGER, jobnum, 1, MPI_COMM_WORLD);
                    MPI_Send(&(lumiout.front()), size, MPI_UNSIGNED, jobnum, 2, MPI_COMM_WORLD);
                    lumiout.clear();
                    size = filesout.size()+1;
                    cout << "Job " << i << " " << jobnum  << " should run on files(simple): " << filesout << endl;
                    MPI_Send(&size, 1, MPI_INTEGER, jobnum, 3, MPI_COMM_WORLD);
                    MPI_Send(const_cast<char *>(filesout.c_str()), size, MPI_CHAR, jobnum, 4, MPI_COMM_WORLD);
                    files.clear();
                    filesout.clear();
                    jobnum++;
                }
            }
            filesout = filelumimap.begin()->first;
            lumiout.clear();
            lumiout.push_back(0);
            lumiout.push_back(0);
            for(int i = jobnum ; i < Batch_NumJobs() ; i++) {
                cout << "Job " << i << ": Warning: running on zero events. " << lumiout.size() << " " << filesout << endl;
                int size = lumiout.size();
                MPI_Send(&size, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD);
                MPI_Send(&(lumiout.front()), size, MPI_UNSIGNED, i, 2, MPI_COMM_WORLD);
                size = filesout.size()+1;
                MPI_Send(&size, 1, MPI_INTEGER, i, 3, MPI_COMM_WORLD);
                MPI_Send(const_cast<char *>(filesout.c_str()), size, MPI_CHAR, i, 4, MPI_COMM_WORLD);
            }

        } else {
            //map<string, vector<pair<UInt_t, UInt_t> >, bool(*)(string, string)> filelumimap(sortfiles);
            map<string, vector<pair<UInt_t, UInt_t> > > filelumimap;
            UInt_t numlumisections = 0;
            for(map< UInt_t, map< UInt_t, Luminosity > >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a) {
                for(map< UInt_t, Luminosity >::iterator b = a->second.begin() ; b != a->second.end() ; ++b) {
                    if((jsonfilter && jsonlist[a->first][b->first]) || !jsonfilter) {
                        numlumisections++;
                        const vector<string> &filesofls = b->second.GetFiles();
                        filelumimap[filesofls[0]].push_back(pair<UInt_t, UInt_t>(a->first, b->first));
                    }
                }
            }
            cerr << "JOB " << Batch_MyId() << ": " << "Total number of lumisections: " << numlumisections << endl;

            vector<UInt_t> lumiout;
            vector<string> files;
            string filesout;
            UInt_t jobnum = 1;
            UInt_t stepval = numlumisections/(Batch_NumJobs()-1);
            UInt_t stepmod = numlumisections % (Batch_NumJobs()-1);
            UInt_t i = 0;
            for(map<string, vector<pair<UInt_t, UInt_t> > >::iterator a = filelumimap.begin() ; a != filelumimap.end() ; ++a) {
                for(vector<pair<UInt_t, UInt_t> >::iterator b = a->second.begin() ; b != a->second.end() ; ++b) {
                    if((jsonfilter && jsonlist[b->first][b->second]) || !jsonfilter) {
                        i++;
                        lumiout.push_back(b->first);
                        lumiout.push_back(b->second);
                        const vector<string> &filesofls = lumilist[b->first][b->second].GetFiles();
                        for(size_t f = 0 ; f < filesofls.size() ; f++) {
                            if(find(files.begin(), files.end(), filesofls[f]) == files.end()) {
                                files.push_back(filesofls[f]);
                            }
                        }
                        if(i == jobnum*stepval + Min(jobnum, stepmod)) {
                            //sort(files.begin(), files.end(), sortfiles);
                            filesout = files[0];
                            for(size_t j = 1 ; j < files.size() ; j++) {
                                //if(files[j] != files[j-1])
                                //{
                                filesout += " " + files[j];
                                //}
                            }
                            int size = lumiout.size();
                            MPI_Send(&size, 1, MPI_INTEGER, jobnum, 1, MPI_COMM_WORLD);
                            MPI_Send(&(lumiout.front()), size, MPI_UNSIGNED, jobnum, 2, MPI_COMM_WORLD);
                            lumiout.clear();
                            size = filesout.size()+1;
                            //cerr << "Job " << i << " " << jobnum  << " should run on files: " << filesout << endl;
                            MPI_Send(&size, 1, MPI_INTEGER, jobnum, 3, MPI_COMM_WORLD);
                            MPI_Send(const_cast<char *>(filesout.c_str()), size, MPI_CHAR, jobnum, 4, MPI_COMM_WORLD);
                            files.clear();
                            filesout.clear();
                            jobnum++;
                        }
                    }
                }
            }
            filesout = filelumimap.begin()->first;
            lumiout.clear();
            lumiout.push_back(0);
            lumiout.push_back(0);
            for(int i = jobnum ; i < Batch_NumJobs() ; i++) {
                cout << "Job " << i << ": Warning: running on zero events. " << lumiout.size() << " " << filesout << endl;
                int size = lumiout.size();
                MPI_Send(&size, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD);
                MPI_Send(&(lumiout.front()), size, MPI_UNSIGNED, i, 2, MPI_COMM_WORLD);
                size = filesout.size()+1;
                MPI_Send(&size, 1, MPI_INTEGER, i, 3, MPI_COMM_WORLD);
                MPI_Send(const_cast<char *>(filesout.c_str()), size, MPI_CHAR, i, 4, MPI_COMM_WORLD);
            }

        }
    } else if(Batch_MyId() > 0) {
        int size;
        MPI_Recv(&size, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &mpistatus);
        UInt_t *intbuff = new UInt_t[size];
        MPI_Recv(intbuff, size, MPI_UNSIGNED, 0, 2, MPI_COMM_WORLD, &mpistatus);
        if(intbuff[0] == 0 && intbuff[1] == 0) {
            cout << "JOB " << Batch_MyId() << ": Warning: there is nothing to do :-( ." << endl;
            batch_emptyjob = true;
        }
        for(int i = 0 ; i < size/2 ; i++) {
            batchselection[intbuff[2*i]][intbuff[2*i+1]] = true;
            //cerr << Batch_MyId() << ": " <<intbuff[2*i] << " " << intbuff[2*i+1] << endl;
        }
        delete[] intbuff;

        MPI_Recv(&size, 1, MPI_INTEGER, 0, 3, MPI_COMM_WORLD, &mpistatus);
        char *buff = new char[size];
        MPI_Recv(buff, size, MPI_CHAR, 0, 4, MPI_COMM_WORLD, &mpistatus);
        string filename;
        vector<string> lumifiles;
        vector<string> inputfiles;
        boost::split(inputfiles, buff, boost::is_any_of(" "));
        for(size_t i = 0 ; i < inputfiles.size() ; i++) {
            AddFile(inputfiles[i]);
            string lumifile = inputfiles[i].substr(0, inputfiles[i].find_last_of("/")) + "/LUMI_INFO.root";
            if(find(lumifiles.begin(), lumifiles.end(), lumifile) == lumifiles.end()) {
                lumifiles.push_back(lumifile);
            }
        }
cerr<<"aww yeah"<<endl;
        for(size_t i = 0 ; i < lumifiles.size() ; i++) {
cerr<<"you gotta"<<endl;
            AddLumiFile(lumifiles[i]);
cerr<<"schwifty"<<endl;
        }
cerr<<"gotta get schwifty in here"<<endl;
        delete[] buff;
    }
#endif
}

//____________________________________________________________________________________
Long64_t Analyse::Loop(Long64_t start, Long64_t end) {
    if(end == -1 || end > GetNumAddedEvents()) {
        end = GetNumAddedEvents();
    }
    if(Batch_MyId() > 0) {
        start = 0;
        end = GetNumAddedEvents();
    }

    if(usepileupinfo && useprimvertexinfo) {
        cerr << "Don't use UsePileUpInfo() and UsePrimVertexInfo() at the same time." << endl;
    }
    GetEvent(0, 1);
    if(IsMC()) {
#ifdef USE_MPI
        if(Batch_MyId() == 0) {
            if(usepileupinfo) {
                for(int i = 1 ; i < Batch_NumJobs() ; i++) {
                    vector<double> buff(pileUpDist.size());
                    cerr << "Master: waiting for Pileup." << i << endl;
                    MPI_Recv(&(buff.front()), buff.size(), MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &mpistatus);
                    for(size_t j = 0 ; j < pileUpDist.size() ; j++) {
                        pileUpDist[j] += buff[j];
                    }
                }
                for(int i = 1 ; i < Batch_NumJobs() ; i++) {
                    MPI_Send(&(pileUpDist.front()), pileUpDist.size(), MPI_DOUBLE, i, 6, MPI_COMM_WORLD);
                    cerr << "Master: sent Pileup to." << i << endl;
                }
            } else if(useprimvertexinfo) {
                for(int i = 1 ; i < Batch_NumJobs() ; i++) {
                    vector<double> buff(primVertexDist.size());
                    cerr << "Master: waiting for Pileup." << i << endl;
                    MPI_Recv(&(buff.front()), buff.size(), MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &mpistatus);
                    for(size_t j = 0 ; j < primVertexDist.size() ; j++) {
                        primVertexDist[j] += buff[j];
                    }
                }
                for(int i = 1 ; i < Batch_NumJobs() ; i++) {
                    MPI_Send(&(primVertexDist.front()), primVertexDist.size(), MPI_DOUBLE, i, 6, MPI_COMM_WORLD);
                    cerr << "Master: sent Pileup to " << i << endl;
                }
            }
        }
#endif
        if(Batch_MyId() != 0) {
            //Pileup Reweighting
            if(usepileupinfo) {
                cerr << "Reading PileUp-Information" << endl;
                if(!batch_emptyjob) {
                    for(Long64_t i = start ; i < end ; i++) {
                        GetEvent(i, 1);
                        if(IsBatchSelected(Run(), LumiBlock())) {
                            pileUpDist[unsigned(pileupscale*2*NumTruePileUpInteractions())]++;
                        }
                    }
                }
#ifdef USE_MPI
                if(Batch_MyId() > 0) {
                    MPI_Send(&(pileUpDist.front()), pileUpDist.size(), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
                    MPI_Recv(&(pileUpDist.front()), pileUpDist.size(), MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &mpistatus);
                }
#endif
                Double_t totnum = accumulate(pileUpDist.begin(), pileUpDist.end(), 0.);
                for(unsigned i = 0 ; i < pileUpDist.size() ; i++) {
                    pileUpDist[i] /= totnum;
                }
            }

            //Vertex Reweighting
            if(useprimvertexinfo) {
                if(!batch_emptyjob) {
                    cerr << "Reading PrimVertex-Information" << endl;
                    for(Long64_t i = start ; i < end ; i++) {
                        GetEvent(i, 2);
                        if(IsBatchSelected(Run(), LumiBlock())) {
                            primVertexDist[unsigned(NumGoodPrimVertices())]++;
                        }
                    }
                }
#ifdef USE_MPI
                if(Batch_MyId() > 0) {
                    MPI_Send(&(primVertexDist.front()), primVertexDist.size(), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
                    cerr << "JOB: " << Batch_MyId() << " PU sent" << endl;
                    MPI_Recv(&(primVertexDist.front()), primVertexDist.size(), MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &mpistatus);
                    cerr << "JOB: " << Batch_MyId() << " PU recv" << endl;
                }
#endif
                Double_t totnum = accumulate(primVertexDist.begin(), primVertexDist.end(), 0.);
                for(unsigned i = 0 ; i < primVertexDist.size() ; i++) {
                    primVertexDist[i] /= totnum;
                }
            }
        }
    } else {
        useprimvertexinfo = false;
        usepileupinfo = false;
    }

    processed = 0;
#ifdef USE_MPI
    if(Batch_MyId() == 0) {
        for(int i = 1 ; i < Batch_NumJobs() ; i++) {
            UInt_t donebyjob = 0;
            int size;
            MPI_Recv(&size, 1, MPI_INTEGER, i, 7, MPI_COMM_WORLD, &mpistatus);
            UInt_t *intbuff = new UInt_t[size];
            MPI_Recv(intbuff, size, MPI_UNSIGNED, i, 8, MPI_COMM_WORLD, &mpistatus);

            for(int j = 0 ; j < size/3 ; j++) {
                UInt_t run = intbuff[3*j];
                UInt_t lumi = intbuff[3*j+1];
                UInt_t numevents = intbuff[3*j+2];
                donebyjob+=numevents;
                processed+=numevents;
                if(run == 0 && lumi == 0) break; //empty job.
                lumilist[run][lumi].EventsProcessed(numevents);
                if(run < minRun) {
                    minRun = run;
                    minLumi = lumi;
                }
                if(run == minRun && lumi < minLumi) {
                    minLumi = lumi;
                }
                if(run > maxRun) {
                    maxRun = run;
                    maxLumi = lumi;
                }
                if(run == maxRun && lumi > maxLumi) {
                    maxLumi = lumi;
                }
            }
            delete[] intbuff;
            cerr << "Master: Job " << i << " has processed " << donebyjob << " events in " << size/3 << " lumisections." << endl;
        }
        cerr << "Master: Total events: " << processed << endl;
        return(processed);
    }
#endif
    BeginLoop();
    if(!batch_emptyjob) {
        if(Batch_MyId() != -1) cerr << "JOB " << Batch_MyId() << ": ";
        if(IsData()) cerr<<"Sample will be treated as DATA."<<endl;
        else cerr<<"Sample will be treated as MC."<<endl;
        cerr << GetNumAddedEvents() << " Events are added." << endl;
        // Looping
        if(Batch_MyId() != -1) cerr << "JOB " << Batch_MyId() << ": ";
        cerr << "Events " << start << " - " << end << " will be processed (" << end - start << " Events)." << endl;
        TStopwatch watch;
        Int_t evres = 0;
        double timeremain;
        watch.Start();
        for(Long64_t i = start ; i < end ; i++) {
            GetEvent(i);

            if(IsBatchSelected(Run(), LumiBlock())) {
                if(jsonfilter) {
                    if(jsonlist[Run()][LumiBlock()]) {
                        evres = AnalyseEvent();
                    } else {
                        evres = 0;
                    }
                } else {
                    evres = AnalyseEvent();
                }

                if(duplicatecheck) {
                    eventlist[Run()][LumiBlock()][Number()]++;
                }

                if(lumicalculation && evres > 0) {
                    ++lumilist[Run()][LumiBlock()];
                    //Setting the Range
                    if(Run() < minRun) {
                        minRun = Run();
                        minLumi = LumiBlock();
                    }
                    if(Run() == minRun && LumiBlock() < minLumi) {
                        minLumi = LumiBlock();
                    }
                    if(Run() > maxRun) {
                        maxRun = Run();
                        maxLumi = LumiBlock();
                    }
                    if(Run() == maxRun && LumiBlock() > maxLumi) {
                        maxLumi = LumiBlock();
                    }
                }
            }
            processed++;

            if(printinfo > 0 && processed%printinfo == 0) {
                //watch.Stop();
                timeremain = watch.RealTime()/processed*(end-start-processed);
                watch.Continue();
                if(Batch_MyId() != -1) cerr << "JOB " << Batch_MyId() << ": ";
                cerr<<"Processed: "<<processed<<"/"<<end-start<<", Event: "<<UInt_t(Number())<<", Run: "<<Run()<<", LumiBlock: "<<LumiBlock()<<", Time: "<<int(timeremain)/3600<<":"<<int(timeremain)%3600/60<<":"<<int(timeremain)%60<<", "<<(100*processed/(end-start))<<"%"<<endl;
            }

            if(evres == -1) {
                break;
            }
        }
    }

#ifdef USE_MPI
    if(Batch_MyId() > 0) {
        vector<UInt_t> lumiinfo;

        for(map< UInt_t, map< UInt_t, bool > >::iterator a = batchselection.begin() ; a != batchselection.end() ; ++a) {
            for(map< UInt_t, bool >::iterator b = a->second.begin() ; b != a->second.end() ; ++b) {
                UInt_t run = a->first;
                UInt_t lumi = b->first;
                UInt_t numprocessed = lumilist[run][lumi].NumEventsProcessed();
                //UInt_t numevents = lumilist[run][lumi].NumEvents();
                //if(numevents != numprocessed)
                //{
                //cerr << "JOB " << Batch_MyId() << ": ";
                //cerr << "Missing Events:" << run << " " << lumi << " " << numprocessed << "/" << numevents << endl;
                //}
                lumiinfo.push_back(run);
                lumiinfo.push_back(lumi);
                lumiinfo.push_back(numprocessed);
            }
        }
        int size = lumiinfo.size();
        cout << "JOB " << Batch_MyId() << ": sending result." << endl;
        MPI_Send(&size, 1, MPI_INTEGER, 0, 7, MPI_COMM_WORLD);
        MPI_Send(&(lumiinfo.front()), size, MPI_UNSIGNED, 0, 8, MPI_COMM_WORLD);
    }
#endif
    EndLoop();
    return(processed);
}

//____________________________________________________________________________________
Int_t Analyse::PrepareSkimming(string filename) {
//    skimfilename = filename;
//    if(skimfile != 0) {
        cerr << "this function shouldn't exist" << endl;
        return(-1);
//    }
}

//____________________________________________________________________________________
Int_t Analyse::SkimEvent() {
    if(skimtree == 0) {
        cerr << "use PrepareSkimming() before SkimEvent()" << endl;
        return(-1);
    }
    selected[Run()][LumiBlock()]++;
    skimtree->Fill();
    return(1);
}

//____________________________________________________________________________________
Int_t Analyse::SetLumi(UInt_t run, UInt_t block, Float_t lumival, Float_t avgpu) {
    if(lumilist.find(run) != lumilist.end() && lumilist.find(run)->second.find(block) != lumilist.find(run)->second.end()) {
        lumilist[run][block].Value(lumival);
        lumilist[run][block].AvgPU(avgpu);
        return(1);
    }
    return(0);
}

//____________________________________________________________________________________
void Analyse::WriteLumiFile(string filename) {
    TFile *file = new TFile(filename.c_str(), "recreate");

    UInt_t run_number;
    UInt_t run_hltcount;
    Char_t run_hltnames[20000];
    Char_t run_hltmunames[10000];
    Char_t run_hltelnames[10000];
    Char_t run_hlttaunames[10000];
    Char_t run_hltphotonnames[10000];
    Char_t run_hltjetnames[10000];
    Char_t run_taudiscriminators[10000];
    UInt_t run_hltprescaletablescount;
    UInt_t run_hltprescaletables[10000];
    UInt_t run_l1algocount;
    UInt_t run_l1algoprescaletablescount;
    UInt_t run_l1algoprescaletables[10000];
    UInt_t run_l1techcount;
    UInt_t run_l1techprescaletablescount;
    UInt_t run_l1techprescaletables[10000];

    TTree *runtree = new TTree("AC1Brun", "AC1Brun", 1);
    runtree->Branch("run_number", &run_number, "run_number/i");
    runtree->Branch("run_hltcount", &run_hltcount, "run_hltcount/i");
    runtree->Branch("run_hltnames", run_hltnames, "run_hltnames[run_hltcount]/C");
    runtree->Branch("run_hltmunames", run_hltmunames, "run_hltmunames/C");
    runtree->Branch("run_hltelnames", run_hltelnames, "run_hltelnames/C");
    runtree->Branch("run_hlttaunames", run_hlttaunames, "run_hlttaunames/C");
    runtree->Branch("run_hltphotonnames", run_hltphotonnames, "run_hltphotonnames/C");
    runtree->Branch("run_hltjetnames", run_hltjetnames, "run_hltjetnames/C");
    runtree->Branch("run_taudiscriminators", run_taudiscriminators, "run_taudiscriminators/C");
    runtree->Branch("run_hltprescaletablescount", &run_hltprescaletablescount, "run_hltprescaletablescount/i");
    runtree->Branch("run_hltprescaletables", run_hltprescaletables, "run_hltprescaletables[run_hltprescaletablescount]/i");
    runtree->Branch("run_l1algoprescaletablescount", &run_l1algoprescaletablescount, "run_l1algoprescaletablescount/i");
    runtree->Branch("run_l1algocount", &run_l1algocount, "run_l1algocount/i");
    runtree->Branch("run_l1algoprescaletables", run_l1algoprescaletables, "run_l1algoprescaletables[run_l1algoprescaletablescount]/i");
    runtree->Branch("run_l1techprescaletablescount", &run_l1techprescaletablescount, "run_l1techprescaletablescount/i");
    runtree->Branch("run_l1techcount", &run_l1techcount, "run_l1techcount/i");
    runtree->Branch("run_l1techprescaletables", run_l1techprescaletables, "run_l1techprescaletables[run_l1techprescaletablescount]/i");

    for(map<UInt_t, RunInfo>::iterator b = runlist.begin() ; b != runlist.end() ; ++b) {
        run_number = b->second.Run();
        run_hltcount = b->second.NumHLT();
        strcpy(run_hltnames, b->second.HLTAllNames().c_str());
        strcpy(run_hltmunames, b->second.HLTMuonAllNames().c_str());
        strcpy(run_hltelnames, b->second.HLTElectronAllNames().c_str());
        strcpy(run_hlttaunames, b->second.HLTTauAllNames().c_str());
        strcpy(run_hltphotonnames, b->second.HLTPhotonAllNames().c_str());
        strcpy(run_hltjetnames, b->second.HLTJetAllNames().c_str());
        strcpy(run_taudiscriminators, b->second.TauDiscriminatorsAllNames().c_str());
        run_hltprescaletablescount = b->second.NumHLT()*b->second.NumHLTTables();
        for(UInt_t i = 0 ; i < b->second.NumHLTTables() ; i++) {
            for(UInt_t j = 0 ; j < b->second.NumHLT() ; j++) {
                run_hltprescaletables[j+b->second.NumHLT()*i] = b->second.HLTPrescale(j,i);
            }
        }

        run_l1algocount = b->second.NumL1Algo();
        run_l1algoprescaletablescount = b->second.NumL1Algo()*b->second.NumL1AlgoTables();
        for(UInt_t i = 0 ; i < b->second.NumL1AlgoTables() ; i++) {
            for(UInt_t j = 0 ; j < b->second.NumL1Algo() ; j++) {
                run_l1algoprescaletables[j+b->second.NumL1Algo()*i] = b->second.L1AlgoPrescale(j,i);
            }
        }

        run_l1techcount = b->second.NumL1Tech();
        run_l1techprescaletablescount = b->second.NumL1Tech()*b->second.NumL1TechTables();
        for(UInt_t i = 0 ; i < b->second.NumL1TechTables() ; i++) {
            for(UInt_t j = 0 ; j < b->second.NumL1Tech() ; j++) {
                run_l1techprescaletables[j+b->second.NumL1Tech()*i] = b->second.L1TechPrescale(j,i);
            }
        }
        runtree->Fill();
    }
    runtree->Write();

    UInt_t lumi_run;
    UInt_t lumi_block;
    Float_t lumi_value;
    Float_t lumi_valueerr;
    Float_t lumi_livefrac;
    Float_t lumi_deadfrac;
    Float_t lumi_avgpu;
    UInt_t lumi_quality;
    UInt_t lumi_eventsprocessed;
    UInt_t lumi_eventsfiltered;
    UInt_t lumi_hltprescaletable;
    UInt_t lumi_l1algoprescaletable;
    UInt_t lumi_l1techprescaletable;
    Char_t lumi_filenames[100000];

    TTree *lumitree = new TTree("AC1Blumi", "AC1Blumi", 1);
    lumitree->Branch("lumi_run", &lumi_run, "lumi_run/i");
    lumitree->Branch("lumi_block", &lumi_block, "lumi_block/i");
    lumitree->Branch("lumi_value", &lumi_value, "lumi_value/F");
    lumitree->Branch("lumi_valueerr", &lumi_valueerr, "lumi_valueerr/F");
    lumitree->Branch("lumi_livefrac", &lumi_livefrac, "lumi_livefrac/F");
    lumitree->Branch("lumi_deadfrac", &lumi_deadfrac, "lumi_deadfrac/F");
    lumitree->Branch("lumi_avgpu", &lumi_avgpu, "lumi_avgpu/F");
    lumitree->Branch("lumi_quality", &lumi_quality, "lumi_quality/i");
    lumitree->Branch("lumi_eventsprocessed", &lumi_eventsprocessed, "lumi_eventsprocessed/i");
    lumitree->Branch("lumi_eventsfiltered", &lumi_eventsfiltered, "lumi_eventsfiltered/i");
    lumitree->Branch("lumi_hltprescaletable", &lumi_hltprescaletable, "lumi_hltprescaletable/i");
    lumitree->Branch("lumi_l1algoprescaletable", &lumi_l1algoprescaletable, "lumi_l1algoprescaletable/i");
    lumitree->Branch("lumi_l1techprescaletable", &lumi_l1techprescaletable, "lumi_l1techprescaletable/i");
    lumitree->Branch("lumi_filenames", lumi_filenames, "lumi_filenames/C");

    for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a) {
        for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; ++b) {
            lumi_run = a->first;
            lumi_block = b->first;
            lumi_value = b->second.LumiValue();
            lumi_valueerr = b->second.LumiError();
            lumi_livefrac = b->second.LiveFraction();
            lumi_deadfrac = b->second.DeadFraction();
            lumi_avgpu = b->second.AvgPU();
            lumi_quality = b->second.Quality();
            lumi_eventsprocessed = b->second.NumEventsOrig();
            if(skimtree != 0) {
                lumi_eventsfiltered = selected[a->first][b->first];
            } else {
                lumi_eventsfiltered = b->second.NumEvents();
            }
            lumi_hltprescaletable = b->second.HLTTable();
            lumi_l1algoprescaletable = b->second.L1AlgoTable();
            lumi_l1techprescaletable = b->second.L1TechTable();
            strcpy(lumi_filenames, b->second.GetFilesString().c_str());
            if(jsonfilter && jsonlist[lumi_run][lumi_block]) lumitree->Fill();
            if(!jsonfilter) lumitree->Fill();
        }

    }
    lumitree->Write();

    file->Write();
    file->Close();
}

//____________________________________________________________________________________
void Analyse::AddLumiFile(string filename, string dir) {

    TFile *lumifile = new TFile(filename.c_str());
    if(lumifile->IsZombie()) {
        cerr << "ERROR AddLumiFile: " << filename << " is not a valid file." << endl;
    }
    TTree *lumitree = 0;
    lumifile->GetObject((dir + string("AC1Blumi")).c_str(), lumitree);
    if(lumitree == 0) {
        cerr << "ERROR AddLumiFile: " << filename << " " << dir + string("AC1Blumi") << " is no TTree." << endl;
    }
    TTree *runtree = 0;
    lumifile->GetObject((dir + string("AC1Brun")).c_str(), runtree);
    if(runtree == 0) {
        cerr << "ERROR AddLumiFile: " << filename << " " << dir + string("AC1Brun") << " is no TTree." << endl;
    }

    UInt_t lumi_run;
    UInt_t lumi_block;
    Float_t lumi_value;
    Float_t lumi_valueerr;
    Float_t lumi_livefrac;
    Float_t lumi_deadfrac;
    Float_t lumi_avgpu = -1.;
    UInt_t lumi_quality;
    UInt_t lumi_eventsprocessed;
    UInt_t lumi_eventsfiltered;
    UInt_t lumi_hltprescaletable;
    UInt_t lumi_l1algoprescaletable;
    UInt_t lumi_l1techprescaletable;
    Char_t lumi_filenames[50000];
    vector<string> filenames;

    lumitree->SetBranchAddress("lumi_run", &lumi_run);
    lumitree->SetBranchAddress("lumi_block", &lumi_block);
    lumitree->SetBranchAddress("lumi_value", &lumi_value);
    lumitree->SetBranchAddress("lumi_valueerr", &lumi_valueerr);
    lumitree->SetBranchAddress("lumi_livefrac", &lumi_livefrac);
    lumitree->SetBranchAddress("lumi_deadfrac", &lumi_deadfrac);
    if(lumitree->GetBranch("lumi_avgpu") != 0) {
        lumitree->SetBranchAddress("lumi_avgpu", &lumi_avgpu);
    }
    lumitree->SetBranchAddress("lumi_quality", &lumi_quality);
    lumitree->SetBranchAddress("lumi_eventsprocessed", &lumi_eventsprocessed);
    lumitree->SetBranchAddress("lumi_eventsfiltered", &lumi_eventsfiltered);
    lumitree->SetBranchAddress("lumi_hltprescaletable", &lumi_hltprescaletable);
    lumitree->SetBranchAddress("lumi_l1algoprescaletable", &lumi_l1algoprescaletable);
    lumitree->SetBranchAddress("lumi_l1techprescaletable", &lumi_l1techprescaletable);
    if(lumitree->GetBranch("lumi_filenames") != 0) {
        lumitree->SetBranchAddress("lumi_filenames", lumi_filenames);
    }

    string filedir = filename.substr(0, filename.find_last_of("/")+1);
    for(UInt_t i = 0 ; i < lumitree->GetEntries() ; i++) {
        lumitree->GetEntry(i);
        if(!IsBatchSelected(lumi_run, lumi_block)) continue;

        filenames.clear();
        if(Batch_MyId() != -1 && lumitree->GetBranch("lumi_filenames") != 0) {
            boost::split(filenames, lumi_filenames, boost::is_any_of(" "));
            for(size_t i = 0 ; i < filenames.size() ; i++) {
                filenames[i] = filedir + filenames[i];
            }
        } else {
            filenames.push_back(filename);
        }

        if(lumilist.find(lumi_run) == lumilist.end() || lumilist.find(lumi_run)->second.find(lumi_block) == lumilist.find(lumi_run)->second.end()) {
            lumilist[lumi_run][lumi_block] = Luminosity(lumi_value, lumi_valueerr, lumi_livefrac, lumi_deadfrac, lumi_avgpu, lumi_quality, lumi_eventsfiltered, lumi_eventsprocessed, lumi_l1algoprescaletable, lumi_l1algoprescaletable, lumi_l1techprescaletable, filenames);
        } else {
            Luminosity newlumi(lumi_value, lumi_valueerr, lumi_livefrac, lumi_deadfrac, lumi_avgpu, lumi_quality, lumi_eventsfiltered, lumi_eventsprocessed, lumi_l1algoprescaletable, lumi_l1algoprescaletable, lumi_l1techprescaletable, filenames);
            lumilist[lumi_run][lumi_block] += newlumi;

            //cout << "Warning: Information of run: " << lumi_run << ", block: " << lumi_block << " has already been filled. For data this means there are duplicated files for MC it might be ok. Check for duplicated events anyway!" << endl;
        }
    }

    UInt_t run_number;
    UInt_t run_hltcount;
    Char_t run_hltnames[20000];
    Char_t run_hltmunames[10000];
    Char_t run_hltelnames[10000];
    Char_t run_hlttaunames[10000];
    Char_t run_hltphotonnames[10000];
    Char_t run_hltjetnames[10000];
    Char_t run_taudiscriminators[10000];
    UInt_t run_hltprescaletablescount;
    UInt_t run_hltprescaletables[10000];
    UInt_t run_l1algocount;
    UInt_t run_l1algoprescaletablescount;
    UInt_t run_l1algoprescaletables[10000];
    UInt_t run_l1techcount;
    UInt_t run_l1techprescaletablescount;
    UInt_t run_l1techprescaletables[10000];

    runtree->SetBranchAddress("run_number", &run_number);
    runtree->SetBranchAddress("run_hltcount", &run_hltcount);
    runtree->SetBranchAddress("run_hltnames", run_hltnames);
    runtree->SetBranchAddress("run_hltmunames", run_hltmunames);
    runtree->SetBranchAddress("run_hltelnames", run_hltelnames);
    runtree->SetBranchAddress("run_hlttaunames", run_hlttaunames);
    runtree->SetBranchAddress("run_hltphotonnames", run_hltphotonnames);
    runtree->SetBranchAddress("run_hltjetnames", run_hltjetnames);
    runtree->SetBranchAddress("run_taudiscriminators", run_taudiscriminators);
    runtree->SetBranchAddress("run_hltprescaletablescount", &run_hltprescaletablescount);
    runtree->SetBranchAddress("run_hltprescaletables", run_hltprescaletables);
    runtree->SetBranchAddress("run_l1algocount", &run_l1algocount);
    runtree->SetBranchAddress("run_l1algoprescaletablescount", &run_l1algoprescaletablescount);
    runtree->SetBranchAddress("run_l1algoprescaletables", run_l1algoprescaletables);
    runtree->SetBranchAddress("run_l1techcount", &run_l1techcount);
    runtree->SetBranchAddress("run_l1techprescaletablescount", &run_l1techprescaletablescount);
    runtree->SetBranchAddress("run_l1techprescaletables", run_l1techprescaletables);

    for(UInt_t i = 0 ; i < runtree->GetEntries() ; i++) {
        runtree->GetEntry(i);
        if(runlist.find(run_number) == runlist.end()) {
            runlist[run_number] = RunInfo(run_number, run_hltcount, run_hltnames, run_hltmunames, run_hltelnames, run_hlttaunames, run_hltphotonnames, run_hltjetnames, run_taudiscriminators, run_hltprescaletablescount, run_hltprescaletables, run_l1algocount, run_l1algoprescaletablescount, run_l1algoprescaletables, run_l1techcount, run_l1techprescaletablescount, run_l1techprescaletables);
        }
    }

    lumicalculation = true;
    lumifile->Close();
}

//____________________________________________________________________________________
Int_t Analyse::IsLumiAvailable() const {
    if(lumilist.find(Run()) != lumilist.end() && lumilist.find(Run())->second.find(LumiBlock()) != lumilist.find(Run())->second.end()) {
        if(lumilist.at(Run()).at(LumiBlock()).LumiValue() != -1) { //We have lumi value or it is MC.
            return(2);
        }
        return(1);
    }
    return(0);
}

//____________________________________________________________________________________
void Analyse::ResetLumiValues() {
    for(map< UInt_t, map< UInt_t, Luminosity > >::iterator a = lumilist.begin() ; a != lumilist.end() ; a++) {
        for(map<UInt_t, Luminosity>::iterator b = (a->second).begin() ; b != (a->second).end() ; b++) {
            b->second.Value(-1.);
        }
    }
}

/*
//____________________________________________________________________________________
Int_t Analyse::GetNumHLTriggers() const {
    if(runlist.find(Run()) != runlist.end()) {
        return(runlist.at(Run()).NumHLT());
    }
    return(-1);
}
//____________________________________________________________________________________
Int_t Analyse::GetHLTriggerIndex(string triggername) const {
    if(runlist.find(Run()) != runlist.end()) {
        return(runlist.at(Run()).HLTIndex(triggername));
    } else {
        cerr << "ERROR GetHLTriggerIndex: Trigger " << triggername << " not found." << endl;
        return(-1);
    }
}

//____________________________________________________________________________________
TriggerSelection *Analyse::AddTriggerSelection(string id, vector<string> triggernames, bool useprescaled) {
    TriggerSelection *triggerselection = new TriggerSelection(this, triggernames, useprescaled);
    triggerselections[id] = triggerselection;
    return(triggerselection);
}

//____________________________________________________________________________________
TriggerSelection *Analyse::GetTriggerSelection(string id) {
    return(triggerselections[id]);
}

//____________________________________________________________________________________
Int_t Analyse::GetHLTrigger(vector<string> triggernames) const {
    Int_t result = -1;
    for(UInt_t i = 0 ; i < triggernames.size() ; i++) {
        Int_t index = GetHLTriggerIndex(triggernames[i]);

        if(index != -1 && GetHLTPrescale(index) == 1) {
            result = 0;
            if(GetHLTrigger(index)) return(1);
        }
    }
    return(result);
}

//____________________________________________________________________________________
string Analyse::GetHLTriggerName(UInt_t index) const {
    if(runlist.find(Run()) != runlist.end()) {
        return(runlist.at(Run()).HLTName(index));
    } else {
        cerr << "ERROR GetHLTriggerName: Could not find name of trigger with index " << index << endl;
        return(string("NotFound"));
    }
}

//____________________________________________________________________________________
Int_t Analyse::GetHLTPrescale(UInt_t triggerindex) const {
    if(IsLumiAvailable()) {
        Int_t triggertable = lumilist.at(Run()).at(LumiBlock()).HLTTable();
        return(runlist.find(Run())->second.HLTPrescale(triggerindex, triggertable));
    } else {
        return(-1);
    }
}
*/

//____________________________________________________________________________________
bool Analyse::EventPassesHLT(std::vector<string> hltnames) const {
    bool matchedAny = false;
    //for (auto &path : hltnames) {// range-based 'for' loops are not allowed in C++98 mode :(
    for (size_t i = 0; i < hltnames.size(); i++) {
        if (hltnames[i] == "IsoMu20")        matchedAny = matchedAny || (bool(IsoMu20Pass));
        else if (hltnames[i] == "IsoTkMu20") matchedAny = matchedAny || (bool(IsoTkMu20Pass));

        else { cerr<<"HLT name "<<hltnames[i]<<" not recognized."<<endl; throw; }
    }

    return(matchedAny);
}

//____________________________________________________________________________________
void Analyse::PrintPrescaleInfo(string triggername) {
    for(map< UInt_t , RunInfo>::iterator a = runlist.begin() ; a != runlist.end() ; ++a) {
        Int_t index = a->second.HLTIndex(triggername);
        if(index < 0) {
            cout << a->first << ": " << triggername << " not available." << endl;
        }
        cout << a->first << ":" << index;
        for(UInt_t i = 0 ; i < a->second.NumHLTTables() ; i++) {
            cout << " " << a->second.HLTPrescale(index, i);
        }
        cout << lumilist.at(a->first).size() <<  endl;
        for(map<UInt_t, Luminosity>::iterator b = lumilist.at(a->first).begin() ; b != lumilist.at(a->first).end() ; b++) {
            Int_t triggertable = b->second.HLTTable();
            Int_t prescale = a->second.HLTPrescale(index, triggertable);
            if(prescale != 1) {
                cout << "(" << a->first << ":" << prescale << "),";
            }
        }
        cout << endl;
    }

}


//____________________________________________________________________________________
void Analyse::PrintPrescaleInfoB(string triggername) {
    triggername += string("_v.*");
    boost::cmatch what;
    boost::regex triggernameregex = boost::regex(triggername.c_str());
    Int_t precur = -2;
    Int_t triggertablecur = -1;
    UInt_t runcur = 0;
    UInt_t lumicur = 0;
    UInt_t runlast = 0;
    UInt_t lumilast = 0;
    Double_t lumi = 0.;
    string triggernamecur;
    Int_t triggerindexcur = -1;
    for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a) {
        string triggernamenow;
        RunInfo &runinfo = runlist.find(a->first)->second;
        for(UInt_t i = 0 ; i < runinfo.NumHLT() ; i++) {
            if(boost::regex_match(runinfo.HLTName(i).c_str(), what, triggernameregex)) {
                triggernamenow = runinfo.HLTName(i);
                break;
            }
        }
        Int_t triggerindex = runinfo.HLTIndex(triggernamenow);
        for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++) {
            if(b->second.LumiValue() == -1) continue;
            lumi+=b->second.LumiValue();
            Int_t triggertablenow = b->second.HLTTable();
            Int_t prenow = -1;
            if(triggerindex != -1) {
                prenow = runinfo.HLTPrescale(triggerindex, triggertablenow);
            }


            if(precur == -2) {
                triggernamecur = triggernamenow;
                triggertablecur = triggertablenow;
                triggerindexcur = triggerindex;
                precur = prenow;
                runcur = a->first;
                lumicur = b->first;
            }

            if(precur != prenow || triggertablecur != triggertablenow || triggernamecur != triggernamenow || triggerindexcur != triggerindex) {
                cout << triggerindexcur << " " << triggernamecur << ", " << runcur << ":" << lumicur << " - " << runlast << ":" << lumilast << ", prescale: " << precur << ", column: "<< triggertablecur  << ", lumi: " << lumi << "/pb"<< endl;
                lumi = 0.;
                triggerindexcur = triggerindex;
                triggernamecur = triggernamenow;
                triggertablecur = triggertablenow;
                precur = prenow;
                runcur = a->first;
                lumicur = b->first;
            }
            runlast = a->first;
            lumilast = b->first;
        }
    }
    cout << triggerindexcur << " " << triggernamecur << ", " << runcur << ":" << lumicur << " - " << runlast << ":" << lumilast << ", prescale: " << precur << ", column: " << triggertablecur  << ", lumi: " << lumi << "/pb"<< endl;
}

Double_t Analyse::GetInstLumi() const {
    map<UInt_t, map<UInt_t, Luminosity> >::const_iterator a = lumilist.find(Run());
    if(a == lumilist.end()) return(-1.);
    map<UInt_t, Luminosity>::const_iterator b = a->second.find(LumiBlock());
    if(b == a->second.end()) return(-1.);
    return(b->second.LumiValue());
}

Double_t Analyse::GetAvgPU() const {
    map<UInt_t, map<UInt_t, Luminosity> >::const_iterator a = lumilist.find(Run());
    if(a == lumilist.end()) return(-1.);
    map<UInt_t, Luminosity>::const_iterator b = a->second.find(LumiBlock());
    if(b == a->second.end()) return(-1.);
    return(b->second.AvgPU());
}

//____________________________________________________________________________________
Double_t Analyse::GetLumi(Int_t format) {
    Double_t lumi = 0.;
    Double_t alllumi = 0.;
    Double_t zerolumi = 0.;
    UInt_t nolumiinfo = 0;
    Double_t events = 0.;
    Double_t eventsprocessed = 0.;
    for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a) {
        Double_t runlumi = 0.;
        Double_t runalllumi = 0.;
        Double_t runzerolumi = 0.;
        Double_t runevents = 0.;
        Double_t runeventsprocessed = 0.;
        UInt_t runnolumiinfo = 0;
        UInt_t numblocks = 0;
        for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++) {
            if(b->second.LumiValue() == -1) {
                continue;
            }
            if(IsInRange(a->first, b->first)) {
                numblocks++;
                //	cout << a->first << " " << b->first << " " << minLumi << " " << maxLumi << endl;
                if(b->second == true) {
                    if(b->second.NumEvents() > 0) {
                        runlumi += b->second.ProcessedLumi();
                        runalllumi += b->second.LumiValue();
                        runevents += b->second.NumEvents();
                        runeventsprocessed += b->second.NumEventsProcessed();
                        if(b->second.ProcessedFraction() > 1) {
                            cerr << "WARNING GetLumi: You ran on " << b->second.ProcessedFraction()*100 <<"% of available events in Run " << a->first << ", Lumiblock " << b->first <<". :-O" << endl;
                        }
                        if(format >= 3) cout << "    Block: " << b->first << ", fraction: " << b->second.ProcessedFraction() << ", lumi: " << b->second.ProcessedLumi() << " pb^-1" << endl;
                    } else {
                        runzerolumi += b->second.LumiValue();
                        runalllumi += b->second.LumiValue();
                    }
                } else if(b->second == false) {
                    runnolumiinfo += b->second.NumEventsProcessed();
                }
            }
        }
        if(a->first <= maxRun && a->first >= minRun) {
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
        cout << "From run " << minRun << ", block " << minLumi << " to run " << maxRun << ", block " << maxLumi << "." << endl;
        cout << "Total lumi: " << lumi+zerolumi << " pb^-1, zero event lumi: " << zerolumi << " pb^-1, (lumi in Range: " << alllumi << " pb^-1), no lumi-info for " << nolumiinfo << " event(s), events: " << eventsprocessed << "(" << events << ")"  << endl;
    }
    return(lumi+zerolumi);
}

//____________________________________________________________________________________
bool Analyse::IsInRange(UInt_t theRun, UInt_t theLumiBlock) {
    if(theRun < maxRun && theRun > minRun) {
        return(true);
    } else if(maxRun != minRun && theRun == maxRun && theLumiBlock <= maxLumi) {
        return(true);
    } else if(maxRun != minRun && theRun == minRun && theLumiBlock >= minLumi) {
        return(true);
    } else if(maxRun == minRun && theRun == minRun && theLumiBlock <= maxLumi && theLumiBlock >= minLumi) {
        return(true);
    } else {
        return(false);
    }
}

//Double_t Analyse::GetLumiBlockLumi()
//{
//}

//____________________________________________________________________________________
void Analyse::PrintLumiOfRuns() {
    Double_t totlumi = 0.;
    Double_t filteredevents = 0;
    Double_t allevents = 0;
    for(map<UInt_t, map<UInt_t, Luminosity> >::iterator a = lumilist.begin() ; a != lumilist.end() ; ++a) {
        Double_t runlumi = 0;
        Double_t runfilteredevents = 0;
        Double_t runallevents = 0;
        for(map<UInt_t, Luminosity>::iterator b = a->second.begin() ; b != a->second.end() ; b++) {
            if(b->second.LumiValue() != -1.) {
                runlumi+=b->second.LumiValue();
                totlumi+=b->second.LumiValue();
            } else {
                cout << "Missing luminosity value for section: " << a->first << " " << b->first << endl;
            }
            filteredevents+=b->second.NumEvents();
            runfilteredevents+=b->second.NumEvents();
            allevents+=b->second.NumEventsOrig();
            runallevents+=b->second.NumEventsOrig();
        }
        cout << "Run: " << a->first << " -> " << runlumi << ", Events: " << (UInt_t)runfilteredevents << "/" << (UInt_t)runallevents <<"="<< runfilteredevents/runallevents << "  sum: " << totlumi << ", Events: " << (UInt_t)filteredevents << "/" << (UInt_t)allevents <<"="<< filteredevents/allevents<< endl;
    }
}

//____________________________________________________________________________________
void Analyse::PrintLumiOfLumiSectionsInRun(UInt_t runnumber) {
    Double_t runlumi = 0.;
    Double_t runfilteredevents = 0.;
    Double_t runallevents = 0.;
    map<UInt_t, map<UInt_t, Luminosity> >::iterator blocklist = lumilist.find(runnumber);
    for(map<UInt_t, Luminosity>::iterator b = blocklist->second.begin() ; b != blocklist->second.end() ; ++b) {
        if(b->second.LumiValue() != -1.) {
            runlumi+=b->second.LumiValue();
            runfilteredevents+=b->second.NumEvents();
            runallevents+=b->second.NumEventsOrig();
            cout << "Section: " << b->first << ", Lumi: " << b->second.LumiValue() << "("<<runlumi<<"), Events: "<< b->second.NumEvents() << "(" << runfilteredevents << ")" << endl;
        } else {
            cout << "Missing lumi info for section: " << b->first << endl;
        }
    }

}

//____________________________________________________________________________________
bool Analyse::LoadJSON(string filename) {
    jsonfilter = true;
    char buf[100];
    string strjson;
    fstream json(filename.c_str());
    while(true) {
        json.read(buf, 100);
        int read = json.gcount();
        strjson += string(buf, read);
        if(read < 100) break;
    }

    UInt_t run;
    UInt_t lumibeg;
    UInt_t lumiend;
    size_t posquota = 0;
    size_t posquotb = 0;
    size_t runbraon = 0;
    size_t runbraoff = 0;
    size_t lumibraon = 0;
    size_t lumibraoff = 0;
    size_t sep;
    bool result = true;
    while(true) {
        posquota = strjson.find("\"", runbraoff+1);
        posquotb = strjson.find("\"", posquota+1);
        if(posquota == string::npos) break;
        run = atoi(strjson.substr(posquota+1, posquotb-1).c_str());
        runbraon = strjson.find("[", posquotb+1);
        lumibraoff = runbraon;
        while(true) {
            lumibraon = strjson.find("[", lumibraoff+1);
            lumibraoff = strjson.find("]", lumibraoff+1);
            if(lumibraoff < lumibraon) break;
            sep = strjson.find(",", lumibraon+1);
            lumibeg = atoi(strjson.substr(lumibraon+1, sep-1).c_str());
            lumiend = atoi(strjson.substr(sep+1, lumibraoff-1).c_str());

            for(UInt_t u = lumibeg ; u <= lumiend ; u++) {
                if(lumilist.find(run) != lumilist.end() && lumilist.find(run)->second.find(u) != lumilist.find(u)->second.end()) {
                    jsonlist[run][u] = true;
                } else {
                    cout << "JSONFilter: Run: " << run << ", LS: " << u << " is not part of your input files." << endl;
                    result = false;
                }
            }

        }
        runbraoff = lumibraoff;
    }
    return(result);
}

//____________________________________________________________________________________
Double_t Analyse::GetPileUpMaxWeight(vector<Double_t> &datadist) const {
    double maxweight = -1.;
    for(size_t nv = 1 ; nv < datadist.size() ; nv++) {
        cout << nv <<  " " << datadist[nv] << " " << pileUpDist[nv] << endl;
        if(pileUpDist[nv] > 0. && datadist[nv]/pileUpDist[nv] > maxweight) {
            maxweight = datadist[nv]/pileUpDist[nv];
        }
    }
    return(maxweight);
}

//____________________________________________________________________________________
Double_t Analyse::GetPileUpWeight(vector<Double_t> &datadist) const {
    if(!usepileupinfo) return(1.);
    unsigned numtrueinteractions = unsigned(2*NumTruePileUpInteractions());
    Double_t dataval = 0.;
    if(datadist.size() > numtrueinteractions) {
        dataval = datadist[numtrueinteractions];
    }
    return(dataval/pileUpDist[numtrueinteractions]);
}


//____________________________________________________________________________________
Double_t Analyse::GetPrimVertexMaxWeight(vector<Double_t> &datadist) const {
    double maxweight = -1.;
    for(size_t nv = 1 ; nv < datadist.size() ; nv++) {
        //cout << nv <<  " " << datadist[nv] << " " << primVertexDist[nv] << endl;
        if(primVertexDist[nv] > 0. && datadist[nv]/primVertexDist[nv] > maxweight) {
            maxweight = datadist[nv]/primVertexDist[nv];
        }
    }
    return(maxweight);
}

//____________________________________________________________________________________
Double_t Analyse::GetPrimVertexWeight(vector<Double_t> &datadist) const {
    if(!useprimvertexinfo) return(1.);
    UInt_t numgoodprimvertices = NumGoodPrimVertices();
    Double_t dataval = 0.;
    if(datadist.size() > numgoodprimvertices) {
        dataval = datadist[numgoodprimvertices];
    }
    return(dataval/primVertexDist[numgoodprimvertices]);
}

//____________________________________________________________________________________
Int_t Analyse::NumGoodPrimVertices() const {
    Int_t numgoodprimvertices = 0;
    for(UInt_t i = 0 ; i < NumPrimVertices() ; i++) {
        if(PrimVertices(i).IsGood()) numgoodprimvertices++;
    }
    return(numgoodprimvertices);
}

Long64_t mem_usage() {
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    Long64_t rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();
    return(rss*sysconf(_SC_PAGE_SIZE));
}
