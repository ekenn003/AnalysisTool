#include "Analyse.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include <TH1D.h>

using namespace std;

class MyAnalysis : public Analyse {
    private:
    public:
        MyAnalysis() {}
        virtual ~MyAnalysis() {}
        virtual Int_t AnalyseEvent() {return(1);}
};

// _____________________________________________________________
int main() {
    TTree::SetMaxTreeSize(4830000000);
    string namebuf;
    string filename;
    string dirname;
    MyAnalysis ana;
    while(true) {
        cin >> namebuf;
        if (namebuf == "END") break;
        else {
            cout << namebuf << endl;
            UInt_t slashpos = namebuf.find_last_of("/");
            if (slashpos == namebuf.size()) {
                filename = namebuf;
            } else {
                filename = namebuf.substr(slashpos+1);
                dirname = namebuf.substr(0,slashpos+1);
            }
            if (filename.find("LUMI_") == 0) continue;
            else ana.AddLumiFile(namebuf.c_str(), "makeroottree/");

        }
    }
    ana.ResetLumiValues();
    cout << "Reading Lumi from text file" << endl;
    fstream file("lumi.cvs");
    char line[1000];
    Float_t lumirecorded;
    Float_t avgpu;
    UInt_t run;
    UInt_t block;
    if (file.is_open()) {
        file.getline(line, 1000);
        while (!file.eof()) {
            file.getline(line, 1000);
            vector<std::string> strs;
//            file.getline(line,1000);
//            vector<std::string> strs;
//            boost::split(strs, line, boost::is_any_of(","));
//            if(strs.size() != 8) continue;
//            run = atoi(strs[0].c_str());
//            vector<std::string> blocks;
//            boost::split(blocks, strs[1], boost::is_any_of(":"));
//            block = atoi(blocks[0].c_str());
//            lumirecorded = atof(strs[6].c_str());
//            avgpu = atof(strs[7].c_str());

            boost::split(strs, line, boost::is_any_of("|"));
            if (strs.size() != 6) continue;
            vector<std::string> runfill;
            boost::split(runfill, strs[0], boost::is_any_of(":"));
            run = atoi(runfill[0].c_str());
            block = atoi(strs[2].c_str()); // not sure if this should be "nls" or "ncms"
            lumirecorded = atof(strs[5].c_str());
            if (lumirecorded == 0.) {
                cout << "WARNING Lumi value is 0. Run: " << run << ", Block: " << block << endl;
                continue;
            }
//            if (ana.SetLumi(run, block, lumirecorded/1000000., avgpu) == 1) {
//                //muhist->Fill(atof(strs[7].c_str())*766./693.);
//                muhist->Fill(avgpu);
//            }
        }
        file.close();
    } else {
        cout << "ERROR: lumi.cvs could not be opened!" << endl;
    }

    TH1D* muhist = new TH1D("mudist", "mudist", 200, 0., 100.);
    TFile *pufile = new TFile("MyDataPileupHistogram.root","READ");
    TH1 *puhist = (TH1D*)(pufile->Get("pileup")->Clone());
    for(int i = 0; i < puhist->GetNbinsX(); ++i) {
        muhist->SetBinContent(i, puhist->GetBinContent(i));
    }
    pufile->Close();


    ana.PrintLumiOfRuns();
    ana.WriteLumiFile(dirname + string("LUMI_INFO.root"));
    TFile* lumfile = new TFile((dirname + string("LUMI_INFO.root")).c_str(), "update");
    muhist->Write();
    lumfile->Write();
    lumfile->Close();
}

