#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session

#INDIR=/afs/cern.ch/user/e/ekennedy/eos/cms/store/user/ekennedy/smh2mu/mc/april2015sets/WH_ZH_HToMuMu/
INDIR=/afs/cern.ch/work/e/ekennedy/work/tuplizer/miniAOD/CMSSW_7_2_4/src/april2015sets/WH_ZH_HToMuMu/rootfiles
WNAME=WH_ZH_HToMuMu_set_test
FSTATE=4mu
OUTDIR=results/april2015set/$WNAME\_$FSTATE

source do_analysis.sh
