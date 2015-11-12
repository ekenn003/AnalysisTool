#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session

INDIR=/afs/cern.ch/user/e/ekennedy/eos/cms/store/user/ekennedy/smh2mu/mc/25ns/ZZTo4L_miniAODv2
WNAME=Run2MC_ZZTo4L_mv2
FSTATE=2mu
OUTDIR=results/25ns/$WNAME\_$FSTATE

source do_analysis.sh
