#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session

# make some files
OUTDIR=results/$WNAME\_$FSTATE
INFILE=input_$WNAME.txt
OUTFILE=output_$WNAME\_$FSTATE.txt
# check if we are in the right place
if [ ! ${PWD##*/} == "AnalysisTool_ZH" ]; then
    echo "ERROR: This script has to be run from AnalysisTool/AnalysisTool_ZH"
    exit 1
fi
# make the working directory if it isn't there
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
else
    echo "WARNING: $OUTDIR already exists"
fi

# set up analysis software
if [ -f "main.cc" ]; then
    rm main.cc
fi
if [ "$FSTATE" == "4mu" ]; then
    ln -s ana_Z2muH2mu_main.cc main.cc
elif [ "$FSTATE" == "2e2mu" ]; then
    ln -s ana_Z2eH2mu_main.cc main.cc
else
    echo "ERROR: final state not found :("
    exit 1
fi
if [ -f "*.o" ]; then
    rm *.o
fi
make

# check that eos is mounted
if [ ! -d "/afs/cern.ch/user/e/ekennedy/eos/cms/store/user/ekennedy/" ]; then
    echo "ERROR: can't find EOS directory - mount it please"
    exit 1
fi
# put the input file names into the input text file
ls $INDIR/*.root > $OUTDIR/$INFILE
#if LUMI_INFO.root doesn't exist for these input files, make it
if [ ! -f "$INDIR/LUMI_INFO.root" ]; then
    echo "no lumi info file... making it"
    echo "END" >> $OUTDIR/$INFILE
    ./../AnalysisTool/merger < $OUTDIR/$INFILE && echo "lumi info file made "
    rm $OUTDIR/$INFILE
    ls $INDIR/*.root > $OUTDIR/$INFILE
fi
echo "END" >> $OUTDIR/$INFILE
# remove the output file if it's there already
if [ -f "$OUTFILE" ]; then
    rm $OUTFILE
fi
# run the analysis
CURDIR=$PWD
cd $OUTDIR
$CURDIR/myanalysis < $INFILE > $OUTFILE && echo "results are in $PWD"
