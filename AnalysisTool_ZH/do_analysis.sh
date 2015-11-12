#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session

# check that eos is mounted
if [ ! -d "/afs/cern.ch/user/e/ekennedy/eos/cms/store/user/ekennedy/" ]; then
    echo "ERROR: can't find EOS directory - mount it please"
    exit 1
fi

# make some files
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
    echo "emptying $OUTDIR in 5 seconds..."
    sleep 1
    echo "emptying $OUTDIR in 4 seconds..."
    sleep 1
    echo "emptying $OUTDIR in 3 seconds..."
    sleep 1
    echo "emptying $OUTDIR in 2 seconds..."
    sleep 1
    echo "emptying $OUTDIR in 1 second..."
    sleep 1
    rm -r $OUTDIR
    mkdir $OUTDIR
fi

# set up analysis software
if [ -f "main.cc" ]; then
    rm main.cc
fi
if [ "$FSTATE" == "4mu" ]; then
    ln -s ana_Z2muH2mu_main.cc main.cc
elif [ "$FSTATE" == "2e2mu" ]; then
    ln -s ana_Z2eH2mu_main.cc main.cc
elif [ "$FSTATE" == "2mu" ]; then
    ln -s ana_2mu_main.cc main.cc
else
    echo "ERROR: final state not found :("
    exit 1
fi
if [ -f "*.o" ]; then
    rm *.o
fi
make


# if LUMI_INFO.root doesn't exist for these input files, make it
if [ ! -f "$INDIR/LUMI_INFO.root" ]; then
    echo "no lumi info file"
    # put the input file names into the input text file
    ls $INDIR/*.root > $OUTDIR/$INFILE
    echo "END" >> $OUTDIR/$INFILE
    ./../AnalysisTool/merger < $OUTDIR/$INFILE && echo "lumi info file made"
    rm $OUTDIR/$INFILE
fi

# put the input file names into the input text file
ls $INDIR/*.root > $OUTDIR/$INFILE
echo "END" >> $OUTDIR/$INFILE

# run the analysis
CURDIR=$PWD
cd $OUTDIR
$CURDIR/myanalysis < $INFILE > $OUTFILE && echo "results are in $PWD"

