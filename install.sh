#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session
# check if we are in the right place
if [ ! ${PWD##*/} == "AnalysisTool" ]; then
    echo "ERROR: This script has to be run from AnalysisTool directory"
    exit 1
fi
CURDIR=$PWD
cd AnalysisTool
./configure --prefix=$CURDIR/AnalysisTool_ZH
make clean
make
make install
source AnalysisToolUseThis
echo ""
echo "Now you can run this tool from the directory AnalysisTool/AnalysisTool_ZH"
echo "You will have to make a softlink to the analysis code that you want to use and then do \"make\" but if you use the run scripts, these steps are already included for you."
cd $CURDIR/AnalysisTool_ZH/
