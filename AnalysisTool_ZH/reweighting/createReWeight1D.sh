#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session lol

INDIR_DYJetsToLL=/afs/cern.ch/user/e/ekennedy/eos/cms/store/user/ekennedy/smh2mu/mc/jan16/DYJetsToLL_jan16/DYJetsToLLtest/AC1B*.root
INDIR_TTJets=/afs/cern.ch/user/e/ekennedy/eos/cms/store/user/ekennedy/smh2mu/mc/jan16/TTJets_jan16/TTJetstest/AC1B*.root

DATAPU=../PileUpData_69000.root

##echo "1. make sure eos is mounted"
##echo "2. make sure there are no typos in dir's"
##echo "3. make sure sample name has no spaces"
##
##
##root -q -b "infoFiles.C(\"$INDIR_DYJetsToLL\",\"DYJetsToLL\")"
##root -q -b "infoFiles.C(\"$INDIR_TTJets\",\"TTJets\")"
##hadd -f nPV_trees.root ReWeightInfo_*.root
##rm ReWeightInfo_*.root

#root -q -b "createReWeightPieces.C(\"DYJetsToLL\",$DATAPU)"
#root -q -b "createReWeightPieces.C(\"TTJets\",$DATAPU)"
root -q -b "createReWeightPieces.C(\"DYJetsToLL\",\"$DATAPU\")"
root -q -b "createReWeightPieces.C(\"TTJets\",\"$DATAPU\")"
#rm nPV_trees.root
hadd -f ReWeight1DWhole.root ReWeightPiece*.root $DATAPU
rm ReWeightPiece*.root

