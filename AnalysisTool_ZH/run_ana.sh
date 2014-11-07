#!/bin/bash
source bin/AnalysisToolUseThis
rm *.o
rm output.txt
make
touch output.txt
#./myanalysis < input_ZH2mu.txt > output.txt

