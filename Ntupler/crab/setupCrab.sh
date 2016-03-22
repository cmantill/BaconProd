#!/bin/bash

eval `scramv1 runtime -sh`
source /cvmfs/cms.cern.ch/crab3/crab.sh
crab --version
cat $HOME/private/$USER.txt | voms-proxy-init -voms cms --valid 168:00 -pwstdin
voms-proxy-info --all