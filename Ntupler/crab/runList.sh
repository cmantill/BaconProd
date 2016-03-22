#!/bin/bash

file=$1
config=makingBacon_MC_25ns_MINIAOD_Gen.py
for x in `cat $file`; do 
    label=`echo $x   | sed "s@/@ @g" | awk '{print $1}'  | sed "s@-madgraphMLM-pythia8@@g"`
    label=`echo $label | sed "s@TuneCUETP8M1_@@g"`
    label=`echo $label | sed "s@_13TeV-powheg@@g"`
    label=`echo $label | sed "s@DMV_NNPDF30_@@g"`
    label=`echo $label | sed "s@DMS_NNPDF30_@@g"`
    label=`echo $label | sed "s@_gSM-1p0_gDM-1p0@@g"`
    #label=`echo $label | sed "s@_gSM-0p25_gDM-1p0@@g"`
    #crab status crab_projects/crab_${label}
    #crab kill crab_projects/crab_${label}
    ./run.sh Spring15_a25ns_DMJets$label $x $config
    #./resubmit.sh Spring15_a25ns_DMJets$label $x $config
done