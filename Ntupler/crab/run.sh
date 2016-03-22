#!/bin/bash

label=$1
dataset=$2
config=$3

#crab status crab_projects/crab_${label}
#crab resubmit crab_projects/crab_${label}
#crab kill crab_projects/crab_${label}
#exit
rm -rf crab_projects/crab_${label}
echo deleting eos folder /store/group/cmst3/group/monojet/production/singleb/${label}
/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select rm -r /store/group/cmst3/group/monojet/production/singleb/${label}
echo creating eos folder /store/group/cmst3/group/monojet/production/singleb/${label}
/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select mkdir /store/group/cmst3/group/monojet/production/singleb/${label}
echo eos folder /store/group/cmst3/group/monojet/production/singleb/${label} created
if [ $config == "makingBacon_Data_25ns_MINIAOD.py" ]; then
    sed "s@XX-LABEL-XX@$label@g" crab_template.py | sed "s@XX-DATASET-XX@$dataset@g" | sed "s@XX-CONFIG-XX@$config@g"   > crab.py
fi
if [ $config == "makingBacon_MC_25ns_MINIAOD.py" ]; then
    sed "s@XX-LABEL-XX@$label@g" crab_template_mc.py | sed "s@XX-DATASET-XX@$dataset@g" | sed "s@XX-CONFIG-XX@$config@g"   > crab.py
fi
if [ $config == "makingBacon_MC_25ns_MINIAOD_Gen.py" ]; then
    sed "s@XX-LABEL-XX@$label@g" crab_template_mc.py | sed "s@XX-DATASET-XX@$dataset@g" | sed "s@XX-CONFIG-XX@$config@g"   > crab.py
fi
crab submit -c crab.py
mv crab.py old/crab_${label}.py
