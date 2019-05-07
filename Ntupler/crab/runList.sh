#!/bin/bash

file=$1
config=$2
isData=`echo $file | grep -c data`
isSig=`echo $file | grep -c sig`
is8X=`echo $file | grep -c 8X`
is9X=`echo $file | grep -c 9X`
is10X=`echo $file | grep -c 10X`
for x in `cat $file | grep -v "#"`; do 
    ext='_ext'
    isExt=`echo $x | grep -c ext`
    if [ $isExt -eq 0 ]; then
	ext=''
    fi
    label=`echo $x   | sed "s@/@ @g" | awk '{print $1_$2}'  | sed "s@-madgraphMLM-pythia8@@g" | sed "s@-@_@g"`
    label=${label}${ext}
    if [ $is8X -ne 0 ]; then
        label=${label}_8X
    fi
    if [ $is10X -ne 0 ]; then
        label=${label}_10X
    fi
    rm -r -f crab_projects/crab_${label}
    echo ${label} $x $config  1
    #./run.sh ${label} $x $config  1
done