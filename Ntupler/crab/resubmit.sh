#!/bin/bash

label=$1
dataset=$2
config=$3

crab status   crab_projects/crab_${label}
crab resubmit crab_projects/crab_${label}
