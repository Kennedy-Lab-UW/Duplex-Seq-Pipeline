#!/bin/bash

inAmbigReads=$1
inNonAmbigReads=$2
outRecoveredReads=$3
basePath=$4

ln -s ${inNonAmbigReads} ${outRecoveredReads}