#!/bin/bash

inAmbigReads=$1
inNonAmbigReads=$2
inWrongSpeciesReads=$3
outPrefix=$4
basePath=$5

cp ${inNonAmbigReads} ${outPrefix}.recovered.temp.bam
cp ${inAmbigReads} ${outPrefix}.ambig.bam
cp ${inWrongSpeciesReads} ${outPrefix}.wrongSpecies.bam