#!/bin/bash

inAmbigReads=$1
inNonAmbigReads=$2
inWrongSpeciesReads=$3
outPrefix=$4
basePath=$5

echo $(pwd)
python3 \
"${basePath}/scripts/RecoveryScripts/SingleChromosomeTargetRecovery.py" \
${inAmbigReads} \
${inAmbigReads}.recovered.temp.bam \
${outPrefix}.ambig.bam \
chrM
samtools merge -c ${outPrefix}.recovered.bam\
 ${inNonAmbigReads} \
 ${inAmbigReads}.recovered.temp.bam
rm ${inAmbigReads}.recovered.temp.bam
cp ${inWrongSpeciesReads} ${outPrefix}.wrongSpecies.bam