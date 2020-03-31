#!/bin/bash

inAmbigReads=$1
inNonAmbigReads=$2
outRecoveredReads=$3
basePath=$4
echo $(pwd)
python3 "${basePath}/scripts/RecoveryScripts/SingleChromosomeTargetRecovery.py" ${inAmbigReads} ${inAmbigReads}.recovered.temp.bam chrM
samtools merge -c ${outRecoveredReads} ${inNonAmbigReads} ${inAmbigReads}.recovered.temp.bam
rm ${inAmbigReads}.recovered.temp.bam