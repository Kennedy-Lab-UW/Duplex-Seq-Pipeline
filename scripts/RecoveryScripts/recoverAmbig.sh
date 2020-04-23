#!/bin/bash

inAmbigReads=$1
inNonAmbigReads=$2
inWrongSpeciesReads=$3
outPrefix=$4
basePath=$5

samtools merge -c \
${outPrefix}.recovered.temp.bam \
${inNonAmbigReads} \
${inAmbigReads}

samtools view -Hb ${outPrefix}.recovered.temp.bam \
> ${outPrefix}.ambig.bam
cp ${inWrongSpeciesReads} \
${outPrefix}.wrongSpecies.bam