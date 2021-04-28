set -e
set -u
set -x

inGenome="$1"
inTaxID="$2"
myPath=${0/blastDbSetup.sh/}

genomeName=$(basename ${inGenome})
baseGenomeName=${genomeName%%.fa?(sta)}

if [ genomeName == baseGenomeName ]; then
    echo "Genome must be a .fa or .fasta file"
    exit 1
fi

makeblastdb \
-dbtype nucl \
-title ${baseGenomeName} \
-out ${baseGenomeName}_db \
-in ${inGenome} \
-taxid ${inTaxID}

blastdb_aliastool -dblist "${baseGenomeName}_db" \
-dbtype nucl -out ${baseGenomeName}_db -title "${baseGenomeName}"
