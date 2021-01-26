set -x

inGenome="$1"
inTaxID="$2"
myPath=${0/blastDbSetup.sh/}

makeblastdb \
-dbtype nucl \
-title ${inGenome/.fa*/} \
-out ${inGenome/.fa*/}_db \
-in ${inGenome} \
-taxid ${inTaxID}

blastdb_aliastool -dblist "${inGenome/.fa*/}_db" \
-dbtype nucl -out ${inGenome/.fa*/}_db -title "${inGenome/.fa*/}"
