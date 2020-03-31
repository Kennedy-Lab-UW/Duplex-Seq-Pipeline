set -x

inGenome="$1"
inTaxID="$2"
myPath=${0/blastDbSetup.sh/}

python3 ${myPath}/AddTaxonID.py ${inGenome} ${inTaxID} ${inGenome/.fa/_taxID.fa}

makeblastdb \
-dbtype nucl \
-title ${inGenome/.fa*/} \
-out ${inGenome/.fa*/}_db \
-in ${inGenome/.fa/_taxID.fa}

echo "#${inGenome/.fa*/}.nal"  > ${inGenome/.fa*/}.nal
echo "TITLE ${inGenome/.fa*/}" >> ${inGenome/.fa*/}.nal
echo "DBLIST ${inGenome/.fa*/}_db" >> ${inGenome/.fa*/}.nal

