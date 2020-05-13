#!/bin/bash

# This script initializes the snakemake envs for Duplex Sequencing.  Before running this script, make sure you have downloaded a copy of GATK 3.8.1.  The GATK3.8.1 jar file is the first argument of this script.  
set -e
set -o pipefail
set -u

#~ #checking that input is a valid gatk 3.8 jar file
#~ jar_version=$(java -jar "${1}" --version | grep -oEi '[0-9]\.[0-9]' | grep -oEim 1 '[0-9]\.[0-9]')
#~ if [[ "$jar_version" != "3.8" ]]; then
    #~ echo "This file is not version 3.8, but $jar_version.  Download GATK 3.8.1 and try again."
    #~ exit 1
#~ fi

snakeDir=$(pwd)
maxCores="${1}"

# Setup test case
echo "Creating test config file"
echo "sample,rglb,rgpl,rgpu,rgsm,reference,target_bed,blast_db,targetTaxonId,baseDir,in1,in2,mqFilt,minMem,maxMem,cutOff,nCutOff,umiLen,spacerLen,locLen,readLen,clipBegin,clipEnd,minClonal,maxClonal,minDepth,maxNs,runSSCS,recovery"> test/testConfig.csv
echo "test1,test,test,test,test,${snakeDir}/test/testRef/testRef.fa,${snakeDir}/test/testTarget/test.bed,${snakeDir}/test/testBlastDb/testBlastDb,9606,testData,testSeq1.fastq.gz,testSeq2.fastq.gz,0,3,200,0.7,0.02,8,1,8,150,7,0,0,0.1,100,1,FALSE,noRecovery.sh" >> test/testConfig.csv
echo "test2,test,test,test,test,${snakeDir}/test/testRef/testRef.fa,${snakeDir}/test/testTarget/test.bed,${snakeDir}/test/testBlastDb/testBlastDb,9606,testData,testSeq1.fastq.gz,testSeq2.fastq.gz,0,3,200,0.7,0.02,8,1,8,150,7,0,0,0.1,100,1,FALSE,recoverAmbig.sh" >> test/testConfig.csv
echo "test3,test,test,test,test,${snakeDir}/test/testRef/testRef.fa,${snakeDir}/test/testTarget/test.bed,${snakeDir}/test/testBlastDb/testBlastDb,9606,testData,testSeq1.fastq.gz,testSeq2.fastq.gz,0,3,200,0.7,0.02,8,1,8,150,7,0,0,0.1,100,1,FALSE,recoverWrongSpecies.sh" >> test/testConfig.csv
echo "test4,test,test,test,test,${snakeDir}/test/testRef/testRef.fa,${snakeDir}/test/testTarget/test.bed,${snakeDir}/test/testBlastDb/testBlastDb,9606,testData,testSeq1.fastq.gz,testSeq2.fastq.gz,0,3,200,0.7,0.02,8,1,8,150,7,0,0,0.1,100,1,FALSE,recoverAll.sh" >> test/testConfig.csv
echo "test5,test,test,test,test,${snakeDir}/test/testRef/testRef.fa,${snakeDir}/test/testTarget/test.bed,none,9606,testData,testSeq1.fastq.gz,testSeq2.fastq.gz,0,3,200,0.7,0.02,8,1,8,150,7,0,0,0.1,100,1,FALSE,noRecovery.sh" >> test/testConfig.csv


# Set up progConfig file
echo "Creating progConfig file"
echo "gatk3:" > DS_progConfig.yaml
echo "samples: test/testConfig.csv" >> DS_progConfig.yaml
echo "maxCores: ${maxCores}" >> DS_progConfig.yaml
echo "vardict_f: \".0000001\"" >> DS_progConfig.yaml
echo "vardict_nmfreq: \".0000001\"" >> DS_progConfig.yaml
echo "vardict_r: \"1\"" >> DS_progConfig.yaml
echo "vardict_V: \"0.00000000001\"" >> DS_progConfig.yaml
echo "vardict_adaptor: GCTCTTCCGATCT,CTCTTCCGATCT,TCTTCCGATCT,CTTCCGATCT,TTCCGATCT,TCCGATCT,CCGATCT,CGATCT" >> DS_progConfig.yaml


echo "Configuring snakemake"
snakemake --cores 1 --use-conda --conda-prefix ${snakeDir}/.snakemake --config gatk3=${1} -- initializeEnvs

echo "Creating run script"
echo "#!/bin/bash" > DS
echo "" >> DS
echo "# This is a run script for the DS snakemake pipeline" >> DS
echo "inConfig=\"\$1\"" >> DS
echo "snakemake -s ${snakeDir}/Snakefile --use-conda -j ${maxCores} --conda-prefix ${snakeDir}/.snakemake --config samples=\"\${inConfig}\"" >> DS
chmod a+x DS

echo "Creating dag script"
echo "#!/bin/bash" > DS-dag
echo "" >> DS-dag
echo "# This is a run script for the DS snakemake pipeline" >> DS-dag
echo "inConfig=\"\$1\"" >> DS-dag
echo "snakemake -s ${snakeDir}/Snakefile --use-conda -j ${maxCores} --dag --conda-prefix ${snakeDir}/.snakemake --config samples=\"\${inConfig}\" -- | dot -Tpdf > \${inConfig}_dag.pdf" >> DS-dag
chmod a+x DS-dag

echo "Creating unlock script"
echo "#!/bin/bash" > DS-unlock
echo "" >> DS-unlock
echo "# This is an unlock script for the DS snakemake pipeline" >> DS-unlock
echo "# Run this if the pipeline gets stuck locked for some reason" >> DS-unlock
echo "inConfig=\"\$1\"" >> DS-unlock
echo "snakemake --unlock -s ${snakeDir}/Snakefile --use-conda -j 1 --conda-prefix ${snakeDir}/.snakemake --config samples=\"\${inConfig}\"" >> DS-unlock
chmod a+x DS-unlock

currentDate=$(date +%F)
echo "Adding DS to path"
if [ -e ~/.bashrc ]
then
cp ~/.bashrc ~/.bashrc_backup_${currentDate}
cat ~/.bashrc \
| sed -e '/^#Duplex\ Sequencing\ Pipeline/,+1 s/^/#/' \
> ~/.bashrc.temp
echo "#Duplex Sequencign Pipeline" >> ~/.bashrc.temp
echo "export PATH=\"${snakeDir}:\$PATH\"" >> ~/.bashrc.temp
echo "" >> ~/.bashrc.temp
mv ~/.bashrc.temp ~/.bashrc
source ~/.bashrc
fi
if [ -e ~/.bash_profile ]
then
cp ~/.bash_profile ~/.bash_profile_backup_${currentDate}
cat ~/.bash_profile \
| sed -e '/^#Duplex\ Sequencing\ Pipeline/,+1 s/^/#/' \
> ~/.bash_profile.temp
echo "#Duplex Sequencign Pipeline" >> ~/.bash_profile.temp
echo "export PATH=\"${snakeDir}:\$PATH\"" >> ~/.bash_profile.temp
echo "" >> ~/.bash_profile
mv ~/.bash_profile.temp ~/.bash_profile
source ~/.bash_profile
fi

echo "Done"
