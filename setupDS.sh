#!/bin/bash

# This script initializes the snakemake envs for Duplex Sequencing.  Before running this script, make sure you have downloaded a copy of GATK 3.8.1.  The GATK3.8.1 jar file is the first argument of this script.  
set -e
set -o pipefail
set -u

#checking that input is a valid gatk 3.8 jar file
jar_version=$(java -jar "${1}" --version | grep -oEi '[0-9]\.[0-9]' | grep -oEim 1 '[0-9]\.[0-9]')
if [[ "$jar_version" != "3.8" ]]; then
    echo "This file is not version 3.8, but $jar_version.  Download GATK 3.8.1 and try again."
    exit 1
fi

snakeDir=$(pwd)
echo "Configuring snakemake"
snakemake --use-conda --conda-prefix ${snakeDir}/.snakemake --config gatk3=${1} -- initializeEnvs

echo "Creating run script"
echo "#!/bin/bash" > DS
echo "" >> DS
echo "# This is a run script for the DS snakemake pipeline" >> DS
echo "inConfig=\"\$1\"" >> DS
echo "snakemake -s ${snakeDir}/Snakefile --use-conda -j 6 --conda-prefix ${snakeDir}/.snakemake --config samples=\"\${inConfig}\"" >> DS
chmod a+x DS

echo "Creating dag script"
echo "#!/bin/bash" > DS-dag
echo "" >> DS-dag
echo "# This is a run script for the DS snakemake pipeline" >> DS-dag
echo "inConfig=\"\$1\"" >> DS-dag
echo "snakemake -s ${snakeDir}/Snakefile --use-conda -j 6 --dag --conda-prefix ${snakeDir}/.snakemake --config samples=\"\${inConfig}\" -- | dot -Tpdf > \${inConfig}_dag.pdf" >> DS-dag
chmod a+x DS-dag

currentDate=$(date +%F)
echo "Adding DS to path"
if [ -e ~/.bashrc ]
then
cp ~/.bashrc ~/.bashrc_backup_${currentDate}
echo "export PATH=\"${snakeDir}:\$PATH\"" >> ~/.bashrc
echo "" >> ~/.bashrc
fi
if [ -e ~/.bash_profile ]
then
cp ~/.bash_profile ~/.bash_profile_backup_${currentDate}
echo "export PATH=\"${snakeDir}:\$PATH\"" >> ~/.bash_profile
echo "" >> ~/.bash_profile
fi

echo "Done"
