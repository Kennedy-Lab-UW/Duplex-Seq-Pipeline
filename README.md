Duplex Sequencing Pipeline
======================

Brendan Kohrn, October 20, 2020  
Duplex Sequencing is copyright Scott Kennedy and Larwrence Loeb, 
University of Washington, Seattle, WA, USA  

## Table of Contents:

1. Glossary
2. Dependencies
3. Setup
4. Genome Setup
5. Contaminant database setup
6. Bed and Interval_list file preparation
7. Configuration File creation
8. Recovery script creation
9. Running pipelines
10. Output file descriptions
11. Extra BAM tags:
12. Testing the pipeline
13. Full and partial reruns
14. Unlocking following a power failure

## 1: Glossary

### Single Stranded Consensus Sequence (SSCS)

A construct created by comparing multiple reads and deciding ambiguities by
simple majority.

### Duplex Consensus Sequence (DCS)

A construct created by comparing two SSCSs.

### Family

A group of reads that shares the same tag sequence.

## 2: Dependencies:
This pipeline is known to work with the following minimum versions of 
the follow required programs:

* Python3.6+
* Snakemake>=5.25.0
* Pandas
* Miniconda/Anaconda=4.7.\*
* bwa=0.7.17.* (for genome setup)
* ncbi-blast=>2.6.0 (installed separately, for contaminant database setup)
* mamba>=0.5.1 (needed for environment setup; install through conda)
* wget (on macOS, install using homebrew; present by default on linux)

Once Python3.6 is installed, snakemake and pandas can be installed using 
pip3 or using whatever package manager you're using.  Blast can be 
downloaded in any of several ways, including some package managers 
(Ubuntu: sudo apt-get install ncbi-blast+).  It can also be installed 
using conda if desired.  

## 3: Setup: 

Find the location where you want the pipeline to be located and clone the 
pipeline using git:

```bash
git clone https://github.com/KennedyLabUW/Duplex-Seq-Pipeline.git
```

change into the directory, and run:

```bash
bash setupDS.sh MAX_CORES
```

where MAX_CORES is the maximum number of cores you want the pipeline to 
be able to use. This will create a test configuration file in the test 
directory (test/testConfig.csv), a basic configuration file (DS_progConfig.yaml), 
a run script (DS), a DAG-creation script (DS-dag), a reset script (DS-clean), 
and an unlock script (DS-unlock). 
After this, with the exception of setting up the Genomes 
(See Section 4) and the (optional) blast contamination database (See 
Section 5), you should be able to run the Duplex-Seq pipeline using 

```bash
bash /path/to/Duplex-Seq-Pipeline/DS CONFIG_CSV.csv
```

where CONFIG_CSV.csv is a configuration CSV file generated as described 
below.  

If you want, you can also add the DS and DS-clean scripts to your path; 
this will simplify invocation if you only have a single release version 
on your system. 

## 4: Genome setup

Put genomes in an easily findable location, such as our references 
directory (i.e. \~/bioinformatics/reference). 

Many genomes can be downloaded from UCSC 
(http://hgdownload.soe.ucsc.edu/downloads.html).  In order to download 
genomes from there, you will need to download the twoBitToFa program from 
the appropriate Utilities directory.  twoBitToFa has the following syntax 
(Copied from UCSC):

```
twoBitToFa - Convert all or part of .2bit file to fasta
usage:
   twoBitToFa input.2bit output.fa
options:
   -seq=name       Restrict this to just one sequence.
   -start=X        Start at given position in sequence (zero-based).
   -end=X          End at given position in sequence (non-inclusive).
   -seqList=file   File containing list of the desired sequence names 
                   in the format seqSpec[:start-end], e.g. chr1 or chr1:0-189
                   where coordinates are half-open zero-based, i.e. [start,end).
   -noMask         Convert sequence to all upper case.
   -bpt=index.bpt  Use bpt index instead of built-in one.
   -bed=input.bed  Grab sequences specified by input.bed. Will exclude introns.
   -bedPos         With -bed, use chrom:start-end as the fasta ID in output.fa.
   -udcDir=/dir/to/cache  Place to put cache for remote bigBed/bigWigs.

Sequence and range may also be specified as part of the input file name using the syntax:

      /path/input.2bit:name
   or
      /path/input.2bit:name
   or
      /path/input.2bit:name:start-end
```

Once you have downloaded a genome and converted it into FASTA format, it 
needs to be indexed.  To do this, open a terminal and navigate to your 
the directory containing your genome.  Note that in the following commands, 
the word "genome.fasta" should be replaced with the file name of your 
genome. These commands **must** be run in the same directory as the 
reference genome

```
    bwa index genome.fasta  
    samtools faidx genome.fasta  
    # /Path/To/PicardTools should be the path to wherever you put your picard tools jar file.    
    # Also, genome.dict should match the name of your fasta (e.g. hg19.fasta to hg19.dict)  
    java -jar /Path/To/PicardTools/picard.jar CreateSequenceDictionary R=genome.fasta O=genome.dict  
```

At the moment, this pipeline does not support compressed genomes.

## 5: Contaminant Database Setup:
The Duplex-Seq pipeline is designed to use a local NCBI Blast instance 
to detect and remove potential contamination from non-target species and 
identify issues arising from pseudogenes. 
*This step is optional, but requires a valid NCBI Blast database file 
name. To run without BLAST, enter "NONE" (with any capitalization) 
in the blast_db field in the config.csv file.*
 
To construct your contaminant database, if desired, first decide on a 
list of species you want to monitor for contaminants.  A suggested 
starting list is:  

 * Human (hg38)  
 * Mouse (mm10)  
 * Rat (rn6)  
 * C. elegens (ce11)  
 * Yeast (sacCer3)  
 * Fruit Fly (dm6)  
 * Cow (bosTau9)  
 * Dog (canFam3)  
 
**It is important that any genome you plan to use for alignment be included in this database, 
in the same version (e.g. if your alignment genome uses UCSC chromosome 
names, your genome in the database cannot use NCBI chromosome names, but 
must also use UCSC chromosome names)**.  

Database setup consists of three steps:

 1. Genome labeling:

Label each record in the genome with the NCBI taxon ID for the species the 
genome is associated with.  TaxIDs can be found using the NCBI taxonomy 
website (https://www.ncbi.nlm.nih.gov/taxonomy).  This is done by running:

```
python3 AddTaxonID.py GENOME.fa TAXID GENOME_taxID.fa
```

where GENOME.fa is the input fasta file with the genome, TAXID is the 
NCBI taxonomy ID for the species associated with the genome, and 
GENOME_taxID.fa is the output labeled genome.  

 2. Sub-database Creation: 

Create the database using:

```
makeblastdb \  
-dbtype nucl \  
-title GENOME \  
-out GENOME_db \  
-in GENOME_taxID.fa  
```

After this, create a .nal file for this database following this template:

```
#GENOME.nal 
TITLE GENOME
DBLIST GENOME_db
```

The blastDbSetup.sh script can be used to automate these steps.  It can 
be run with:

```
bash /path/to/pipeline/setupBlastDb/blastDbSetup.sh GENOME.fa TAXON_ID
```

 3. Full database creation:

Once you've created all your sub-databases (e.g. GENOME1, GENOME2, 
GENOME3, ..., GENOME_N), create a .nal file to represent the full 
database following this template:

```
#contaminantDb.nal 
TITLE contaminantDb
DBLIST GENOME1_db GENOME2_db GENOME3_db ... GENOME_N_db
```

If, at a later time, you need to change your contaminant database, you 
then only need to rebuild the modified portion, and then update this 
file to reflect that.  

## 6: Bed file preparation

A bed file is a file which details regions of the genome in which you 
are interested.  The syntax for bed files is described at 
https://genome.ucsc.edu/FAQ/FAQformat.html#format1.  

If you know which genes you are targeting and are using a common 
published genome, the target bed file can be downloaded from the UCSC Table 
Browser (https://genome.ucsc.edu/cgi-bin/hgTables).  Otherwise, it can 
be created using results from a BLAST search or any other method you 
like.  **This pipeline does not currently support overlapping intervals, but
does support blocks as described in the bed spec.**  

This pipeline can also optionally use a masking bed file, for if there are 
genomic regions where you know there will be a high number of artifactual 
variants, or which you would like to ignore for other reasons. The masking bed 
file only needs to contain the first three bed columns (chrom, start, and stop), 
and will ignore any other columns provided.  

## 7: Configuration file creation:

Use the ConfigTemplate to create a new file with the appropriate headers. 
For each row, fill in the information about a particular sample:

| Header           | Required or Default    | Information |
| ---------------- | ---------------------- | ----------- |
| sample           | Required               | A unique identifier for a sample; this will be used to name all output files for this sample |
| rglb             | Required               | Read Group Library Identifier |
| rgpl             | Required               | Read Group Platform; usually Illumina |
| rgpu             | Required               | Read Group Platform Unit |
| rgsm             | Required               | Read Group Sample |
| reference        | Required               | The path to the prepared reference genome to use with this sample.  |
| target_bed       | Required               | A bed file showing where the targets are for this particular sample |  
| maskBed | NONE | A bed file to use for masking variants. |  
| blast_db         | Required               | The blast database to use for contaminant filtering; must include your target genome.  |
| targetTaxonId    | Required               | The taxon ID of the species you are expecting to be present in the sample.  |
| baseDir          | Required               | The directory the input files are in, and where the output files will be created. |
| in1              | Required               | The read1 fastq (or fastq.gz, or fq.gz, or fq) file for this sample. Note that this is just the name of the file, and not the full path.  |
| in2              | Required               | The read2 fastq (or fastq.gz, or fq.gz, or fq) file for this sample. Note that this is just the name of the file, and not the full path.   |
| mqFilt           | 0                      | A threshold for mapping quality filtering, if desired. |
| minMem           | 0                      | The minimum number of reads that must be in a family for consensus making |
| maxMem           | 200                    | The maximum number of reads in a family the consensus maker should consider. |
| cutOff           | 0.9                    | The threshold for consensus making; the consensus maker will require at least this much agreement on a per base pair level. |
| nCutOff          | 1                      | The maximum proportion of N bases in an output consensus sequence. |
| umiLen           | 8                      | The length of the UMI in this sample |
| spacerLen        | 1                      | The length of the spacer sequence in this sample |
| locLen           | 10                     | The localization length to use for this sample |
| readLen          | 101                    | The length of a read for this sample |  
| adapterSeq       | "ANNNNNNNNAGATCGGAAGAG" | The adapter sequence used in library preperation, with UMI bases as Ns, and spacer sequence included.  Used by cutadapt for adapter clipping |  
| clipBegin        | 7                      | How many bases to clip off the 5' end of the read |
| clipEnd          | 0                      | How many bases to clip off the 3' end of the read |
| minClonal        | 0                      | The minimum clonality to use for count_muts generation |
| maxClonal        | 0.1                    | The maximum clonality to use for count_muts generation |
| minDepth         | 100                    | The minimum depth to use for count_muts generation |
| maxNs            | 1                      | The maximum proportion of N bases to use for count_muts generation |
| recovery         | "noRecovery.sh"        | The recovery script to use in attempting to recover ambiguously mapped reads (as determine by blast alignment vs bwa alignment).  Recovery script creation is discussed in 8; below.  |  
| cm_outputs       | "GB"                   | Select which sections of the countmuts to output, in addition to 'OVERALL'.  String of one or more of 'G', 'B', and 'N'.  G -> output GENE sections for each bed line; B -> output 'BLOCK' sections for each block in the bed line (if present); 'N' -> Only output overall frequencies.  Overrides all other options. |  
| cm_sumTypes      | "GT"                   | How to calculate OVERALL and GENE blocks for countmuts output. The first character controls summing for overall: G -> OVERALL = sum(GENEs); B -> OVERALL = sum(BLOCKs).  In sum(GENEs) mode, this will ignore BLOCKs for the purposes of calculating OVERALL.  The second character controls summing for each GENE: T -> GENE = Whole gene, ignoring BLOCKs; B -> GENE = sum(BLOCKs).  |
| cm_filters | "near_indel:clustered" | Select which filters to apply during frequency calculation. These filters will also be applied during muts_per_cycle calculation. |  
| runSSCS          | false                  | true or false; whether to do full analysis for SSCS data.  |  
| rerun_type       | 0 (Required for rerun) | What type of rerun you want to do.  0 -> no rerun;  1 -> rerun variant caller;  2 -> rerun postBlastRecovery; 3 -> rerun BLAST and alignment;  4 -> rerun consensus maker.  |  

Save the file as a .csv file with unix line endings (LF).

## 8: Recovery script creation
As part of its operation, this pipeline filters out correct-species reads 
where the blast mapping and bwa mapping positions disagree or where blast 
is unable to determine conclusively where the read maps (i.e.E-scores are 
the same). There is a step which provides an option to recover those reads 
by using a user-generated bash script.  Currently, we use a bash script to 
call a python script which will actually accomplish the recovery, but this 
functionality may be changed in the future.  In general, these scripts must:

1. accept ambiguous reads as $1
2. accept non-ambiguous reads as $2
3. accept reads labeled as incorrect species as $3
4. take a name base for output files as $4
5. take a basePath for location of script files as $5

The script must create the following output files: 

 * ${4}.recovered.temp.bam
 * ${4}.ambig.bam
 * ${4}.wrongSpecies.bam

All script files must be stored in scripts/RecoveryScripts.  

Any packages required by your recovery script can be added to the envs/DS_env_recovery.yaml file.  By default, this contains the conda packages (from bioconda or conda-forge):  

 * Python 3.6
 * samtools
 * pysam
 * snakemake
 * pandas
 * regex
 * biopython
 * pandoc

Note that if you decide to write your own recovery scripts, you are 
responsible for ensuring that your recovery scripts actually work as 
intended, and do not break your computer.  We are unable to guarantee 
that any sample run with a recovery script other than one of the ones 
provided with the pipeline will produce accurate data.  

## 9: Running pipelines:

The pipeline can be run using the DS command created by the setup script:

```bash
/path/to/Duplex-Seq-Pipeline/DS CONFIG_CSV.csv
```

## 10: Output file descriptions:

The pipeline will create a set of summary files covering all samples, as 
well as a file directory structure for each sample.  The summary files are:

| File Name | Description |  
| --------- | ----------- | 
| summary.csv | A csv file with summary metrics for all samples. |  
| summaryDepth.pdf | A pdf file containing depth per target plots for all samples.  |  
| summaryFamilySize.pdf | A pdf file containing family size plots for all samples.  | 
| summaryInsertSize.pdf | A pdf file containing insert size plots for all samples.  |
| summaryMutsByCycle.pdf | A pdf file containing non-SNP mutations per cycle for all samples.  |  


This directory structure looks like this:

```bash
.
└──SAMP_DIR
    ├── Final
    │   ├── dcs
    │   │   ├── FilteredReads
    │   │   │   ├── SAMPLE_dcs.postRecovery.ambig.bam
    │   │   │   ├── SAMPLE_dcs.postRecovery.ambig.bam.bai
    │   │   │   ├── SAMPLE_dcs.postRecovery.wrongSpecies.bam
    │   │   │   └── SAMPLE_dcs.postRecovery.wrongSpecies.bam.bai
    │   │   ├── SAMPLE.dcs.countmuts.csv
    │   │   ├── SAMPLE.dcs.final.bam
    │   │   ├── SAMPLE.dcs.final.bam.bai
    │   │   ├── SAMPLE.dcs.mutated.bam
    │   │   ├── SAMPLE.dcs.mutated.bam.bai
    │   │   ├── SAMPLE.dcs.snps.vcf
    │   │   └── SAMPLE.dcs.vcf
    │   ├── SAMPLE.report.html
    │   └── sscs
    │       ├── SAMPLE.sscs.countmuts.csv
    │       ├── SAMPLE.sscs.final.bam
    │       ├── SAMPLE.sscs.final.bam.bai
    │       ├── SAMPLE.sscs.mutated.bam
    │       ├── SAMPLE.sscs.mutated.bam.bai
    │       ├── SAMPLE.sscs.snps.vcf
    │       └── SAMPLE.sscs.vcf
    ├── Intermediate
    │   ├── ConsensusMakerOutputs
    │   │   ├── SAMPLE_aln_seq1.fq.gz
    │   │   ├── SAMPLE_aln_seq2.fq.gz
    │   │   ├── SAMPLE_read1_dcs.fq.gz
    │   │   ├── SAMPLE_read1_sscs.fq.gz
    │   │   ├── SAMPLE_read2_dcs.fq.gz
    │   │   └── SAMPLE_read2_sscs.fq.gz
    │   └── postBlast
    │       ├── FilteredReads
    │       │   ├── SAMPLE_dcs.ambig.sort.bam
    │       │   ├── SAMPLE_dcs.ambig.sort.bam.bai
    │       │   ├── SAMPLE_dcs.wrongSpecies.sort.bam
    │       │   └── SAMPLE_dcs.wrongSpecies.sort.bam.bai
    │       ├── SAMPLE_dcs.blast.xml
    │       ├── SAMPLE_dcs.preBlast.mutated.bam
    │       └── SAMPLE_dcs.preBlast.unmutated.bam
    ├── logs
    │   └── Log files
    ├── testSeq1.fastq.gz
    ├── testSeq2.fastq.gz
    └── Stats
        ├── data
        │   ├── SAMPLE_cmStats.txt
        │   ├── SAMPLE.dcs_ambiguity_counts.txt
        │   ├── SAMPLE.dcs.iSize_Metrics.txt
        │   ├── SAMPLE.dcs_MutsPerCycle.dat.csv
        │   ├── SAMPLE.dcs.mutsPerRead.txt
        │   ├── SAMPLE.sscs_MutsPerCycle.dat.csv
        │   ├── SAMPLE.sscs.mutsPerRead.txt
        │   ├── SAMPLE.dcs.depth.txt
        │   ├── SAMPLE.dcs.depth.summary.csv
        │   ├── SAMPLE.sscs.depth.txt
        │   ├── SAMPLE.sscs.depth.summary.csv
        │   ├── SAMPLE_dcs.speciesComp.txt
        │   ├── SAMPLE_mem.dcs.sort.flagstats.txt
        │   ├── SAMPLE_mem.sscs.sort.flagstats.txt
        │   ├── SAMPLE.dcs.endClip.metrics.txt
        │   ├── SAMPLE.dcs.overlapClip.metrics.txt
        │   ├── SAMPLE.sscs.endClip.metrics.txt
        │   ├── SAMPLE.sscs.overlapClip.metrics.txt
        │   ├── SAMPLE_onTargetCount.txt
        │   ├── SAMPLE.tagstats.txt
        │   └── SAMPLE.temp.sort.flagstats.txt
        ├── SAMPLE.report.ipynb
        └── plots
            ├── SAMPLE.dcs.iSize_Histogram.png
            ├── SAMPLE.dcs.mutsPerRead.png
            ├── SAMPLE.dcs.targetCoverage.png
            ├── SAMPLE.dcs_BasePerPosInclNs.png
            ├── SAMPLE.dcs_BasePerPosWithoutNs.png
            ├── SAMPLE.sscs.mutsPerRead.png
            ├── SAMPLE.sscs_BasePerPosInclNs.png
            ├── SAMPLE.sscs_BasePerPosWithoutNs.png
            ├── SAMPLE_fam_size_relation.png
            └── SAMPLE_family_size.png
```

File descriptions are as follows:  

| Directory | File name | Description | When Generated |  
| --------- | ------------ | -------------- | -------- |  
| SAMP_DIR | SAMPLE_seq1.fastq.gz | Input read 1 file | Input |  
| SAMP_DIR | SAMPLE_seq2.fastq.gz | Input read 2 file | Input |  
| SAMP_DIR | Final | Directory containing final bam and vcf files | Always |  
| SAMP_DIR/Final | dcs | Directory containing final dcs files | Always |  
| SAMP_DIR/Final/dcs | FilteredReads | Directory containing reads still filtered after postBlastRecovery | blast_db!=NONE |  
| SAMP_DIR/Final/dcs/FilteredReads | SAMPLE.dcs.postRecovery.ambig.bam | File containing reads still considered ambiguous after postBlastRecovery.  May be empty.  | blast_db!=NONE |  
| SAMP_DIR/Final/dcs/FilteredReads | SAMPLE.dcs.postRecovery.ambig.bam.bai | Index for SAMPLE.dcs.postRecovery.ambig.bam | blast_db!=NONE |  
| SAMP_DIR/Final/dcs/FilteredReads | SAMPLE.dcs.postRecovery.wrongSpecies.bam | File containing reads still considered ambiguous after postBlastRecovery. May be empty. | blast_db!=NONE |  
| SAMP_DIR/Final/dcs/FilteredReads | SAMPLE.dcs.postRecovery.wrongSpecies.bam.bai | Index for SAMPLE.dcs.postRecovery.wrongSpecies.bam | blast_db!=NONE |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.countmuts.csv |  Countmuts file, showing a summary of mutation data for DCS reads.   | Always |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.final.bam | Final file for DCS reads, including all reads that overlap the bed file.   | Always |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.final.bam.bai | Index for final DCS reads.   | Always |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.mutated.bam | File containing DCS reads with non-SNP mutations | Always |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.mutated.bam.bai | Index for DCS mutated reads file.   | Always |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.snps.vcf | VCF file containing SNPs overlapping bed file in DCS. | Always |  
| SAMP_DIR/Final/dcs | SAMPLE.dcs.vcf | VCF file containing all variants overlapping bed file in DCS. | Always |  
| SAMP_DIR/Final | SAMPLE.report.html | Summary report for this sample | Always |  
| SAMP_DIR/Final | sscs | Directory containing final SSCS files | Always |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.countmuts.csv |  Countmuts file, showing a summary of mutation data for SSCS reads.   | runSscs=True |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.final.bam | Final file for SSCS reads, including all reads that overlap the bed file.   | Always |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.final.bam.bai | Index for final SSCS reads.   | Always |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.mutated.bam | File containing SSCS reads with non-SNP mutations | runSscs=True |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.mutated.bam.bai | Index for SSCS mutated reads file.   | runSscs=True |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.snps.vcf | VCF file containing SNPs overlapping bed file in SSCS. | runSscs=True |  
| SAMP_DIR/Final/sscs | SAMPLE.sscs.vcf | VCF file containing all variants overlapping bed file in SSCS. | runSscs=True |  
| SAMP_DIR | Intermediate | Directory containing intermediate checkpointing files | Always |  
| SAMP_DIR/Intermediate | ConsensusMakerOutputs | Directory for post-consensus maker checkpoint files | Always |  
| SAMP_DIR/Intermediate/ConsensusMakerOutputs | SAMPLE_aln_seq1.fq.gz | Read 1 file for raw on-target determination | Always |  
| SAMP_DIR/Intermediate/ConsensusMakerOutputs | SAMPLE_aln_seq2.fq.gz | Read 2 file for raw on-target determination | Always |  
| SAMP_DIR/Intermediate/ConsensusMakerOutputs | SAMPLE_read1_dcs.fq.gz |  DCS Read 1 File  | Always |  
| SAMP_DIR/Intermediate/ConsensusMakerOutputs | SAMPLE_read1_sscs.fq.gz |  SSCS Read 1 file  | Always |  
| SAMP_DIR/Intermediate/ConsensusMakerOutputs | SAMPLE_read2_dcs.fq.gz |  DCS Read 2 File  | Always |  
| SAMP_DIR/Intermediate/ConsensusMakerOutputs | SAMPLE_read2_sscs.fq.gz |  SSCS Read 2 file  | Always |  
| SAMP_DIR/Intermediate | postBlast | Directory for post-BLAST checkpoint files.  Only affects DCS.   | Always |  
| SAMP_DIR/Intermediate/postBlast | FilteredReads | Directory for reads that got filtered out of DCS processing due to BLAST analysis indicating that they were either the wrong species or ambiguously mapped.   | Always |  
| SAMP_DIR/Intermediate/postBlast/FilteredReads | SAMPLE_dcs.ambig.sort.bam | Reads that were filtered out due to ambiguous mapping according to BLAST alignment.   | Always |  
| SAMP_DIR/Intermediate/postBlast/FilteredReads | SAMPLE_dcs.ambig.sort.bam.bai | Index for ambiguous reads file | Always |  
| SAMP_DIR/Intermediate/postBlast/FilteredReads | SAMPLE_dcs.wrongSpecies.sort.bam | Reads that were filtered out due to BLAST alignment indicating that they were from the wrong species, or where species of origin could not be determined.   | Always |  
| SAMP_DIR/Intermediate/postBlast/FilteredReads | SAMPLE_dcs.wrongSpecies.sort.bam.bai | Index for wrong-species file | Always |  
| SAMP_DIR/Intermediate/postBlast | SAMPLE_dcs.blast.xml | BLAST xml output | Always |  
| SAMP_DIR/Intermediate/postBlast | SAMPLE_dcs.preBlast.mutated.bam | DCS with potential non-SNP variants that were submitted to BLAST. | Always |  
| SAMP_DIR/Intermediate/postBlast | SAMPLE_dcs.preBlast.unmutated.bam | DCS reads without non-SNP variants.   | Always |  
| SAMP_DIR | logs | Directory containing log files for this sample.   | Always |  
| SAMP_DIR | Stats | Directory containing statistics files | Always |  
| SAMP_DIR/Stats | data | Directory containing statistics data files.   | Always |  
| SAMP_DIR/Stats/data | SAMPLE_cmStats.txt | Statistics from the Consensus Maker | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs_ambiguity_counts.txt | Statistics on ambiguity counts | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs.iSize_Metrics.txt | Statistics on insert size in DCS | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs_MutsPerCycle.dat.csv | Statistics file for non-SNP mutations per cycle in DCS reads | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs.mutsPerRead.txt | Statistics file for non-SNP mutations per read in DCS reads | Always |  
| SAMP_DIR/Stats/data | SAMPLE.sscs_MutsPerCycle.dat.csv |  Text data of error rate per cycle in unclipped SSCS  | runSscs=True |  
| SAMP_DIR/Stats/data | SAMPLE.sscs.mutsPerRead.txt | Statistics file for non-SNP mutations per read in SSCS reads | runSscs=True |  
| SAMP_DIR/Stats/data | SAMPLE.dcs.depth.txt | Per-base coverage and N counts for final DCS  | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs.depth.summary.csv | Per-bed region min, mean, median, and max non-zero depth for final DCS.  | Always |  
| SAMP_DIR/Stats/data | SAMPLE.sscs.depth.txt | Per-base coverage and N counts for final SSCS  | runSscs=True |  
| SAMP_DIR/Stats/data | SAMPLE.sscs.depth.summary.csv | Per-bed region min, mean, median, and max non-zero depth for final SSCS.  | runSscs=True |  
| SAMP_DIR/Stats/data | SAMPLE_dcs.speciesComp.txt | File containing species assignment data for DCS reads | Always |  
| SAMP_DIR/Stats/data | SAMPLE_mem.dcs.sort.flagstats.txt | Initial alignment statistics for DCS reads | Always |  
| SAMP_DIR/Stats/data | SAMPLE_mem.sscs.sort.flagstats.txt | Initial alignment statistics for SSCS reads | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs.endClip.metrics.txt | Statistics on fixed end clipping in DCS | Always |  
| SAMP_DIR/Stats/data | SAMPLE.dcs.overlapClip.metrics.txt | Statistics on overlap clipping in DCS | Always |  
| SAMP_DIR/Stats/data | SAMPLE.sscs.endClip.metrics.txt | Statistics on fixed end clipping in SSCS | runSscs=True |  
| SAMP_DIR/Stats/data | SAMPLE.sscs.overlapClip.metrics.txt | Statistics on overlap clipping in SSCS | runSscs=True |  
| SAMP_DIR/Stats/data | SAMPLE_onTargetCount.txt | Raw on target statistics | Always |  
| SAMP_DIR/Stats/data | SAMPLE.tagstats.txt |  Family size data (in text form)  | Always |  
| SAMP_DIR/Stats/data | SAMPLE.temp.sort.flagstats.txt | Statistics on initial read counts | Always |  
| SAMP_DIR/Stats | SAMPLE.report.ipynb | iPython notebook for the HTML report | Always |  
| SAMP_DIR/Stats | plots | Directory containing statistics plots. | Always |  
| SAMP_DIR/Stats/plots | SAMPLE.dcs.iSize_Histogram.png |  Histogram of insert size metrics for un-clipped DCS  | Always |  
| SAMP_DIR/Stats/plots | SAMPLE.dcs.mutsPerRead.png | Plot of mutations per read in DCS | Always |  
| SAMP_DIR/Stats/plots | SAMPLE.dcs.targetCoverage.png | Plot of per-target coverage in DCS | Always |  
| SAMP_DIR/Stats/plots | SAMPLE.dcs_BasePerPosInclNs.png |  Plot of error rate per cycle for DCS, including Ns | Always |  
| SAMP_DIR/Stats/plots | SAMPLE.dcs_BasePerPosWithoutNs.png |  Plot of error rate per cycle for DCS  | Always |  
| SAMP_DIR/Stats/plots | SAMPLE.sscs.mutsPerRead.png | Plot of mutations per read in SSCS | runSscs=True |  
| SAMP_DIR/Stats/plots | SAMPLE.sscs_BasePerPosInclNs.png |  Plot of error rate per cycle SSCS, including Ns | runSscs=True |  
| SAMP_DIR/Stats/plots | SAMPLE.sscs_BasePerPosWithoutNs.png |  Plot of error rate per cycle for SSCS  | runSscs=True |  
| SAMP_DIR/Stats/plots | SAMPLE_fam_size_relation.png |  Plot of relationship between a:b and b:a families  | Always |  
| SAMP_DIR/Stats/plots | SAMPLE_family_size.png | Plot of family size distribution | Always |  

## 11: Extra BAM tags:

| Tag | Type | Meaning |  
| --- | ---- | ------- |  
| XF | String | In SSCS, represents family size.  In DCS, colon-deliniated list of family sizes for ab1:ba2 or ba1:ab2. |  
| t0 | Integer | With BLAST, TaxID of the species the read matched most closely. Note that negative numbers have special meanings; see below. |  
| t# | Integer | With BLAST, TaxID of the #th BLAST hit. |  
| c# | String | With BLAST, the chromosome of the #th BLAST hit. |  
| p# | Integer | With BLAST, the position of the #th BLAST hit. |  
| l# | Integer | With BLAST, the length of the #th BLAST hit. |  
| YB | String | With BLAST, True if a read was BLASTed. |  
| am | integer | With BLAST, the ambiguity code for the read. See below for meanings. |  

Some negative values for tag t0 have special meanings: 

| t0 value | meaning |  
| -------- | ------- |  
| -1 | Between species BLAST tie. |  
| -3 | No BLAST results. |  
| -4 | Read was submitted to BLAST, but BLAST did not attempt alignment. |  

Ambiguity codes have manings: 

| am value | meaning |  
| 0 | BLAST has 1 best match, matches bwa position |  
| 1 | BLAST has 1 best match, does not match bwa position |  
| 2 | BLAST has 2+ best matches, correct species |  
| 3 | BLAST has 2+ best matches, at least one incorrect species |  
| 4 | BLAST did not attempt alignment, or no blast matches |  
| 5 | Not BLASTed |  


## 12: Testing the pipeline

The newly setup pipeline can be tested using provided data and files 
located in the 'test' directory. To test the pipeline, change into the 
'test' directory and invoking at the command prompt: 

```bash
DS testConfig.csv
```
A final set of output reports can be found in the testData/Final directory and 
be compared to the reports in the ExpectedReports directory located in 
the parent test directory.

## 13: Full and partial reruns

Sometimes it may be necessary to rerun all or part of the pipeline for 
various reasons.  In order to facilitate this, we have provided a script 
(DS-clean) which will prepare samples to rerun based on the "rerun_type" 
column in the config file.  

The following table lists some of the reasons you might want to rerun all 
or part of the pipeline, and how much of the pipeline you want to rerun 
in those cases.  

| Issue | Amount to rerun | rerun_type |  
| ----- | --------------- | ----------------- |  
| <ul><li>Wrong target bed file used</li><li>Wrong clipping parameters used</li></ul> | From pre-variant calling | 1 |  
| Wrong target taxon ID used | From post-blast | 2 |  
| <ul><li>Wrong contaminant db used</li><li>Wrong reference genome used</li><li>Wrong adapter sequence used</li></ul> | From post-Consensus Maker | 3 |  
| Wrong consensus making parameters used | From beginning | 4 |  

To finish preparing for and executing a rerun, run:

```bash
/path/to/Duplex-Seq-Pipeline/DS-clean CONFIG_CSV.csv
/path/to/Duplex-Seq-Pipeline/DS CONFIG_CSV.csv
```

## 14: Unlocking following a power failure
In the event that pipeline execution is interrupted (such as by a power failure), 
the directory can be unlocked in order to restart using the provided 
DS-unlock script:

```bash 
/path/to/Duplex-Seq-Pipeline/DS-unlock CONFIG_CSV.csv
```
