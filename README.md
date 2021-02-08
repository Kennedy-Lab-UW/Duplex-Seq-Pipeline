Duplex Sequencing Pipeline
======================

Brendan Kohrn, March 30, 2020  
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
10. Output file description
11. Testing the pipeline
12. Full and partial reruns

## 1: Glossary

### Single Stranded Consensus Sequence (SSCS)

A construct created by comparing multiple reads and deciding ambiguities by
simple majority.

### Duplex Consensus Sequence (DCS)

A construct created by comparing two SSCSs.

### Family

A group of reads that shares the same tag sequence.

## 2: Dependencies:
This pipeline is known to work with the following minimum versions of the follow required programs:

* Python3.6+
* Snakemake=5.25.0
* Pandas
* Miniconda/Anaconda=4.7.\*
* A GATK3.8.1 .jar file
* bwa=0.7.17.* (for genome setup)
* ncbi-blast=>2.6.0 (installed separately, for contaminant database setup)
* wget (on macOS, install using homebrew; present by default on linux)

Once Python3.6 is installed, snakemake and pandas can be installed using 
pip3 or using whatever package manager you're using.  GATK3.8.1 can be 
downloaded from 
https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2.  
**We need to use GATK3.8.1 since it is the last version of GATK that 
includes the IndelRealigner functionality; this part was removed in 
GATK4**.  Blast can be downloaded in any of several ways, including some 
package managers (Ubuntu: sudo apt-get install ncbi-blast+).  It can 
also be installed using conda if desired.  

## 3: Setup: 

Find the location where you want the pipeline to be located and clone the pipeline using git:

```bash
git clone https://github.com/KennedyLabUW/Duplex-Seq-Pipeline.git
```

change into the directory, and run:

```bash
bash setupDS.sh GATK_JAR_PATH MAX_CORES
```

where GATK_JAR_PATH is the path to your pre-downloaded GATK jar, and MAX_CORES is the maximum number of cores you want the pipeline to be able to use. After this, with the exception of setting up the Genomes (See Section 4) and the (optional) blast contamination database (See Section 5), you should be able to run the Duplex-Seq pipeline using 

```bash
DS CONFIG_CSV.csv
```

where CONFIG_CSV.csv is a configuration CSV file generated as described below.  

## 4: Genome setup

Put genomes in an easily findable location, such as our references directory (i.e. ~/bioinformatics/reference). 

Many genomes can be downloaded from UCSC (http://hgdownload.soe.ucsc.edu/downloads.html).  In order to download genomes from there, you will need to download the twoBitToFa program from the appropriate Utilities directory.  twoBitToFa has the following syntax (Copied from UCSC):

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

Once you have downloaded a genome and converted it into FASTA format, it needs to be indexed.  To do this, open a terminal and navigate to your the directory containing your genome.  Note that in the following commands, the word "genome.fasta" should be replaced with the file name of your genome. These commands **must** be run in the same directory as the reference genome

```
    bwa index genome.fasta  
    samtools faidx genome.fasta  
    # /Path/To/PicardTools should be the path to wherever you put your picard tools jar file.    
    # Also, genome.dict should match the name of your fasta (e.g. hg19.fasta to hg19.dict)  
    java -jar /Path/To/PicardTools/picard.jar CreateSequenceDictionary R=genome.fasta O=genome.dict  
```

At the moment, this pipeline does not support compressed genomes.

## 5: Contaminant Database Setup:
The Duplex-Seq pipeline is designed to use a local NCBI Blast instance to detect and remove potential contamination from non-target species and identify issues arising from pseudogenes.
*This step is optional, but requires a non valid NCBI Blast database file name. We recommend either a dummy path and file name or '.' to be present in the blast_db field in the config.csv file. Blank entries will result in pipeline failure.*
 
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
names, your genome in the database cannot use NCBI chromosme names, but 
must also use UCSC chromosome names)**.  

Database setup consists of three steps:

 1. Genome labeling:

Label each record in the genome with the ncbi taxID for the species the 
genome is associated with.  TaxIDs can be found using the NCBI taxonomy 
website (https://www.ncbi.nlm.nih.gov/taxonomy).  This is done by running:

```
python3 AddTaxonID.py GENOME.fa TAXID GENOME_taxID.fa
```

where GENOME.fa is the input fasta file with the genome, TAXID is the NCBI taxonomy ID for the species associated with the genome, and GENOME_taxID.fa is the output labeled genome.  

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

## 6: Bed  file preparation

A bed file is a file which details regions of the genome in which you 
are interested.  The syntax for bed files is described at 
https://genome.ucsc.edu/FAQ/FAQformat.html#format1.  

If you know which genes you are targeting and are using a common 
published genome, the bed file can be downloaded from the UCSC Table 
Browser (https://genome.ucsc.edu/cgi-bin/hgTables).  Otherwise, it can 
be created using results from a BLAST search or any other method you 
like.  **This pipeline does not currently support overlapping intervals, 
does support blocks as described in the bed spec.**  

## 7: Configuration file creation:

Use the ConfigTemplate to create a new file with the appropriate headers.  For each row, fill in the information about a particular sample:

| Header           | Required or Default | Information |
| ---------------- | ------------------- | ----------- |
| sample           | Required        | A unique identifier for a sample; this will be used to name all output files for this sample |
| rglb             | Required        | Read Group Library Identifier |
| rgpl             | Required        | Read Group Platform; usually illumina |
| rgpu             | Required        | Read Group Platform Unit |
| rgsm             | Required        | Read Group Sample |
| reference        | Required        | The path to the prepared reference genome to use with this sample.  |
| target_bed       | Required        | A bed file showing where the targets are for this particular sample |  
| blast_db         | Required        | The blast database to use for contaminant filtering; must include your target genome.  |
| targetTaxonId    | Required        | The taxon ID of the species you are expecting to be present in the sample.  |
| baseDir          | Required        | The directory the input files are in, and where the output files will be created. |
| in1              | Required        | The read1 fastq (or fastq.gz, or fq.gz, or fq) file for this sample. Note that this is just the name of the file, and not the full path.  |
| in2              | Required        | The read2 fastq (or fastq.gz, or fq.gz, or fq) file for this sample. Note that this is just the name of the file, and not the full path.   |
| mqFilt           | 0               | A threshold for mapping quality filtering, if desired. |
| minMem           | 0               | The minimum number of reads that must be in a family for consensus making |
| maxMem           | 200             | The maximum number of reads in a family the consensus maker should consider. |
| cutOff           | 0.9             | The threshold for consensus making; the consensus maker will require at least this much agreement on a per base pair level. |
| nCutOff          | 1               | The maximum proportion of N bases in an output consensus sequence. |
| umiLen           | 8               | The length of the UMI in this sample |
| spacerLen        | 1               | The length of the spacer sequence in this sample |
| locLen           | 10              | The localization length to use for this sample |
| readLen          | 101             | The length of a read for this sample |
| clipBegin        | 7               | How many bases to clip off the 5' end of the read |
| clipEnd          | 0               | How many bases to clip off the 3' end of the read |
| minClonal        | 0               | The minimum clonality to use for count_muts generation |
| maxClonal        | 0.1             | The maximum clonality to use for count_muts generation |
| minDepth         | 100             | The minimum depth to use for count_muts generation |
| maxNs            | 1               | The maximum proportion of N bases to use for count_muts generation |
| cm_outputs       | "GB"            | Select which sections of the countmuts to output, in addition to 'OVERALL'.  String of one or more of 'G', 'B', and 'N'.  G -> output GENE sections for each bed line; B -> output 'BLOCK' sections for each block in the bed line (if present); 'N' -> Only output overall frequencies.  Overrides all other options. |  
| cm_sumTypes      | "GT"            | How to calculate OVERALL and GENE blocks for countmuts output. The first character controls summing for overall: G -> OVERALL = sum(GENEs); B -> OVERALL = sum(BLOCKs).  In sum(GENEs) mode, this will ignore BLOCKs for the purposes of calculating OVERALL.  The second character controls summing for each GENE: T -> GENE = Whole gene, ignoring BLOCKs; B -> GENE = sum(BLOCKs).  |
| runSSCS          | false           | true or false; whether to do full analysis for SSCS data.  | 
| recovery         | "noRecovery_noSynLink.sh" | The recovery script to use in attempting to recover ambiguously mapped reads (as determine by blast alignment vs bwa alignment).  Recovery script creation is discussed in 8; below.  |  

Save the file as a .csv file with unix line endings (LF).

## 8: Recovery script creation
As part of its operation, this pipeline filters out correct-species reads where the blast mapping and bwa mapping positions disagree or where blast is unable to determine conclusively where the read maps (i.e.E-scores are the same). There is a step which provides an option to recover those reads by using a user-generated bash script.  Currently, we use a bash script to call a python script which will actually accomplish the recovery, but this functionality may be changed in the future.  In general, these scripts must:

1. accept ambiguous reads as $1
2. accept non-ambiguous reads as $2
3. take a file name for output reads as $3
4. take a basePath for location of script files as $4

All script files must be stored in scripts/RecoveryScripts.  

## 9: Running pipelines:

The pipeline can be run using the DS command created by the setup script:

```bash
DS CONFIG_CSV.csv
```

## 10: Output file descriptions:

The pipeline will create a set of summary files covering all samples, as well as a file directory structure for each sample.  The summary files are:

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
└── SAMP_DIR
   ├── Final
   │   ├── dcs
   │   │   ├── SAMPLE.dcs.countmuts.txt
   │   │   ├── SAMPLE.dcs.final.bam
   │   │   ├── SAMPLE.dcs.final.bam.bai
   │   │   ├── SAMPLE.dcs.mutated.bam
   │   │   ├── SAMPLE.dcs.mutated.bam.bai
   │   │   ├── SAMPLE.dcs.snps.vcf
   │   │   └── SAMPLE.dcs.vcf
   │   ├── SAMPLE.report.html
   │   └── sscs
   │       ├── SAMPLE.sscs.final.bam
   │       └── SAMPLE.sscs.final.bam.bai
   ├── Intermediate
   │   ├── ConsensusMakerOutputs
   │   │   ├── SAMPLE_aln_seq1.fq.gz
   │   │   ├── SAMPLE_aln_seq2.fq.gz
   │   │   ├── SAMPLE_read1_dcs.fq.gz
   │   │   ├── SAMPLE_read1_sscs.fq.gz
   │   │   ├── SAMPLE_read2_dcs.fq.gz
   │   │   └── SAMPLE_read2_sscs.fq.gz
   │   ├── postBlast
   │   │   ├── FilteredReads
   │   │   │   ├── SAMPLE_dcs.ambig.sort.bam
   │   │   │   ├── SAMPLE_dcs.ambig.sort.bam.bai
   │   │   │   ├── SAMPLE_dcs.wrongSpecies.sort.bam
   │   │   │   └── SAMPLE_dcs.wrongSpecies.sort.bam.bai
   │   │   ├── SAMPLE_dcs.blast.xml
   │   │   ├── SAMPLE_dcs.preBlast.mutated.bam
   │   │   └── SAMPLE_dcs.preBlast.unmutated.bam
   │   └── PreVariantCallsCp
   │       ├── dcs
   │       │   ├── SAMPLE.dcs.clipped.bai
   │       │   └── SAMPLE.dcs.clipped.bam
   │       └── sscs
   │           ├── SAMPLE.sscs.clipped.bai
   │           └── SAMPLE.sscs.clipped.bam
   ├── logs
   │   └── Log Files
   ├── SAMPLE_config.sh
   ├── SAMPLE_seq1.fastq.gz
   ├── SAMPLE_seq2.fastq.gz
   └── Stats
       ├── data
       │   ├── SAMPLE_cmStats.txt
       │   ├── SAMPLE.dcs_ambiguity_counts.txt
       │   ├── SAMPLE.dcs.clipped.metrics.txt
       │   ├── SAMPLE.dcs.filt.clipped.metrics.txt
       │   ├── SAMPLE.dcs.filt.no_overlap.metrics.txt
       │   ├── SAMPLE.dcs.iSize_Metrics.txt
       │   ├── SAMPLE.dcs_MutsPerCycle.dat.csv
       │   ├── SAMPLE.dcs.mutsPerRead.txt
       │   ├── SAMPLE.dcs.region.mutpos.vcf_depth.txt
       │   ├── SAMPLE.dcs.region.snpFiltered.mutsPerRead.txt
       │   ├── SAMPLE_dcs.speciesComp.txt
       │   ├── SAMPLE_mem.dcs.sort.flagstats.txt
       │   ├── SAMPLE_mem.sscs.sort.flagstats.txt
       │   ├── SAMPLE_onTargetCount.txt
       │   ├── SAMPLE.sscs.clipped.metrics.txt
       │   ├── SAMPLE.sscs.filt.clipped.metrics.txt
       │   ├── SAMPLE.sscs.filt.no_overlap.metrics.txt
       │   ├── SAMPLE.tagstats.txt
       │   └── SAMPLE.temp.sort.flagstats.txt
       ├── SAMPLE.report.ipynb
       └── plots
           ├── SAMPLE.dcs_BasePerPosInclNs.png
           ├── SAMPLE.dcs_BasePerPosWithoutNs.png
           ├── SAMPLE.dcs.iSize_Histogram.png
           ├── SAMPLE.dcs.mutsPerRead.png
           ├── SAMPLE.dcs.region.snpFiltered.mutsPerRead.png
           ├── SAMPLE.dcs.targetCoverage.png
           ├── SAMPLE_family_size.png
           └── SAMPLE_fam_size_relation.png
```

File descriptions are as follows:  

| Directory | File name | Description | When Generated |  
| --------- | ------------ | -------------- | -------- |  
| . | SAMPLE_config.sh | Per-sample config file, showing run parameters for this sample.   | Always |  
| . | SAMPLE_seq1.fastq.gz | Input read 1 file | Input |  
| . | SAMPLE_seq2.fastq.gz | Input read 2 file | Input |  
| . | Final | Directory containing final bam and vcf files | Always |  
| Final | dcs | Directory containing final dcs files | Always |  
| Final/dcs | SAMPLE.dcs.countmuts.txt |  Countmuts file, showing a summary of mutation data for DCS reads.   | Always |  
| Final/dcs | SAMPLE.dcs.final.bam | Final file for DCS reads, including all reads that overlap the bed file.   | Always |  
| Final/dcs | SAMPLE.dcs.final.bam.bai | Index for final DCS reads.   | Always |  
| Final/dcs | SAMPLE.dcs.mutated.bam | File containing DCS reads with non-SNP mutations | Always |  
| Final/dcs | SAMPLE.dcs.mutated.bam.bai | Index for DCS mutated reads file.   | Always |  
| Final/dcs | SAMPLE.dcs.snps.vcf | VCF file containing SNPs overlapping bed file in DCS. | Always |  
| Final/dcs | SAMPLE.dcs.vcf | VCF file containing all variants overlapping bed file in DCS. | Always |  
| Final | SAMPLE.report.html | Summary report for this sample | Always |  
| Final | sscs | Directory containing final SSCS files | Always |  
|  | SAMPLE.sscs.countmuts.txt |  Countmuts file, showing a summary of mutation data for SSCS reads.   |  |  
| Final/sscs | SAMPLE.sscs.final.bam | Final file for SSCS reads, including all reads that overlap the bed file.   | Always |  
| Final/sscs | SAMPLE.sscs.final.bam.bai | Index for final SSCS reads.   | Always |  
| Final/sscs | SAMPLE.sscs.mutated.bam | File containing SSCS reads with non-SNP mutations | runSscs=True |  
| Final/sscs | SAMPLE.sscs.mutated.bam.bai | Index for SSCS mutated reads file.   | runSscs=True |  
| Final/sscs | SAMPLE.sscs.snps.vcf | VCF file containing SNPs overlapping bed file in SSCS. | runSscs=True |  
| Final/sscs | SAMPLE.sscs.vcf | VCF file containing all variants overlapping bed file in SSCS. | runSscs=True |  
| . | Intermediate | Directory containing intermediate checkpointing files | Always |  
| Intermediate | ConsensusMakerOutputs | Directory for post-consensus maker checkpoint files | Always |  
| Intermediate/ConsensusMakerOutputs | SAMPLE_aln_seq1.fq.gz | Read 1 file for raw on-target determination | Always |  
| Intermediate/ConsensusMakerOutputs | SAMPLE_aln_seq2.fq.gz | Read 2 file for raw on-target determination | Always |  
| Intermediate/ConsensusMakerOutputs | SAMPLE_read1_dcs.fq.gz |  DCS Read 1 File  | Always |  
| Intermediate/ConsensusMakerOutputs | SAMPLE_read1_sscs.fq.gz |  SSCS Read 1 file  | Always |  
| Intermediate/ConsensusMakerOutputs | SAMPLE_read2_dcs.fq.gz |  DCS Read 2 File  | Always |  
| Intermediate/ConsensusMakerOutputs | SAMPLE_read2_sscs.fq.gz |  SSCS Read 2 file  | Always |  
| Intermediate | postBlast | Directory for post-BLAST checkpoint files.  Only affects DCS.   | Always |  
| Intermediate/postBlast | FilteredReads | Directory for reads that got filtered out of DCS processing due to BLAST analysis indicating that they were either the wrong species or ambiguously mapped.   | Always |  
| Intermediate/postBlast/FilteredReads | SAMPLE_dcs.ambig.sort.bam | Reads that were filtered out due to ambiguous mapping according to BLAST alignment.   | Always |  
| Intermediate/postBlast/FilteredReads | SAMPLE_dcs.ambig.sort.bam.bai | Index for ambiguous reads file | Always |  
| Intermediate/postBlast/FilteredReads | SAMPLE_dcs.wrongSpecies.sort.bam | Reads that were filtered out due to BLAST alignment indicating that they were from the wrong species, or where species of origin could not be determined.   | Always |  
| Intermediate/postBlast/FilteredReads | SAMPLE_dcs.wrongSpecies.sort.bam.bai | Index for wrong-species file | Always |  
| Intermediate/postBlast | SAMPLE_dcs.blast.xml | BLAST xml output | Always |  
| Intermediate/postBlast | SAMPLE_dcs.preBlast.mutated.bam | DCS with potential non-SNP variants that were submitted to BLAST. | Always |  
| Intermediate/postBlast | SAMPLE_dcs.preBlast.unmutated.bam | DCS reads without non-SNP variants.   | Always |  
| Intermediate | PreVariantCallsCp | Directory for final checkpoint, immediately before variant calling.   | Always |  
| Intermediate/PreVariantCallsCp | dcs | Directory containing final checkpoint DCS bam files | Always |  
| Intermediate/PreVariantCallsCp/dcs | SAMPLE.dcs.clipped.bai | Final file containing all DCS reads, including those that don’t overlap the bed file.   | Always |  
| Intermediate/PreVariantCallsCp/dcs | SAMPLE.dcs.clipped.bam | Final file containing all DCS reads, including those that don’t overlap the bed file.   | Always |  
| Intermediate/PreVariantCallsCp | sscs | Directory containing final checkpoint SSCS bam files | Always |  
| Intermediate/PreVariantCallsCp/sscs | SAMPLE.sscs.clipped.bai | Index for final file containing all SSCS reads, including those that don’t overlap the bed file.   | Always |  
| Intermediate/PreVariantCallsCp/sscs | SAMPLE.sscs.clipped.bam | Final file containing all SSCS reads, including those that don’t overlap the bed file.   | Always |  
| . | logs | Directory containing log files for this sample.   | Always |  
| . | Stats | Directory containing statistics files | Always |  
| Stats | data | Directory containing statistics data files.   | Always |  
| Stats/data | SAMPLE_cmStats.txt | Statistics from the Consensus Maker | Always |  
| Stats/data | SAMPLE.dcs_ambiguity_counts.txt | Statistics on ambiguity counts | Always |  
| Stats/data | SAMPLE.dcs.clipped.metrics.txt | Statistics on fixed end clipping in DCS | Always |  
| Stats/data | SAMPLE.dcs.filt.clipped.metrics.txt | Statistics on overlap clipping in DCS | Always |  
| Stats/data | SAMPLE.dcs.iSize_Metrics.txt | Statistics on insert size in DCS | Always |  
| Stats/data | SAMPLE.dcs_MutsPerCycle.dat.csv | Statistics file for non-SNP mutations per cycle in DCS reads | Always |  
| Stats/data | SAMPLE.dcs.mutsPerRead.txt | Statistics file for non-SNP mutations per read in DCS reads | Always |  
| Stats/data | SAMPLE.dcs.region.mutpos.vcf_depth.txt | Per-base coverage and N counts for final DCS  | Always |  
| Stats/data | SAMPLE_dcs.speciesComp.txt | File containing species assignment data for DCS reads | Always |  
| Stats/data | SAMPLE_mem.dcs.sort.flagstats.txt | Initial alignment statistics for DCS reads | Always |  
| Stats/data | SAMPLE_mem.sscs.sort.flagstats.txt | Initial alignment statistics for SSCS reads | Always |  
| Stats/data | SAMPLE_onTargetCount.txt | Raw on target statistics | Always |  
| Stats/data | SAMPLE.sscs.clipped.metrics.txt | Statistics on fixed end clipping in SSCS | Always |  
| Stats/data | SAMPLE.sscs.filt.clipped.metrics.txt | Statistics on overlap clipping in SSCS | Always |  
| Stats/data | SAMPLE.sscs_MutsPerCycle.dat.csv |  Text data of error rate per cycle in unclipped SSCS  | runSscs=True |  
| Stats/data | SAMPLE.sscs.mutsPerRead.txt | Statistics file for non-SNP mutations per read in SSCS reads | runSscs=True |  
| Stats/data | SAMPLE.sscs.region.mutpos.vcf_depth.txt | Per-base coverage and N counts for final SSCS  | runSscs=True |  
| Stats/data | SAMPLE.tagstats.txt |  Family size data (in text form)  | Always |  
| Stats/data | SAMPLE.temp.sort.flagstats.txt | Statistics on initial read counts | Always |  
| Stats | SAMPLE.report.ipynb | iPython notebook for the HTML report | Always |  
| Stats | plots | Directory containing statistics plots. | Always |  
| Stats/plots | SAMPLE.dcs_BasePerPosInclNs.png |  Plot of error rate per cycle for DCS, including Ns | Always |  
| Stats/plots | SAMPLE.dcs_BasePerPosWithoutNs.png |  Plot of error rate per cycle for DCS  | Always |  
| Stats/plots | SAMPLE.dcs.iSize_Histogram.png |  Histogram of insert size metrics for un-clipped DCS  | Always |  
| Stats/plots | SAMPLE.dcs.mutsPerRead.png | Plot of mutations per read in DCS | Always |  
| Stats/plots | SAMPLE.dcs.targetCoverage.png | Plot of per-target coverage in DCS | Always |  
| Stats/plots | SAMPLE_family_size.png | Plot of family size distribution | Always |  
| Stats/plots | SAMPLE_fam_size_relation.png |  Plot of relationship between a:b and b:a families  | Always |  
| Stats/plots | SAMPLE.sscs_BasePerPosInclNs.png |  Plot of error rate per cycle SSCS, including Ns | runSscs=True |  
| Stats/plots | SAMPLE.sscs_BasePerPosWithoutNs.png |  Plot of error rate per cycle for SSCS  | runSscs=True |  
| Stats/plots | SAMPLE.sscs.mutsPerRead.png | Plot of mutations per read in SSCS | runSscs=True |  

## 11: Testing the pipeline

The newly setup pipeline can be tested using provided data and files located in the 'test' directory. 
To test the pipeline, change into the 'test' directory and invoking at the command prompt: 

```bash
DS testConfig.csv
```
A final expected output report can be found in the testData/Final directory and be compared to the expected_report.html file
located in the parent test directory.

## 12: Full and partial reruns

Sometimes it may be necessary to rerun all or part of the pipeline for various reasons.  The following table lists some of the reasons you might want to rerun all or part of the pipeline, how much of the pipeline you want to rerun in those cases, and how to carry out the rerun.  

| Issue | Amount to rerun | Preparation steps |  
| ----- | --------------- | ----------------- |  
| Wrong bed file used | From pre-variant calling | Remove the following files: <ul><li>summary.csv</li><li>summaryMutsByCycle.pdf</li><li>summaryDepth.pdf</li><li>summaryFamilySize.pdf</li><li>summaryInsertSize.pdf</li><li>All _config.sh files for samples that you want to force a rerun of</li></ul>Remove the following directories: <ul><li>All Final directories in samples that you want to force a rerun of.  </li></ul> |
| <ul><li>Wrong clipping parameters used</li><li>Wrong target taxon ID used</li></ul> | From post-blast | Remove the following files: <ul><li>summary.csv</li><li>summaryMutsByCycle.pdf</li><li>summaryDepth.pdf</li><li>summaryFamilySize.pdf</li><li>summaryInsertSize.pdf</li><li>All _config.sh files for samples that you want to force a rerun of</li></ul>Remove the following directories: <ul><li>All Final directories in samples that you want to force a rerun of.  </li><li>All Intermediate/PreVariantCallsCp directories in samples that you want to force a rerun of</li></ul> |
| <ul><li>Wrong contaminant db used</li><li>Wrong reference genome used</li></ul> | From post-Consensus Maker | Remove the following files: <ul><li>summary.csv</li><li>summaryMutsByCycle.pdf</li><li>summaryDepth.pdf</li><li>summaryFamilySize.pdf</li><li>summaryInsertSize.pdf</li><li>All _config.sh files for samples that you want to force a rerun of</li></ul>Remove the following directories: <ul><li>All Final directories in samples that you want to force a rerun of.  </li><li>All Intermediate/PreVariantCallsCp directories in samples that you want to force a rerun of</li><li>All Intermediate/postBlast directories in samples that you want to force a rerun of</li></ul> |
| Wrong consensus making parameters used | From beginning | Remove the following files: <ul><li>summary.csv</li><li>summaryMutsByCycle.pdf</li><li>summaryDepth.pdf</li><li>summaryFamilySize.pdf</li><li>summaryInsertSize.pdf</li><li>All _config.sh files for samples that you want to force a rerun of</li></ul>Remove the following directories: <ul><li>All Final directories in samples that you want to force a rerun of.  </li><li>All Intermediate directories in samples that you want to force a rerun of</li><li>All Stats directories in samples that you want to force a rerun of</li></ul>
