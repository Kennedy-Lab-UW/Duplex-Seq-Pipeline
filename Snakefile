import pandas as pd
from snakemake.utils import validate
import multiprocessing as mp
import sys

configfile: f"{sys.path[0]}/DS_progConfig.yaml"
# ~ validate(config, "config.schema.yaml")

samples = pd.read_csv(config["samples"]).set_index("sample", drop=False)
validate(samples, f"{sys.path[0]}/DS_baseSchema.yaml")


config["maxCores"] = min(config["maxCores"], mp.cpu_count())

# ~ print(samples.index)

#print(config["maxCores"])

def get_sample(wildcards):
    return samples.loc[wildcards.sample, "sample"]
def get_rglb(wildcards):
    return samples.loc[wildcards.sample, "rglb"]
def get_rgpl(wildcards):
    return samples.loc[wildcards.sample, "rgpl"]
def get_rgpu(wildcards):
    return samples.loc[wildcards.sample, "rgpu"]
def get_rgsm(wildcards):
    return samples.loc[wildcards.sample, "rgsm"]
def get_reference(wildcards):
    return samples.loc[wildcards.sample, "reference"]
def get_target_bed(wildcards):
    return samples.loc[wildcards.sample, "target_bed"]
def get_blast_db(wildcards):
    return samples.loc[wildcards.sample, "blast_db"]
def get_blast_db_path(wildcards):
    if samples.loc[wildcards.sample, "blast_db"].upper() == "NONE":
        return "NONE.nal"
    else:
        return f'{samples.loc[wildcards.sample, "blast_db"]}.nal'
def get_target_taxon(wildcards):
    return samples.loc[wildcards.sample, "targetTaxonId"]
def get_baseDir(wildcards):
    return samples.loc[wildcards.sample, "baseDir"]
def get_in1(wildcards):
    return f'{samples.loc[wildcards.sample, "baseDir"]}/{samples.loc[wildcards.sample, "in1"]}'
def get_in2(wildcards):
    return f'{samples.loc[wildcards.sample, "baseDir"]}/{samples.loc[wildcards.sample, "in2"]}'
def get_mqFilt(wildcards):
    return samples.loc[wildcards.sample, "mqFilt"]
def get_minMem(wildcards):
    return samples.loc[wildcards.sample, "minMem"]
def get_maxMem(wildcards):
    return samples.loc[wildcards.sample, "maxMem"]
def get_cutOff(wildcards):
    return samples.loc[wildcards.sample, "cutOff"]
def get_nCutOff(wildcards):
    return samples.loc[wildcards.sample, "nCutOff"]
def get_umiLen(wildcards):
    return samples.loc[wildcards.sample, "umiLen"]
def get_spacerLen(wildcards):
    return samples.loc[wildcards.sample, "spacerLen"]
def get_locLen(wildcards):
    return samples.loc[wildcards.sample, "locLen"]
def get_readLen(wildcards):
    return samples.loc[wildcards.sample, "readLen"]
def get_clipBegin(wildcards):
    return samples.loc[wildcards.sample, "clipBegin"]
def get_clipEnd(wildcards):
    return samples.loc[wildcards.sample, "clipEnd"]
def get_runSSCS(wildcards):
    return samples.loc[wildcards.sample, "runSSCS"]
def get_runDCS(wildcards):
    return samples.loc[wildcards.sample, "runDCS"]
def get_makeDCS(wildcards):
    return samples.loc[wildcards.sample, "makeDCS"]
def get_cm_outputs(wildcards):
    return samples.loc[wildcards.sample, "cm_outputs"]
def get_cm_sumTypes(wildcards):
    return samples.loc[wildcards.sample, "cm_sumTypes"]
def get_minClonal(wildcards):
    return samples.loc[wildcards.sample, "minClonal"]
def get_maxClonal(wildcards):
    return samples.loc[wildcards.sample, "maxClonal"]
def get_minDepth(wildcards):
    return samples.loc[wildcards.sample, "minDepth"]
def get_maxNs(wildcards):
    return samples.loc[wildcards.sample, "maxNs"]
def get_cleanup(wildcards):
    return samples.loc[wildcards.sample, "cleanup"]
def get_recovery(wildcards):
    return f'{sys.path[0]}/scripts/RecoveryScripts/{samples.loc[wildcards.sample, "recovery"]}'

def getVarDictBam(wildcards):
    if wildcards.sampType == 'dcs':
        if get_blast_db_path(wildcards) == "NONE.nal":
            output=(f"{wildcards.runPath}/"
                f"{wildcards.sample}_mem"
                f".dcs.sort.bam"
                )
        else:
            output=(f"{wildcards.runPath}/"
                f"{wildcards.sample}_dcs"
                f".postRecovery.recovered.temp.bam"
                )
    elif wildcards.sampType == 'sscs':
        output=(f"{wildcards.runPath}/"
                f"{wildcards.sample}_mem"
                f".sscs.sort.bam"
                )
    else:
        raise Exception("Wrong sampType")
    return output

def get_outFiles(prefix="", sampType="dcs", suffix=".clipped.bam"):
    outList = []
    for sampIter in samples.index:
        if sampType == "dcs":
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"{prefix}"
                f"{sampIter}.{sampType}{suffix}"
                )
        elif sampType == "":
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"{prefix}"
                f"{sampIter}{suffix}"
                )
        elif samples.loc[sampIter,'runSSCS']:
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"{prefix}"
                f"{sampIter}.{sampType}{suffix}"
                )
    return outList
    
def get_outCountMuts(sampType="dcs"):
    outList = []
    for sampIter in samples.index:
        if sampType != "sscs":
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"Final/{sampType}/"
                f"{sampIter}.{sampType}.countmuts.csv"
                )
        elif samples.loc[sampIter,'runSSCS']:
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"Final/{sampType}/"
                f"{sampIter}.{sampType}.countmuts.csv"
                )
    return outList

def get_outMems(sampType="dcs"):
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/"
            f"Stats/data/"
            f"{sampIter}_mem.{sampType}.sort.flagstats.txt"
            )
    return outList

def getOutConfig(wildcards):
    return f'{samples.loc[wildcards.sample, "baseDir"]}/{samples.loc[wildcards.sample, "baseDir"]}_config.sh'

def get_outMutPos(sampType="dcs"):
    outList = []
    for sampIter in samples.index:
        if sampType != "sscs":
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"{sampIter}.{sampType}.mutpos"
                )
        elif samples.loc[sampIter,'runSSCS']:
            # ~ print("RunSSCS True")
            outList.append(
                f"{samples.loc[sampIter,'baseDir']}/"
                f"{sampIter}.{sampType}.mutpos"
                )
        # ~ else:
            # ~ print("RunSSCS False")
    return outList

def get_outOnTarget():
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/"
            f"Stats/data/"
            f"{sampIter}_onTargetCount.txt"
            )
    return outList

def getDepthFiles():
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/"
            f"Stats/plots/"
            f"{sampIter}.dcs.targetCoverage.png"
            )
    return outList

def getInsertFiles():
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/"
            f"Stats/plots/"
            f"{sampIter}.dcs.iSize_Histogram.png"
            )
    return outList

def getMutsByCycFiles():
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/"
            f"Stats/data/"
            f"{sampIter}.dcs_MutsPerCycle.dat.csv"
            )
    return outList

def getFamSizeFiles():
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/"
            f"Stats/data/{sampIter}.tagstats.txt"
            )
    return outList
    
def getReportInput(wildcards):
    # BLAST switch:
    if get_blast_db_path(wildcards) == "NONE.nal":
        withBlast = False
    else:
        withBlast = True
    outArgs = []
    # Post-CM checkpoint files
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/ConsensusMakerOutputs/{wildcards.sample}_read1_sscs.fq.gz')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/ConsensusMakerOutputs/{wildcards.sample}_read2_sscs.fq.gz')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/ConsensusMakerOutputs/{wildcards.sample}_read1_dcs.fq.gz')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/ConsensusMakerOutputs/{wildcards.sample}_read2_dcs.fq.gz')
    if withBlast:
        # post-BLAST checkpoint files:
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/{wildcards.sample}_dcs.blast.xml')
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.mutated.bam')
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.unmutated.bam')
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.ambig.sort.bam')
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.wrongSpecies.sort.bam')
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.ambig.sort.bam.bai')
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.wrongSpecies.sort.bam.bai')
        #'BLAST Filtering Stats File
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}_dcs.speciesComp.txt')
        # Ambiguity stats file
        outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}.dcs_ambiguity_counts.txt')
    # Raw read stats files
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}.temp.sort.flagstats.txt')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}_onTargetCount.txt')
    #'Consensus Maker Stats File
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}.tagstats.txt')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}_cmStats.txt')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}_family_size.png')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}_fam_size_relation.png')
    #'SSCS Alignment Stats File
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}_mem.sscs.sort.flagstats.txt')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/data/{wildcards.sample}_mem.dcs.sort.flagstats.txt')
    # Final BAM files
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/sscs/{wildcards.sample}.sscs.final.bam')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/sscs/{wildcards.sample}.sscs.final.bam.bai')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/dcs/{wildcards.sample}.dcs.final.bam')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/dcs/{wildcards.sample}.dcs.final.bam.bai')
    #'Clipping Stats files
    #'Mutations Stats Files
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/dcs/{wildcards.sample}.dcs.countmuts.csv')
    #'Final stats files
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}.dcs.iSize_Histogram.png')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/dcs/{wildcards.sample}.dcs.mutated.bam')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Final/dcs/{wildcards.sample}.dcs.mutated.bam.bai')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}.dcs_BasePerPosInclNs.png')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}.dcs_BasePerPosWithoutNs.png')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}.dcs.mutsPerRead.png')
    outArgs.append(f'{samples.loc[wildcards.sample, "baseDir"]}/Stats/plots/{wildcards.sample}.dcs.targetCoverage.png')
    
    return outArgs

def get_reports():
    outList = []
    for sampIter in samples.index:
        outList.append(
            f"{samples.loc[sampIter,'baseDir']}/Final/"
            f"{sampIter}.report.html"
            )
    return outList

def getSummaryInput():
    outFiles = []
    outFiles.extend(get_outFiles(prefix="Stats/data/", sampType="", suffix=".temp.sort.flagstats.txt"))
    outFiles.extend(get_outFiles(prefix="Stats/data/", sampType="", suffix="_mem.sscs.sort.flagstats.txt"))
    outFiles.extend(get_outFiles(prefix="Stats/data/", sampType="", suffix="_mem.dcs.sort.flagstats.txt"))
    outFiles.extend(get_outFiles(prefix="Stats/data/", sampType="", suffix=".tagstats.txt"))
    outFiles.extend(get_outCountMuts(sampType="sscs"))
    outFiles.extend(get_outCountMuts(sampType="dcs"))
    outFiles.extend(get_outMems(sampType="sscs"))
    outFiles.extend(get_outMems(sampType="dcs"))
    # ~ outFiles.extend(get_outMutPos(sampType="sscs"))
    # ~ outFiles.extend(get_outMutPos(sampType="dcs"))
    outFiles.extend(get_outFiles(prefix="Final/sscs/", sampType="sscs", suffix=".mutated.bam"))
    outFiles.extend(get_outFiles(prefix="Final/dcs/", sampType="dcs", suffix=".mutated.bam"))
    outFiles.extend(get_outFiles(prefix="Final/sscs/", sampType="sscs", suffix=".vcf"))
    outFiles.extend(get_outFiles(prefix="Final/dcs/", sampType="dcs", suffix=".vcf"))
    outFiles.extend(get_outFiles(prefix="Stats/data/", sampType="dcs", suffix=".depth.txt"))
    outFiles.extend(get_reports())
    return outFiles

wildcard_constraints:
    sampType="(sscs)|(dcs)",
    sample="|".join([f"({x})" for x in samples.index]),
    runPath="|".join([f"({samples.loc[x, 'baseDir']})" for x in samples.index])

# Run all steps of the pipeline
rule all:
    input:
        f"{config['samples']}.summary.csv", 
        f"{config['samples']}.summaryDepth.pdf", 
        f"{config['samples']}.summaryInsertSize.pdf", 
        f"{config['samples']}.summaryMutsByCycle.pdf", 
        f"{config['samples']}.summaryFamilySize.pdf"
        # ~ get_outCountMuts("sscs"),
        # ~ get_outCountMuts("dcs")
        # ~ get_outFiles(sampType="dcs", suffix=".filt.no_overlap.bam")
    output:
        temp(touch(".ruleAllFinished"))
    shell:
        """
        
        rm -rf */picardTempDir
        rm -f Rplots.pdf
        """
        

rule rerun:
    input:
        ".rerunPrepDone",
        ".ruleAllFinished"
    shell:
        """
        rm .rerunPrepDone
        """

# Environment set up rules
# Implemented to save time at run-time on environment setup
rule initializeEnvs:
    input:
        "full.initialized"
        #~ "gatk.initialized",
        #~ "DS.initialized",
        #~ "bwa.initialized",
        #~ "fgbio.initialized"
    params:
        basePath = sys.path[0]
    output:
        ".env_initialized",
    shell:
        """
        touch "{output}"
        """

rule initializeFullEnv:
    output:
        temp("full.initialized")
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        echo "Initializing envs"
        touch full.initialized
        """

rule initializeDS_env:
    output:
        temp("DS.initialized")
    conda:
        "envs/DS_env.yaml"
    shell:
        """
        touch DS.initialized
        """
rule initialize bwa_env:
    output:
        temp("bwa.initialized")
    conda:
        "envs/bwa_env.yaml"
    shell:
        """
        touch bwa.initialized
        """
rule initialize fgbio_env:
    output:
        temp("fgbio.initialized")
    conda:
        "envs/fgbio_env.yaml"
    shell:
        """
        touch fgbio.initialized
        """

# Setup output directories
rule makeDirs:
    output:
        outConfigTmp = temp(touch("{runPath}/.{sample}_dirsMade")), 
    shell:
        """
        cd {wildcards.runPath}
        mkdir -p Intermediate/ConsensusMakerOutputs \
        Intermediate/postBlast/FilteredReads \
        Final/dcs Final/sscs Final/dcs/FilteredReads \
        Stats/data Stats/plots \
        Intermediate/PreVariantCallsCp/sscs \
        Intermediate/PreVariantCallsCp/dcs
        cd ../
        """

# Make consensus sequences from input fastqs
rule makeConsensus:
    params:
        sample = get_sample,
        cutOff = get_cutOff,
        nCut = get_nCutOff,
        umiLen = get_umiLen,
        spacerLen = get_spacerLen,
        locLen = get_locLen,
        basePath = sys.path[0],
        runPath = get_baseDir, 
        minmem = get_minMem, 
        maxmem = get_maxMem
    input:
        in1=get_in1,
        in2=get_in2,
        dirsMade="{runPath}/.{sample}_dirsMade"
    output:
        r1_sscs = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_read1_sscs.fq.gz",
        r2_sscs = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_read2_sscs.fq.gz",
        r1_dcs = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_read1_dcs.fq.gz",
        r2_dcs = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_read2_dcs.fq.gz",
        umi_proc = temp("{runPath}/{sample}.temp.sort.bam"),
        tagstats = "{runPath}/Stats/data/{sample}.tagstats.txt",
        tagPlot1 = "{runPath}/Stats/plots/{sample}_family_size.png",
        tagPlot2 = "{runPath}/Stats/plots/{sample}_fam_size_relation.png",
        outAln1 = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_aln_seq1.fq.gz",
        outAln2 = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_aln_seq2.fq.gz",
        outStats = "{runPath}/Stats/data/{sample}_cmStats.txt"
    priority: 50
    conda:
       "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_makeConsensus.log"
    threads: min(max(int(config["maxCores"]/2), 1),4)
    shell:
        """
        set -x
        cd {params.runPath}
        pwd
        ls
        picard FastqToSam \
        F1=../{input.in1} \
        F2=../{input.in2} \
        O=/dev/stdout \
        SM={params.sample} \
        TMP_DIR=picardTempDir \
        SORT_ORDER=unsorted \
        | python3 {params.basePath}/scripts/UnifiedConsensusMaker.py \
        --input /dev/stdin \
        --taglen {params.umiLen} \
        --spacerlen {params.spacerLen} \
        --loclen {params.locLen} \
        --write-sscs \
        --prefix {params.sample} \
        --tagstats \
        --cutoff {params.cutOff} \
        --Ncutoff {params.nCut} \
        --numCores {threads} \
        --minmem {params.minmem} \
        --maxmem {params.maxmem}
        
        mv {wildcards.sample}_*.fq.gz Intermediate/ConsensusMakerOutputs/
        mv {wildcards.sample}_cmStats.txt Stats/data/
        mv {wildcards.sample}.tagstats.txt Stats/data/
        mv {wildcards.sample}_fam*.png Stats/plots/
        cd ../
        """

# Calculate # on target raw reads based on outAln1 and outAln2 from makeConsensus
rule getOnTarget:
    params:
        basePath = sys.path[0],
        sample = get_sample,
        runPath = get_baseDir,
        rgpu = get_rgpu,
        rgpl = get_rgpl,
        rgsm = get_rgsm,
        rglb = get_rglb
    input:
        in1 = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_aln_seq1.fq.gz",
        in2 = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_aln_seq2.fq.gz",
        inRef = get_reference,
        inBed = get_target_bed
    output:
        outBam = temp("{runPath}/{sample}_mem.aln.sort.bam"),
        outBai = temp("{runPath}/{sample}_mem.aln.sort.bam.bai"),
        outOnTarget = "{runPath}/Stats/data/{sample}_onTargetCount.txt"
    conda:
       "envs/DS_env_full.yaml"
    threads: min(max(int(config["maxCores"]/2), 1),4)
    shell:
        """
        set -x
        cd {params.runPath}
        bwa mem \
        -L 1,1 \
        -C \
        -R "@RG\\tID:{params.sample}\\tLB:{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tSM:{params.rgsm}" \
        {input.inRef} \
        ../{input.in1} \
        ../{input.in2} \
        | samtools sort  -o ../{output.outBam} -
        samtools index ../{output.outBam}
        echo "$(samtools view -c -L {input.inBed} ../{output.outBam}) reads on target" > ../{output.outOnTarget}
        echo "$(samtools view -c ../{output.outBam}) total reads" >> ../{output.outOnTarget}
        cd ../
        """

# Align SSCS and DCS to provided reference genome, and sort based on position
rule alignReads:
    params:
        sample = get_sample,
        rgpu = get_rgpu,
        rgpl = get_rgpl,
        rgsm = get_rgsm,
        rglb = get_rglb,
        basePath = sys.path[0],
        runPath = get_baseDir
    priority: 45
    input:
        in1 = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_read1_{sampType}.fq.gz",
        in2 = "{runPath}/Intermediate/ConsensusMakerOutputs/{sample}_read2_{sampType}.fq.gz",
        inRef = get_reference
    output:
        outBam = temp("{runPath}/{sample}_mem.{sampType}.sort.bam"),
        outBai = temp("{runPath}/{sample}_mem.{sampType}.sort.bam.bai"),
        tempDir = temp(touch(directory("{runPath}/{sample}.{sampType}.alignReads.samtoolsTemp")))
    conda:
       "envs/DS_env_full.yaml"
    threads: min(max(int(config["maxCores"]/2), 1),4)
    log:
         "{runPath}/logs/{sample}_bwa_{sampType}.log"
    shell:
        """
        set -x
        mkdir {wildcards.runPath}/{wildcards.sample}.{wildcards.sampType}.alignReads.samtoolsTemp
        cd {params.runPath}
        bwa mem \
        -L 1,1 \
        -C \
        -R "@RG\\tID:{params.sample}\\tLB:{params.rglb}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\\tSM:{params.rgsm}" \
        {input.inRef} \
        ../{input.in1} \
        ../{input.in2} \
        | samtools sort \
        -T {wildcards.sample}.{wildcards.sampType}.alignReads.samtoolsTemp \
        -o {wildcards.sample}_mem.{wildcards.sampType}.sort.bam -
        samtools index {wildcards.sample}_mem.{wildcards.sampType}.sort.bam
        cd ../
        """

# Filter out secondary and supplamentery alignments prior to BLAST filters.
rule PreBlastFilter:
    input:
        inBam="{runPath}/{sample}_mem.dcs.sort.bam",
        inBai="{runPath}/{sample}_mem.dcs.sort.bam.bai",
        inRef = get_reference, 
        inBlastDb = get_blast_db_path
    output:
        outBam1 = temp("{runPath}/{sample}_mem.dcs.nonSecSup.bam"),
        outBai1 = temp("{runPath}/{sample}_mem.dcs.nonSecSup.bam.bai"),
        outBam2 = temp("{runPath}/{sample}_mem.dcs.SecSup.bam")
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_PreBlastFilter_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}
        samtools view -b -F 0x900 \
        -U {wildcards.sample}_mem.dcs.SecSup.bam \
        -o {wildcards.sample}_mem.dcs.nonSecSup.bam \
        {wildcards.sample}_mem.dcs.sort.bam
        samtools index {wildcards.sample}_mem.dcs.nonSecSup.bam
        cd ../
        """

# Make bam index for temporary bam files
rule makeTempBai:
    input:
        inBam = "{runPath}/{fileBase}.temp.bam"
    output:
        outBai = temp("{runPath}/{fileBase}.temp.bam.bai")
    conda:
       "envs/DS_env_full.yaml"
    shell:
        """
        set -x
        samtools index {input.inBam} {output.outBai}
        """

# make bam index for non-temporary bam files
rule makeBai:
    input:
        inBam = "{runPath}/{fileBase}.bam"
    output:
        outBai = "{runPath}/{fileBase}.bam.bai"
    conda:
       "envs/DS_env_full.yaml"
    shell:
        """
        set -x
        samtools index {input.inBam} {output.outBai}
        """

# Get samtools flagstat output for a bam file
rule getFlagstats:
    input:
        inBam = "{runPath}/{fileBase}.bam"
    output:
        outBai = "{runPath}/Stats/data/{fileBase}.flagstats.txt"
    conda:
       "envs/DS_env_full.yaml"
    shell:
        """
        set -x
        cd {wildcards.runPath}
        samtools flagstat {wildcards.fileBase}.bam \
        > Stats/data/{wildcards.fileBase}.flagstats.txt
        cd ../
        """

# Do naive variant calling to allow us to not include reads that just have 
# SNPs in them.  
rule PreBlastProcessing1:
    params:
        basePath = sys.path[0],
        rgsm = get_rgsm
    priority: 43
    input:
        inBam="{runPath}/{sample}_mem.dcs.nonSecSup.bam",
        inBai="{runPath}/{sample}_mem.dcs.nonSecSup.bam.bai",
        inRef = get_reference
    output:
        tempVars = temp("{runPath}/{sample}_dcs.vars.vcf"), 
        tempDepth = temp("{runPath}/{sample}_dcs.vars.vcf_depth.txt")
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_preBlast1_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}
        python3 {params.basePath}/scripts/BamToMutposVCF.py \
        --inBam {wildcards.sample}_mem.dcs.nonSecSup.bam \
        --inFasta {input.inRef} -o {wildcards.sample}_dcs.vars.vcf \
        --samp_name {params.rgsm}
        cd ../
        """

# Get SNPs from naive variant calling
rule PreBlastProcessing2:
    params:
        basePath = sys.path[0]
    priority: 43
    input:
        inVCF="{runPath}/{sample}_dcs.vars.vcf",
    output:
        tempMarked = temp("{runPath}/{sample}_dcs.Markedvars.vcf"),
        tempSnps = temp("{runPath}/{sample}_dcs.snps.vcf")
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_preBlast2_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}
        python3 {params.basePath}/scripts/SNP_finder.py \
        --in_file {wildcards.sample}_dcs.vars.vcf \
        -o {wildcards.sample}_dcs.Markedvars.vcf \
        -s {wildcards.sample}_dcs.snps.vcf
        cd ../
        """

# Separate out reads with non-SNP variants for examination with BLAST
rule PreBlastProcessing3:
    params:
        basePath = sys.path[0],
        rgsm = get_rgsm,
        readLength = get_readLen
    input:
        inBam="{runPath}/{sample}_mem.dcs.nonSecSup.bam",
        inBai="{runPath}/{sample}_mem.dcs.nonSecSup.bam.bai",
        inSnps = "{runPath}/{sample}_dcs.snps.vcf"
    output:
        tempBam1 = "{runPath}/Intermediate/postBlast/{sample}_dcs.preBlast.unmutated.bam",
        tempBam2 = "{runPath}/Intermediate/postBlast/{sample}_dcs.preBlast.mutated.bam", 
        tempMPC = temp("{runPath}/Intermediate/postBlast/{sample}_dcs.snpFiltered_MutsPerCycle.dat.csv")
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_preBlast3_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}
        python3 {params.basePath}/scripts/countMutsPerCycle.py  \
        --inFile {wildcards.sample}_mem.dcs.nonSecSup.bam \
        --inSnps {wildcards.sample}_dcs.snps.vcf \
        -o Intermediate/postBlast/{wildcards.sample}_dcs.snpFiltered \
        -l {params.readLength} -g -b -t 0
        mv Intermediate/postBlast/{wildcards.sample}_dcs.snpFiltered.badReads.t0.bam \
        Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.mutated.bam
        mv Intermediate/postBlast/{wildcards.sample}_dcs.snpFiltered.goodReads.t0.bam \
        Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.unmutated.bam
        
        cd ../
        """

# BLAST reads with non-SNP variants
rule BLAST:
    params:
        basePath = sys.path[0],
        db = get_blast_db
    priority: 42
    threads: config["maxCores"]
    input:
        inBam1="{runPath}/Intermediate/postBlast/{sample}_dcs.preBlast.mutated.bam",
        inBlastDbPath = get_blast_db_path
    output:
        outXML = "{runPath}/Intermediate/postBlast/{sample}_dcs.blast.xml",
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_blast_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}

        samtools fasta Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.mutated.bam | \
        blastn -task blastn \
        -db {params.db} \
        -outfmt 5 \
        -max_hsps 2 \
        -max_target_seqs 2 \
        -num_threads {threads} \
        | python3 {params.basePath}/scripts/blastMonitor.py \
        > Intermediate/postBlast/{wildcards.sample}_dcs.blast.xml

        cd ../
        """

# Label reads with taxIDs based on BLAST results
rule PostBlastProcessing1:
    params:
        basePath = sys.path[0]
    priority: 41
    input:
        inBam1="{runPath}/Intermediate/postBlast/{sample}_dcs.preBlast.mutated.bam",
        inXML = ancient("{runPath}/Intermediate/postBlast/{sample}_dcs.blast.xml"),
    output:
        tempBam3 = temp("{runPath}/{sample}_dcs.speciesLabeled.bam"),
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_postBlast1_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}
        python3 {params.basePath}/scripts/blastFilter.py \
        Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.mutated.bam \
        Intermediate/postBlast/{wildcards.sample}_dcs.blast.xml \
        {wildcards.sample}_dcs
        cd ../
        """
# Multi-step process.  
# 1. Merge species-labeled and unblasted reads into a single bam file
# 2. Sort reads based on read name
# 3. Filter out incorrect species reads
# 4. Sort all outputs into position sorting
rule PostBlastProcessing2:
    params:
        basePath = sys.path[0],
        taxID = get_target_taxon
    priority: 41
    input:
        inBam1="{runPath}/{sample}_dcs.speciesLabeled.bam",
        inBam2="{runPath}/Intermediate/postBlast/{sample}_dcs.preBlast.unmutated.bam",
    output:
        tempBam4 = temp("{runPath}/{sample}_dcs.wrongSpecies.bam"),
        tempAmbigBam = temp("{runPath}/{sample}_dcs.ambig.bam"),
        outBadBam = "{runPath}/Intermediate/postBlast/FilteredReads/{sample}_dcs.wrongSpecies.sort.bam",
        outBam = temp("{runPath}/{sample}_dcs.speciesFilt.sort.temp.bam"),
        outAmbBam = "{runPath}/Intermediate/postBlast/FilteredReads/{sample}_dcs.ambig.sort.bam",
        outSpecComp = "{runPath}/Stats/data/{sample}_dcs.speciesComp.txt",
        tempDir1 = temp(touch(directory("{runPath}/{sample}.dcs.postBlast1.samtoolsTemp"))),
        tempDir2 = temp(touch(directory("{runPath}/{sample}.dcs.postBlast2.samtoolsTemp"))),
        tempDir3 = temp(touch(directory("{runPath}/{sample}.dcs.postBlast3.samtoolsTemp"))),
        tempDir4 = temp(touch(directory("{runPath}/{sample}.dcs.postBlast4.samtoolsTemp")))
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_postBlast2_dcs.log"
    shell:
        """
        set -x
        cd {wildcards.runPath}

        mkdir {wildcards.sample}.dcs.postBlast1.samtoolsTemp
        mkdir {wildcards.sample}.dcs.postBlast2.samtoolsTemp
        mkdir {wildcards.sample}.dcs.postBlast3.samtoolsTemp

        samtools merge -c - \
        Intermediate/postBlast/{wildcards.sample}_dcs.preBlast.unmutated.bam \
        {wildcards.sample}_dcs.speciesLabeled.bam \
        | samtools sort -n \
        -T {wildcards.sample}.dcs.postBlast1.samtoolsTemp \
        | python3 {params.basePath}/scripts/postBlastFilter.py \
        {wildcards.sample}_dcs \
        {params.taxID} \
        | samtools sort \
        -T {wildcards.sample}.dcs.postBlast2.samtoolsTemp \
        -o {wildcards.sample}_dcs.speciesFilt.sort.temp.bam
        samtools sort -o Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.wrongSpecies.sort.bam \
        -T {wildcards.sample}.dcs.postBlast3.samtoolsTemp \
        {wildcards.sample}_dcs.wrongSpecies.bam
        samtools sort -o Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.ambig.sort.bam \
        -T {wildcards.sample}.dcs.postBlast4.samtoolsTemp \
        {wildcards.sample}_dcs.ambig.bam
        mv {wildcards.sample}_dcs.speciesComp.txt Stats/data/{wildcards.sample}_dcs.speciesComp.txt
        cd ../
        """

# Run user-defined recovery script
rule postBlastRecovery:
    params:
        basePath = sys.path[0],
        taxID = get_target_taxon
    priority: 40
    input:
        inAmbigBam="{runPath}/Intermediate/postBlast/FilteredReads/{sample}_dcs.ambig.sort.bam",
        inNonAmbigBam="{runPath}/{sample}_dcs.speciesFilt.sort.temp.bam",
        inWrongSpeciesBam="{runPath}/Intermediate/postBlast/FilteredReads/{sample}_dcs.wrongSpecies.sort.bam",
        inRecoveryScript=get_recovery
    output:
        outBam = temp("{runPath}/{sample}_dcs.postRecovery.recovered.temp.bam"),
        outAmbigBam = "{runPath}/Final/dcs/FilteredReads/{sample}_dcs.postRecovery.ambig.bam",
        outWrongSpeciesBam = "{runPath}/Final/dcs/FilteredReads/{sample}_dcs.postRecovery.wrongSpecies.bam"
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_postBlast2_dcs.log"
    shell:
        """
        cd {wildcards.runPath}
        bash {input.inRecoveryScript} \
        Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.ambig.sort.bam \
        {wildcards.sample}_dcs.speciesFilt.sort.temp.bam \
        Intermediate/postBlast/FilteredReads/{wildcards.sample}_dcs.wrongSpecies.sort.bam \
        {wildcards.sample}_dcs.postRecovery \
        "{params.basePath}"
        mv {wildcards.sample}_dcs.postRecovery.ambig.bam \
        Final/dcs/FilteredReads/{wildcards.sample}_dcs.postRecovery.ambig.bam
        mv {wildcards.sample}_dcs.postRecovery.wrongSpecies.bam \
        Final/dcs/FilteredReads/{wildcards.sample}_dcs.postRecovery.wrongSpecies.bam
        cd ../
        """

# Count number of reads in files after post-blast recovery with different 
# ambiguity scores
rule CountAmbig:
    params:
        basePath = sys.path[0],
    input:
        inNonAmbigFile = "{runPath}/{sample}_dcs.postRecovery.recovered.temp.bam",
        inAmbigFile = "{runPath}/Final/dcs/FilteredReads/{sample}_dcs.postRecovery.ambig.bam"
    output:
        "{runPath}/Stats/data/{sample}.dcs_ambiguity_counts.txt"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {wildcards.runPath}
        python3 {params.basePath}/scripts/countAmbiguityClasses.py \
        Stats/data/{wildcards.sample}.dcs \
        {wildcards.sample}_dcs.postRecovery.recovered.temp.bam \
        Final/dcs/FilteredReads/{wildcards.sample}_dcs.postRecovery.ambig.bam
        cd ../
        """

rule makePreVariantCallCp:
    params:
    input:
        inBam = getVarDictBam,
    output:
        outBam = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cp {input.inBam} {output.outBam}
        """

rule varDict:
    params:
        clip5 = get_clipBegin,
        clip3 = get_clipEnd,
        readLength = get_readLen,
        umiLen = get_umiLen,
        spacerLen = get_spacerLen,
        sampName = get_rgsm, 
        basePath = sys.path[0],
        vardict_f = config["vardict_f"],
        vardict_nmfreq = config["vardict_nmfreq"],
        vardict_r = config["vardict_r"],
        vardict_V = config["vardict_V"],
        vardict_adaptor=config["vardict_adaptor"]
    input:
        inBam = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam", 
        inBai = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam.bai",
        inBed = get_target_bed,
        inRef = get_reference
    output:
        outVars = temp("{runPath}/{sample}.{sampType}.varDict.txt"), 
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {wildcards.runPath}
        vardict-java -b Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.final.bam \
        -UN -P {params.clip5} \
        -T $(( {params.readLength}-{params.umiLen}-{params.spacerLen}-{params.clip3} )) \
        -f {params.vardict_f} -p \
        -G {input.inRef} \
        -nmfreq {params.vardict_nmfreq} \
        -N {params.sampName} \
        -r {params.vardict_r} \
        -v -c 1 -S 2 -E 3 -U -g 4 \
        {input.inBed} \
        -h -V {params.vardict_V} \
        --adaptor {params.vardict_adaptor} \
        > {wildcards.sample}.{wildcards.sampType}.varDict.txt
        cd ..
        """

rule varDict_Ns:
    params:
        clip5 = get_clipBegin,
        clip3 = get_clipEnd,
        readLength = get_readLen,
        umiLen = get_umiLen,
        spacerLen = get_spacerLen,
        sampName = get_rgsm, 
        basePath = sys.path[0],
        vardict_f = config["vardict_f"],
        vardict_nmfreq = config["vardict_nmfreq"],
        vardict_r = config["vardict_r"],
        vardict_V = config["vardict_V"],
        vardict_adaptor=config["vardict_adaptor"]
    input:
        inBam = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam", 
        inBai = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam.bai",
        inBed = get_target_bed,
        inRef = get_reference
    output:
        outVars = temp("{runPath}/{sample}.{sampType}.varDict.Ns.txt"), 
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {wildcards.runPath}
        vardict-java -b Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.final.bam \
        -UN -P {params.clip5} \
        -T $(( {params.readLength}-{params.umiLen}-{params.spacerLen}-{params.clip3} )) \
        -f {params.vardict_f} -p -K \
        -G {input.inRef} \
        -nmfreq {params.vardict_nmfreq} \
        -N {params.sampName} \
        -r {params.vardict_r} \
        -v -c 1 -S 2 -E 3 -U -g 4 \
        {input.inBed} \
        -h -V {params.vardict_V} \
        --adaptor {params.vardict_adaptor} \
        > {wildcards.sample}.{wildcards.sampType}.varDict.Ns.txt
        cd ..
        """

rule varDict2VCF:
    params:
        basePath = sys.path[0],
        sampName = get_rgsm, 
        minDepth = get_minDepth, 
    input:
        inVarDict = "{runPath}/{sample}.{sampType}.varDict.txt", 
        inVarDict_Ns = "{runPath}/{sample}.{sampType}.varDict.Ns.txt", 
        inBam = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam", 
        inBai = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam.bai",
    output:
        outVCF = "{runPath}/Final/{sampType}/{sample}.{sampType}.vcf",
        outSNPs = "{runPath}/Final/{sampType}/{sample}.{sampType}.snps.vcf"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {wildcards.runPath}
        python {params.basePath}/scripts/varDictToVCF.py \
        -i {wildcards.sample}.{wildcards.sampType}.varDict.txt \
        -n {wildcards.sample}.{wildcards.sampType}.varDict.Ns.txt \
        -b Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.final.bam \
        -o Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.vcf \
        -s Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.snps.vcf \
        --samp_name {params.sampName} \
        -d {params.minDepth}
        cd ..
        """

rule makeDepth:
    params:
        basePath = sys.path[0],
        runPath = get_baseDir,
        sampName = get_rgsm,
    input:
        inBam = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam",
        inBai = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam.bai",
        inRef = get_reference
    output:
        outDepth = "{runPath}/Stats/data/{sample}.{sampType}.depth.txt"
    conda:
        "envs/DS_env_full.yaml"
    log:
         "{runPath}/logs/{sample}_makeDepth_{sampType}.log"
    shell:
        """
        cd {params.runPath}
        python {params.basePath}/scripts/makeDepthFile.py \
        -i ../{input.inBam} \
        -f {input.inRef} \
        -o {wildcards.sample}.{wildcards.sampType}.depth.txt \
        --samp_name {params.sampName}
        
        mv {wildcards.sample}.{wildcards.sampType}.depth.txt Stats/data/{wildcards.sample}.{wildcards.sampType}.depth.txt
        cd ..
        """

# Make countMuts file from VCF file
rule makeCountMuts:
    params:
        basePath = sys.path[0],
        runPath = get_baseDir,
        minClonal = get_minClonal,
        maxClonal = get_maxClonal,
        minDepth = get_minDepth,
        maxNs = get_maxNs,
        sampName = get_rgsm,
        cm_outputs = get_cm_outputs,
        cm_sums = get_cm_sumTypes
    input:
        inVCF = "{runPath}/Final/{sampType}/{sample}.{sampType}.vcf",
        inBam = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam",
        inBai = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam.bai",
        inRef = get_reference, 
        inBed = get_target_bed
    output:
        outCountMuts = "{runPath}/Final/{sampType}/{sample}.{sampType}.countmuts.csv", 
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {params.runPath}
        # BamToCountMuts:
        python {params.basePath}/scripts/bamToCountMuts.py \
        --samp_name {params.sampName} \
        -i ../{input.inBam} \
        -f {input.inRef} \
        -b {input.inBed} \
        -d {params.minDepth} \
        -v Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.vcf \
        --outputType {params.cm_outputs} \
        --sumType {params.cm_sums} \
        -c {params.minClonal} \
        -C {params.maxClonal} \
        -n {params.maxNs} \
        -u \
        -o Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.countmuts.csv
        cd ..
        """

# Collect insert size data
rule InsertSize:
    params:
        sample = get_sample,
        basePath = sys.path[0],
        runPath = get_baseDir,
        readLength = get_readLen
    input:
        inRef = get_reference,
        inAligned = "{runPath}/{sample}_mem.{sampType}.sort.bam",
        inAlignedBai = "{runPath}/{sample}_mem.{sampType}.sort.bam.bai"
    output:
        outMetrics = touch("{runPath}/Stats/data/{sample}.{sampType}.iSize_Metrics.txt"),
        out_iSizeHist = temp(touch("{runPath}/{sample}.{sampType}.iSize_Histogram.pdf")),
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_stats_{sampType}.log"
    shell:
        """
        cd {params.runPath}
        # Plot insert-size histogram (using unfiltered and unclipped data)
        echo "# Dummy insert size file" > Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "## HISTOGRAM" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "insert_size	All_Reads.fr_count" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        echo "0	0" >> Stats/data/{wildcards.sample}.{wildcards.sampType}.iSize_Metrics.txt
        
        picard CollectInsertSizeMetrics \
        I=../{input.inAligned} \
        O=../{output.outMetrics} \
        H=../{output.out_iSizeHist} M=0.5 \
        TMP_DIR=picardTempDir
        cd ../
        """

# Create insert size plot
rule PlotInsertSize:
    params:
        basePath = sys.path[0],
    input:
        "{runPath}/Stats/data/{sample}.{sampType}.iSize_Metrics.txt",
    output:
        temp("{runPath}/Stats/plots/{sample}.{sampType}.iSize_Histogram.tiff")
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {wildcards.runPath}
        Rscript {params.basePath}/scripts/plotInsertSize.R \
        {wildcards.sample}.{wildcards.sampType}
        cd ..
        """

# Create coverage plot
rule PlotCoverage:
    params:
        basePath = sys.path[0],
    input:
        inBed = get_target_bed,
        inDepth = "{runPath}/Stats/data/{sample}.{sampType}.depth.txt",
        inVCF = "{runPath}/Final/{sampType}/{sample}.{sampType}.vcf"
    output:
        temp("{runPath}/Stats/plots/{sample}.{sampType}.targetCoverage.tiff")
    conda:
        "envs/DS_env_full.yaml"
    log:
        "{runPath}/logs/{sample}_stats_{sampType}.log"
    shell:
        """
        cd {wildcards.runPath}
        Rscript {params.basePath}/scripts/MakeDepthPlot.R \
        {input.inBed} \
        {wildcards.sample}.{wildcards.sampType} \
        {wildcards.sampType}
        cd ..
        """

# Plot number of non-SNP variants per sequencing cycle
rule MutsPerCycle:
    params:
        sample = get_sample,
        basePath = sys.path[0],
        runPath = get_baseDir,
        readLength = get_readLen
    input:
        inRef = get_reference,
        inFinal = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam",
        inFinalBai = "{runPath}/Final/{sampType}/{sample}.{sampType}.final.bam.bai",
        inSnps = "{runPath}/Final/{sampType}/{sample}.{sampType}.snps.vcf"
    output:
        outBam2 = "{runPath}/Final/{sampType}/{sample}.{sampType}.mutated.bam",
        outErrPerCyc2_WN = "{runPath}/Stats/plots/{sample}.{sampType}_BasePerPosInclNs.png",
        outErrPerCyc2_WoN = "{runPath}/Stats/plots/{sample}.{sampType}_BasePerPosWithoutNs.png",
        outErrPerCycDat2 = "{runPath}/Stats/data/{sample}.{sampType}_MutsPerCycle.dat.csv",
        outMutsPerRead = "{runPath}/Stats/data/{sample}.{sampType}.mutsPerRead.txt",
        outMutsPerReadPlot = "{runPath}/Stats/plots/{sample}.{sampType}.mutsPerRead.png"
    conda:
         "envs/DS_env_full.yaml"
    log:
         "{runPath}/logs/{sample}_stats_{sampType}.log"
    shell:
        """
        cd {params.runPath}
        python3 {params.basePath}/scripts/countMutsPerCycle.py  \
        --inFile Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.final.bam \
        --inSnps Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.snps.vcf \
        -o {wildcards.sample}.{wildcards.sampType} \
        -l {params.readLength} -b -t 0 -c --text_file
        
        mv {wildcards.sample}.{wildcards.sampType}*.png Stats/plots/
        mv {wildcards.sample}.{wildcards.sampType}_MutsPerCycle.dat.csv Stats/data/
        mv {wildcards.sample}.{wildcards.sampType}.mutsPerRead.txt Stats/data/
        mv {wildcards.sample}.{wildcards.sampType}.badReads.t0.bam Final/{wildcards.sampType}/{wildcards.sample}.{wildcards.sampType}.mutated.bam
        cd ../
        """

# Convert TIFF figures created by R into PNG figures
# We implemented this due to artifacts in the raw PNG output produced by
# R, and a desire to embed the resulting figures in the report html file.
# This issue may be revisited in the future, if we come accross a better 
# PNG device for R.  
rule tiff2png:
    params:
        basePath = sys.path[0],
    input:
        "{prefix}.tiff"
    output:
        "{prefix}.png"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        python {params.basePath}/scripts/tiff2png.py \
        {wildcards.prefix}.tiff \
        {wildcards.prefix}.png
        """

# make summary files
rule makeSummaryCSV:
    params:
        basePath = sys.path[0],
        configPath = config["samples"]
    input:
        getSummaryInput()
    output:
        outSum = f"{config['samples']}.summary.csv"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        python {params.basePath}/scripts/retrieveSummary.py \
        --config {params.configPath}
        """

rule makeSummaryDepth:
    params:
        basePath = sys.path[0],
        configPath = config["samples"]
    input:
        getDepthFiles()
    output:
        outSum = f"{config['samples']}.summaryDepth.pdf"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        python {params.basePath}/scripts/plotSummaryDepth.py \
        {params.configPath}
        """

rule makeSummaryInsertSize:
    params:
        basePath = sys.path[0],
        configPath = config["samples"]
    input:
        getInsertFiles()
    output:
        outSum = f"{config['samples']}.summaryInsertSize.pdf"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        python {params.basePath}/scripts/plotSummaryInsertSize.py \
        {params.configPath}
        """
        
rule makeSummaryMutsByCycle:
    params:
        basePath = sys.path[0],
        configPath = config["samples"]
    input:
        getMutsByCycFiles()
    output:
        outSum = f"{config['samples']}.summaryMutsByCycle.pdf"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        Rscript {params.basePath}/scripts/plotSummaryMutsByCycle.R \
        {params.configPath}
        """

rule makeSummaryFamilySize:
    params:
        basePath = sys.path[0],
        configPath = config["samples"]
    input:
        getFamSizeFiles()
    output:
        outSum = f"{config['samples']}.summaryFamilySize.pdf"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        Rscript {params.basePath}/scripts/plotSummaryFamilySize.R \
        {params.configPath}
        """

# make per-sample report file
rule makeReport:
    params:
        sample = get_sample,
        rglb = get_rglb,
        rgpl = get_rgpl,
        rgsm = get_rgsm,
        rgpu = get_rgpu,
        reference = get_reference,
        target_bed = get_target_bed,
        baseDir = get_baseDir,
        in1 = get_in1,
        in2 = get_in2,
        mqFilt = get_mqFilt,
        minMem = get_minMem,
        maxMem = get_maxMem,
        cutOff = get_cutOff,
        nCutOff = get_nCutOff,
        umiLen = get_umiLen,
        spacerLen = get_spacerLen,
        locLen = get_locLen,
        clipBegin = get_clipBegin,
        clipEnd = get_clipEnd,
        minClonal = get_minClonal,
        maxClonal = get_maxClonal,
        minDepth = get_minDepth,
        maxNs = get_maxNs,
        outConfig = getOutConfig,
        rLen = get_readLen,
        contaminantDb = get_blast_db, 
        version = next(open(f"{sys.path[0]}/VERSION",'r')).strip()
    input:
        getReportInput
    output:
        touch("{runPath}/Stats/{sample}.report.ipynb")
    run:
        import nbformat as nbf
        #from matplotlib import pyplot as plt
        import numpy as np
        nb = nbf.v4.new_notebook()
        myCells = []
        myCells.append(nbf.v4.new_code_cell("""\
%pylab inline
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
from IPython.display import Image
import numpy as np
"""))
        myCells.append(nbf.v4.new_markdown_cell("""\
#Duplex Sequencing Summary"""
            ))

        # ToC
        myCells.append(nbf.v4.new_markdown_cell("""\
##Table of Contents:
[Top](#Duplex-Sequencing-Summary)  
1. [Table of Contents](#Table-of-Contents:)
2. [Glosary](#Glosary:)  
3. [Parameters](#Parameters:)
4. [Consensus Maker Statistics](#Consensus-Maker-Statistics:)
5. [Family Size Plots](#Family-Size-Plots:)
6. [Read Statistics](#Read-Statistics:)
7. [Alignment Statistics](#Alignment-Statistics:)
8. [Consensus Making Ratios](#Consensus-Making-Ratios:)
9. [BLAST Statistics](#BLAST-Statistics:)
10. [Insert Size Graph](#Insert-size-graph:)
11. [Depth per Target](#Depth-per-Target:)
12. [Muts per Cycle](#Muts-per-Cycle:)
13. [Countmuts output](#Countmuts-output:)
"""
            ))
        # Glossary
        myCells.append(nbf.v4.new_markdown_cell(
            f'##Glosary:  \n'
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f'###Read:  \n'
            f'A single DNA sequence; one half of an Illumina paired-end read.  \n'
            f'###Paired-end read:  \n'
            f'A pair of DNA sequences from the same molecule in the final library; analogous to cluster. Most sequencing facilities will refer to paired-end reads as reads.  \n'
            f'###Family:  \n'
            f'A group of reads originating from the same end of the same strand of the same original (pre-library preperation) DNA molecule.  \n'
            f'###SSCS:  \n'
            f'A consensus sequence made by comparing all reads in a family at each base, and selecting the most common base at each position for the consensus base if it matches the stringency requirement of what proportion of reads match the base.  \n'
            f'###DCS:  \n'
            f'A consensus made by comparing two SSCS from the same end of the same molecule; each SSCS molecule represents one strand of the molecule.  \n'
            f'###Depth:  \n'
            f'The number of DCS reads, present at a given position in the reference genome; equivalent to the number of original molecules sequenced.  \n'
            ))

        # Parameters
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Parameters:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f"  \n"
            f"| Parameter | Value |  \n"
            f"| --------- | ----- |  \n"
            f"|Version|{params.version}|  \n"
            f"|RunID |{wildcards.sample}|  \n"
            f"|Run Directory | {wildcards.runPath} |  \n"
            f"|Inputs| |  \n"
            f"|R1| {params.in1}|  \n"
            f"|R2| {params.in2}|  \n"
            f"|Read Length| {params.rLen}|  \n"
            f"|UMI Length| {params.umiLen}|  \n"
            f"|Spacer Length| {params.spacerLen}|  \n"
            f"|Localization Length| {params.locLen}|  \n"
            f"|Minimum family size| {params.minMem}|  \n"
            f"|Maximum family size| {params.maxMem}|  \n"
            f"|Consensus Making Stringency| {params.cutOff}|  \n"
            f"|Max proportion Ns per DCS| {params.nCutOff}|  \n"
            f"|Genome| {params.reference}|  \n"
            f"|Bed File| {params.target_bed}|  \n"
            f"|Contaminant_filter DB| {params.contaminantDb}|  \n"
            f"|RGSM| {params.rgsm}|  \n"
            f"|RGLB| {params.rglb}|  \n"
            f"|RGPL| {params.rgpl}|  \n"
            f"|RGPU| {params.rgpu}|  \n"
            f"|Clipping |  \n"
            f"|5' | {params.clipBegin}|  \n"
            f"|3' | {params.clipEnd}|"
            ))
        # Read Consensus Maker statistics:
        cmStats = '  \n'.join([x.strip() for x in open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}_cmStats.txt", 'r').readlines()[1:]])
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Consensus Maker Statistics:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f"```\n{cmStats}```"
            ))
        # import family size plots
        myCells.append(nbf.v4.new_markdown_cell(
            f"###Family Size Plots:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            ))
        myCells.append(nbf.v4.new_code_cell(
            f'Image(f"plots/{wildcards.sample}_family_size.png")'
            ))
        myCells.append(nbf.v4.new_code_cell(
            f'Image(f"plots/{wildcards.sample}_fam_size_relation.png")'
            ))
        # Raw Read Counts:
        rawFlagstats = open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}.temp.sort.flagstats.txt", 'r').readlines()
        sscsFlagstats = open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}_mem.sscs.sort.flagstats.txt", 'r').readlines()
        dcsFlagstats = open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}_mem.dcs.sort.flagstats.txt", 'r').readlines()
        rawTarget = open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}_onTargetCount.txt", 'r').readlines()
        # Alignment Statistics:
        rawReads = int(rawFlagstats[0].split()[0])
        sscsReads=int(sscsFlagstats[0].split()[0])
        mappedSscs=int(sscsFlagstats[4].split()[0])
        dcsReads=int(dcsFlagstats[0].split()[0])
        mappedDcs=int(dcsFlagstats[4].split()[0])
        if int(rawTarget[1].split()[0]) == 0:
            rawOnTarget=0
        else:
            rawOnTarget=int(rawTarget[0].split()[0])/int(rawTarget[1].split()[0])
        if sscsReads == 0:
            sscsMapped=0
            raw_sscs=0
        else:
            sscsMapped=mappedSscs/sscsReads
            raw_sscs=rawReads/sscsReads
        if dcsReads == 0:
            dcsMapped=0
            sscs_dcs=0
        else:
            dcsMapped=mappedDcs/dcsReads
            sscs_dcs=sscsReads/dcsReads
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Read Statistics:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f"  \n"
            f"| | |  \n"
            f"| --- | --- |  \n"
            f"| Raw Reads | {rawReads} |  \n"
            f"| SSCS Reads: | {sscsReads} |  \n"
            f"| Mapped SSCS Reads: | {mappedSscs} |  \n"
            f"| DCS Reads: | {dcsReads} |  \n"
            f"| Mapped DCS Reads: | {mappedDcs} |  \n"
            f"##Alignment Statistics:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f"  \n"
            f"| | |  \n"
            f"| --- | --- |  \n"
            f"| Raw on target | {round(rawOnTarget,4)*100}% |  \n"
            f"| SSCS Mapped | {round(sscsMapped, 4)*100}% |  \n"
            f"| DCS Mapped | {round(dcsMapped, 4)*100}% |  \n"
            f"##Consensus Making Ratios:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f"  \n"
            f"| | |  \n"
            f"| --- | --- |  \n"
            f"| Raw/SSCS | {round(raw_sscs, 2)} |  \n"
            f"| SSCS/DCS | {round(sscs_dcs, 2)} |  \n"
            ))
        if params.contaminantDb.upper() != "NONE":
            # BLAST Statitics:
            # Table of species (taxIDs)
            blastSpecSrc = open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}_dcs.speciesComp.txt", 'r').readlines()
            blastSpecTab = [
                "| TaxID | Number Reads |  ",
                "| ----- | ------------ |  "
                ]
            for line in blastSpecSrc:
                linebins = line.strip().split()
                blastSpecTab.append(f"| {linebins[0]} | {linebins[1]} |  ")
            blastSpec = "\n".join(blastSpecTab)
            # Ambiguity Composition
            inAmbigs = "  \n".join([x.strip() for x in open(f"{wildcards.runPath}/Stats/data/{wildcards.sample}.dcs_ambiguity_counts.txt", 'r').readlines()])
            myCells.append(nbf.v4.new_markdown_cell(
                f"##BLAST Statistics:  \n"
                f"[Top](#Duplex-Sequencing-Summary)  \n"
                f"###Species Composition:  \n"
                f"To look up taxIDs, use [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)  \n"
                f" * TaxIDs of -1 are assigned to between species BLAST ties. \n"
                f" * TaxIDs of -3 are assigned to records with no BLAST results. \n"
                f" * TaxIDs of -4 are assigned to records which were submitted to BLAST, for which BLAST did not attempt alignment.  \n"
                f"  \n{blastSpec}\n"
                f"###Ambiguity Composition:  \n"
                f"```\n{inAmbigs}```"
                ))
        else:
            # BLAST Statitics:
            myCells.append(nbf.v4.new_markdown_cell(
                f"##BLAST Statistics:  \n"
                f"[Top](#Duplex-Sequencing-Summary)  \n"
                f"No blast was run on this sample.  \n"
                ))
        # Insert size graph
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Insert size graph:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            ))
        # ~ myCells.append(nbf.v4.new_code_cell(
            # ~ f'from IPython.display import IFrame\n'
            # ~ f'IFrame(f"{wildcards.sample}.dcs.iSize_Histogram.pdf",600,600)\n'
            # ~ ))
        myCells.append(nbf.v4.new_code_cell(
            f'Image(f"plots/{wildcards.sample}.dcs.iSize_Histogram.png")'
            ))

        # Depth statistics
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Depth per Target:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            ))
        myCells.append(nbf.v4.new_code_cell(
            f'Image(f"plots/{wildcards.sample}.dcs.targetCoverage.png")'
            ))

        # Muts per read pos statistics
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Muts per Cycle:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            ))
        myCells.append(nbf.v4.new_code_cell(
            f'Image(f"plots/{wildcards.sample}.dcs_BasePerPosWithoutNs.png")'
            ))

        # Mutation statistics
        cmTable1 = [
            "| Parameter | Value |  \n", 
            "| --------- | ----- |  \n"
            ]
        cmTable2  = [
            "| MUTATION_TYPE | MUTATION_CLASS | COUNT | DENOMINATOR | FREQUENCY |  \n",
            "| ------------- | -------------- | ----- | ----------- | --------- |  \n"
            ]
        cmFile = open(f"{wildcards.runPath}/Final/dcs/{wildcards.sample}.dcs.countmuts.csv", 'r')
        for line in cmFile:
            if "##" in line:
                linebins = line.strip().split('\t')
                if "Input file:" in line:
                    cmTable1.append(
                        f"| Input file | {linebins[1]} |  \n"
                        )
                elif "Input reference:" in line:
                    cmTable1.append(
                        f"| Input reference | {linebins[1]} |  \n"
                        )
                elif "Input bed:" in line:
                    cmTable1.append(
                        f"| Input bed | {linebins[1]} |  \n"
                        )
                elif "Minimum Depth" in line:
                    cmTable1.append(
                        f"| Minimum Depth | {linebins[1]} |  \n"
                        )
                elif "Clonality" in line: 
                    cmTable1.append(
                        f"| Clonality | {linebins[1]} |  \n"
                        )
                elif "Unique" in line:
                    cmTable1.append(
                        f"| Unique | True |  \n"
                        )
            elif "#" not in line:
                if "OVERALL" in line:
                    cmTable2.append(
                        f"| {' | '.join([x for x in line.strip().split(',')[2:]])} |  \n"
                        )
        cmFile.close()
        
        myCells.append(nbf.v4.new_markdown_cell(
            f"##Countmuts output:  \n"
            f"[Top](#Duplex-Sequencing-Summary)  \n"
            f"###Parameters:  \n"
            f"{''.join(cmTable1)}"
            f"  \n"
            f"###Overall Mutation Counts:  \n"
            f"{''.join(cmTable2)}"
            ))

        nb['cells'] = myCells
        nbf.write(nb, f"{wildcards.runPath}/Stats/{wildcards.sample}.report.ipynb")

# Compile report ipython notebook
rule compileReport:
    input:
        "{runPath}/Stats/{sample}.report.ipynb"
    output:
        "{runPath}/Final/{sample}.report.html"
    conda:
        "envs/DS_env_full.yaml"
    shell:
        """
        cd {wildcards.runPath}/Stats
        jupyter nbconvert --to notebook --execute --inplace {wildcards.sample}.report.ipynb
        jupyter nbconvert {wildcards.sample}.report.ipynb --no-input --stdout > ../Final/{wildcards.sample}.report.html
        cd ../../
        """

# Ruleorders
ruleorder: alignReads > makeBai
ruleorder: PreBlastFilter > makeBai
ruleorder: PreBlastProcessing1 > makeBai
ruleorder: PreBlastProcessing3 > makeBai
ruleorder: makeTempBai > makeBai
