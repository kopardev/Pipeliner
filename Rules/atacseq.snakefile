from os.path import join
import os


configfile: "run.json"
WORKDIR = config['project']['workpath']
os.system("touch "+join(WORKDIR,"dummy"))
SINGULARITY_SIF = config['bin'][pfamily]['SINGULARITY_SIF']
# SAMPLES=["test_input"]
GENOME = config['references'][pfamily]['GENOME']
# SCRIPTS_DIR="/opt2"
SCRIPTS_DIR = config['bin'][pfamily]['SCRIPTS_DIR']
N_REPLICATES_MAX= config['bin'][pfamily]['N_REPLICATES_MAX']
INDEX_DIR = config['references'][pfamily]['INDEX_DIR']

# ATACSeq only works with paired end data
if config['project']['nends'] != 2 :
    exit("ATACseq pipeline currently only supports Paired-End data")

# define functions
def check_existence(filename):
  if not os.path.exists(filename):
    exit("File: %s does not exists!"%(filename))

def check_readaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.R_OK):
    exit("File: %s exists, but cannot be read!"%(filename))

def check_writeaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.W_OK):
    exit("File: %s exists, but cannot be read!"%(filename))

# read groups and samples
check_readaccess(join(WORKDIR,"groups.tab"))
group2samples=dict()
sample2label=dict()
with open("groups.tab") as g:
    for line in g:
        line=line.strip().split("\t")
        if len(line)!=3:
            exit("File: groups.tab should have 3 tab-delimited columns!")
        if not line[1] in group2samples:
            group2samples[line[1]]=list()
        group2samples[line[1]].append(line[0])
        sample2label[line[0]]=line[2]
GROUPS=list(group2samples.keys())
SAMPLES=list(sample2label.keys())

# OTHER FUNCTIONS

def get_sample_tagAlignFiles_for_group(wildcards):
    d=dict()
    for i in range(N_REPLICATES_MAX):
        d["tagAlign"+str(i+1)]=join(WORKDIR,"dummy")
        d["genomefile"+str(i+1)]=join(WORKDIR,"dummy")
    for i,f in enumerate(group2samples[wildcards.grp]):
        d["tagAlign"+str(i+1)]=join(WORKDIR,"tagAlign",f+".tagAlign.gz")
        d["genomefile"+str(i+1)]=join(WORKDIR,"bam",f+".genome")
    return d

def get_sample_qsortedbamfiles_for_group(wildcards):
    d=dict()
    for i in range(N_REPLICATES_MAX):
        d["qsortedbam"+str(i+1)]=join(WORKDIR,"dummy")
    for i,f in enumerate(group2samples[wildcards.grp]):
        d["qsortedbam"+str(i+1)]=join(WORKDIR,"bam",f+".qsorted.bam")
    return d


def get_frip_input(wildcards):
    d=dict()
    for i in range(N_REPLICATES_MAX):
        d["tagAlign"+str(i+1)]=join(WORKDIR,"dummy")
        d["macs_sample_peaks"+str(i+1)]=join(WORKDIR,"dummy")
        d["genrich_sample_peaks"+str(i+1)]=join(WORKDIR,"dummy")
        d["macs_group_peaks"+str(i+1)]=join(WORKDIR,"dummy")
        d["genrich_group_peaks"+str(i+1)]=join(WORKDIR,"dummy")
    for i,f in enumerate(group2samples[wildcards.grp]):
        d["tagAlign"+str(i+1)]=join(WORKDIR,"tagAlign",f+".tagAlign.gz")
        d["macs_sample_peaks"+str(i+1)]=join(WORKDIR,"peaks","macs2",wildcards.grp+".sample.macs2.peakfiles")
        d["genrich_sample_peaks"+str(i+1)]=join(WORKDIR,"peaks","genrich",wildcards.grp+".sample.genrich.peakfiles")
        d["macs_group_peaks"+str(i+1)]=join(WORKDIR,"peaks","macs2",wildcards.grp+".group.macs2.peakfiles")
        d["genrich_group_peaks"+str(i+1)]=join(WORKDIR,"peaks","genrich",wildcards.grp+".group.genrich.peakfiles")
    return d


# RULES

rule all:
    input:
        expand(join(WORKDIR,"tagAlign","{sample}.tagAlign.gz"),sample=SAMPLES),
        expand(join(WORKDIR,"bam","{sample}.dedup.bam"),sample=SAMPLES),
        expand(join(WORKDIR,"bam","{sample}.genome"),sample=SAMPLES),
        expand(join(WORKDIR,"bam","{sample}.qsorted.bam"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R1_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R2_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","{sample}.nreads.txt"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R1.noBL_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R2.noBL_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","{sample}.dupmetric"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","preseq","{sample}.nrf"),sample=SAMPLES),
        expand(join(WORKDIR,"peaks","macs2","{grp}.group.macs2.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","macs2","{grp}.sample.macs2.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","macs2","{grp}.macs2.tn5knicksbedfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.group.genrich.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.sample.genrich.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.genrich.tn5knicksbedfiles"),grp=GROUPS),
        expand(join(WORKDIR,"QC","fld","{sample}.fld.txt"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","tss","{sample}.tss.txt"),sample=SAMPLES),
        join(WORKDIR,"QC","jaccard","macs2.sample.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","macs2.group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","macs2.sample_group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","genrich.sample.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","genrich.group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","genrich.sample_group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","allmethods.sample.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","allmethods.group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","allmethods.sample_group.jaccard.pca.html"),
        expand(join(WORKDIR,"QC","frip","{grp}.frip"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","{method}","motifs","{grp}.{method}.motiffiles"),grp=GROUPS,method=["macs2","genrich"]),
        join(WORKDIR,"QC","multiqc_report.html"),
        join(WORKDIR,"QC","QCStats.txt")


rule atac_tss:
    input:
        tagalign=join(WORKDIR,"tagAlign","{sample}.tagAlign.gz")
    output:
        tss=join(WORKDIR,"QC","tss","{sample}.tss.txt")
    params:
        rname='pl:tss',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        sample="{sample}"
    threads: 56
    shell:"""
set -e -x -o pipefail
module load singularity
rsync -Laz --progress {input.tagalign} /lscratch/$SLURM_JOBID/
cd /lscratch/$SLURM_JOBID
tagAlign=$(echo {input.tagalign}|awk -F"/" "{{print \$NF}}")
outfn=$(echo {output.tss}|awk -F"/" "{{print \$NF}}")

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/ \
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_tagAlign2TSS.bash \
--tagaligngz $tagAlign \
--genome {params.genome} \
--tsstxt $outfn \
--ncpus {threads} \
--scriptsfolder {params.scriptsdir}

rsync -az --progress $outfn {output.tss}

"""

rule atac_motifs:
    input:
        join(WORKDIR,"peaks","{method}","{grp}.group.{method}.peakfiles"),
        join(WORKDIR,"peaks","{method}","{grp}.sample.{method}.peakfiles")
    output:
        join(WORKDIR,"peaks","{method}","motifs","{grp}.{method}.motiffiles")
    params:
        rname='pl:motifs',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        grp="{grp}",
        method="{method}"
    threads: 56
    shell:"""
set -e -x -o pipefail
module load singularity
# export SINGULARITYENV_PYTHONNOUSERSITE="1"
for i in {input};do
    while read sample group filepath;do
        rsync -Laz --progress $filepath /lscratch/$SLURM_JOBID/
        cd /lscratch/$SLURM_JOBID
        narrowpeak_fn=$(echo $filepath|awk -F"/" "{{print \$NF}}")
        PYTHONNOUSERSITE=1 singularity exec --cleanenv \
        -B {params.workdir}/:/workdir/,{params.index_dir}/:/index,/lscratch/$SLURM_JOBID/:/data2 \
        {params.singularity_sif} \
        bash {params.scriptsdir}/ccbr_atac_motif_enrichment.bash \
        --narrowpeak $narrowpeak_fn \
        --genome {params.genome} \
        --threads {threads}
        rsync -Laz --progress /lscratch/$SLURM_JOBID/motif_enrichment ${{filepath}}_motif_enrichment
        echo -ne "$sample\t$group\t${{filepath}}_motif_enrichment\n" >> {output}
        rm -rf /lscratch/$SLURM_JOBID/*
        cd {params.workdir}
    done < $i
done
"""    


rule atac_frip_init:
    input:
        unpack(get_frip_input)
    output:
        outfile=join(WORKDIR,"QC","frip","{grp}.init.frip")
    params:
        rname='pl:fripinit',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        grp="{grp}"
    shell:"""
set -e -x -o pipefail
nsamples=0
for i in {input};do
    if [ $(basename $i) != "dummy" ]; then
    nsamples=$(echo "$nsamples+1"|bc)
    fi
done
nsamples=$(echo "$nsamples/5"|bc)

# the groups has only one samples ..  no replicates
if [ "$nsamples" -eq "1" ]; then
    echo "NO REPLICATES...THIS CODE IS NOT WRITTEN YET!!!"
    exit 1
fi

if [ "$nsamples" -gt "4" ]; then
    echo "5 or more REPLICATES...THIS CODE IS NOT WRITTEN YET!!!"
    exit 1
fi

# the groups has 2 replicates
if [ "$nsamples" -ge "2" ]; then

    tagAlign=$(echo {input.tagAlign1}|awk -F"/" "{{print \$NF}}")
    samplename=$(echo $tagAlign|awk -F".tagAlign" "{{print \$1}}")
    while read s g n;do
        if [ "$samplename" == "$s" ];then
            echo -ne "$s\t{input.tagAlign1}\tmacs2\tindividual_replicate\t$n\n" >> {output.outfile}
        fi
    done < {input.macs_sample_peaks1}
    while read s g n;do
        echo -ne "$samplename\t{input.tagAlign1}\tmacs2\treplicates_consensus\t$n\n" >> {output.outfile}
    done < {input.macs_group_peaks1}
    while read s g n;do
        if [ "$samplename" == "$s" ];then
            echo -ne "$s\t{input.tagAlign1}\tgenrich\tindividual_replicate\t$n\n" >> {output.outfile}
        fi
    done < {input.genrich_sample_peaks1}
    while read s g n;do
        echo -ne "$samplename\t{input.tagAlign1}\tgenrich\treplicates_consensus\t$n\n" >> {output.outfile}
    done < {input.genrich_group_peaks1}

    tagAlign=$(echo {input.tagAlign2}|awk -F"/" "{{print \$NF}}")
    samplename=$(echo $tagAlign|awk -F".tagAlign" "{{print \$1}}")
    while read s g n;do
        if [ "$samplename" == "$s" ];then
            echo -ne "$s\t{input.tagAlign2}\tmacs2\tindividual_replicate\t$n\n" >> {output.outfile}
        fi
    done < {input.macs_sample_peaks2}
    while read s g n;do
        echo -ne "$samplename\t{input.tagAlign2}\tmacs2\treplicates_consensus\t$n\n" >> {output.outfile}
    done < {input.macs_group_peaks2}
    while read s g n;do
        if [ "$samplename" == "$s" ];then
            echo -ne "$s\t{input.tagAlign2}\tgenrich\tindividual_replicate\t$n\n" >> {output.outfile}
        fi
    done < {input.genrich_sample_peaks2}
    while read s g n;do
        echo -ne "$samplename\t{input.tagAlign2}\tgenrich\treplicates_consensus\t$n\n" >> {output.outfile}
    done < {input.genrich_group_peaks2}

    if [ "$nsamples" -ge "3" ];then
        tagAlign=$(echo {input.tagAlign3}|awk -F"/" "{{print \$NF}}")
        samplename=$(echo $tagAlign|awk -F".tagAlign" "{{print \$1}}")
        while read s g n;do
            if [ "$samplename" == "$s" ];then
                echo -ne "$s\t{input.tagAlign3}\tmacs2\tindividual_replicate\t$n\n" >> {output.outfile}
            fi
        done < {input.macs_sample_peaks3}
        while read s g n;do
            echo -ne "$samplename\t{input.tagAlign3}\tmacs2\treplicates_consensus\t$n\n" >> {output.outfile}
        done < {input.macs_group_peaks3}
        while read s g n;do
            if [ "$samplename" == "$s" ];then
                echo -ne "$s\t{input.tagAlign3}\tgenrich\tindividual_replicate\t$n\n" >> {output.outfile}
            fi
        done < {input.genrich_sample_peaks3}
        while read s g n;do
            echo -ne "$samplename\t{input.tagAlign3}\tgenrich\treplicates_consensus\t$n\n" >> {output.outfile}
        done < {input.genrich_group_peaks3}
    fi
fi


"""

rule atac_frip:
    input:
        join(WORKDIR,"QC","frip","{grp}.init.frip")
    output:
        join(WORKDIR,"QC","frip","{grp}.frip")
    params:
        rname='pl:frip',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        grp="{grp}"
    shell:"""
set -e -x -o pipefail
cd /lscratch/$SLURM_JOBID
module load singularity
while read samplename tagalign method sample_or_group narrowpeak;do
    tagalign_fn=$(echo $tagalign|awk -F"/" "{{print \$NF}}")
    samplename=$(echo $tagalign_fn|awk -F".tagAlign" "{{print \$1}}")
    narrowpeak_fn=$(echo $narrowpeak|awk -F"/" "{{print \$NF}}")
    echo $tagalign
    echo $narrowpeak
    rsync -Laz --progress $tagalign /lscratch/$SLURM_JOBID/
    rsync -Laz --progress $narrowpeak /lscratch/$SLURM_JOBID/
    ls -alrth
    PYTHONNOUSERSITE=1 singularity exec --cleanenv \
    -B {params.workdir}/:/data2/,{params.index_dir}/:/index,/lscratch/$SLURM_JOBID:/data3/ \
    {params.singularity_sif} \
    bash {params.scriptsdir}/ccbr_frip.bash \
    --narrowpeak /data3/$narrowpeak_fn \
    --tagalign /data3/$tagalign_fn \
    --samplename $samplename \
    --genome {params.genome} \
    --out /data3/${{samplename}}.${{method}}.${{sample_or_group}}.frip
    frip=$(cat ${{samplename}}.${{method}}.${{sample_or_group}}.frip|grep -m1 "FRiP"|awk "{{print \$NF}}")
    fridhs=$(cat ${{samplename}}.${{method}}.${{sample_or_group}}.frip|grep -m1 "FRiDHS"|awk "{{print \$NF}}")
    fripro=$(cat ${{samplename}}.${{method}}.${{sample_or_group}}.frip|grep -m1 "FRiPro"|awk "{{print \$NF}}")
    frienh=$(cat ${{samplename}}.${{method}}.${{sample_or_group}}.frip|grep -m1 "FRiEnh"|awk "{{print \$NF}}")
    echo -ne "${{samplename}}\t${{tagalign}}\t${{method}}\t${{sample_or_group}}\t${{narrowpeak}}\t${{frip}}\t${{fridhs}}\t${{fripro}}\t${{frienh}}\n" >> {params.grp}.frip
done < {input}
rsync -Laz --progress {params.grp}.frip {output}
"""        


rule atac_fqscreen:
    input:
        expand(join(WORKDIR,"trim","{sample}.{r1r2}.trim.fastq.gz"),sample=SAMPLES,r1r2=["R1","R2"])
    params:
        rname='pl:fqscreen',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC")
    threads: 56
    output:
        fqscreendir=join(WORKDIR,"QC","FQscreen")
shell:"""
set -e -x -o pipefail

for f in {input};do
    rsync -Laz --progress $f /lscratch/$SLURM_JOBID/
done
cd /lscratch/$SLURM_JOBID

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_fqscreen.bash \
--thread {threads} \
--infastq $(ls *fastq.gz) \
--confurl https://hpc.nih.gov/~CCBR/sbg_reference_bundle/fastq_screen.conf \
--dblisturl https://hpc.nih.gov/~CCBR/sbg_reference_bundle/fastq_screen_databases_list.txt

for f in $(ls *_screen.*);do
rsync -az --progress $f {output.fqscreendir}/
done

"""


rule atac_macs_peakcalling:
    input:
        unpack(get_sample_tagAlignFiles_for_group)
    params:
        rname='pl:macs',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        grp="{grp}"
    output:
        groupPeakFileList=join(WORKDIR,"peaks","macs2","{grp}.group.macs2.peakfiles"),
        samplePeakFileList=join(WORKDIR,"peaks","macs2","{grp}.sample.macs2.peakfiles"),
        tn5knicksFileList=join(WORKDIR,"peaks","macs2","{grp}.macs2.tn5knicksbedfiles")
    shell:"""
set -e -x -o pipefail
module load singularity
for f in {input};do
rsync -Laz --progress $f /lscratch/$SLURM_JOBID/
done
if [ ! -d {params.workdir}/peaks/macs2/bigwig ];then mkdir -p {params.workdir}/peaks/macs2/bigwig;fi
if [ ! -d {params.workdir}/peaks/macs2/tn5knicks ];then mkdir -p {params.workdir}/peaks/macs2/tn5knicks;fi
if [ ! -d {params.qcdir}/peak_annotation ];then mkdir -p {params.qcdir}/peak_annotation;fi

cd /lscratch/$SLURM_JOBID

nsamples=0
for i in {input};do
    if [ $(basename $i) != "dummy" ];then
    nsamples=$(echo "$nsamples+1"|bc)
    fi
done
nsamples=$(echo "$nsamples/2"|bc)
rep1name=$(echo {input.tagAlign1}|awk -F"/" "{{print \$NF}}"|awk -F".tagAlign" "{{print \$1}}")
tagAlign1=$(echo {input.tagAlign1}|awk -F"/" "{{print \$NF}}")
genomefile1=$(echo {input.genomefile1}|awk -F"/" "{{print \$NF}}")

# the groups has only one samples ..  no replicates
if [ "$nsamples" -eq "1" ];then
echo "NO REPLICATES...THIS CODE IS NOT WRITTEN YET!!!"
exit 1
fi

# the groups has 2 replicates
if [ "$nsamples" -eq "2" ];then
rep2name=$(echo {input.tagAlign2}|awk -F"/" "{{print \$NF}}"|awk -F".tagAlign" "{{print \$1}}")
tagAlign2=$(echo {input.tagAlign2}|awk -F"/" "{{print \$NF}}")

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/,{params.index_dir}/:/index \
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_macs2_peak_calling.bash \
--tagalign1 $tagAlign1 \
--tagalign2 $tagAlign2 \
--rep1name $rep1name \
--rep2name $rep2name \
--samplename {params.grp} \
--genomefile $genomefile1 \
--genome {params.genome} \
--consensusbedfile {params.grp}.macs2.consensus.bed \
--scriptsfolder {params.scriptsdir}

ls -larth

for f in $rep1name $rep2name {params.grp}
do
    rsync -az --progress ${{f}}.macs2_summits.bed {params.workdir}/peaks/macs2/
    rsync -az --progress ${{f}}.macs2.narrowPeak {params.workdir}/peaks/macs2/
    rsync -az --progress ${{f}}.macs2.qfilter.narrowPeak {params.workdir}/peaks/macs2/

    for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
    do
        rsync -az --progress ${{f}}.macs2.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
        rsync -az --progress ${{f}}.macs2.qfilter.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
    done

    rsync -az --progress ${{f}}.macs2.bw {params.workdir}/peaks/macs2/bigwig/
    rsync -az --progress ${{f}}.macs2.tn5knicks.bed.gz {params.workdir}/peaks/macs2/tn5knicks/
    rsync -az --progress ${{f}}.macs2.tn5knicks.bam {params.workdir}/peaks/macs2/tn5knicks/
    rsync -az --progress ${{f}}.macs2.tn5knicks.bam.bai {params.workdir}/peaks/macs2/tn5knicks/

    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/macs2/tn5knicks/${{f}}.macs2.tn5knicks.bam\n" >> {output.tn5knicksFileList}
done

for f in $rep1name $rep2name
do
    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/macs2/${{f}}.macs2.narrowPeak\n" >> {output.samplePeakFileList}
done
fi

# the groups has 3 replicates
if [ "$nsamples" -eq "3" ];then
rep2name=$(echo {input.tagAlign2}|awk -F"/" "{{print \$NF}}"|awk -F".tagAlign" "{{print \$1}}")
tagAlign2=$(echo {input.tagAlign2}|awk -F"/" "{{print \$NF}}")
rep3name=$(echo {input.tagAlign3}|awk -F"/" "{{print \$NF}}"|awk -F".tagAlign" "{{print \$1}}")
tagAlign3=$(echo {input.tagAlign3}|awk -F"/" "{{print \$NF}}")

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/,{params.index_dir}/:/index \
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_macs2_peak_calling.bash \
--tagalign1 $tagAlign1 \
--tagalign2 $tagAlign2 \
--tagalign3 $tagAlign3 \
--rep1name $rep1name \
--rep2name $rep2name \
--rep3name $rep3name \
--samplename {params.grp} \
--genome {params.genome} \
--genomefile $genomefile1 \
--consensusbedfile {params.grp}.macs2.consensus.bed \
--scriptsfolder {params.scriptsdir}

ls -larth

for f in $rep1name $rep2name $rep3name {params.grp}
do
    rsync -az --progress ${{f}}.macs2_summits.bed {params.workdir}/peaks/macs2/
    rsync -az --progress ${{f}}.macs2.narrowPeak {params.workdir}/peaks/macs2/
    rsync -az --progress ${{f}}.macs2.qfilter.narrowPeak {params.workdir}/peaks/macs2/

    for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
    do
        rsync -az --progress ${{f}}.macs2.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
        rsync -az --progress ${{f}}.macs2.qfilter.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
    done

    rsync -az --progress ${{f}}.macs2.bw {params.workdir}/peaks/macs2/bigwig/
    rsync -az --progress ${{f}}.macs2.tn5knicks.bed.gz {params.workdir}/peaks/macs2/tn5knicks/
    rsync -az --progress ${{f}}.macs2.tn5knicks.bam {params.workdir}/peaks/macs2/tn5knicks/
    rsync -az --progress ${{f}}.macs2.tn5knicks.bam.bai {params.workdir}/peaks/macs2/tn5knicks/


    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/macs2/tn5knicks/${{f}}.macs2.tn5knicks.bam\n" >> {output.tn5knicksFileList}
done

for f in $rep1name $rep2name $rep3name
do
    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/macs2/${{f}}.macs2.narrowPeak\n" >> {output.samplePeakFileList}
done

fi

# for all replicates there will be consensus peaks
if [ "$nsamples" -ge "2" ];then

rsync -az --progress {params.grp}.macs2.consensus.bed {params.workdir}/peaks/macs2/
for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
do
    rsync -az --progress {params.grp}.macs2.consensus.bed.${{g}} {params.qcdir}/peak_annotation/
done

echo -ne "{params.grp}\t{params.grp}\t{params.workdir}/peaks/macs2/{params.grp}.macs2.consensus.bed\n" > {output.groupPeakFileList}

else
 touch {output.groupPeakFileList}
fi

# the groups as more than 3 replicates ... rare
if [ "$nsamples" -gt "3" ];then
echo "THIS CODE IS NOT WRITTEN YET!!!"
exit 1
fi

"""


rule atac_fld:
    input:
        dedupbam=join(WORKDIR,"bam","{sample}.dedup.bam")
    output:
        fld=join(WORKDIR,"QC","fld","{sample}.fld.txt")
    params:
        rname='pl:fld',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        sample="{sample}"
    shell:"""
set -e -x -o pipefail
module load singularity
for f in {input};do
rsync -Laz --progress $f /lscratch/$SLURM_JOBID/
done
bamfile=$(echo "{input.dedupbam}"|awk -F"/" "{{print \$NF}}")
fldfile=$(echo "{output.fld}"|awk -F"/" "{{print \$NF}}")
cd /lscratch/$SLURM_JOBID

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/ \
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_bam2FLD.bash \
--dedupbam $bamfile \
--fldout $fldfile \
--scriptsfolder {params.scriptsdir}

rsync -az --progress $fldfile {output.fld}

"""


rule atac_genrich_peakcalling:
    input:
        unpack(get_sample_qsortedbamfiles_for_group)
    params:
        rname='pl:genrich',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC"),
        grp="{grp}"
    output:
        groupPeakFileList=join(WORKDIR,"peaks","genrich","{grp}.group.genrich.peakfiles"),
        samplePeakFileList=join(WORKDIR,"peaks","genrich","{grp}.sample.genrich.peakfiles"),
        tn5knicksFileList=join(WORKDIR,"peaks","genrich","{grp}.genrich.tn5knicksbedfiles")
    shell:"""
set -e -x -o pipefail
module load singularity
for f in {input};do
rsync -Laz --progress $f /lscratch/$SLURM_JOBID/
done
cd /lscratch/$SLURM_JOBID

nsamples=0
for i in {input};do
    if [ $(basename $i) != "dummy" ]; then
    nsamples=$(echo "$nsamples+1"|bc)
    fi
done

if [ ! -d {params.workdir}/peaks/genrich/bigwig ];then mkdir -p {params.workdir}/peaks/genrich/bigwig;fi
if [ ! -d {params.workdir}/peaks/genrich/tn5knicks ];then mkdir -p {params.workdir}/peaks/genrich/tn5knicks;fi
if [ ! -d {params.qcdir}/peak_annotation ];then mkdir -p {params.qcdir}/peak_annotation;fi

rep1name=$(echo {input.qsortedbam1}|awk -F"/" "{{print \$NF}}"|awk -F".qsorted" "{{print \$1}}")
qsortedbam1=$(echo {input.qsortedbam1}|awk -F"/" "{{print \$NF}}")
peakfile1="${{rep1name}}.genrich.narrowPeak"

# the groups has only one samples ..  no replicates
if [ "$nsamples" -eq "1" ];then
echo "NO REPLICATES...THIS CODE IS NOT WRITTEN YET!!!"
exit 1
fi

# the groups as more than 2 replicates ... rare
if [ "$nsamples" -gt "3" ];then
echo "THIS CODE IS NOT WRITTEN YET!!!"
exit 1
fi

# the groups has 2 replicates
if [ "$nsamples" -eq "2" ];then
rep2name=$(echo {input.qsortedbam2}|awk -F"/" "{{print \$NF}}"|awk -F".qsorted" "{{print \$1}}")
qsortedbam2=$(echo {input.qsortedbam2}|awk -F"/" "{{print \$NF}}")
peakfile2="${{rep2name}}.genrich.narrowPeak"

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/,{params.index_dir}/:/index \
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_genrich_peak_calling.bash \
--bamrep1 $qsortedbam1 \
--bamrep2 $qsortedbam2 \
--peakfile1 $peakfile1 \
--peakfile2 $peakfile2 \
--genome {params.genome} \
--pooledpeakfile {params.grp}.genrich.narrowPeak \
--consensusbedfile {params.grp}.genrich.consensus.bed \
--scriptsfolder {params.scriptsdir}

ls -alrth


for f in $rep1name $rep2name {params.grp}
do
    rsync -az --progress ${{f}}.genrich.narrowPeak {params.workdir}/peaks/genrich/
    rsync -az --progress ${{f}}.genrich.qfilter.narrowPeak {params.workdir}/peaks/genrich/

    for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
    do
        rsync -az --progress ${{f}}.genrich.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
        rsync -az --progress ${{f}}.genrich.qfilter.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
    done

    rsync -az --progress ${{f}}.genrich.reads.bw {params.workdir}/peaks/genrich/bigwig/
    rsync -az --progress ${{f}}.genrich.tn5knicks.bed.gz {params.workdir}/peaks/genrich/tn5knicks/
    rsync -az --progress ${{f}}.genrich.tn5knicks.bam {params.workdir}/peaks/genrich/tn5knicks/
    rsync -az --progress ${{f}}.genrich.tn5knicks.bam.bai {params.workdir}/peaks/genrich/tn5knicks/

    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/genrich/tn5knicks/${{f}}.genrich.tn5knicks.bam\n" >> {output.tn5knicksFileList}
done

for f in $rep1name $rep2name
do
    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/genrich/${{f}}.genrich.narrowPeak\n" >> {output.samplePeakFileList}
done

fi

# the groups has 3 replicates
if [ "$nsamples" -eq "3" ];then
rep2name=$(echo {input.qsortedbam2}|awk -F"/" "{{print \$NF}}"|awk -F".qsorted" "{{print \$1}}")
qsortedbam2=$(echo {input.qsortedbam2}|awk -F"/" "{{print \$NF}}")
peakfile2="${{rep2name}}.genrich.narrowPeak"

rep3name=$(echo {input.qsortedbam3}|awk -F"/" "{{print \$NF}}"|awk -F".qsorted" "{{print \$1}}")
qsortedbam3=$(echo {input.qsortedbam3}|awk -F"/" "{{print \$NF}}")
peakfile3="${{rep3name}}.genrich.narrowPeak"

PYTHONNOUSERSITE=1 singularity exec --cleanenv \
-B {params.workdir}/:/data2/,{params.index_dir}/:/index \
{params.singularity_sif} \
bash {params.scriptsdir}/ccbr_genrich_peak_calling.bash \
--bamrep1 $qsortedbam1 \
--bamrep2 $qsortedbam2 \
--bamrep3 $qsortedbam3 \
--peakfile1 $peakfile1 \
--peakfile2 $peakfile2 \
--peakfile3 $peakfile3 \
--genome {params.genome} \
--pooledpeakfile {params.grp}.genrich.narrowPeak \
--consensusbedfile {params.grp}.genrich.consensus.bed \
--scriptsfolder {params.scriptsdir}

ls -alrth

for f in $rep1name $rep2name $rep3name {params.grp}
do
    rsync -az --progress ${{f}}.genrich.narrowPeak {params.workdir}/peaks/genrich/
    rsync -az --progress ${{f}}.genrich.qfilter.narrowPeak {params.workdir}/peaks/genrich/

    for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
    do
        rsync -az --progress ${{f}}.genrich.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
        rsync -az --progress ${{f}}.genrich.qfilter.narrowPeak.${{g}} {params.qcdir}/peak_annotation/
    done

    rsync -az --progress ${{f}}.genrich.reads.bw {params.workdir}/peaks/genrich/bigwig/
    rsync -az --progress ${{f}}.genrich.tn5knicks.bed.gz {params.workdir}/peaks/genrich/tn5knicks/
    rsync -az --progress ${{f}}.genrich.tn5knicks.bam {params.workdir}/peaks/genrich/tn5knicks/
    rsync -az --progress ${{f}}.genrich.tn5knicks.bam.bai {params.workdir}/peaks/genrich/tn5knicks/

    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/genrich/tn5knicks/${{f}}.genrich.tn5knicks.bam\n" >> {output.tn5knicksFileList}
done

for f in $rep1name $rep2name $rep3name
do
    echo -ne "${{f}}\t{params.grp}\t{params.workdir}/peaks/genrich/${{f}}.genrich.narrowPeak\n" >> {output.samplePeakFileList}
done

fi

# if you have replicates .. you will have consensus peaks
if [ "$nsamples" -ge "2" ];then

rsync -az --progress {params.grp}.genrich.consensus.bed {params.workdir}/peaks/genrich/
for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
do
    rsync -az --progress {params.grp}.genrich.consensus.bed.${{g}} {params.qcdir}/peak_annotation/
done

echo -ne "{params.grp}\t{params.grp}\t{params.workdir}/peaks/genrich/{params.grp}.genrich.consensus.bed\n" > {output.groupPeakFileList}

else
touch {output.groupPeakFileList}
fi
"""


rule jaccard:
    input:
        expand(join(WORKDIR,"peaks","macs2","{grp}.group.macs2.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","macs2","{grp}.sample.macs2.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.group.genrich.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.sample.genrich.peakfiles"),grp=GROUPS)
    output:
        macs2persamplejaccardpca=join(WORKDIR,"QC","jaccard","macs2.sample.jaccard.pca.html"),
        macs2pergroupjaccardpca=join(WORKDIR,"QC","jaccard","macs2.group.jaccard.pca.html"),
        macs2persamplegroupjaccardpca=join(WORKDIR,"QC","jaccard","macs2.sample_group.jaccard.pca.html"),
        genrichpersamplejaccardpca=join(WORKDIR,"QC","jaccard","genrich.sample.jaccard.pca.html"),
        genrichpergroupjaccardpca=join(WORKDIR,"QC","jaccard","genrich.group.jaccard.pca.html"),
        genrichpersamplegroupjaccardpca=join(WORKDIR,"QC","jaccard","genrich.sample_group.jaccard.pca.html"),
        allmethodspersamplejaccardpca=join(WORKDIR,"QC","jaccard","allmethods.sample.jaccard.pca.html"),
        allmethodspergroupjaccardpca=join(WORKDIR,"QC","jaccard","allmethods.group.jaccard.pca.html"),
        allmethodspersamplegroupjaccardpca=join(WORKDIR,"QC","jaccard","allmethods.sample_group.jaccard.pca.html")
    params:
        rname='pl:jaccard',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC")
    shell:"""
set -e -x -o pipefail
module load singularity
for f in {input};do
    rsync -Laz --progress $f /lscratch/$SLURM_JOBID/
    while read label glabel file;do
        rsync -Laz --progress $file /lscratch/$SLURM_JOBID/
    done < $f
done
cd /lscratch/$SLURM_JOBID

for f in $(ls *narrowPeak);do
        PYTHONNOUSERSITE=1 singularity exec --cleanenv \
        -B {params.workdir}/:/data2/,{params.index_dir}/:/index \
        {params.singularity_sif} bedSort $f $f
done


awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *group.macs2.peakfiles > macs2.group.peakfiles.tmp
awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *sample.macs2.peakfiles > macs2.sample.peakfiles.tmp
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "1000" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < macs2.sample.peakfiles.tmp > macs2.sample.peakfiles
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "1000" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < macs2.group.peakfiles.tmp > macs2.group.peakfiles
cat macs2.group.peakfiles macs2.sample.peakfiles > macs2.sample_group.peakfiles
awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *group.genrich.peakfiles > genrich.group.peakfiles.tmp
awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *sample.genrich.peakfiles > genrich.sample.peakfiles.tmp
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "1000" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < genrich.sample.peakfiles.tmp > genrich.sample.peakfiles
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "1000" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < genrich.group.peakfiles.tmp > genrich.group.peakfiles
cat genrich.group.peakfiles genrich.sample.peakfiles > genrich.sample_group.peakfiles

for m in "macs2" "genrich";do
    while read sample group file;do
        echo -ne "${{sample}}_${{m}}\t${{group}}_${{m}}\t${{file}}\n"
    done < ${{m}}.sample.peakfiles
done > allmethods.sample.peakfiles
for m in "macs2" "genrich";do
    while read sample group file;do
        echo -ne "${{sample}}_${{m}}\t${{group}}_${{m}}\t${{file}}\n"
    done < ${{m}}.group.peakfiles
done > allmethods.group.peakfiles
for m in "macs2" "genrich";do
    while read sample group file;do
        echo -ne "${{sample}}_${{m}}\t${{group}}_${{m}}\t${{file}}\n"
    done < ${{m}}.sample_group.peakfiles
done > allmethods.sample_group.peakfiles

for m in "macs2" "genrich" "allmethods";do
    for f in "group" "sample" "sample_group";do
        PYTHONNOUSERSITE=1 singularity exec --cleanenv \
        -B {params.workdir}/:/data2/,{params.index_dir}/:/index \
        {params.singularity_sif} \
        bash {params.scriptsdir}/ccbr_jaccard_pca.bash \
        --inputfilelist ${{m}}.${{f}}.peakfiles \
        --pairwise ${{m}}.${{f}}.jaccard.pairwise.txt \
        --pcahtml ${{m}}.${{f}}.jaccard.pca.html \
        --scriptsfolder {params.scriptsdir}
        rsync -az --progress ${{m}}.${{f}}.jaccard.pairwise.txt {params.qcdir}/jaccard/
        rsync -az --progress ${{m}}.${{f}}.jaccard.pca.html {params.qcdir}/jaccard/
    done
done



"""


rule atac_multiqc:
    input:
        expand(join(WORKDIR,"QC","fastqc","{sample}.R1_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R2_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","{sample}.nreads.txt"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R1.noBL_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","fastqc","{sample}.R2.noBL_fastqc.zip"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","{sample}.dupmetric"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","preseq","{sample}.nrf"),sample=SAMPLES),
        expand(join(WORKDIR,"peaks","macs2","{grp}.group.macs2.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","macs2","{grp}.sample.macs2.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.group.genrich.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"peaks","genrich","{grp}.sample.genrich.peakfiles"),grp=GROUPS),
        expand(join(WORKDIR,"QC","fld","{sample}.fld.txt"),sample=SAMPLES),
        expand(join(WORKDIR,"QC","tss","{sample}.tss.txt"),sample=SAMPLES),
        join(WORKDIR,"QC","jaccard","macs2.sample.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","macs2.group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","macs2.sample_group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","genrich.sample.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","genrich.group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","genrich.sample_group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","allmethods.sample.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","allmethods.group.jaccard.pca.html"),
        join(WORKDIR,"QC","jaccard","allmethods.sample_group.jaccard.pca.html"),
        expand(join(WORKDIR,"QC","frip","{grp}.frip"),grp=GROUPS)
    output:
        multiqchtml=join(WORKDIR,"QC","multiqc_report.html"),
        qcstatstable=join(WORKDIR,"QC","QCStats.txt")
    params:
        rname='pl:qc',
        genome=GENOME,
        singularity_sif=SINGULARITY_SIF,
        index_dir=INDEX_DIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTS_DIR,
        qcdir=join(WORKDIR,"QC")
    shell:"""
set -e -x -o pipefail
module load singularity
cd {params.qcdir}
PYTHONNOUSERSITE=1 singularity exec --cleanenv \
 -B {params.workdir}/:/data2/,{params.index_dir}/:/index \
 {params.singularity_sif} \
 bash {params.scriptsdir}/ccbr_atac_qc.bash \
 --qcfolder /data2/QC \
 --scriptsfolder {params.scriptsdir}
"""