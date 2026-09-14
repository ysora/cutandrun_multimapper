import os
import re
import pandas as pd
from snakemake.utils import validate, min_version

shell.prefix("set -eo pipefail; ")

configfile: "config.yaml"

localrules: all

samples_df = pd.read_csv(config['samples']).set_index("cond", drop=False)

samples_dict = samples_df.to_dict(orient='index')
conditions = [samples_dict[x]['cond'] for x in samples_df.index]

ALL_SAMPLES = samples_df.index.to_list()
ALL_FASTQC  = expand(["01.fastqc/{sample}_R1_001_fastqc.html", "01.fastqc/{sample}_R2_001_fastqc.html"], sample = ALL_SAMPLES)
ALL_TRIM_ADAPTER_FASTQC = expand(["01.fastqc_trim_galore/{sample}_R1_001_val_1_fastqc.html", "01.fastqc_trim_galore/{sample}_R2_001_val_2_fastqc.html"], sample = ALL_SAMPLES)
ALL_TRIM_FASTQ = expand(["01.trimmed/{sample}_R1_001_val_1.fq.gz", "01.trimmed/{sample}_R2_001_val_2.fq.gz"], sample = ALL_SAMPLES)
ALL_BAM = expand("02.bam/{sample}.sorted.bam", sample=ALL_SAMPLES)
ALL_BAM_INDEX = expand("02.bam/{sample}.sorted.bam.bai", sample=ALL_SAMPLES)
ALL_SELECT = expand("03.record_select/{sample}.one_random_in_multi_best_hits.bed.sorted", sample=ALL_SAMPLES)
ALL_EXT = expand("03.record_select/{sample}.one_random_in_multi_best_hits.ext200.bed.sorted", sample=ALL_SAMPLES)
ALL_WINDOW = expand("04.analysis_window/{sample}.bins30.sorted.bed", sample=ALL_SAMPLES)
ALL_BEDGRAPH=expand("05.bedGraph/{sample}.bin30.ext200.count.bedGraph", sample=ALL_SAMPLES)
ALL_BEDGRAPH_CPM=expand("05.bedGraph/{sample}.bin30.ext200.count.cpm.bedGraph", sample=ALL_SAMPLES)
ALL_BEDGRAPH_CPM_S=expand("05.bedGraph/{sample}.bin30.ext200.count.cpm.smooth300.bedGraph", sample=ALL_SAMPLES)
ALL_BigWig_CPM_S = expand("06.bigwig_sliding/{sample}.cpm.bin30.ext200.smooth300.bw", sample=ALL_SAMPLES)

rule all:
    input: ALL_FASTQC + ALL_TRIM_ADAPTER_FASTQC + ALL_TRIM_FASTQ + ALL_BAM + ALL_BAM_INDEX + ALL_SELECT + ALL_EXT + ALL_WINDOW + ALL_BEDGRAPH + ALL_BEDGRAPH_CPM + ALL_BEDGRAPH_CPM_S + ALL_BigWig_CPM_S
 
rule fastqc:
    input:  fq1 = lambda wildcards: samples_dict[wildcards.sample]['fq1'],
            fq2 = lambda wildcards: samples_dict[wildcards.sample]['fq2']
    output: "01.fastqc/{sample}_R1_001_fastqc.html",
            temp("01.fastqc/{sample}_R1_001_fastqc.zip"),
            "01.fastqc/{sample}_R2_001_fastqc.html",
            temp("01.fastqc/{sample}_R2_001_fastqc.zip")
    log:    "00.log/{sample}.fastqc"
    threads: samples_df.shape[0]*2
    resources:
        mem  = 2 ,
        time = 20
    message: "fastqc {input}: {threads} / {resources.mem}"
    shell:
        """
        echo $(python --version)
        fastqc -o 01.fastqc -f fastq --noextract {input.fq1} 2>{log}
        fastqc -o 01.fastqc -f fastq --noextract {input.fq2} 2>{log}
        """
        
rule trim_adapter:
    input:  fq1 = lambda wildcards: samples_dict[wildcards.sample]['fq1'],
            fq2 = lambda wildcards: samples_dict[wildcards.sample]['fq2']
    output: 
        fq1_trimmed = "01.trimmed/{sample}_R1_001_val_1.fq.gz",
        fq2_trimmed = "01.trimmed/{sample}_R2_001_val_2.fq.gz",
        fastqc_r1 = "01.fastqc_trim_galore/{sample}_R1_001_val_1_fastqc.html",
        fastqc_r2 = "01.fastqc_trim_galore/{sample}_R2_001_val_2_fastqc.html"
    log: "00.log/{sample}.trim_galore.log"
    threads: 4
    resources:
        mem  = 4,
        time = 30
    message: "Trimming adapters for {input.fq1} and {input.fq2}"
    shell:
        """
        /usr/bin/trim_galore --paired --gzip --fastqc --fastqc_args "--outdir 01.fastqc_trim_galore" \
        --output_dir 01.trimmed --length 20 --cores {threads} \
        {input.fq1} {input.fq2} &> {log}
        """

rule align:
    input:
        fq1_trimmed = "01.trimmed/{sample}_R1_001_val_1.fq.gz",
        fq2_trimmed = "01.trimmed/{sample}_R2_001_val_2.fq.gz"
    output: bam="02.bwa_{sample}/{sample}.sorted.bam",
            sai1=temp("02.bwa_{sample}/{sample}.R1.sai"),
            sai2=temp("02.bwa_{sample}/{sample}.R2.sai")
    threads: 20
    resources:
        mem    = 80,
        time   = 105
    message: "aligning {input}: {threads} threads / {resources.mem}"
    log: "00.log/{sample}.bwa"
    shell:
        """
        set -euo pipefail
        bwa aln -k 0 -t {threads} -M 100 -O 100 -E 100 -R 150 {config[Reference]} {input.fq1_trimmed} > {output.sai1}
        bwa aln -k 0 -t {threads} -M 100 -O 100 -E 100 -R 150 {config[Reference]} {input.fq2_trimmed} > {output.sai2}
        bwa sampe -n 5000 -N 0 {config[Reference]} {output.sai1} {output.sai2} {input.fq1_trimmed} {input.fq2_trimmed} \
        | samtools view -b -@ {threads} -F 4 -f 2 -h - \
        | samtools sort -m 2G -@ {threads} -O bam -T {output.bam}.tmp -o {output.bam}
        """                      

rule AddorReplace_read_groups:
    input: "02.bwa_{sample}/{sample}.sorted.bam"
    output: temp("02.bwa_{sample}/{sample}.RG.bam")
    log:    "00.log/{sample}.picard_RGgroup"
    message: "Add or replace read groups using picard : {input}"
    params:
        name = "02.bwa_{sample}"
    shell:
        """
        IFS='_' read -ra parts <<< "{wildcards.sample}"
        RGID="${{parts[1]}}"
        RGSM="${{parts[@]:2}}"
    RGSM=$(IFS=_; echo "${{parts[*]:2}}")
        set +o pipefail
        first_read=$(samtools view {input} | head -n 1 | cut -f1)
        set -o pipefail
        flowcell=$(echo "$first_read" | cut -d':' -f3)
        lane=$(echo "$first_read" | cut -d':' -f4)
        RGPU="${{flowcell}}.${{lane}}"
        RGPL="illumina"
        RGLB="{wildcards.sample}_lib1"

        picard AddOrReplaceReadGroups \
                I={input} \
                O={output} \
                RGID="$RGID" \
                RGLB="$RGLB" \
                RGPL="$RGPL" \
                RGPU="$RGPU" \
                RGSM="$RGSM"
    
        """
rule remove_duplicates:
    input:  "02.bwa_{sample}/{sample}.RG.bam"
    output: "02.bam/{sample}.sorted.bam",
            temp(directory("02.bam/temp/{sample}"))
    log:    "00.log/{sample}.picard"
    message: "Remove duplicates using picard: {input}"
    params:
        name = "02.bwa_{sample}"

    shell:
        """
        picard MarkDuplicates \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\
            METRICS_FILE=out.metrics \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true  \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR={output[1]} \
            INPUT={input} \
            OUTPUT={output[0]} 2> {log}
        if [ -d {params.name} ]
        then
            rm -rf {params.name}
        fi
        """

rule index_bam:
    input:  "02.bam/{sample}.sorted.bam"
    output: "02.bam/{sample}.sorted.bam.bai"
    log:    "00.log/{sample}.index_bam"
    threads: 1
    resources:
        mem   = 500,
        time  = 10
    message: "index_bam {input}: {threads} threads / {resources.mem}"
    shell:
        """
        module load samtools
        samtools index {input} 2> {log}
        """

rule picard_insert_size_histogram:
    input:  "02.bam/{sample}.sorted.bam"
    output: a="01.fastqc/{sample}_insertSizes.txt",
            b="01.fastqc/{sample}_insertSizes.pdf"
    log:    "00.log/{sample}.picard_insert_size_histogram"
    threads: 1
    message: "Picard insert size histogram for {input}"
    shell:
        """
         picard CollectInsertSizeMetrics \
         I={input} \
         O={output.a} \
         H={output.b} 2> {log}
        """
        
rule assign_to_one_random_best_hit_site:
    input: "02.bam/{sample}.sorted.bam"
    output: bed=temp("03.record_select/{sample}.one_random_in_multi_best_hits.bed")
    log:    "00.log/{sample}.score_adjustment"
    threads: 20
    message: "Score adjustment for {input}"
    params:
        perlcodepath = "./take_one_random_read.pl"
    shell:
        """
        samtools view {input} | perl {params.perlcodepath} > {output.bed} 2> {log}
        """

rule bed_sort:
    input:
        bed = "03.record_select/{sample}.one_random_in_multi_best_hits.bed"
    output:
        bedcopy=temp("03.record_select/{sample}.one_random_in_multi_best_hits.bed.sorted.copy.bed"),
        sorted = "03.record_select/{sample}.one_random_in_multi_best_hits.bed.sorted"
    log:
        "00.log/{sample}.sort_bed"
    threads: 20
    message: "Sort genomic coordinates using BEDOPS sort-bed for {input}"
    shell:
        r"""
        set -euo pipefail
        cp {input.bed} {output.bedcopy}
        sort-bed {output.bedcopy} > {output.sorted} 2> {log}
        """


rule extend_reads_to_200:
    input:
        bed = "03.record_select/{sample}.one_random_in_multi_best_hits.bed.sorted",
        chromsizes={config["chromsize"]}
    output:
        ext = "03.record_select/{sample}.one_random_in_multi_best_hits.ext200.bed.sorted"
    log:
        "00.log/{sample}.extend200"
    threads: 20
    message: "Reads were strand-aware expanded to 200 bp and clipped to chromosome boundaries. {input}"
    shell:
        r"""
        set -euo pipefail
        # chromsizes: chr \t length
        awk 'BEGIN{{OFS="\t"}}
             FNR==NR {{len[$1]=$2; next}}
             {{
               chr=$1; s=$2; e=$3; strand=$6;
               if (strand=="+")    {{ ns=s; ne=s+200; }}
               else if (strand=="-") {{ ns=e-200; ne=e; }}
               else                 {{ ns=s; ne=e; }}  

               if (ns<0) ns=0;
               if (ne>len[chr]) ne=len[chr];

               # 4th has 1-based values so that bedmap --count is available in next step.
               print chr, ns, ne, 1;
             }}' {input.chromsizes} {input.bed} \
          | sort-bed - > {output.ext} 2> {log}
        """


rule make_bins30:
    input:
        chromsizes={config["chromsize"]}
    output:
        bounds = "04.analysis_window/{sample}.mm.genome.bounds.bed",
        bins   = "04.analysis_window/{sample}.bins30.sorted.bed"
    log:
        "00.log/{sample}.make_bins30"
    params:
        BIN = 30
    message: "Create fixed windows of 30bp."
    shell:
        r"""
        set -euo pipefail
        # genome bounds (chr \t 0 \t length)
        awk 'BEGIN{{OFS="\t"}} {{print $1, 0, $2}}' {input.chromsizes} > {output.bounds}
        bedops --chop {params.BIN} {output.bounds} \
          | sort-bed - > {output.bins} 2> {log}
        """

rule bedGraph_bin30_count:
    input:
        frame="03.record_select/{sample}.one_random_in_multi_best_hits.ext200.bed.sorted",
        bins="04.analysis_window/{sample}.bins30.sorted.bed"
    output:
        bdg="05.bedGraph/{sample}.bin30.ext200.count.bedGraph"
    log:
        "00.log/{sample}.bin30_count"
    threads: 20
    message: "Create a bedgraph counting extended reads mapped to each 30bp bin: {input.frame}"
    shell:
        r"""
        set -euo pipefail
        bedmap --echo --count --delim "\t" --unmapped-val 0 \
          {input.bins} {input.frame} > {output.bdg} 2> {log}
        """


rule normalized_bedGraph_bin30_cpm:
    input:
        frame="03.record_select/{sample}.one_random_in_multi_best_hits.bed.sorted",
        bdg="05.bedGraph/{sample}.bin30.ext200.count.bedGraph"
    output:
        cpm="05.bedGraph/{sample}.bin30.ext200.count.cpm.bedGraph"
    log:
        "00.log/{sample}.bin30_cpm"
    threads: 20
    message: "CPM normalization on: {input.frame}"
    shell:
        r"""
        set -euo pipefail
        total=$(wc -l < {input.frame})
        awk -v total="${{total}}" 'BEGIN{{OFS="\t"}}
            {{
              if($3 > $2) {{
                cpm = ($4 / total) * 1000000;
                print $1, $2, $3, cpm;
              }}
            }}' {input.bdg} > {output.cpm} 2> {log}
        """


rule smooth_cpm_bedGraph_300:
    input:  "05.bedGraph/{sample}.bin30.ext200.count.cpm.bedGraph"
    output: "05.bedGraph/{sample}.bin30.ext200.count.cpm.smooth300.bedGraph"
    log:    "00.log/{sample}.smooth300"
    params:
        BIN=30,
        SMOOTH=300,
        MODE="mean"   # "mean" or "gaussian"
    message: "Smoothing on: {input}"
    shell:
        """
        set -euo pipefail
        python - << 'PY'
import sys, math, collections
bin_size = {params.BIN}
smooth_len = {params.SMOOTH}
mode = "{params.MODE}"
half_window_bins = int((smooth_len / bin_size) // 2)

by_chr = collections.OrderedDict()
with open("{input}", "r") as f:
    for line in f:
        if not line.strip():
            continue
        chrom, s, e, val = line.strip().split("\\t")
        s = int(s); e = int(e); v = float(val)
        by_chr.setdefault(chrom, []).append((s, e, v))

win = 2 * half_window_bins + 1
if mode == "gaussian":
    sigma = smooth_len / 6.0
    sigma_bins = max(sigma / bin_size, 1e-6)
    kernel = [math.exp(-0.5 * ((i - half_window_bins)**2) / (sigma_bins**2)) for i in range(win)]
else:
    kernel = [1.0] * win
norm = sum(kernel)
kernel = [k / norm for k in kernel]

with open("{output}", "w") as out:
    for chrom, arr in by_chr.items():
        vals = [v for _,_,v in arr]
        starts = [s for s,_,_ in arr]
        ends = [e for _,e,_ in arr]
        n = len(vals)
        for i in range(n):
            acc = 0.0; wsum = 0.0
            for k in range(win):
                j = i + k - half_window_bins
                if 0 <= j < n:
                    acc += kernel[k] * vals[j]
                    wsum += kernel[k]
            sm = acc / wsum if wsum > 0 else vals[i]
            out.write("%s\\t%d\\t%d\\t%f\\n" % (chrom, starts[i], ends[i], sm))
PY
        """

rule make_normalized_BigWig_bin30_smooth:
    input: "05.bedGraph/{sample}.bin30.ext200.count.cpm.smooth300.bedGraph"
    output: "06.bigwig_sliding/{sample}.cpm.bin30.ext200.smooth300.bw"
    log:    "00.log/{sample}.cpm.bin30.ext200.smooth300.BigWig"
    threads: 20
    message: "Convert to a BigWig format: {input}"
    shell:
        r"""
        set -euo pipefail
        sort -k1,1 -k2,2n {input} > {input}.sorted 2> {log}
        bedGraphToBigWig {input}.sorted {config[chromsize]} {output} 2>> {log}
        rm -f {input}.sorted
        """
