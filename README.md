# cutandrun_multimapper
Snakemake Pipeline for Multi-mapper CUT&RUN Analysis

## Overview

This Snakemake pipeline aligns multi-mapping reads by randomly assigning them to one of their exact-matching sites in the reference genome. By default, reads mapping to up to alternative 5,000 locations are retained, with no mismatches allowed. These parameters can be customized as needed.

## Required software
- Snakemake
- Python 3
- Perl
- Bash
- FastQC
- Trim Galore
- Cutadapt
- BWA
- SAMtools
- Picard
- BEDOPS
- UCSC bedGraphToBigWig

## Required files in a working driectory
- A folder including paired-end FASTQ files (e.g., fastq/)
- Sample sheet (samples.csv)
- Configuration file (config.yaml)
- take_one_random_read.pl

## Required reference files (entries in config.yaml)
- Reference: Reference genome FASTA indexed with BWA
- sample: Sample sheet (samples.csv)
- chromsize: Chromosome sizes file including two columns (chromosome and length (in bp).

## Run
cd /my/working/directory

snakemake --cores <number of cores>

## Workflow

This Snakemake pipeline processes paired-end CUT&RUN sequencing data from raw FASTQ files to normalized BigWig tracks. The workflow consists of read quality control, adapter trimming, alignment, duplicate removal, multi-mapping read selection, read extension, genomic binning, signal quantification, normalization, smoothing, and BigWig generation.

### 1. Read quality control

Raw paired-end FASTQ files are first assessed using **FastQC** to evaluate sequencing quality and identify potential issues with the input data.

### 2. Adapter trimming

Sequencing adapters and low-quality terminal bases are removed using **Trim Galore** with paired-end mode. FastQC is subsequently run on the trimmed reads to assess the post-trimming read quality.

### 3. Read alignment

Trimmed reads are aligned to the reference genome using **BWA aln**. The alignment is performed with no mismatches allowed (`-k 0`), and up to 5,000 alternative alignments are retained for each read pair.

Only properly paired alignments are retained for downstream processing.

### 4. Read group assignment and duplicate removal

Read groups are added using **Picard AddOrReplaceReadGroups**. PCR or optical duplicates are subsequently identified and removed using **Picard MarkDuplicates**.

### 5. Multi-mapping read selection

For reads with multiple equally best alignment sites, one alignment site is randomly selected using the custom `take_one_random_read.pl` script. This step assigns each multi-mapping read to a single genomic location while retaining reads that have multiple exact or equally scoring alignment sites.

The resulting BED files are sorted using **BEDOPS `sort-bed`**.

### 6. Strand-aware read extension

Selected reads are extended to 200 bp in a strand-aware manner to represent the estimated fragment span. For reads on the `+` strand, the interval is extended downstream from the read start, whereas reads on the `−` strand are extended upstream from the read end.

Extended intervals are clipped to chromosome boundaries.

### 7. Genomic bin generation

The reference genome is divided into fixed, non-overlapping **30-bp genomic bins** using BEDOPS. These bins serve as the common genomic framework for signal quantification.

### 8. Signal quantification

The number of extended reads overlapping each 30-bp genomic bin is calculated using **BEDOPS `bedmap --count`**.

This produces a genome-wide bedGraph containing the raw read count for each 30-bp bin.

### 9. CPM normalization

Raw bin-level read counts are normalized to **counts per million (CPM)** using the total number of selected reads for each sample.

### 10. Signal smoothing

The CPM-normalized signal is smoothed using a **300-bp window**. By default, a mean-based smoothing approach is used, with an option to use Gaussian smoothing.

### 11. BigWig generation

The smoothed bedGraph files are converted to indexed **BigWig** files using the UCSC `bedGraphToBigWig` utility.

The resulting BigWig tracks can be used for genome-browser visualization and downstream CUT&RUN signal analysis.


#### Citation

This pipeline was developed and used for the analysis presented in:

Yang, Yoon, Jaiswal et al. Recruitment of TEX15 by HUSH2 drives repression of multicopy gene families in fertility regulation. _Nature Communications_. Accepted for publication.

#### - Contact: Sora Yoon, Ph.D. <sorayoon@upenn.edu>
