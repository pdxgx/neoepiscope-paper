#!/usr/bin/bash

## Benchmarking of BWA for neoepiscope (BWA v0.7.9a-r786)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_bwa.sh {HOME_DIR} {SAMPLE_ID} {FASTQ1} {FASTQ2} {REFERENCE} {BWA} {SAMTOOLS}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## SAMPLE_ID is your tumor sample ID or normal sample ID (must be run for both)
## FASTQ1 is your first fastq file; FASTQ2 is your second
## REFERENCE is the path to your hg38 reference fasta
## BWA is the path to your BWA executable
## SAMTOOLS is the path to your samtools 1.3.1 executable

HOME_DIR=$1
SAMPLE_ID=$2
FASTQ1=$3
FASTQ2=$4
BAM=${HOME_DIR}/${SAMPLE_ID}.bam
REFERENCE=$5
TIMELOG=${HOME_DIR}/${SAMPLE_ID}.bwa.time_log
CPUINFO=${HOME_DIR}/${SAMPLE_ID}.bwa.cpu_data
BWA=$6
SAMTOOLS=$7

cat /proc/cpuinfo > $CPUINFO

echo 'BWA RUN TIME' > $TIMELOG

time ($BWA mem -t 4 $REFERENCE $FASTQ1 $FASTQ2 | $SAMTOOLS sort -@4 -O BAM -o $BAM -) 2>> $TIMELOG
