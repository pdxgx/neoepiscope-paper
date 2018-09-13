#!/usr/bin/bash

## Benchmarking of MuTect 2 for neoepiscope (GATK v3.5)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_filter_mutect.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.filtermutect.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.filtermutect.cpu_data
RAW_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.raw.vcf
FILTERED_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf

cat /proc/cpuinfo > $CPUINFO

echo 'FILTER MUTECT RUN TIME' > $TIMELOG

time (grep '^#' $RAW_VCF > $FILTERED_VCF && grep -v '^#' $RAW_VCF | awk '$7 == "PASS"' \
		>> $FILTERED_VCF 2>&1) 2>> $TIMELOG
