#!/usr/bin/bash

## Benchmarking of neoepitope enumeration for neoepiscope (neoepiscope v0.1.0)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_neoepiscope.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {HLA_ALLELES}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## HLA_ALLELES is your comma-separated list of HLA alleles

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
HLA_ALLELES=$4
NEOEPISCOPE=${HOME_DIR}/neoepiscope/neoepiscope
HAPLOTYPES=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.haplotypes
PREPPED_HAPLOTYPES=${HAPLOTYPES}.prepped
VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.merged.vcf
OUTPUT=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.neoepiscope.out
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.neoepiscope.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.neoepiscope.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'NEOEPISCOPE RUN TIME' > $TIMELOG

time (neoepiscope prep -v $VCF -c $HAPLOTYPES -o $PREPPED_HAPLOTYPES \
		&& neoepiscope call -v $VCF -c $PREPPED_HAPLOTYPES -k 9 -b GRCh38 -o $OUTPUT \
    -p netMHCpan 4 affinity,rank -a $HLA_ALLELES \
	        2>&1) 2>> $TIMELOG || exit 1
