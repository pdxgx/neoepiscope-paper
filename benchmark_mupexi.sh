#!/usr/bin/bash

## Benchmarking of neoepitope enumeration for neoepiscope (MuPeXI v1.2.0)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_mupexi.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {HLA_ALLELES} {MUPEXI}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## HLA_ALLELES is your comma-separated list of HLA alleles
## MUPEXI is the path to your MuPeXI python script

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
HLA_ALLELES=$4
MUPEXI=$5
VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mupexi.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mupexi.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'MUPEXI RUN TIME' > $TIMELOG

time (python $MUPEXI -v $VCF -a $HLA_ALLELES -l 9 -d $HOME_DIR --netmhc-full-anal \
		-p ${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID} -F 4 \
		-c ${HOME_DIR}MuPeXI/config.ini 2>&1) 2>> $TIMELOG
