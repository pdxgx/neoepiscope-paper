#!/usr/bin/bash

## Benchmarking of neoepitope enumeration for neoepiscope
## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_tsnad.sh {HOME_DIR} {PATIENT} {TSNAD} {CONFIG}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## PATIENT is the patient ID
## TSNAD is the path to your TSNAD antigen_predicting_pipeline.py
## CONFIG is the path to your TSNAD config file for the patient

HOME_DIR=$1
PATIENT=$2
TUMOR_SAMPLE_ID=${2}_tumor
NORMAL_SAMPLE_ID=${2}_normal
TSNAD=$3
CONFIG=$4
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.tsnad.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.tsnad.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'TSNAD RUN TIME' > $TIMELOG

time (mkdir ${HOME_DIR}${PATIENT}_tsnad/ && python $TSNAD -c $CONFIG 2>&1) 2>> $TIMELOG
