#!/usr/bin/bash

## Benchmarking of HLA typing for neoepiscope (Optitype v1.0)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_optitype.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {FASTQ1} {FASTQ2} {OPTITYPE} {CONFIG}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## FASTQ1 is your first fastq file; FASTQ2 is your second
## OPTITYPE is the path to your Optitype python script
## CONFIG is the path to your Optitype config file

HOME_DIR=$1
SAMPLE_ID=$2
FASTQ1=$3
FASTQ2=$4
OPTITYPE=$5
CONFIG=$6
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}.optitype.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}.optitype.cpu_data

#export PATH=/mnt/lustre1/CompBio/bin:$PATH
export HDF5_USE_FILE_LOCKING=FALSE

cat /proc/cpuinfo > $CPUINFO

echo 'OPTITYPE RUN TIME' > $TIMELOG

time (python $OPTITYPE --i $FASTQ1 $FASTQ2 --dna -e 1 -o $HOME_DIR --config $CONFIG \
		--prefix $SAMPLE_ID --verbose 2>&1) 2>> $TIMELOG
