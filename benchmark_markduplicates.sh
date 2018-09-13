#!/usr/bin/bash

## Benchmarking of Mark Duplicates for neoepiscope (GATK v4.0.0.0)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_markduplicates.sh {HOME_DIR} {SAMPLE_ID} {GATK}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## SAMPLE_ID is your tumor sample ID or normal sample ID (must be run for both)
## GATK is the path to your GATK v4 jar file

HOME_DIR=$1
SAMPLE_ID=$2
BAM=${HOME_DIR}/${SAMPLE_ID}.bam
TIMELOG=${HOME_DIR}/${SAMPLE_ID}.markduplicates.time_log
CPUINFO=${HOME_DIR}/${SAMPLE_ID}.markduplicates.cpu_data
GATK=$3

cat /proc/cpuinfo > $CPUINFO

echo 'MARK DUPLICATES RUN TIME' > $TIMELOG

time ($GATK MarkDuplicates --INPUT $BAM --METRICS_FILE ${HOME_DIR}${SAMPLE_ID}.duplicate_metrics \
		--OUTPUT ${HOME_DIR}${SAMPLE_ID}.marked.bam 2>&1) 2>> $TIMELOG || exit 1


echo 'SORT SAM RUN TIME' >> $TIMELOG

time ($GATK SortSam --INPUT ${HOME_DIR}${SAMPLE_ID}.marked.bam \
		--OUTPUT ${HOME_DIR}${SAMPLE_ID}.marked.sorted.bam --SORT_ORDER coordinate 2>&1) 2>> $TIMELOG
