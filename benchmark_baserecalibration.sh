#!/usr/bin/bash

## Benchmarking of Base Recalibration for neoepiscope (GATK v4.0.0.0)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_baserecalibration.sh {HOME_DIR} {SAMPLE_ID} {REFERENCE} {DBSNP} {GATK} {PICARD}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## SAMPLE_ID is your tumor sample ID or normal sample ID (must be run for both)
## REFERENCE is the path to your hg38 reference fasta
## DBSNP is the path to your hg38 DBSNP VCF
## GATK is the path to your GATK v4 jar file
## PICARD is the path to your picard tools AddOrReplaceReadGroups.jar file

HOME_DIR=$1
SAMPLE_ID=$2
BAM=${HOME_DIR}/${SAMPLE_ID}.marked.sorted.bam
REFERENCE=$3
DBSNP=$4
TIMELOG=${HOME_DIR}/${SAMPLE_ID}.baserecalibration.time_log
CPUINFO=${HOME_DIR}/${SAMPLE_ID}.baserecalibration.cpu_data
GATK=$5
PICARD=$6

cat /proc/cpuinfo > $CPUINFO

echo 'READ GROUP RUN TIME' > $TIMELOG

time (java -jar $PICARD I=$BAM O=${HOME_DIR}${SAMPLE_ID}.marked.sorted.reheadered.bam RGLB=lib1 RGPL=ILLUMINA \
    RGSM=${SAMPLE_ID} RGPU=foo 2>&1) 2>> $TIMELOG || exit 1

echo 'BASE RECALIBRATOR RUN TIME' >> $TIMELOG

time ($GATK BaseRecalibrator --input ${HOME_DIR}${SAMPLE_ID}.marked.sorted.reheadered.bam --reference $REFERENCE --known-sites $DBSNP \
		--output ${HOME_DIR}${SAMPLE_ID}.recalibration_table 2>&1) 2>> $TIMELOG || exit 1


echo 'APPLY RECALIBRATION RUN TIME' >> $TIMELOG

time ($GATK ApplyBQSR --input ${HOME_DIR}${SAMPLE_ID}.marked.sorted.reheadered.bam \
      --bqsr-recal-file ${HOME_DIR}${SAMPLE_ID}.recalibration_table \
		--output ${HOME_DIR}${SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam 2>&1) 2>> $TIMELOG
