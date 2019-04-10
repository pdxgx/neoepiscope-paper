#!/usr/bin/bash

## Benchmarking of MuTect 2 for neoepiscope (GATK v3.5)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_mutect.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {REFERENCE} {DBSNP} {GATK}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## REFERENCE is the path to your hg38 reference fasta
## DBSNP is the path to your hg38 DBSNP VCF
## GATK is the path to your GATK v3.5 jar file

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
TUMOR_BAM=${HOME_DIR}/${TUMOR_SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
NORMAL_BAM=${HOME_DIR}/${NORMAL_SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
REFERENCE=$4
DBSNP=$5
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.cpu_data
GATK=$6
RAW_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.raw.vcf
FILTERED_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf

cat /proc/cpuinfo > $CPUINFO

echo 'MUTECT RUN TIME' > $TIMELOG

time (java -jar $GATK --analysis_type MuTect2 -nct 4 --input_file:tumor $TUMOR_BAM --input_file:normal $NORMAL_BAM \
    --reference_sequence $REFERENCE --dbsnp $DBSNP --out $RAW_VCF 2>&1) 2>> $TIMELOG || exit 1

echo 'FILTER MUTECT RUN TIME' > $TIMELOG

time (grep '^#' $RAW_VCF > $FILTERED_VCF && grep -v '^#' $RAW_VCF | awk '$7 == "PASS"' \
		>> $FILTERED_VCF 2>&1) 2>> $TIMELOG
