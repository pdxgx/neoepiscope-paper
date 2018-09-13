#!/usr/bin/bash

## Benchmarking of haplotype phasing for neoepiscope (HapCUT2, no version specified)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_hapcut2.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {HAPCUT2}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## BUILD_DIR is the path to your HapCUT2 build directory

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
TUMOR_BAM=${HOME_DIR}/${TUMOR_SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
SOMATIC_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf
GERMLINE_VCF=${HOME_DIR}/${NORMAL_SAMPLE_ID}.germline.filtered.vcf
MERGED_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.merged.vcf
FRAGMENTS=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.fragments
HAPLOTYPES=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.haplotypes
HAPCUT2=$4
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.hapcut2.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.hapcut2.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'MERGE RUN TIME' > $TIMELOG

time (neoepiscope merge -s $SOMATIC_VCF -g $GERMLINE_VCF -o $MERGED_VCF \
		2>&1) 2>> $TIMELOG || exit 1

echo 'EXTRACT HAIRS RUNTIME' >> $TIMELOG

time (${HAPCUT2}/extractHAIRS --indels 1 --bam $TUMOR_BAM --VCF $MERGED_VCF \
		--out $FRAGMENTS 2>&1) 2>> $TIMELOG || exit 1

echo 'HAPCUT RUNTIME' >> $TIMELOG

time (${HAPCUT2}/HAPCUT2 --fragments $FRAGMENTS --vcf $MERGED_VCF --output \
		$HAPLOTYPES 2>&1) 2>> $TIMELOG || exit 1
