#!/usr/bin/bash

## Benchmarking of haplotype phasing for neoepiscope (HapCUT2, no version specified)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_nogermline_hapcut2.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {HAPCUT2}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## BUILD_DIR is the path to your HapCUT2 build directory

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
TUMOR_BAM=${HOME_DIR}/${TUMOR_SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
VCF=${HOME_DIR}/old_vcfs/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.vcf
FRAGMENTS=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.tumor.fragments
HAPLOTYPES=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.tumor.haplotypes
HAPCUT2=$4
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.hapcut2_tumor.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.hapcut2_tumor.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'EXTRACT HAIRS RUNTIME' >> $TIMELOG

time (${HAPCUT2}/extractHAIRS --indels 1 --bam $TUMOR_BAM --VCF $VCF \
		--out $FRAGMENTS 2>&1) 2>> $TIMELOG || exit 1

echo 'HAPCUT RUNTIME' >> $TIMELOG

time (${HAPCUT2}/HAPCUT2 --fragments $FRAGMENTS --vcf $VCF --output \
		$HAPLOTYPES 2>&1) 2>> $TIMELOG || exit 1
