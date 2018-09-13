#!/usr/bin/bash

## Benchmarking of Variant Effect Predictor for neoepiscope (VEP v91.3)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l
## add --nodelist=exanode-[#]-[#] to request a specific node

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_vep.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {VEP_DIR}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## VEP_DIR contains your VEP v91.3 executable, cache, plugins directory, and text file of contig synonyms (called 'chr_synonyms.txt')

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf
TUMOR_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.vcf
ANNOTATED_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.annotated.vcf
VEP_DIR=$4
VEP=${VEP_DIR}/vep
VEP_CACHE=${VEP_DIR}/cache/
SYNONYMS=${VEP_DIR}/chr_synonyms.txt
PLUGINS=${VEP_DIR}/plugins/
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.vep.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.vep.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'ISOLATE TUMOR RUN TIME' > $TIMELOG

time (grep '^##' $VCF > $TUMOR_VCF && grep -v '^##' $VCF | cut -f1-10 >> $TUMOR_VCF \
		2>&1) 2>> $TIMELOG || exit 1

echo 'VEP RUN TIME' >> $TIMELOG

time($VEP --offline --fork 4 --cache --dir_cache $VEP_CACHE --cache_version 91 --assembly GRCh38 \
		--input_file $TUMOR_VCF --format vcf --output_file $ANNOTATED_VCF \
		--vcf --symbol --synonyms $SYNONYMS --terms SO --plugin Downstream \
		--plugin Wildtype --dir_plugins $PLUGINS 2>&1) 2>> $TIMELOG
