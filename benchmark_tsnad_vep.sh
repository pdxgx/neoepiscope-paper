#!/usr/bin/bash

## Benchmarking of Variant Effect Predictor for neoepiscope (VEP v91.3)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l
## add --nodelist=exanode-[#]-[#] to request a specific node

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_tsnad_vep.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {VEP_DIR}
## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## VEP_DIR contains your VEP v91.3 executable, cache, plugins directory, and text file of contig synonyms (called 'chr_synonyms.txt')

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf
TUMOR_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor_tsnad.vcf
ANNOTATED_DATA=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.annotated.filtered.txt
TUMOR_BAM=${HOME_DIR}${TUMOR_SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
VEP_DIR=$4
VEP=${VEP_DIR}/vep
VEP_CACHE=${VEP_DIR}/cache/
SYNONYMS=${VEP_DIR}/chr_synonyms.txt
PLUGINS=${VEP_DIR}/plugins/
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.vep_tsnad.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.vep_tsnad.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'ISOLATE TUMOR RUN TIME' > $TIMELOG

time (grep '^##' $VCF > $TUMOR_VCF && grep -v "^##" $VCF | grep "^#" | \
	  sed -E "s/(.+)TUMOR(.+)/\1$TUMOR_SAMPLE_ID/" >> $TUMOR_VCF \
	  && grep -v '^#' $VCF | cut -f1-10 >> $TUMOR_VCF 2>&1) 2>> $TIMELOG || exit 1

echo 'VEP RUN TIME' >> $TIMELOG

time(vep --offline --no_stats --fork 4 --cache --dir_cache $VEP_CACHE --cache_version 94 --assembly GRCh38 \
	 --input_file $TUMOR_VCF --format vcf --output_file STDOUT --tab --canonical \
	 --symbol --synonyms $SYNONYMS --terms SO --plugin Downstream \
	 --plugin Wildtype --dir_plugins $PLUGINS | filter_vep --output_file STDOUT \
	 -filter "CANONICAL is YES and Consequence is missense_variant" | grep -E '#|ENSG' > \
	 $ANNOTATED_DATA  2>&1) 2>> $TIMELOG
