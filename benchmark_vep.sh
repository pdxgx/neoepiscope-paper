#!/usr/bin/bash

## Benchmarking of Variant Effect Predictor for neoepiscope (VEP v91.3)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l
## add --nodelist=exanode-[#]-[#] to request a specific node

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_vep.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {VEP_DIR} {GATK} \
##										 {REFERENCE} {REFERENCE_DICT} {PICARD} {BGZIP} {TABIX}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## VEP_DIR contains your VEP v91.3 executable, cache, plugins directory, and text file of contig synonyms (called 'chr_synonyms.txt')
## GATK is the path to your GATK jar file
## REFERENCE is the path to your hg38 reference fasta
## REFERENCE_DICT is the path to your hg38 reference dictionary file
## PICARD is the path to your picard tools jar file
## BGZIP is the path to your bgzip install
## TABIX is the path to your tabix install

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.filtered.vcf
TUMOR_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.vcf
GERMLINE_VCF=${HOME_DIR}${NORMAL_SAMPLE_ID}.germline.filtered.vcf
ADJ_GERMLINE_VCF=${HOME_DIR}${NORMAL_SAMPLE_ID}.germline.filtered.adj.vcf
ANNOTATED_VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.annotated.vcf
COMBINED_VCF=${HOME_DIR}${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.combine_variants.vcf
SORTED_VCF=${HOME_DIR}${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.combine_variants.sorted.vcf
PHASED_VCF=${HOME_DIR}${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.combine_variants.phased.vcf
ANNOTATED_PHASED_VCF=${HOME_DIR}${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.combine_variants.phased.annotated.vcf
TUMOR_BAM=${HOME_DIR}${TUMOR_SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
VEP_DIR=$4
VEP=${VEP_DIR}/vep
VEP_CACHE=${VEP_DIR}/cache/
SYNONYMS=${VEP_DIR}/chr_synonyms.txt
PLUGINS=${VEP_DIR}/plugins/
GATK=$5
REFERENCE=$6
REFERENCE_DICT=$7
PICARD=$8
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.vep.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.vep.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'ISOLATE TUMOR RUN TIME' > $TIMELOG

time (grep '^##' $VCF > $TUMOR_VCF && grep -v "^##" $VCF | grep "^#" | \
	  sed -E "s/(.+)TUMOR(.+)/\1$TUMOR_SAMPLE_ID/" >> $TUMOR_VCF \
	  && grep -v '^#' $VCF | cut -f1-10 >> $TUMOR_VCF 2>&1) 2>> $TIMELOG || exit 1

echo 'RENAME GERMLINE RUN TIME' >> $TIMELOG

time (grep '^##' $GERMLINE_VCF > $ADJ_GERMLINE_VCF && grep -v "^##" $GERMLINE_VCF | grep "^#" | \
       sed -E "s/(.+)$NORMAL_SAMPLE_ID/\1$TUMOR_SAMPLE_ID/" >> $ADJ_GERMLINE_VCF && grep -v "^#" $GERMLINE_VCF \
       >> $ADJ_GERMLINE_VCF 2>&1) 2>> $TIMELOG || exit 1

echo 'COMBINEVARIANTS RUN TIME' >> $TIMELOG

time(java -jar $GATK -T CombineVariants -R $REFERENCE --variant $ADJ_GERMLINE_VCF --variant $TUMOR_VCF \
	 -o $COMBINED_VCF --assumeIdenticalSamples 2>&1) 2>> $TIMELOG

echo 'PHASING RUN TIME' >> $TIMELOG

time(java -jar $GATK -T ReadBackedPhasing -R $REFERENCE -I $TUMOR_BAM --variant $COMBINED_VCF -L $COMBINED_VCF \
	 -o $PHASED_VCF 2>&1) 2>> $TIMELOG

echo 'PHASED VEP RUNTIME' >> $TIMELOG

time(vep --no_stats --offline --fork 4 --cache --dir_cache $VEP_CACHE --cache_version 94 --assembly GRCh38 \
	 --input_file $PHASED_VCF --format vcf --output_file $ANNOTATED_PHASED_VCF \
	 --vcf --symbol --synonyms $SYNONYMS --terms SO --plugin Downstream \
	 --plugin Wildtype --dir_plugins $PLUGINS && $BGZIP \
	 -c $ANNOTATED_PHASED_VCF > ${ANNOTATED_PHASED_VCF}.gz && $TABIX -p vcf ${ANNOTATED_PHASED_VCF}.gz 2>&1) 2>> $TIMELOG
