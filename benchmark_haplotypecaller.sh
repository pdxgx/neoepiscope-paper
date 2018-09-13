#!/usr/bin/bash

## Benchmarking of HaplotypeCaller for neoepiscope (GATK v3.5)

## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_haplotypecaller.sh {HOME_DIR} {SAMPLE_ID} {REFERENCE} {DBSNP} {GATK}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## REFERENCE is the path to your hg38 reference fasta
## DBSNP is the path to your hg38 DBSNP VCF
## GATK is the path to your GATK v3.5 jar file

HOME_DIR=$1
SAMPLE_ID=$2
BAM=${HOME_DIR}/${SAMPLE_ID}.marked.sorted.reheadered.recalibrated.bam
REFERENCE=$3
DBSNP=$4
TIMELOG=${HOME_DIR}/${SAMPLE_ID}.haplotypecaller.time_log
CPUINFO=${HOME_DIR}/${SAMPLE_ID}.haplotypecaller.cpu_data
GATK=$5
RAW_VCF=${HOME_DIR}/${SAMPLE_ID}.germline.raw.vcf
POST_FILTER_VCF=${HOME_DIR}/${SAMPLE_ID}.germline.postfilter.vcf
FILTERED_VCF=${HOME_DIR}/${SAMPLE_ID}.germline.filtered.vcf

cat /proc/cpuinfo > $CPUINFO

echo 'HAPLOTYPE CALLER RUN TIME' > $TIMELOG

time (java -jar $GATK --analysis_type HaplotypeCaller -nct 4 --input_file $BAM --reference_sequence $REFERENCE --dbsnp $DBSNP \
		--out $RAW_VCF 2>&1) 2>> $TIMELOG || exit 1

echo 'FILTER VARIANTS RUN TIME' >> $TIMELOG

time (java -jar $GATK --analysis_type VariantFiltration --variant $RAW_VCF --out $POST_FILTER_VCF \
    --reference_sequence $REFERENCE \
		--clusterSize 3 --clusterWindowSize 15 --missingValuesInExpressionsShouldEvaluateAsFailing \
		--filterName 'QDFilter' --filterExpression 'QD < 2.0' \
		--filterName 'QUALFilter' --filterExpression 'QUAL < 100.0' \
		--filterName DPFilter --filterExpression 'DP < 10.0' \
	2>&1) 2>> $TIMELOG || exit 1

echo 'FINAL FILTER RUN TIME' >> $TIMELOG

time (grep '^#' $POST_FILTER_VCF > $FILTERED_VCF && grep -v '^#' $POST_FILTER_VCF | awk '$7 == "PASS"' \
		>> $FILTERED_VCF 2>&1) 2>> $TIMELOG
