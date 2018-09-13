#!/usr/bin/bash

## Benchmarking of neoepitope enumeration for neoepiscope
## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_pvacseq.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {HLA_ALLELES}

## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## HLA_ALLELES is your comma-separated list of HLA alleles

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
HLA_ALLELES=$4
VCF=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.mutect.tumor.annotated.vcf
IEDB=${HOME_DIR}/iedb_data/
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.pvacseq.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.pvacseq.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'PVACSEQ RUN TIME' > $TIMELOG

time (source activate benchmark_env && pvacseq run $VCF ${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID} \
		 $HLA_ALLELES NetMHCpan $HOME_DIR -e 9 -l 17 -b 1000000 \
		--iedb-install-directory $IEDB && source deactivate benchmark_env \
		&& mv ${HOME_DIR}/MHC_Class_I/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.final.tsv ${HOME_DIR} && \
		rm -r ${HOME_DIR}/MHC_Class_I/ 2>&1) 2>> $TIMELOG || exit 1
