#!/usr/bin/bash

## Benchmarking of neoepitope enumeration for neoepiscope
## Request an exclusive interactive node:
## srun --partition=exacloud --time=1-12:00:00 --exclusive --pty /usr/bin/bash -l

## Run script w/ taskset:
## taskset -c 0,1,2,3 ./benchmark_neopredpipe.sh {HOME_DIR} {TUMOR_SAMPLE_ID} {NORMAL_SAMPLE_ID} {NEOPREDPIPE} {HLA_FILE}
## HOME_DIR is the path to the directory containing your benchmarking scripts and result files
## TUMOR_SAMPLE_ID is your tumor sample ID
## NORMAL_SAMPLE_ID is your normal sample ID
## HLA_ALLELES is your comma-separated list of HLA alleles
## NEOPREDPIPE is the path to your NeoPredPipe main_netMHCpan_pipe.py script
## HLA_FILE is the path to your NeoPredPipe HLA input file

HOME_DIR=$1
TUMOR_SAMPLE_ID=$2
NORMAL_SAMPLE_ID=$3
NEOPREDPIPE=$4
HLA_FILE=$5
TIMELOG=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.neopredpipe.time_log
CPUINFO=${HOME_DIR}/${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID}.neopredpipe.cpu_data

cat /proc/cpuinfo > $CPUINFO

echo 'NEOPREDPIPE RUN TIME' > $TIMELOG

#mkdir $VCF_DIR && cp $VCF $VCF_DIR                                                                            

time (mkdir $VCF_DIR && cp $VCF $VCF_DIR && python $NEOPREDPIPE -I $VCF_DIR -H $HLA_FILE \
	  -o $HOME_DIR -n ${TUMOR_SAMPLE_ID}_v_${NORMAL_SAMPLE_ID} -E 9 -l -a 2>&1) 2>> $TIMELOG
