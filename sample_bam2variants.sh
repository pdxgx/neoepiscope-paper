#!/bin/bash

set -e

mkdir -p /PATH/TO/VARIANT_OUTPUT_DIRECTORY/TUMOR_SAMPLE_NAME/tmp
mkdir -p /PATH/TO/VARIANT_OUTPUT_DIRECTORY/TUMOR_SAMPLE_NAME/cache


cwltool \
--outdir /PATH/TO/VARIANT_OUTPUT_DIRECTORY/TUMOR_SAMPLE_NAME \
--tmpdir-prefix /PATH/TO/VARIANT_OUTPUT_DIRECTORY/TUMOR_SAMPLE_NAME/tmp \
--cachedir /PATH/TO/VARIANT_OUTPUT_DIRECTORY/TUMOR_SAMPLE_NAME/cache \
/PATH/TO/mc3/mc3_variant.cwl \
/PATH/TO/sample_bam2variants.json