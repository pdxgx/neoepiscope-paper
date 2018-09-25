Reproducing results from neoepiscope paper
-----

This repository holds scripts and instructions for reproducing the benchmarking for [neoepiscope](https://github.com/ohsu-comp-bio/neoepiscope) as described in our [manuscript](https://www.biorxiv.org/content/early/2018/09/14/418129).

----

Requirements:
-----

##### Software:

[Anaconda2](https://www.anaconda.com/download/)

[BWA v0.7.9a-r786](https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2/download)

[GATK v3.5](https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.5-0-g36282e4)

[GATK v4.0.0.0](https://github.com/broadinstitute/gatk/releases/download/4.0.0.0/gatk-4.0.0.0.zip)

[HapCUT2](https://github.com/vibansal/HapCUT2) (we used commit 0f4ef9a0ea07d7eb32b306debee96c4ad00a13c2)

[MuPeXI v1.2.0](https://github.com/ambj/MuPeXI)

[neoepiscope v0.1.0](https://github.com/ohsu-comp-bio/neoepiscope)

[Picard Tools v1.110](https://github.com/broadinstitute/picard/releases/tag/1.110)

[pvactools v1.0.8](https://github.com/griffithlab/pVACtools)

[Samtools v1.3.1](https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2/download)

[Variant Effect Predictor (VEP) v91.3](https://github.com/Ensembl/ensembl-vep)


##### Reference files:

Hg38 reference fasta (available from the Broad Institute [resource bundle](https://software.broadinstitute.org/gatk/download/bundle))

Hg38 DBSNP VCF (available from the Broad Institute [resource bundle](https://software.broadinstitute.org/gatk/download/bundle))

VEP GRCh38 cache (available from [Ensembl](https://uswest.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache))

[Downstream VEP plugin](https://github.com/Ensembl/VEP_plugins/blob/release/93/Downstream.pm)


##### Data files:

We used WES paired-end fastq files from Bassani-Sternberg et al. [1] for benchmarking. The files are available at the [European Genome-phenome Archive](https://www.ebi.ac.uk/ega/home) (EGA) under accession number [EGAS00001002050](https://www.ebi.ac.uk/ega/studies/EGAS00001002050). We used matched tumor and normal samples from patients Mel5, Mel8, and Mel12.

----

Before you start:
-----

1) Create a python3 virtual environment w/ conda called benchmark\_env:

```conda create -n benchmark_env python=3.4```

2) Move the VEP cache data to a subdirectory called 'cache' in your VEP main directory

3) Copy the 'chr_synonyms.txt' file in the VEP cache data to your VEP main directory

4) Move the VEP Downstream plugin to a subdirectory called 'plugins' in your VEP main directory

----

Running benchmarking scripts (Benchmarking section of MATERIALS AND METHODS)
-----

We ran our benchmarking on an exclusive node of our institution's computer cluster, using the node's first four processors for multithreading (when possible). Each script benchmarks a different step in the neoepitope calling pipeline, and produces relevant output files from these steps along with files summarizing the CPU information and run time information (we used real time in our assessments). As an example, below are the commands to run the benchmarking for patient Mel5 from within the repository:


**1) Benchmark BWA performance for tumor and normal samples:**

_SCRIPTS:_

[benchmark_bwa.sh](benchmark_bwa.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_bwa.sh HOME_DIRECTORY Mel5_tumor TUMOR_FASTQ1 TUMOR_FASTQ2 REFERENCE_FASTA BWA SAMTOOLS```

```taskset -c 0,1,2,3 ./benchmark_bwa.sh HOME_DIRECTORY Mel5_normal NORMAL_FASTQ1 NORMAL_FASTQ2 REFERENCE_FASTA BWA SAMTOOLS```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

TUMOR\_FASTQ1/2 and NORMAL\_FASTQ1/2 are the paired end tumor and normal fastqs

REFERENCE\_FASTA is the path to your reference fasta

BWA is the path to your BWA executable

SAMTOOLS is the path to your samtools executable

_OUTPUTS:_

CPU info is output in Mel5\_tumor.bwa.cpu\_data and Mel5\_normal.bwa.cpu\_data 

Run time info is output in Mel5\_tumor.bwa.time_log and Mel5\_normal.bwa.time\_log

----

**2) Benchmark BAM processing for tumor and normal samples**

_SCRIPTS:_

[benchmark_markduplicates.sh](benchmark_markduplicates.sh)

[benchmark_baserecalibration.sh](benchmark_baserecalibration.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_markduplicates.sh HOME_DIRECTORY Mel5_tumor GATK```

```taskset -c 0,1,2,3 ./benchmark_markduplicates.sh HOME_DIRECTORY Mel5_normal GATK```

```taskset -c 0,1,2,3 ./benchmark_baserecalibration.sh HOME_DIRECTORY Mel5_tumor REFERENCE_FASTA DBSNP GATK PICARD```

```taskset -c 0,1,2,3 ./benchmark_baserecalibration.sh HOME_DIRECTORY Mel5_normal REFERENCE_FASTA DBSNP GATK PICARD```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

REFERENCE\_FASTA is the path to your reference fasta

DBSNP is the path to your DBSNP VCF

GATK is the path to your GATK jar file

PICARD is the path to your PICARD executable

_OUTPUTS:_

CPU info is output in Mel5\_tumor.markduplicates.cpu\_data, Mel5\_normal.markduplicates.cpu\_data, Mel5\_tumor.baserecalibration.cpu\_data, and Mel5\_normal.baserecalibration.cpu\_data

Run time info is output in Mel5\_tumor.markduplicates.time\_log, Mel5\_normal.markduplicates.time\_log, Mel5\_tumor.baserecalibration.time\_log, and Mel5\_normal.baserecalibration.time\_log

----

**3) Benchmark somatic variant calling**

_SCRIPTS:_

[benchmark_mutect.sh](benchmark_mutect.sh)

[benchmark_filter_mutect.sh](benchmark_filter_mutect.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_mutect.sh HOME_DIRECTORY Mel5_tumor Mel5_normal REFERENCE_FASTA DBSNP GATK```

```taskset -c 0,1,2,3 ./benchmark_filter_mutect.sh HOME_DIRECTORY Mel5_tumor Mel5_normal```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

REFERENCE\_FASTA is the path to your reference fasta

DBSNP is the path to your DBSNP VCF

GATK is the path to your GATK jar file

_OUTPUTS:_

CPU info is output in Mel5\_tumor\_v\_Mel5\_normal.mutect.cpu\_data, Mel5\_tumor\_v\_Mel5\_normal.mutect.cpu\_data, Mel5\_tumor\_v\_Mel5\_normal.filtermutect.cpu\_data and Mel5\_tumor\_v\_Mel5\_normal.filtermutect.cpu\_data, 

Run time info is output in Mel5\_tumor\_v\_Mel5\_normal.mutect.time\_log, Mel5\_tumor\_v\_Mel5\_normal.mutect.time\_log, Mel5\_tumor\_v\_Mel5\_normal.filtermutect.time\_log, and Mel5\_tumor\_v\_Mel5\_normal.filtermutect.time\_log

----

**4) Benchmark germline variant calling**

_SCRIPTS:_

[benchmark_haplotypecaller.sh](benchmark_haplotypecaller.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_haplotypecaller.sh HOME_DIRECTORY Mel5_normal REFERENCE_FASTA DBSNP GATK```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

REFERENCE\_FASTA is the path to your reference fasta

DBSNP is the path to your DBSNP VCF

GATK is the path to your GATK jar file

_OUTPUTS:_

CPU info is output in Mel5\_normal.haplotypecaller.cpu\_data

Run time info is output in Mel5\_normal.haplotypecaller.time\_log

----

**5) Benchmark haplotype phasing**

_SCRIPTS:_

[benchmark_hapcut2.sh](benchmark_hapcut2.sh)

[benchmark_nogermline_hapcut2.sh](benchmark_nogermline_hapcut2.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_hapcut2.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HAPCUT2```

```taskset -c 0,1,2,3 ./benchmark_nogermline_hapcut2.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HAPCUT2```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

HAPCUT2 is the path to your HapCUT2 build directory

_OUTPUTS:_

CPU info is output in Mel5\_tumor\_v\_Mel5\_normal.hapcut2.cpu\_data and Mel5\_tumor\_v\_Mel5\_normal.hapcut2\_tumor.cpu\_data

Run time info is output in Mel5\_tumor\_v\_Mel5\_normal.hapcut2.time\_log and Mel5\_tumor\_v\_Mel5\_normal.hapcut2.time\_log

----

**6) Benchmark HLA typing**

_SCRIPTS:_

[benchmark_optitype.sh](benchmark_optitype.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_optitype.sh HOME_DIRECTORY Mel5_tumor TUMOR_FASTQ1 TUMOR_FASTQ2 OPTITYPE CONFIG```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

TUMOR\_FASTQ1/2 are the paired end tumor fastqs

OPTITYPE is the path to your Optitype python script

CONFIG is the path to your Optitype config file

_OUTPUTS:_

CPU info is output in Mel5\_tumor.optitype.cpu\_data

Run time info is output in Mel5\_tumor.optitype.time\_log

----

**7) Benchmark pVACseq**

_SCRIPTS:_

[benchmark_vep.sh](benchmark_vep.sh)

[benchmark_pvacseq.sh](benchmark_pvacseq.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_vep.sh HOME_DIRECTORY Mel5_tumor Mel5_normal VEP_DIRECTORY```

```taskset -c 0,1,2,3 ./benchmark_pvacseq.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HLA-A*01:01,HLA-B*08:01,HLA-C*07:01```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

VEP\_DIRECTORY is the path to your VEP directory

_OUTPUTS:_

CPU info is output in Mel5\_tumor\_v\_Mel5\_normal.vep.cpu\_data and Mel5\_tumor\_v\_Mel5\_normal.pvacseq.cpu\_data

Run time info is output in Mel5\_tumor\_v\_Mel5\_normal.vep.time\_log and Mel5\_tumor\_v\_Mel5\_normal.pvacseq.time\_log

pVACseq neoepitopes are output in Mel5\_tumor\_v\_Mel5\_normal.final.tsv

----

**8) Benchmark MuPeXI**

_SCRIPTS:_

[benchmark_mupexi.sh](benchmark_mupexi.sh)

_COMMANDS:_

```taskset -c 0,1,2,3 ./benchmark_mupexi.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HLA-A01:01,HLA-B08:01,HLA-C07:01 MUPEXI```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

MUPEXI is the path to your MuPeXI python script

_OUTPUTS:_

CPU info is output in Mel5\_tumor\_v\_Mel5\_normal.mupexi.cpu\_data

Run time info is output in Mel5\_tumor\_v\_Mel5\_normal.mupexi.time\_log

MuPeXI neoepitopes are output in Mel5\_tumor\_v\_Mel5\_normal.mupexi

----

**9) Benchmark neoepiscope**

_SCRIPTS_:

[benchmark_neoepiscope.sh](benchmark_neoepiscope.sh)

[benchmark_nogermline_neoepiscope.sh](benchmark_nogermline_neoepiscope.sh)

[benchmark_comprehensive_neoepiscope.sh](benchmark_comprehensive_neoepiscope.sh)

_COMMANDS_:

```taskset -c 0,1,2,3 ./benchmark_neoepiscope.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HLA-A*01:01,HLA-B*08:01,HLA-C*07:01```

```taskset -c 0,1,2,3 ./benchmark_nogermline_neoepiscope.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HLA-A*01:01,HLA-B*08:01,HLA-C*07:01```

```taskset -c 0,1,2,3 ./benchmark_comprehensive_neoepiscope.sh HOME_DIRECTORY Mel5_tumor Mel5_normal HLA-A*01:01,HLA-B*08:01,HLA-C*07:01```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

_OUTPUTS:_

CPU info is output in Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.cpu\_data, Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.tumor.cpu\_data and Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.comprehensive.cpu\_data

Run time info is output in Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.time\_log, Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.tumor.time\_log and Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.comprehensive.time\_log

neoepiscope neoepitopes are output in Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.out, Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.tumor.out and Mel5\_tumor\_v\_Mel5\_normal.neoepiscope.comprehensive.out

----

**Post-processing**

The same commands above can be run for patients Mel8 and Mel12 by replacing 'Mel5' with 'Mel8' or 'Mel12' in the commands.

To compile which epitope sequences were enumerated by each caller and for each patient, you can run our python script to parse the output of each tool for each patient.

```python epitope_comparison.py -d HOME_DIRECTORY```

_INPUTS:_

HOME\_DIRECTORY is the path to this repository

_OUTPUTS:_

Mel5.peptide\_overlap.out, Mel8.peptide\_overlap.out, Mel12.peptide\_overlap.out and combined.peptide\_overlap.out are tab-delimited text files summarizing the neoepitopes predicted by all tools and which tools predicted them (on a per-patient basis, or for all patients combined).

----

Identifying phased variants (Variant identification and phasing section of MATERIALS AND METHODS)
-----

**Data availability**

We used paired tumor-normal WES data of melanoma patients from Carreno et al., Gao et al., Hugo et al., Roh et al., Snyder et al., Van Allen et al., and Zaretsky et al. (3-9); NSCLC patients from Rizvi et al. (10); and colon, endometrial, and thyroid cancer patients from Le et al. [11]. 

**Read alignment and BAM processing**

We aligned WES reads to the GRCh37d5 genome and generated genome-coordinate sorted alignments with duplicates marked using the [Sanger cgpmap workflow](https://github.com/cancerit/dockstore-cgpmap) (commit 0bacb0bee2e5c04b268c629d589ff1c551d34745). We realigned around indels and perform base recalibration using [gatk-cocleaning-tool](https://github.com/OpenGenomics/gatk-cocleaning-tool) (commit d2bafc23221f6a8dceedd45a534163e0e1bf5c68). The relevant reference bundle is available [here](https://github.com/cancerit/dockstore-cgpwxs/wiki) (see "Reference bundle"), and the necessary variant files can be found in the [Broad Institute Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).

The file [fastq2bam.cwl.yaml](fastq2bam.cwl.yaml) can be used with [sample_fastq2bam.sh](sample_fastq2bam.sh) and [sample_fastq2bam.json](sample_fastq2bam.json) to run the workflow. In the sample shell script and json file, change the paths to match the relevant paths for your computer.

**Somatic variant calling**

We used the [mc3 workflow](https://github.com/OpenGenomics/mc3) to call somatic variants. The relevant reference bundle is available [here](https://github.com/cancerit/dockstore-cgpwxs/wiki) (see "Reference bundle"), and the necessary variant files can be found in the [Broad Institute Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).

The file [mc3_variant.cwl](https://github.com/OpenGenomics/mc3/blob/master/mc3_variant.cwl) (commit 72a24b55544e3011ede1c46b13d531a7d05ef4e0) can be used with [sample_bam2variants.sh](sample_bam2variants.sh) and [sample_bam2variants.json](sample_bam2variants.json) to run the workflow. In the sample shell script and json file, change the paths to match the relevant paths for your computer.

**Somatic variant processing/consensus calling**

Following variant calling with mc3, we processed VCFs produced by MuSE, MuTect, Pindel, RADIA, SomaticSniper, and VarScan2 using [vt](https://genome.sph.umich.edu/wiki/Vt):

First, we normalized each VCF:

```vt normalize INPUT_VCF -n -r /PATH/TO/core_ref_GRCh37d5/genome.fa -o OUTPUT_VCF```

Then we decomposed block substitutions:

```vt decompose_blocksub -a -p INPUT_VCF -o OUTPUT_VCF```

Then we decomposed multi-allelic variants:

```vt decompose -s INPUT_VCF -o OUTPUT_VCF```

Then we sorted variants:

```vt sort INPUT_VCF -o OUTPUT_VCF```

Then we took uniq variants:

```vt uniq INPUT_VCF -o OUTPUT_VCF```

After processing VCFs with vt, we ran [VCF_parse.py](VCF_parse.py) to produce a consensus call set of variants produced by at least 2 callers (or called by Pindel and overlapping at least one other variant):

```python2.7 VCF_parse.py -v MuSE.uniq.vcf,MuTect.uniq.vcf,Pindel.uniq.vcf,RADIA.uniq.vcf,SomaticSniper.uniq.vcf,VarScan2.uniq.vcf -c MuSE,MuTect,Pindel,RADIA,SomaticSniper,VarScan2 -o /PATH/TO/OUTPUT_DIRECTORY -n 2 -s SAMPLE_NAME```

This outputs a file called `SAMPLE_NAME.consensus.vcf` in the `OUTPUT_DIRECTORY` used for downstream analyses.

**Germline variant calling and processing**

We used [GATK v3.7](https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.7-0-gcfedb67) to run HaplotypeCaller for calling germline variants, and to run VariantFiltration for filtering the results of HaplotypeCaller. The relevant reference bundle is available [here](https://github.com/cancerit/dockstore-cgpwxs/wiki) (see "Reference bundle"), and the necessary variant files can be found in the [Broad Institute Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).

HaplotypeCaller was run using default options:

```java -jar GenomeAnalysisTK.jar -R /PATH/TO/core_ref_GRCh37d5/genome.fa -T HaplotypeCaller -I /PATH/TO/NORMAL.realigned.cleaned.bam -o /PATH/TO/germline.vcf```

VariantFiltration was run as below:

```java -jar GenomeAnalysisTK.jar -R /PATH/TO/core_ref_GRCh37d5/genome.fa -T VariantFiltration --variant /PATH/TO/germline.vcf -o /PATH/TO/germline.filtered.vcf --clusterSize 3 --clusterWindowSize 15 --missingValuesInExpressionsShouldEvaluateAsFailing --filterName 'QDFilter' --filterExpression 'QD < 2.0' --filterName 'QUALFilter' --filterExpression 'QUAL < 100.0' --filterName DPFilter --filterExpression 'DP < 10.0'```

Only variants passing all filters were retained for downstream analyses.

**Haplotype phasing**

Before performing haplotype phasing, we swapped the sample columns in our consensus somatic variants VCF and merged our filtered germline and consensus somatic variants using [`neoepiscope`](https://github.com/pdxgx/neoepiscope):

```neoepiscope swap -i SOMATIC_VCF -o SWAPPED_SOMATIC_VCF```

```neoepiscope merge -g GERMLINE_VCF -s SWAPPED_SOMATIC_VCF -o MERGED_VCF```

Then, we predicted haplotypes using [HapCUT2](https://github.com/vibansal/HapCUT2):

```extractHAIRS --indels 1 --bam TUMOR_SAMPLE_NAME.realigned.cleaned.bam --VCF MERGED_VCF --out FRAGMENT_FILE```

```HAPCUT2 --fragments FRAGMENT_FILE --vcf MERGED_VCF --output HAPLOTYPES```

Finally, we prepared our haplotype predictions for `neoepiscope` neoepitope prediction using `neoepiscope`:

```neoepiscope prep -v MERGED_VCF -c HAPLOTYPES -o PREPPED_HAPLOTYPES```

**Neoepitope prediction**

We predicted neoepitopes of 8-24 amino acids in length with `neoepiscope`:

```neoepiscope call -b hg19 -c PREPPED_HAPLOTYPES -o OUTPUT_FILE -k 8,24```

**Phasing prevalence**

We used a combination of two scripts to analyze phasing. First, we ran [phasing_analysis.py](phasing_analysis.py) to identify instances of variant phasing within 33 or 72 bp across all patients and tumors:

```python phasing_analysis.py -o OUTPUT_DIR -n NEOEPISCOPE_DATA_DIR -c HAPLOTYPE_DIR -d 33```

```python phasing_analysis.py -o OUTPUT_DIR -n NEOEPISCOPE_DATA_DIR -c HAPLOTYPE_DIR -d 72```

The HAPLOTYPE_DIR is the directory containing the prepped haplotype files from HapCUT2/`neoepiscope merge` used to predict neoepitopes. For consistency with the script, the naming convention on these files should be `PATIENT_ID.TUMOR_SAMPLE_ID.hapcut.out.prepped`. The NEOEPISCOPE_DATA_DIR is the directory into which `neoepiscope download` saves it's data - make sure that you have run the downloader and selected 'yes' to download the hg19 GTF and bowtie index files.

To produce summary statistics and figures, we used the R script [phasing_stats.R](phasing_stats.R).

----

References:
-----

1. Wood MA, Nguyen A, Struck AJ, Ellrott K, Nellore A, Thompson RF. [neoepiscope improves neoepitope prediction with multi-variant phasing](https://www.biorxiv.org/content/early/2018/09/14/418129). Preprint.

2. Bassani-Sternberg M, Bräunlein E, Klar R, Engleitner T, Sinitcyn P, Audehm S, et al. [Direct identification of clinically relevant neoepitopes presented on native human melanoma tissue by mass spectrometry](https://www.nature.com/articles/ncomms13404). Nat Commun. 2016;7: 13404.

3. Carreno BM, Magrini V, Becker-Hapak M, Kaabinejadian S, Hundal J, Petti AA, et al. [A dendritic cell vaccine increases the breadth and diversity of melanoma neoantigen-specific T cells](http://science.sciencemag.org/content/348/6236/803). Science. 2015;348: 803–808.

4. Gao J, Shi LZ, Zhao H, Chen J, Xiong L, He Q, et al. [Loss of IFN-γ Pathway Genes in Tumor Cells as a Mechanism of Resistance to Anti-CTLA-4 Therapy](https://www.sciencedirect.com/science/article/pii/S0092867416311679?via%3Dihub). Cell. 2016;167: 397–404.e9.

5. Hugo W, Zaretsky JM, Sun L, Song C, Moreno BH, Hu-Lieskovan S, et al. [Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma](https://www.sciencedirect.com/science/article/pii/S009286741630215X?via%3Dihub). Cell. 2017;168: 542.

6. Roh W, Chen P-L, Reuben A, Spencer CN, Prieto PA, Miller JP, et al. [Integrated molecular analysis of tumor biopsies on sequential CTLA-4 and PD-1 blockade reveals markers of response and resistance](http://stm.sciencemag.org/content/9/379/eaah3560.short). Sci Transl Med. 2017;9.

7. Snyder A, Makarov V, Merghoub T, Yuan J, Zaretsky JM, Desrichard A, et al. [Genetic basis for clinical response to CTLA-4 blockade in melanoma](https://www.nejm.org/doi/full/10.1056/NEJMoa1406498). N Engl J Med. 2014;371: 2189–2199.

8. Van Allen EM, Miao D, Schilling B, Shukla SA, Blank C, Zimmer L, et al. [Genomic correlates of response to CTLA-4 blockade in metastatic melanoma](http://science.sciencemag.org/content/350/6257/207.long). Science. 2015;350: 207–211.

9. Zaretsky JM, Garcia-Diaz A, Shin DS, Escuin-Ordinas H, Hugo W, Hu-Lieskovan S, et al. [Mutations Associated with Acquired Resistance to PD-1 Blockade in Melanoma](https://www.nejm.org/doi/full/10.1056/NEJMoa1604958). N Engl J Med. 2016;375: 819–829.

10. Rizvi NA, Hellmann MD, Snyder A, Kvistborg P, Makarov V, Havel JJ, et al. Cancer immunology. [Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer](http://science.sciencemag.org/content/348/6230/124.long). Science. 2015;348: 124–128.

11. Le DT, Durham JN, Smith KN, Wang H, Bartlett BR, Aulakh LK, et al. [Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade](http://science.sciencemag.org/content/357/6349/409.long). Science. 2017;357: 409–413.
