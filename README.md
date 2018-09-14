Reproducing results from neoepiscope paper
-----

This repository holds scripts and instructions for reproducing the benchmarking for [neoepiscope](https://github.com/ohsu-comp-bio/neoepiscope) as described in our [manuscript]().

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

Reproducing Variant Identification and Phasing
-----

**Data availability**

We used paired tumor-normal WES data of melanoma patients from Carreno et al., Gao et al., Hugo et al., Roh et al., Snyder et al., Van Allen et al., and Zaretsky et al. (3-9); NSCLC patients from Rizvi et al. (10); and colon, endometrial, and thyroid cancer patients from Le et al. [11].

**Read alignment and BAM processing**

We aligned WES reads to the GRCh37d5 genome using the [Sanger cgpmap workflow](https://github.com/cancerit/dockstore-cgpmap). The relevant reference bundle is available [here](https://github.com/cancerit/dockstore-cgpwxs/wiki) (see "Reference bundle"), and the necessary variant files can be found in the [Broad Institute Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).

The file [fastq2bam.cwl.yaml](fastq2bam.cwl.yaml) can be used with [sample_fastq2bam.sh](sample_fastq2bam.sh) and [sample_fastq2bam.json](sample_fastq2bam.json) to run the workflow. In the sample shell script and json file, change the paths to match the relevant paths for your computer.

**MORE COMING SOON**

----

References:
-----

1. Wood MA, Nguyen A, Struck AJ, Ellrott K, Nellore A, Thompson RF. neoepiscope improves neoepitope prediction with multi-variant phasing. Preprint.

2. Bassani-Sternberg M, Bräunlein E, Klar R, Engleitner T, Sinitcyn P, Audehm S, et al. [Direct identification of clinically relevant neoepitopes presented on native human melanoma tissue by mass spectrometry](https://www.nature.com/articles/ncomms13404). Nat Commun. 2016;7: 13404.

3. Carreno BM, Magrini V, Becker-Hapak M, Kaabinejadian S, Hundal J, Petti AA, et al. A dendritic cell vaccine increases the breadth and diversity of melanoma neoantigen-specific T cells. Science. 2015;348: 803–808.

4. Gao J, Shi LZ, Zhao H, Chen J, Xiong L, He Q, et al. Loss of IFN-γ Pathway Genes in Tumor Cells as a Mechanism of Resistance to Anti-CTLA-4 Therapy. Cell. 2016;167: 397–404.e9.

5. Hugo W, Zaretsky JM, Sun L, Song C, Moreno BH, Hu-Lieskovan S, et al. Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell. 2017;168: 542.

6. Roh W, Chen P-L, Reuben A, Spencer CN, Prieto PA, Miller JP, et al. Integrated molecular analysis of tumor biopsies on sequential CTLA-4 and PD-1 blockade reveals markers of response and resistance. Sci Transl Med. 2017;9.

7. Snyder A, Makarov V, Merghoub T, Yuan J, Zaretsky JM, Desrichard A, et al. Genetic basis for clinical response to CTLA-4 blockade in melanoma. N Engl J Med. 2014;371: 2189–2199.

8. Van Allen EM, Miao D, Schilling B, Shukla SA, Blank C, Zimmer L, et al. Genomic correlates of response to CTLA-4 blockade in metastatic melanoma. Science. 2015;350: 207–211.

9. Zaretsky JM, Garcia-Diaz A, Shin DS, Escuin-Ordinas H, Hugo W, Hu-Lieskovan S, et al. Mutations Associated with Acquired Resistance to PD-1 Blockade in Melanoma. N Engl J Med. 2016;375: 819–829.

10. Rizvi NA, Hellmann MD, Snyder A, Kvistborg P, Makarov V, Havel JJ, et al. Cancer immunology. Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer. Science. 2015;348: 124–128.

11. Le DT, Durham JN, Smith KN, Wang H, Bartlett BR, Aulakh LK, et al. Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. Science. 2017;357: 409–413.
