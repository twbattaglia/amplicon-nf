# amplicon-nf

*amplicon-nf* is a bioinformatics pipeline that can be used to process DNA data obtained from HLA-unbiased genetic screens.

This repository contains the code pertinent to our publication Nature Biotechnology. Please refer to materials and methods section of the article for details.

>Chiara M. Cattaneo, Thomas Battaglia, Jos Urbanus, Ziva Moravec, Rhianne Voogd, Rosa de Groot, Koen J. Hartemink, John B.A.G. Haanen, Emile E. Voest, Ton N. Schumacher, and Wouter Scheper, Identification of patient-specific CD4+ and CD8+ T cell neoantigens through HLA-unbiased genetic screens, _Nature Biotechnology_ (Accepted)

## Installation
The pipeline can be install by pulling the repo (below) or by using the command: `nextflow run twbattaglia/amplicon-nf`
```
git clone https://github.com/twbattaglia/amplicon-nf
```

## Usage

```
nextflow run amplicon-nf/main.nf \
--fastq 'library/*.fq' \
--library 'library.txt' \
--length 63 \
--mode 'library' \
--outdir 'results' 
```

## Required inputs

#### Demultiplexed FASTQ files
The input data must be unmodified demultiplexed FASTQ files (quality filtering will be performed within the pipeline). Depending upon the set up, the data maybe need to be reverse complemented before downstream processing.

#### Library file
This file refers a TSV table (NO HEADERS) that contains (1) A unique oligo indentifier, (2) the genomic context of the neoantigen and (3) the unique barcode that refers to the sequence. See below for an example:

```
oligo_1	CTGAGCTGGGTGAAACACGACGTGGACGCCCGCAGGCAGCATGTCTCACGGCTCATGAAGTGTGTGCGGCTGCCCTTGCTGAGCCGCGACTTC	TGGTAGTCGCAT
oligo_2	GACGGCATCGTGACCGCCGAGGAGCTGGAGAACGTGCCCACACTCTCGCTGCAGCCAATAGGCACCTTAAATAGCCACTTCGTGCGGCTGGCC	CTCAAGTCGCAT
oligo_3	GGTCAGCGGTGGACCTTCAGCCCCTCCTGCCTGGTGGCCTACCGGCTTGAGGAGGATGCCAACCTGGACGTGGCCGAGCGCGCCCGGGAGAAC	AGTCGCTCGCAT
oligo_4	CCCCCCTTCTTCAGCAAGGAGCAGCCACAGGCCTTGAACTTTGGATGCATTGGGATGGTGATCGGGCACGAGATCACGCACGGCTTTGACGAC	TGAACCTCGCAT
```


## Profiling modes
There are two modes available for processessing and quantifying neoantigen oligo sequences (*Library* and *Barcode*). The `Library` mode is meant for data that is generatate from a single sequence that refers to the genomics context of the neoantigen with the mutation centered (See _Methods_ section of the manuscript). The `Barcode` mode is meant for data that contains replicates of both the mutant and wildtype neoantigen sequence and contains a unique 12 nt barcode sequence.


### Mode 1: Library
```
nextflow run amplicon-nf/main.nf \
--fastq 'demulti_fastq/*.fq.gz' \
--library 'ITO66-oligos.txt' \
--length 63 \
--mode 'library' \
--mapper 'bbmap' \
--outdir 'results'
```

### Mode 2: Barcode

```
nextflow run amplicon-nf/main.nf \
--fastq 'demulti_fastq/*.fq.gz' \
--library 'TAA-oligo.txt' \
--length 12 \
--reverse \
--mode 'barcode' \
--mapper 'bowtie2' \
--outdir 'results' 
```

## Parameters
```
=============================================================
 Amplicon-seq Nextflow Pipeline vnull
 Author: Thomas W. Battaglia
============================================================
Usage:
The typical command for running the pipeline is as follows:

Required arguments:
  --fastq [file]                        Folder of FASTQ files that will be used for downstream analysis.
  --library [file]                      Tab-separated file (Column 1 = ID, Column 2 = Library) describing the barcode sequences associated with each entry in the library.
  --mode [str]                          Processing modes:
                                            --mode 'library' | Specifies that the reads will be mapped directly to the library. Assumes a perfect match for quantification.
                                            --mode 'barcode' | Specifies that the reads have barcodes at the end of the read that are used for quantification.
  --length [int]                        Will be the length of the oligo or length of the barcode attached.
  --outdir [file]                       Local directory to store results

Optional arguments:
  --print [bool]                        Option to print channel of files without running full pipeline (Default = False)
  --check [bool]                        Only run the first few steps to get a sense of the sequence data distribution. (Default = False)
  --reverse [bool]                      Option to reverse complement the initial reads before analysis (Default = False)
  --cut5 [int]                          Parameter to cut 5' end in Bowtie2 alignment (Default = 0)
  --cut3 [int]                          Parameter to cut 3' end in Bowtie2 alignment (Default = 0)
  --quality [int]                       Minimum base quality score required when filtering with Cutadapt. (Default = 10)
  --univec [file]                       Path to UniVec database in FASTQ format. (Default is taken from within pipeline)
  --dedeup [bool]                       Option to enable removal of optical/PCR duplicates (Default = False)
  --perfect [bool]                      Option to enable perfect matching with BBMap (Default = False)
  --semiperfect [bool]                  Option to enable semi-perfect matching with BBMap (Default = False)
  --mapper [str]                        Select an alignment tool when mapping reads to barcodes. Can be of :
                                            --mapper 'bbmap' | Default
                                            --mapper 'bowtie2'
```


