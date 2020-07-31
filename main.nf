#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    =============================================================
     Amplicon-seq Nextflow Pipeline v${workflow.manifest.version}
     Author: ${workflow.manifest.author}
    ============================================================
    Usage:
    The typical command for running the pipeline is as follows:

    Required arguments:
      --fastq [file]                        Folder of FASTQ files that will be used for downstream analysis.
      --library [file]                      Tab-separated file (Column 1 = ID, Column 2 = Library) describing the barcode sequences associated with each entry in the library.
      --mode [str]                          Processing modes:
                                                --mode 'library' | Specifies that the reads will be mapped directly to the library. Assumes a perfect match for quantification.
                                                --mode 'barcode' | Specifies that the reads have barcodes at the end of the read that are used for quantification.
      --length [int]                        Will be the length of the oligo or legnth of the barcode attached.
      --outdir [file]                       Can be a local folder (will have egress fees) or GCP bucket folder.

    Optional arguments:
      --print [bool]                        Option to print channel of files without running full pipeline (Default = False)
      --check [bool]                        Only run the first few steps to get a sense of the sequence data distribution. (Default = False)
      --quality [int]                       Minimum base quality score required when filtering with Cutadapt. (Default = 20)
      --univec [file]                       Path to UniVec database in FASTQ format. (Default is taken from within pipeline)
      --mapper [str]                        Select an alignment tool when mapping reads to library. Can be of :
                                                --mapper 'bowtie2' | Default
                                                --mapper 'bbmap'
      --perfect [bool]                      Option to set BBmap to only select perfect matches (Default = true)
      --mismatch [int]                      Number of mismatches allowed. Filtered after bowtie2 alignment
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Create a channel of input files and gather sample identifiers
if(params.print){
  Channel
      .fromPath(params.fastq)
      .map { file -> tuple(file.simpleName, file) }
      .ifEmpty { error "Cannot find any reads matching: ${params.fastq}" }
      .println()
} else {
  Channel
      .fromPath(params.fastq)
      .map { file -> tuple(file.simpleName, file) }
      .ifEmpty { error "Cannot find any reads matching: ${params.fastq}" }
      .into { fastq_ch1; fastq_ch2 ; fastq_ch3}
}

// Check if required arguments are given
if (params.length == 0) {
    exit 1, "Invalid aligner option: ${params.length}. Must supply the length of the amplicon."
}

// Gather the adapters/vector file
univec_ch = file("$baseDir/database/UniVec.fa")
if (univec_ch.isEmpty()) {exit 1, "File ${univec_ch.getName()} is empty!"}

// Check if library is given properly
library_ch = file(params.library)
if (library_ch.isEmpty()) {exit 1, "File ${library_ch.getName()} is empty!"}

// Run FASTQC on each sample for initial quality control
process fastqc {
  publishDir "$params.outdir/01_quality/$sample_id", mode: 'copy'
  cpus 2

  input:
    set sample_id, file(reads) from fastq_ch1

  output:
    file("*") into fastqc_raw

  script:
    """
    fastqc \
    --outdir . \
    --format fastq \
    --quiet \
    --threads ${task.cpus} \
    ${reads}
    """
}

// Decontaminate reads of vector sequences
process remove_vector {
  publishDir "$params.outdir/02_vector/$sample_id", mode: 'copy'
  cpus 2

  input:
    set sample_id, file(reads) from fastq_ch2

  output:
    set val(sample_id), file("${sample_id}-decontam.fq.gz") into decontam_fastq
    file("${sample_id}*.{html,zip}") into fastqc_decontam

  script:
    """
    # Remove adapters and vectors
    fastq-mcf \
    -o ${sample_id}-decontam.fq \
    ${univec_ch} \
    ${reads}

    # Gzip to save space
    gzip ${sample_id}-decontam.fq

    # Run FastQC on decontaminated data
    fastqc \
    --outdir . \
    --quiet \
    --threads ${task.cpus} \
    ${sample_id}-decontam.fq.gz
    """
}

// Remove low quality reads from the decontaminated sequences
process trim_filter {
  publishDir "$params.outdir/03_filter/$sample_id", mode: 'copy'
  cpus 2

  input:
    set sample_id, file(reads) from decontam_fastq

  output:
    set val(sample_id), file("${sample_id}_trimmed-right.fq.gz") into filter_fastq
    file("${sample_id}*.{html,zip}") into fastqc_filter

  when:
    params.check == false & params.mode == "library"

  script:
    min_len = ( params.length - 3 )
    trim_len = ( params.length - 1 )
    """
    trim_galore \
    --basename ${sample_id} \
    --output_dir . \
    --length ${min_len} \
    --quality ${params.quality} \
    --fastqc \
    ${reads}

    # Add maximum read cutoff
    reformat.sh \
    in=${sample_id}_trimmed.fq.gz \
    out=${sample_id}_trimmed-right.fq.gz \
    forcetrimright=${trim_len}
    """
}

// Create library and check Hamming distance
process check_library {
  publishDir "$params.outdir/04_mapping/library", mode: 'copy'
  cpus 2

  input:
    file(library) from library_ch

  output:
    file "library-oligos.fa" into library_fa_ch1,library_fa_ch2
    file("library-oligos.txt") into library_txt
    file("library-oligos.aln") into library_aln
    file("*.png") into library_figs

  script:
    if( params.mode == 'library' )
      """
      # Verify the input table is correctly formatted (TODO)

      # Keep only the 1st and 2nd column
      cut -f 1,2 ${library} > library-oligos.txt

      # Convert table to FASTA
      seqkit tab2fx library-oligos.txt > library-oligos.fa

      # Run MUSCLE to create MSA
      muscle -in library-oligos.fa -out library-oligos.aln

      # Check the fasta library for distance issues
      check_library.py -i library-oligos.fa

      """
    else if ( params.mode == 'barcode' )
      """
      # Verify the input table is correctly formatted (TODO)

      # Keep only the 1st and 3rd column
      cut -f 1,3 ${library} > library-oligos.txt

      # Convert table to FASTA
      seqkit tab2fx library-oligos.txt > library-oligos.fa

      # Run MUSCLE to create MSA
      muscle -in library-oligos.fa -out library-oligos.aln

      # Check the fasta library
      check_library.py -i library-oligos.fa
      """
    else
        error "Invalid mode: ${params.mode}"
}

// Map reads directly to library
process mapping {
  publishDir "$params.outdir/04_mapping/bam", mode: 'copy', pattern: '*.bam'
  publishDir "$params.outdir/04_mapping/counts", mode: 'copy', pattern: '*-counts.txt'
  publishDir "$params.outdir/04_mapping/report", mode: 'copy', pattern: '*-report.txt'
  cpus 4

  input:
    file(library) from library_fa_ch1
    set sample_id, file(reads) from filter_fastq

  output:
    file("${sample_id}-counts.txt") into map_counts
    file("${sample_id}-mapped.bam") into map_bam
    file("${sample_id}-report.txt") into map_report
    file("${sample_id}-rpkm.txt") optional true into map_rpkm

  when:
    params.check == false | params.mode == "library"

  script:
    if( params.mapper == 'bowtie2' )
        """
        # Make small bowtie2 index
        bowtie2-build \
        --threads ${task.cpus} \
        ${library} \
        genome.index

        # Perform alignment
        bowtie2 \
        -x genome.index \
        -U ${reads} \
        -S ${sample_id}-mapped.sam \
        --very-sensitive \
        --norc 2> ${sample_id}-report.txt

        # Convert to BAM
        samtools view -S -f bam \
        -o ${sample_id}-mapped.bam \
        ${sample_id}-mapped.sam

        # Get counts from alignment
        pileup.sh \
        in=${sample_id}-mapped.bam \
        out=${sample_id}-counts.txt
        """
    else if( params.mapper == 'bbmap'  )
        """
        bbmap.sh \
        in=${reads} \
        out=${sample_id}-mapped.sam \
        outu=${sample_id}-unmapped.sam \
        ref=${library} \
        threads=${task.cpus} \
        statsfile=${sample_id}-report.txt \
        vslow \
        perfectmode \
        ambiguous=toss \
        -Xmx6g \
        nodisk

        # Convert to BAM
        samtools view -S -f bam \
        -o ${sample_id}-mapped.bam \
        ${sample_id}-mapped.sam

        # Get counts from alignment
        pileup.sh \
        in=${sample_id}-mapped.bam \
        out=${sample_id}-counts.txt
        """
    else
        error "Invalid alignment mode: ${params.mapper}"
}

// Merge the tables for downstream analysis
process extract_barcodes {
  publishDir "$params.outdir/03_extract/$sample_id", mode: 'copy'
  cpus 2

  input:
    file(library) from library_fa_ch2
    set sample_id, file(reads) from fastq_ch3

  output:
    set val(sample_id), file("${sample_id}-split.fq") into extract_fastq
    file("${sample_id}*.{html,zip}") into fastqc_extract

  when:
    params.check == false && params.mode == "barcode"

  script:
    """
    # Extract the last N sequences (barcodes)
    cat ${reads} | seqkit subseq -r -${params.length}:-1 > ${sample_id}-split.fq

    # Run FastQC on decontaminated data
    fastqc \
    --outdir . \
    --quiet \
    --threads ${task.cpus} \
    ${sample_id}-split.fq
    """
}

// Merge the tables for downstream analysis
process map_barcodes {
  publishDir "$params.outdir/04_mapping/bam", mode: 'copy', pattern: '*.bam'
  publishDir "$params.outdir/04_mapping/counts", mode: 'copy', pattern: '*-counts.txt'
  publishDir "$params.outdir/04_mapping/report", mode: 'copy', pattern: '*-report.txt'
  cpus 2

  input:
    file(library) from library_fa_ch2
    set sample_id, file(reads) from extract_fastq

  output:
    file("${sample_id}-counts.txt") into barcode_counts
    file("${sample_id}-mapped.bam") into barcode_bam
    file("${sample_id}-report.txt") into barcode_report

  when:
    params.check == false && params.mode == "barcode"

  script:
  if( params.mapper == 'bowtie2' )
      """
      # Make small bowtie2 index
      bowtie2-build \
      --threads ${task.cpus} \
      ${library} \
      genome.index

      # Perform alignment
      bowtie2 \
      -x genome.index \
      -U ${reads} \
      -S ${sample_id}-mapped.sam \
      --very-sensitive \
      --norc 2> ${sample_id}-report.txt

      # Convert to BAM
      samtools view -S -f bam \
      -o ${sample_id}-mapped.bam \
      ${sample_id}-mapped.sam

      # Get counts from alignment
      pileup.sh \
      in=${sample_id}-mapped.bam \
      out=${sample_id}-counts.txt
      """
  else if( params.mapper == 'bbmap'  )
      """
      bbmap.sh \
      in=${reads} \
      out=${sample_id}-mapped.sam \
      outu=${sample_id}-unmapped.sam \
      ref=${library} \
      threads=${task.cpus} \
      statsfile=${sample_id}-report.txt \
      perfectmode=${params.perfect} \
      ambiguous=toss \
      -Xmx6g \
      nodisk

      # Convert to BAM
      samtools view -S -f bam \
      -o ${sample_id}-mapped.bam \
      ${sample_id}-mapped.sam

      # Get counts from alignment
      pileup.sh \
      in=${sample_id}-mapped.bam \
      out=${sample_id}-counts.txt
      """
  else
      error "Invalid alignment mode: ${params.mapper}"
}

// Merge the tables for downstream analysis
process merge_tables_barcodes {
  publishDir "$params.outdir/04_mapping/", mode: 'copy'
  cpus 2

  input:
    file(tables) from barcode_counts.collect()

  output:
    file("*")

  when:
    params.check == false && params.mode == "barcode"

  script:
    """
    merge_tables.py \
    -i $tables
    """
}

// Merge the tables for downstream analysis
process merge_tables_library {
  publishDir "$params.outdir/04_mapping/", mode: 'copy'
  cpus 2

  input:
    file(tables) from map_counts.collect()

  output:
    file("*")

  when:
    params.check == false && params.mode == "library"

  script:
    """
    merge_tables.py \
    -i $tables
    """
}
