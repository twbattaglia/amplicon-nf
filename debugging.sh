# Amplicon-seq: A package for quanitifying amplicon runs (e.g antigen drop out or TCR profiling)
# @ Thomas W. Battaglia
# 29-07-2029

# Run check to get initial quality results
nextflow run main.nf \
--fastq 'testing/fastq/*.fq' \
--library 'testing/chiara-library.txt' \
--length 63 \
-resume \
--check

# Run full pipeline with library-oriented
nextflow run main.nf \
--fastq 'testing/fastq/*.fq' \
--length 63 \
-resume
