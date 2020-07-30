# Amplicon-seq: A package for quanitifying amplicon runs (e.g antigen drop out or TCR profiling)
# @ Thomas W. Battaglia
# 29-07-2029

# Run check to get initial quality results
nextflow run main.nf \
--fastq 'testing/library/*.fq' \
--library 'testing/library/library.txt' \
--length 63 \
--mode 'library' \
--outdir 'res-library' \
-resume \
--check

# Run check to get initial quality results
nextflow run main.nf \
--fastq 'testing/barcode/*.fq' \
--library 'testing/barcode/library.txt' \
--length 17 \
--mismatch 2 \
--mode 'barcode' \
--outdir 'res-barcode' \
-resume
