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
-resume

# Run check to get initial quality results (barcode 17)
nextflow run main.nf \
--fastq 'testing/barcode-n17/*.fq' \
--library 'testing/barcode-n17/library.txt' \
--length 17 \
--mode 'barcode' \
--outdir 'res-barcode-17' \
-resume

# Run check to get initial quality results (barcode 17)
nextflow run main.nf \
--fastq 'testing/barcode-n12/*.fq' \
--library 'testing/barcode-n12/library.txt' \
--length 12 \
--mode 'barcode' \
--outdir 'res-barcode-12' \
-resume
