# Programs and scripts needed
* [ccs](https://github.com/PacificBiosciences/ccs)
* [lima](https://github.com/PacificBiosciences/barcoding)
* assign_sample.py


# Processing PacBio CCS reads
* Make CCS reads from PacBio subread bam file
```bash
ccs --report-file ccs_report_p5_rq0.999.txt --minLength=100 --maxLength=2000 --num-threads=12 --min-passes=5 --min-rq=0.999 m54089_180212_172613.subreads.bam pilot_run_ccs_p5_rq0.999.bam 
```
* Demultiplex the data into separate bam files (each bam file now corresponds to a sample)
```bash
mkdir split_bams
cd split_bams
lima --different --ccs --num-threads 44 --min-length 100 --min-score 26 --split-bam-named ../pilot_run_ccs_p5_rq0.999.bam ../barcodes.fasta pilot_run_ccs_p5_rq0.999_demux.bam
```
* Remove the bam files given the barcode combinations
```
assign_sample.py
```
