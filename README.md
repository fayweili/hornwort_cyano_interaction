# Programs and scripts needed
* [ccs](https://github.com/PacificBiosciences/ccs)
* [lima](https://github.com/PacificBiosciences/barcoding)
* [bam2fastq](https://gsl.hudsonalpha.org/information/software/bam2fastq)
* assign_sample.py


# Processing PacBio CCS reads
Make CCS reads from PacBio subread bam file, requiring 5 minimum passes and expected accuracy of >0.999
```bash
ccs --report-file ccs_report_p5_rq0.999.txt --minLength=100 --maxLength=2000 --num-threads=12 --min-passes=5 --min-rq=0.999 m54089_180212_172613.subreads.bam pilot_run_ccs_p5_rq0.999.bam 
```
Demultiplex the big bam file into separate bam files (each bam file now corresponds to a sample)
```bash
mkdir split_bams
cd split_bams
lima --different --ccs --num-threads 44 --min-length 100 --min-score 26 --split-bam-named ../pilot_run_ccs_p5_rq0.999.bam ../barcodes.fasta pilot_run_ccs_p5_rq0.999_demux.bam
```
Rename the bam files given the barcode combinations
```
assign_sample.py barcode-map.txt
```
* the `barcode-map.txt` file should look like this:
```
BCF1	BCR1	mock_community_5taxa
BCF2	BCR1	negative_control
BCF3	BCR1	G1_1_1_Notothylas
BCF4	BCR1	G1_1_5_soil
BCF5	BCR1	G1_2_1_Notothylas
BCF6	BCR1	G1_2_5_soil
BCF7	BCR1	P1_1_1_Phaeoceros
BCF8	BCR1	P1_1_5_Anthoceros
```
* also note that you might need to edit `assign_sample.py` so that it can find the correct bam files in your situation

