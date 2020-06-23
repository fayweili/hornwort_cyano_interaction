# Programs and scripts needed
* [ccs](https://github.com/PacificBiosciences/ccs)
* [lima](https://github.com/PacificBiosciences/barcoding)
* [bam2fastq](https://gsl.hudsonalpha.org/information/software/bam2fastq)
* [assign_sample.py](https://github.com/fayweili/hornwort_cyano_interaction/blob/master/scripts/assign_sample.py)


# Processing PacBio CCS reads
**Make CCS reads from PacBio subread bam file, requiring 5 minimum passes and expected accuracy of >0.999**
```bash
ccs --report-file ccs_report_p5_rq0.999.txt --minLength=100 --maxLength=2000 --num-threads=12 --min-passes=5 --min-rq=0.999 m54089_180212_172613.subreads.bam pilot_run_ccs_p5_rq0.999.bam 
```
**Demultiplex the big bam file into separate bam files (each bam file now corresponds to a unique sample)**
```bash
mkdir split_bams
cd split_bams
lima --different --ccs --num-threads 44 --min-length 100 --min-score 26 --split-bam-named ../pilot_run_ccs_p5_rq0.999.bam ../barcodes.fasta pilot_run_ccs_p5_rq0.999_demux.bam
```
**Rename the bam files given the barcode combinations**
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
**Remove primer sequences and filter by quality and length using dada2**
```R
## Loading libraries ##
library(phyloseq)
library(dada2)
library(ggplot2)
library(Biostrings)

## Reading files ##
path <- "/Users/fay-weili/Desktop/Hornwort_amplicon/dada2/sample_fastq" # CHANGE ME to location of the fastq file
fns <- list.files(path, pattern="fastq.gz", full.names=TRUE)
sample.names <- sapply(strsplit(basename(fns), ".fastq.gz"), '[', 1) # get sample names from fastq file names
cw <- "CGTAGCTTCCGGTGGTATCCACGT" # primer 1
cx <- "GGGGCAGGTAAGAAAGGGTTTCGTA" # primer 2
rc <- dada2:::rc 

## Remove primers ##
nop <- file.path(paste(path,'_deprimers', sep=''), basename(fns))
prim <- removePrimers(fns, nop, primer.fwd=cw, primer.rev=dada2:::rc(cx), orient=TRUE, verbose=TRUE)
lens.fn <- lapply(nop, function(fn) nchar(getSequences(fn))) # plot len distribution
lens <- do.call(c, lens.fn)
hist(lens, 100)

## Filter ##
filt <- file.path(paste(path,'_deprimers_lenfiltered', sep=''), basename(fns))
track <- filterAndTrim(nop, filt, minQ=3, minLen=400, maxLen=1200, maxN=0, rm.phix=FALSE, maxEE=2, verbose=TRUE)
```
