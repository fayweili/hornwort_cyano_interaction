library(phyloseq)
library(dada2)
library(ggplot2)
library(Biostrings)

## Reading files ##
path <- "/Users/fay-weili/Desktop/Hornwort_amplicon/dada2/sample_fastq" # CHANGE ME to location of the fastq file
fns <- list.files(path, pattern="fastq.gz", full.names=TRUE)
#fn <- file.path(path, "JNP2_6_F01_ccs_demux.BCF1--BCR1.fastq.gz")
sample.names <- sapply(strsplit(basename(fns), ".fastq.gz"), '[', 1) # get sample names from fastq file names
cw <- "CGTAGCTTCCGGTGGTATCCACGT" # primer 1
cx <- "GGGGCAGGTAAGAAAGGGTTTCGTA" # primer 2
rc <- dada2:::rc 
theme_set(theme_bw())

## Remove primers ##
nop <- file.path(paste(path,'_deprimers', sep=''), basename(fns))
prim <- removePrimers(fns, nop, primer.fwd=cw, primer.rev=dada2:::rc(cx), orient=TRUE, verbose=TRUE)
lens.fn <- lapply(nop, function(fn) nchar(getSequences(fn))) # plot len distribution
lens <- do.call(c, lens.fn)
hist(lens, 100)

## Filter ##
filt <- file.path(paste(path,'_deprimers_lenfiltered', sep=''), basename(fns))
track <- filterAndTrim(nop, filt, minQ=3, minLen=400, maxLen=1200, maxN=0, rm.phix=FALSE, maxEE=2, verbose=TRUE)

## DADA2 ##
drp <- derepFastq(filt, verbose=TRUE, qualityType="FastqQuality") # dereplicate
err <- learnErrors(drp, BAND_SIZE=32, multithread=TRUE, errorEstimationFunction=PacBioErrfun) # learn error
plotErrors(err)
dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sapply(dd, function(x) sum(x$denoised)))
seqtab <- makeSequenceTable(dd); dim(seqtab)

## Save seqtab file ##
saveRDS(seqtab, 'rds')

## Check chimera ##
bim <- isBimeraDenovo(seqtab, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## Read metadata ##
#sample_df <- read.table("/Users/fayweili/Desktop/Hornwort_microbiome/dada2_test/TimeSeriesMeta_lib1.csv", header=TRUE, row.names=1, sep=",", stringsAsFactors=FALSE)

## Make phyloseq object ##
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE)) #,sample_data(sample_df))
ps <- prune_samples(sample_names(ps) != "mock_community_5taxa_JNP1_5_E01.fastq.gz", ps)
ps <- prune_samples(sample_names(ps) != "negative_control_JNP1_5_E01.fastq.gz", ps)

## Rename taxa names from sequences to ASV ##
dna <- DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
writeXStringSet((refseq(ps)), 'ASV.fa') # write sequences to fasta

## Filter out non-rbcLX ASVs based on usearch (ran elsewhere) ##
ASV_off_target <- c("ASV1182", "ASV1542", "ASV1614", "ASV1859", "ASV1858", "ASV1278", "ASV1915", "ASV1761", "ASV1760", "ASV1833", "ASV1911", "ASV1912", "ASV1609", "ASV1919", "ASV1252", "ASV1850", "ASV1525", "ASV1527", "ASV1335", "ASV1905", "ASV1409", "ASV1490", "ASV1781", "ASV1393", "ASV1783", "ASV853", "ASV1866", "ASV1541", "ASV1502", "ASV1687", "ASV1568", "ASV1457", "ASV1383", "ASV1907", "ASV1693", "ASV1826", "ASV1605", "ASV1638", "ASV1340", "ASV1779", "ASV1432", "ASV1940", "ASV1828", "ASV1829", "ASV1909", "ASV1867", "ASV1906", "ASV1458", "ASV1825", "ASV1279", "ASV1047", "ASV1606", "ASV1455", "ASV1712", "ASV1910", "ASV1341", "ASV1692", "ASV1218")
ps_ontarget <- prune_taxa(!(taxa_names(ps) %in% ASV_off_target), ps)
OTU = as(otu_table(ps_ontarget), "matrix")
OTUdf = as.data.frame(OTU)
write.csv(OTUdf, "ASV_on_target_table.csv") # write OTU table

## Transform otu table
#ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
#ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
#plot_ordination(ps.prop, ord.nmds.bray, color="Site", title="Bray NMDS")
#plot_ordination(ps.prop, ord.nmds.bray, color="Quadrat", title="Bray NMDS")
