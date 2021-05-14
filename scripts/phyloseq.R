library(phyloseq)
library(ape)
library(ggplot2)
library(plyr)
library(DESeq2)
library(ggbeeswarm)
library(patchwork)
library(vegan)
library(dplyr)
library(tidyr)
library(microbiome)
library(RVAideMemoire)
library(exactRankTests)
library(compositions)
library(ranacapa)
library(picante)
library(UpSetR)
library(cowplot)
library(ggplotify)
library(reshape2)

# Use these when plotting the trees #
library(ggtree)
library(treeio)
library(naniar)
library(phytools)

source("/Users/fay-weili/bin/ANCOM/scripts/ancom_v2.1.R")
theme_set(theme_bw())
cbp1 <- c("#E69F00", "#56B4E9", "#009E73","#999999", 
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_ASV_deseq <- function(gene_to_plot) {
  geneCounts <- plotCounts(dds, gene = "ASV1", intgroup = 'Type', returnData = TRUE)
  ggplot(geneCounts, aes(x = Type, y = count, color = Type)) + geom_boxplot() + 
    geom_beeswarm(cex = 1, size=3) + labs(title = gene_to_plot) + 
    scale_colour_manual(values = cbp1) 
}

plot_ASV <- function(ASV_to_plot, phyloseq_obj, violin=FALSE) {
  #ASV_to_plot <- "ASV1"
  #phyloseq_obj <- potato3_physeq_transformed_plantsonly
  ASV_t_subset <- t(otu_table(phyloseq_obj))[,ASV_to_plot]
  ASV_count <- cbind(ASV_t_subset, sample_data(phyloseq_obj))
  if (violin==FALSE) {
    ggplot(ASV_count, aes_string(x = "Type", y = ASV_to_plot, color = "Type")) + geom_boxplot() + 
      geom_beeswarm(cex = 1, size=3) + labs(title = ASV_to_plot, y="Relative abundance") + 
      scale_colour_manual(values = cbp1) }
  else {
    ggplot(ASV_count, aes_string(x = "Type", y = ASV_to_plot)) + 
      geom_violin(aes_string(fill="Type"), scale="width") + geom_boxplot(width=0.1) + 
      labs(title = ASV_to_plot, y="Relative abundance") + 
      scale_colour_manual(values = cbp1) + scale_fill_manual(values = cbp1) }
}

plot_ASV_soil <- function(ASV_to_plot, phyloseq_obj, violin=FALSE) {
  #ASV_to_plot <- "ASV8"
  #phyloseq_obj <- potato_physeq_transformed
  ASV_t_subset <- t(otu_table(phyloseq_obj))[,ASV_to_plot]
  ASV_count <- cbind(ASV_t_subset, sample_data(phyloseq_obj))
  ASV_count <- filter(ASV_count, Type=="soil")
  if (violin==FALSE) {
    ggplot(ASV_count, aes_string(x = "Quadrat", y = ASV_to_plot)) + geom_boxplot(color = "#999999") + 
      geom_beeswarm(cex = 1, size=3, color = "#999999") + labs(title = ASV_to_plot, y="Relative abundance") }
  else {  
    ggplot(ASV_count, aes_string(x = "Quadrat", y = ASV_to_plot)) + 
      geom_violin(fill="#999999", scale="width") + geom_boxplot(width=0.1) + 
      labs(title = ASV_to_plot, y="Relative abundance") }
}

phyloseq2upsetr <- function(phyloseq_obj, category){
  #phyloseq_obj <- time_series_physeq_plant_filtered
  #category <- "Quadrat"
  sampleda <- sample_data(phyloseq_obj)[,match(category, colnames(sample_data(phyloseq_obj))),drop=FALSE]
  asv_meta <- merge(t(otu_table(phyloseq_obj)), sampleda, by=0)
  asv_meta <- tibble::column_to_rownames(asv_meta, var="Row.names")
  ### group by Type
  asv_meta <- data.frame(plyr::ddply(asv_meta, category, plyr::numcolwise(sum)), 
                         check.names=FALSE, stringsAsFactors=FALSE)   
  rownames(asv_meta) <- as.vector(asv_meta[[category]])
  asv_meta[[category]] <- NULL
  ### turn to absence/present
  upset_data <- apply(asv_meta, 1, 
                      function(x){unlist(lapply(x, function(x){if(x>0){1}else{0}}))})
  upset_data <- data.frame(upset_data, check.names=FALSE)
  return(upset_data)
  }

phyloseq2core <- function(phyloseq_obj, type, quadrat, time){
  #phyloseq_obj <- potato_physeq_plantsonly_transformed
  #quadrat <- 'Grossman1'
  #time <- 'T1'
  #type <- 'Noto'
  #phyloseq_obj <- subset_samples(grossman2_physeq_transformed_Noto, Time=='T1')
  tb <- otu_table(phyloseq_obj) %>% as("matrix") %>% tibble::as_tibble(rownames = "ASV")
  nonzero <- function(x) sum(x != 0)
  overhalf <- function(x) sum(x > 0.5)
  countall <- function(x) sum(x >= 0)
  rtb <- tb %>% rowwise(ASV)
  table <- rtb %>% mutate(colno = countall(c_across(where(is.numeric))), total = (nonzero(c_across(where(is.numeric))) - 1)/colno) %>% 
    select(total) %>% arrange(desc(total)) %>% filter(total >= 0.5) 
  table <- table %>% tibble::add_column(Type=type, Quadrat=quadrat, Time=time)
  if (dim(table)[1] != 0) {
    asv_list <- pull(table,ASV)
    abun_tb <- tb %>% filter(ASV %in% asv_list) %>% tibble::column_to_rownames(var='ASV')
    #pheatmap(abun_tb,cluster_rows=FALSE, cluster_cols=FALSE)
    colnames(abun_tb) <- c(1,2,3,4)
    rownames(abun_tb) <- paste(rownames(abun_tb), quadrat, type, time, sep = "_")
    abun_tb_melt <- reshape::melt(as.matrix(abun_tb),value.name="abundance",varnames=c("ASV","replicate"))
    order <- rownames(abun_tb)
    abun_tb_melt$ASV <- factor(abun_tb_melt$ASV, levels=order)
    abun_tb_melt <- abun_tb_melt %>% tibble::add_column(Type=type, Quadrat=quadrat, Time=time)
    #ggplot(data=abun_tb_melt, aes(x=ASV, y=replicate, fill=value)) + geom_tile()
    return(abun_tb_melt)
    }
  else {
    print('null')
    }
  #ggplot(data=abun_tb_melt) + geom_point(aes(y=ASV, x=abundance))
  
}

tbl2hmap= function(tbl.long, rowVar, colVar, valueVar,
                   colAnnVars = NULL, rowAnnVars = NULL) {
  Mat0 = dplyr::select(tbl.long, one_of(c(rowVar, colVar,valueVar)))%>%
    tidyr::spread_(key = colVar, value=valueVar)
  Mat = select(Mat0, -one_of(rowVar)) %>% as.matrix()
  rownames(Mat) = unlist(select(Mat0,one_of(rowVar)))
  
  colAnn = rowAnn= NULL
  
  if(!is.null(colAnnVars)) {
    colAnn0 = dplyr::select(tbl.long,one_of( c(colVar,colAnnVars))) %>%
      unique()
    colAnn = select(colAnn0, - one_of(colVar)) %>% as.data.frame()
    rownames(colAnn) = unlist(select( colAnn0,one_of(colVar)))
  }
  
  if(!is.null(rowAnnVars)) {
    rowAnn0 = dplyr::select(tbl.long,one_of(c( rowVar,rowAnnVars))) %>%
      unique()
    rowAnn = select(rowAnn0, -one_of( rowVar)) %>% as.data.frame()
    rownames(rowAnn) = unlist(select(rowAnn0,one_of(rowVar)))
  }
  
  list(mat = Mat, rowAnn =rowAnn, colAnn = colAnn)
}

filter5perc_count = function(x){
  x[(x / sum(x)) < (0.05)] <- 0
  return(x) }  
filter3perc_count = function(x){
  x[(x / sum(x)) < (0.03)] <- 0
  return(x) } 

plotDistances = function(p = grossman_physeq_transformed, m = "wunifrac", var1 = "Time", var2 = "Quadrat", var3 = "Type") {
  # calc distances
  wu = phyloseq::distance(p, m)
  wu.m = melt(as.matrix(wu))
  
  # remove self-comparisons
  wu.m = wu.m %>%
    filter(as.character(Var1) != as.character(Var2)) %>%
    mutate_if(is.factor,as.character)
  colnames(wu.m) = c("sample1", "sample2", "value")
  
  # get sample data (S4 error OK and expected)
  sd = as_tibble(sample_data(p),rownames="sample") %>% 
    select(sample, var1, var2, var3) %>% mutate_if(is.factor,as.character)

  # combined distances with sample data
  colnames(sd) = c("sample1", "Var1a", "Var2a", "Var3a")
  wu.sd = left_join(wu.m, sd, by = "sample1")
  
  colnames(sd) = c("sample2", "Var1b", "Var2b", "Var3b")
  wu.sd = left_join(wu.sd, sd, by = "sample2")
  
  wu.sd = wu.sd %>% filter(as.character(Var1a) == as.character(Var1b)) %>% 
    filter(as.character(Var2a) == as.character(Var2b)) %>%
    filter(as.character(Var3a) == as.character(Var3b)) %>%
    mutate_if(is.character, stringr::str_replace_all, pattern = "T", replacement = "")
  
  # plot 
  ggplot(wu.sd, aes(x=Var1a, y=value, colour=Var3a)) + 
    geom_point() + 
    geom_smooth(aes(x=as.numeric(Var1a), y=value),method = "loess",se=T) + 
    facet_grid(Var3a ~ Var2a) + scale_colour_manual(values = cbp1) +
    xlab("Time points") + ylab("Pairwise weighted unifrac distance") 
  #linearMod <- lm(value ~ as.numeric(Var1a), data=wu.sd)
  #summary(linearMod)
 }
grossman_time_variance <- plotDistances(p = grossman_physeq_transformed, m = "wunifrac", var1 = "Time", var2 = "Quadrat", var3 = "Type")
potato_time_variance <- plotDistances(p = potato_physeq_transformed, m = "wunifrac", var1 = "Time", var2 = "Quadrat", var3 = "Type")
(grossman_time_variance / plot_spacer()) | potato_time_variance
grossman_time_variance | potato_time_variance
ggsave('unifrac_througthT.pdf', device = "pdf", width = 10, height = 5)
ggsave('unifrac_througthT.svg', device = "svg", width = 10, height = 5)

setwd("/Users/fay-weili/Box/hornwort_amplicon/dada2/")

# Preparing PhyloSeq  -----------------------------------------------------
  ## Read ASV table 
    ASV_file <- as.matrix(read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/ASV_Table_timeseries_cleaned.txt", header=TRUE, sep = "\t", row.names = 1))
  ## Read sample metadata
    sample_meta_file <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/TimeSeriesMeta.csv", header=TRUE, sep = ",", row.names = 1)
  ## Read the tree file
    tree <- read_tree("/Users/fay-weili/Box/hornwort_amplicon/dada2/ASV_on_target_timeseries_pasta_rooted.tre")
  ## Make phyloseq object
    ASV_table <- otu_table(ASV_file, taxa_are_rows = TRUE)
    sample_meta <- sample_data(sample_meta_file)
    time_series_physeq = phyloseq(ASV_table, sample_meta, tree)
    time_series_physeq_transformed <- transform_sample_counts(time_series_physeq, function(x) x/sum(x))
    time_series_physeq_soil_transformed <- subset_samples(time_series_physeq_transformed, Type=="soil")
    time_series_physeq_plant_transformed <- subset_samples(time_series_physeq_transformed, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
    time_series_physeq_soil <- subset_samples(time_series_physeq, Type=="soil")
    time_series_physeq_soil_filtered <- transform_sample_counts(time_series_physeq_soil, fun = filter3perc_count)  
    time_series_physeq_plant <- subset_samples(time_series_physeq, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
    time_series_physeq_plant_filtered <- transform_sample_counts(time_series_physeq_plant, fun = filter3perc_count)  
    time_series_physeq_filtered <- merge_phyloseq(time_series_physeq_soil_filtered, time_series_physeq_plant_filtered)
    time_series_physeq_filtered_transformed <- transform_sample_counts(time_series_physeq_filtered, function(x) x/sum(x))
    
    time_series_physeq_soil_filtered_transformed <- transform_sample_counts(time_series_physeq_soil_filtered, function(x) x/sum(x))
    tb <- otu_table(time_series_physeq_soil_filtered_transformed) %>% as("matrix") %>% tibble::as_tibble(rownames = "ASV")
    nonzero <- function(x) sum(x != 0)
    x <- numcolwise(nonzero)(tb)
    mean(t(x))
    
    time_series_physeq_plant_filtered_transformed <- transform_sample_counts(time_series_physeq_plant_filtered, function(x) x/sum(x))  
    time_series_physeq_filtered_all <- transform_sample_counts(time_series_physeq, fun = filter3perc_count)  
    time_series_physeq_filtered_all_transformed <- transform_sample_counts(time_series_physeq_filtered_all, function(x) x/sum(x))  
    #write.table(otu_table(time_series_physeq_plant_filtered_transformed), quote=F, file="time_series_physeq_plant_filtered_transformed.txt", sep="\t")
    #write.table(otu_table(time_series_physeq_filtered_all_transformed), quote=F, file="time_series_physeq_filtered_all_transformed.txt", sep="\t")
  
  ## UpsetR
    p_upset_quadrat <- upset(phyloseq2upsetr(time_series_physeq_filtered, "Quadrat"),order.by = "freq")
    P_upset_potato_type <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Potato1"|Quadrat=="Potato2"|Quadrat=="Potato3"), "Type"),order.by = "freq")
    P_upset_grossman_type <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Grossman1"|Quadrat=="Grossman2"), "Type"),order.by = "freq")
    P_upset_all_type <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Grossman1"|Quadrat=="Grossman2"|Quadrat=="Potato1"|Quadrat=="Potato2"|Quadrat=="Potato3"), "Type"),order.by = "freq")
    upset_df <- phyloseq2upsetr(time_series_physeq_filtered, "Type")
    write.table(phyloseq2upsetr(time_series_physeq_filtered, "Quadrat"), "ASV_by_Q.txt", quote=F, sep="\t")
    write.table(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Potato1"|Quadrat=="Potato2"|Quadrat=="Potato3"), "Type"), "ASV_by_Type_PH.txt", quote=F, sep="\t")
    write.table(otu_table(time_series_physeq_filtered_transformed), "ASV_table_for_Jessica.txt", quote=F, sep="\t")
    
    #no singleton i.e. all ASVs appear in more than one sample
    time_series_physeq_filtered_nosingleton <- filter_taxa(time_series_physeq_filtered, function(x){sum(x > 0) > 1}, prune = TRUE)
    p_upset_quadrat <- upset(phyloseq2upsetr(time_series_physeq_filtered_nosingleton, "Quadrat"),order.by = "freq",
                  sets = c("Grossman2","Grossman1","Potato3","Potato2","Potato1"),keep.order = TRUE)
    P_upset_potato_type <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered_nosingleton, Quadrat=="Potato1"|Quadrat=="Potato2"|Quadrat=="Potato3"), "Type"),order.by = "freq",
                  sets = c("soil","Notothylas","Phaeoceros","Anthoceros"),keep.order = TRUE)
    P_upset_grossman_type <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered_nosingleton, Quadrat=="Grossman1"|Quadrat=="Grossman2"), "Type"),order.by = "freq")
    
    ### plots
    library(grid)
    library(svglite)
    p_upset_quadrat 
    grid.edit('arrange',name='arrange2')
    vp1 = grid.grab()
    P_upset_potato_type 
    grid.edit('arrange',name='arrange2')
    vp2 = grid.grab()
    #P_upset_grossman_type
    #grid.edit('arrange',name='arrange2')
    #vp3 = grid.grab()
    svglite("upset_plot_nosingleton.svg", width = 9, height = 3.5)
    plot_grid(vp1,vp2, ncol = 2, labels = c('A', 'B'))
    dev.off()
    
    ### Time points
    P_upset_potato_Q1_all_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Potato1"), "Time"),order.by = "freq",
                                        sets = c("T4","T3","T2","T1"),keep.order = TRUE)
    P_upset_potato_Q1_plant_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_plant_filtered, Quadrat=="Potato1"), "Time"),order.by = "freq",
                                          sets = c("T4","T3","T2","T1"),keep.order = TRUE)
    P_upset_potato_Q1_soil_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_soil, Quadrat=="Potato1"), "Time"),order.by = "freq",
                                         sets = c("T4","T3","T2","T1"),keep.order = TRUE)
    P_upset_potato_Q2_all_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Potato2"), "Time"),order.by = "freq",
                                        sets = c("T4","T3","T2","T1"),keep.order = TRUE)
    P_upset_potato_Q3_all_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Potato3"), "Time"),order.by = "freq",
                                        sets = c("T4","T3","T2","T1"),keep.order = TRUE)
    P_upset_grossman_Q1_all_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Grossman1"), "Time"),order.by = "freq",
                                          sets = c("T4","T3","T2","T1"),keep.order = TRUE)
    P_upset_grossman_Q2_all_time <- upset(phyloseq2upsetr(subset_samples(time_series_physeq_filtered, Quadrat=="Grossman2"), "Time"),order.by = "freq",
                                          sets = c("T3","T2","T1"),keep.order = TRUE)
    
    P_upset_grossman_Q1_all_time 
    grid.edit('arrange',name='arrange2')
    vp1 = grid.grab()
    P_upset_grossman_Q2_all_time 
    grid.edit('arrange',name='arrange2')
    vp2 = grid.grab()
    svglite("upset_plot_GrossmanQ_by_time_filtered.svg", width = 9, height = 5)
    plot_grid(vp1,vp2, ncol = 2, labels = c('GrossmanQ1', 'GrossmanQ2'))
    dev.off()
 
    P_upset_potato_Q1_all_time 
    grid.edit('arrange',name='arrange2')
    vp12 = grid.grab()
    P_upset_potato_Q2_all_time 
    grid.edit('arrange',name='arrange2')
    vp22 = grid.grab()
    P_upset_potato_Q3_all_time 
    grid.edit('arrange',name='arrange2')
    vp32 = grid.grab()
    svglite("upset_plot_PotatoQ_by_time_filtered.svg", width = 9, height = 9)
    plot_grid(vp12,vp22,vp32,vp1,vp2, ncol = 3)
    dev.off()   
    
  ## PCoA of everything
    ### Color by Type
      time_series_physeq_transformed.ord <- ordinate(time_series_physeq_filtered_transformed, "PCoA", "bray")
      time_series_pcoa_all <- plot_ordination(time_series_physeq_filtered_transformed, time_series_physeq_transformed.ord, type="samples", color="Quadrat") 
      time_series_pcoa_all_DF <- plot_ordination(time_series_physeq_filtered_transformed, time_series_physeq_transformed.ord, type="samples", color="Quadrat", justDF=TRUE) 
      time_series_pcoa_all_soil <- ggplot(data = subset(time_series_pcoa_all_DF, Type=='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, shape=17, aes(color=Quadrat)) + geom_point(shape = 2,size = 1.8,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
        coord_fixed() + xlim(-0.52, 0.52)+ ylim(-0.4, 0.5) + ggtitle("Soil samples")
      time_series_pcoa_all_plant <- ggplot(data = subset(time_series_pcoa_all_DF, Type!='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, aes(color=Quadrat)) + geom_point(shape = 1,size = 1.8,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
        coord_fixed()+ xlim(-0.52, 0.52) + ylim(-0.4, 0.5) + ggtitle("Plant samples")
      time_series_pcoa_all_soil + time_series_pcoa_all_plant + plot_layout(guides = 'collect')
      ggsave("potato_grossman_pcoa_plants_soils_allinone.pdf", device = "pdf", width = 10, height = 5)
      ggsave("potato_grossman_pcoa_plants_soils_allinone.svg", device = "svg", width = 10, height = 5)
      #### unifrac
      time_series_physeq_transformed.ord <- ordinate(time_series_physeq_filtered_transformed, "PCoA", "wunifrac")
      time_series_pcoa_all <- plot_ordination(time_series_physeq_filtered_transformed, time_series_physeq_transformed.ord, type="samples", color="Quadrat") 
      time_series_pcoa_all_DF <- plot_ordination(time_series_physeq_filtered_transformed, time_series_physeq_transformed.ord, type="samples", color="Quadrat", justDF=TRUE) 
      time_series_pcoa_all_soil <- ggplot(data = subset(time_series_pcoa_all_DF, Type=='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, shape=17, aes(color=Quadrat)) + geom_point(shape = 2,size = 1.8,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
        xlim(-0.13, 0.175)+ ylim(-0.125, 0.08) + ggtitle("Soil samples")
      time_series_pcoa_all_plant <- ggplot(data = subset(time_series_pcoa_all_DF, Type!='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, aes(color=Quadrat)) + geom_point(shape = 1,size = 1.8,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
        xlim(-0.13, 0.175)+ ylim(-0.125, 0.08) + ggtitle("Plant samples")
      time_series_pcoa_all_soil + time_series_pcoa_all_plant + plot_layout(guides = 'collect')
      ggsave("potato_grossman_pcoa_plants_soils_allinone_unifrac.pdf", device = "pdf", width = 12, height = 3)
      ggsave("potato_grossman_pcoa_plants_soils_allinone_unifrac.svg", device = "svg", width = 12, height = 3)
      
      time_series_pcoa_all <- ggplot(data = time_series_pcoa_all_DF, mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, aes(color=Quadrat, shape=)) + geom_point(shape = 1,size = 1.8,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
        coord_fixed()+ xlim(-0.13, 0.175)+ ylim(-0.125, 0.08) + ggtitle("Plant samples")
      
    ### Color by Time
      time_series_physeq_transformed.ord <- ordinate(time_series_physeq_transformed, "PCoA", "bray")
      time_series_pcoa_all_DF <- plot_ordination(time_series_physeq_transformed, time_series_physeq_transformed.ord, type="samples", color="Time", justDF=TRUE) 
      time_series_pcoa_all_soil <- ggplot(data = subset(time_series_pcoa_all_DF, Type=='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, shape=17, aes(color=Time)) + geom_point(shape = 2,size = 1.8,colour = "black") + scale_color_brewer(palette="YlGnBu") +
        coord_fixed() + xlim(-0.52, 0.52)+ ylim(-0.4, 0.5) + ggtitle("Soil samples")
      time_series_pcoa_all_plant <- ggplot(data = subset(time_series_pcoa_all_DF, Type!='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
        geom_point(size=1.8, aes(color=Time)) + geom_point(shape = 1,size = 1.8,colour = "black") + scale_color_brewer(palette="YlGnBu")+
        coord_fixed()+ xlim(-0.52, 0.52) + ylim(-0.4, 0.5) + ggtitle("Plant samples")
      time_series_pcoa_all_soil + time_series_pcoa_all_plant + plot_layout(guides = 'collect')
      
      time_series_bray_dist <- phyloseq::distance(time_series_physeq_transformed, method="bray")
      pairwise_dist <- dist2list(time_series_bray_dist)
      ggplot(pairwise_dist) + geom_histogram(aes(value))
  ## PERMANOVA 
    ## All 
      time_series_wunifrac_dist <- phyloseq::distance(time_series_physeq_transformed, method="wunifrac")
      
      simplified_type <- sample_data(time_series_physeq_transformed)$Type
      levels(simplified_type)<-c(levels(simplified_type),"plant")
      simplified_type[simplified_type=="Notothylas"]<-"plant"
      simplified_type[simplified_type=="Anthoceros"]<-"plant"
      simplified_type[simplified_type=="Phaeoceros"]<-"plant"
      
      site <- sample_data(time_series_physeq_transformed)$Quadrat
      levels(site)<-c(levels(site),"Potato","Grossman")
      site[site=="Potato1"]<-"Potato"
      site[site=="Potato2"]<-"Potato"
      site[site=="Potato3"]<-"Potato"
      site[site=="Grossman1"]<-"Grossman"
      site[site=="Grossman2"]<-"Grossman"
      
      adonis(time_series_wunifrac_dist ~ sample_data(time_series_physeq_transformed)$Quadrat+simplified_type, permutations = 10000)
      adonis(time_series_wunifrac_dist ~ sample_data(time_series_physeq_transformed)$Quadrat+sample_data(time_series_physeq_transformed)$Time+simplified_type, permutations = 10000)
      adonis(time_series_wunifrac_dist ~ sample_data(time_series_physeq_transformed)$Time, permutations = 10000)
      adonis(time_series_wunifrac_dist ~ site, permutations = 10000)
      adonis(time_series_wunifrac_dist ~ site + sample_data(time_series_physeq_transformed)$Quadrat + sample_data(time_series_physeq_transformed)$Time + simplified_type, permutations = 10000)
      
      beta <- betadisper(time_series_bray_dist, sample_data(time_series_physeq_transformed)$Type)
      permutest(beta)
    
    ## Potato
      potato_wunifrac_dist <- phyloseq::distance(potato_physeq_transformed, method="wunifrac")
      simplified_type <- sample_data(potato_physeq_transformed)$Type
      levels(simplified_type)<-c(levels(simplified_type),"plant")
      simplified_type[simplified_type=="Notothylas"]<-"plant"
      simplified_type[simplified_type=="Anthoceros"]<-"plant"
      simplified_type[simplified_type=="Phaeoceros"]<-"plant"
      adonis(potato_wunifrac_dist ~ sample_data(potato_physeq_transformed)$Quadrat+sample_data(potato_physeq_transformed)$Time+simplified_type, permutations = 10000)
  
    ## Grossman
      grossman_wunifrac_dist <- phyloseq::distance(grossman_physeq_transformed, method="wunifrac")
      simplified_type <- sample_data(grossman_physeq_transformed)$Type
      levels(simplified_type)<-c(levels(simplified_type),"plant")
      simplified_type[simplified_type=="Notothylas"]<-"plant"
      adonis(grossman_wunifrac_dist ~ sample_data(grossman_physeq_transformed)$Quadrat+sample_data(grossman_physeq_transformed)$Time+simplified_type, permutations = 10000)
      
  ### test plot
    time_series_pcoa_all_Grossman_soil <- ggplot(data = subset(time_series_pcoa_all_DF, Type=='soil'&(Quadrat=='Grossman1'|Quadrat=='Grossman2')), mapping = aes(x = Axis.1, y = Axis.2)) +
      geom_point(size=3, shape=17, aes(color=Quadrat)) + geom_point(shape = 2,size = 3,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
      coord_fixed()+ xlim(-0.25, -0.05) + ylim(0, 0.26) + ggtitle("Soil samples")
    time_series_pcoa_all_Grossman_plant <- ggplot(data = subset(time_series_pcoa_all_DF, Type!='soil'&(Quadrat=='Grossman1'|Quadrat=='Grossman2')), mapping = aes(x = Axis.1, y = Axis.2)) +
      geom_point(size=3, aes(color=Quadrat)) + geom_point(shape = 1,size = 3,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3","steelblue","gold2","darkseagreen")) +
      coord_fixed()+ xlim(-0.25, -0.05) + ylim(0, 0.26) + ggtitle("Plant samples")
    time_series_pcoa_all_Grossman_soil + time_series_pcoa_all_Grossman_plant + plot_layout(guides='collect')
  ## Network
    ig <- make_network(time_series_physeq_transformed, distance="bray", max.dist=0.4)
    net_all <- plot_network(ig, time_series_physeq_transformed, color = "Time", label=NULL, point_size=4) +
      geom_point(shape = 1,size = 4,colour = "black") #+ scale_colour_manual(values = cbp1) #+ coord_fixed()
    
    
# Plot Tree ---------------------------------------------------------------
  ## Get dataset ready
  ## for plants, remove ASV relative abundance <5%
  ## for soils, keep as is
    time_series_physeq_plantsonly <- subset_samples(time_series_physeq, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
    time_series_physeq_plantsonly_filtered <- transform_sample_counts(time_series_physeq_plantsonly, fun = filter5perc_count)  
    time_series_physeq_plantsonly_filtered_transformed <- transform_sample_counts(time_series_physeq_plantsonly_filtered, function(x) x / sum(x) )
    time_series_physeq_filtered_transformed <- merge_phyloseq(time_series_physeq_soil_transformed, time_series_physeq_plantsonly_filtered_transformed)
  ## Make a dummy ASV table
  ## columns as soil, Anthoceros, Notothylas, Phaeoceros, rows as ASVs
  ## count as number of times this ASV found in soil, Anthoceros, Notothylas, or Phaeoceros samples
    all_ASV_count <- apply(otu_table(time_series_physeq_filtered_transformed) > 0, 1, sum)
    all_ASV_count <- as.data.frame(all_ASV_count)
    colnames(all_ASV_count) <- "a_all_count"
    all_ASV_count <- all_ASV_count[sort(row.names(all_ASV_count)),,drop=F]
  ## count ASV in soil
    soil_ASV_count <- apply(otu_table(subset_samples(time_series_physeq_filtered_transformed, Type=="soil")) > 0, 1, sum)
    soil_ASV_count <- as.data.frame(soil_ASV_count)
    colnames(soil_ASV_count) <- "b_soil_count"
    soil_ASV_count <- soil_ASV_count[sort(row.names(soil_ASV_count)),,drop=F]
  ## count ASV in Anthoceros
    Antho_ASV_count <- apply(otu_table(subset_samples(time_series_physeq_filtered_transformed, Type=="Anthoceros")) > 0, 1, sum)
    Antho_ASV_count <- as.data.frame(Antho_ASV_count)
    colnames(Antho_ASV_count) <- "c_Antho_count"
    Antho_ASV_count <- Antho_ASV_count[sort(row.names(Antho_ASV_count)),,drop=F]
  ## count ASV in Notothylas
    Noto_ASV_count <- apply(otu_table(subset_samples(time_series_physeq_filtered_transformed, Type=="Notothylas")) > 0, 1, sum)
    Noto_ASV_count <- as.data.frame(Noto_ASV_count)
    colnames(Noto_ASV_count) <- "d_Noto_count"
    Noto_ASV_count <- Noto_ASV_count[sort(row.names(Noto_ASV_count)),,drop=F]
  ## count ASV in Notothylas_grossman
    sub_grossman <- subset_samples(time_series_physeq_filtered_transformed, Quadrat=="Grossman1" | Quadrat=="Grossman2")
    sub_grossman_Noto <- subset_samples(sub_grossman, Type=="Notothylas")
    Noto_ASV_count_grossman <- apply(otu_table(sub_grossman_Noto) > 0, 1, sum)
    Noto_ASV_count_grossman <- as.data.frame(Noto_ASV_count_grossman)
    colnames(Noto_ASV_count_grossman) <- "dd_Noto_count"
    Noto_ASV_count_grossman <- Noto_ASV_count_grossman[sort(row.names(Noto_ASV_count_grossman)),,drop=F]
  ## count ASV in Notothylas_potato
    sub_potato <- subset_samples(time_series_physeq_filtered_transformed, Quadrat=="Potato1" | Quadrat=="Potato2"| Quadrat=="Potato3")
    sub_potato_Noto <- subset_samples(sub_potato, Type=="Notothylas")
    Noto_ASV_count_potato <- apply(otu_table(sub_potato_Noto) > 0, 1, sum)
    Noto_ASV_count_potato <- as.data.frame(Noto_ASV_count_potato)
    colnames(Noto_ASV_count_potato) <- "ddd_Noto_count"
    Noto_ASV_count_potato <- Noto_ASV_count_potato[sort(row.names(Noto_ASV_count_potato)),,drop=F]
  ## count ASV in Phaeoceros
    Phae_ASV_count <- apply(otu_table(subset_samples(time_series_physeq_filtered_transformed, Type=="Phaeoceros")) > 0, 1, sum)
    Phae_ASV_count <- as.data.frame(Phae_ASV_count)
    colnames(Phae_ASV_count) <- "e_Phaeo_count"
    Phae_ASV_count <- Phae_ASV_count[sort(row.names(Phae_ASV_count)),,drop=F]
  ## Put all together and make a new phyloseq object
    combined_ASV_count <- cbind(all_ASV_count, soil_ASV_count, Antho_ASV_count, Noto_ASV_count, Noto_ASV_count_grossman, Noto_ASV_count_potato, Phae_ASV_count)
    combined_ASV_count_meta <- sample_data(data.frame(row.names=c("a_all_count","b_soil_count","c_Antho_count","d_Noto_count","dd_Noto_count","ddd_Noto_count","e_Phaeo_count"), "Type"=c("0_all","1_soil","2_Antho","3_Noto","3_Noto_grossman","3_Noto_potato","4_Phaeo"))) #dummy meta
    combined_ASV_count_phyloseq = phyloseq(otu_table(combined_ASV_count, taxa_are_rows = TRUE), combined_ASV_count_meta, tree)
    combined_ASV_count_phyloseq_filtered = filter_taxa(combined_ASV_count_phyloseq, function(x) sum(x) > 0, TRUE) # remove ASV of zero count in ALL categories
  ## Use phyloseq to plot tree
    #plot_tree(combined_ASV_count_phyloseq_filtered, size="abundance", color="Type", label.tips="taxa_names", ladderize="left", base.spacing=0.02, plot.margin=0.1, sizebase=5, justify="left", text.size=1.5) +
    #  scale_colour_manual(values = c("#F0E442","#999999","#E69F00","#56B4E9","#009E73")) 
    #ggsave("ASV_on_target_timeseries_pasta_rooted_annotated.pdf", device = "pdf", width = 15, height = 20)
    #plot_tree(combined_ASV_count_phyloseq_filtered, size="abundance", color="Type", label.tips="taxa_names", ladderize="left", base.spacing=0.02, plot.margin=0.1, sizebase=5, text.size=1.5) + coord_polar(theta="y") +
    #  scale_colour_manual(values = c("#F0E442","#999999","#E69F00","#56B4E9","#009E73")) 
    #ggsave("ASV_on_target_timeseries_pasta_rooted_annotated_cir.pdf", device = "pdf", width = 20, height = 20)
  ## Use ggtree to plot tree
    tree <- phy_tree(combined_ASV_count_phyloseq_filtered)
    ASV_data <- as.data.frame(otu_table(combined_ASV_count_phyloseq_filtered))
    ASV_data <- tibble::rownames_to_column(ASV_data)
    ASV_data <- replace_with_na(ASV_data, replace = list(b_soil_count=0, c_Antho_count=0, d_Noto_count=0,dd_Noto_count=0,ddd_Noto_count=0, e_Phaeo_count=0)) # turn zero count into NA so that ggtree will not plot zeros
    #ASV_data_clade <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/ASV_on_target_timeseries_tax.txt", header=TRUE, sep = "\t")
    ### regular
      d1 <- data.frame(x = seq(1.65, 1.85, length.out = 5)) # control the x position/spread of circles
      p <- ggtree(tree) %<+% ASV_data #%<+% ASV_data_clade # attach ASV_data to the tree/graph object
      p <- p  + 
        geom_tippoint(aes(size=b_soil_count), x = d1$x[1], shape = 16,na.rm=TRUE,color="#999999") + 
        geom_tippoint(aes(size=c_Antho_count), x = d1$x[2], shape = 16,na.rm=TRUE,color="#E69F00") + 
        geom_tippoint(aes(size=d_Noto_count), x = d1$x[3], shape = 16,na.rm=TRUE,color="#56B4E9") + 
        geom_tippoint(aes(size=e_Phaeo_count), x = d1$x[4], shape = 16,na.rm=TRUE,color="#009E73") + 
        geom_tiplab(offset = .25, align=TRUE, linesize=0.1, size=1) 
      #groupClade(p, c(findMRCA(tree, tips=c('ASV1089', 'ASV1634')), findMRCA(tree, tips=c('ASV165', 'ASV1410'))),"clades") + aes(color=clades)
      ggsave("ASV_on_target_timeseries_pasta_rooted_annotated_gg.pdf", device = "pdf", width = 15, height = 20)
      ggsave("ASV_on_target_timeseries_pasta_rooted_annotated_gg.svg", device = "svg", width = 15, height = 20)
    ### circular
      d1 <- tibble(x = seq(1.65, 2.05, length.out = 4)) # control the x position of circles
      p <- ggtree(tree,layout="circular",size=0.2) %<+% ASV_data # attach ASV_data to the tree/graph object
      p <- p+ xlim_tree(4)+
        geom_tippoint(aes(size=b_soil_count), x = d1$x[1], shape = 19,na.rm=TRUE,color="#999999") + 
        geom_tippoint(aes(size=c_Antho_count), x = d1$x[2], shape = 19,na.rm=TRUE,color="#E69F00") + 
        geom_tippoint(aes(size=d_Noto_count), x = d1$x[3], shape = 19,na.rm=TRUE,color="#56B4E9") + 
        #geom_tippoint(aes(size=ddd_Noto_count), x = d1$x[4], shape = 1,na.rm=TRUE,color="#56B4E9") + #potato only
        #geom_tippoint(aes(size=dd_Noto_count), x = d1$x[5], shape = 19,na.rm=TRUE,color="#56B4E9") + #grossman only
        geom_tippoint(aes(size=e_Phaeo_count), x = d1$x[4], shape = 19,na.rm=TRUE,color="#009E73") + 
        geom_tiplab(offset = 0.55, align=TRUE, linesize=0.1, size=2) +
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV1410', 'ASV1763')), label='Nostoc lichens Petigera/Leptogium/Collema/Lobaria/Warshan KVJ10Blasia/KVJ20/KVJ11/KVJ2/Nmoss2/ Nostoc punctiforme PCC73102 / Nostoc ATCC53789', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV521', 'ASV503')), label='Nostoc lichens Petigera', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV178', 'ASV209')), label='Nostoc sphaeroides ASM344365', fontsize=3, align=T, offset=0.85, barsize=0.2) + 
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV50', 'ASV15')), label='Nostoc lichens Peltigera/Massalonga/Leptogium/ Nostoc commune NIES4072', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV253', 'ASV1633')), label='Nostoc lichens Peltigera/Nephroma/Gunnera', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV236', 'ASV1685')), label='Trichormus variabilis ATCC29413 / Nostoc PCC7120 / Nostoc Warshan Nmoss5 / Nostoc Warshan Nmoss6', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV71', 'ASV232')), label='Nostoc cycadae WK1 / Anabaenopsis circularis NIES21', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV275', 'ASV123')), label='Nostoc CENA543 / Nostoc CENA21', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV91', 'ASV480')), label='Nostoc NIES2111 / Nostoc NIES3756', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_strip('ASV1459', 'ASV1459', label='Nodularia CCY9414', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV1581', 'ASV179')), label='Sytonema NIES4073', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV23', 'ASV246')), label='Nostoc lichens Leptogium/Parmotrema/Peltigera', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV26', 'ASV553')), label='Cylindrospermum PCC7417', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV825', 'ASV417')), label='Anabaena/Trichormus/Dolichospermum/Apanizomenon', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV812', 'ASV1082')), label='Nostoc_ellipsosporum_AJ632066', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV155', 'ASV165')), label='Scytonema_UTEX_2349 / Hassallia_KY417069', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_cladelabel(node=findMRCA(tree, tips=c('ASV1634', 'ASV1089')), label='Microcoleus PCC 7113', fontsize=3, align=T, offset=0.85, barsize=0.2)+
        geom_strip('ASV1489', 'ASV1489', label='Kamptonema PCC6407/Ocillatoria PCC7112', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_strip('ASV225', 'ASV225', label='Oscillatoria PCC6304 ', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_strip('ASV1425', 'ASV1425', label='Oculatella ', fontsize=3, align=T, offset=0.85, barsize=0.2) +
        geom_strip('ASV1708', 'ASV634', label='Psedanabaena', fontsize=3, align=T, offset=0.85, barsize=0.2) + 
        geom_hilight(node=findMRCA(tree, tips=c('ASV155', 'ASV165')), fill="gray",alpha=0.3,extend=1.195) + 
        geom_hilight(node=findMRCA(tree, tips=c('ASV23', 'ASV246')), fill="gray",alpha=0.3,extend=1.16) +
        geom_hilight(node=findMRCA(tree, tips=c('ASV40', 'ASV21')), fill="gray",alpha=0.3,extend=1.23) +
        #geom_hilight(node=findMRCA(tree, tips=c('ASV563', 'ASV97')), fill="gray",alpha=0.3,extend=0.972) +
        #geom_hilight(node=findMRCA(tree, tips=c('ASV1320', 'ASV1233')), fill="gray",alpha=0.3,extend=0.88) +
        geom_hilight(node=findMRCA(tree, tips=c('ASV451', 'ASV1233')), fill="gray",alpha=0.3,extend=0.88) +
        geom_hilight(node=findMRCA(tree, tips=c('ASV1354', 'ASV605')), fill="gray",alpha=0.3,extend=1.04)
        #geom_hilight(node=345,extend=0.5) + geom_hilight(node=399,extend=0.5) + geom_hilight(node=406,extend=0.5) + geom_hilight(node=481,extend=0.5) +
        # geom_hilight(node=findMRCA(tree, tips=c('ASV1410', 'ASV1763')),extend=0.5) #+ 
        #+ geom_label2(aes(subset=!isTip, label=node), size=1, color="darkred", alpha=0.5)
        #+ scale_size_continuous()
      groupClade(p, c(findMRCA(tree, tips=c('ASV1089', 'ASV1634')), findMRCA(tree, tips=c('ASV165', 'ASV1410'))),"clades") + 
        aes(color=clades) + scale_color_manual(values=c("#1B9E77", "#D95F02","#252525")) #+ scale_color_brewer(palette="Dark2")
      ggsave("ASV_on_target_timeseries_pasta_rooted_annotated_gg_cir_test.pdf", device = "pdf", width = 15, height = 20)
      ggsave("ASV_on_target_timeseries_pasta_rooted_annotated_gg_cir_test.svg", device = "svg", width = 15, height = 20)
      ggsave("ASV_on_target_timeseries_pasta_rooted_annotated_gg_cir.svg", device = "svg", width = 15, height = 20)
      
# Distribution of reads and ASV -----------------------------------------------
  ## Read number distribution
    read_count_file <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/sample_fastq_deprimers_lenfiltered_JNP1-3/PG_read_number_distribution")
    read_count_file$V1 <- read_count_file$V1 /4
    mu <- data.frame(Type="all", grp.mean=mean(read_count_file$V1))
    p_readdistr <- ggplot(read_count_file, aes(x=V1)) + geom_histogram(color="black", fill="lightblue", binwidth=500, alpha=0.5) +
      xlab("Number of reads") + ylab("Sample count") + ggtitle("Distribution of read number") + 
      geom_vline(data=mu, aes(xintercept=grp.mean), linetype="dashed") 
  ## ASV no. distribution
    ASV_distribution_func <- function(observationThreshold){
      ASV_count_plant <- apply(otu_table(time_series_physeq_plant_transformed) > observationThreshold, 2, sum)
      ASV_count_plant_df <- as.data.frame(ASV_count_plant)
      ASV_count_plant_df$Type = rep("plant", length(ASV_count_plant))
      colnames(ASV_count_plant_df) <- c("ASV_count","Type")
      ASV_count_soil <- apply(otu_table(time_series_physeq_soil_transformed) > observationThreshold, 2, sum)
      ASV_count_soil_df <- as.data.frame(ASV_count_soil)
      ASV_count_soil_df$Type = rep("soil", length(ASV_count_soil))
      colnames(ASV_count_soil_df) <- c("ASV_count","Type")
      ASV_count <- rbind(ASV_count_plant_df, ASV_count_soil_df)
      return(ASV_count)
      }
    ### When observationThreshold = 0.03
      ASV_count <- ASV_distribution_func(0.03)
      mu <- ddply(ASV_count, "Type", summarise, grp.mean=mean(ASV_count)) # calculate mean
      p_ASVdistr_3 <- ggplot(as.data.frame(subset(ASV_count,Type=='plant')), aes(x=ASV_count, fill=Type, color=Type)) + 
        geom_histogram(binwidth=1, alpha=0.5) + scale_x_continuous(breaks=seq(1,10,by=2)) +
        xlab("Number of ASV") + ylab("Sample count") + ggtitle("Distribution of ASV") +
        scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
        geom_vline(data=subset(mu,Type=='plant'), aes(xintercept=grp.mean, color=Type), linetype="dashed") 
    ### When observationThreshold = 0.05
      ASV_count <- ASV_distribution_func(0.05)
      mu <- ddply(ASV_count, "Type", summarise, grp.mean=mean(ASV_count)) # calculate mean
      p_ASVdistr_5 <- ggplot(as.data.frame(subset(ASV_count,Type=='plant')), aes(x=ASV_count, fill=Type, color=Type)) + 
        geom_histogram(binwidth=1, alpha=0.5) +
        xlab("Number of ASV") + ylab("Sample count") + ggtitle("Distribution of ASV") +
        scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
        geom_vline(data=subset(mu,Type=='plant'), aes(xintercept=grp.mean, color=Type), linetype="dashed") 
    ### When observationThreshold = 0
      ASV_count <- ASV_distribution_func(0)
      mu <- ddply(ASV_count, "Type", summarise, grp.mean=mean(ASV_count)) # calculate mean
      p_ASVdistr_0 <- ggplot(as.data.frame(ASV_count), aes(x=ASV_count, fill=Type, color=Type)) + 
        geom_histogram(binwidth=2, alpha=0.5, position="identity") +
        xlab("Number of ASV") + ylab("Sample count") + ggtitle("Distribution of ASV") +
        scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
        geom_vline(data=mu, aes(xintercept=grp.mean, color=Type), linetype="dashed") 
    ### Plot
      p_readdistr + p_ASVdistr_0 + p_ASVdistr_3 + plot_layout(guides = 'collect')
      p_distr <- p_readdistr + p_ASVdistr_0 + p_ASVdistr_3 + plot_layout(guides = 'collect')
      ggsave("Read_ASV_distribution.pdf", device = "pdf", width = 18, height = 6)
    
    ### Read count vs ASV count
      read_count_file <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/sample_fastq_deprimers_lenfiltered_JNP1-3/PG_read_number_distribution_wFilename")
      read_count_file <- tibble::column_to_rownames(read_count_file, var="V1")
      read_count_file$V2 <- read_count_file$V2 /4
      read_count_ASV_count <- merge(read_count_file, ASV_count, by="row.names")
      p_read_count_ASV_count <- ggplot(data=read_count_ASV_count, mapping=aes(V2,ASV_count)) + 
        geom_point(aes(color=Type)) + geom_smooth(method = lm) + 
        xlab("Read count") + ylab("Number of ASVs") 
      ggsave("read_count_ASV_count.svg", device = "svg", width = 8, height = 4)
      linearMod <- lm(ASV_count ~ V2, data=read_count_ASV_count)
      summary(linearMod)
      
# Rarefraction curve ------------------------------------------------------------
  ggrare(time_series_physeq, color="Type", step=100, parallel=T) + 
      scale_fill_manual(values = cbp1) + scale_color_manual(values = cbp1) + facet_wrap(~Type)
  ggsave("rarefraction.pdf", device = "pdf", width = 8, height = 6)

# Species accumulation curve ------------------------------------------------------------
  data_obj <- subset_samples(time_series_physeq, Quadrat=="Grossman1")
  data_obj <- time_series_physeq
  sp1 <- specaccum(t(otu_table(data_obj)))
  sp2 <- specaccum(t(otu_table(data_obj)), "random")
  p1 <- plot(sp2, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
  boxplot(sp2, col="yellow", add=TRUE, pch="+")
  par(mfrow=c(5,2))
  par(mar=c(2,2,2,2))
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_plant_filtered, Quadrat=="Grossman1") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Grossman1 plants")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_soil, Quadrat=="Grossman1") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Grossman1 soil")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_plant_filtered, Quadrat=="Grossman2") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Grossman2 plants")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_soil, Quadrat=="Grossman2") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Grossman2 soil")
  
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_plant_filtered, Quadrat=="Potato1") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Potato1 plants")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_soil, Quadrat=="Potato1") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Potato2 soil")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_plant_filtered, Quadrat=="Potato2") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Potato2 plants")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_soil, Quadrat=="Potato2") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Potato2 soil")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_plant_filtered, Quadrat=="Potato3") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Potato3 plants")
  plot(specaccum(t(otu_table( subset_samples(time_series_physeq_soil, Quadrat=="Potato3") )), "random"), ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Potato3 soil")
  
# GrossmanPond ------------------------------------------------------------
  ## Subset GrossmanPond samples
    grossman_physeq <- subset_samples(time_series_physeq, Quadrat=="Grossman1"|Quadrat=="Grossman2")
  ## Transform
    ## Two ways to transform: 1. simply turn into relative abundance, 
    ##                        2. zero out ASV with <5% abundance first before turn into relative abundance
    ## 1. Transform counts into relative abundance
    grossman_physeq_transformed <- transform_sample_counts(grossman_physeq, function(x) x/sum(x))
    ## 2. a) Transform plant counts into proportion AND remove low abundant ASV (5%) while keeping the sum=1
    ##    b) transform soil sounts into proportion 
    ##    c) merge the two
    grossman_physeq_plantsonly <- subset_samples(grossman_physeq, Type=="Notothylas")
    grossman_physeq_plantsonly_filtered <- transform_sample_counts(grossman_physeq_plantsonly, fun = filter3perc_count)  
    grossman_physeq_plantsonly_transformed <- transform_sample_counts(grossman_physeq_plantsonly_filtered, function(x) x / sum(x) )
    grossman_physeq_soil <- subset_samples(grossman_physeq, Type=="soil")
    grossman_physeq_soil_transformed <- transform_sample_counts(grossman_physeq_soil, function(x) x/sum(x))  
    grossman_physeq_transformed <- merge_phyloseq(grossman_physeq_plantsonly_transformed, grossman_physeq_soil_transformed)
    grossman_physeq_filtered <- merge_phyloseq(grossman_physeq_plantsonly_filtered, grossman_physeq_soil)
  ## ALL ordination plot
    grossman_physeq_transformed.ord <- ordinate(grossman_physeq_transformed, "PCoA", "bray")
    grossman_pcoa_all <- plot_ordination(grossman_physeq_transformed, grossman_physeq_transformed.ord, title='GrossmanPond all', type="samples", color="Quadrat", shape="Type") +
      geom_point(size=4) + scale_color_manual(values=c("lightcoral","#7570B3")) + 
      coord_fixed()
    grossman_pcoa_DF_all <- plot_ordination(grossman_physeq_transformed, grossman_physeq_transformed.ord, type="samples", color="Quadrat", justDF=TRUE) 
    grossman_pcoa_soil <- ggplot(data = subset(grossman_pcoa_DF_all, Type=='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
      geom_point(size=4, shape=17, aes(color=Quadrat)) + geom_point(shape = 2,size = 4,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3")) + 
      coord_fixed() + xlim(-0.52, 0.52)+ ylim(-0.55, 0.35) + ggtitle("PotatoHill Soils")
    grossman_pcoa_plant <- ggplot(data = subset(grossman_pcoa_DF_all, Type!='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
      geom_point(size=4, aes(color=Quadrat)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("lightcoral","#7570B3")) + 
      coord_fixed() + xlim(-0.52, 0.52)+ ylim(-0.55, 0.35) + ggtitle("PotatoHill Soils")
    ig <- make_network(grossman_physeq_transformed, distance="bray", type="samples", max.dist=0.6)
    grossman_net <- plot_network(ig, grossman_physeq_transformed, color = "Type", shape = "Quadrat", label=NULL)
  ## PERMANOVA test
    grossman_bray_dist <- phyloseq::distance(grossman_physeq_transformed, method="bray")
    adonis(grossman_bray_dist ~ sample_data(grossman_physeq_transformed)$Type, permutations=10000)
    adonis(grossman_bray_dist ~ sample_data(grossman_physeq_transformed)$Quadrat, permutations=10000)
    adonis(grossman_bray_dist ~ sample_data(grossman_physeq_transformed)$Quadrat+sample_data(grossman_physeq_transformed)$Type, permutations = 10000)
  ## UpSet 
    upset(phyloseq2upsetr(grossman_physeq_filtered, "Type"),order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(grossman_physeq_filtered, "Quadrat"),order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(subset_samples(grossman_physeq_filtered, Quadrat=="Grossman1"), "Type"),order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(subset_samples(grossman_physeq_filtered, Quadrat=="Grossman2"), "Type"),order.by = "freq", empty.intersections = "on")
    
  ## ^Quadrat 1 ====
    ### subset
      grossman1_physeq_transformed <- subset_samples(grossman_physeq_transformed, Quadrat=="Grossman1")
    ### ordination plot
      grossman1_physeq_transformed.ord <- ordinate(grossman1_physeq_transformed, "PCoA", "wunifrac")
      grossman1_pcoa_all <- plot_ordination(grossman1_physeq_transformed, grossman1_physeq_transformed.ord, type="samples", color="Time") +
        geom_point(size=4) + scale_color_brewer(palette="YlGnBu") + labs(title="Grossman Q1 by Time") + facet_wrap(~Type, nrow=1) #+ coord_fixed()  
    ### PERMANOVA
      grossman1_bray_dist <- phyloseq::distance(grossman1_physeq_transformed, method="wunifrac")
      permanova <- adonis(grossman1_bray_dist ~ sample_data(grossman1_physeq_transformed)$Time, permutations = 10000)
    ### network
      ig <- make_network(grossman1_physeq_transformed, distance="bray", type="samples", max.dist=0.4)
      grossman1_net <- plot_network(ig, grossman1_physeq_transformed, color = "Type", label=NULL)
    ### Notothylas by time
      grossman1_physeq_transformed_Noto <- subset_samples(grossman1_physeq_transformed, Type=="Notothylas")
      grossman1_physeq_transformed_Noto.ord <- ordinate(grossman1_physeq_transformed_Noto, "PCoA", "bray")
      grossman1_bray_dist <- phyloseq::distance(grossman1_physeq_transformed_Noto, method="bray")
      permanova <- adonis(grossman1_bray_dist ~ sample_data(grossman1_physeq_transformed_Noto)$Time, permutations = 10000)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      grossman1_pcoa_Noto <- plot_ordination(grossman1_physeq_transformed_Noto, grossman1_physeq_transformed_Noto.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Grossman Q1 Notothylas",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(grossman1_physeq_transformed_Noto, distance="bray", max.dist=0.4)
      grossman1_Noto_net <- plot_network(ig, grossman1_physeq_transformed_Noto, color = "Time", label=NULL)
    ### soil by time
      grossman1_physeq_transformed_soil <- subset_samples(grossman1_physeq_transformed, Type=="soil")
      grossman1_physeq_transformed_soil.ord <- ordinate(grossman1_physeq_transformed_soil, "PCoA", "bray")
      grossman1_bray_dist <- phyloseq::distance(grossman1_physeq_transformed_soil, method="bray")
      permanova <- adonis(grossman1_bray_dist ~ sample_data(grossman1_physeq_transformed_soil)$Time, permutations = 10000)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      grossman1_pcoa_soil <- plot_ordination(grossman1_physeq_transformed_soil, grossman1_physeq_transformed_soil.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Grossman Q1 soil",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(grossman1_physeq_transformed_soil, distance="bray", max.dist=0.4)
      grossman1_soil_net <- plot_network(ig, grossman1_physeq_transformed_soil, color = "Time", label=NULL)
    
  ## ^Quadrat 2 ====
    ### transform counts into proportion
      grossman2_physeq_transformed <- subset_samples(grossman_physeq_transformed, Quadrat=="Grossman2")
    ### ordination plot
      grossman2_physeq_transformed.ord <- ordinate(grossman2_physeq_transformed, "PCoA", "wunifrac")
      grossman2_pcoa_all <- plot_ordination(grossman2_physeq_transformed, grossman2_physeq_transformed.ord, type="samples", color="Time") +
        geom_point(size=4) + scale_color_brewer(palette="YlGnBu") + labs(title="Grossman Q2 by Time") + facet_wrap(~Type, nrow=1) #+ coord_fixed()  
    ### PERMANOVA
      grossman2_bray_dist <- phyloseq::distance(grossman2_physeq_transformed, method="wunifrac")
      permanova <- adonis(grossman2_bray_dist ~ sample_data(grossman2_physeq_transformed)$Time, permutations = 10000)
    ### network
      ig <- make_network(grossman2_physeq_transformed, distance="bray", max.dist=0.8)
      grossman2_net <- plot_network(ig, grossman2_physeq_transformed, color = "Type", label=NULL)
    ### Notothylas by time
      grossman2_physeq_transformed_Noto <- subset_samples(grossman2_physeq_transformed, Type=="Notothylas")
      grossman2_physeq_transformed_Noto.ord <- ordinate(grossman2_physeq_transformed_Noto, "PCoA", "bray")
      grossman2_bray_dist <- phyloseq::distance(grossman2_physeq_transformed_Noto, method="bray")
      permanova <- adonis(grossman2_bray_dist ~ sample_data(grossman2_physeq_transformed_Noto)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      grossman2_pcoa_Noto <- plot_ordination(grossman2_physeq_transformed_Noto, grossman2_physeq_transformed_Noto.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Grossman Q2 Notothylas",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(grossman2_physeq_transformed_Noto, distance="bray", max.dist=0.4)
      grossman2_Noto_net <- plot_network(ig, grossman2_physeq_transformed_Noto, color = "Time", label=NULL)
    ### soil by time
      grossman2_physeq_transformed_soil <- subset_samples(grossman2_physeq_transformed, Type=="soil")
      grossman2_physeq_transformed_soil.ord <- ordinate(grossman2_physeq_transformed_soil, "PCoA", "bray")
      grossman2_bray_dist <- phyloseq::distance(grossman2_physeq_transformed_soil, method="bray")
      permanova <- adonis(grossman2_bray_dist ~ sample_data(grossman2_physeq_transformed_soil)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      grossman2_pcoa_soil <- plot_ordination(grossman2_physeq_transformed_soil, grossman2_physeq_transformed_soil.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Grossman Q2 soil",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(grossman2_physeq_transformed_soil, distance="bray", max.dist=0.4)
      grossman2_soil_net <- plot_network(ig, grossman2_physeq_transformed_soil, color = "Time", label=NULL)
      
  ## Make plots  
    #(grossman_pcoa_all | (grossman1_pcoa_all / grossman2_pcoa_all)) + plot_layout(guides = 'collect')
    grossman_pcoa_all
    ggsave("grossman_pcoa_all.pdf", device = "pdf", width = 18, height = 12)
    (grossman1_pcoa_soil | grossman1_pcoa_Noto | grossman2_pcoa_soil | grossman2_pcoa_Noto ) + plot_layout(guides = 'collect')
    ggsave("grossman_pcoa_byQbyT.pdf", device = "pdf", width = 18, height = 5)

  ## Core ASV
    grossman1_T1 <- phyloseq2core(subset_samples(grossman1_physeq_transformed_Noto, Time=='T1'), 'Notothylas', 'Grossman1', 'T1')
    grossman1_T2 <- phyloseq2core(subset_samples(grossman1_physeq_transformed_Noto, Time=='T2'), 'Notothylas', 'Grossman1', 'T2')
    grossman1_T3 <- phyloseq2core(subset_samples(grossman1_physeq_transformed_Noto, Time=='T3'), 'Notothylas', 'Grossman1', 'T3')
    grossman1_T4 <- phyloseq2core(subset_samples(grossman1_physeq_transformed_Noto, Time=='T4'), 'Notothylas', 'Grossman1', 'T4')
    grossman2_T1 <- phyloseq2core(subset_samples(grossman2_physeq_transformed_Noto, Time=='T1'), 'Notothylas', 'Grossman2', 'T1')
    grossman2_T2 <- phyloseq2core(subset_samples(grossman2_physeq_transformed_Noto, Time=='T2'), 'Notothylas', 'Grossman2', 'T2')
    grossman2_T3 <- phyloseq2core(subset_samples(grossman2_physeq_transformed_Noto, Time=='T3'), 'Notothylas', 'Grossman2', 'T3')
    
    grossman1_T1_soil <- phyloseq2core(subset_samples(grossman1_physeq_transformed_soil, Time=='T1'), 'soil', 'Grossman1', 'T1')
    grossman1_T2_soil <- phyloseq2core(subset_samples(grossman1_physeq_transformed_soil, Time=='T2'), 'soil', 'Grossman1', 'T2')
    grossman1_T3_soil <- phyloseq2core(subset_samples(grossman1_physeq_transformed_soil, Time=='T3'), 'soil', 'Grossman1', 'T3')
    grossman1_T4_soil <- phyloseq2core(subset_samples(grossman1_physeq_transformed_soil, Time=='T4'), 'soil', 'Grossman1', 'T4')
    grossman2_T1_soil <- phyloseq2core(subset_samples(grossman2_physeq_transformed_soil, Time=='T1'), 'soil', 'Grossman2', 'T1')
    grossman2_T2_soil <- phyloseq2core(subset_samples(grossman2_physeq_transformed_soil, Time=='T2'), 'soil', 'Grossman2', 'T2')
    grossman2_T3_soil <- phyloseq2core(subset_samples(grossman2_physeq_transformed_soil, Time=='T3'), 'soil', 'Grossman2', 'T3')
    
    grossman_Tall <- bind_rows(grossman1_T1, grossman1_T2, grossman1_T3, grossman1_T4, grossman2_T1, grossman2_T3)
    grossman_core_stickgraph <- ggplot(data=grossman_Tall) + geom_point(aes(y = reorder(ASV, desc(ASV)), x=abundance, shape=Time, color=Type)) +
      xlab('Relative abundance') + ylab('') + scale_color_manual(values=c("steelblue","gold2","darkseagreen"))
    #grossman_heatmap <- ggplot(data=grossman_Tall, aes(x=ASV, y=replicate, fill=value)) + geom_tile()
    tbl2hmap(tbl.long = grossman_Tall, rowVar = 'ASV', colVar = 'replicate', valueVar = 'value')$mat
    grossman_heatmap <- as.ggplot( pheatmap(tbl2hmap(tbl.long = grossman_Tall, rowVar = 'ASV', colVar = 'replicate', valueVar = 'value')$mat ,cluster_rows=FALSE, cluster_cols=FALSE) )
    
    
# Potato Hill ------------------------------------------------------------
  ## Subset PotatoHill samples
    potato_physeq <- subset_samples(time_series_physeq, Quadrat=="Potato1"|Quadrat=="Potato2"|Quadrat=="Potato3")
  ## Transform
    ## Two ways to transform: 1. simply turn into relative abundance, 
    ##                        2. zero out ASV with <5% abundance first before turn into relative abundance
    ## 1. Transform counts into relative abundance
    #potato_physeq_transformed <- transform_sample_counts(potato_physeq, function(x) x/sum(x))
    ## 2. a) Transform plant counts into proportion AND remove low abundant ASV (5%) while keeping the sum=1
    ##    b) transform soil sounts into proportion 
    ##    c) merge the two
    potato_physeq_plantsonly <- subset_samples(potato_physeq, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
    potato_physeq_plantsonly_filtered <- transform_sample_counts(potato_physeq_plantsonly, fun = filter3perc_count)  
    potato_physeq_plantsonly_transformed <- transform_sample_counts(potato_physeq_plantsonly_filtered, function(x) x / sum(x) )
    potato_physeq_soil <- subset_samples(potato_physeq, Type=="soil")
    potato_physeq_soil_transformed <- transform_sample_counts(potato_physeq_soil, function(x) x/sum(x))  
    potato_physeq_transformed <- merge_phyloseq(potato_physeq_plantsonly_transformed, potato_physeq_soil_transformed)
    potato_physeq_filtered <- merge_phyloseq(potato_physeq_plantsonly_filtered, potato_physeq_soil)
  ## ALL ordination plot
    potato_physeq_transformed.ord <- ordinate(potato_physeq_transformed, "PCoA", "wunifrac")
    potato_pcoa_all <- plot_ordination(potato_physeq_transformed, potato_physeq_transformed.ord, type="samples", color="Quadrat", shape="Type") #+ geom_polygon(aes(fill=Quadrat))
    potato_pcoa_DF_all <- plot_ordination(potato_physeq_transformed, potato_physeq_transformed.ord, type="samples", color="Quadrat", justDF=TRUE) 
  ## ALL ordination but soil and plant separately
    potato_pcoa_all_soil <- ggplot(data = subset(potato_pcoa_DF_all, Type=='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
      geom_point(size=4, shape=17, aes(color=Quadrat)) + geom_point(shape = 2,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) + 
      coord_fixed() + xlim(-0.52, 0.52)+ ylim(-0.55, 0.35) + ggtitle("PotatoHill Soils")
    potato_pcoa_all_plantsonly <- ggplot(data = subset(potato_pcoa_DF_all, Type!='soil'), mapping = aes(x = Axis.1, y = Axis.2)) +
      geom_point(size=4, aes(color=Quadrat)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) + 
      coord_fixed()+ xlim(-0.52, 0.52) + ylim(-0.55, 0.35) + ggtitle("PotatoHill Plants")
  ## Plant samples, unweighted unifrac
    potato_physeq_plantsonly_transformed.ord <- ordinate(potato_physeq_plantsonly_transformed, "PCoA", "unifrac")
    potato_pcoa_all_plantsonly <- plot_ordination(potato_physeq_plantsonly_transformed, potato_physeq_plantsonly_transformed.ord, color="Quadrat") +
      geom_point(size=4, aes(color=Quadrat)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) +
      ggtitle("Plants unweighted unifrac")
    potato_pcoa_all_plantsonly <- plot_ordination(potato_physeq_plantsonly_transformed, potato_physeq_plantsonly_transformed.ord, color="Type") +
      geom_point(size=4, aes(color=Type)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) +
      ggtitle("Plants unweighted unifrac")
  ## Soil samples, unweighted unifrac
    potato_physeq_soil_transformed.ord <- ordinate(potato_physeq_soil_transformed, "PCoA", "unifrac")
    potato_pcoa_all_soil <- plot_ordination(potato_physeq_soil_transformed, potato_physeq_soil_transformed.ord, color="Quadrat") +
      geom_point(size=4, aes(color=Quadrat)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) +
      ggtitle("Soils unweighted unifrac")
  ## Plant samples, weighted unifrac
    potato_physeq_plantsonly_transformed_w.ord <- ordinate(potato_physeq_plantsonly_transformed, "PCoA", "wunifrac")
    potato_pcoa_all_plantsonly_w <- plot_ordination(potato_physeq_plantsonly_transformed, potato_physeq_plantsonly_transformed_w.ord, color="Quadrat") +
      geom_point(size=4, aes(color=Quadrat)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) +
      ggtitle("Plants weighted unifrac")
    potato_pcoa_all_plantsonly_w <- plot_ordination(potato_physeq_plantsonly_transformed, potato_physeq_plantsonly_transformed_w.ord, color="Type") +
      geom_point(size=4, aes(color=Type)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) +
      ggtitle("Plants weighted unifrac")
  ## Soil samples, weighted unifrac
    potato_physeq_soil_transformed_w.ord <- ordinate(potato_physeq_soil_transformed, "PCoA", "wunifrac", )
    potato_pcoa_all_soil_w <- plot_ordination(potato_physeq_soil_transformed, potato_physeq_soil_transformed_w.ord, color="Quadrat") +
      geom_point(size=4, aes(color=Quadrat)) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_manual(values=c("steelblue","gold2","darkseagreen")) +
      ggtitle("Soils weighted unifrac")
    
  (potato_pcoa_all_soil_w | potato_pcoa_all_soil) / (potato_pcoa_all_plantsonly_w | potato_pcoa_all_plantsonly) + plot_layout(guides = 'collect')
  ggsave("potato_pcoa_plants_soils_weighted_unweighted_unifrac.svg", device = "svg", width = 12, height = 10)
    
  ## PERMANOVA test
    potato_bray_dist <- phyloseq::distance(potato_physeq_transformed, method="bray")
    adonis(potato_bray_dist ~ sample_data(potato_physeq_transformed)$Type, permutations = 10000)
    adonis(potato_bray_dist ~ sample_data(potato_physeq_transformed)$Quadrat, permutations = 10000)
    adonis(potato_bray_dist ~ sample_data(potato_physeq_transformed)$Quadrat+sample_data(potato_physeq_transformed)$Type, permutations = 10000)
    simplified_type <- sample_data(potato_physeq_transformed)$Type
    levels(simplified_type)<-c(levels(simplified_type),"plant")
    simplified_type[simplified_type=="Notothylas"]<-"plant"
    simplified_type[simplified_type=="Anthoceros"]<-"plant"
    simplified_type[simplified_type=="Phaeoceros"]<-"plant"
    adonis(potato_bray_dist ~ sample_data(potato_physeq_transformed)$Quadrat+simplified_type, permutations = 10000)
  ## UpSet 
    upset(phyloseq2upsetr(potato_physeq_filtered, "Type"), sets = c("soil", "Phaeoceros", "Notothylas", "Anthoceros"), keep.order = T, order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(potato_physeq_filtered, "QuadratxType"), sets = c("soil", "Phaeoceros", "Notothylas", "Anthoceros"), keep.order = T, nsets=15, order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(subset_samples(potato_physeq_filtered, Quadrat=="Potato1"), "Type"), sets = c("soil", "Phaeoceros", "Notothylas", "Anthoceros"), keep.order = T, order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(subset_samples(potato_physeq_filtered, Quadrat=="Potato2"), "Type"), sets = c("soil", "Phaeoceros", "Notothylas", "Anthoceros"), keep.order = T, order.by = "freq", empty.intersections = "on")
    upset(phyloseq2upsetr(subset_samples(potato_physeq_filtered, Quadrat=="Potato3"), "Type"), sets = c("soil", "Phaeoceros", "Notothylas", "Anthoceros"), keep.order = T, order.by = "freq", empty.intersections = "on")
    
    intersect_df <- phyloseq2upsetr(potato_physeq_filtered, "Type")
    intersect_df <- tibble::rownames_to_column(intersect_df, "ASV")
    uniqueASV <- as.vector(filter(intersect_df, (Anthoceros==0&Notothylas==0&soil==0)|(Anthoceros==0&Phaeoceros==0&soil==0)|(Notothylas==0&Phaeoceros==0&soil==0))$ASV)
    uniqueASV_melt <- reshape2::melt(otu_table( prune_taxa(uniqueASV, potato_physeq_transformed) ),value.name="abundance",varnames=c("ASV","Sample"))
    uniqueASV_counts <- filter(uniqueASV_melt, abundance > 0) %>% mutate(UniqShared="unqiue")
    sharedASV <- as.vector(filter(intersect_df, !((Anthoceros==0&Notothylas==0&soil==0)|(Anthoceros==0&Phaeoceros==0&soil==0)|(Notothylas==0&Phaeoceros==0&soil==0)))$ASV)
    sharedASV_melt <- reshape2::melt(otu_table( prune_taxa(sharedASV, potato_physeq_transformed)   ),value.name="abundance",varnames=c("ASV","Sample"))
    sharedASV_count <- filter(sharedASV_melt, abundance > 0) %>% mutate(UniqShared="shared")
    ggplot(data=rbind(uniqueASV_counts,sharedASV_count)) + geom_violin(aes(y=abundance,x=UniqShared)) + geom_boxplot(aes(y=abundance,x=UniqShared),width=0.1)
  ## Network
    ig <- make_network(potato_physeq_transformed, distance="bray", max.dist=0.4)
    potato_net_all <- plot_network(ig, potato_physeq_transformed, color = "Type", label=NULL, point_size=4) +
      geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    
  ## ^Quadrat 3 ====
    ### transform counts into proportion
      potato3_physeq_transformed <- subset_samples(potato_physeq_transformed, Quadrat=="Potato3")
    ### PERMANOVA: soil + plants by type
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed)$Type, permutations = 10000)
      permanova_des_all <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      pairwise.perm.manova(potato3_bray_dist,sample_data(potato3_physeq_transformed)$Type,nperm=999)
      anova(betadisper(potato3_bray_dist, sample_data(potato3_physeq_transformed)$Type))  
    ### PERMANOVA: plants only by type
      potato3_physeq_transformed_plantsonly <- subset_samples(potato3_physeq_transformed, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_plantsonly, method="bray")
      #adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_plantsonly)$Time)
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_plantsonly)$Type, permutations = 10000)
      permanova_des_plant <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      pairwise.perm.manova(potato3_bray_dist,sample_data(potato3_physeq_transformed_plantsonly)$Type,nperm=999)
      anova(betadisper(potato3_bray_dist, sample_data(potato3_physeq_transformed_plantsonly)$Type))  
    ### PERMANOVA: plants only by type Unifrac
      potato3_unifrac_dist <- phyloseq::distance(potato3_physeq_transformed_plantsonly, method="wunifrac")
      permanova <- adonis(potato3_unifrac_dist ~ sample_data(potato3_physeq_transformed_plantsonly)$Type, permutations = 10000)
      permanova_des_plant_unifrac <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      beta <- betadisper(potato3_unifrac_dist, sample_data(potato3_physeq_transformed_plantsonly)$Type)
      permutest(beta)
    ### PERMANOVA: soil + plants by time
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed, method="wunifrac")
      adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed)$Type + sample_data(potato3_physeq_transformed)$Time,permutations = 10000)
      adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed)$Time,permutations = 10000)
    ### ordination plot: soil + plants by type
      potato3_physeq_transformed.ord <- ordinate(potato3_physeq_transformed, "PCoA", "bray")
      potato3_pcoa_all <- plot_ordination(potato3_physeq_transformed, potato3_physeq_transformed.ord, type="samples", color="Time") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 plant+soil",caption=permanova_des_all) #+ coord_fixed() 
    ### ordination plot: plants only by type
      potato3_physeq_transformed_plantsonly.ord <- ordinate(potato3_physeq_transformed_plantsonly, "PCoA")
      potato3_pcoa <- plot_ordination(potato3_physeq_transformed_plantsonly, potato3_physeq_transformed_plantsonly.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 plant",caption=permanova_des_plant) #+ coord_fixed()  
    ### ordination plot: plants only by type Unifrac
      potato3_physeq_transformed_plantsonly_unifrac.ord <- ordinate(potato3_physeq_transformed_plantsonly, "PCoA", "wunifrac")
      potato3_pcoa_unifrac <- plot_ordination(potato3_physeq_transformed_plantsonly, potato3_physeq_transformed_plantsonly_unifrac.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 plant",caption=permanova_des_plant_unifrac) #+ coord_fixed()  
    ### ordination plot: soil + plants by time Unifrac
      potato3_physeq_transformed_unifrac.ord <- ordinate(potato3_physeq_transformed, "PCoA", "wunifrac")
      potato3_pcoa_byT <- plot_ordination(potato3_physeq_transformed, potato3_physeq_transformed_unifrac.ord, type="samples", color="Time") +
        geom_point(size=4) + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q3 by Time") + facet_wrap(~Type, nrow=1) #+ coord_fixed()  
    ### plot top influencer taxa among plant hosts
      permanova <- adonis(t(otu_table(potato3_physeq_transformed_plantsonly)) ~ sample_data(potato3_physeq_transformed_plantsonly)$Type)
      coef <- coefficients(permanova)["sample_data(potato3_physeq_transformed_plantsonly)$Type1",]
      top.coef <- coef[rev(order(abs(coef)))[1:10]]
      df <- as.data.frame(top.coef)
      colnames(df) <- "coef"
      top_q3_plant_influencer <- ggplot(df, aes(x=coef, y=reorder(rownames(df),coef))) + geom_bar(stat="identity") + 
        ggtitle("Top influencer ASV in PotatoHill Q3 plants") + ylab("ASV")
    ### network: soil + plants
      plot_net(potato3_physeq_transformed, distance = "bray", maxdist = 0.4, color = "Type", laymeth="fruchterman.reingold")
      ig <- make_network(potato3_physeq_transformed, distance="bray", max.dist=0.4)
      potato3_net_all <- plot_network(ig, potato3_physeq_transformed, color = "Type", label=NULL, point_size=4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### network: plants only
      plot_net(potato3_physeq_transformed_plantsonly, distance = "bray", maxdist = 0.4, color = "Type", laymeth="fruchterman.reingold")
      ig <- make_network(potato3_physeq_transformed_plantsonly, distance="bray", max.dist=0.4)
      potato3_net <- plot_network(ig, potato3_physeq_transformed_plantsonly, color = "Type", label=NULL, point_size=4) + 
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### network: plants only Unifrac
      ig <- make_network(potato3_physeq_transformed_plantsonly, distance="wunifrac", max.dist=0.1)
      potato3_net_unifrac <- plot_network(ig, potato3_physeq_transformed_plantsonly, color = "Type", label=NULL, point_size=4) + 
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### T1
      potato3_physeq_transformed_T1 <- subset_samples(potato3_physeq_transformed_plantsonly, Time=="T1")
      potato3_physeq_transformed_T1.ord <- ordinate(potato3_physeq_transformed_T1, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_T1, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_T1)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_T1 <- plot_ordination(potato3_physeq_transformed_T1, potato3_physeq_transformed_T1.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 T1",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato3_physeq_transformed_T1, distance="bray", max.dist=0.4)
      potato3_T1_net <- plot_network(ig, potato3_physeq_transformed_T1, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T2
      potato3_physeq_transformed_T2 <- subset_samples(potato3_physeq_transformed_plantsonly, Time=="T2")
      potato3_physeq_transformed_T2.ord <- ordinate(potato3_physeq_transformed_T2, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_T2, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_T2)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_T2 <- plot_ordination(potato3_physeq_transformed_T2, potato3_physeq_transformed_T2.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 T2",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato3_physeq_transformed_T2, distance="bray", max.dist=0.4)
      potato3_T2_net <- plot_network(ig, potato3_physeq_transformed_T2, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T3
      potato3_physeq_transformed_T3 <- subset_samples(potato3_physeq_transformed_plantsonly, Time=="T3")
      potato3_physeq_transformed_T3.ord <- ordinate(potato3_physeq_transformed_T3, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_T3, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_T3)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_T3 <- plot_ordination(potato3_physeq_transformed_T3, potato3_physeq_transformed_T3.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 T3",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato3_physeq_transformed_T3, distance="bray", max.dist=0.4)
      potato3_T3_net <- plot_network(ig, potato3_physeq_transformed_T3, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T4
      potato3_physeq_transformed_T4 <- subset_samples(potato3_physeq_transformed_plantsonly, Time=="T4")
      potato3_physeq_transformed_T4.ord <- ordinate(potato3_physeq_transformed_T4, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_T4, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_T4)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_T4 <- plot_ordination(potato3_physeq_transformed_T4, potato3_physeq_transformed_T4.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q3 T4",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato3_physeq_transformed_T4, distance="bray", max.dist=0.4)
      potato3_T4_net <- plot_network(ig, potato3_physeq_transformed_T4, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### Anthoceros by time
      potato3_physeq_transformed_Antho <- subset_samples(potato3_physeq_transformed_plantsonly, Type=="Anthoceros")
      potato3_physeq_transformed_Antho.ord <- ordinate(potato3_physeq_transformed_Antho, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_Antho, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_Antho)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_Antho <- plot_ordination(potato3_physeq_transformed_Antho, potato3_physeq_transformed_Antho.ord, type="samples", color="Time") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q3 Anthoceros",caption=permanova_des) #+ coord_fixed() 
    ### Notothylas by time
      potato3_physeq_transformed_Noto <- subset_samples(potato3_physeq_transformed_plantsonly, Type=="Notothylas")
      potato3_physeq_transformed_Noto.ord <- ordinate(potato3_physeq_transformed_Noto, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_Noto, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_Noto)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_Noto <- plot_ordination(potato3_physeq_transformed_Noto, potato3_physeq_transformed_Noto.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q3 Notothylas",caption=permanova_des) #+ coord_fixed()  
    ### Phaeoceros by time
      potato3_physeq_transformed_Phaeo <- subset_samples(potato3_physeq_transformed_plantsonly, Type=="Phaeoceros")
      potato3_physeq_transformed_Phaeo.ord <- ordinate(potato3_physeq_transformed_Phaeo, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_Phaeo, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_soil)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_Phaeo <- plot_ordination(potato3_physeq_transformed_Phaeo, potato3_physeq_transformed_Phaeo.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q3 Phaeoceros",caption=permanova_des) #+ coord_fixed() 
    ### soil by time
      potato3_physeq_transformed_soil <- subset_samples(potato3_physeq_transformed, Type=="soil")
      potato3_physeq_transformed_soil.ord <- ordinate(potato3_physeq_transformed_soil, "PCoA", "bray")
      potato3_bray_dist <- phyloseq::distance(potato3_physeq_transformed_soil, method="bray")
      permanova <- adonis(potato3_bray_dist ~ sample_data(potato3_physeq_transformed_soil)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato3_pcoa_soil <- plot_ordination(potato3_physeq_transformed_soil, potato3_physeq_transformed_soil.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q3 soil",caption=permanova_des) #+ coord_fixed() 
    
    (potato3_pcoa_T1 | potato3_pcoa_T2 | potato3_pcoa_T3 | potato3_pcoa_T4) + plot_layout(guides = 'collect')
    (potato3_pcoa_Antho | potato3_pcoa_Noto | potato3_pcoa_Phaeo | potato3_pcoa_soil)+ plot_layout(guides = 'collect')
    
    ### core abundance 
      potato3_Tall <- NULL
      for (time in c('T1','T2','T3','T4')) {
        table1 <- phyloseq2core(subset_samples(potato3_physeq_transformed_Antho, Time==time), 'Anthoceros', 'Potato3', time)
        table2 <- phyloseq2core(subset_samples(potato3_physeq_transformed_Noto, Time==time), 'Notothylas', 'Potato3', time)
        table3 <- phyloseq2core(subset_samples(potato3_physeq_transformed_Phaeo, Time==time), 'Phaeoceros', 'Potato3', time)
        potato3_Tall <- bind_rows(potato3_Tall, table1, table2, table3)
      }
      #potato3_core_heatmap <- as.ggplot(pheatmap(potato3_Tall,cluster_rows=FALSE, cluster_cols=FALSE, fontsize=8))
      potato3_core_stickgraph <- ggplot(data=potato3_Tall) + geom_point(aes(y = reorder(ASV, desc(ASV)), x=abundance, shape=Time, color=Type)) +
        xlab('Relative abundance') + ylab('') + scale_color_manual(values=c("steelblue","gold2","darkseagreen"))
      potato3_core_heatmap <- as.ggplot( pheatmap(tbl2hmap(tbl.long = potato3_Tall, rowVar = 'ASV', colVar = 'replicate', valueVar = 'value')$mat ,cluster_rows=FALSE, cluster_cols=FALSE) )
      
  ## ^Quadrat 1 ====
    ### transform counts into proportion
      potato1_physeq_transformed <- subset_samples(potato_physeq_transformed, Quadrat=="Potato1")
    ### PERMANOVA: soil + plants
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed)$Type)
      permanova_des_all <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      pairwise.perm.manova(potato1_bray_dist,sample_data(potato1_physeq_transformed)$Type,nperm=999)
    ### PERMANOVA: plants only
      potato1_physeq_transformed_plantsonly <- subset_samples(potato1_physeq_transformed, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_plantsonly, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_plantsonly)$Type, permutation=10000)
      permanova_des_plants <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      pairwise.perm.manova(potato1_bray_dist,sample_data(potato1_physeq_transformed_plantsonly)$Type,nperm=999)
    ### PERMANOVA: plants only Unifrac
      potato1_unifrac_dist <- phyloseq::distance(potato1_physeq_transformed_plantsonly, method="wunifrac")
      permanova <- adonis(potato1_unifrac_dist ~ sample_data(potato1_physeq_transformed_plantsonly)$Type, permutation=10000)
      permanova_des_plants_unifrac <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
    ### PERMANOVA: soil + plants by time
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed, method="wunifrac")
      adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed)$Type + sample_data(potato1_physeq_transformed)$Time, permutations = 10000)
      adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed)$Time, permutations = 10000)
    ### ordination plot: soil + plants
      potato1_physeq_transformed.ord <- ordinate(potato1_physeq_transformed, "PCoA", "bray")
      potato1_pcoa_all <- plot_ordination(potato1_physeq_transformed, potato1_physeq_transformed.ord, type="samples", color="Type") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 plant+soil",caption=permanova_des_all) #+ coord_fixed() 
    ### ordination plot: plants only
      potato1_physeq_transformed_plantsonly.ord <- ordinate(potato1_physeq_transformed_plantsonly, "PCoA", "bray")
      potato1_pcoa <- plot_ordination(potato1_physeq_transformed_plantsonly, potato1_physeq_transformed_plantsonly.ord, type="samples", color="Type") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 plant",caption=permanova_des_plants) #+ coord_fixed() 
    ### ordination plot: plants only Unifrac
      potato1_physeq_transformed_plantsonly_unifrac.ord <- ordinate(potato1_physeq_transformed_plantsonly, "PCoA", "wunifrac")
      potato1_pcoa_unifrac <- plot_ordination(potato1_physeq_transformed_plantsonly, potato1_physeq_transformed_plantsonly_unifrac.ord, type="samples", color="Type") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 plant",caption=permanova_des_plants_unifrac) #+ coord_fixed() 
    ### ordination plot: soil + plants by time
      potato1_physeq_transformed.ord <- ordinate(potato1_physeq_transformed, "PCoA", "wunifrac")
      potato1_pcoa_byT <- plot_ordination(potato1_physeq_transformed, potato1_physeq_transformed.ord, type="samples", color="Time") +
        geom_point(size=4) + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q1 by Time") + facet_wrap(~Type, nrow=1) #+ coord_fixed()  
      plot_ordination(potato1_physeq_transformed, potato1_physeq_transformed.ord, type="samples", color="Time") +
        geom_point(size=4) + scale_color_brewer(palette="YlGnBu")
    ### network: soil + plants
      plot_net(potato1_physeq_transformed, distance = "bray", maxdist = 0.2, color = "Type", laymeth="fruchterman.reingold")
      ig <- make_network(potato1_physeq_transformed, distance="bray", max.dist=0.4)
      potato1_net_all <- plot_network(ig, potato1_physeq_transformed, color = "Type", label=NULL, point_size=4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### network: plants only
      plot_net(potato1_physeq_transformed_plantsonly, distance = "bray", maxdist = 0.2, color = "Type", laymeth="fruchterman.reingold")
      ig <- make_network(potato1_physeq_transformed_plantsonly, distance="bray", max.dist=0.4)
      potato1_net <- plot_network(ig, potato1_physeq_transformed_plantsonly, color = "Type", label=NULL, point_size = 4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### network: plants only Unifract
      ig <- make_network(potato1_physeq_transformed_plantsonly, distance="wunifrac", max.dist=0.1)
      potato1_net_unifrac <- plot_network(ig, potato1_physeq_transformed_plantsonly, color = "Type", label=NULL, point_size = 4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### T1
      potato1_physeq_transformed_T1 <- subset_samples(potato1_physeq_transformed_plantsonly, Time=="T1")
      potato1_physeq_transformed_T1.ord <- ordinate(potato1_physeq_transformed_T1, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_T1, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_T1)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_T1 <- plot_ordination(potato1_physeq_transformed_T1, potato1_physeq_transformed_T1.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 T1",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato1_physeq_transformed_T1, distance="bray", max.dist=0.4)
      potato1_T1_net <- plot_network(ig, potato1_physeq_transformed_T1, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T2
      potato1_physeq_transformed_T2 <- subset_samples(potato1_physeq_transformed_plantsonly, Time=="T2")
      potato1_physeq_transformed_T2.ord <- ordinate(potato1_physeq_transformed_T2, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_T2, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_T2)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_T2 <- plot_ordination(potato1_physeq_transformed_T2, potato1_physeq_transformed_T2.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 T2",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato1_physeq_transformed_T2, distance="bray", max.dist=0.4)
      potato1_T2_net <- plot_network(ig, potato1_physeq_transformed_T2, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T3
      potato1_physeq_transformed_T3 <- subset_samples(potato1_physeq_transformed_plantsonly, Time=="T3")
      potato1_physeq_transformed_T3.ord <- ordinate(potato1_physeq_transformed_T3, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_T3, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_T3)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_T3 <- plot_ordination(potato1_physeq_transformed_T3, potato1_physeq_transformed_T3.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 T3",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato1_physeq_transformed_T3, distance="bray", max.dist=0.4)
      potato1_T3_net <- plot_network(ig, potato1_physeq_transformed_T3, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T4
      potato1_physeq_transformed_T4 <- subset_samples(potato1_physeq_transformed_plantsonly, Time=="T4")
      potato1_physeq_transformed_T4.ord <- ordinate(potato1_physeq_transformed_T4, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_T4, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_T4)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_T4 <- plot_ordination(potato1_physeq_transformed_T4, potato1_physeq_transformed_T4.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q1 T4",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato1_physeq_transformed_T4, distance="bray", max.dist=0.4)
      potato1_T4_net <- plot_network(ig, potato1_physeq_transformed_T4, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()      
    ### Anthoceros by time
      potato1_physeq_transformed_Antho <- subset_samples(potato1_physeq_transformed_plantsonly, Type=="Anthoceros")
      potato1_physeq_transformed_Antho.ord <- ordinate(potato1_physeq_transformed_Antho, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_Antho, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_Antho)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_Antho <- plot_ordination(potato1_physeq_transformed_Antho, potato1_physeq_transformed_Antho.ord, type="samples", color="Time") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q1 Anthoceros",caption=permanova_des) #+ coord_fixed() 
    ### Notothylas by time
      potato1_physeq_transformed_Noto <- subset_samples(potato1_physeq_transformed_plantsonly, Type=="Notothylas")
      potato1_physeq_transformed_Noto.ord <- ordinate(potato1_physeq_transformed_Noto, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_Noto, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_Noto)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_Noto <- plot_ordination(potato1_physeq_transformed_Noto, potato1_physeq_transformed_Noto.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q1 Notothylas",caption=permanova_des) #+ coord_fixed()  
    ### Phaeoceros by time
      potato1_physeq_transformed_Phaeo <- subset_samples(potato1_physeq_transformed_plantsonly, Type=="Phaeoceros")
      potato1_physeq_transformed_Phaeo.ord <- ordinate(potato1_physeq_transformed_Phaeo, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_Phaeo, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_Phaeo)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_Phaeo <- plot_ordination(potato1_physeq_transformed_Phaeo, potato1_physeq_transformed_Phaeo.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q1 Phaeoceros",caption=permanova_des) #+ coord_fixed() 
    ### soil by time
      potato1_physeq_transformed_soil <- subset_samples(potato1_physeq_transformed, Type=="soil")
      potato1_physeq_transformed_soil.ord <- ordinate(potato1_physeq_transformed_soil, "PCoA", "bray")
      potato1_bray_dist <- phyloseq::distance(potato1_physeq_transformed_soil, method="bray")
      permanova <- adonis(potato1_bray_dist ~ sample_data(potato1_physeq_transformed_soil)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato1_pcoa_soil <- plot_ordination(potato1_physeq_transformed_soil, potato1_physeq_transformed_soil.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q1 soil",caption=permanova_des) #+ coord_fixed() 
      
    (potato1_pcoa_T1 | potato1_pcoa_T2 | potato1_pcoa_T3 | potato1_pcoa_T4) + plot_layout(guides = 'collect')
    (potato1_pcoa_Antho | potato1_pcoa_Noto | potato1_pcoa_Phaeo | potato1_pcoa_soil)+ plot_layout(guides = 'collect')
     
    ### core abundance 
      potato1_Tall <- NULL
      for (time in c('T1','T2','T3','T4')) {
        table1 <- phyloseq2core(subset_samples(potato1_physeq_transformed_Antho, Time==time), 'Anthoceros', 'Potato1', time)
        if (time != 'T4'){
        table2 <- phyloseq2core(subset_samples(potato1_physeq_transformed_Noto, Time==time), 'Notothylas', 'Potato1', time)
        }
        else {table2 <- NULL}
        table3 <- phyloseq2core(subset_samples(potato1_physeq_transformed_Phaeo, Time==time), 'Phaeoceros', 'Potato1', time)
        potato1_Tall <- bind_rows(potato1_Tall, table1, table2, table3)
      }
      #potato1_core_heatmap <- as.ggplot(pheatmap(potato1_Tall,cluster_rows=FALSE, cluster_cols=FALSE, fontsize=8))
      potato1_core_stickgraph <- ggplot(data=potato1_Tall) + geom_point(aes(y = reorder(ASV, desc(ASV)), x=abundance, shape=Time, color=Type)) +
        xlab('Relative abundance') + ylab('') + scale_color_manual(values=c("steelblue","gold2","darkseagreen"))
      potato1_core_heatmap <- as.ggplot( pheatmap(tbl2hmap(tbl.long = potato1_Tall, rowVar = 'ASV', colVar = 'replicate', valueVar = 'value')$mat ,cluster_rows=FALSE, cluster_cols=FALSE) )
      
  ## ^Quadrat 2 ====
    ### transform counts into proportion
      potato2_physeq_transformed <- subset_samples(potato_physeq_transformed, Quadrat=="Potato2")
    ### PERMANOVA: soil + plants
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed)$Type, permutations = 10000)
      permanova_des_all <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      pairwise.perm.manova(potato2_bray_dist,sample_data(potato2_physeq_transformed)$Type,nperm=10000)
    ### PERMANOVA: plants only
      potato2_physeq_transformed_plantsonly <- subset_samples(potato2_physeq_transformed, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_plantsonly, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_plantsonly)$Type, permutations = 10000)
      permanova_des_plants <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      pairwise.perm.manova(potato2_bray_dist,sample_data(potato2_physeq_transformed_plantsonly)$Type,nperm=999)
    ### PERMANOVA: plants only Unifrac
      potato2_unifrac_dist <- phyloseq::distance(potato2_physeq_transformed_plantsonly, method="wunifrac")
      permanova <- adonis(potato2_unifrac_dist ~ sample_data(potato2_physeq_transformed_plantsonly)$Type, permutations = 10000)
      permanova_des_plants_unifrac <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
    ### PERMANOVA: soil + plants by time 
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed, method="wunifrac")
      adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed)$Type + sample_data(potato2_physeq_transformed)$Time, permutations = 10000)
      adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed)$Time, permutations = 10000)
    ### ANOSIM: test
      var_group <- get_variable(potato2_physeq_transformed_plantsonly, "Type")
      anosim(distance(potato2_physeq_transformed_plantsonly, "bray"), var_group, permutations = 10000)    
    ### ordination plot: soil + plants
      potato2_physeq_transformed.ord <- ordinate(potato2_physeq_transformed, "PCoA", "bray")
      potato2_pcoa_all <- plot_ordination(potato2_physeq_transformed, potato2_physeq_transformed.ord, type="samples", color="Type") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 plant+soil",caption=permanova_des_all) #+ coord_fixed() 
    ### ordination plot: plants only
      #potato2_physeq_transformed_plantsonly.ord <- ordinate(potato2_physeq_transformed_plantsonly, "PCoA", "unifrac", weighted=TRUE)
      potato2_physeq_transformed_plantsonly.ord <- ordinate(potato2_physeq_transformed_plantsonly, "PCoA", "bray")
      potato2_pcoa <- plot_ordination(potato2_physeq_transformed_plantsonly, potato2_physeq_transformed_plantsonly.ord, type="samples", color="Type") + #+ geom_polygon(aes(fill=Type))
        geom_point(size=4) + scale_colour_manual(values = cbp1) + coord_fixed()
      potato2_pcoa <- plot_ordination(potato2_physeq_transformed_plantsonly, potato2_physeq_transformed_plantsonly.ord, type="samples", color="Type") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 plant",caption=permanova_des_plants) #+ coord_fixed() 
    ### ordination plot: plants only Unifrac
      potato2_physeq_transformed_plantsonly_unifrac.ord <- ordinate(potato2_physeq_transformed_plantsonly, "PCoA", "wunifrac")
      potato2_pcoa_unifrac <- plot_ordination(potato2_physeq_transformed_plantsonly, potato2_physeq_transformed_plantsonly_unifrac.ord, type="samples", color="Type") + 
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 plant",caption=permanova_des_plants_unifrac) #+ coord_fixed() 
    ### ordination plot: soil + plants by time
      potato2_physeq_transformed_unifrac.ord <- ordinate(potato2_physeq_transformed, "PCoA", "wunifrac")
      potato2_pcoa_byT <- plot_ordination(potato2_physeq_transformed, potato2_physeq_transformed_unifrac.ord, type="samples", color="Time") +
        geom_point(size=4) + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q2 by time") + facet_wrap(~Type, nrow=1) #+ coord_fixed()  
    ### network: soil + plants
      plot_net(potato2_physeq_transformed, distance = "bray", maxdist = 0.2, color = "Type", laymeth="fruchterman.reingold")
      ig <- make_network(potato2_physeq_transformed, distance="bray", max.dist=0.4)
      potato2_net_all <- plot_network(ig, potato2_physeq_transformed, color = "Type", label=NULL, point_size=4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### network: plants only
      plot_net(potato2_physeq_transformed_plantsonly, distance = "bray", maxdist = 0.2, color = "Type", laymeth="fruchterman.reingold")
      ig <- make_network(potato2_physeq_transformed_plantsonly, distance="bray", max.dist=0.4)
      potato2_net <- plot_network(ig, potato2_physeq_transformed_plantsonly, color = "Type", label=NULL, point_size = 4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### network: plants only Unifrac
      ig <- make_network(potato2_physeq_transformed_plantsonly, distance="wunifrac", max.dist=0.1)
      potato2_net_unifrac <- plot_network(ig, potato2_physeq_transformed_plantsonly, color = "Type", label=NULL, point_size = 4) +
        geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) #+ coord_fixed()
    ### T1
      potato2_physeq_transformed_T1 <- subset_samples(potato2_physeq_transformed_plantsonly, Time=="T1")
      potato2_physeq_transformed_T1.ord <- ordinate(potato2_physeq_transformed_T1, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_T1, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_T1)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_T1 <- plot_ordination(potato2_physeq_transformed_T1, potato2_physeq_transformed_T1.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 T1",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato2_physeq_transformed_T1, distance="bray", max.dist=0.4)
      potato2_T1_net <- plot_network(ig, potato2_physeq_transformed_T1, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T2
      potato2_physeq_transformed_T2 <- subset_samples(potato2_physeq_transformed_plantsonly, Time=="T2")
      potato2_physeq_transformed_T2.ord <- ordinate(potato2_physeq_transformed_T2, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_T2, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_T2)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_T2 <- plot_ordination(potato2_physeq_transformed_T2, potato2_physeq_transformed_T2.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 T2",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato2_physeq_transformed_T2, distance="bray", max.dist=0.4)
      potato2_T2_net <- plot_network(ig, potato2_physeq_transformed_T2, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T3
      potato2_physeq_transformed_T3 <- subset_samples(potato2_physeq_transformed_plantsonly, Time=="T3")
      potato2_physeq_transformed_T3.ord <- ordinate(potato2_physeq_transformed_T3, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_T3, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_T3)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_T3 <- plot_ordination(potato2_physeq_transformed_T3, potato2_physeq_transformed_T3.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 T3",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato2_physeq_transformed_T3, distance="bray", max.dist=0.4)
      potato2_T3_net <- plot_network(ig, potato2_physeq_transformed_T3, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()
    ### T4
      potato2_physeq_transformed_T4 <- subset_samples(potato2_physeq_transformed_plantsonly, Time=="T4")
      potato2_physeq_transformed_T4.ord <- ordinate(potato2_physeq_transformed_T4, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_T4, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_T4)$Type)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_T4 <- plot_ordination(potato2_physeq_transformed_T4, potato2_physeq_transformed_T4.ord, type="samples", color="Type") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_colour_manual(values = cbp1) + labs(title="Potato Q2 T4",caption=permanova_des) #+ coord_fixed() 
      ig <- make_network(potato2_physeq_transformed_T4, distance="bray", max.dist=0.4)
      potato2_T4_net <- plot_network(ig, potato2_physeq_transformed_T4, color = "Type", label=NULL, point_size=4) +
        scale_colour_manual(values = cbp1) + coord_fixed()      
      
    ### Anthoceros by time
      potato2_physeq_transformed_Antho <- subset_samples(potato2_physeq_transformed_plantsonly, Type=="Anthoceros")
      potato2_physeq_transformed_Antho.ord <- ordinate(potato2_physeq_transformed_Antho, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_Antho, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_Antho)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_Antho <- plot_ordination(potato2_physeq_transformed_Antho, potato2_physeq_transformed_Antho.ord, type="samples", color="Time") +
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q2 Anthoceros",caption=permanova_des) #+ coord_fixed() 
    ### Notothylas by time
      potato2_physeq_transformed_Noto <- subset_samples(potato2_physeq_transformed_plantsonly, Type=="Notothylas")
      potato2_physeq_transformed_Noto.ord <- ordinate(potato2_physeq_transformed_Noto, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_Noto, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_Noto)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_Noto <- plot_ordination(potato2_physeq_transformed_Noto, potato2_physeq_transformed_Noto.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q2 Notothylas",caption=permanova_des) #+ coord_fixed()  
    ### Phaeoceros by time
      potato2_physeq_transformed_Phaeo <- subset_samples(potato2_physeq_transformed_plantsonly, Type=="Phaeoceros")
      potato2_physeq_transformed_Phaeo.ord <- ordinate(potato2_physeq_transformed_Phaeo, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_Phaeo, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_Phaeo)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_Phaeo <- plot_ordination(potato2_physeq_transformed_Phaeo, potato2_physeq_transformed_Phaeo.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q2 Phaeoceros",caption=permanova_des) #+ coord_fixed() 
    ### soil by time
      potato2_physeq_transformed_soil <- subset_samples(potato2_physeq_transformed, Type=="soil")
      potato2_physeq_transformed_soil.ord <- ordinate(potato2_physeq_transformed_soil, "PCoA", "bray")
      potato2_bray_dist <- phyloseq::distance(potato2_physeq_transformed_soil, method="bray")
      permanova <- adonis(potato2_bray_dist ~ sample_data(potato2_physeq_transformed_soil)$Time)
      permanova_des <- paste("permanova ","R2=",substr(as.character(permanova$aov.tab$R2[1]),1,5)," ","p=",as.character(permanova$aov.tab$`Pr(>F)`[1]), sep = "")
      potato2_pcoa_soil <- plot_ordination(potato2_physeq_transformed_soil, potato2_physeq_transformed_soil.ord, type="samples", color="Time")+
        geom_point(size=4) + geom_point(shape = 1,size = 4,colour = "black") + scale_color_brewer(palette="YlGnBu") + labs(title="Potato Q2 soil",caption=permanova_des) #+ coord_fixed() 
    
    (potato2_pcoa_T1 | potato2_pcoa_T2 | potato2_pcoa_T3 | potato2_pcoa_T4) + plot_layout(guides = 'collect')
    (potato2_pcoa_Antho | potato2_pcoa_Noto | potato2_pcoa_Phaeo | potato2_pcoa_soil)+ plot_layout(guides = 'collect')
    
    ### core abundance 
      potato2_Tall <- NULL
      for (time in c('T1','T2','T3','T4')) {
        table1 <- phyloseq2core(subset_samples(potato2_physeq_transformed_Antho, Time==time), 'Anthoceros', 'Potato2', time)
        table2 <- phyloseq2core(subset_samples(potato2_physeq_transformed_Noto, Time==time), 'Notothylas', 'Potato2', time)
        table3 <- phyloseq2core(subset_samples(potato2_physeq_transformed_Phaeo, Time==time), 'Phaeoceros', 'Potato2', time)
        potato2_Tall <- bind_rows(potato2_Tall, table1, table2, table3)
      }
      #potato2_core_heatmap <- as.ggplot(pheatmap(potato2_Tall,cluster_rows=FALSE, cluster_cols=FALSE, fontsize=8))
      potato2_core_stickgraph <- ggplot(data=potato2_Tall) + geom_point(aes(y = reorder(ASV, desc(ASV)), x=abundance, shape=Time, color=Type)) +
        xlab('Relative abundance') + ylab('') + scale_color_manual(values=c("steelblue","gold2","darkseagreen"))
      potato2_core_heatmap <- as.ggplot( pheatmap(tbl2hmap(tbl.long = potato2_Tall, rowVar = 'ASV', colVar = 'replicate', valueVar = 'value')$mat ,cluster_rows=FALSE, cluster_cols=FALSE) )
      
# Save Plots --------------------------------------------------------------
  potato_pcoa_all_soil + potato_pcoa_all_plantsonly + plot_layout(guides = 'collect')
  ggsave("potato_pcoa_plants_soils.pdf", device = "pdf", width = 18, height = 12)
    
  grossman_pcoa_all + potato_pcoa_all_soil + potato_pcoa_all_plantsonly + plot_layout(guides = 'collect')
  ggsave("potato_grossman_pcoa_plants_soils.pdf", device = "pdf", width = 18, height = 12)
  ggsave("potato_grossman_pcoa_plants_soils.svg", device = "svg", width = 18, height = 12)
  
  (potato1_pcoa | potato2_pcoa | potato3_pcoa)/(potato1_net | potato2_net | potato3_net) + plot_layout(guides = 'collect')
  ggsave("potato_pcoa_plants_byQ.pdf", device = "pdf", width = 18, height = 12)
  ggsave("potato_pcoa_plants_byQ.svg", device = "svg", width = 18, height = 12)

  (potato1_pcoa_unifrac | potato2_pcoa_unifrac | potato3_pcoa_unifrac) + plot_layout(guides = 'collect')
  ggsave("potato_pcoa_plants_byQ_unifrac.pdf", device = "pdf", width = 18, height = 6)
  ggsave("potato_pcoa_plants_byQ_unifrac.svg", device = "svg", width = 18, height = 6)
  
  (potato1_pcoa_all | potato2_pcoa_all | potato3_pcoa_all)/(potato1_net_all | potato2_net_all | potato3_net_all) + plot_layout(guides = 'collect')
  ggsave("potato_pcoa_plants_soils_byQ.pdf", device = "pdf", width = 18, height = 12)
  
  (potato1_pcoa_T1 | potato1_pcoa_T2 | potato1_pcoa_T3 | potato1_pcoa_T4) / 
    (potato2_pcoa_T1 | potato2_pcoa_T2 | potato2_pcoa_T3 | potato2_pcoa_T4) /
    (potato3_pcoa_T1 | potato3_pcoa_T2 | potato3_pcoa_T3 | potato3_pcoa_T4) + plot_layout(guides = 'collect')
  ggsave("potato_pcoa_plants_soils_byQ_byT.pdf", device = "pdf", width = 18, height = 12)
  
  (grossman1_pcoa_Noto | grossman1_pcoa_soil | grossman2_pcoa_Noto | grossman2_pcoa_soil) /
    (potato1_pcoa_Antho | potato1_pcoa_Noto | potato1_pcoa_Phaeo | potato1_pcoa_soil) /
    (potato2_pcoa_Antho | potato2_pcoa_Noto | potato2_pcoa_Phaeo | potato2_pcoa_soil) /
    (potato3_pcoa_Antho | potato3_pcoa_Noto | potato3_pcoa_Phaeo | potato3_pcoa_soil) + plot_layout(guides = 'collect')
  ggsave("potato_grossman_pcoa_plants_soils_byQ_byType.pdf", device = "pdf", width = 18, height = 17)
  ggsave("potato_grossman_pcoa_plants_soils_byQ_byType.svg", device = "svg", width = 18, height = 17)

  ((grossman1_pcoa_all | grossman2_pcoa_all) + plot_layout(guides = 'collect')) / potato1_pcoa_byT / potato2_pcoa_byT / potato3_pcoa_byT
  ggsave("potato_grossman_pcoa_plants_soils_byQ_byTime.pdf", device = "pdf", width = 18, height = 17)
  ggsave("potato_grossman_pcoa_plants_soils_byQ_byTime.svg", device = "svg", width = 18, height = 17)
  
  ((grossman_core_stickgraph / guide_area() ) | potato1_core_stickgraph | potato2_core_stickgraph | potato3_core_stickgraph) + plot_layout(guides = 'collect')
  core_Tall <- bind_rows(grossman_Tall, potato1_Tall, potato2_Tall, potato3_Tall)
  core_Tall_heatmap <- as.ggplot( pheatmap(tbl2hmap(tbl.long = core_Tall, rowVar = 'ASV', colVar = 'replicate', valueVar = 'value')$mat ,cluster_rows=FALSE, cluster_cols=FALSE) )
  ggsave("core_Tall_heatmap.pdf", device = "pdf", width = 4, height = 18)
  
# DESeq2 on Q3 ------------------------------------------------------------
  fdr <- 0.05
  fold2changethres <- 0
  potato3_physeq <- subset_samples(potato_physeq, Quadrat=="Potato3")
  potato3_physeq_plantsonly <- subset_samples(potato3_physeq, Type!="soil")
  #potato3_physeq_plantsonly <- subset_samples(potato_physeq_plantsonly, Quadrat=="Potato3")
  dds <- phyloseq_to_deseq2(potato3_physeq_plantsonly, ~ Type)
  dds <- DESeq(dds, test="Wald", fitType="parametric", sfType="poscounts")
  res <- results(dds, contrast = c("Type","Anthoceros","Notothylas"),alpha=fdr)
  resDF <- cbind(counts(dds, normalized=TRUE), as.data.frame(res))
  resDFSig <- subset(resDF, padj<fdr & abs(log2FoldChange)>=fold2changethres)
  p1 <- plot_ASV('ASV1', potato3_physeq_transformed)
  p2 <- plot_ASV('ASV36', potato3_physeq_transformed)
  p3 <- plot_ASV('ASV42', potato3_physeq_transformed)
  p4 <- plot_ASV('ASV45', potato3_physeq_transformed)
  p5 <- plot_ASV('ASV52', potato3_physeq_transformed)
  p8 <- plot_ASV_deseq('ASV1')
  p9 <- plot_ASV_deseq('ASV36')
  p10 <- plot_ASV_deseq('ASV42')
  p11 <- plot_ASV_deseq('ASV45')
  p12 <- plot_ASV_deseq('ASV52')

  (p1 | p2 | p3 | p4| p5) /(p8 | p9 | p10 | p11| p12)+ plot_layout(guides = 'collect')
  
  ### MDS plot from DESeq2
  sampleDists <- dist(t(assay(dds)))
  sampleDistMatrix <- as.matrix(sampleDists)
  mds <- as.data.frame(colData(dds)) %>% cbind(cmdscale(sampleDistMatrix))
  ggplot(mds, aes(x = `1`, y = `2`, color=Type)) +
    geom_point(size = 3) + coord_fixed()

# ANCOM ------------------------------------------------------------
  dataset <- subset_samples(potato_physeq, Type!="soil")
  dataset <- subset_samples(potato_physeq_filtered, Type!="soil")
  ## ^ANCOM on Q3 ====
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato3")), 
              meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato3")), 
              out_cut = 0.05, zero_cut = 0.90, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Type", alpha = 0.05)
    ANCOM_res_time <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Time", alpha = 0.05)
    p <- guide_area() 
    for (x in as.vector(subset(ANCOM_res$out, detected_0.7==TRUE)$taxa_id)) {
      print(x)
      print(plot_ASV(x, potato3_physeq_transformed))
      p = p + plot_ASV(x, potato3_physeq_transformed)
    }
    p = p + plot_layout(guides = 'collect')
    ## make better plot
    p_ASV1 <- plot_ASV("ASV1", potato3_physeq_transformed_plantsonly, violin = T)
    p_ASV8 <- plot_ASV("ASV8", potato3_physeq_transformed_plantsonly, violin = T)
    p_ASV1_soil <- plot_ASV_soil("ASV1", potato_physeq_transformed, violin = T)
    p_ASV8_soil <- plot_ASV_soil("ASV8", potato_physeq_transformed, violin = T)
    (p_ASV1 | p_ASV8) / (p_ASV1_soil | p_ASV8_soil) + plot_layout(guides = 'collect')
    (p_ASV1 | p_ASV1_soil) / (p_ASV8 | p_ASV8_soil) + plot_layout(guides = 'collect')
    ggsave("potato_Q3_plantsonly_ANCOM_0.90_v2.pdf", device = "pdf", width = 6, height = 5)
    ggsave("potato_Q3_plantsonly_ANCOM_0.90_v2.svg", device = "svg", width = 9, height = 4)
    
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato3")), 
              meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato3")), 
              out_cut = 0.05, zero_cut = 0.95, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Type", alpha = 0.05)
    p2 <- guide_area() 
    for (x in as.vector(subset(ANCOM_res$out, detected_0.7==TRUE)$taxa_id)) {
      print(x)
      print(plot_ASV(x, potato3_physeq_transformed))
      p2 = p2 + plot_ASV(x, potato3_physeq_transformed)
    }
    p2 = p2 + plot_layout(guides = 'collect')
    ggsave("potato_Q3_plantsonly_ANCOM_0.95.pdf", device = "pdf", width = 18, height = 12)
    
  ## ^ANCOM on Q2 ====
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato2")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato2")), 
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Type", alpha = 0.05)
    ANCOM_res_time <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Time", alpha = 0.05)
    p <- guide_area() 
    for (x in as.vector(subset(ANCOM_res$out, detected_0.7==TRUE)$taxa_id)) {
      print(x)
      print(plot_ASV(x, potato2_physeq_transformed))
      p = p + plot_ASV(x, potato2_physeq_transformed)
    }
    p = p + plot_layout(guides = 'collect')
    ggsave("potato_Q2_plantsonly_ANCOM_0.90.pdf", device = "pdf", width = 18, height = 12)
    
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato2")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato2")), 
                                       out_cut = 0.05, zero_cut = 0.95, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Type", alpha = 0.05)
    p2 <- guide_area() 
    for (x in as.vector(subset(ANCOM_res$out, detected_0.7==TRUE)$taxa_id)) {
      print(x)
      print(plot_ASV(x, potato2_physeq_transformed))
      p2 = p2 + plot_ASV(x, potato2_physeq_transformed)
    }
    p2 = p2 + plot_layout(guides = 'collect')
    ggsave("potato_Q2_plantsonly_ANCOM_0.95.pdf", device = "pdf", width = 18, height = 12)
  
  ## ^ANCOM on Q1 ====
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato1")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato1")), 
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Type", alpha = 0.05)
    ANCOM_res_time <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Time", alpha = 0.05)
    p <- guide_area() 
    for (x in as.vector(subset(ANCOM_res$out, detected_0.7==TRUE)$taxa_id)) {
      print(x)
      print(plot_ASV(x, potato1_physeq_transformed))
      p = p + plot_ASV(x, potato1_physeq_transformed)
    }
    p = p + plot_layout(guides = 'collect')
    ggsave("potato_Q1_ANCOM_0.90.pdf", device = "pdf", width = 18, height = 12)
    
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato1")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato1")), 
                                       out_cut = 0.05, zero_cut = 0.95, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Type", alpha = 0.05)
    p2 <- guide_area() 
    for (x in as.vector(subset(ANCOM_res$out, detected_0.7==TRUE)$taxa_id)) {
      print(x)
      print(plot_ASV(x, potato1_physeq_transformed))
      p2 = p2 + plot_ASV(x, potato1_physeq_transformed)
    }
    p2 = p2 + plot_layout(guides = 'collect')
    ggsave("potato_Q1_ANCOM_0.95.pdf", device = "pdf", width = 18, height = 12)
    
    dataset <- subset_samples(potato_physeq, Type=="soil")
    ## ^ANCOM on Q3 soil by time ====
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato3")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato3")), 
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Time", alpha = 0.05)
    ## ^ANCOM on Q2 soil by time ====
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato2")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato2")), 
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Time", alpha = 0.05)
    ## ^ANCOM on Q1 soil by time ====
    prepro = feature_table_pre_process(feature_table=otu_table(subset_samples(dataset,Quadrat=="Potato1")), 
                                       meta_data=sample_data(subset_samples(dataset,Quadrat=="Potato1")), 
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut = 500, neg_lb=FALSE)
    ANCOM_res <- ANCOM(feature_table=prepro$feature_table, meta_data=prepro$meta_data, main_var="Time", alpha = 0.05)
    
# Mock ------------------------------------------------------------
  # Read ASV table
    ASV_mock_file <- as.matrix(read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/mock/ASV_Table_mock.csv", header=TRUE, sep = ",", row.names = 1))
  # Read sample metadata
    sample_mock_meta_file <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/mock/MockMeta.csv", header=TRUE, sep = ",", row.names = 1)
  # Make phyloseq object
    ASV_table_mock <- otu_table(ASV_mock_file, taxa_are_rows = TRUE)
    sample_meta_mock <- sample_data(sample_mock_meta_file)
    mock_physeq = phyloseq(ASV_table_mock, sample_meta_mock) #, ASV_tree_file)
    mock_physeq_transformed  = transform_sample_counts(mock_physeq, function(x) x / sum(x) )
    #mock_physeq_transformed3 <- transform_sample_counts(mock_physeq, fun = filter5perc_count)  
    #mock_physeq_transformed3 <- transform_sample_counts(mock_physeq_transformed3, function(x) x / sum(x) )
    mock_physeq_transformed_filtered = filter_taxa(mock_physeq_transformed, function(x) sum(x) > 0, TRUE)  
    mock_physeq_transformed_filtered_sub = subset_samples(mock_physeq_transformed_filtered, Run=="JNP1"|Run=="JNP2"|Run=="JNP3")
    mock_physeq_transformed_filtered_sub = filter_taxa(mock_physeq_transformed_filtered_sub, function(x) sum(x) > 0, TRUE) 
  # Make a tax table with just ASV names 
    taxmat <- matrix(rownames(otu_table(mock_physeq_transformed_filtered_sub)))
    rownames(taxmat) <- rownames(otu_table(mock_physeq_transformed_filtered_sub))
    colnames(taxmat) <- "Species"
    mock_physeq_transformed_filtered_sub_tax = merge_phyloseq(mock_physeq_transformed_filtered_sub, tax_table(taxmat))
    p_ASV <- plot_bar(mock_physeq_transformed_filtered_sub_tax, fill="Species") + theme(legend.position = "none") 
    #p_ASV$data$Sample <- factor(p_ASV$data$Sample, levels = list("pos_mock_community_pilot", "mock_community_3taxa_JNP4", "mock_community_5taxa_JNP1", "mock_community_5taxa_JNP2", "mock_community_5taxa_JNP3", "mock_community_5taxa_JNP4"))
    p_ASV
    mock_tree <- read.newick("/Users/fay-weili/Box/hornwort_amplicon/dada2/mock/ASV_mock_ref_namechang_Lib1-3.fasttree.tre")
    p_tree <- ggtree(mock_tree)+ geom_tiplab()
    p_tree <- rotate_tree(p_tree, 90)
  # Chi-square test
    mock_physeq_subset <- subset_samples(mock_physeq, Run=="JNP1"|Run=="JNP2"|Run=="JNP3")
    mock_physeq_subset_filtered <- filter_taxa(mock_physeq_subset, function(x) sum(x) > 0, TRUE) 
    mock_physeq_subset_filtered_tax <- merge_phyloseq(mock_physeq_subset_filtered, tax_table(taxmat))
    mock_physeq_subset_filtered_tax_prune <- subset_taxa(mock_physeq_subset_filtered_tax, Species!="ASV25")
    ## JNP1
      mock_count <- as.data.frame(otu_table(mock_physeq_subset_filtered_tax_prune))$mock_community_5taxa_JNP1
      chi_df <- cbind(mock_count, rep(sum(mock_count) * 0.25, 5))
      chisq.test(chi_df)
    ## JNP2
      mock_count <- as.data.frame(otu_table(mock_physeq_subset_filtered_tax_prune))$mock_community_5taxa_JNP2
      chi_df <- cbind(mock_count, rep(sum(mock_count) * 0.25, 5))
      chisq.test(chi_df)
    ## JNP3
      mock_count <- as.data.frame(otu_table(mock_physeq_subset_filtered_tax_prune))$mock_community_5taxa_JNP3
      chi_df <- cbind(mock_count, rep(sum(mock_count) * 0.25, 5))
      chisq.test(chi_df)
    
  # Read OTU 0.97 table
    otu97 <- as.matrix(read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/sample_fastq_deprimers_lenfiltered_JNP1-3_vsearch_2/all_OTUtable97_mocks.txt", header=TRUE, sep = "\t", row.names = 1))
    otu97_no_singleton <- replace(otu97, otu97 == 1, 0)
  # Read sample metadata
    sample_mock_meta_file_otu <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/sample_fastq_deprimers_lenfiltered_JNP1-3_vsearch/MockMeta.csv", header=TRUE, sep = ",", row.names = 1)
  # Make phyloseq object
    otu97_table_mock <- otu_table(otu97_no_singleton, taxa_are_rows = TRUE)
    sample_meta_mock_otu <- sample_data(sample_mock_meta_file_otu)
    otu97_mock_physeq = phyloseq(otu97_table_mock, sample_meta_mock_otu) #, ASV_tree_file)
    otu97_mock_physeq_transformed  = transform_sample_counts(otu97_mock_physeq, function(x) x / sum(x) )
    otu97_mock_physeq_transformed_sub = subset_samples(otu97_mock_physeq_transformed, Run=="JNP1"|Run=="JNP2"|Run=="JNP3")
  # Make a tax table with just ASV names 
    taxmat <- matrix(rownames(otu_table(otu97_mock_physeq_transformed_sub)))
    rownames(taxmat) <- rownames(otu_table(otu97_mock_physeq_transformed_sub))
    colnames(taxmat) <- "Species"
    otu97_mock_physeq_transformed_tax = merge_phyloseq(otu97_mock_physeq_transformed_sub, tax_table(taxmat))
    otu97_mock_physeq_transformed_tax = filter_taxa(otu97_mock_physeq_transformed_tax, function(x) sum(x) > 0, TRUE)
    p_OTU97 <- plot_bar(otu97_mock_physeq_transformed_tax, fill="Species") + theme(legend.position = "none")
    #p_OTU97$data$Sample <- factor(p_OTU97$data$Sample, levels = list("pos_mock_community_pilot_run", "mock_community_3taxa_JNP4_8_H01", "mock_community_5taxa_JNP1_5_E01", "mock_community_5taxa_JNP2_6_F01", "mock_community_5taxa_JNP3_7_G01", "mock_community_5taxa_JNP4_8_H01"))
    p_OTU97
    colSums(otu_table(otu97_mock_physeq_transformed_tax) != 0)
    
  # Read OTU 0.95 table
    otu95 <- as.matrix(read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/sample_fastq_deprimers_lenfiltered_JNP1-3_vsearch_2/all_OTUtable95_mocks.txt", header=TRUE, sep = "\t", row.names = 1))
    otu95_no_singleton <- replace(otu95, otu95 == 1, 0)
  # Read sample metadata
    sample_mock_meta_file_otu <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/sample_fastq_deprimers_lenfiltered_JNP1-3_vsearch/MockMeta.csv", header=TRUE, sep = ",", row.names = 1)
  # Make phyloseq object
    otu95_table_mock <- otu_table(otu95_no_singleton, taxa_are_rows = TRUE)
    sample_meta_mock_otu <- sample_data(sample_mock_meta_file_otu)
    otu95_mock_physeq = phyloseq(otu95_table_mock, sample_meta_mock_otu) #, ASV_tree_file)
    otu95_mock_physeq_transformed  = transform_sample_counts(otu95_mock_physeq, function(x) x / sum(x) )
    otu95_mock_physeq_transformed_sub = subset_samples(otu95_mock_physeq_transformed, Run=="JNP1"|Run=="JNP2"|Run=="JNP3")
  # Make a tax table with just ASV names 
    taxmat <- matrix(rownames(otu_table(otu95_mock_physeq_transformed_sub)))
    rownames(taxmat) <- rownames(otu_table(otu95_mock_physeq_transformed_sub))
    colnames(taxmat) <- "Species"
    otu95_mock_physeq_transformed_tax = merge_phyloseq(otu95_mock_physeq_transformed_sub, tax_table(taxmat))
    otu95_mock_physeq_transformed_tax = filter_taxa(otu95_mock_physeq_transformed_tax, function(x) sum(x) > 0, TRUE)
    p_OTU95 <- plot_bar(otu95_mock_physeq_transformed_tax, fill="Species") + theme(legend.position = "none")
    #p_OTU95$data$Sample <- factor(p_OTU95$data$Sample, levels = list("pos_mock_community_pilot_run", "mock_community_3taxa_JNP4_8_H01", "mock_community_5taxa_JNP1_5_E01", "mock_community_5taxa_JNP2_6_F01", "mock_community_5taxa_JNP3_7_G01", "mock_community_5taxa_JNP4_8_H01"))
    p_OTU95
    colSums(otu_table(otu95_mock_physeq_transformed_tax) != 0)
    
    p_ASV | p_tree / (p_OTU97 | p_OTU95) 
    p_tree | p_ASV | p_OTU97 | p_OTU95 
    ggsave("mock_comparison2.pdf", device = "pdf", width = 18, height = 12)
  
# Alpha diversity ------------------------------------------------------------
  potato_physeq <- subset_samples(time_series_physeq, Quadrat=="Potato1"|Quadrat=="Potato2"|Quadrat=="Potato3")
  potato_physeq_plantsonly <- subset_samples(potato_physeq, Type=="Anthoceros"|Type=="Notothylas"|Type=="Phaeoceros")
  potato_physeq_plantsonly_filtered <- transform_sample_counts(potato_physeq_plantsonly, fun = filter3perc_count)  
  potato_physeq_filtered <- merge_phyloseq(potato_physeq_plantsonly_filtered, subset_samples(potato_physeq, Type=="soil"))
  grossman_physeq <- subset_samples(time_series_physeq, Quadrat=="Grossman1"|Quadrat=="Grossman2")
  grossman_physeq_plantsonly <- subset_samples(grossman_physeq, Type=="Notothylas")
  grossman_physeq_plantsonly_filtered <- transform_sample_counts(grossman_physeq_plantsonly, fun = filter3perc_count)  
  grossman_physeq_filtered <- merge_phyloseq(grossman_physeq_plantsonly_filtered, subset_samples(grossman_physeq, Type=="soil"))
  phydist <- cophenetic(tree)
  ## ^Plant samples no filtering ====
    ### Grossman
      alphadiv <- estimate_richness(grossman_physeq, measures=c("Shannon","Chao1"))
      pd_df <- pd(t(as.data.frame(otu_table(grossman_physeq))), tree, include.root=T)
      mpd_df <- mpd(t(as.data.frame(otu_table(grossman_physeq))), phydist, abundance.weighted = T)
      mntd_df <- mntd(t(as.data.frame(otu_table(grossman_physeq))), phydist, abundance.weighted = T)
      alphadiv$PD <- pd_df$PD
      alphadiv$MPD <- mpd_df
      alphadiv$MNTD <- mntd_df
      grossman_alphadiv <- merge(alphadiv, sample_data(grossman_physeq), by=0)
      row.names(grossman_alphadiv) <- grossman_alphadiv$Row.names
      grossman_alphadiv[1] <- NULL
    ### Potato Hill
      alphadiv <- estimate_richness(potato_physeq, measures=c("Shannon","Chao1"))
      pd_df <- pd(t(as.data.frame(otu_table(potato_physeq))), tree, include.root=T)
      mpd_df <- mpd(t(as.data.frame(otu_table(potato_physeq))), phydist, abundance.weighted = T)
      mntd_df <- mntd(t(as.data.frame(otu_table(potato_physeq))), phydist, abundance.weighted = T)
      alphadiv$PD <- pd_df$PD
      alphadiv$MPD <- mpd_df
      alphadiv$MNTD <- mntd_df
      potato_alphadiv <- merge(alphadiv, sample_data(potato_physeq), by=0)
      row.names(potato_alphadiv) <- potato_alphadiv$Row.names
      potato_alphadiv[1] <- NULL
    ### Combine the two
      all_alphadiv <- rbind(potato_alphadiv, grossman_alphadiv)
    ### Shannon
      shannon_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = Shannon, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_shannon <- aov(Shannon ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_shannon)
      TukeyHSD(anov_shannon, which = "Type")
    ### Chao1
      chao_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = Chao1, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_chao <- aov(Chao1 ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_chao)
      TukeyHSD(anov_chao, which = "Type")
    ### PD
      PD_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = PD, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_PD <- aov(PD ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_PD)
      TukeyHSD(anov_PD, which = "Type")
    ### MPD
      MPD_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = MPD, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_MPD <- aov(MPD ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_MPD)
      TukeyHSD(anov_MPD, which = "Type")
    ### MPD
      MNTD_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = MNTD, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_MNTD <- aov(MNTD ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_MNTD)
      TukeyHSD(anov_MNTD, which = "Type")
    ### Make plot object
      p_alpha_nofilt <- (chao_all | shannon_all | PD_all) + plot_layout(guides = 'collect')
  ## ^Plant samples filtered at 3 percent ====
    ### Grossman
      alphadiv <- estimate_richness(grossman_physeq_filtered, measures=c("Shannon","Chao1"))
      pd_df <- pd(t(as.data.frame(otu_table(grossman_physeq_filtered))), tree, include.root=T)
      alphadiv$PD <- pd_df$PD
      grossman_alphadiv <- merge(alphadiv, sample_data(grossman_physeq_filtered), by=0)
      row.names(grossman_alphadiv) <- grossman_alphadiv$Row.names
      grossman_alphadiv[1] <- NULL
    ### Potato Hill
      alphadiv <- estimate_richness(potato_physeq_filtered, measures=c("Shannon","Chao1"))
      pd_df <- pd(t(as.data.frame(otu_table(potato_physeq_filtered))), tree, include.root=T)
      alphadiv$PD <- pd_df$PD
      potato_alphadiv <- merge(alphadiv, sample_data(potato_physeq_filtered), by=0)
      row.names(potato_alphadiv) <- potato_alphadiv$Row.names
      potato_alphadiv[1] <- NULL
    ### Combine the two
      all_alphadiv <- rbind(potato_alphadiv, grossman_alphadiv)
    ### Shannon
      shannon_all <- ggplot(potato_alphadiv, aes(x = Quadrat, y = Shannon, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_shannon <- aov(Shannon ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_shannon)
      TukeyHSD(anov_shannon, which = "Type")
    ### Chao1
      chao_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = Chao1, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_chao <- aov(Chao1 ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_chao)
      TukeyHSD(anov_chao, which = "Type")
    ### PD
      PD_all <- ggplot(all_alphadiv, aes(x = Quadrat, y = PD, fill = Type)) + geom_boxplot(position = position_dodge(preserve = "single")) + 
        scale_fill_manual(values = cbp1) #+ facet_grid(.~Quadrat)
      anov_PD <- aov(PD ~ Type + Quadrat + Time, data = potato_alphadiv)
      summary(anov_PD)
      TukeyHSD(anov_PD, which = "Type")
    ### Make plot object
      p_alpha_3percfilt <- (chao_all | shannon_all | PD_all) + plot_layout(guides = 'collect')  
  ## ^Plant samples filtered at 3 percent; no soil ====
    ### Grossman
      alphadiv <- estimate_richness(subset_samples(grossman_physeq_filtered, Type!="soil"), measures=c("Shannon","Chao1","Simpson"))
      pd_df <- pd(t(as.data.frame(otu_table(subset_samples(grossman_physeq_filtered, Type!="soil")))), tree, include.root=T)
      alphadiv$PD <- pd_df$PD
      mpd_df <- mpd(t(as.data.frame(otu_table(subset_samples(grossman_physeq_filtered, Type!="soil")))), phydist, abundance.weighted = T)
      alphadiv$MPD <- mpd_df
      grossman_alphadiv <- merge(alphadiv, sample_data(subset_samples(grossman_physeq_filtered, Type!="soil")), by=0)
      row.names(grossman_alphadiv) <- grossman_alphadiv$Row.names
      grossman_alphadiv[1] <- NULL
      grossman_alphadiv$Site <- rep("Grossman", nrow(grossman_alphadiv))
    ### Potato Hill
      alphadiv <- estimate_richness(subset_samples(potato_physeq_filtered, Type!="soil"), measures=c("Shannon","Chao1", "Simpson"))
      pd_df <- pd(t(as.data.frame(otu_table(subset_samples(potato_physeq_filtered, Type!="soil")))), tree, include.root=T)
      alphadiv$PD <- pd_df$PD
      mpd_df <- mpd(t(as.data.frame(otu_table(subset_samples(potato_physeq_filtered, Type!="soil")))), phydist, abundance.weighted = T)
      alphadiv$MPD <- mpd_df
      potato_alphadiv <- merge(alphadiv, sample_data(subset_samples(potato_physeq_filtered, Type!="soil")), by=0)
      row.names(potato_alphadiv) <- potato_alphadiv$Row.names
      potato_alphadiv[1] <- NULL
      potato_alphadiv$Site <- rep("Potato", nrow(potato_alphadiv))
    ### Combine the two
      all_alphadiv <- rbind(potato_alphadiv, grossman_alphadiv)
    ### Shannon
      shannon_all <- ggplot(all_alphadiv, aes(x = Type, y = Shannon)) + #ggplot(all_alphadiv, aes(x = QuadratxType, y = Shannon, label = Quadrat)) + 
        geom_violin(aes(fill=Type), scale="width") + #, position="dodge") + 
        geom_boxplot(width=0.1) + #, position="dodge") + 
        #geom_text(position = position_dodge(width = 1), aes(x=Quadrat, y=1, angle=90)) + 
        scale_fill_manual(values = cbp1) + scale_color_manual(values = cbp1)  #+ facet_grid(.~Quadrat)
      anov_shannon <- aov(Shannon ~ Type * Quadrat * Time, data = all_alphadiv)
      summary(anov_shannon)
      TukeyHSD(anov_shannon, which = "Type")
    ### Chao1
      chao_all <- ggplot(all_alphadiv, aes(x = Type, y = Chao1)) + 
        geom_violin(aes(fill=Type), scale="width") + 
        geom_boxplot(width=0.1) + 
        scale_fill_manual(values = cbp1) + scale_color_manual(values = cbp1)  #+ facet_grid(.~Quadrat)
      anov_chao <- aov(Chao1 ~ Type * Quadrat * Time, data = all_alphadiv)
      summary(anov_chao)
      TukeyHSD(anov_chao, which = "Type")
    ### Simpson
      simpson_all <- ggplot(all_alphadiv, aes(x = Type, y = Simpson)) + 
        geom_violin(aes(fill=Type), scale="width") + 
        geom_boxplot(width=0.1) + 
        scale_fill_manual(values = cbp1) + scale_color_manual(values = cbp1)  #+ facet_grid(.~Quadrat)
      anov_simpson <- aov(Chao1 ~ Type * Quadrat * Time, data = all_alphadiv)
      summary(anov_simpson)
      TukeyHSD(anov_simpson, which = "Type")
    ### PD
      PD_all <- ggplot(all_alphadiv, aes(x = Type, y = PD)) + 
        geom_violin(aes(fill=Type), scale="width") + 
        geom_boxplot(width=0.1) + 
        scale_fill_manual(values = cbp1) + scale_color_manual(values = cbp1)  #+ facet_grid(.~Quadrat)
      anov_PD <- aov(PD ~ Type * Quadrat * Time, data = potato_alphadiv)
      summary(anov_PD)
      TukeyHSD(anov_PD, which = "Type")
    ### MPD
      MPD_all <- ggplot(all_alphadiv, aes(x = Type, y = MPD)) + 
        geom_violin(aes(fill=Type), scale="width") + 
        geom_boxplot(width=0.1) + 
        scale_fill_manual(values = cbp1) + scale_color_manual(values = cbp1)  #+ facet_grid(.~Quadrat)
      anov_MPD <- aov(MPD ~ Type + Quadrat + Time, data = all_alphadiv)
      summary(anov_MPD)
      TukeyHSD(anov_MPD, which = "Type")
      
    ### Make plot object
     p_alpha_3percfilt_nosoil <- (chao_all | simpson_all | shannon_all | PD_all | MPD_all) + plot_layout(guides = 'collect')
  ## ^Plot ====
    p_alpha_3percfilt_nosoil
    ggsave("alpha_diversity_filt_nosoil.pdf", device = "pdf", width = 12, height = 5)
    
    p_distr / p_alpha_3percfilt_nosoil
    ggsave("Read_ASV_distribution_alpha_diversity.pdf", device = "pdf", width = 12, height = 6)
    
  ## Alpha Diversity by species
    alphadiv_by_cat <- function(phyloseq_obj, category){
      #phyloseq_obj <- potato_physeq_plantsonly_filtered
      #category <- "Type"
      sampleda <- sample_data(phyloseq_obj)[,match(category, colnames(sample_data(phyloseq_obj))),drop=FALSE]
      asv_meta <- merge(t(otu_table(phyloseq_obj)), sampleda, by=0)
      asv_meta <- tibble::column_to_rownames(asv_meta, var="Row.names")
      asv_summed_by_type <- data.frame(plyr::ddply(asv_meta, category, plyr::numcolwise(sum)), 
                             check.names=FALSE, stringsAsFactors=FALSE)   
      asv_summed_by_type <- tibble::column_to_rownames(asv_summed_by_type, var=category)
      asv_summed_by_type_table <- otu_table(asv_summed_by_type, taxa_are_rows = F)
      
      alphadiv <- estimate_richness(asv_summed_by_type_table, measures=c("Shannon","Chao1", "Simpson"))
      pd_df <- pd(asv_summed_by_type, tree, include.root=T)
      alphadiv$PD <- pd_df$PD
      mpd_df <- mpd(asv_summed_by_type, phydist, abundance.weighted = T)
      alphadiv$MPD <- mpd_df
      return(alphadiv)
    }
    alphadiv_by_cat(grossman_physeq_plantsonly_filtered,"Type")  
    alphadiv_by_cat(potato_physeq_plantsonly_filtered,"QuadratxType")  
    alphadiv_by_cat(time_series_physeq_filtered,"Type")

# Top ASV ------------------------------------------------------------
    ASV_melt <- reshape2::melt(otu_table(time_series_physeq_filtered_tranformed),value.name="abundance",varnames=c("ASV","Sample"))
    top_ASV_count <- ASV_melt %>% group_by(Sample) %>% summarize(topASV=max(abundance))
    meta <- as_tibble(sample_data(time_series_physeq_filtered_tranformed), rownames="Sample")
    top_ASV_count_meta <- left_join(top_ASV_count, meta, by = "Sample")
    count_50 <- top_ASV_count_meta %>% filter(topASV>0.5) %>% group_by(Type) %>% summarize(count=n())
    top_ASV_count_meta %>% group_by(Type) %>% summarize(count=n())
    
    ggplot(top_ASV_count_meta) + geom_histogram(aes(x=topASV, fill=Type),binwidth = 0.075) + 
      geom_vline(data=top_ASV_count_meta, aes(xintercept=0.5), linetype="dashed") + 
      labs(title = "", y="Number of sample", x="Relative abundance of the most abundant ASV") + 
      scale_fill_manual(values = cbp1) + facet_wrap(~Type,ncol=1) 
    ggsave("top_ASV_distri.pdf", device = "pdf", width = 5, height = 8)

    ggplot(top_ASV_count_meta, aes(x = Type, y = topASV)) + 
      geom_violin(aes(fill=Type), scale="width") + geom_boxplot(width=0.1) +
      labs(title = "", y="Relative abundance") + 
      scale_colour_manual(values = cbp1) + scale_fill_manual(values = cbp1)
    anov_topASV <- aov(topASV ~ Type + Quadrat + Time, data = top_ASV_count_meta)
    summary(anov_topASV)
    TukeyHSD(anov_topASV, which = "Type")
 
    x_coor <- vector()
    y_coor <- vector()
    sample_name_list <- vector()
    phyloseq_obj <- time_series_physeq_transformed
    for ( sample in 1:ncol(as.data.frame(otu_table(phyloseq_obj))) ) {
      sample_name <- colnames(as.data.frame(otu_table(phyloseq_obj)))[sample]
      list <- as.data.frame(otu_table(phyloseq_obj))[,sample]
      ASV_freq_list <- sort(list, decreasing = F)
      ASV_freq_list <- ASV_freq_list[ ASV_freq_list > 0 ]
      counter = 1
      for (freq in ASV_freq_list) { 
        y_coor <- append(y_coor, freq)
        sample_name_list <- append(sample_name_list, sample_name)
        if (counter==1) {
          x_coor <- append(x_coor, 0) 
          x_accumulator <- freq }
        else {
          x_coor <- append(x_coor, x_accumulator)
          x_accumulator <- x_accumulator + freq }
        #print(counter)
        if (counter==length(ASV_freq_list)) {
          y_coor <- append(y_coor, freq)
          x_coor <- append(x_coor, x_accumulator) 
          sample_name_list <- append(sample_name_list, sample_name)}
        counter = counter + 1
        }
    }
    table <- tibble(a=x_coor, b=y_coor, Sample=sample_name_list)
    ASV_accum_meta <- left_join(table, meta, by = "Sample")
    ggplot(ASV_accum_meta) + geom_step(aes(x=a, y=b, group=Sample, color=Type)) +
      labs(title = "", x="Proportion of the community", y="Relative abundance of each ASV") + 
      scale_color_manual(values = cbp1) + facet_wrap(~Type,ncol=2) 
    ggsave("all_ASV_distri.pdf", device = "pdf", width = 12, height = 8)
    
### BSL_SHH #### 
  # Read ASV table
  ASV_BSLSHH_file <- as.matrix(read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/BSL_SHH/ASV_Table_SHHBSL.csv", header=TRUE, sep = ",", row.names = 1))
  # Read sample metadata
  sample_BSLSHH_meta_file <- read.table("/Users/fay-weili/Box/hornwort_amplicon/dada2/BSL_SHH/SHHBSLMeta.csv", header=TRUE, sep = ",", row.names = 1)
  # Make phyloseq object
  ASV_table_BSLSHH <- otu_table(ASV_BSLSHH_file, taxa_are_rows = TRUE)
  sample_meta_BSLSHH <- sample_data(sample_BSLSHH_meta_file)
  BSLSHH_physeq = phyloseq(ASV_table_BSLSHH, sample_meta_BSLSHH) #, ASV_tree_file)
  BSLSHH_physeq_transformed  = transform_sample_counts(BSLSHH_physeq, function(x) x / sum(x) )
  BSLSHH_physeq_transformed_filtered = filter_taxa(BSLSHH_physeq_transformed, function(x) sum(x) > 0, TRUE)  
  ## ALL ordination plot
  BSLSHH_physeq_transformed.ord <- ordinate(BSLSHH_physeq_transformed, "PCoA", "bray")
  plot_ordination(BSLSHH_physeq_transformed, BSLSHH_physeq_transformed.ord, type="samples", color="Site") #+ geom_polygon(aes(fill=Quadrat))
  ## BSL
  BSL_physeq_transformed <- subset_samples(BSLSHH_physeq_transformed, Site=="BSL")
  BSL_physeq_transformed.ord <- ordinate(BSL_physeq_transformed, "PCoA", "bray")
  BSL_pcoa <- plot_ordination(BSL_physeq_transformed, BSL_physeq_transformed.ord, type="samples", color="Type", title="BSL") #+ geom_polygon(aes(fill=Quadrat))
  ## SHH
  SHH_physeq_transformed <- subset_samples(BSLSHH_physeq_transformed, Site=="SHH")
  SHH_physeq_transformed.ord <- ordinate(SHH_physeq_transformed, "PCoA", "bray")
  SHH_pcoa <- plot_ordination(SHH_physeq_transformed, SHH_physeq_transformed.ord, type="samples", color="Type", title="SHH") #+ geom_polygon(aes(fill=Quadrat))
  
  BSL_pcoa + SHH_pcoa

  

  
  
  
  
  
  
  
  
  