library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(topGO)
library(pheatmap)
library(VennDiagram)
library(grDevices)
library(Pigengene)
library(GO.db)
library(reshape2)
library(egg)
library(eulerr)
library(openxlsx)
library(SuperExactTest)
library(cowplot)
library(ggpubr)
source("dactylorhiza_functions.R")


############################################################
#       make supplementary table of read count stats       #
############################################################

# the file here is generated throughb a simple script found in /scratch/botany/katie/orchid/collate_mapping_stats.sh

star_summary_stats<-read.table("rna_seq_summary_statistics.txt") %>% set_colnames(c("Library", "Input Reads", "Uniquely Mapped", "Percent Multimapped"))
colnames(star_summary_stats)
star_summary_stats$Library<-str_remove(star_summary_stats$Library, "Log.final.out")
star_summary_stats$`Percent Multimapped` <- str_remove(star_summary_stats$`Percent Multimapped`, "%") %>% as.numeric()
star_summary_stats$Library<-str_split_fixed(star_summary_stats$Library, "/", 3)[,3]

star_summary_stats %<>% mutate(Tissue=case_when(endsWith(star_summary_stats$Library, "r") == TRUE ~ "Root",
                                                !endsWith(star_summary_stats$Library, "r") == TRUE ~ "Leaf")) %>% 
  mutate(species=case_when(substr(Library,1,1) == "m" ~ "majalis", substr(Library,1,1) == "t" ~ "traunsteineri")) %>%
  mutate(locality=case_when(substr(Library,3,3) == "S" ~ "St Ulrich", substr(Library,3,3) == "K" ~ "Kitzbuhl")) 


star_summary_stats %>% summary()

write.xlsx(star_summary_stats, file = "SupplementaryTable1_read_mapping_summary.xlsx")


aggregate(`Uniquely Mapped` ~ locality, data = star_summary_stats, mean)


#############################################################################
#        load in the GO ID mappings to transcripts for use later on         #
#############################################################################

# for the GO term enrichment tests
mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")

#names(mp) <- names(mp) %>% 
#  str_remove("-RA") %>%
#  str_remove("-RB") %>%
#  str_remove("-RC")

#############################################################################
#       make a function fo annotating samples sing naming convention        #
#############################################################################

annotate_samples<-function(sample_tibble){
  sample_tibble %<>% 
    mutate(species=case_when(substr(sample_id,1,1) == "m" ~ "majalis", substr(sample_id,1,1) == "t" ~ "traunsteineri")) %>%
    mutate(locality=case_when(substr(sample_id,3,3) == "S" ~ "St Ulrich", substr(sample_id,3,3) == "K" ~ "Kitzbuhl")) %>%
    mutate(environment=environment<-case_when(substr(sample_id,2,2) == "M" ~ "majalis", substr(sample_id,2,2) == "T" ~ "traunsteineri")) %>%
    mutate(replicate=substr(sample_id,4,4)) %>%
    mutate(treatment=case_when(species == environment ~ "native", species != environment ~ "transplant"))
  sample_tibble %<>% column_to_rownames("sample_id")
  return(sample_tibble)
}


##################################################################################################
#      make a function to read in featurecounts and remove unecessary strings in colnames        #
##################################################################################################

read_in_featurecounts<-function(input_path, strings_to_remove){
  # read in each featurecounts file and join them together
  df <- list.files(path=input_path, full.names = TRUE) %>% 
    lapply(read_tsv, skip=1) %>% 
    purrr::reduce(left_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
  # get and store the first standard featurecounts column names
  fc_cols<-colnames(df)[1:6]
  #get and store the egbne lengths to return for DE later
  gene_lengths<-df$Length
  # take the 7th until the last colnames of the data frame
  sample_names<-colnames(df)[7:length(colnames(df))]
  # remove any prefixes that need to be removed
  for (strings in strings_to_remove){
    sample_names<-str_replace(sample_names, strings, "")
  }
  # now set the df colnames to the shortened sample names for ease of reading
  colnames(df)<-c(fc_cols, sample_names)
  # get the relevant columns (geneid and sample counts) and set the geneid column to be rownames
  df_counts<-df %>% 
    dplyr::select(c(1,7:length(colnames(df)))) %>% 
    column_to_rownames("Geneid")
  return(list(counts=df_counts, lengths=gene_lengths))
}

#############################################################################
#          read in samples, edit out irrelevant samples from root           #
#############################################################################

leaf_samples<-read.table("Leaf_samples.txt", header=FALSE, col.names = c("sample_id"))
root_samples<-read.table("Root_samples.txt", header=FALSE, col.names = c("sample_id")) %>% 
  mutate(assay=case_when(endsWith(sample_id, "_s") == TRUE ~ "sRNA",
                         endsWith(sample_id, "_s") == FALSE ~ "RNAseq")) %>%
  mutate(tissue=case_when(substr(sample_id,5,5) == "l" ~ "leaf",
                          substr(sample_id,5,5) == "r" ~ "root")) %>%
  filter(tissue == "root" & assay == "RNAseq") %>% dplyr::select(sample_id)


#################################################################
#    specify samples to exclude from sample design and counts   #
#################################################################
samples_to_exclude<-c("mTS2", "mTS2r", 
                      "tTS4", "tTS4r",
                      "tTK1", "tTK1r",
                      "mTK1", "mTK1r",
                      "tMK4", "tMK4r",
                      "tMS4", "tMS4r")



###########################################################
#          annotate samples using above function          #
###########################################################
leaf_samples<-annotate_samples(leaf_samples)
root_samples<-annotate_samples(root_samples)

leaf_samples <- leaf_samples[!(row.names(leaf_samples) %in% samples_to_exclude),]
root_samples <- root_samples[!(row.names(root_samples) %in% samples_to_exclude),]


####################################################################################################
#  define unwanted strings and read in featurecounts, joining and removing unwanted string         #
####################################################################################################

strings_to_remove_root<-c("root_samples/", "Aligned.sortedByCoord.out.bam")
df_root<-read_in_featurecounts('dactylorhiza_root_featurecounts', strings_to_remove_root)
df_root$counts %<>% dplyr::select(-contains(samples_to_exclude))
df_counts_root<-df_root$counts
df_lengths_root<-df_root$lengths

strings_to_remove_leaf<-c("StUlrich/", "Kitzbuhel/", "Aligned.sortedByCoord.out.bam")
df_leaf<-read_in_featurecounts('dactylorhiza_leaf_featurecounts', strings_to_remove_leaf)
df_leaf$counts %<>% dplyr::select(-contains(samples_to_exclude))
df_counts_leaf<-df_leaf$counts
df_lengths_leaf<-df_leaf$lengths



###########################################################
#       get number of genes with zero reads in a sample  #
###########################################################


number_zero_leaf_genes<-apply(df_counts_leaf, 2, function(x) length(which(x == 0))) %>% data.frame() %>% rownames_to_column(var=c("name")) %>% set_colnames(c("name", "number_zero_genes"))
number_zero_root_genes<-apply(df_counts_root, 2, function(x) length(which(x == 0))) %>% data.frame() %>% rownames_to_column(var=c("name")) %>% set_colnames(c("name", "number_zero_genes"))


###############################
#        Gather PCA data      #
###############################
cpm_threshold = 1
min_count_per_sample = 3

# make a dds object from the total root samples (no subsetting)
root_dds<-specify_comparison(root_samples, df_counts_root, "1 == 1")
root_dds <- DESeqDataSetFromMatrix(countData = root_dds[["counts"]],
                                   colData = root_dds[["samples"]],
                                   design = ~ species + locality)
mcols(root_dds)$basepairs<-df_lengths_root
root_dds <- estimateSizeFactors(root_dds)
idx <- rowSums( counts(root_dds, normalized=TRUE) >= cpm_threshold ) >= min_count_per_sample
root_dds <- root_dds[idx,]

# make a dds object from the total leaf samples (no subsetting)
leaf_dds<-specify_comparison(leaf_samples, df_counts_leaf, "1 == 1")
leaf_dds <- DESeqDataSetFromMatrix(countData = leaf_dds[["counts"]],
                                   colData = leaf_dds[["samples"]],
                                   design = ~ species + locality)
mcols(leaf_dds)$basepairs<-df_lengths_leaf
leaf_dds <- estimateSizeFactors(leaf_dds)
idx <- rowSums( counts(leaf_dds, normalized=TRUE) >= cpm_threshold ) >= min_count_per_sample
leaf_dds <- leaf_dds[idx,]


# perform variance stabilizig transformation
root_vst<-varianceStabilizingTransformation(root_dds)
leaf_vst<-varianceStabilizingTransformation(leaf_dds)


# plot the PCA 
pcaRootData <-plotPCA(root_vst, intgroup=c("treatment", "species", "locality"), ntop = 2400, returnData = TRUE)
percentRootVar <- round(100 * attr(pcaRootData, "percentVar"))

pcaLeafData <-plotPCA(leaf_vst, intgroup=c("treatment", "species", "locality"), ntop = 2400, returnData = TRUE)
percentLeafVar <- round(100 * attr(pcaLeafData, "percentVar"))

pcaRootData %<>% mutate("Tissue" = "Root")
pcaLeafData %<>% mutate("Tissue" = "Leaf")

pca_plot_leaf<-ggplot(pcaLeafData, aes(PC1, PC2, color=locality, fill=locality, label=pcaLeafData$name, shape=interaction(species, treatment))) +
  geom_point(size=5, stroke = 1.5) +
  ggtitle("Leaf")

pca_plot_root<-ggplot(pcaRootData, aes(PC1, PC2, color=locality, fill=locality, label=pcaRootData$name, shape=interaction(species, treatment))) +
  geom_point(size=5, stroke = 1.5) + 
  ggtitle("Root")

pca_plot_leaf + geom_text_repel()
pca_plot_root + geom_text_repel()

# show PCA plots
pca_plot_leaf
pca_plot_root


##########################################
#        Get thresholds of zero reads    #
##########################################

number_zero_leaf_genes_PCA <- inner_join(pcaLeafData, number_zero_leaf_genes)
number_zero_leaf_genes_PCA %<>% mutate(colour_by=ifelse(number_zero_leaf_genes_PCA$number_zero_genes > 40000, "to exclude", "to keep"))

number_zero_root_genes_PCA <- inner_join(pcaRootData, number_zero_root_genes)
number_zero_root_genes_PCA %<>% mutate(colour_by=ifelse(number_zero_root_genes_PCA$number_zero_genes > 31500, "to exclude", "to keep"))

number_zero_genes_PCA <- rbind(number_zero_leaf_genes_PCA,
                               number_zero_root_genes_PCA)



png("dactylorhizaQC_PC1_vs_zeroCountGenes.png", width=1200, height=900)
ggplot(number_zero_genes_PCA, aes(x=PC1, y=number_zero_genes, colour=colour_by)) +
  geom_point(alpha=0.5, size = 10) + 
  facet_wrap(~ Tissue, ncol=2) +
  theme(text = element_text(size = 20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        legend.title=element_blank(),
        strip.text.x = element_text(size = 27), 
        #strip.text.x = element_blank(), 
        plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) +
  xlab("PC1") +
  ylab("Number of genes with zero counts")
dev.off()


png("dactylorhizaQC_PC2_vs_zeroCountGenes.png", width=1200, height=900)
ggplot(number_zero_genes_PCA, aes(x=PC2, y=number_zero_genes, colour=colour_by)) +
  geom_point(alpha=0.5, size = 10) + 
  facet_wrap(~ Tissue + treatment, ncol=2) +
  theme(text = element_text(size = 20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        legend.title=element_blank(),
        strip.text.x = element_text(size = 27), 
        #strip.text.x = element_blank(), 
        plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) +
  xlab("PC2") +
  ylab("")
dev.off()



###############################################
#       plot PCAs without excluded samples    #
###############################################

# print the excluded samples
number_zero_leaf_genes_PCA %>% filter(colour_by == "to exclude")
number_zero_root_genes_PCA %>% filter(colour_by == "to exclude")

# excliude to exlude samples and plot PCAs
number_zero_leaf_genes_PCA_keep <- number_zero_leaf_genes_PCA %>% filter(colour_by == "to keep")
number_zero_root_genes_PCA_keep <- number_zero_root_genes_PCA %>% filter(colour_by == "to keep")

pcaData <- rbind(number_zero_leaf_genes_PCA_keep, number_zero_root_genes_PCA_keep)

png("dactylorhizaQC_PCA_without_excluded_samples.png", width=1000, height=400)
ggplot(pcaData, aes(PC1, PC2, color=locality, fill=locality, shape=interaction(species, treatment))) +
  geom_point(size=5, stroke = 1.5) +
  #xlab(paste0("PC1: ",percentRootVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentRootVar[2],"% variance")) + 
  #coord_fixed() + 
  theme(legend.title = element_blank(),
        text = element_text(size = 20), 
        legend.text=element_text(size=20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        plot.margin = unit(c(0,0,0,0), "cm"),
        strip.text.x = element_text(size = 20),) +
  scale_shape_manual(labels=c("D. majalis native", "D. traunstaineri native", "D. majalis transplant", "D. traunstaineri transplant"), values = c(16, 15, 10, 7)) +
  facet_wrap(~ Tissue, ncol=2)
dev.off()


exclude_samples <- c("mTK1", "tMK4", "tMS4")
pcaData <- rbind(pcaRootData, pcaLeafData)
pcaLeafData_filt <- pcaLeafData %>% filter(!(name %in% exclude_samples))

ggplot(pcaLeafData_filt, aes(PC1, PC2, color=locality, fill=locality, shape=interaction(species, treatment))) +
  geom_point(size=5, stroke = 1.5) +
  scale_shape_manual(labels=c("D. majalis native", "D. traunstaineri native", "D. majalis transplant", "D. traunstaineri transplant"), values = c(16, 15, 10, 7)) +
  facet_wrap(~ Tissue, ncol=2)
dev.off()






