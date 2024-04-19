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


samples_to_exclude_leaf<-c("mTS2", 
                      "tTS4",
                      "tTK1",
                      "mTK1",
                      "tMK4",
                      "tMS4",
                      "tTK3")

samples_to_exclude_root<-c("mTS2r", 
                           "tTS4r",
                           "tTK1r",
                           "mTK2r", 
                           "tMK4r")



###########################################################
#          annotate samples using above function          #
###########################################################
leaf_samples<-annotate_samples(leaf_samples)
root_samples<-annotate_samples(root_samples)

leaf_samples <- leaf_samples[!(row.names(leaf_samples) %in% samples_to_exclude_leaf),]
root_samples <- root_samples[!(row.names(root_samples) %in% samples_to_exclude_root),]


####################################################################################################
#  define unwanted strings and read in featurecounts, joining and removing unwanted string         #
####################################################################################################

strings_to_remove_root<-c("root_samples/", "Aligned.sortedByCoord.out.bam")
df_root<-read_in_featurecounts('dactylorhiza_root_featurecounts', strings_to_remove_root)
df_root$counts %<>% dplyr::select(-contains(samples_to_exclude_root))
df_counts_root<-df_root$counts
df_lengths_root<-df_root$lengths

strings_to_remove_leaf<-c("StUlrich/", "Kitzbuhel/", "Aligned.sortedByCoord.out.bam")
df_leaf<-read_in_featurecounts('dactylorhiza_leaf_featurecounts', strings_to_remove_leaf)
df_leaf$counts %<>% dplyr::select(-contains(samples_to_exclude_leaf))
df_counts_leaf<-df_leaf$counts
df_lengths_leaf<-df_leaf$lengths


######################################
#      Figure 2 Draw PCA plots       #
######################################

# make a dds object from the total root samples (no subsetting)
root_kitzbuhel_dds<-specify_comparison(root_samples, df_counts_root, "locality == 'Kitzbuhl'")
root_stulrich_dds<-specify_comparison(root_samples, df_counts_root, "locality == 'St Ulrich'")


root_kitzbuhel_dds <- DESeqDataSetFromMatrix(countData = root_kitzbuhel_dds[["counts"]],
                                   colData = root_kitzbuhel_dds[["samples"]],
                                   design = ~ species + treatment)

root_stulrich_dds <- DESeqDataSetFromMatrix(countData = root_stulrich_dds[["counts"]],
                                             colData = root_stulrich_dds[["samples"]],
                                             design = ~ species + treatment)



# make a dds object from the total leaf samples (no subsetting)
leaf_kitzbuhel_dds<-specify_comparison(leaf_samples, df_counts_leaf, "locality == 'Kitzbuhl'")
leaf_stulrich_dds<-specify_comparison(leaf_samples, df_counts_leaf, "locality == 'St Ulrich'")

leaf_kitzbuhel_dds <- DESeqDataSetFromMatrix(countData = leaf_kitzbuhel_dds[["counts"]],
                                   colData = leaf_kitzbuhel_dds[["samples"]],
                                   design = ~ species + treatment)

leaf_stulrich_dds <- DESeqDataSetFromMatrix(countData = leaf_stulrich_dds[["counts"]],
                                   colData = leaf_stulrich_dds[["samples"]],
                                   design = ~ species + treatment)


# perform variance stabilizig transformation
root_kitzbuhel_vst<-varianceStabilizingTransformation(root_kitzbuhel_dds)
root_stulrich_vst<-varianceStabilizingTransformation(root_stulrich_dds)
leaf_kitzbuhel_vst<-varianceStabilizingTransformation(leaf_kitzbuhel_dds)
leaf_stulrich_vst<-varianceStabilizingTransformation(leaf_stulrich_dds)


# plot the PCA 
pcaRootKitzbuhelData <-plotPCA(root_kitzbuhel_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentRootKitzbuhelVar <- round(100 * attr(pcaRootKitzbuhelData, "percentVar"))
pcaRootStUlrichData <-plotPCA(root_stulrich_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentRootStUlrichVar <- round(100 * attr(pcaRootStUlrichData, "percentVar"))

pcaLeafKitzbuhelData <-plotPCA(leaf_kitzbuhel_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentLeafKitzbuhelVar <- round(100 * attr(pcaLeafKitzbuhelData, "percentVar"))
pcaLeafStUlrichData <-plotPCA(leaf_stulrich_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentLeafStUlrichVar <- round(100 * attr(pcaLeafStUlrichData, "percentVar"))

pcaRootKitzbuhelData %<>% mutate("Tissue" = "Root", "Locality" = "Kitzbuhel")
pcaRootStUlrichData %<>% mutate("Tissue" = "Root", "Locality" = "St Ulrich")

pcaLeafKitzbuhelData %<>% mutate("Tissue" = "Leaf", "Locality" = "Kitzbuhel")
pcaLeafStUlrichData %<>% mutate("Tissue" = "Leaf", "Locality" = "St Ulrich")

pcaData <- rbind(pcaRootKitzbuhelData, pcaRootStUlrichData, pcaLeafKitzbuhelData, pcaLeafStUlrichData)

pca_plot<-ggplot(pcaData, aes(PC1, PC2, color=Locality, fill=Locality, shape=interaction(species, treatment))) +
  geom_point(size=15, stroke = 1.5) +
  #xlab(paste0("PC1: ",percentRootVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentRootVar[2],"% variance")) + 
  #coord_fixed() + 
  theme(legend.title = element_blank(),
        text = element_text(size = 40), 
        legend.text=element_text(size=35),
        axis.text=element_text(size=45),
        axis.title=element_text(size=45),
        plot.margin = unit(c(1,1,1,1), "cm"),
        strip.text.x = element_text(size = 45)) +
  scale_shape_manual(labels=c("D. majalis native", "D. traunstaineri native", "D. majalis transplant", "D. traunstaineri transplant"), values = c(16, 15, 10, 7)) +
  facet_wrap(~ Tissue + Locality, ncol=2)

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2.png", width = 2000, height = 1600)
pca_plot
dev.off()

percentRootKitzbuhelVar[1]
percentRootKitzbuhelVar[2]
percentRootStUlrichVar[1]
percentRootStUlrichVar[2]
percentLeafKitzbuhelVar[1]
percentLeafKitzbuhelVar[2]
percentLeafStUlrichVar[1]
percentLeafStUlrichVar[2]

###############################################################
#      Count how many genes are expressed total per tissue    #
###############################################################
leaf_kitzbuhel_names<-leaf_samples %>% filter(locality == "Kitzbuhl") %>% rownames()
leaf_stulrich_names<-leaf_samples %>% filter(locality == "St Ulrich") %>% rownames()

root_kitzbuhel_names<-root_samples %>% filter(locality == "Kitzbuhl") %>% rownames()
root_stulrich_names<-root_samples %>% filter(locality == "St Ulrich") %>% rownames()


mcols(leaf_dds)$basepairs<-df_leaf$lengths
leaf_dds <- estimateSizeFactors(leaf_dds)
length(which(rowSums(counts(leaf_dds, normalized=TRUE)[,leaf_kitzbuhel_names] >= 1 ) >= 3))
length(which(rowSums(counts(leaf_dds, normalized=TRUE)[,leaf_stulrich_names] >= 1 ) >= 3))


mcols(root_dds)$basepairs<-df_root$lengths
root_dds <- estimateSizeFactors(root_dds)
length(which(rowSums( counts(root_dds, normalized=TRUE)[,root_kitzbuhel_names] >= 1 ) >= 3))
length(which(rowSums( counts(root_dds, normalized=TRUE)[,root_stulrich_names] >= 1 ) >= 3))

#############################################################
#      Constitutively differentially expressed genes        #
#############################################################

# these are genes which are differentially expressed between majalis and traunsteineri regardless of environemnt
# we compare traunsteineri and majalis in the traunsteineri environmernt
# the we compare traunsteineri and majalis in the majalis environmernt
# the genes in common between these two groups are genes which are differentially expressed irrespective of the environment
# we do this for Kitzbuhl and St Ulrich seprarately

# format of comparison variable names:
# speciesA_speciesB_tissue_environment_locality

####################
#     KITZBUHL     #
#################### 

# root 
traunsteineri_majalis_root_M_kitzbuhl<-specify_comparison(root_samples, df_counts_root, 'environment == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("species", df_lengths_root)
traunsteineri_majalis_root_T_kitzbuhl<-specify_comparison(root_samples, df_counts_root, 'environment == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("species", df_lengths_root)


# leaf 
traunsteineri_majalis_leaf_M_kitzbuhl<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("species", df_lengths_leaf)
traunsteineri_majalis_leaf_T_kitzbuhl<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("species", df_lengths_leaf)


####################
#    ST ULRICH     #
#################### 

# root 
traunsteineri_majalis_root_M_stulrich<-specify_comparison(root_samples, df_counts_root, 'environment == "majalis" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_root)
traunsteineri_majalis_root_T_stulrich<-specify_comparison(root_samples, df_counts_root, 'environment == "traunsteineri" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_root)

# leaf 
traunsteineri_majalis_leaf_M_stulrich<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "majalis" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_leaf)
traunsteineri_majalis_leaf_T_stulrich<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "traunsteineri" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_leaf)

##############################################
#    Figure 3 Plotting constitutively DEGs   #
##############################################

#pvalue_threshold <- 0.05
#logfc_threshold<-1.5

tmM_stu_root<-traunsteineri_majalis_root_M_stulrich$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")


tmT_stu_root<-traunsteineri_majalis_root_T_stulrich$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

tmM_ktz_root<-traunsteineri_majalis_root_M_kitzbuhl$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

tmT_ktz_root<-traunsteineri_majalis_root_T_kitzbuhl$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

tmM_stu_leaf<-traunsteineri_majalis_leaf_M_stulrich$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

tmT_stu_leaf<-traunsteineri_majalis_leaf_T_stulrich$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

tmM_ktz_leaf<-traunsteineri_majalis_leaf_M_kitzbuhl$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

tmT_ktz_leaf<-traunsteineri_majalis_leaf_T_kitzbuhl$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id")

colname_string<-c("gene_id", "majalis_env", "majalis_env_padj", "traunst_env", "traunst_env_padj", "locality", "tissue")

stu_root<-inner_join(tmM_stu_root, 
                     tmT_stu_root, 
                     by="gene_id") %>% mutate(locality="St. Ulrich", tissue="Root") %>%
  set_colnames(colname_string)

ktz_root<-inner_join(tmM_ktz_root, 
                     tmT_ktz_root, 
                     by="gene_id") %>% mutate(locality="Kitzbuhl", tissue="Root") %>%
  set_colnames(colname_string)

stu_leaf<-inner_join(tmM_stu_leaf, 
                     tmT_stu_leaf, 
                     by="gene_id") %>% mutate(locality="St. Ulrich", tissue="Leaf") %>%
  set_colnames(colname_string)

ktz_leaf<-inner_join(tmM_ktz_leaf, 
                     tmT_ktz_leaf, 
                     by="gene_id") %>% mutate(locality="Kitzbuhl", tissue="Leaf") %>%
  set_colnames(colname_string)

all_bound<-rbind(stu_root,
      ktz_root,
      stu_leaf,
      ktz_leaf)



all_bound %<>% mutate(status=case_when(majalis_env_padj < 0.05 &
                                         traunst_env_padj < 0.05 ~ "Constitutively DE",
                                         majalis_env_padj < 0.05 &
                                         traunst_env_padj > 0.05 ~ "DE in D. majalis environment only",
                                         majalis_env_padj > 0.05 &
                                         traunst_env_padj < 0.05 ~ "DE in D. traunsteineri environment only"))


# designate all unlabelled status values Not Significant         
all_bound %<>% replace_na(list(status = "Not significant"))


# arrange the data so that the DE points are plotted last and on top of the grey non significant points
all_bound %<>% arrange(factor(status, levels = c("Not significant", 
                                                 "DE in D. majalis environment only", 
                                                 "DE in D. traunsteineri environment only", 
                                                 "Constitutively DE")))



#################################
#        Plot  Constitutive     #
#################################

bar_dataset<-all_bound %>% dplyr::select(status, locality, tissue) %>% filter(status != "Not significant") %>% melt()
bar <- ggplot(bar_dataset, aes(x=locality, fill=status)) +
  geom_bar(position="stack") + facet_wrap(~ tissue, ncol=2) +
  scale_fill_manual(values = c("DE in D. majalis environment only" = "gold", 
                               "DE in D. traunsteineri environment only" = "deeppink",
                               "Constitutively DE" = "dodgerblue")) +
  ylab("Number of Genes") +
  #xlab("Locality") +
  theme(text = element_text(size = 20), 
        legend.position = "none",
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 27),
        axis.text.x = element_text(angle = 45, size = 20, vjust = 1, hjust=1))

scatt <-ggplot(all_bound, aes(x=majalis_env, y=traunst_env, colour=status)) +
  geom_point(alpha=0.5, size = 4) + 
  ylim(-10, 10) +
  xlim(-10, 10) +
  ylab("Fold change in traunsteineri environment") +
  xlab("Fold change in majalis environment") +
  scale_color_manual(values = c("DE in D. majalis environment only" = "gold", 
                                "DE in D. traunsteineri environment only" = "deeppink",
                                "Not significant" = "grey",
                                "Constitutively DE" = "dodgerblue")) + 
  facet_wrap(~ tissue + locality, ncol=2) +
  theme(text = element_text(size = 20), 
        legend.text=element_text(size=20),
        #legend.title=element_text(size=23),
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        legend.title=element_blank(),
        strip.text.x = element_text(size = 27), 
        plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) +
  guides(colour = guide_legend(override.aes = list(size=10)))

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure3.png", width = 1400, height = 580)
egg::ggarrange(bar, scatt, ncol=2, widths = c(0.5, 1), labels = c('A', 'B'), label.args = list(gp = grid::gpar(font = 2, cex = 2.6)))
dev.off()

# how many DE genes are there in each category?
all_bound %>% dplyr::select("status") %>% pull() %>% table()

# how many genes are upregulated in traunsteineri relative to majalis?
all_bound %>% filter(majalis_env > 2 & traunst_env > 2 & majalis_env_padj < 0.05 & traunst_env_padj < 0.05)

# how many genes are upregulated in majalis relative to traunsteineri?
all_bound %>% filter(majalis_env > 2 & traunst_env > 2 & majalis_env_padj < 0.05 & traunst_env_padj < 0.05)



###################################################
#        Plot GO term enrichment Constitutive     #
###################################################


de_majalis_only_root<-all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Root") %>% dplyr::select(gene_id) %>% pull()
de_traunsteineri_only_root<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Root") %>% dplyr::select(gene_id) %>% pull()
de_majalis_only_leaf<-all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Leaf") %>% dplyr::select(gene_id) %>% pull()
de_traunsteineri_only_leaf<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Leaf") %>% dplyr::select(gene_id) %>% pull()
de_constitutive_leaf <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Leaf") %>% dplyr::select(gene_id) %>% pull()
de_constitutive_root <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Root") %>% dplyr::select(gene_id) %>% pull()

de_majalis_only_root_GO<-get_enriched_terms(de_majalis_only_root, mp)
de_traunsteineri_only_root_GO<-get_enriched_terms(de_traunsteineri_only_root, mp)
de_majalis_only_leaf_GO<-get_enriched_terms(de_majalis_only_leaf, mp)
de_traunsteineri_only_leaf_GO<-get_enriched_terms(de_traunsteineri_only_leaf, mp)
de_constitutive_leaf_GO<-get_enriched_terms(de_constitutive_leaf, mp)
de_constitutive_root_GO<-get_enriched_terms(de_constitutive_root, mp)

de_majalis_only_root_GO_df <- de_majalis_only_root_GO %>% prepare_go_df() %>% mutate(Environment="D. majalis Environment Only", tissue="Root", placeholder="placeholder")
de_traunsteineri_only_root_GO_df <- de_traunsteineri_only_root_GO %>% prepare_go_df()  %>% mutate(Environment="D. traunsteineri Environment Only", tissue="Root", placeholder="placeholder")
de_majalis_only_leaf_GO_df <- de_majalis_only_leaf_GO %>% prepare_go_df()  %>% mutate(Environment="D. majalis Environment Only", tissue="Leaf", placeholder="placeholder")
de_traunsteineri_only_leaf_GO_df <- de_traunsteineri_only_leaf_GO %>% prepare_go_df()  %>% mutate(Environment="D. traunsteineri Environment Only", tissue="Leaf", placeholder="placeholder")
de_constitutive_leaf_GO_df <- de_constitutive_leaf_GO %>% prepare_go_df() %>% mutate(Environment="Constitutively DE", tissue="Leaf", placeholder="placeholder")
de_constitutive_root_GO_df<- de_constitutive_root_GO %>% prepare_go_df() %>% mutate(Environment="Constitutively DE", tissue="Root", placeholder="placeholder")

root_go_bound <- rbind(de_majalis_only_root_GO_df, de_traunsteineri_only_root_GO_df, de_constitutive_root_GO_df)
root_go_bound$Term[root_go_bound$Term == "regulation of transcription from RNA polymerase II promoter in response to stress"] <- "reg of transcription from RNApolII promoter in resp. to stress"
leaf_go_bound <- rbind(de_majalis_only_leaf_GO_df, de_traunsteineri_only_leaf_GO_df, de_constitutive_leaf_GO_df)

#root_go_bound[root_go_bound$Term == "regulation of transcription from RNA polymerase II promoter in response to stress",]$Term = "reg. of transcription from RNApolII promoter in response to stress"
#leaf_go_bound[leaf_go_bound$Term == "pectin catabolic process" & leaf_go_bound$Environment == "D. traunsteineri Environment Only",]$Term = "pectin catabolic process "

root_go_bound %<>% mutate(Term = fct_reorder(Term, Environment))
# update duplicated value so it doesnt cause problems in factor redordering
leaf_go_bound[leaf_go_bound$Term == "sodium ion import across plasma membrane" & leaf_go_bound$Environment == "Constitutively DE",]$Term <- "sodium ion import across plasma membrane "
leaf_go_bound %<>% mutate(Term = fct_reorder(Term, Environment))




a<-ggplot(root_go_bound, aes(x=placeholder, y=Term, color = Environment, size=Rich.score)) + 
  geom_point() + facet_grid(scales="free", space= "free", switch = "y", cols=vars(tissue)) + 
  theme(text = element_text(size = 80), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 70),
        strip.text.x = element_text(size = 70),
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        panel.spacing=unit(1, "lines"),
        plot.margin = margin(0, 0, 0, 0, "cm")) + 
  scale_size_continuous(range = c(18, 60)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("dodgerblue", "gold", "deeppink")) +
  theme(legend.position = "none") 

b <- ggplot(leaf_go_bound, aes(x=placeholder, y=Term, color = Environment, size=Rich.score)) + 
  geom_point() + facet_grid(scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 80), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 70),
        strip.text.x = element_text(size = 70),
        legend.text=element_text(size=60),
        legend.title=element_blank(),
        legend.position = "none",
        panel.spacing=unit(1, "lines"),
        plot.margin = margin(0, 0, 0, 0, "cm")) + 
  scale_size_continuous(range = c(18, 60)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("dodgerblue", "gold", "deeppink")) +
  scale_y_discrete(position = "right") 


b_legend<-ggplot(leaf_go_bound, aes(x=placeholder, y=Term, color = Environment, size=Rich.score)) + 
  geom_point() + facet_grid(scales="free", space= "free", cols=vars(tissue)) + 
  theme(legend.margin=margin(c(0,0,0,0)),
        legend.text=element_text(size=70),
        legend.title=element_text(size=70)) + 
  scale_color_manual(values=c("dodgerblue", "gold", "deeppink")) +
  guides(color = guide_legend(override.aes = list(size = 30))) +
  scale_size_continuous(range = c(18, 60))


leg <- cowplot::get_legend(b_legend)
leg<-as_ggplot(leg)

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure4.png", height=2000, width=5500)
plot_grid(a, b, leg,
          ncol = 3,
          rel_widths = c(2, 2, 1.1),
          rel_heights = c(2, 2, 1.1))
dev.off()


######################################################################
#    Test for whether more constitutively DEGs are found by chance   #
######################################################################




rootK_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Root" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
rootK_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Root" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
rootK_DEG_M <- c(rootK_DEG_M, rootK_consititutive)
rootK_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Root" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
rootK_DEG_T <- c(rootK_DEG_T, rootK_consititutive)

rootS_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Root" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
rootS_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Root" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
rootS_DEG_M <- c(rootS_DEG_M, rootS_consititutive)
rootS_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Root" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
rootS_DEG_T <- c(rootS_DEG_T, rootS_consititutive)


leafK_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Leaf" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
leafK_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Leaf" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
leafK_DEG_M <- c(leafK_DEG_M, leafK_consititutive)
leafK_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Leaf" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
leafK_DEG_T <- c(leafK_DEG_T, leafK_consititutive)

leafS_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Leaf" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
leafS_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Leaf" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
leafS_DEG_M <- c(leafS_DEG_M, leafS_consititutive)
leafS_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Leaf" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
leafS_DEG_T <- c(leafS_DEG_T, leafS_consititutive)




total_genes_tested_root<-unique(c(ktz_root$gene_id, stu_root$gene_id)) %>% length()
total_genes_tested_leaf<-unique(c(ktz_leaf$gene_id, stu_leaf$gene_id)) %>% length()

ResultRoot=supertest(list(rootK_DEG_M=rootK_DEG_M, 
                          rootK_DEG_T=rootK_DEG_T,
                          rootS_DEG_T=rootS_DEG_T,
                          rootS_DEG_M=rootS_DEG_M), n=total_genes_tested_root) 

ResultLeaf=supertest(list(leafK_DEG_M=leafK_DEG_M, 
                          leafK_DEG_T=leafK_DEG_T, 
                          leafS_DEG_T=leafS_DEG_T, 
                          leafS_DEG_M=leafS_DEG_M), n=total_genes_tested_leaf) 


summary(ResultRoot)$Table %>% write.xlsx(file = "SupplementaryTable2_Root_Exact_Test.xlsx")
summary(ResultLeaf)$Table %>% write.xlsx(file = "SupplementaryTable3_Leaf_Exact_Test.xlsx")

#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Supplemetary_Fig1.png", width = 1000, height = 1000)
#plot(ResultRoot, Layout="landscape", degree=2:7, sort.by="size")
#dev.off()
#
#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Supplemetary_Fig2.png", width = 1000, height = 1000)
#plot(ResultLeaf, Layout="landscape", degree=2:7, sort.by="size")
#dev.off()



#############################################
#    Table 1 Constitutively DEG GO terms    #
#############################################

# get a table of GO terms mapped to their description
godb_table<-toTable(GOTERM)

# make an empty data frame to add info to
constitutive_annotation<-data.frame(go_id=c(),
                                    Term=c(),
                                    Ontology=c(),
                                    expression_pattern=c(),
                                    gene_id=c(),
                                    locality=c(),
                                    tissue=c())

# get info for constitutive genes, and add expression pattern info
all_bound_constitutive <- all_bound %>% filter(status == "Constitutively DE")
all_bound_constitutive %<>% mutate(pattern=case_when(majalis_env > 2 & traunst_env > 2 ~ "traunsteineri > majalis", 
                                                    majalis_env < -2 & traunst_env < -2 ~ "majalis > traunsteineri",
                                                    majalis_env > 2 & traunst_env < -2 ~ "opposite", 
                                                    majalis_env < -2 & traunst_env > 2 ~ "opposite")) %>% na.omit()

# loop through the data and add GO term info into a table, write to excel output
for (i in 1:nrow(all_bound_constitutive)){
  gene_name<-(all_bound_constitutive[i,]$gene_id)
  go_ids<-mp[gene_name][[1]]
  print(go_ids)
  if (length(go_ids) >=1){
    go<-godb_table[godb_table$go_id %in% go_ids,]
    go<-godb_table[godb_table$go_id %in% go_ids,]
    new_rows<-data.frame(go[,c("go_id", "Term", "Ontology")] %>% unique())
    new_rows$expression_pattern=all_bound_constitutive[i,]$pattern
    new_rows$gene_id=gene_name
    new_rows$locality=all_bound_constitutive[i,]$locality
    new_rows$tissue=all_bound_constitutive[i,]$tissue
    constitutive_annotation <- rbind(constitutive_annotation, new_rows)
  }
}

constitutive_annotation %>% filter(!(Ontology == "CC")) %>% arrange(expression_pattern, Ontology) %>% write.xlsx(file = "Table_S1_consitutive_go_terms.xlsx")



#############################################################
#             Majalis effect of transplantation             #
#############################################################

###############################
#            LEAF             #
###############################
transplant_majalis_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$lengths, cpm_threshold=1, min_count_per_sample=5)
transplant_majalis_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$lengths, cpm_threshold=1, min_count_per_sample=5)

###############################
#            ROOT             #
###############################
transplant_majalis_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$lengths)
transplant_majalis_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$lengths)


#############################################################
#         Traunsteineri effect of transplantation           #
#############################################################

###############################
#            LEAF             #
###############################
transplant_traunsteineri_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$lengths)
transplant_traunsteineri_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$lengths)

###############################
#            ROOT             #
###############################
transplant_traunsteineri_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$lengths)
transplant_traunsteineri_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$lengths)


#############################################
#    Venn Diagram of number of DE genes     #
#############################################


a<-get_significant_genes(transplant_majalis_kitzbuhl_leaf)
b<-get_significant_genes(transplant_majalis_stulrich_leaf)
c<-get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf)
d<-get_significant_genes(transplant_traunsteineri_stulrich_leaf)


library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")


myCol<-c("deeppink", "yellowgreen", "deepskyblue", "orange1")
# Chart
venn.diagram(
  x = list(a, b, c, d),
  category.names = c("mK", "mS", "tK", "tS"),
  filename = 'Figure5A_Leaf.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1080 , 
  resolution = 1000,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)


e<-get_significant_genes(transplant_majalis_kitzbuhl_root)
f<-get_significant_genes(transplant_majalis_stulrich_root)
g<-get_significant_genes(transplant_traunsteineri_kitzbuhl_root)
h<-get_significant_genes(transplant_traunsteineri_stulrich_root)

venn.diagram(
  x = list(e, f, g, h),
  category.names = c("mK", "mS", "tK", "tS"),
  filename = 'Figure5A_Root.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1080 , 
  width = 1080 , 
  resolution = 1000,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  cat.fontface = "bold",
  cat.fontfamily = "sans"
)



#####################################
#        DE gene count plots        #
#####################################

# make a function to label up or down differential expression (significant)
# mainly so you dont have to repeat the code
label_expression_direction<-function(results_object){
  results_object %>% data.frame() %>% 
    mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                   log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>%
    replace_na(list(diffexpressed = "Not significant"))
}

mKLeaf<-transplant_majalis_kitzbuhl_leaf$results %>% label_expression_direction()
mSLeaf<-transplant_majalis_stulrich_leaf$results %>% label_expression_direction()
mKRoot<-transplant_majalis_kitzbuhl_root$results %>% label_expression_direction()
mSRoot<-transplant_majalis_stulrich_root$results %>% label_expression_direction()
mSRoot %>% filter(diffexpressed == "downregulated")

tKLeaf<-transplant_traunsteineri_kitzbuhl_leaf$results %>% label_expression_direction()
tSLeaf<-transplant_traunsteineri_stulrich_leaf$results %>% label_expression_direction()
tKRoot<-transplant_traunsteineri_kitzbuhl_root$results %>% label_expression_direction()
tSRoot<-transplant_traunsteineri_stulrich_root$results %>% label_expression_direction()

species=c(rep("D. majalis", 4), rep("D traunsteineri", 4))
tissue=c(rep("Leaf", 2), rep("Root", 2), rep("Leaf", 2), rep("Root", 2))
locality=c(rep(c("Kitzbuhel", "St. Ulrich"), 4))
individual<-c("mKLeaf", "mSLeaf", "mKRoot", "mSRoot", "tKLeaf", "tSLeaf", "tKRoot", "tSRoot")
upregulated=c(mKLeaf %>% filter(diffexpressed == "upregulated") %>% nrow(), 
              mSLeaf %>% filter(diffexpressed == "upregulated") %>% nrow(),
              mKRoot %>% filter(diffexpressed == "upregulated") %>% nrow(),
              mSRoot %>% filter(diffexpressed == "upregulated") %>% nrow(),
              tKLeaf %>% filter(diffexpressed == "upregulated") %>% nrow(),
              tSLeaf %>% filter(diffexpressed == "upregulated") %>% nrow(),
              tKRoot %>% filter(diffexpressed == "upregulated") %>% nrow(),
              tSRoot %>% filter(diffexpressed == "upregulated") %>% nrow())
downregulated=c(mKLeaf %>% filter(diffexpressed == "downregulated") %>% nrow(), 
               mSLeaf %>% filter(diffexpressed == "downregulated") %>% nrow(),
               mKRoot %>% filter(diffexpressed == "downregulated") %>% nrow(),
               mSRoot %>% filter(diffexpressed == "downregulated") %>% nrow(),
               tKLeaf %>% filter(diffexpressed == "downregulated") %>% nrow(),
               tSLeaf %>% filter(diffexpressed == "downregulated") %>% nrow(),
               tKRoot %>% filter(diffexpressed == "downregulated") %>% nrow(),
               tSRoot %>% filter(diffexpressed == "downregulated") %>% nrow())

de_counts<-data.frame(species=species,
           tissue=tissue,
           locality=locality,
           individual=individual,
           upregulated=upregulated,
           downregulated=downregulated)


# make the downregulated values negative for mirrorred barplot
de_counts$downregulated <- (-de_counts$downregulated)

de_counts %<>% melt()

de_counts$Tissue<-paste(de_counts$Tissue, "Plastic")

colnames(de_counts)<-c("Species", "Tissue", "Locality", "Individual", "Direction", "Number of genes")

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure5B.png", width = 900, height = 1200)
ggplot(de_counts, aes(x=Individual, y=`Number of genes`, fill=Direction)) + 
  geom_bar(stat="identity", position="identity") +
  facet_wrap(~ Tissue + Locality, scales = "free", ncol=2) +
  theme(text = element_text(size = 32), 
        #axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.x = element_text(angle = 45, size = 20, vjust = 1, hjust=1),
        #axis.title.y=element_blank(),
        strip.text.y = element_text(size = 32),
        legend.text=element_text(size=27),
        legend.title=element_text(size=27)) +
  scale_fill_manual(values=c("brown1", "deepskyblue1")) +
  ylab("Number of DE genes") + scale_x_discrete(labels=c("mKLeaf" = "D. majalis", 
                                                         "mSLeaf" = "D. majalis", 
                                                         "mKRoot" = "D. majalis", 
                                                         "mSRoot" = "D. majalis", 
                                                         "tKLeaf" = "D. traunsteineri", 
                                                         "tSLeaf" = "D. traunsteineri", 
                                                         "tKRoot" = "D. traunsteineri", 
                                                         "tSRoot" = "D. traunsteineri")) +
  ylim(-400, 500)
dev.off()


######**** NB after this I took 5A and 5B figureds and collated them in a pptx to create final composite figure ****##########



#############################################
#            GO term enrichment             #
#############################################

#######################################
#        majalis leaf kitzbuhl        #
#######################################
transplant_majalis_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = FALSE)$up, mp) 
transplant_majalis_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = FALSE)$down, mp) 

#######################################
#        majalis leaf st ulrich        #
#######################################
transplant_majalis_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = FALSE)$up, mp) 
transplant_majalis_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = FALSE)$down, mp) 



#######################################
#        majalis root kitzbuhl        #
#######################################
transplant_majalis_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE)$down, mp) 

#######################################
#        majalis root st ulrich       #
#######################################
transplant_majalis_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE)$up, mp) 
transplant_majalis_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE)$down, mp) 

#######################################
#     traunsteineri leaf kitzbuhl     #
#######################################
transplant_traunsteineri_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE)$up, mp) 
transplant_traunsteineri_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE)$down, mp) 


#######################################
#     traunsteineri leaf st ulrich    #
#######################################
transplant_traunsteineri_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE)$up, mp) 
transplant_traunsteineri_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE)$down, mp) 

#######################################
#     traunsteineri root kitzbuhl     #
#######################################
transplant_traunsteineri_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE)$up, mp) 
transplant_traunsteineri_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE)$down, mp) 

#######################################
#     traunsteineri root st ulrich    #
#######################################
transplant_traunsteineri_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE)$up, mp) 
transplant_traunsteineri_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE)$down, mp) 

#######################################
#        Plot GO term enrichment      #
#######################################


prepare_go_df<-function(topgo_object){
  topgo_object %<>% 
    mutate(`Rich score`=Significant/Annotated) %>% 
    dplyr::select(Term, GO.ID, Significant, Annotated, classicFisher, `Rich score`) %>% 
    #filter(classicFisher < 0.05) %>% 
    filter(classicFisher < 0.05 & `Rich score` >= 0.07) %>% 
    head(5) %>% 
    data.frame()
  return(topgo_object)
}

mKLeafUp<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_up) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. majalis")
mSLeafUp<-prepare_go_df(transplant_majalis_stulrich_leaf_up) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. majalis")
tKLeafUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_up) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSLeafUp<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_up) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. traunsteineri")
mKLeafDown<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_down) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. majalis")
mSLeafDown<-prepare_go_df(transplant_majalis_stulrich_leaf_down) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. majalis")
tKLeafDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_down) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSLeafDown<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_down) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. traunsteineri")

mKRootUp<-prepare_go_df(transplant_majalis_kitzbuhl_root_up) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. majalis")
mSRootUp<-prepare_go_df(transplant_majalis_stulrich_root_up) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. majalis")
#tKRootUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up) %>% mutate(comparison="D.t.Kitz", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
tKRootUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")

tSRootUp<-prepare_go_df(transplant_traunsteineri_stulrich_root_up) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. traunsteineri")
mKRootDown<-prepare_go_df(transplant_majalis_kitzbuhl_root_down) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. majalis")
mSRootDown<-prepare_go_df(transplant_majalis_stulrich_root_down) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. majalis")
#tKRootDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down) %>% mutate(comparison="D.t.Kitz", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
tKRootDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSRootDown<-prepare_go_df(transplant_traunsteineri_stulrich_root_down) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. traunsteineri")



leaf_go_bound<-rbind(mKLeafUp,
                     mSLeafUp,
                     tKLeafUp,
                     tSLeafUp,
                     mKLeafDown,
                     mSLeafDown,
                     tKLeafDown,
                     tSLeafDown)

root_go_bound<-rbind(mKRootUp,
                    mSRootUp,
                    tKRootUp,
                    tSRootUp,
                    mKRootDown,
                    mSRootDown,
                    tKRootDown,
                    tSRootDown)




colnames(leaf_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Rich factor", "comparison", "Direction", "Locality", "Species")
colnames(root_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Rich factor", "comparison", "Direction", "Locality", "Species")

# I am adding a new column with a placeholder the same across all rows
# I will use this column as my X axis for the GO term plots
# if I use the comparison column, it staggers the plots, and it looks better if the points are stacked
leaf_go_bound_newcol <- leaf_go_bound %>% mutate(newcol="placeholder", tissue="Leaf")
root_go_bound_newcol <- root_go_bound %>% mutate(newcol="placeholder", tissue="Root")


# I had to shorten the name of one of the facet grids (D.t. Kitz due to soace
# this alphabetically rearranges the orderr so locking in order of localitoies this way
#root_go_bound_newcol$comparison <- factor(root_go_bound$comparison, levels = unique(root_go_bound$comparison))
root_go_bound_newcol[root_go_bound_newcol$Term == "positive regulation of transcription from RNA polymerase II promoter in response to heat stress",]$Term <- "+ve reg. transcription from RNApolII prmtr in response to heat stress"
root_go_bound_newcol[root_go_bound_newcol$Term == "regulation of transcription from RNA polymerase II promoter in response to stress",]$Term <- "reg. transcription from RNApolII promoter in response to stress"
leaf_go_bound_newcol[leaf_go_bound_newcol$Term == "peptidyl-diphthamide biosynthetic process from peptidyl-histidine",]$Term <- "peptidyl-diphthamide biosyn. proc. from peptidyl-histidine"

leaf_go_bound_newcol[leaf_go_bound_newcol$Term == "negative regulation of phosphoprotein phosphatase activity",]$Term <- "-ve regulation of phosphoprotein phosphatase activity"


leaf_go_bound_newcol$Term <- leaf_go_bound_newcol$Term %>% as.character() %>% make.unique(sep = " ")
root_go_bound_newcol$Term <- root_go_bound_newcol$Term %>% as.character() %>% make.unique(sep = " ")



leaf_go_bound_newcol_order <- leaf_go_bound_newcol %>%
  arrange(comparison) %>%   
  mutate(Term = factor(Term, unique(Term)),
         Direction = factor(Direction, c("Plastic Upregulated", "Plastic Downregulated")))

root_go_bound_newcol_order <- root_go_bound_newcol %>%
  arrange(comparison) %>%   
  mutate(Term = factor(Term, unique(Term)),
         Direction = factor(Direction, c("Plastic Upregulated", "Plastic Downregulated")))

a<-ggplot(leaf_go_bound_newcol_order, aes(x=newcol, y=Term, color = Species, shape=Locality, size=`Rich factor`)) + 
  geom_point() + facet_grid(rows=vars(Direction), scales="free_y", space= "free", switch = "y", cols=vars(tissue)) + 
  theme(text = element_text(size = 71), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 65),
        strip.text.x = element_text(size = 70),
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        panel.spacing=unit(1, "lines")) + 
  scale_size_continuous(range = c(13, 33)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("mediumorchid1", "mediumseagreen")) +
  theme(legend.position = "none")



# I had to shorten the name of one of the facet grids (D.t. Kitz due to soace
# this alphabetically rearranges the orderr so locking in order of localitoies this way
#root_go_bound_newcol$comparison <- factor(root_go_bound$comparison, levels = unique(root_go_bound$comparison))


#root_go_bound_newcol[root_go_bound_newcol$Term == "positive regulation of transcription from RNA polymerase II promoter in response to heat stress",]$Term = "+ve reg. of transcription from RNApolII promoter in response to heat stress"
#root_go_bound_newcol[root_go_bound_newcol$Term == "regulation of transcription from RNA polymerase II promoter in response to stress" ,]$Term = "reg. of transcription from RNApolII promoter in response to stress" 


b<-ggplot(root_go_bound_newcol_order, aes(x=newcol, y=Term, color = Species, shape=Locality, size=`Rich factor`)) + 
  geom_point() + facet_grid(rows=vars(Direction), scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 71), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 70),
        strip.text.x = element_text(size = 70),
        legend.text=element_text(size=45),
        legend.title=element_text(size=45),
        panel.spacing=unit(1, "lines"),
        legend.key.size = unit(5, 'cm'), 
        legend.key.height = unit(3, 'cm'), 
        legend.key.width = unit(3, 'cm')) + 
  scale_size_continuous(range = c(13, 33)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("mediumorchid1", "mediumseagreen")) +
  scale_y_discrete(position = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 29)),
         shape = guide_legend(override.aes = list(size = 29)))

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure6.png", height=2800, width=4200)
egg::ggarrange(a, b, ncol=2)
dev.off()






all_env_volcano <- all_env %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.0005 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.0005 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = interaction(locality, diffexpressed))) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.0005), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "red", "grey", "grey", "yellow", "pink")) +
  #ggtitle("transplant stulrich kitzbuhl root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)
































##### comapre metal ion trauns up constitutive with majalis up metal plstic up


# define GO terms which we want to pull out
# traunsteineri constitutively upregualted related to metal ion binding
metal_genes_trans_constitut_up<-c("")




####################################################################################
#     get the genes in the GO term that topGO assigned in majalis St Ulrich upreg  #
####################################################################################

mS_up<-get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$up
geneID2GO<-mp
geneSel<-mS_up
geneSel<-factor(as.integer(names(geneID2GO) %in% geneSel))
names(geneSel)<-names(geneID2GO)

# set up the topGO object with upregulated majalis root genes
sampleGOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = geneSel, 
                    nodeSize = 10,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# get the Dinc genes in term
allGO = genesInTerm(sampleGOdata)

# now get the constitutive traunsteineri UP gene annotations and try and see if there is any overlap
traunst_const_up<-constitutive_annotation %>% 
  filter(!(Ontology == "CC") & 
           expression_pattern == "T > M") %>% arrange(expression_pattern, Ontology) %>% 
  dplyr:: select(gene_id) %>% pull()

tfg<-c()
for (i in names(allGO[transplant_majalis_stulrich_root_up$GO.ID])){
  my_intersect<-(intersect(allGO[transplant_majalis_stulrich_root_up$GO.ID][[i]], traunst_const_up))
  if (!is_empty(my_intersect)){
    print(i)
    print(my_intersect)
    tfg<-c(tfg, my_intersect)
  }
}

tfg %>% unique()

tes2t<-specify_comparison(root_samples, df_counts_root, '1 == 1')
dds <- DESeqDataSetFromMatrix(countData = tes2t$counts,
                              colData = tes2t$samples,
                              design = ~ treatment)
dds<-DESeq(dds)
dds_counts<-counts(dds, normalized=TRUE)
  
#interesting<-dds_counts[tfg %>% unique(),]
interesting<-dds_counts[mS_up %>% unique(),]


intersect(tfg %>% unique(), mS_up %>% unique())

heatmap_data <- t(apply(interesting, 1, function(x) x/mean(x)))

# log normalise, remove NAs
heatmap_data <- log2(heatmap_data)
# get data range for breaks
max_data <- max(heatmap_data, na.rm = TRUE)
min_data <- -min(heatmap_data, na.rm = TRUE)
range <- min(max_data, min_data)
heatmap_data[is.infinite(heatmap_data)] <- NA
heatmap_data[is.nan(heatmap_data)] <- NA
heatmap_data %<>% data.frame() %>% drop_na()

new_column_order<-heatmap_data %>% colnames %>% sort()
heatmap_data <- heatmap_data %>% data.frame() %>% dplyr::select(new_column_order)

annotation_col<-data.frame(c(rep("majalis native kitzbuhl", 5),
             rep("majalis native st ulrich", 4),
             rep("majalis transplant kitzbuhl", 5),
             rep("majalis transplant st ulrich", 4),
             rep("traunsteineri transplant kitzbuhl", 5),
             rep("traunsteineri transplant st ulrich", 4),
             rep("traunsteineri native kitzbuhl", 5),
             rep("traunsteineri native st ulrich", 4)), row.names=colnames(heatmap_data)
             )


colnames(annotation_col)<-"treatment"


# plot heatmap
pheatmap(heatmap_data,
         breaks = seq(-range, range, length.out = 100),
         cluster_rows = TRUE, cluster_cols = FALSE,
         treeheight_row = 0, treeheight_col = 0,
         show_rownames = F, show_colnames = T, scale="none",
         annotation_col = annotation_col)



###### repeat above code for majalis up traunst up

mS_up<-get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$up
geneID2GO<-mp
geneSel<-mS_up
geneSel<-factor(as.integer(names(geneID2GO) %in% geneSel))
names(geneSel)<-names(geneID2GO)

# set up the topGO object with upregulated majalis root genes
sampleGOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = geneSel, 
                    nodeSize = 10,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)

# get the Dinc genes in term
allGO = genesInTerm(sampleGOdata)

# also can get the majalis constitutiveky up
majalis_const_up<-constitutive_annotation %>% 
  filter(!(Ontology == "CC") & 
           expression_pattern == "M > T") %>% arrange(expression_pattern, Ontology) %>% 
  dplyr:: select(gene_id) %>% pull()





transplant_majalis_kitzbuhl_leaf_up
transplant_majalis_kitzbuhl_leaf_down
transplant_majalis_stulrich_leaf_up
transplant_majalis_stulrich_leaf_down
transplant_majalis_kitzbuhl_root_up
transplant_majalis_kitzbuhl_root_down
transplant_majalis_stulrich_root_up
transplant_majalis_stulrich_root_down
transplant_traunsteineri_kitzbuhl_leaf_up
transplant_traunsteineri_kitzbuhl_leaf_down
transplant_traunsteineri_stulrich_leaf_up
transplant_traunsteineri_stulrich_leaf_down
transplant_traunsteineri_kitzbuhl_root_up
transplant_traunsteineri_kitzbuhl_root_down
transplant_traunsteineri_stulrich_root_up
transplant_traunsteineri_stulrich_root_down










##### experimenting with manual DEseq to see what the heatmap problem is, and why its so variable between reps

test<-specify_comparison(root_samples, df_counts_root, 'environment == "majalis" & locality == "St Ulrich"')

tes2t<-specify_comparison(root_samples, df_counts_root, '1 == 1')

colSums(tes2t$counts) %>% data.frame()

dds <- DESeqDataSetFromMatrix(countData = tes2t$counts,
                              colData = tes2t$samples,
                              design = ~ treatment)
test$counts

colSums(test$counts)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]


#mcols(dds)$basepairs<-gene_lengths
dds <- DESeq(dds)
res <- results(dds, contrast=c("treatment", "native", "transplant"))

################


de_genes<-dds %>% results() %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% rownames() 
gene_counts<-counts(dds, normalized=TRUE)
de_gene_counts<-gene_counts[de_genes,]

max_data <- max(de_gene_counts, na.rm = TRUE)
min_data <- -min(de_gene_counts, na.rm = TRUE)
range <- min(max_data, min_data)


draw_heatmap2(dds)











plotCounts(dds, gene=really_de_genes[1], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[2], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[3], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[4], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[5], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[6], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[7], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[8], intgroup="treatment")
plotCounts(dds, gene=really_de_genes[9], intgroup="treatment")







#################################################
#################################################
##                                             ##
##            SUPPLEMENTARY FIGURES.           ##
##                                             ##
#################################################
#################################################



#############################################
#        Volcano plots majalis        #
#############################################

transplant_majalis_kitzbuhl_leaf_volcano <- transplant_majalis_kitzbuhl_leaf$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich kitzbuhl leaf") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)

transplant_majalis_stulrich_leaf_volcano <- transplant_majalis_stulrich_leaf$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich st ulrich leaf") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)

transplant_majalis_kitzbuhl_root_volcano <- transplant_majalis_kitzbuhl_root$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich kitzbuhl root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)

transplant_majalis_stulrich_root_volcano <- transplant_majalis_stulrich_root$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich st ulrich root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)


cowplot::plot_grid(transplant_majalis_kitzbuhl_leaf_volcano,
                   transplant_majalis_stulrich_leaf_volcano,
                   transplant_majalis_kitzbuhl_root_volcano,
                   transplant_majalis_stulrich_root_volcano,
                   labels = c('A', 
                              'B', 
                              "C", 
                              "D"))

#############################################
#        Volcano plots traunsteineri        #
#############################################


transplant_traunsteineri_kitzbuhl_leaf_volcano <- transplant_traunsteineri_kitzbuhl_leaf$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich kitzbuhl root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)

transplant_traunsteineri_stulrich_leaf_volcano <- transplant_traunsteineri_stulrich_leaf$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich st ulrich root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)




transplant_traunsteineri_kitzbuhl_root_volcano <- transplant_traunsteineri_kitzbuhl_root$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich kitzbuhl root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)

transplant_traunsteineri_stulrich_root_volcano <- transplant_traunsteineri_stulrich_root$results %>% 
  data.frame() %>% 
  mutate(diffexpressed=case_when(log2FoldChange >= 2 & padj <= 0.05 ~ "upregulated",
                                 log2FoldChange <= -2 & padj <= 0.05 ~ "downregulated")) %>% 
  replace_na(list(diffexpressed = "Not significant")) %>%
  ggplot(data = ., aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
  geom_point(size = 2) + 
  geom_vline(xintercept = c(-2, 2), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00")) +
  #ggtitle("transplant stulrich st ulrich root") +
  #theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-5, 5) + 
  ylim(0, 15)


cowplot::plot_grid(transplant_traunsteineri_kitzbuhl_leaf_volcano,
                   transplant_traunsteineri_stulrich_leaf_volcano,
                   transplant_traunsteineri_kitzbuhl_root_volcano,
                   transplant_traunsteineri_stulrich_root_volcano,
                   labels = c('A', 
                              'B', 
                              "C", 
                              "D"))



#library(GSEABase)
#fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
#äslim <- getOBOCollection(fl)
#myIDs<-transplant_majalis_kitzbuhl_root_up$GO.ID
#myCollection <- GOCollection(myIDs)
#goSlim(myCollection, slim, "MF")
