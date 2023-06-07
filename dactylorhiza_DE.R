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
source("dactylorhiza_functions.R")


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


###########################################################
#          annotate samples using above function          #
###########################################################
leaf_samples<-annotate_samples(leaf_samples)
root_samples<-annotate_samples(root_samples)



####################################################################################################
#  define unwanted strings and read in featurecounts, joining and removing unwanted string         #
####################################################################################################

strings_to_remove_root<-c("root_samples/", "Aligned.sortedByCoord.out.bam")
df_root<-read_in_featurecounts('dactylorhiza_root_featurecounts', strings_to_remove_root)
df_counts_root<-df_root$counts
df_lengths_root<-df_root$lengths

strings_to_remove_leaf<-c("StUlrich/", "Kitzbuhel/", "Aligned.sortedByCoord.out.bam")
df_leaf<-read_in_featurecounts('dactylorhiza_leaf_featurecounts', strings_to_remove_leaf)
df_counts_leaf<-df_leaf$counts
df_lengths_leaf<-df_leaf$lengths








###############################
#       Draw PCA plots        #
###############################

# make a dds object from the total root samples (no subsetting)
root_dds<-specify_comparison(root_samples, df_counts_root, "1 == 1")
root_dds <- DESeqDataSetFromMatrix(countData = root_dds[["counts"]],
                                   colData = root_dds[["samples"]],
                                   design = ~ species + locality)

# make a dds object from the total leaf samples (no subsetting)
leaf_dds<-specify_comparison(leaf_samples, df_counts_leaf, "1 == 1")
leaf_dds <- DESeqDataSetFromMatrix(countData = leaf_dds[["counts"]],
                                   colData = leaf_dds[["samples"]],
                                   design = ~ species + locality)


# perform variance stabilizig transformation
root_vst<-varianceStabilizingTransformation(root_dds)
leaf_vst<-varianceStabilizingTransformation(leaf_dds)


# remove outlier sample
root_vst<-root_vst[,-30]

# plot the PCA 
pcaData <-plotPCA(root_vst, intgroup=c("treatment", "species", "locality"), ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=locality, shape=interaction(species, treatment))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  theme(legend.title = element_blank()) 



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

draw_heatmap2(traunsteineri_majalis_root_M_kitzbuhl$dds)
draw_heatmap2(traunsteineri_majalis_root_T_kitzbuhl$dds)

# leaf 
traunsteineri_majalis_leaf_M_kitzbuhl<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("species", df_lengths_leaf)
traunsteineri_majalis_leaf_T_kitzbuhl<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("species", df_lengths_leaf)

draw_heatmap2(traunsteineri_majalis_leaf_M_kitzbuhl$dds)
draw_heatmap2(traunsteineri_majalis_leaf_T_kitzbuhl$dds)



# getv the intersection of the two lists of DE genes
# n.b. I need to change this to take into account direction of change
intersect(get_significant_genes(traunsteineri_majalis_leaf_M_kitzbuhl), get_significant_genes(traunsteineri_majalis_leaf_T_kitzbuhl))
intersect(get_significant_genes(traunsteineri_majalis_root_M_kitzbuhl), get_significant_genes(traunsteineri_majalis_root_T_kitzbuhl))


####################
#    ST ULRICH     #
#################### 


# root 
traunsteineri_majalis_root_M_stulrich<-specify_comparison(root_samples, df_counts_root, 'environment == "majalis" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_root)
traunsteineri_majalis_root_T_stulrich<-specify_comparison(root_samples, df_counts_root, 'environment == "traunsteineri" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_root)

# leaf 
traunsteineri_majalis_leaf_M_stulrich<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "majalis" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_leaf)
traunsteineri_majalis_leaf_T_stulrich<-specify_comparison(leaf_samples, df_counts_leaf, 'environment == "traunsteineri" & locality == "St Ulrich"') %>% run_diffexp("species", df_lengths_leaf)


draw_heatmap2(traunsteineri_majalis_root_M_stulrich$dds)
draw_heatmap2(traunsteineri_majalis_root_T_stulrich$dds)

draw_heatmap2(traunsteineri_majalis_leaf_M_stulrich$dds)
draw_heatmap2(traunsteineri_majalis_leaf_T_stulrich$dds)




results(traunsteineri_majalis_root_M_stulrich$dds, contrast=c("treatment", "native", "transplant"))
(traunsteineri_majalis_root_M_stulrich$dds)

traunsteineri_majalis_root_M_stulrich$dds$
###########################################
#      Plotting constitutively DEGs       #
###########################################

pvalue_threshold <- 0.05
#logfc_threshold<-1

tmM_stu_root<-traunsteineri_majalis_root_M_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")


tmT_stu_root<-traunsteineri_majalis_root_T_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

tmM_ktz_root<-traunsteineri_majalis_root_M_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

tmT_ktz_root<-traunsteineri_majalis_root_T_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

tmM_stu_leaf<-traunsteineri_majalis_leaf_M_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

tmT_stu_leaf<-traunsteineri_majalis_leaf_T_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

tmM_ktz_leaf<-traunsteineri_majalis_leaf_M_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

tmT_ktz_leaf<-traunsteineri_majalis_leaf_T_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  rownames_to_column(var="gene_id")

stu_root<-inner_join(tmM_stu_root, 
                     tmT_stu_root, 
                     by="gene_id") %>% mutate(comparison="St. Ulrich root") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

ktz_root<-inner_join(tmM_ktz_root, 
                     tmT_ktz_root, 
                     by="gene_id") %>% mutate(comparison="Kitzbuhl root") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

stu_leaf<-inner_join(tmM_stu_leaf, 
                     tmT_stu_leaf, 
                     by="gene_id") %>% mutate(comparison="St. Ulrich leaf") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

ktz_leaf<-inner_join(tmM_ktz_leaf, 
                     tmT_ktz_leaf, 
                     by="gene_id") %>% mutate(comparison="Kitzbuhl leaf") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

all_bound<-rbind(stu_root,
      ktz_root,
      stu_leaf,
      ktz_leaf)


all_bound %<>% mutate(status=case_when(majalis_env > 2 & traunst_env > 2 | majalis_env < -2 & traunst_env < -2 ~ "DE"))
all_bound %<>% replace_na(list(status = "Not DE"))
                     
ggplot(all_bound, aes(x=majalis_env, y=traunst_env, colour=comparison)) +
  geom_point() + 
  ylab("Fold change in traunsteineri environment (p < 0.05)") +
  xlab("Fold change in majalis environment (p < 0.05)")
















# for the GO term enrichment tests
mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")

names(mp) <- names(mp) %>% 
  str_remove("-RA") %>%
  str_remove("-RB") %>%
  str_remove("-RC")
  
  


#############################################################
#             Majalis effect of transplantation             #
#############################################################

###############################
#            LEAF             #
###############################

# diffexp
transplant_majalis_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$lengths, cpm_threshold=1, min_count_per_sample=5)
transplant_majalis_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$lengths, cpm_threshold=1, min_count_per_sample=5)

# heatmap
draw_heatmap2(transplant_majalis_kitzbuhl_leaf$dds)
draw_heatmap2(transplant_majalis_stulrich_leaf$dds)


#draw_heatmap(transplant_majalis_kitzbuhl_leaf)
#draw_heatmap(transplant_majalis_stulrich_leaf)



plotMA(transplant_majalis_kitzbuhl_leaf$results, ylim=c(-10,10))
plotMA(transplant_majalis_stulrich_leaf$results, ylim=c(-10,10))



transplant_majalis_kitzbuhl_leaf$results %>% data.frame() %>% filter(log2FoldChange > 5 & padj < 0.05)

plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc060018", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc069484", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc089174", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc30703", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc31917", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc41693", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc127790", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc135798", intgroup="treatment")
plotCounts(transplant_majalis_kitzbuhl_leaf$dds, gene="Dinc134292", intgroup="treatment")




weird<-transplant_majalis_kitzbuhl_leaf$results %>% data.frame() %>% filter(abs(log2FoldChange) > 5) %>% rownames()

df_counts_leaf %>% dplyr::select("mTK1") %>% rownames_to_column(var="gene_name") %>% filter(gene_name %in% weird & mTK1 > 1000)

# go enrichment kitzbuhl
transplant_majalis_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_kitzbuhl_leaf_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_kitzbuhl_leaf_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

# go enrichment stulrich
transplant_majalis_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_stulrich_leaf_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_stulrich_leaf_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

intersect(transplant_majalis_kitzbuhl_leaf_up$Term, transplant_majalis_stulrich_leaf_up$Term)
intersect(transplant_majalis_kitzbuhl_leaf_down$Term, transplant_majalis_stulrich_leaf_down$Term)


###############################
#            ROOT             #
###############################

# diffexp
transplant_majalis_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$lengths)
transplant_majalis_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$lengths)

#heatmap 
#draw_heatmap(transplant_majalis_kitzbuhl_root)
#draw_heatmap(transplant_majalis_stulrich_root)

draw_heatmap2(transplant_majalis_kitzbuhl_leaf$dds)
draw_heatmap2(transplant_majalis_stulrich_leaf$dds)


dds_fpkm<-fpkm(transplant_majalis_kitzbuhl_root$dds)
new_column_order<-dds_fpkm %>% colnames %>% sort()
dds_fpkm <- dds_fpkm %>% data.frame() %>% dplyr::select(new_column_order)
dds_significant<-transplant_majalis_kitzbuhl_root$results %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% rownames()
dds_toplot<-dds_fpkm[rownames(dds_fpkm) %in% dds_significant,]



sapply(c("treatment"), function(x) transplant_majalis_kitzbuhl_root[["dds"]][[x]]) %>% data.frame()

transplant_majalis_kitzbuhl_root[["dds"]]$treatment

# go enrichment kitzbuhl
transplant_majalis_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_kitzbuhl_root_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_kitzbuhl_root_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

# go enrichment stulrich
transplant_majalis_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_stulrich_root_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_stulrich_root_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

intersect(transplant_majalis_stulrich_root_up$Term, transplant_majalis_kitzbuhl_root_up$Term)
intersect(transplant_majalis_stulrich_root_down$Term, transplant_majalis_kitzbuhl_root_down$Term)




#############################################################
#         Traunsteineri effect of transplantation           #
#############################################################

###############################
#            LEAF             #
###############################

# diffexp
transplant_traunsteineri_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$lengths)
transplant_traunsteineri_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$lengths)

# heatmap
draw_heatmap(transplant_traunsteineri_kitzbuhl_leaf)
draw_heatmap(transplant_traunsteineri_stulrich_leaf)
draw_heatmap2(transplant_traunsteineri_kitzbuhl_leaf$dds)
draw_heatmap2(transplant_traunsteineri_stulrich_leaf$dds)


read.table("~/Desktop/Diospyros/impolita_seqlengths.txt") %>% dplyr::select(V2) %>% pull() %>% summary()

# go enrichment kitzbuhl
transplant_traunsteineri_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_traunsteineri_kitzbuhl_leaf_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_traunsteineri_kitzbuhl_leaf_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

# go enrichment stulrich
transplant_traunsteineri_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_traunsteineri_stulrich_leaf_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_traunsteineri_stulrich_leaf_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)


intersect(transplant_traunsteineri_kitzbuhl_leaf_up$Term, transplant_traunsteineri_stulrich_leaf_up$Term)
intersect(transplant_traunsteineri_stulrich_leaf_down$Term, transplant_traunsteineri_kitzbuhl_leaf_down$Term)



###############################
#            ROOT             #
###############################

# diffexp
transplant_traunsteineri_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$lengths)
transplant_traunsteineri_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$lengths)

# heatmap
draw_heatmap(transplant_traunsteineri_kitzbuhl_root)
draw_heatmap(transplant_traunsteineri_stulrich_root)

# go enrichment kitzbuhl
transplant_traunsteineri_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_traunsteineri_kitzbuhl_root_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_traunsteineri_kitzbuhl_root_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

# go enrichment stulrich
transplant_traunsteineri_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_traunsteineri_stulrich_root_up %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_traunsteineri_stulrich_root_down %<>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

intersect(transplant_traunsteineri_kitzbuhl_root_up$Term, transplant_traunsteineri_stulrich_root_up$Term)
intersect(transplant_traunsteineri_stulrich_root_down$Term, transplant_traunsteineri_kitzbuhl_root_down$Term)




#############################################
#         Count number of DE genes          #
#############################################

fc=5
pv=0.005

get_significant_genes(transplant_majalis_kitzbuhl_leaf, fold_change=fc, pvalue=pv) %>% length()
get_significant_genes(transplant_majalis_stulrich_leaf, fold_change=fc, pvalue=pv) %>% length()
get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, fold_change=fc, pvalue=pv) %>% length()
get_significant_genes(transplant_traunsteineri_stulrich_leaf, fold_change=fc, pvalue=pv) %>% length()




get_significant_genes(transplant_majalis_kitzbuhl_root, fold_change=fc, pvalue=pv) %>% length()
get_significant_genes(transplant_majalis_stulrich_root, fold_change=fc, pvalue=pv) %>% length()
get_significant_genes(transplant_traunsteineri_kitzbuhl_root, fold_change=fc, pvalue=pv) %>% length()
get_significant_genes(transplant_traunsteineri_stulrich_root, fold_change=fc, pvalue=pv) %>% length()


intersect(get_significant_genes(transplant_majalis_kitzbuhl_root), get_significant_genes(transplant_majalis_stulrich_root))
intersect(get_significant_genes(transplant_traunsteineri_kitzbuhl_root), get_significant_genes(transplant_traunsteineri_stulrich_root))






a<-get_significant_genes(transplant_majalis_kitzbuhl_leaf)
b<-get_significant_genes(transplant_majalis_stulrich_leaf)
c<-get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf)
d<-get_significant_genes(transplant_traunsteineri_stulrich_leaf)



library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

# Chart
venn.diagram(
  x = list(a, b, c, d),
  category.names = c("mK", "mS", "tK", "tS"),
  filename = 'venv_orchid_de_leaf',
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
  filename = 'venv_orchid_de_root',
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



df_root$lengths



##### experimenting with manual DEseq to see what the heatmap problem is, and why its so variable between reps

test<-specify_comparison(root_samples, df_counts_root, 'environment == "majalis" & locality == "St Ulrich"')

tes2t<-specify_comparison(root_samples, df_counts_root, '1 == 1')

colSums(tes2t$counts) %>% data.frame()

dds <- DESeqDataSetFromMatrix(countData = test$counts,
                              colData = test$samples,
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

