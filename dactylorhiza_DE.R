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




###########################################
#      Plotting constitutively DEGs       #
###########################################

pvalue_threshold <- 0.1
logfc_threshold<-0.5

tmM_stu_root<-traunsteineri_majalis_root_M_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")


tmT_stu_root<-traunsteineri_majalis_root_T_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")


tmM_ktz_root<-traunsteineri_majalis_root_M_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")


tmT_ktz_root<-traunsteineri_majalis_root_T_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")

tmM_stu_leaf<-traunsteineri_majalis_leaf_M_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")

tmT_stu_leaf<-traunsteineri_majalis_leaf_T_stulrich$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")

tmM_ktz_leaf<-traunsteineri_majalis_leaf_M_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")

tmT_ktz_leaf<-traunsteineri_majalis_leaf_T_kitzbuhl$results %>% 
  data.frame() %>% 
  filter(padj < pvalue_threshold & abs(log2FoldChange) > logfc_threshold) %>% 
  dplyr::select(log2FoldChange) %>% 
  filter(log2FoldChange < 10) %>% 
  rownames_to_column(var="gene_id")

stu_root<-inner_join(tmM_stu_root, 
                     tmT_stu_root, 
                     by="gene_id") %>% mutate(comparison="stu_root") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

ktz_root<-inner_join(tmM_ktz_root, 
                     tmT_ktz_root, 
                     by="gene_id") %>% mutate(comparison="ktz_root") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

stu_leaf<-inner_join(tmM_stu_leaf, 
                     tmT_stu_leaf, 
                     by="gene_id") %>% mutate(comparison="stu_leaf") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

ktz_leaf<-inner_join(tmM_ktz_leaf, 
                     tmT_ktz_leaf, 
                     by="gene_id") %>% mutate(comparison="ktz_leaf") %>%
  set_colnames(c("gene_id", "majalis_env", "traunst_env", "comparison"))

all_bound<-rbind(stu_root,
      ktz_root,
      stu_leaf,
      ktz_leaf)

ggplot(all_bound, aes(x=majalis_env, y=traunst_env, colour=comparison)) +
  geom_point()


#### maybe for later, splitting up the groups of genes to annotate separately 
#### but not sure if its worth it due to now sign DEGs
plot(all_bound$majalis_env, all_bound$traunst_env)
abline(a=0,b=2.5)
abline(a=0,b=-3)
abline(a=0.5,b=0)





















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
transplant_majalis_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$lengths)
transplant_majalis_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$lengths)


# heatmap
draw_heatmap(transplant_majalis_kitzbuhl_leaf)
draw_heatmap(transplant_majalis_stulrich_leaf)

# go enrichment kitzbuhl
transplant_majalis_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_kitzbuhl_leaf_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_kitzbuhl_leaf_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

# go enrichment stulrich
transplant_majalis_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_stulrich_leaf_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_stulrich_leaf_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

###############################
#            ROOT             #
###############################

# diffexp
transplant_majalis_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$lengths)
transplant_majalis_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$lengths)

#heatmap 
draw_heatmap(transplant_majalis_kitzbuhl_root)
draw_heatmap(transplant_majalis_stulrich_root)

# go enrichment kitzbuhl
transplant_majalis_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_kitzbuhl_root_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_kitzbuhl_root_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)

# go enrichment stulrich
transplant_majalis_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$down, mp) 
transplant_majalis_stulrich_root_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected, classicFisher)
transplant_majalis_stulrich_root_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(Term, classicFisher)


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



get_significant_genes(transplant_majalis_kitzbuhl_leaf) %>% length()
get_significant_genes(transplant_majalis_stulrich_leaf) %>% length()
get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf) %>% length()
get_significant_genes(transplant_traunsteineri_stulrich_leaf) %>% length()


get_significant_genes(transplant_majalis_kitzbuhl_root) %>% length()
get_significant_genes(transplant_majalis_stulrich_root) %>% length()
get_significant_genes(transplant_traunsteineri_kitzbuhl_root) %>% length()
get_significant_genes(transplant_traunsteineri_stulrich_root) %>% length()






