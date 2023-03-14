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






x<-traunsteineri_majalis_root_M_stulrich$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")
y<-traunsteineri_majalis_root_T_stulrich$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")

x<-traunsteineri_majalis_root_M_kitzbuhl$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")
y<-traunsteineri_majalis_root_T_kitzbuhl$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")


x<-traunsteineri_majalis_leaf_M_stulrich$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")
y<-traunsteineri_majalis_leaf_T_stulrich$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")

x<-traunsteineri_majalis_leaf_M_kitzbuhl$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")
y<-traunsteineri_majalis_leaf_T_kitzbuhl$results %>% data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 0.5) %>% dplyr::select(log2FoldChange) %>% filter(log2FoldChange < 10) %>% rownames_to_column(var="gene_id")



z<-inner_join(x, y, by="gene_id")
colnames(z)<-c("gene_id", "majalis_env", "traunst_env")


ggplot(z, aes(x=majalis_env, y=traunst_env)) +
  geom_point()


# getv the intersection of the two lists of DE genes
# n.b. I need to change this to take into account direction of change
inter1<-intersect(get_significant_genes(traunsteineri_majalis_root_M_stulrich), get_significant_genes(traunsteineri_majalis_root_T_stulrich))
inter2<-intersect(get_significant_genes(traunsteineri_majalis_leaf_M_stulrich), get_significant_genes(traunsteineri_majalis_leaf_T_stulrich))
inter<-c(inter1, inter2)


test<-get_dds_object(specify_comparison(root_samples, df_counts_root, 'locality == "St Ulrich"'), 
               "~treatment",
               df_lengths_root,
               5,
               10)


pheatmap::pheatmap(fpkm(test)[inter,])


draw_heatmap(fpkm(test)[inter,], custom = TRUE)













transplant_majalis_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$Length)
transplant_majalis_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$Length)
transplant_traunsteineri_kitzbuhl_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_leaf$Length)
transplant_traunsteineri_stulrich_leaf<-specify_comparison(leaf_samples, df_counts_leaf, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df_leaf$Length)

transplant_majalis_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$Length)
transplant_majalis_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$Length)
transplant_traunsteineri_kitzbuhl_root<-specify_comparison(root_samples, df_counts_root, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df_root$Length)
transplant_traunsteineri_stulrich_root<-specify_comparison(root_samples, df_counts_root, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df_root$Length)


get_significant_genes(transplant_majalis_kitzbuhl_leaf) %>% length()
get_significant_genes(transplant_majalis_stulrich_leaf) %>% length()
get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf) %>% length()
get_significant_genes(transplant_traunsteineri_stulrich_leaf) %>% length()


get_significant_genes(transplant_majalis_kitzbuhl_root) %>% length()
get_significant_genes(transplant_majalis_stulrich_root) %>% length()
get_significant_genes(transplant_traunsteineri_kitzbuhl_root) %>% length()
get_significant_genes(transplant_traunsteineri_stulrich_root) %>% length()


