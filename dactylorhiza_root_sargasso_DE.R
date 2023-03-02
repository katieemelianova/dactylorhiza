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

###########################################################
#  get root sample data into a data frame using filenames #
###########################################################

root_samples_sargasso<-list.files(path='dactylorhiza_root_sargasso_featurecounts', full.names = TRUE) %>%
  str_replace("dactylorhiza_root_sargasso_featurecounts/", "") %>%
  str_replace(".featurecounts", "") %>%
  str_replace("___filtered", "") %>%
  tibble() %>%
  set_colnames("sample_id") %>%
  mutate(species=case_when(substr(sample_id,1,1) == "m" ~ "majalis",
                   substr(sample_id,1,1) == "t" ~ "traunsteineri")) %>%
  mutate(environment=case_when(substr(sample_id,2,2) == "M" ~ "majalis",
                           substr(sample_id,2,2) == "T" ~ "traunsteineri")) %>%
  mutate(locality=case_when(substr(sample_id,3,3) == "S" ~ "St Ulrich",
                               substr(sample_id,3,3) == "K" ~ "Kitzbuhl")) %>%
  mutate(subgenome=case_when(substr(sample_id,9,9) == "f" ~ "fuchsii",
                            substr(sample_id,9,9) == "i" ~ "incarnata")) %>%
  mutate(treatment=case_when(species == environment ~ "native",
                             species != environment ~ "transplant")) %>%
  column_to_rownames("sample_id")
  
  ########################################################################
#  read in and join featurecounts files, make colnames easier to read  #
########################################################################

df <- list.files(path='dactylorhiza_root_sargasso_featurecounts', full.names = TRUE) %>% 
  lapply(read_tsv, skip=1) %>% 
  purrr::reduce(left_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))

new_df_colnames<-colnames(df) %>% 
  str_replace("sargasso_root_best/filtered_reads/", "") %>%
  str_replace("___filteredAligned.sortedByCoord.out.bam", "")

colnames(df) <- new_df_colnames

#############################################################
#  remove unnecessary colnames and set geneids as rownames  #
#############################################################

df_counts <- df %>% 
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>%
  column_to_rownames("Geneid")

########################
#    draw a PCA plot   #
########################

root_sargasso_dds<-specify_comparison(root_samples_sargasso, df_counts, "1 == 1")
root_sargasso_dds <- DESeqDataSetFromMatrix(countData = root_sargasso_dds[["counts"]],
                                   colData = root_sargasso_dds[["samples"]],
                                   design = ~ species + locality + subgenome)

root_sargasso_dds<-varianceStabilizingTransformation(root_sargasso_dds)

plotPCA(root_sargasso_dds, intgroup=c("species", "subgenome"), ntop = 1000, returnData = FALSE)











