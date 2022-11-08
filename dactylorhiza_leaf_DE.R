library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(topGO)
library(pheatmap)
source("/Users/katieemelianova/Desktop/Dactylorhiza/dactylorhiza/dactylorhiza_functions.R")


###########################################################
#       get leaf sample data into a data frame            #
###########################################################

leaf_samples<-read.table("/Users/katieemelianova/Desktop/Dactylorhiza/Leaf_samples.txt", header=FALSE, col.names = c("sample_id"))
# create new columns detailing  the species, environemt, tissue, assay and location
leaf_samples$species<-case_when(substr(leaf_samples$sample_id,1,1) == "m" ~ "majalis",
                                substr(leaf_samples$sample_id,1,1) == "t" ~ "traunsteineri")
leaf_samples$locality<-case_when(substr(leaf_samples$sample_id,3,3) == "S" ~ "St Ulrich",
                                 substr(leaf_samples$sample_id,3,3) == "K" ~ "Kitzbuhl")
leaf_samples$environment<-case_when(substr(leaf_samples$sample_id,2,2) == "M" ~ "majalis",
                                    substr(leaf_samples$sample_id,2,2) == "T" ~ "traunsteineri")
leaf_samples$replicate<-substr(leaf_samples$sample_id,4,4)




####################################################################
#               Now reading in reads for leaf                      #
####################################################################


# read in each featurecounts file and join them together
setwd("/Users/katieemelianova/Desktop/Dactylorhiza/dactylorhiza_leaf_featurecounts")
df_leaf <- list.files(path='/Users/katieemelianova/Desktop/Dactylorhiza/dactylorhiza_leaf_featurecounts') %>% 
  lapply(read_tsv, skip=1) %>% 
  purrr::reduce(left_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))

# get and store the first standard featurecounts column names
fc_cols<-colnames(df_leaf)[1:6]

# take the 7th until the last colnames of the data frame
sample_names<-colnames(df_leaf)[7:length(colnames(df_leaf))]
sample_names<-str_replace(sample_names, "StUlrich/", "")
sample_names<-str_replace(sample_names, "Kitzbuhel/", "")
sample_names<-str_replace(sample_names, "Aligned.sortedByCoord.out.bam", "")

# now set the df colnames to the shortened sample names for ease of reading
colnames(df_leaf)<-c(fc_cols, sample_names)

# get the relevant columns (geneid and sample counts) and set the geneid column to be rownames
df_counts_leaf<-df_leaf %>% 
  dplyr::select(c(1,7:length(colnames(df_leaf)))) %>% 
  column_to_rownames("Geneid")

#make the sample table, setting the sample id as the rowname
leaf_samples %<>% column_to_rownames("sample_id")

# remove irrelevant columns and add a transplant column to be used in the design
leaf_samples %<>% mutate(treatment=case_when(species == environment ~ "native",
                                             species != environment ~ "transplant"))



