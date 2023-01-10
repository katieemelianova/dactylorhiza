
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
#       get root sample data into a data frame            #
###########################################################


# read in table: I got the same names from the RootRNAseq multiplex1 and multiplex2 dirs
root_samples<-read.table("/Users/katieemelianova/Desktop/Dactylorhiza/Root_samples.txt", header=FALSE, col.names = c("sample_id"))

# create new columns detailing  the species, environemt, tissue, assay and location
root_samples$species<-case_when(substr(root_samples$sample_id,1,1) == "m" ~ "majalis",
          substr(root_samples$sample_id,1,1) == "t" ~ "traunsteineri")

root_samples$environment<-case_when(substr(root_samples$sample_id,2,2) == "M" ~ "majalis",
                                    substr(root_samples$sample_id,2,2) == "T" ~ "traunsteineri")

root_samples$locality<-case_when(substr(root_samples$sample_id,3,3) == "S" ~ "St Ulrich",
                                 substr(root_samples$sample_id,3,3) == "K" ~ "Kitzbuhl")

root_samples$replicate<-substr(root_samples$sample_id,4,4)

root_samples$tissue<-case_when(substr(root_samples$sample_id,5,5) == "l" ~ "leaf",
                                 substr(root_samples$sample_id,5,5) == "r" ~ "root")

root_samples$assay<-case_when(endsWith(root_samples$sample_id, "_s") == TRUE ~ "sRNA",
          endsWith(root_samples$sample_id, "_s") == FALSE ~ "RNAseq")

# now we can filter by the tissue and assay: I think a few repeat libraries from leaf were sequenced along with the root samples
root_samples %>% filter(tissue == "root" & assay == "RNAseq") %>% dplyr::select(sample_id)


#################################################################################################
#   At this point the sample names are extracted and used to process the fastqs on lisc         #
#                         STAR mapping, featurecounts and fastqc                                #
#   Now reading the featurecounts back in here for DEseq differential expression                #
#################################################################################################

# read in each featurecounts file and join them together
setwd("~/Desktop/Dactylorhiza/dactylorhiza_root_featurecounts/")
df <- list.files(path='/Users/katieemelianova/Desktop/Dactylorhiza/dactylorhiza_root_featurecounts') %>% 
  lapply(read_tsv, skip=1) %>% 
  purrr::reduce(left_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))

# get and store the first standard featurecounts column names
fc_cols<-colnames(df)[1:6]

# take the 7th until the last colnames of the data frame
sample_names<-colnames(df)[7:length(colnames(df))]
sample_names<-str_replace(sample_names, "root_samples/", "")
sample_names<-str_replace(sample_names, "Aligned.sortedByCoord.out.bam", "")

# now set the df colnames to the shortened sample names for ease of reading
colnames(df)<-c(fc_cols, sample_names)

# get the relevant columns (geneid and sample counts) and set the geneid column to be rownames
df_counts<-df %>% 
  dplyr::select(c(1,7:length(colnames(df)))) %>% 
  column_to_rownames("Geneid")


#make the sample table, setting the sample id as the rowname
# first need to filter on assay and tissue as a few leaf samples made their way into this library, and exlude the smRNA assay
root_samples %<>% filter(tissue == "root" & assay == "RNAseq") %>% column_to_rownames("sample_id")

# remove irrelevant columns and add a transplant column to be used in the design
root_samples %<>% mutate(treatment=case_when(species == environment ~ "native",
                                            species != environment ~ "transplant")) %>% 
  dplyr::select(-c(tissue, "assay"))




######################################################################################
#   baseline DE genes between native majalis and native traunsteineri per locality   #
######################################################################################

# we want to find out what genes are DE between majalis and traunsteineri in their native areas per locality
# one comparison for native majalis and traunsteineri in St Ulrich 
# one comparison for the same but in Kitzbuhl

native_kitzbuhl<-specify_comparison(root_samples, df_counts, "treatment == 'native' & locality == 'Kitzbuhl'") %>% run_diffexp("species", df$Length)
native_stulrich<-specify_comparison(root_samples, df_counts, "treatment == 'native' & locality == 'St Ulrich'") %>% run_diffexp("species", df$Length)


############################################################################################
#        DE genes between native and transplanted species, regardless of locality       #
############################################################################################

transplant_majalis<-specify_comparison(root_samples, df_counts, 'species == "majalis"') %>% run_diffexp("treatment", df$Length)
transplant_traunsteineri<-specify_comparison(root_samples, df_counts, 'species == "traunsteineri"') %>% run_diffexp("treatment", df$Length)


############################################################################################
#        DE genes in each species when transplanted compared to native, per locality       #
############################################################################################


transplant_majalis_kitzbuhl<-specify_comparison(root_samples, df_counts, 'species == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df$Length)
transplant_majalis_stulrich<-specify_comparison(root_samples, df_counts, 'species == "majalis" & locality == "St Ulrich"') %>% run_diffexp("treatment", df$Length)
transplant_traunsteineri_kitzbuhl<-specify_comparison(root_samples, df_counts, 'species == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("treatment", df$Length)
transplant_traunsteineri_stulrich<-specify_comparison(root_samples, df_counts, 'species == "traunsteineri"& locality == "St Ulrich"') %>% run_diffexp("treatment", df$Length)

############################################################################################
#        DE genes in each environment between native species and transplanted species      #
############################################################################################

majalis_vs_traunsteineri_majalis_kitzbuhl<-specify_comparison(root_samples, df_counts, 'environment == "majalis" & locality == "Kitzbuhl"') %>% run_diffexp("species", df$Length)
majalis_vs_traunsteineri_traunsteineri_kitzbuhl<-specify_comparison(root_samples, df_counts, 'environment == "traunsteineri" & locality == "Kitzbuhl"') %>% run_diffexp("species", df$Length)
majalis_vs_traunsteineri_majalis_stulrich<-specify_comparison(root_samples, df_counts, 'environment == "majalis" & locality == "St Ulrich"') %>% run_diffexp("species", df$Length)
majalis_vs_traunsteineri_traunsteineri_stulrich<-specify_comparison(root_samples, df_counts, 'environment == "traunsteineri" & locality == "St Ulrich"') %>% run_diffexp("species", df$Length)





transplant_majalis_kitzbuhl$results %>% data.frame() %>% filter(padj < 0.05) %>% nrow()
transplant_majalis_stulrich$results %>% data.frame() %>% filter(padj < 0.05) %>% nrow()
transplant_traunsteineri_kitzbuhl$results %>% data.frame() %>% filter(padj < 0.05) %>% nrow()
transplant_traunsteineri_stulrich$results %>% data.frame() %>% filter(padj < 0.05) %>% nrow()


pval<-0.05
fold_change<-2
transplant_traunsteineri_stulrich$results %>% data.frame() %>% filter(log2FoldChange > fold_change & padj < pval) %>% rownames()

get_significant_genes(transplant_traunsteineri_stulrich, fold_change = 2, pvalue = 0.05)

results_object<-transplant_traunsteineri_stulrich

results_object$results %>% data.frame() %>% filter(log2FoldChange > integer(fold_change) & padj < integer(pval)) %>% rownames()


################################################################
#               GO term enrichment analysis                    #
################################################################

mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")




one<-transplant_majalis_kitzbuhl$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()
two<-transplant_majalis_stulrich$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()
three<-transplant_traunsteineri_kitzbuhl$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()
four<-transplant_traunsteineri_stulrich$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()

# remove the :cds bt from the ends of the gene names so they match the mappings
one<-sapply(one, function(x) str_split(x, ":")[[1]][1]) %>% unname()
two<-sapply(two, function(x) str_split(x, ":")[[1]][1]) %>% unname()
three<-sapply(three, function(x) str_split(x, ":")[[1]][1]) %>% unname()
four<-sapply(four, function(x) str_split(x, ":")[[1]][1]) %>% unname()

five<-majalis_vs_traunsteineri_majalis_kitzbuhl$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()
five<-sapply(five, function(x) str_split(x, ":")[[1]][1]) %>% unname()
majalis_vs_traunsteineri_majalis_kitzbuhl_go<-get_enriched_terms(five, mp)



majalis_kitzbuhl_trans<-get_enriched_terms(one, mp)
majalis_stulrich_trans<-get_enriched_terms(two, mp)

traunsteineri_kitzbuhl_trans<-get_enriched_terms(three, mp)
traunsteineri_sturich_trans<-get_enriched_terms(four, mp)


################################################################
#                          plot heatmaps                       #
################################################################


# native vs transplanted majalis, all localities
draw_heatmap(transplant_majalis)

# native vs transplanted traunsteineri, all localities
draw_heatmap(transplant_traunsteineri)

# native vs transplanted majalis, kitzbuhl
draw_heatmap(transplant_majalis_kitzbuhl)

# native vs transplanted majalis, St Ulrich
draw_heatmap(transplant_majalis_stulrich)

# native vs transplanted traunsteineri, kitzbuhl
draw_heatmap(transplant_traunsteineri_kitzbuhl)

# native vs transplanted traunsteineri, St Ulrich
draw_heatmap(transplant_traunsteineri_stulrich)


################################################################
#                         draw PCA plots                       #
################################################################


# make a dds object from the total root samples (no subsetting)
root_dds<-specify_comparison(root_samples, df_counts, "1 == 1")
root_dds <- DESeqDataSetFromMatrix(countData = root_dds[["counts"]],
                                   colData = root_dds[["samples"]],
                                   design = ~ treatment)

test<-varianceStabilizingTransformation(root_dds)

# remove outlier sample
test<-test[,-30]

plotPCA(test, intgroup=c("species", "locality"), ntop = 2000, returnData = FALSE)
plotPCA(test, intgroup="treatment", ntop = 1000, returnData = FALSE)
plotPCA(test, intgroup="locality", ntop = 1000, returnData = FALSE)
plotPCA(test, intgroup="environment", ntop = 1000, returnData = FALSE)

Use 52K gff




test<-fpkm(transplant_majalis_kitzbuhl$dds) %>% colnames %>% sort()
fpkm(transplant_majalis_kitzbuhl$dds) %>% data.frame() %>% dplyr::select(test) %>% colnames()


test<-fpkm(transplant_majalis_kitzbuhl$dds)
test_res<-transplant_majalis_kitzbuhl$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()

stress<-c("Dinc068413-RA",
          "Dinc122846-RA")
pheatmap(test[rownames(test) %in% stress,], scale = "row", cluster_cols = FALSE)


mcols(root_dds)$basepairs<-df$Length

root_dds <- estimateSizeFactors(root_dds)
idx <- rowSums( counts(root_dds, normalized=TRUE) >= 5 ) >= 3
root_dds <- root_dds[idx,]

# add in a gene length column for FPKM calculation

root_fpkms<-fpkm(root_dds)

root_fpkms[rownames(root_fpkms) %in% three,]

root_fpkms %>% filter(rownames(root_fpkms) %in% two)
pheatmap(root_fpkms[rownames(root_fpkms) %in% one,], scale = "row", cluster_cols = FALSE)

head(root_fpkms, 1000)



read_tsv("~/Desktop/Dactylorhiza/fungal_blast_pid_len_eval") %>%
  set_colnames(c("qseqid", "pid", "length", "evalue")) %>%
  filter(length > 100 & pid > 95) %>%
  dplyr::select(qseqid) %>%
  unique() %>%
  nrow()


#######
sarg<-read_csv("sargasso_stats.txt", col_names = c("strategy", "test", "Assigned_Hits_incarnata", "Assigned_Reads_incarnata", "Rejected_Hits_incarnata", "Rejected_Reads_incarnata", "Ambiguous_Hits_incarnata", "Ambiguous_Reads_incarnata", "Assigned_Hits_fuchsii", "Assigned_Reads_fuchsii", "Rejected_Hits_fuchsii", "Rejected_Reads_fuchsii", "Ambiguous_Hits_fuchsii", "Ambiguous_Reads_fuchsii"))

sarg %>% dplyr::select(Assigned_Reads_incarnata, Ambiguous_Reads_incarnata, Rejected_Reads_incarnata)

sarg %>% colnames()

sarg %>% mutate(incarnata_pct_accepted_reads=100*(Assigned_Reads_incarnata/(Assigned_Reads_incarnata + Ambiguous_Reads_incarnata + Rejected_Reads_incarnata))) %>% 
  mutate(fuchsii_pct_accepted_reads=100*(Assigned_Reads_fuchsii/(Assigned_Reads_fuchsii + Rejected_Reads_fuchsii + Ambiguous_Reads_fuchsii))) %>% 
  dplyr::select(strategy, incarnata_pct_accepted_reads, fuchsii_pct_accepted_reads, Assigned_Reads_incarnata, Assigned_Reads_fuchsii) 



mutate(fuchsii_pct_accepted_reads=100*(Assigned_Reads_fuchsii/Assigned_Reads_fuchsii + Rejected_Reads_fuchsii + Ambiguous_Reads_fuchsii)) %>% 
  