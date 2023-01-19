
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


#### testing all together multi factor with interaction

all_samples<-specify_comparison(root_samples, df_counts, "1 == 1")
all_samples$samples
design_term<-"species + treatment + species:treatment"
formula_parsed<-paste("~", design_term)
dds <- DESeqDataSetFromMatrix(countData = all_samples[["counts"]],
                              colData = all_samples[["samples"]],
                              design = as.formula(formula_parsed))
mcols(dds)$basepairs<-df$Length
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 8 ) >= 5
dds <- dds[idx,]
dds <- DESeq(dds)
res <- results(dds)




### now trying to specify fancier designs using the functions defined in dact_funxtions.R
all_samples_treatment<-specify_comparison(root_samples, df_counts, "1 == 1") %>% run_diffexp("treatment", df$Length)
all_samples_species_treatment<-specify_comparison(root_samples, df_counts, "1 == 1") %>% run_diffexp("species + treatment", df$Length)
all_samples_species_locality_treatment<-specify_comparison(root_samples, df_counts, "1 == 1") %>% run_diffexp("species + locality + treatment", df$Length)
all_samples_species_treatment_interaction<-specify_comparison(root_samples, df_counts, "1 == 1") %>% run_diffexp("species + treatment + species:treatment", df$Length)

# one way of getting DE genes for a partcular effect
results(all_samples_species_treatment$dds, contrast = c("species", "majalis", "traunsteineri")) %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% nrow()
results(all_samples_species_treatment$dds, contrast = c("treatment", "native", "transplant")) %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% nrow()

# list the names of different effects in the design
resultsNames(all_samples_species_treatment$dds)  

# get results for one of the effects only 
results(all_samples_species_treatment$dds, name = "species_traunsteineri_vs_majalis")

# draw a heatmap of each comaprison
draw_heatmap(all_samples_treatment)
draw_heatmap(all_samples_species_treatment)
draw_heatmap(all_samples_species_locality_treatment)
draw_heatmap(all_samples_species_treatment_interaction)

# get sig DE genes for each model fitted
a<-get_significant_genes(all_samples_species)
b<-get_significant_genes(all_samples_species_treatment)
c<-get_significant_genes(all_samples_species_locality_treatment)
d<-get_significant_genes(all_samples_species_treatment_interaction)

length(a)
length(b)
length(c)
length(d)

# plot a venn diagram of the genes significantly DE in each fitted model
venn<-venn.diagram(
  x = list(a, b, c, d),
  category.names = c("trt" , "sp+trt" , "sp+loc+trt", "sp+trt+sp:trt"),
  filename = NULL)

pdf(file="venn.pdf")
grid.draw(venn)
dev.off()

intersect(a,b) %>% length()
intersect(a,c) %>% length()
intersect(b,c) %>% length()

mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")
test_species<-get_significant_genes(all_samples_species, mappings_format = TRUE)
test_species_treatment<-get_significant_genes(all_samples_species_treatment, mappings_format=TRUE)
test_species_treatment_interaction<-get_significant_genes(all_samples_species_treatment_interaction, mappings_format=TRUE)

test_enrich_species<-get_enriched_terms(test_species, mp) %>% data.frame() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected)
test_enrich_species_treatment<-get_enriched_terms(test_species_treatment, mp) %>% data.frame() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected)
test_enrich_species_treatment_interaction<-get_enriched_terms(test_species_treatment_interaction, mp) %>% data.frame() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected)


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


################################################################
#               GO term enrichment analysis                    #
################################################################

mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")
test<-get_significant_genes(transplant_traunsteineri, mappings_format=TRUE)

#transplant_traunsteineri
#transplant_majalis_kitzbuhl
#transplant_majalis_stulrich
#transplant_traunsteineri_kitzbuhl
#transplant_traunsteineri_stulrich


test_enrich<-get_enriched_terms(test, mp) %>% data.frame() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected)

test_enrich %>% data.frame() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term, Annotated, Significant, Expected)

majalis_kitzbuhl_trans<-get_enriched_terms(one, mp)
majalis_stulrich_trans<-get_enriched_terms(two, mp)

traunsteineri_kitzbuhl_trans<-get_enriched_terms(three, mp)
traunsteineri_sturich_trans<-get_enriched_terms(four, mp)


################################################################
#                          plot heatmaps                       #
################################################################


# native vs transplanted majalis, all localities
#draw_heatmap(transplant_majalis)
transplant_majalis_enrich<-return_enrichment_table(transplant_majalis)

# native vs transplanted traunsteineri, all localities
#draw_heatmap(transplant_traunsteineri)
transplant_traunsteineri_enrich<-return_enrichment_table(transplant_traunsteineri)

# native vs transplanted majalis, kitzbuhl
#draw_heatmap(transplant_majalis_kitzbuhl)
transplant_majalis_kitzbuhl_enrich<-return_enrichment_table(transplant_majalis_kitzbuhl)

# native vs transplanted majalis, St Ulrich
#draw_heatmap(transplant_majalis_stulrich)
transplant_majalis_stulrich_enrich<-return_enrichment_table(transplant_majalis_stulrich)

# native vs transplanted traunsteineri, kitzbuhl
#draw_heatmap(transplant_traunsteineri_kitzbuhl)
transplant_traunsteineri_kitzbuhl_enrich<-return_enrichment_table(transplant_traunsteineri_kitzbuhl)

# native vs transplanted traunsteineri, St Ulrich
#draw_heatmap(transplant_traunsteineri_stulrich)
transplant_traunsteineri_stulrich_enrich<-return_enrichment_table(transplant_traunsteineri_stulrich)



################################################################
#                         draw PCA plots                       #
################################################################

model.matrix(~test$species + test$locality)

# make a dds object from the total root samples (no subsetting)
root_dds<-specify_comparison(root_samples, df_counts, "1 == 1")
root_dds <- DESeqDataSetFromMatrix(countData = root_dds[["counts"]],
                                   colData = root_dds[["samples"]],
                                   design = ~ species + locality)

test<-varianceStabilizingTransformation(root_dds)
test$locality
test_limma<-limma::removeBatchEffect(assay(test), test$locality, design=model.matrix(~test$species + test$locality))




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
  