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

root_sargasso_dds <- estimateSizeFactors(root_sargasso_dds)

root_sargasso_dds_vst<-varianceStabilizingTransformation(root_sargasso_dds)

plotPCA(root_sargasso_dds_vst, intgroup=c("species", "subgenome"), ntop = 1000, returnData = FALSE)







###########################################################
#    get normalised counts and filter for low expression  #
###########################################################

# get normalised counts (default median ratio)
root_sargasso_dds_norm<-counts(root_sargasso_dds, normalized=TRUE) 



#index the genes you want to keep by minimum expression of X in N samples
idx <- rowSums( root_sargasso_dds_norm >= 5 ) >= 30 

# use that index to select out the actual rows in the normalised expression marix
root_sargasso_dds_norm <- root_sargasso_dds_norm[idx,] %>% data.frame()


# remove genes where combinaed expression is less than 5 (normalised)
#root_sargasso_dds_norm<-root_sargasso_dds_norm[rowSums(root_sargasso_dds_norm) > 5,] %>% data.frame()
tibble(root_sargasso_dds_norm) %>% colnames()



##################################################################
#     traunsteineri incarnata native vs transplant ratio plot    #
##################################################################

# get traunsteineri native incarnata and fuchsii nornalised expression values
traunsteineri_incarnata_native<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "incarnata" & treatment == "native"')$counts
traunsteineri_fuchsii_native<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "fuchsii" & treatment == "native"')$counts

# get traunsteineri transplant incarnata and fuchsii nornalised expression values
traunsteineri_incarnata_transplant<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "incarnata" & treatment == "transplant"')$counts
traunsteineri_fuchsii_transplant<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "fuchsii" & treatment == "transplant"')$counts




traunsteineri_meanexp<-data.frame(traunsteineri_incarnata_native=rowMeans(traunsteineri_incarnata_native),
           traunsteineri_fuchsii_native=rowMeans(traunsteineri_fuchsii_native),
           traunsteineri_incarnata_transplant=rowMeans(traunsteineri_incarnata_transplant),
           traunsteineri_fuchsii_transplant=rowMeans(traunsteineri_fuchsii_transplant))

# add a small number to everything in the combined dataframe to ensure non NAN result from division including 0
traunsteineri_meanexp <- traunsteineri_meanexp + 0.00001

#get the rownames to set them later
rows<-traunsteineri_meanexp %>% rownames()

# calculate the homeolog expression ratio
# h1h2native = expression of ncarnata homeolog / expression of fuchsii homeolog in native habitat
h1h2native<-log10(traunsteineri_meanexp$traunsteineri_incarnata_native/traunsteineri_meanexp$traunsteineri_fuchsii_native)
h1h2transplant<-log10(traunsteineri_meanexp$traunsteineri_incarnata_transplant/traunsteineri_meanexp$traunsteineri_fuchsii_transplant)

# I will use log transformation to nornalise the distribution of homeolog differences
# an example of how using log on different combinations of H1/H2 values gives the same result (just pos and neg)
# both homeologs same
#log(1/1)

# H1 = 10 and H2 = 50 vs H1 = 50 and H2 = 10
log(10/50)
log(50/10)

# H1 = 10 and H2 = 100 vs H1 = 100 and H2 = 10
log(10/100)
log(100/10)








h1h2_all<-data.frame(h1h2native=h1h2native,
           h1h2transplant=h1h2transplant,
           row.names = rows) %>% set_colnames(c("h1h2native", "h1h2transplant"))

#h1h2_all %<>% filter(h1h2native < 300 & h1h2transplant < 300)

h1h2_all %<>% mutate(test=case_when(h1h2native > 5 & h1h2transplant < 5 ~ "native homeolog shift",
                              h1h2transplant > 5 & h1h2native < 5 ~ "transplant homeolog shift",
                              h1h2transplant > 5 & h1h2native > 5 ~ "both homeolog shift",
                              h1h2transplant < 5 & h1h2native < 5 ~ "none homeolog shift"))


ggplot(h1h2_all, aes(x=h1h2native, y=h1h2transplant)) + 
  geom_point() +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")






# You want all points where intercept + x is more than y: 
# df$is_below <- 25 + df$x > df$y. 
# To clarify: you want all points (x, y) where the line (25 + 1*x) is above (>) the y-value. – 

plot(h1h2_all$h1h2native, h1h2_all$h1h2transplant, xlab="incarnata:fuchsii homeolog ratio native", ylab="incarnata:fuchsii homeolog ratio transplant")
abline(a=6.5,b=1)
abline(a=-7,b=1)

abline(a=-6.5,b=-1)
abline(a=6.5,b=-1)
abline(a=2.5,b=1)
abline(a=-2.3,b=1)

abline(a=-2.3,b=-1)
abline(a=-18.3,b=-1)
abline(a=17.3,b=-1)





h1h2_all %<>% mutate(cluster=case_when((6.5 + 1 * h1h2native) < h1h2transplant ~ "cluster1",
                                            (-7 + 1 * h1h2native) > h1h2transplant ~ "cluster2",
                                            (-6.5 + -1 * h1h2native) > h1h2transplant & (2.2 + 1 * h1h2native) > h1h2transplant & (-2.3 + 1 * h1h2native) < h1h2transplant ~ "cluster3",
                                            (-6.5 + -1 * h1h2native) < h1h2transplant & (6.5 + -1 * h1h2native) > h1h2transplant & (2.5 + 1 * h1h2native) > h1h2transplant & (-2.3 + 1 * h1h2native) < h1h2transplant  ~ "cluster4",
                                            (2.5 + 1 * h1h2native) < h1h2transplant & (-2.3 + -1 * h1h2native) < h1h2transplant ~ "cluster5",
                                            (2.5 + 1 * h1h2native) < h1h2transplant & (-2.3 + -1 * h1h2native) > h1h2transplant ~ "cluster6",
                                            (-2.3 + 1 * h1h2native) > h1h2transplant & (-2.3 + -1 * h1h2native) > h1h2transplant ~ "cluster7",
                                            (-2.3 + 1 * h1h2native) > h1h2transplant & (-2.3 + -1 * h1h2native) < h1h2transplant ~ "cluster8",
                                            (-2.3 + 1 * h1h2native) < h1h2transplant & (6.5 + -1 * h1h2native) < h1h2transplant & (2.5 + 1 * h1h2native) > h1h2transplant ~ "cluster9",
                                            ))



ggplot(h1h2_all, aes(x=h1h2native, y=h1h2transplant, color=cluster)) + 
  geom_point() +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")

clust1<-h1h2_all %>% filter(cluster == "cluster1") %>% rownames()
clust2<-h1h2_all %>% filter(cluster == "cluster2") %>% rownames()
clust3<-h1h2_all %>% filter(cluster == "cluster3") %>% rownames()
clust4<-h1h2_all %>% filter(cluster == "cluster4") %>% rownames()
clust5<-h1h2_all %>% filter(cluster == "cluster5") %>% rownames()
clust6<-h1h2_all %>% filter(cluster == "cluster6") %>% rownames()
clust7<-h1h2_all %>% filter(cluster == "cluster7") %>% rownames()
clust8<-h1h2_all %>% filter(cluster == "cluster8") %>% rownames()
clust9<-h1h2_all %>% filter(cluster == "cluster9") %>% rownames()



clust1_enrich<-get_enriched_terms(clust1, mp) 
clust2_enrich<-get_enriched_terms(clust2, mp) 
clust3_enrich<-get_enriched_terms(clust3, mp) 
clust4_enrich<-get_enriched_terms(clust4, mp) 
clust5_enrich<-get_enriched_terms(clust5, mp) 
clust6_enrich<-get_enriched_terms(clust6, mp)
clust7_enrich<-get_enriched_terms(clust7, mp) 
clust8_enrich<-get_enriched_terms(clust8, mp) 
clust9_enrich<-get_enriched_terms(clust9, mp) 


clust1_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust2_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust3_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust4_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust5_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust6_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust7_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust8_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()
clust9_enrich %>% tibble() %>% filter(classicFisher < 0.05) %>% dplyr::select(Term) %>% pull()

##################################################################
#           GO term enrichment of homeolog shift genes           #
##################################################################

mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")
names(mp) %<>% str_replace("-RA", "")
#subset the genes annotated with go to to only those genes in our dataset
mp<-mp[rownames(root_sargasso_dds_norm)]



transplant_shift_genes<-h1h2_all %>% filter(test == "transplant homeolog shift") %>% rownames()
get_enriched_terms(transplant_shift_genes, mp) 


test<-h1h2_all %>% filter(cluster == "cluster1") %>% rownames()


test 



test2<-root_sargasso_dds_norm[test,] %>% dplyr::select(contains("fuchsii")) %>% rownames_to_column(var="gene") %>% melt() %>% dplyr::select(-c("variable")) %>% mutate(species="fucshii")
test3<-root_sargasso_dds_norm[test,] %>% dplyr::select(contains("incarnata")) %>% rownames_to_column(var="gene") %>% melt() %>% dplyr::select(-c("variable")) %>% mutate(species="incarnata")
test4<-rbind(test2, test3)

ggplot(data = test4, aes(x=gene, y=log(value))) + geom_boxplot(aes(fill=species))













