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
library(qgraph)
library(dodgr)
library(reshape2)
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

#index the genes you want to keep by minimum expression of 5 in 20 samples
idx_keep <- rowSums( root_sargasso_dds_norm >= 5 ) >= 20 

#index the genes you want to remove due to crazy high expression (more than 1000 in at least one sample)
idx_remove <- rowSums( root_sargasso_dds_norm >= 100000 ) < 2


# use that index to select out the actual rows in the normalised expression marix
root_sargasso_dds_norm <- root_sargasso_dds_norm[idx_keep,] %>% data.frame()
#root_sargasso_dds_norm<-root_sargasso_dds_norm[!idx_remove,] %>% data.frame()

# remove genes where combinaed expression is less than 5 (normalised)
#root_sargasso_dds_norm<-root_sargasso_dds_norm[rowSums(root_sargasso_dds_norm) > 5,] %>% data.frame()


##################################################################
#     traunsteineri incarnata native vs transplant ratio plot    #
##################################################################

# get traunsteineri native incarnata and fuchsii nornalised expression values
traunsteineri_incarnata_native<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "incarnata" & treatment == "native"')$counts
traunsteineri_fuchsii_native<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "fuchsii" & treatment == "native"')$counts

# get traunsteineri transplant incarnata and fuchsii nornalised expression values
traunsteineri_incarnata_transplant<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "incarnata" & treatment == "transplant"')$counts
traunsteineri_fuchsii_transplant<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "traunsteineri" & subgenome == "fuchsii" & treatment == "transplant"')$counts


# quickly changing it to majalis to see what the plot looks like, delete and rewrite later
majalis_incarnata_native<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "majalis" & subgenome == "incarnata" & treatment == "native"')$counts
majalis_fuchsii_native<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "majalis" & subgenome == "fuchsii" & treatment == "native"')$counts

# get traunsteineri transplant incarnata and fuchsii nornalised expression values
majalis_incarnata_transplant<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "majalis" & subgenome == "incarnata" & treatment == "transplant"')$counts
majalis_fuchsii_transplant<-specify_comparison(root_samples_sargasso, root_sargasso_dds_norm, 'species == "majalis" & subgenome == "fuchsii" & treatment == "transplant"')$counts


meanexp<-data.frame(
  traunsteineri_incarnata_native = rowMeans(traunsteineri_incarnata_native),
  traunsteineri_fuchsii_native = rowMeans(traunsteineri_fuchsii_native),
  traunsteineri_incarnata_transplant = rowMeans(traunsteineri_incarnata_transplant),
  traunsteineri_fuchsii_transplant = rowMeans(traunsteineri_fuchsii_transplant),
  majalis_incarnata_native = rowMeans(majalis_incarnata_native),
  majalis_fuchsii_native = rowMeans(majalis_fuchsii_native),
  majalis_incarnata_transplant = rowMeans(majalis_incarnata_transplant),
  majalis_fuchsii_transplant = rowMeans(majalis_fuchsii_transplant)
)

# add a small number to everything in the combined dataframe to ensure non NAN result from division including 0
meanexp <- meanexp + 0.00001

#get the rownames to set them later
rows<-meanexp %>% rownames()

# calculate the homeolog expression ratio
# h1h2native = expression of ncarnata homeolog / expression of fuchsii homeolog in native habitat
h1h2native_traunsteineri<-log10(meanexp$traunsteineri_incarnata_native/meanexp$traunsteineri_fuchsii_native)
h1h2transplant_traunsteineri<-log10(meanexp$traunsteineri_incarnata_transplant/meanexp$traunsteineri_fuchsii_transplant)

h1h2native_majalis<-log10(meanexp$majalis_incarnata_native/meanexp$majalis_fuchsii_native)
h1h2transplant_majalis<-log10(meanexp$majalis_incarnata_transplant/meanexp$majalis_fuchsii_transplant)




# I will use log transformation to nornalise the distribution of homeolog differences
# an example of how using log on different combinations of H1/H2 values gives the same result (just pos and neg)
# both homeologs same
#log(1/1)

# H1 = 10 and H2 = 50 vs H1 = 50 and H2 = 10
#log(10/50)
#log(50/10)

# H1 = 10 and H2 = 100 vs H1 = 100 and H2 = 10
#log(10/100)
#log(100/10)








h1h2_all<-data.frame(h1h2native_traunsteineri=h1h2native_traunsteineri,
                     h1h2transplant_traunsteineri=h1h2transplant_traunsteineri,
                     h1h2native_majalis=h1h2native_majalis,
                     h1h2transplant_majalis=h1h2transplant_majalis,
           row.names = rows) %>% set_colnames(c("h1h2native_traunsteineri", 
                                                "h1h2transplant_traunsteineri", 
                                                "h1h2native_majalis",
                                                "h1h2transplant_majalis"))

#h1h2_all %<>% filter(h1h2native < 300 & h1h2transplant < 300)

#h1h2_all %<>% mutate(test=case_when(h1h2native > 5 & h1h2transplant < 5 ~ "native homeolog shift",
#                              h1h2transplant > 5 & h1h2native < 5 ~ "transplant homeolog shift",
#                              h1h2transplant > 5 & h1h2native > 5 ~ "both homeolog shift",
#                              h1h2transplant < 5 & h1h2native < 5 ~ "none homeolog shift"))


ggplot(h1h2_all, aes(x=h1h2native_traunsteineri, y=h1h2transplant_traunsteineri)) + 
  geom_point() +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")

ggplot(h1h2_all, aes(x=h1h2native_majalis, y=h1h2transplant_majalis)) + 
  geom_point() +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")





# You want all points where intercept + x is more than y: 
# df$is_below <- 25 + df$x > df$y. 
# To clarify: you want all points (x, y) where the line (25 + 1*x) is above (>) the y-value. – 

plot(h1h2_all$h1h2native, h1h2_all$h1h2transplant, xlab="incarnata:fuchsii homeolog ratio native", ylab="incarnata:fuchsii homeolog ratio transplant")
abline(a=6.5,b=1)
abline(a=-7,b=1)

abline(a=-8.5,b=-1)
abline(a=8.5,b=-1)
abline(a=2.5,b=1)
abline(a=-3.1,b=1)

abline(a=-2.3,b=-1)
abline(a=-18.3,b=-1)
abline(a=17.3,b=-1)





h1h2_all %<>% mutate(traunsteineri_cluster=case_when((-8.5 + -1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri & (2.2 + 1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri & (-2.3 + 1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri ~ "cluster1",
                                        (-10.5 + -1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri & (10.5 + -1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri & (2.5 + 1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri & (-3 + 1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri  ~ "cluster2",
                                        (2.5 + 1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri & (-2.3 + -1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri ~ "cluster3",
                                        (2.5 + 1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri & (-2.3 + -1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri ~ "cluster4",
                                        (-2.3 + 1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri & (-2.3 + -1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri ~ "cluster5",
                                        (-2.3 + 1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri & (-2.3 + -1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri ~ "cluster6",
                                        (-2.3 + 1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri & (6.5 + -1 * h1h2native_traunsteineri) < h1h2transplant_traunsteineri & (2.5 + 1 * h1h2native_traunsteineri) > h1h2transplant_traunsteineri ~ "cluster7",
                                        ))

h1h2_all %<>% mutate(majalis_cluster=case_when((-8.5 + -1 * h1h2native_majalis) > h1h2transplant_majalis & (2.2 + 1 * h1h2native_majalis) > h1h2transplant_majalis & (-2.3 + 1 * h1h2native_majalis) < h1h2transplant_majalis ~ "cluster1",
                                                     (-10.5 + -1 * h1h2native_majalis) < h1h2transplant_majalis & (10.5 + -1 * h1h2native_majalis) > h1h2transplant_majalis & (2.5 + 1 * h1h2native_majalis) > h1h2transplant_majalis & (-3 + 1 * h1h2native_majalis) < h1h2transplant_majalis  ~ "cluster2",
                                                     (2.5 + 1 * h1h2native_majalis) < h1h2transplant_majalis & (-2.3 + -1 * h1h2native_majalis) < h1h2transplant_majalis ~ "cluster3",
                                                     (2.5 + 1 * h1h2native_majalis) < h1h2transplant_majalis & (-2.3 + -1 * h1h2native_majalis) > h1h2transplant_majalis ~ "cluster4",
                                                     (-2.3 + 1 * h1h2native_majalis) > h1h2transplant_majalis & (-2.3 + -1 * h1h2native_majalis) > h1h2transplant_majalis ~ "cluster5",
                                                     (-2.3 + 1 * h1h2native_majalis) > h1h2transplant_majalis & (-2.3 + -1 * h1h2native_majalis) < h1h2transplant_majalis ~ "cluster6",
                                                     (-2.3 + 1 * h1h2native_majalis) < h1h2transplant_majalis & (6.5 + -1 * h1h2native_majalis) < h1h2transplant_majalis & (2.5 + 1 * h1h2native_majalis) > h1h2transplant_majalis ~ "cluster7",
))


ggplot(h1h2_all, aes(x=h1h2native_traunsteineri, y=h1h2transplant_traunsteineri, color=traunsteineri_cluster)) + 
  geom_point() +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")



ggplot(h1h2_all, aes(x=h1h2native_majalis, y=h1h2transplant_majalis, color=majalis_cluster)) + 
  geom_point() +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")


# make a table of the number of genes in each cluster per species
h1h2_all %>% dplyr::select(traunsteineri_cluster) %>% pull() %>% table()
h1h2_all %>% dplyr::select(majalis_cluster) %>% pull() %>% table()


h1h2_all %<>% mutate(membership=case_when(traunsteineri_cluster == majalis_cluster ~ "same",
                              traunsteineri_cluster != majalis_cluster ~"different"))


ggplot(h1h2_all, aes(x=h1h2native_majalis, y=h1h2transplant_majalis, color=membership)) + 
  geom_point(size=0.5) +
  ylab("incarnata:fuchsii homeolog ratio transplant") + 
  xlab("incarnata:fuchsii homeolog ratio native")

# traunsteineri  clusters
h1h2_all %>% 
  dplyr::select(traunsteineri_cluster, majalis_cluster) %>% 
  table() %>%
  rowSums()/8980 * 100

# majalis  clusters
h1h2_all %>% 
  dplyr::select(traunsteineri_cluster, majalis_cluster) %>% 
  table() %>%
  colSums()/8980 * 100



h1h2_all %>% 
  dplyr::select(traunsteineri_cluster, majalis_cluster) %>% 
  table() %>% 
  pheatmap(scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)

test <- h1h2_all %>% 
  dplyr::select(traunsteineri_cluster, majalis_cluster) %>% 
  filter(traunsteineri_cluster != majalis_cluster) %>%
  data.frame() %>%
  dcast(traunsteineri_cluster~majalis_cluster) %>%
  column_to_rownames("traunsteineri_cluster") %>% 
  as.matrix() %>% 
  melt() %>% 
  set_colnames(c("traunsteineri", "majalis", "count"))








# this clunky bit of code is to colour the edges of the network graph by the direction of comparison
# i.e. compare cluster 1 to cluster 2 vs cluster 2 to cluster 1
bin1=list()
bin2=list()
cols=c()
for (i in 1:nrow(mg)){
  # turn each row of two columns into a vector of two
  row_vector<-(as.numeric(mg[i,]))

  if (!(list(row_vector) %in% bin1) & !(list(rev(row_vector)) %in% bin1)) {
    # turn each vector of two into a list, append it to a prepared list outside of loop
    bin1<-append(bin1, c(list(row_vector)))
    # put the reverse of the row vector (opposite direction comparison) into the other bin
    bin2<-append(bin2, c(list(rev(row_vector))))
  }
  # after putting the vectors and their inverse into different bins
  # check which bin the vector is in
  # if in the first bin, its the first way around, pick it a colour
  #if in the second bin, pick it a different coloir
  if ((list(row_vector) %in% bin1)) {
    cols=c(cols, "purple")
  } else if ((list(row_vector) %in% bin2)){
    cols=c(cols, "red") 
  }
}

test_qgraph<-qgraph(test, 
                    edge.labels=TRUE, 
                    edge.label.color='black', 
                    #mode = "direct", 
                    edge.color=cols,
                    esize=10)








3 = transplant incarnata biased
5 = transplant fucshii biased

4 = native fuchsii biased 
6 = native incarnata biased


clust1<-h1h2_all %>% filter(cluster == "cluster1") %>% rownames()
clust2<-h1h2_all %>% filter(cluster == "cluster2") %>% rownames()
clust3<-h1h2_all %>% filter(cluster == "cluster3") %>% rownames()
clust4<-h1h2_all %>% filter(cluster == "cluster4") %>% rownames()
clust5<-h1h2_all %>% filter(cluster == "cluster5") %>% rownames()
clust6<-h1h2_all %>% filter(cluster == "cluster6") %>% rownames()
clust7<-h1h2_all %>% filter(cluster == "cluster7") %>% rownames()








clust1_enrich<-get_enriched_terms(clust1, mp) 
clust2_enrich<-get_enriched_terms(clust2, mp) 
clust3_enrich<-get_enriched_terms(clust3, mp) 
clust4_enrich<-get_enriched_terms(clust4, mp) 
clust5_enrich<-get_enriched_terms(clust5, mp) 
clust6_enrich<-get_enriched_terms(clust6, mp)
clust7_enrich<-get_enriched_terms(clust7, mp) 



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

#cluster 3 and 5



test2<-root_sargasso_dds_norm[test,] %>% dplyr::select(contains("fuchsii")) %>% rownames_to_column(var="gene") %>% melt() %>% dplyr::select(-c("variable")) %>% mutate(species="fucshii")
test3<-root_sargasso_dds_norm[test,] %>% dplyr::select(contains("incarnata")) %>% rownames_to_column(var="gene") %>% melt() %>% dplyr::select(-c("variable")) %>% mutate(species="incarnata")
test4<-rbind(test2, test3)

ggplot(data = test4, aes(x=gene, y=log(value))) + geom_boxplot(aes(fill=species))

plot_homeolog_heatmap<-function(cluster_name){
  test<-h1h2_all %>% filter(cluster == cluster_name) %>% rownames()
  test2<-root_sargasso_dds_norm[test,] %>% dplyr::select(contains("fuchsii"))
  test3<-root_sargasso_dds_norm[test,] %>% dplyr::select(contains("incarnata"))
  test4<-cbind(test2, test3)
  pheatmap(test4, 
           #scale = "row", 
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           treeheight_row = 0, treeheight_col = 0,
  )
}


3 = transplant incarnata biased
5 = transplant fucshii biased

4 = native fuchsii biased 
6 = native incarnata biased

plot_homeolog_heatmap("cluster6")













