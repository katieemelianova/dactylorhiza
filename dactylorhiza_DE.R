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
source("dactylorhiza_functions.R")


#############################################################################
#        load in the GO ID mappings to transcripts for use later on         #
#############################################################################

# for the GO term enrichment tests
mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")

names(mp) <- names(mp) %>% 
  str_remove("-RA") %>%
  str_remove("-RB") %>%
  str_remove("-RC")

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
samples_to_exclude<-c("mTS2", "mTS2r", 
                      "tTS4", "tTS4r",
                      "tTK1", "tTK1r")

###########################################################
#          annotate samples using above function          #
###########################################################
leaf_samples<-annotate_samples(leaf_samples)
root_samples<-annotate_samples(root_samples)

leaf_samples <- leaf_samples[!(row.names(leaf_samples) %in% samples_to_exclude),]
root_samples <- root_samples[!(row.names(root_samples) %in% samples_to_exclude),]


####################################################################################################
#  define unwanted strings and read in featurecounts, joining and removing unwanted string         #
####################################################################################################

strings_to_remove_root<-c("root_samples/", "Aligned.sortedByCoord.out.bam")
df_root<-read_in_featurecounts('dactylorhiza_root_featurecounts', strings_to_remove_root)
df_root$counts %<>% dplyr::select(-contains(samples_to_exclude))
df_counts_root<-df_root$counts
df_lengths_root<-df_root$lengths

strings_to_remove_leaf<-c("StUlrich/", "Kitzbuhel/", "Aligned.sortedByCoord.out.bam")
df_leaf<-read_in_featurecounts('dactylorhiza_leaf_featurecounts', strings_to_remove_leaf)
df_leaf$counts %<>% dplyr::select(-contains(samples_to_exclude))
df_counts_leaf<-df_leaf$counts
df_lengths_leaf<-df_leaf$lengths

######################################
#      Figure 2 Draw PCA plots       #
######################################

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
root_vst<-root_vst[,-c(28)]

# plot the PCA 
pcaRootData <-plotPCA(root_vst, intgroup=c("treatment", "species", "locality"), ntop = 1000, returnData = TRUE)
percentRootVar <- round(100 * attr(pcaRootData, "percentVar"))

pcaLeafData <-plotPCA(leaf_vst, intgroup=c("treatment", "species", "locality"), ntop = 1000, returnData = TRUE)
percentLeafVar <- round(100 * attr(pcaLeafData, "percentVar"))

pcaRootData %<>% mutate("Tissue" = "Root")
pcaLeafData %<>% mutate("Tissue" = "Leaf")
pcaData <- rbind(pcaRootData, pcaLeafData)

pca_plot<-ggplot(pcaData, aes(PC1, PC2, color=locality, fill=locality, shape=interaction(species, treatment))) +
  geom_point(size=7, stroke = 1.5) +
  #xlab(paste0("PC1: ",percentRootVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentRootVar[2],"% variance")) + 
  #coord_fixed() + 
  theme(legend.title = element_blank(),
        text = element_text(size = 20), 
        legend.text=element_text(size=25),
        axis.text=element_text(size=24),
        axis.title=element_text(size=28),
        plot.margin = unit(c(0,0,0,0), "cm"),
        strip.text.x = element_text(size = 27),) +
  scale_shape_manual(values = c(16, 15, 10, 7)) +
  facet_wrap(~ Tissue, ncol=1)

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2.png", width = 1000, height = 1000)
pca_plot
dev.off()

percentRootVar[1]
percentRootVar[2]
percentLeafVar[1]
percentLeafVar[2]


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

pvalue_threshold <- 0.05
logfc_threshold<-2

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
                                         #(majalis_env > 2 | majalis_env < -2) &
                                         #(traunst_env < 2 | traunst_env > -2) &
                                         traunst_env_padj > 0.05 ~ "DE in D. majalis environment only",
                                         #(majalis_env < 2 | majalis_env > -2) &
                                         #(traunst_env > 2 | traunst_env < -2) &
                                         majalis_env_padj > 0.05 &
                                         traunst_env_padj < 0.05 ~ "DE in D. traunsteineri environment only"))


# designate all unlabelled status values Not Significant         
all_bound %<>% replace_na(list(status = "Not significant"))


# arrange the data so that the DE points are plotted last and on top of the grey non significant points
all_bound %<>% arrange(factor(status, levels = c("Not significant", 
                                                 "DE in D. majalis environment only", 
                                                 "DE in D. traunsteineri environment only", 
                                                 "Constitutively DE")))



all_bound %>% dplyr::select("status", "locality", "tissue") %>% table()


all_bound %>% filter(status == "DE in D. majalis environment only" & traunst_env > 0 & tissue == "Root") %>% nrow()
all_bound %>% filter(status == "DE in D. majalis environment only" & traunst_env > 0 & tissue == "Leaf") %>% nrow()

all_bound %>% filter(status == "DE in D. majalis environment only" & traunst_env < 0 & tissue == "Root") %>% nrow()
all_bound %>% filter(status == "DE in D. majalis environment only" & traunst_env < 0 & tissue == "Leaf") %>% nrow()


all_bound %>% filter(status == "DE in D. traunsteineri environment only" & traunst_env > 0) %>% nrow()
all_bound %>% filter(status == "DE in D. traunsteineri environment only" & traunst_env < 0) %>% nrow()

test <- ggplot(all_bound, aes(x=majalis_env, y=traunst_env, colour=status)) +
  geom_point(alpha=0.5, size = 4)

ggplot(all_bound, aes(x=status, fill=status)) +
  geom_bar(position="stack") + facet_wrap(~ tissue + locality, ncol=2)

test<-all_bound %>% dplyr::select(status, locality, tissue) %>% filter(status != "Not significant") %>% melt()
ggplot(test, aes(x=locality, fill=status)) +
  geom_bar(position="stack") + facet_wrap(~ tissue, ncol=2)

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure3.png", width = 1300, height = 1000)
ggplot(all_bound, aes(x=majalis_env, y=traunst_env, colour=status)) +
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
        strip.text.x = element_text(size = 27)) +
  guides(colour = guide_legend(override.aes = list(size=10)))
dev.off()


# how many DE genes are there in each category?
all_bound %>% dplyr::select("status") %>% pull() %>% table()

# how many genes are upregulated in traunsteineri relative to majalis?
all_bound %>% filter(majalis_env > 2 & traunst_env > 2 & majalis_env_padj < 0.05 & traunst_env_padj < 0.05)

# how many genes are upregulated in majalis relative to traunsteineri?
all_bound %>% filter(majalis_env > 2 & traunst_env > 2 & majalis_env_padj < 0.05 & traunst_env_padj < 0.05)


#a<-get_significant_genes(transplant_majalis_kitzbuhl_leaf)
#b<-get_significant_genes(transplant_majalis_stulrich_leaf)
#c<-get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf)
#d<-get_significant_genes(transplant_traunsteineri_stulrich_leaf)
#
#intersect(a, c)
#
#mylist=list(mkL=a,
#            msL=b,
#            tkL=c,
#            tsL=d,
#            mkR=e,
#            msR=f,
#            tkR=g,
#            tsR=h,
#            constitutive=const)
#
#
#upset(fromList(mylist), order.by = "freq", nsets = 10)


#e<-get_significant_genes(transplant_majalis_kitzbuhl_root)
#f<-get_significant_genes(transplant_majalis_stulrich_root)
#g<-get_significant_genes(transplant_traunsteineri_kitzbuhl_root)
#h<-get_significant_genes(transplant_traunsteineri_stulrich_root)



#const_traunst_up <- all_bound %>% filter(status == "Constitutively DE" & majalis_env < 0 & traunst_env < 0) %>% dplyr::select(gene_id) %>% pull()
#const_majalis_up <- all_bound %>% filter(status == "Constitutively DE" & majalis_env > 0 & traunst_env > 0) %>% dplyr::select(gene_id) %>% pull()
#const <- all_bound %>% filter(status == "Constitutively DE") %>% dplyr::select(gene_id) %>% pull()


#const %<>% str_remove(".t1:cds") %>% str_remove("cds") %>% str_remove("-RA:") %>% str_remove(".t2")
#const_traunst_up %<>% str_remove(".t1:cds") %>% str_remove("cds") %>% str_remove("-RA:") %>% str_remove(".t2")
#const_majalis_up %<>% str_remove(".t1:cds") %>% str_remove("cds") %>% str_remove("-RA:") %>% str_remove(".t2")
#
#consitutive_trauns_up_goenrich<-get_enriched_terms(const_traunst_up, mp) 
#consitutive_majalis_up_goenrich<-get_enriched_terms(const_majalis_up, mp) 
#consitutive_goenrich<-get_enriched_terms(const, mp) 
#
#
#consitutive_trauns_up_goenrich %>% filter(Significant > 1)
#consitutive_majalis_up_goenrich %>% filter(Significant > 1)
#consitutive_goenrich %>% filter(Significant > 1 & classicFisher < 0.05)
#
#intersect(a, const)
#intersect(a, const)



#############################################
#    Table 1 Constitutively DEG GO terms    #
#############################################

godb_table<-toTable(GOTERM)

constitutive_annotation<-data.frame(Term=c(),
                                    Ontology=c(),
                                    expression_pattern=c(),
                                    gene_id=c(),
                                    comparison=c())

for (i in 1:nrow(gene_ids_constitutive_traunst_up)){
  go_ids<-mp[gene_ids_constitutive_traunst_up[i,]$gene_id][[1]]
  if (length(go_ids) >=1){
    print(gene_ids_constitutive_traunst_up[i,])
    go<-godb_table[godb_table$go_id %in% go_ids,]
    
    new_rows<-data.frame(go[,c("Term", "Ontology")] %>% unique())
    new_rows$expression_pattern<-"T > M"
    new_rows$gene_id=gene_ids_constitutive_traunst_up[i,]$gene_id
    new_rows$comparison=gene_ids_constitutive_traunst_up[i,]$comparison
    constitutive_annotation <- rbind(constitutive_annotation, new_rows)
  }
}

for (i in 1:nrow(gene_ids_constitutive_majalis_up)){
  go_ids<-mp[gene_ids_constitutive_majalis_up[i,]$gene_id][[1]]
  if (length(go_ids) >=1){
    print(gene_ids_constitutive_majalis_up[i,])
    go<-godb_table[godb_table$go_id %in% go_ids,]
    new_rows<-data.frame(go[,c("Term", "Ontology")] %>% unique())
    new_rows$expression_pattern<-"M > T"
    new_rows$gene_id=gene_ids_constitutive_majalis_up[i,]$gene_id
    new_rows$comparison=gene_ids_constitutive_majalis_up[i,]$comparison
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


######################################################
#    Figure 4 Venn Diagram of number of DE genes     #
######################################################


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
  filename = 'Figure4a.png',
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
  filename = 'Figure4b.png',
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



#############################################
#       Figure 5 DE gene count plots        #
#############################################

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

colnames(de_counts)<-c("Species", "Tissue", "Locality", "Individual", "Direction", "Number of genes")

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure5.png", width = 900, height = 1200)
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
  

de_counts %>% filter(Tissue == "Leaf")
de_counts %>% filter(Tissue == "Leaf" & Locality == "St. Ulrich")

de_counts %>% filter(Tissue == "Root")



#############################################
#            GO term enrichment             #
#############################################



#######################################
#        majalis leaf kitzbuhl        #
#######################################
transplant_majalis_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 

#######################################
#        majalis leaf st ulrich        #
#######################################
transplant_majalis_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 



#######################################
#        majalis root kitzbuhl        #
#######################################
transplant_majalis_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$down, mp) 

#######################################
#        majalis root st ulrich       #
#######################################
transplant_majalis_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_majalis_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE, mappings_format = TRUE)$down, mp) 

#######################################
#     traunsteineri leaf kitzbuhl     #
#######################################
transplant_traunsteineri_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 


#######################################
#     traunsteineri leaf st ulrich    #
#######################################
transplant_traunsteineri_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE, mappings_format = TRUE)$down, mp) 

#######################################
#     traunsteineri root kitzbuhl     #
#######################################
transplant_traunsteineri_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE, mappings_format = TRUE)$down, mp) 

#######################################
#     traunsteineri root st ulrich    #
#######################################
transplant_traunsteineri_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE, mappings_format = TRUE)$up, mp) 
transplant_traunsteineri_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE, mappings_format = TRUE)$down, mp) 




#######################################
#        Plot GO term enrichment      #
#######################################


prepare_go_df<-function(topgo_object){
  topgo_object %<>% 
    mutate(`Rich score`=Significant/Annotated) %>% 
    dplyr::select(Term, GO.ID, Significant, Annotated, classicFisher, `Rich score`) %>% 
    #filter(classicFisher < 0.05) %>% 
    filter(classicFisher < 0.05 & `Rich score` >= 0.01) %>% 
    head(10) %>% 
    data.frame()
  return(topgo_object)
}



mKLeafUp<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_up) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="up")
mSLeafUp<-prepare_go_df(transplant_majalis_stulrich_leaf_up) %>% mutate(comparison="D. majalis St. Ulrich", Direction="up")
tKLeafUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_up) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="up")
tSLeafUp<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_up) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="up")
mKLeafDown<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_down) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="down")
mSLeafDown<-prepare_go_df(transplant_majalis_stulrich_leaf_down) %>% mutate(comparison="D. majalis St. Ulrich", Direction="down")
tKLeafDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_down) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="down")
tSLeafDown<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_down) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="down")

mKRootUp<-prepare_go_df(transplant_majalis_kitzbuhl_root_up) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="up")
mSRootUp<-prepare_go_df(transplant_majalis_stulrich_root_up) %>% mutate(comparison="D. majalis St. Ulrich", Direction="up")
#tKRootUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up) %>% mutate(comparison="D.t.Kitz", Direction="up")
tKRootUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="up")

tSRootUp<-prepare_go_df(transplant_traunsteineri_stulrich_root_up) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="up")
mKRootDown<-prepare_go_df(transplant_majalis_kitzbuhl_root_down) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="down")
mSRootDown<-prepare_go_df(transplant_majalis_stulrich_root_down) %>% mutate(comparison="D. majalis St. Ulrich", Direction="down")
#tKRootDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down) %>% mutate(comparison="D.t.Kitz", Direction="down")
tKRootDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="down")
tSRootDown<-prepare_go_df(transplant_traunsteineri_stulrich_root_down) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="down")


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

colnames(leaf_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Rich factor", "comparison", "Direction")
colnames(root_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Rich factor", "comparison", "Direction")

root_go_bound %>% filter(comparison == "D. majalis St. Ulrich" & `Rich factor` > 0.01 & classicFisher < 0.05)
root_go_bound %>% filter(comparison == "D. traunsteineri Kitzbuhl" & `Rich factor` > 0.0001 & classicFisher < 0.05)
root_go_bound %>% filter(comparison == "D. majalis Kitzbuhl" & classicFisher < 0.05)

# I am adding a new column with a placeholder the same across all rows
# I will use this column as my X axis for the GO term plots
# if I use the comparison column, it staggers the plots, and it looks better if the points are stacked
leaf_go_bound_newcol <- leaf_go_bound %>% mutate(newcol="placeholder", tissue="Leaf")
root_go_bound_newcol <- root_go_bound %>% mutate(newcol="placeholder", tissue="Root")


a<-ggplot(leaf_go_bound_newcol, aes(x=newcol, y=Term, color = Direction, size=`Rich factor`)) + 
  geom_point() + facet_grid(rows=vars(comparison), scales="free", space= "free", switch = "y", cols=vars(tissue)) + 
  theme(text = element_text(size = 60), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 49),
        legend.text=element_text(size=36),
        legend.title=element_text(size=40),
        panel.spacing=unit(1, "lines")) + 
  scale_size_continuous(range = c(9, 18)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("blue", "red")) +
  theme(legend.position = "none")

# I had to shorten the name of one of the facet grids (D.t. Kitz due to soace
# this alphabetically rearranges the orderr so locking in order of localitoies this way
root_go_bound_newcol$comparison <- factor(root_go_bound$comparison, levels = unique(root_go_bound$comparison))


b<-ggplot(root_go_bound_newcol, aes(x=newcol, y=Term, color = Direction, size=`Rich factor`)) + 
  geom_point() + facet_grid(rows=vars(comparison), scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 60), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 49),
        legend.text=element_text(size=36),
        legend.title=element_text(size=40),
        panel.spacing=unit(1, "lines")) + 
  scale_size_continuous(range = c(9, 18)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("blue", "red")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none") 

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure6.png", height=3000, width=4200)
egg::ggarrange(a, b, ncol=2)
dev.off()


########################################################
#      plot heatmaps of reciprocally plastic genes     #
########################################################

st_ulrich_root_effect_of_environment<-specify_comparison(root_samples, df_counts_root, 'locality == "St Ulrich"') %>% run_diffexp("environment", df_root$lengths)
kitzbuhl_root_effect_of_environment<-specify_comparison(root_samples, df_counts_root, 'locality == "Kitzbuhl"') %>% run_diffexp("environment", df_root$lengths)
st_ulrich_leaf_effect_of_environment<-specify_comparison(leaf_samples, df_counts_leaf, 'locality == "St Ulrich"') %>% run_diffexp("environment", df_leaf$lengths)
kitzbuhl_leaf_effect_of_environment<-specify_comparison(leaf_samples, df_counts_leaf, 'locality == "Kitzbuhl"') %>% run_diffexp("environment", df_leaf$lengths)

normalise_data_for_heatmap<-function(dds){
  de_genes<-dds %>% results() %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.0005) %>% rownames() 
  gene_counts<-counts(dds, normalized=TRUE)
  de_gene_counts<-gene_counts[de_genes,]
  new_column_order<-de_gene_counts %>% colnames %>% sort()
  de_gene_counts <- de_gene_counts %>% data.frame() %>% dplyr::select(new_column_order)
  heatmap_data <- t(apply(de_gene_counts, 1, function(x) x/mean(x)))
  #heatmap_data <- t(apply(de_gene_counts, 1, function(x) x/sd(x)))
  # log normalise, remove NAs
  heatmap_data <- log2(heatmap_data)
  # get data range for breaks
  max_data <- max(heatmap_data, na.rm = TRUE)
  min_data <- -min(heatmap_data, na.rm = TRUE)
  range <- min(max_data, min_data)
  heatmap_data[is.infinite(heatmap_data)] <- NA
  heatmap_data[is.nan(heatmap_data)] <- NA
  heatmap_data %<>% data.frame() %>% drop_na()
  return(heatmap_data)
}

st_ulrich_root_effect_of_environment_heatmap<-normalise_data_for_heatmap(st_ulrich_root_effect_of_environment$dds)
kitzbuhl_root_effect_of_environment_heatmap<-normalise_data_for_heatmap(kitzbuhl_root_effect_of_environment$dds)
st_ulrich_leaf_effect_of_environment_heatmap<-normalise_data_for_heatmap(st_ulrich_leaf_effect_of_environment$dds)
kitzbuhl_leaf_effect_of_environment_heatmap<-normalise_data_for_heatmap(kitzbuhl_leaf_effect_of_environment$dds)



# theres an "r" on the end of the root samples so need to remove these to be able to rbind by column with leaf
colnames(st_ulrich_root_effect_of_environment_heatmap) %<>% str_remove("r")
colnames(kitzbuhl_root_effect_of_environment_heatmap) %<>% str_remove("r")


# creatig the grouping to enable me to make the annotation DF for the grouped heatmap
st_ulrich_leaf_effect_of_environment_heatmap %<>% mutate(grouping="Leaf")
st_ulrich_root_effect_of_environment_heatmap %<>% mutate(grouping="Root")
kitzbuhl_leaf_effect_of_environment_heatmap %<>% mutate(grouping="Leaf")
kitzbuhl_root_effect_of_environment_heatmap %<>% mutate(grouping="Root")

# make stUlrich annotation
st_ulrich_annotation<-data.frame(c(st_ulrich_leaf_effect_of_environment_heatmap$grouping, 
                                   st_ulrich_root_effect_of_environment_heatmap$grouping))

colnames(st_ulrich_annotation) <-"Tissue"
# use make.names to rename duplicate genes and not break data.frame function
rownames(st_ulrich_annotation) <- 
  make.names(c(rownames(st_ulrich_leaf_effect_of_environment_heatmap),
               rownames(st_ulrich_root_effect_of_environment_heatmap)), 
             unique = TRUE)


# make Kitzbuhl annotation
kitzbuhl_annotation<-data.frame(c(kitzbuhl_leaf_effect_of_environment_heatmap$grouping, 
                                   kitzbuhl_root_effect_of_environment_heatmap$grouping))

colnames(kitzbuhl_annotation) <-"Tissue"
# use make.names to rename duplicate genes and not break data.frame function
rownames(kitzbuhl_annotation) <- 
  make.names(c(rownames(kitzbuhl_leaf_effect_of_environment_heatmap), 
               rownames(kitzbuhl_root_effect_of_environment_heatmap)), 
             unique = TRUE)

# remove the grouping column you made to get the annotation dataframe
st_ulrich_leaf_effect_of_environment_heatmap %<>% dplyr::select(-grouping) %>% set_rownames(st_ulrich_annotation %>% filter(Tissue == "Leaf") %>% rownames())
st_ulrich_root_effect_of_environment_heatmap %<>% dplyr::select(-grouping) %>% set_rownames(st_ulrich_annotation %>% filter(Tissue == "Root") %>% rownames())
kitzbuhl_leaf_effect_of_environment_heatmap  %<>% dplyr::select(-grouping) %>% set_rownames(kitzbuhl_annotation %>% filter(Tissue == "Leaf") %>% rownames())
kitzbuhl_root_effect_of_environment_heatmap  %<>% dplyr::select(-grouping) %>% set_rownames(kitzbuhl_annotation %>% filter(Tissue == "Root") %>% rownames())


st_ulrich_effect_of_environment_heatmap<-rbind(st_ulrich_leaf_effect_of_environment_heatmap, st_ulrich_root_effect_of_environment_heatmap)
kitzbuhl_effect_of_environment_heatmap<-rbind(kitzbuhl_leaf_effect_of_environment_heatmap, kitzbuhl_root_effect_of_environment_heatmap)

stulrich_annotation_col = data.frame(Species = c(rep("D. majalis", 7), rep("D. traunsteineri", 7)), Environment=c(rep("M", 4), rep("T", 3), rep("M", 4), rep("T", 3)))
kitzbuhl_annotation_col = data.frame(Species = c(rep("D. majalis", 10), rep("D. traunsteineri", 9)), Environment=c(rep("M", 5), rep("T", 5), rep("M", 5), rep("T", 4)))


rownames(stulrich_annotation_col)<-st_ulrich_effect_of_environment_heatmap %>% colnames()  
rownames(kitzbuhl_annotation_col)<-kitzbuhl_effect_of_environment_heatmap %>% colnames()  

# define colours of annotation to be used
ann_colors = list(Environment = c(M="deepskyblue", T="goldenrod1"),
                  Species = c(`D. majalis`="lightpink", `D. traunsteineri`="lightgreen"))

nrow(st_ulrich_annotation)
nrow(kitzbuhl_annotation)
kitzbuhl_annotation$Tissue %>% table()
st_ulrich_annotation$Tissue %>% table()


png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure7a.png", height=2000, width=2300)
pheatmap.type(Data=st_ulrich_effect_of_environment_heatmap, 
              annRow=st_ulrich_annotation, 
              show_rownames=FALSE,
              cluster_cols=FALSE,
              treeheight_row = 0, treeheight_col = 0,
              show_colnames = F, scale="none",
              border_color = NA,
              annotation_col=stulrich_annotation_col,
              annotation_colors=ann_colors,
              annotation_names_row=F,
              annotation_names_col=F,
              fontsize = 45,
              legend=FALSE,
              gaps_row=c(19),
              gaps_col=c(4,8,12))
dev.off()



png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure7b.png", height=2000, width=2300)
pheatmap.type(Data=kitzbuhl_effect_of_environment_heatmap, 
              annRow=kitzbuhl_annotation, 
              show_rownames=FALSE,
              cluster_cols=FALSE,
              treeheight_row = 0, treeheight_col = 0,
              show_colnames = F, scale="none",
              border_color = NA,
              annotation_col=kitzbuhl_annotation_col,
              annotation_colors=ann_colors,
              annotation_names_row=F,
              annotation_names_col=F,
              fontsize = 45,
              legend=FALSE,
              gaps_row=c(6),
              gaps_col=c(5,10,15))
dev.off()




#######################################################################################################
#      Figure 9 plot venn diagram of effect of transplant (both species) and shared plastic genes     #
#######################################################################################################


shared_plastic_stulrich<-st_ulrich_effect_of_environment_heatmap %>% rownames()
shared_plastic_kitzbuhl<-kitzbuhl_effect_of_environment_heatmap %>% rownames()

stulrich_transplant_genes<-c(get_significant_genes(transplant_traunsteineri_stulrich_leaf),
                             get_significant_genes(transplant_traunsteineri_stulrich_root),
                             get_significant_genes(transplant_majalis_stulrich_leaf),
                             get_significant_genes(transplant_majalis_stulrich_root)) %>% unique()

kitzbuhl_transplant_genes<-c(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf),
                             get_significant_genes(transplant_traunsteineri_kitzbuhl_root),
                             get_significant_genes(transplant_majalis_kitzbuhl_leaf),
                             get_significant_genes(transplant_majalis_kitzbuhl_root)) %>% unique()



length(kitzbuhl_transplant_genes)
kitzbuhl_effect_of_environment_heatmap %>% rownames() %>% length()
intersect(kitzbuhl_transplant_genes, kitzbuhl_effect_of_environment_heatmap %>% rownames()) %>% length()

length(stulrich_transplant_genes)
st_ulrich_effect_of_environment_heatmap %>% rownames() %>% length()
intersect(stulrich_transplant_genes, st_ulrich_effect_of_environment_heatmap %>% rownames()) %>% length()



stu_intersect_euler <- intersect(stulrich_transplant_genes, shared_plastic_stulrich) %>% length()
stu_trans_euler <- (stulrich_transplant_genes %>% length()) - stu_intersect_euler
stu_plastic_euler <- (shared_plastic_stulrich %>% length()) - stu_intersect_euler

ktz_intersect_euler <- intersect(kitzbuhl_transplant_genes, shared_plastic_kitzbuhl) %>% length()
ktz_trans_euler <- (kitzbuhl_transplant_genes %>% length()) - ktz_intersect_euler
ktz_plastic_euler <- (shared_plastic_kitzbuhl %>% length()) - ktz_intersect_euler


set.seed(1)
stulrich_euler <- c("St Ulrich transplant" = stu_trans_euler,
                    "St Ulrich shared plastic" = stu_plastic_euler, 
                    "St Ulrich transplant&St Ulrich shared plastic" = stu_intersect_euler)

kitzbuhl_euler <- c("Kitzbuhl transplant" = ktz_trans_euler,
                    "Kitzbuhl shared plastic" = ktz_plastic_euler, 
                    "Kitzbuhl transplant&Kitzbuhl shared plastic"  = ktz_intersect_euler)



p1 <- plot(euler(stulrich_euler), 
           #quantities = TRUE,
           fills = c("maroon2", "tan1"),
           quantities = list(cex = 4),
           legend=list(fontsize=28))

p2 <- plot(euler(kitzbuhl_euler), 
           #quantities = TRUE,
           quantities = list(cex = 4),
           fills = c("seagreen3", "dodgerblue"), 
           legend=list(fontsize=28))

png(file="Figure8.png", width = 1200, height = 1200)
gridExtra::grid.arrange(p1, p2)
dev.off()



one<-transplant_traunsteineri_stulrich_leaf_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()
two<-transplant_traunsteineri_stulrich_leaf_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()
three<-transplant_traunsteineri_stulrich_root_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()
four<-transplant_traunsteineri_stulrich_root_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()

five<-transplant_majalis_stulrich_leaf_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()
six<-transplant_majalis_stulrich_leaf_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()
seven<-transplant_majalis_stulrich_root_up %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()
eight<-transplant_majalis_stulrich_root_down %>% filter(as.numeric(classicFisher) < 0.05) %>% dplyr::select(GO.ID) %>% pull()


intersect(test_go, c(one, two, three, four, five, six, seven, eight)) %>% length()

test_go %>% length()

####### getting annotations of the plastic reciprocal genes

gene_to_go_table<-function(genes){
  go_id_list<-c()
  for (gene in genes){
    go_ids<-mp[gene]
    if (length(go_ids) >=1){
      go_id_list <-c(go_id_list, go_ids[[1]])
    }
  }
  go_id_list<-unique(go_id_list)
  return(go_id_list)
}

st_ulrich_leaf_effect_of_environment_go<-st_ulrich_leaf_effect_of_environment_heatmap %>% rownames() %>% gene_to_go_table()
st_ulrich_root_effect_of_environment_go<-st_ulrich_root_effect_of_environment_heatmap %>% rownames() %>% gene_to_go_table()
kitzbuhl_leaf_effect_of_environment_go<-kitzbuhl_leaf_effect_of_environment_heatmap %>% rownames() %>% gene_to_go_table()
kitzbuhl_root_effect_of_environment_go<-kitzbuhl_root_effect_of_environment_heatmap %>% rownames() %>% gene_to_go_table()




transplant_majalis_kitzbuhl_leaf_genes<-get_significant_genes(transplant_majalis_kitzbuhl_leaf)
transplant_majalis_stulrich_leaf_genes<-get_significant_genes(transplant_majalis_stulrich_leaf)
transplant_traunsteineri_kitzbuhl_leaf_genes<-get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf)
transplant_traunsteineri_stulrich_leaf_genes<-get_significant_genes(transplant_traunsteineri_stulrich_leaf)


st_ulrich_leaf_effect_of_environment_genes<-get_significant_genes(st_ulrich_leaf_effect_of_environment)
st_ulrich_root_effect_of_environment_genes<-get_significant_genes(st_ulrich_root_effect_of_environment)
kitzbuhl_leaf_effect_of_environment_genes<-get_significant_genes(kitzbuhl_leaf_effect_of_environment)
kitzbuhl_root_effect_of_environment_genes<-get_significant_genes(kitzbuhl_root_effect_of_environment)




listInputKitzbuhl<-list(transplant_majalis_kitzbuhl_leaf_genes=transplant_majalis_kitzbuhl_leaf_genes,
                transplant_traunsteineri_kitzbuhl_leaf_genes=transplant_traunsteineri_kitzbuhl_leaf_genes,
                kitzbuhl_leaf_effect_of_environment_genes=kitzbuhl_leaf_effect_of_environment_genes,
                kitzbuhl_root_effect_of_environment_genes=kitzbuhl_root_effect_of_environment_genes)

listInputStUlrich<-list(
                        transplant_majalis_stulrich_leaf_genes=transplant_majalis_stulrich_leaf_genes,
                        transplant_traunsteineri_stulrich_leaf_genes=transplant_traunsteineri_stulrich_leaf_genes,
                        st_ulrich_leaf_effect_of_environment_genes=st_ulrich_leaf_effect_of_environment_genes,
                        st_ulrich_root_effect_of_environment_genes=st_ulrich_root_effect_of_environment_genes)

listInputAll<-list(transplant_majalis_kitzbuhl_leaf_genes=transplant_majalis_kitzbuhl_leaf_genes,
                        transplant_majalis_stulrich_leaf_genes=transplant_majalis_stulrich_leaf_genes,
                        transplant_traunsteineri_kitzbuhl_leaf_genes=transplant_traunsteineri_kitzbuhl_leaf_genes,
                        transplant_traunsteineri_stulrich_leaf_genes=transplant_traunsteineri_stulrich_leaf_genes,
                        st_ulrich_leaf_effect_of_environment_genes=st_ulrich_leaf_effect_of_environment_genes,
                        st_ulrich_root_effect_of_environment_genes=st_ulrich_root_effect_of_environment_genes,
                        kitzbuhl_leaf_effect_of_environment_genes=kitzbuhl_leaf_effect_of_environment_genes,
                        kitzbuhl_root_effect_of_environment_genes=kitzbuhl_root_effect_of_environment_genes)


st_ulrich_leaf_effect_of_environment_go
st_ulrich_root_effect_of_environment_go
kitzbuhl_leaf_effect_of_environment_go
kitzbuhl_root_effect_of_environment_go
  
listInputGO<-list(st_ulrich_leaf_effect_of_environment_go=st_ulrich_leaf_effect_of_environment_go,
                  st_ulrich_root_effect_of_environment_go=st_ulrich_root_effect_of_environment_go,
                  kitzbuhl_leaf_effect_of_environment_go=kitzbuhl_leaf_effect_of_environment_go,
                  kitzbuhl_root_effect_of_environment_go=kitzbuhl_root_effect_of_environment_go)

  


upset(fromList(listInputAll), order.by = "freq", nsets = 10)



stu_root<-st_ulrich_root_effect_of_environment$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id") %>%
  mutate("locality" = "stulrich")

kitz_root<-kitzbuhl_root_effect_of_environment$results %>% 
  data.frame() %>% 
  dplyr::select(log2FoldChange, padj) %>% 
  rownames_to_column(var="gene_id") %>%
  mutate("locality" = "kitzzbul")

all_env<-rbind(stu_root, kitz_root)




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

