
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(magrittr)
library(DESeq2)


####################################################################
#      keep a record of species and source iof ITS sequence        #
####################################################################

# Ceratobasidium: https://unite.ut.ee/bl_forw_sh.php?sh_name=SH0110961.10FU#fndtn-panel3
# Tulasnella: https://unite.ut.ee/bl_forw_sh.php?sh_name=SH1010020.10FU#fndtn-panel1
# Serendipita: https://unite.ut.ee/bl_forw_sh.php?sh_name=SH0188460.10FU#fndtn-panel3
# Sebacina: https://unite.ut.ee/bl_forw_sh.php?sh_name=SH0136516.10FU#fndtn-panel1
# Russula: https://unite.ut.ee/bl_forw_sh.php?sh_name=SH0197565.10FU#fndtn-panel1


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


######################################################
#                    read in counts                  #
######################################################

strings_to_remove_root<-c("root_samples/", "_mycorrhizal.sorted.bam")
mycorrhizal_counts<-read_in_featurecounts("mycorrhizal_sequences_featurecounts", strings_to_remove_root)


#######################################################
#   normalise counts by ITS length and library size   #
#######################################################

library_size<-read_delim("library_size.txt", col_names = FALSE) %>% 
  set_colnames(c("sample", "library_size")) %>% 
  mutate(sample=str_replace(sample, "r.1", "r")) %>%
  mutate(scaling_factor=library_size/1000000)


# divide counts by per-million scaling factor
mycorrhizal_counts$libnorm_counts<-apply(mycorrhizal_counts$counts, 1, function(x) x/library_size$scaling_factor) %>% t()

# divide libsize normalised counts by gene length to get fpkm
mycorrhizal_counts$libnorm_counts<-apply(mycorrhizal_counts$libnorm_counts, 2, function(x) x/(mycorrhizal_counts$lengths/1000))


######################################################
#      get root sample df for annotating heatmap     #
######################################################

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

root_samples<-read.table("Root_samples.txt", header=FALSE, col.names = c("sample_id")) %>% 
  mutate(assay=case_when(endsWith(sample_id, "_s") == TRUE ~ "sRNA",
                         endsWith(sample_id, "_s") == FALSE ~ "RNAseq")) %>%
  mutate(tissue=case_when(substr(sample_id,5,5) == "l" ~ "leaf",
                          substr(sample_id,5,5) == "r" ~ "root")) %>%
  filter(tissue == "root" & assay == "RNAseq") %>% dplyr::select(sample_id)

root_samples<-annotate_samples(root_samples) %>% dplyr::select(c("species", "treatment"))


####################################################
#         make  histogram of read counts           #
####################################################

pdf("SupplementaryFigure3_mycorrhizal_reads_histogram.pdf", height=15, width=15)
melt(mycorrhizal_counts$libnorm_counts) %>%
  set_colnames(c("Taxon", "sample_id", "reads_mapped")) %>%
  mutate(species=case_when(substr(sample_id,1,1) == "m" ~ "majalis", substr(sample_id,1,1) == "t" ~ "traunsteineri")) %>%
  mutate(locality=case_when(substr(sample_id,3,3) == "S" ~ "St Ulrich", substr(sample_id,3,3) == "K" ~ "Kitzbuhl")) %>%
  mutate(environment=environment<-case_when(substr(sample_id,2,2) == "M" ~ "majalis", substr(sample_id,2,2) == "T" ~ "traunsteineri")) %>%
  mutate(treatment=case_when(species == environment ~ "Native", species != environment ~ "Transplant")) %>%
  mutate(Taxon=str_split_i(Taxon, "_", 1)) %>%
  ggplot(aes(x=Taxon, y=log(reads_mapped), fill=species), colour=black) + 
  geom_boxplot(lwd=0.7) +
  scale_fill_manual(values=c("gold", "dodgerblue")) +
  scale_colour_manual(values=c("dodgerblue", "gold")) +
  facet_wrap(~treatment, ncol=1) +
  ylab("log(FPKM)") +
  theme(legend.title = element_blank(),
        plot.background = element_rect(color = "black"),
        text = element_text(size = 45),
        axis.text=element_text(size=35),
        axis.title.y=element_text(size=45),
        axis.title.x=element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "grey", fill=NA, size=2),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  guides(fill="none")
dev.off()

####################################################
#         make heatmap of read counts           #
####################################################



heatmap_data<-mycorrhizal_counts$libnorm_counts
rownames(heatmap_data) <- str_split_i(rownames(heatmap_data), "_", 1)
ann_colors = list(treatment = c(native="black", transplant="white"),
                  species=c(majalis="gold", traunsteineri="dodgerblue"))
pdf("SupplementaryFigure4_mycorrhizal_reads_heatmap.pdf", height=10, width=15)
pheatmap::pheatmap(heatmap_data, 
                   scale = "row", 
                   cluster_cols = FALSE,
                   cluster_rows = FALSE,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   treeheight_row = 0, 
                   treeheight_col = 0,
                   annotation_col = root_samples,
                   gaps_col=c(9, 18, 27),
                   fontsize = 22,
                   annotation_colors=ann_colors)
dev.off()




