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
library(SuperExactTest)
library(ggvenn)
library(cowplot)
library(ggpubr)
library(GSEABase)
source("dactylorhiza_functions.R")


############################################################
#       make supplementary table of read count stats       #
############################################################

# the file here is generated throughb a simple script found in /scratch/botany/katie/orchid/collate_mapping_stats.sh

star_summary_stats<-read.table("rna_seq_summary_statistics.txt") %>% set_colnames(c("Library", "Input Reads", "Uniquely Mapped", "Percent Multimapped"))
colnames(star_summary_stats)
star_summary_stats$Library<-str_remove(star_summary_stats$Library, "Log.final.out")
star_summary_stats$`Percent Multimapped` <- str_remove(star_summary_stats$`Percent Multimapped`, "%") %>% as.numeric()
star_summary_stats$Library<-str_split_fixed(star_summary_stats$Library, "/", 3)[,3]

star_summary_stats %<>% mutate(Tissue=case_when(endsWith(star_summary_stats$Library, "r") == TRUE ~ "Root",
                                                !endsWith(star_summary_stats$Library, "r") == TRUE ~ "Leaf")) %>% 
  mutate(species=case_when(substr(Library,1,1) == "m" ~ "majalis", substr(Library,1,1) == "t" ~ "traunsteineri")) %>%
  mutate(locality=case_when(substr(Library,3,3) == "S" ~ "St Ulrich", substr(Library,3,3) == "K" ~ "Kitzbuhl")) 

  
star_summary_stats %>% summary()




aggregate(`Uniquely Mapped` ~ locality, data = star_summary_stats, mean)


#############################################################################
#        load in the GO ID mappings to transcripts for use later on         #
#############################################################################

# for the GO term enrichment tests
mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")

#names(mp) <- names(mp) %>% 
#  str_remove("-RA") %>%
#  str_remove("-RB") %>%
#  str_remove("-RC")

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


samples_to_exclude_leaf<-c("mTS2", 
                      "tTS4",
                      "tTK1",
                      "mTK1",
                      "tMK4",
                      "tMS4",
                      "tTK3")

samples_to_exclude_root<-c("mTS2r", 
                           "tTS4r",
                           "tTK1r",
                           "mTK2r", 
                           "tMK4r")



###########################################################
#          annotate samples using above function          #
###########################################################
leaf_samples<-annotate_samples(leaf_samples)
root_samples<-annotate_samples(root_samples)

leaf_samples <- leaf_samples[!(row.names(leaf_samples) %in% samples_to_exclude_leaf),]
root_samples <- root_samples[!(row.names(root_samples) %in% samples_to_exclude_root),]


####################################################################################################
#  define unwanted strings and read in featurecounts, joining and removing unwanted string         #
####################################################################################################

strings_to_remove_root<-c("root_samples/", "Aligned.sortedByCoord.out.bam")
df_root<-read_in_featurecounts('dactylorhiza_root_featurecounts', strings_to_remove_root)
df_root$counts %<>% dplyr::select(-contains(samples_to_exclude_root))
df_counts_root<-df_root$counts
df_lengths_root<-df_root$lengths

strings_to_remove_leaf<-c("StUlrich/", "Kitzbuhel/", "Aligned.sortedByCoord.out.bam")
df_leaf<-read_in_featurecounts('dactylorhiza_leaf_featurecounts', strings_to_remove_leaf)
df_leaf$counts %<>% dplyr::select(-contains(samples_to_exclude_leaf))
df_counts_leaf<-df_leaf$counts
df_lengths_leaf<-df_leaf$lengths


######################################
#      Figure 2 Draw PCA plots       #
######################################

# make a dds object from the total root samples (no subsetting)
root_kitzbuhel_dds<-specify_comparison(root_samples, df_counts_root, "locality == 'Kitzbuhl'")
root_stulrich_dds<-specify_comparison(root_samples, df_counts_root, "locality == 'St Ulrich'")


root_kitzbuhel_dds <- DESeqDataSetFromMatrix(countData = root_kitzbuhel_dds[["counts"]],
                                   colData = root_kitzbuhel_dds[["samples"]],
                                   design = ~ species + treatment)

root_stulrich_dds <- DESeqDataSetFromMatrix(countData = root_stulrich_dds[["counts"]],
                                             colData = root_stulrich_dds[["samples"]],
                                             design = ~ species + treatment)



# make a dds object from the total leaf samples
leaf_kitzbuhel_dds<-specify_comparison(leaf_samples, df_counts_leaf, "locality == 'Kitzbuhl'")
leaf_stulrich_dds<-specify_comparison(leaf_samples, df_counts_leaf, "locality == 'St Ulrich'")

leaf_kitzbuhel_dds <- DESeqDataSetFromMatrix(countData = leaf_kitzbuhel_dds[["counts"]],
                                   colData = leaf_kitzbuhel_dds[["samples"]],
                                   design = ~ species + treatment)

leaf_stulrich_dds <- DESeqDataSetFromMatrix(countData = leaf_stulrich_dds[["counts"]],
                                   colData = leaf_stulrich_dds[["samples"]],
                                   design = ~ species + treatment)


# perform variance stabilizig transformation
root_kitzbuhel_vst<-varianceStabilizingTransformation(root_kitzbuhel_dds)
root_stulrich_vst<-varianceStabilizingTransformation(root_stulrich_dds)
leaf_kitzbuhel_vst<-varianceStabilizingTransformation(leaf_kitzbuhel_dds)
leaf_stulrich_vst<-varianceStabilizingTransformation(leaf_stulrich_dds)


# plot the PCA 
pcaRootKitzbuhelData <-plotPCA(root_kitzbuhel_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentRootKitzbuhelVar <- round(100 * attr(pcaRootKitzbuhelData, "percentVar"))
pcaRootStUlrichData <-plotPCA(root_stulrich_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentRootStUlrichVar <- round(100 * attr(pcaRootStUlrichData, "percentVar"))

pcaLeafKitzbuhelData <-plotPCA(leaf_kitzbuhel_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentLeafKitzbuhelVar <- round(100 * attr(pcaLeafKitzbuhelData, "percentVar"))
pcaLeafStUlrichData <-plotPCA(leaf_stulrich_vst, intgroup=c("species", "treatment"), ntop = 1000, returnData = TRUE)
percentLeafStUlrichVar <- round(100 * attr(pcaLeafStUlrichData, "percentVar"))

pcaRootKitzbuhelData %<>% mutate("Tissue" = "Root", "Locality" = "Kitzbuhel")
pcaRootStUlrichData %<>% mutate("Tissue" = "Root", "Locality" = "St Ulrich")

pcaLeafKitzbuhelData %<>% mutate("Tissue" = "Leaf", "Locality" = "Kitzbuhel")
pcaLeafStUlrichData %<>% mutate("Tissue" = "Leaf", "Locality" = "St Ulrich")

pcaData <- rbind(pcaRootKitzbuhelData, pcaRootStUlrichData, pcaLeafKitzbuhelData, pcaLeafStUlrichData)


pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2A_RootKitzbuhel_ntop1000.pdf", width = 19, height = 19)
ggplot(pcaRootKitzbuhelData, aes(PC1, PC2)) +          
  geom_point(size = 35, 
             stroke = 13, aes(fill = species, colour=group, shape=species)) +
  xlab(paste0("PC1: ", percentRootKitzbuhelVar[1],"%")) +
  ylab(paste0("PC2: ", percentRootKitzbuhelVar[2],"%")) + 
  scale_colour_manual(values = c("gold", "dodgerblue", "dodgerblue", "gold")) +
  scale_fill_manual(values = c("gold", "dodgerblue")) +
  scale_shape_manual(values=c(22, 21)) +
  theme(legend.title = element_blank(),
        text = element_text(size = 85),
        axis.text=element_text(size=85),
        axis.title=element_text(size=85),
        legend.position="none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank())
dev.off()


pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2B_RootStUlrich_ntop1000.pdf", width = 19, height = 19)
ggplot(pcaRootStUlrichData, aes(PC1, PC2)) +          
  geom_point(size = 35, 
             stroke = 13, aes(fill = species, colour=group, shape=species)) +
  xlab(paste0("PC1: ", percentRootStUlrichVar[1],"%")) +
  ylab(paste0("PC2: ", percentRootStUlrichVar[2],"%")) + 
  scale_colour_manual(values = c("gold", "dodgerblue", "dodgerblue", "gold")) +
  scale_fill_manual(values = c("gold", "dodgerblue")) +
  scale_shape_manual(values=c(22, 21)) +
  theme(legend.title = element_blank(),
        text = element_text(size = 85),
        axis.text=element_text(size=85),
        axis.title=element_text(size=85),
        legend.position="none",
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank())
dev.off()


pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2C_LeafKitzbuhel_ntop1000.pdf", width = 19, height = 19)
ggplot(pcaLeafKitzbuhelData, aes(PC1, PC2)) +          
  geom_point(size = 35, 
             stroke = 13, aes(fill = species, colour=group, shape=species)) +
  xlab(paste0("PC1: ", percentLeafKitzbuhelVar[1],"%")) +
  ylab(paste0("PC2: ", percentLeafKitzbuhelVar[2],"%")) + 
  scale_colour_manual(values = c("gold", "dodgerblue", "dodgerblue", "gold")) +
  scale_fill_manual(values = c("gold", "dodgerblue")) +
  scale_shape_manual(values=c(22, 21)) +
  theme(legend.title = element_blank(),
        text = element_text(size = 85),
        axis.text=element_text(size=85),
        axis.title=element_text(size=85),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position="none",
        panel.background = element_blank())
dev.off()




pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2D_LeafStUlrich_ntop1000.pdf", width = 19, height = 19)
ggplot(pcaLeafStUlrichData, aes(PC1, PC2)) +          
  geom_point(size = 35, 
             stroke = 13, aes(fill = species, colour=group, shape=species)) +
  xlab(paste0("PC1: ", percentLeafStUlrichVar[1],"%")) +
  ylab(paste0("PC2: ", percentLeafStUlrichVar[2],"%")) + 
  scale_colour_manual(values = c("gold", "dodgerblue", "dodgerblue", "gold")) +
  scale_fill_manual(values = c("gold", "dodgerblue")) +
  scale_shape_manual(values=c(22, 21)) +
  theme(legend.title = element_blank(),
        text = element_text(size = 85),
        axis.text=element_text(size=85),
        axis.title=element_text(size=85),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.position="none",
        panel.background = element_blank())
dev.off()


########################################################################################
#  here I generate the barplot for around the schematic on Figure 2.                   #
# Need to run below result table ciode to get the reslts table and get the numbers.    #
########################################################################################

result_table


mtML <- result_table %>% filter(comparison %in% c("traunsteineri vs majalis (M) up", "traunsteineri vs majalis (M) down")) %>% dplyr::select(c("comparison", "Leaf.Kitzbuhel", "Leaf.St..Ulrich"))
mtTL <- result_table %>% filter(comparison %in% c("traunsteineri vs majalis (T) up", "traunsteineri vs majalis (T) down")) %>% dplyr::select(c("comparison", "Leaf.Kitzbuhel", "Leaf.St..Ulrich"))

mMmTL <- result_table %>% filter(comparison %in% c("majalis (M-T) up", "majalis (M-T) down")) %>% dplyr::select(c("comparison", "Leaf.Kitzbuhel", "Leaf.St..Ulrich"))
tTtML <- result_table %>% filter(comparison %in% c("traunsteineri (T-M) up", "traunsteineri (T-M) down")) %>% dplyr::select(c("comparison", "Leaf.Kitzbuhel", "Leaf.St..Ulrich"))



mtMR <- result_table %>% filter(comparison %in% c("traunsteineri vs majalis (M) up", "traunsteineri vs majalis (M) down")) %>% dplyr::select(c("comparison", "Root.Kitzbuhel", "Root.St..Ulrich"))
mtTR <- result_table %>% filter(comparison %in% c("traunsteineri vs majalis (T) up", "traunsteineri vs majalis (T) down")) %>% dplyr::select(c("comparison", "Root.Kitzbuhel", "Root.St..Ulrich"))

mMmTR <- result_table %>% filter(comparison %in% c("majalis (M-T) up", "majalis (M-T) down")) %>% dplyr::select(c("comparison", "Root.Kitzbuhel", "Root.St..Ulrich"))
tTtMR <- result_table %>% filter(comparison %in% c("traunsteineri (T-M) up", "traunsteineri (T-M) down")) %>% dplyr::select(c("comparison", "Root.Kitzbuhel", "Root.St..Ulrich"))


########################################################
########################################################
####                     LEAF                      #####
########################################################
########################################################

########################################
#               mtML                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_mtML.pdf", width = 12, height = 13)
ggplot(data=mtML %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,1100) +
  theme(text = element_text(size = 65), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(3, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.2, size = 50, position = "stack") 
dev.off()

########################################
#               mtTL                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_mtTL.pdf", width = 12, height = 13)
ggplot(data=mtTL %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,1050) +
  theme(text = element_text(size = 45), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.4, size = 50, position = "stack") 
dev.off()



########################################
#               mMmTL                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_mMmTL.pdf", width = 12, height = 13)
ggplot(data=mMmTL %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,1050) +
  theme(text = element_text(size = 45), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.4, size = 50, position = "stack") 
dev.off()



########################################
#               tTtML                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_tTtML.pdf", width = 12, height = 13)
ggplot(data=tTtML %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,1050) +
  theme(text = element_text(size = 45), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.4, size = 50, position = "stack") 
dev.off()




########################################################
########################################################
####                     ROOT                      #####
########################################################
########################################################
 
 



########################################
#               mtMR                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_mtMR.pdf", width = 12, height = 13)
ggplot(data=mtMR %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,2350) +
  theme(text = element_text(size = 65), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(3, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.15, size = 50, position = "stack") 
dev.off()

########################################
#               mtTR                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_mtTR.pdf", width = 12, height = 13)
ggplot(data=mtTR %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,2350) +
  theme(text = element_text(size = 45), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.4, size = 50, position = "stack") 
dev.off()



########################################
#               mMmTR                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_mMmTR.pdf", width = 12, height = 13)
ggplot(data=mMmTR %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,2350) +
  theme(text = element_text(size = 45), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.4, size = 50, position = "stack") 
dev.off()



########################################
#               tTtMR                  #
########################################

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2_tTtMR.pdf", width = 12, height = 13)
ggplot(data=tTtMR %>% melt() %>% mutate(Locality=c("", "K", "", "S")), aes(x=variable, y=value, fill=comparison)) +
  geom_bar(stat="identity", colour="black") +
  ylim(0,2350) +
  theme(text = element_text(size = 45), 
        legend.position="none",
        #legend.title=element_text(size=23),
        axis.text.y =element_text(size=120),
        axis.text.x = element_blank(),
        axis.title=element_blank(),
        legend.title=element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values=c("darkseagreen3", "burlywood2")) +
  geom_text(aes(label = Locality), group="term",
            vjust = -0.4, size = 50, position = "stack") 
dev.off()




#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2.png", width = 2000, height = 1600)
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure2.pdf", width = 28, height = 22)
pca_plot
dev.off()

percentRootKitzbuhelVar[1]
percentRootKitzbuhelVar[2]
percentRootStUlrichVar[1]
percentRootStUlrichVar[2]
percentLeafKitzbuhelVar[1]
percentLeafKitzbuhelVar[2]
percentLeafStUlrichVar[1]
percentLeafStUlrichVar[2]

###############################################################
#      Count how many genes are expressed total per tissue    #
###############################################################


# make a dds object from the total root samples (no subsetting)
root_dds<-specify_comparison(root_samples, df_counts_root, "1 == 1")
leaf_dds<-specify_comparison(leaf_samples, df_counts_leaf, "1 == 1")


root_dds <- DESeqDataSetFromMatrix(countData = root_dds[["counts"]],
                                             colData = root_dds[["samples"]],
                                             design = ~ species + treatment)
leaf_dds <- DESeqDataSetFromMatrix(countData = leaf_dds[["counts"]],
                                             colData = leaf_dds[["samples"]],
                                             design = ~ species + treatment)


leaf_kitzbuhel_names<-leaf_samples %>% filter(locality == "Kitzbuhl") %>% rownames()
leaf_stulrich_names<-leaf_samples %>% filter(locality == "St Ulrich") %>% rownames()

root_kitzbuhel_names<-root_samples %>% filter(locality == "Kitzbuhl") %>% rownames()
root_stulrich_names<-root_samples %>% filter(locality == "St Ulrich") %>% rownames()


mcols(leaf_dds)$basepairs<-df_leaf$lengths
leaf_dds <- estimateSizeFactors(leaf_dds)
length(which(rowSums(counts(leaf_dds, normalized=TRUE)[,leaf_kitzbuhel_names] >= 1 ) >= 3))
length(which(rowSums(counts(leaf_dds, normalized=TRUE)[,leaf_stulrich_names] >= 1 ) >= 3))


mcols(root_dds)$basepairs<-df_root$lengths
root_dds <- estimateSizeFactors(root_dds)
length(which(rowSums( counts(root_dds, normalized=TRUE)[,root_kitzbuhel_names] >= 1 ) >= 3))
length(which(rowSums( counts(root_dds, normalized=TRUE)[,root_stulrich_names] >= 1 ) >= 3))

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

#pvalue_threshold <- 0.05
#logfc_threshold<-1.5

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
                                         traunst_env_padj < 0.05 & majalis_env <= 0 & traunst_env <= 0 ~ "Constitutively DE",
                                      majalis_env_padj < 0.05 &
                                        traunst_env_padj < 0.05 & majalis_env >= 0 & traunst_env >= 0 ~ "Constitutively DE",
                                         majalis_env_padj < 0.05 &
                                         traunst_env_padj > 0.05 ~ "DE in D. majalis environment",
                                         majalis_env_padj > 0.05 &
                                         traunst_env_padj < 0.05 ~ "DE in D. traunsteineri environment"))


# designate all unlabelled status values Not Significant         
all_bound %<>% replace_na(list(status = "Not significant"))


# arrange the data so that the DE points are plotted last and on top of the grey non significant points
all_bound %<>% arrange(factor(status, levels = c("Not significant", 
                                                 "DE in D. majalis environment", 
                                                 "DE in D. traunsteineri environment", 
                                                 "Constitutively DE")))



#################################
#        Plot  Constitutive     #
#################################

bar_dataset<-all_bound %>% dplyr::select(status, locality, tissue) %>% filter(status != "Not significant") %>% melt()
bar <- ggplot(bar_dataset, aes(x=locality, fill=status)) +
  geom_bar(position="stack") + facet_wrap(~ tissue, ncol=2) +
  scale_fill_manual(values = c("DE in D. majalis environment" = "gold", 
                               "DE in D. traunsteineri environment" = "dodgerblue",
                               "Constitutively DE" = "deeppink")) +
  ylab("Number of Genes") +
  #xlab("Locality") +
  theme(text = element_text(size = 20), 
        legend.position = "none",
        axis.text=element_text(size=24),
        axis.title=element_text(size=42),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 40),
        axis.text.x = element_text(angle = 45, size = 30, vjust = 1, hjust=1))

scatt <-ggplot(all_bound, aes(x=majalis_env, y=traunst_env, colour=status)) +
  geom_point(alpha=0.5, size = 5) + 
  ylim(-10, 10) +
  xlim(-10, 10) +
  ylab("Fold change in traunsteineri environment") +
  xlab("Fold change in majalis environment") +
  scale_color_manual(values = c("DE in D. majalis environment" = "gold", 
                                "DE in D. traunsteineri environment" = "dodgerblue",
                                "Not significant" = "grey",
                                "Constitutively DE" = "deeppink")) + 
  facet_wrap(~ tissue + locality, ncol=2) +
  theme(text = element_text(size = 21), 
        legend.text=element_text(size=35),
        #legend.title=element_text(size=23),
        axis.text=element_text(size=30),
        axis.title=element_text(size=42),
        legend.title=element_blank(),
        strip.text.x = element_text(size = 40), 
        plot.margin = margin(0.5, 0.5, 0.5, 2, "cm")) +
  guides(colour = guide_legend(override.aes = list(size=10)))

#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure3.png", width = 1800, height = 900)
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure3.pdf", width = 26, height =12.5, onefile=FALSE)
egg::ggarrange(bar, scatt, ncol=2, widths = c(0.45, 1.1), labels = c('A', 'B'), label.args = list(gp = grid::gpar(font = 2, cex = 3.5)))
dev.off()

# how many DE genes are there in each category?
all_bound %>% dplyr::select("status") %>% pull() %>% table()

# how many genes are upregulated in traunsteineri relative to majalis?
all_bound %>% filter(majalis_env > 2 & traunst_env > 2 & majalis_env_padj < 0.05 & traunst_env_padj < 0.05)

# how many genes are upregulated in majalis relative to traunsteineri?
all_bound %>% filter(majalis_env > 2 & traunst_env > 2 & majalis_env_padj < 0.05 & traunst_env_padj < 0.05)



###################################################
#        Plot GO term enrichment Constitutive     #
###################################################


de_majalis_only_root<-all_bound %>% filter(status == "DE in D. majalis environment" & tissue == "Root") %>% dplyr::select(gene_id) %>% pull()
de_traunsteineri_only_root<-all_bound %>% filter(status == "DE in D. traunsteineri environment" & tissue == "Root") %>% dplyr::select(gene_id) %>% pull()
de_majalis_only_leaf<-all_bound %>% filter(status == "DE in D. majalis environment" & tissue == "Leaf") %>% dplyr::select(gene_id) %>% pull()
de_traunsteineri_only_leaf<-all_bound %>% filter(status == "DE in D. traunsteineri environment" & tissue == "Leaf") %>% dplyr::select(gene_id) %>% pull()
de_constitutive_leaf <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Leaf") %>% dplyr::select(gene_id) %>% pull()
de_constitutive_root <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Root") %>% dplyr::select(gene_id) %>% pull()

de_majalis_only_root_GO<-get_enriched_terms(de_majalis_only_root, mp)
de_traunsteineri_only_root_GO<-get_enriched_terms(de_traunsteineri_only_root, mp)
de_majalis_only_leaf_GO<-get_enriched_terms(de_majalis_only_leaf, mp)
de_traunsteineri_only_leaf_GO<-get_enriched_terms(de_traunsteineri_only_leaf, mp)
de_constitutive_leaf_GO<-get_enriched_terms(de_constitutive_leaf, mp)
de_constitutive_root_GO<-get_enriched_terms(de_constitutive_root, mp)

de_majalis_only_root_GO_df <- de_majalis_only_root_GO %>% prepare_go_df() %>% mutate(Environment="D. majalis Environment Only", tissue="Root", placeholder="placeholder")
de_traunsteineri_only_root_GO_df <- de_traunsteineri_only_root_GO %>% prepare_go_df()  %>% mutate(Environment="D. traunsteineri Environment Only", tissue="Root", placeholder="placeholder")
de_majalis_only_leaf_GO_df <- de_majalis_only_leaf_GO %>% prepare_go_df()  %>% mutate(Environment="D. majalis Environment Only", tissue="Leaf", placeholder="placeholder")
de_traunsteineri_only_leaf_GO_df <- de_traunsteineri_only_leaf_GO %>% prepare_go_df()  %>% mutate(Environment="D. traunsteineri Environment Only", tissue="Leaf", placeholder="placeholder")
de_constitutive_leaf_GO_df <- de_constitutive_leaf_GO %>% prepare_go_df() %>% mutate(Environment="Constitutively DE", tissue="Leaf", placeholder="placeholder")
de_constitutive_root_GO_df<- de_constitutive_root_GO %>% prepare_go_df() %>% mutate(Environment="Constitutively DE", tissue="Root", placeholder="placeholder")

root_go_bound <- rbind(de_majalis_only_root_GO_df, de_traunsteineri_only_root_GO_df, de_constitutive_root_GO_df)
root_go_bound$Term[root_go_bound$Term == "regulation of transcription from RNA polymerase II promoter in response to stress"] <- "reg of transcription from RNApolII promoter in resp. to stress"
leaf_go_bound <- rbind(de_majalis_only_leaf_GO_df, de_traunsteineri_only_leaf_GO_df, de_constitutive_leaf_GO_df)

#root_go_bound[root_go_bound$Term == "regulation of transcription from RNA polymerase II promoter in response to stress",]$Term = "reg. of transcription from RNApolII promoter in response to stress"
#leaf_go_bound[leaf_go_bound$Term == "pectin catabolic process" & leaf_go_bound$Environment == "D. traunsteineri Environment Only",]$Term = "pectin catabolic process "

root_go_bound %<>% mutate(Term = fct_reorder(Term, Environment))
# update duplicated value so it doesnt cause problems in factor redordering
#leaf_go_bound[leaf_go_bound$Term == "sodium ion import across plasma membrane" & leaf_go_bound$Environment == "Constitutively DE",]$Term <- "sodium ion import across plasma membrane "
leaf_go_bound %<>% mutate(Term = fct_reorder(Term, Environment))



a<-ggplot(leaf_go_bound, aes(x=placeholder, y=Term, color = Environment)) + 
  geom_point(size=47) + facet_grid(scales="free", space= "free", switch = "y", cols=vars(tissue)) + 
  theme(text = element_text(size = 102), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.x = element_blank(),
        #strip.text.y = element_text(size = 90),
        #strip.text.x = element_text(size = 90),
        legend.title=element_blank(),
        panel.spacing=unit(1, "lines"),
        plot.margin = margin(0, 0, 0, 0, "cm")) + 
  scale_color_manual(values=c("deeppink", "gold", "dodgerblue")) +
  theme(legend.position = "none") 


b <- ggplot(root_go_bound, aes(x=placeholder, y=Term, color = Environment)) + 
  geom_point(size=47) + facet_grid(scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 102), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        #strip.text.y = element_text(size = 90),
        #strip.text.x = element_text(size = 90),
        strip.text.x = element_blank(),
        legend.title=element_blank(),
        legend.position = "none",
        #legend.position = c(0.1,0.2),
        panel.spacing=unit(1, "lines"),
        plot.margin = margin(0, 11, 0, 0, "cm")) + 
  scale_color_manual(values=c("deeppink", "gold", "dodgerblue")) +
  scale_y_discrete(position = "right") 


b_legend<-ggplot(leaf_go_bound, aes(x=placeholder, y=Term, color = Environment)) + 
  geom_point() + facet_grid(scales="free", space= "free", cols=vars(tissue)) + 
  theme(legend.margin=margin(c(0,0,0,0)),
        #legend.position = c(1,1),
        legend.text=element_text(size=83),
        legend.title=element_text(size=83)) + 
  scale_color_manual(values=c("deeppink", "gold", "dodgerblue")) +
  guides(color = guide_legend(override.aes = list(size = 30)))


leg <- cowplot::get_legend(b_legend)
leg<-as_ggplot(leg)

#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure4.png", height=2000, width=7000)
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure4.pdf", height=35, width=86)
plot_grid(a, b, NULL, leg,
          ncol = 4,
          rel_widths = c(2, 2, -0.8, 1.2),
          rel_heights = c(2, 2, -0.8, 1.2))
dev.off()


#############################################################################
#    Make a table of constitutive GO terms total for supplementary table4  #
############################################################################


de_majalis_only_root_GO_df <- de_majalis_only_root_GO %>% prepare_go_df(top_terms = 1000) %>% mutate(Environment="D. majalis Environment Only", tissue="Root")
de_traunsteineri_only_root_GO_df <- de_traunsteineri_only_root_GO %>% prepare_go_df(top_terms = 1000)  %>% mutate(Environment="D. traunsteineri Environment Only", tissue="Root")
de_majalis_only_leaf_GO_df <- de_majalis_only_leaf_GO %>% prepare_go_df(top_terms = 1000)  %>% mutate(Environment="D. majalis Environment Only", tissue="Leaf")
de_traunsteineri_only_leaf_GO_df <- de_traunsteineri_only_leaf_GO %>% prepare_go_df(top_terms = 1000)  %>% mutate(Environment="D. traunsteineri Environment Only", tissue="Leaf")
de_constitutive_leaf_GO_df <- de_constitutive_leaf_GO %>% prepare_go_df(top_terms = 1000) %>% mutate(Environment="Constitutively DE", tissue="Leaf")
de_constitutive_root_GO_df<- de_constitutive_root_GO %>% prepare_go_df(top_terms = 1000) %>% mutate(Environment="Constitutively DE", tissue="Root")

rbind(de_majalis_only_root_GO_df,
      de_traunsteineri_only_root_GO_df,
      de_majalis_only_leaf_GO_df,
      de_traunsteineri_only_leaf_GO_df,
      de_constitutive_leaf_GO_df,
      de_constitutive_root_GO_df) %>% 
  write.xlsx(file = "SupplementaryTable5_constitutive_and_specific_GOterms.xlsx")




######################################################################
#    Test for whether more constitutively DEGs are found by chance   #
######################################################################

rootK_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Root" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
rootK_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Root" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
rootK_DEG_M <- c(rootK_DEG_M, rootK_consititutive)
rootK_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Root" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
rootK_DEG_T <- c(rootK_DEG_T, rootK_consititutive)

rootS_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Root" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
rootS_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Root" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
rootS_DEG_M <- c(rootS_DEG_M, rootS_consititutive)
rootS_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Root" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
rootS_DEG_T <- c(rootS_DEG_T, rootS_consititutive)


leafK_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Leaf" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
leafK_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Leaf" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
leafK_DEG_M <- c(leafK_DEG_M, leafK_consititutive)
leafK_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Leaf" & locality == "Kitzbuhl") %>% dplyr::select(gene_id) %>% pull()
leafK_DEG_T <- c(leafK_DEG_T, leafK_consititutive)

leafS_consititutive <- all_bound %>% filter(status == "Constitutively DE" & tissue == "Leaf" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
leafS_DEG_M <- all_bound %>% filter(status == "DE in D. majalis environment only" & tissue == "Leaf" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
leafS_DEG_M <- c(leafS_DEG_M, leafS_consititutive)
leafS_DEG_T<-all_bound %>% filter(status == "DE in D. traunsteineri environment only" & tissue == "Leaf" & locality == "St. Ulrich") %>% dplyr::select(gene_id) %>% pull()
leafS_DEG_T <- c(leafS_DEG_T, leafS_consititutive)




total_genes_tested_root<-unique(c(ktz_root$gene_id, stu_root$gene_id)) %>% length()
total_genes_tested_leaf<-unique(c(ktz_leaf$gene_id, stu_leaf$gene_id)) %>% length()

ResultRoot=supertest(list(rootK_DEG_M=rootK_DEG_M, 
                          rootK_DEG_T=rootK_DEG_T,
                          rootS_DEG_T=rootS_DEG_T,
                          rootS_DEG_M=rootS_DEG_M), n=total_genes_tested_root) 

ResultLeaf=supertest(list(leafK_DEG_M=leafK_DEG_M, 
                          leafK_DEG_T=leafK_DEG_T, 
                          leafS_DEG_T=leafS_DEG_T, 
                          leafS_DEG_M=leafS_DEG_M), n=total_genes_tested_leaf) 


summary(ResultRoot)$Table %>% dplyr::select(-"Elements") %>% write.xlsx(file = "SupplementaryTable2_Root_Exact_Test.xlsx")
summary(ResultLeaf)$Table %>% dplyr::select(-"Elements") %>% write.xlsx(file = "SupplementaryTable3_Leaf_Exact_Test.xlsx")

#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Supplemetary_Fig1.png", width = 1000, height = 1000)
plot(ResultRoot, Layout="landscape", degree=2:7, sort.by="size")
#dev.off()
#
#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Supplemetary_Fig2.png", width = 1000, height = 1000)
plot(ResultLeaf, Layout="landscape", degree=2:7, sort.by="size")
#dev.off()




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


#############################################
#    Venn Diagram of number of DE genes     #
#############################################


a<-get_significant_genes(transplant_majalis_kitzbuhl_leaf)
b<-get_significant_genes(transplant_majalis_stulrich_leaf)
c<-get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf)
d<-get_significant_genes(transplant_traunsteineri_stulrich_leaf)


length(intersect(a, b))/length(unique(c(a, b))) * 100
length(intersect(c, d))/length(unique(c(c, d))) * 100

length(intersect(d, e))/length(unique(c(d, e))) * 100
length(intersect(f, g))/length(unique(c(f, g))) * 100


library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

myCol<-c("deeppink", "yellowgreen", "deepskyblue", "orange1")


#pdf("Figure5A_Leaf.pdf", height = 25, width = 25)
#png("Figure5A_Leaf.png", height = 2000, width = 2000)
leaf_venn<-ggvenn(
  list(" "=a, "  "=b, "   "=c, "    "=d), 
  fill_color = myCol,
  stroke_size = 0.9, 
  set_name_size = 21.5, 
  text_size = 28, 
  fill_alpha = 0.2,
  show_percentage = FALSE
) + theme(plot.margin = unit(c(2,2,2,2), "cm"))
#dev.off()

e<-get_significant_genes(transplant_majalis_kitzbuhl_root)
f<-get_significant_genes(transplant_majalis_stulrich_root)
g<-get_significant_genes(transplant_traunsteineri_kitzbuhl_root)
h<-get_significant_genes(transplant_traunsteineri_stulrich_root)

#pdf("Figure5A_Root.pdf", height = 25, width = 25)
#png("Figure5A_Root.png", height = 2000, width = 2000)
root_venn<-ggvenn(
  list(" "=e, "  "=f, "   "=g, "    "=h), 
  fill_color = myCol,
  stroke_size = 0.9, 
  set_name_size = 21.5, 
  text_size = 28, 
  fill_alpha = 0.2,
  show_percentage = FALSE
) + theme(plot.margin = unit(c(2,2,2,2), "cm"))
#dev.off()


pdf("Figure5.pdf", height = 25, width = 45)
cowplot::plot_grid(leaf_venn, root_venn)
                   #labels = c('Leaf', 
                   #           'Root'),
                   #label_size = 80,
                   #label_x = c(0.4, 0.4))
dev.off()


#####################################
#        DE gene count plots        #
#####################################

# make a function to label up or down differential expression (significant)
# mainly so you dont have to repeat the code
label_expression_direction<-function(results_object){
  results_object %>% data.frame() %>% 
    mutate(diffexpressed=case_when(log2FoldChange >= 1.5 & padj <= 0.05 ~ "upregulated",
                                   log2FoldChange <= -1.5 & padj <= 0.05 ~ "downregulated")) %>%
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
de_counts$Tissue<-paste(de_counts$Tissue, "Plastic")
colnames(de_counts)<-c("Species", "Tissue", "Locality", "Individual", "Direction", "Number of genes")

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure5B.png", width = 900, height = 1200)
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
  ylim(-580, 580)
dev.off()


######**** NB after this I took 5A and 5B figureds and collated them in a pptx to create final composite figure ****##########



#########################################
#    Gather all DEGs and make a table   #
#########################################


# plastic genes (within species comparison)

mTK_leaf_degs<-get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE)
mTS_leaf_degs<-get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE)
mTK_root_degs<-get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE)
mTS_root_degs<-get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE)


tMK_leaf_degs<-get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE)
tMS_leaf_degs<-get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE)
tMK_root_degs<-get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE)
tMS_root_degs<-get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE)


`majalis (M-T Kitzbuhel) leaf up` <- mTK_leaf_degs$up %>% length()
`majalis (M-T Kitzbuhel) leaf down` <- mTK_leaf_degs$down %>% length()
`majalis (M-T Kitzbuhel) root up` <- mTK_root_degs$up %>% length()
`majalis (M-T Kitzbuhel) root down` <- mTK_root_degs$down %>% length()
`traunsteineri (T-M Kitzbuhel) leaf up` <- tMK_leaf_degs$up %>% length()
`traunsteineri (T-M Kitzbuhel) leaf down` <- tMK_leaf_degs$down %>% length()
`traunsteineri (T-M Kitzbuhel) root up` <- tMK_root_degs$up %>% length()
`traunsteineri (T-M Kitzbuhel) root down` <- tMK_root_degs$down %>% length()

`majalis (M-T St. Ulrich) leaf up` <- mTS_leaf_degs$up %>% length()
`majalis (M-T St. Ulrich) leaf down` <- mTS_leaf_degs$down %>% length()
`majalis (M-T St. Ulrich) root up` <- mTS_root_degs$up %>% length()
`majalis (M-T St. Ulrich) root down` <- mTS_root_degs$down %>% length()
`traunsteineri (T-M St. Ulrich) leaf up` <- tMS_leaf_degs$up %>% length()
`traunsteineri (T-M St. Ulrich) leaf down` <- tMS_leaf_degs$down %>% length()
`traunsteineri (T-M St. Ulrich) root up` <- tMS_root_degs$up %>% length()
`traunsteineri (T-M St. Ulrich) root down` <- tMS_root_degs$down %>% length()


# consitutively differentially expressed genes (between species comparisons)

tmMK_leaf_degs<-get_significant_genes(traunsteineri_majalis_leaf_M_kitzbuhl, directional = TRUE)
tmTK_leaf_degs<-get_significant_genes(traunsteineri_majalis_leaf_T_kitzbuhl, directional = TRUE)
tmMK_root_degs<-get_significant_genes(traunsteineri_majalis_root_M_kitzbuhl, directional = TRUE)
tmTK_root_degs<-get_significant_genes(traunsteineri_majalis_root_T_kitzbuhl, directional = TRUE)

tmMS_leaf_degs<-get_significant_genes(traunsteineri_majalis_leaf_M_stulrich, directional = TRUE)
tmTS_leaf_degs<-get_significant_genes(traunsteineri_majalis_leaf_T_stulrich, directional = TRUE)
tmMS_root_degs<-get_significant_genes(traunsteineri_majalis_root_M_stulrich, directional = TRUE)
tmTS_root_degs<-get_significant_genes(traunsteineri_majalis_root_T_stulrich, directional = TRUE)


`traunsteineri vs majalis (M Kitzbuhel) leaf up` <- tmMK_leaf_degs$up %>% length()
`traunsteineri vs majalis (M Kitzbuhel) leaf down` <- tmMK_leaf_degs$down %>% length()
`traunsteineri vs majalis (T Kitzbuhel) leaf up` <- tmTK_leaf_degs$up %>% length()
`traunsteineri vs majalis (T Kitzbuhel) leaf down` <- tmTK_leaf_degs$down %>% length()
`traunsteineri vs majalis (M Kitzbuhel) root up` <- tmMK_root_degs$up %>% length()
`traunsteineri vs majalis (M Kitzbuhel) root down` <- tmMK_root_degs$down %>% length()
`traunsteineri vs majalis (T Kitzbuhel) root up` <- tmTK_root_degs$up %>% length()
`traunsteineri vs majalis (T Kitzbuhel) root down` <- tmTK_root_degs$down %>% length()

`traunsteineri vs majalis (M St. Ulrich) leaf up` <- tmMS_leaf_degs$up %>% length()
`traunsteineri vs majalis (M St. Ulrich) leaf down` <- tmMS_leaf_degs$down %>% length()
`traunsteineri vs majalis (T St. Ulrich) leaf up` <- tmTS_leaf_degs$up %>% length()
`traunsteineri vs majalis (T St. Ulrich) leaf down` <- tmTS_leaf_degs$down %>% length()
`traunsteineri vs majalis (M St. Ulrich) root up` <- tmMS_root_degs$up %>% length()
`traunsteineri vs majalis (M St. Ulrich) root down` <- tmMS_root_degs$down %>% length()
`traunsteineri vs majalis (T St. Ulrich) root up` <- tmTS_root_degs$up %>% length()
`traunsteineri vs majalis (T St. Ulrich) root down` <- tmTS_root_degs$down %>% length()



`traunsteineri vs majalis (Kitzbuhel) Leaf Mup Tup` <- intersect(tmMK_leaf_degs$up, tmTK_leaf_degs$up) %>% length()
`traunsteineri vs majalis (Kitzbuhel) Leaf Mdown Tdown` <- intersect(tmMK_leaf_degs$down, tmTK_leaf_degs$down) %>% length()
`traunsteineri vs majalis (Kitzbuhel) Leaf Mup Tdown` <- intersect(tmMK_leaf_degs$up, tmTK_leaf_degs$down) %>% length()
`traunsteineri vs majalis (Kitzbuhel) Leaf Mdown Tup` <- intersect(tmMK_leaf_degs$down, tmTK_leaf_degs$up) %>% length()

`traunsteineri vs majalis (St. Ulrich) Leaf Mup Tup` <- intersect(tmMS_leaf_degs$up, tmTS_leaf_degs$up) %>% length()
`traunsteineri vs majalis (St. Ulrich) Leaf Mdown Tdown` <- intersect(tmMS_leaf_degs$down, tmTS_leaf_degs$down) %>% length()
`traunsteineri vs majalis (St. Ulrich) Leaf Mup Tdown` <- intersect(tmMS_leaf_degs$up, tmTS_leaf_degs$down) %>% length()
`traunsteineri vs majalis (St. Ulrich) Leaf Mdown Tup` <- intersect(tmMS_leaf_degs$down, tmTS_leaf_degs$up) %>% length()

`traunsteineri vs majalis (Kitzbuhel) Root Mup_Tup` <- intersect(tmMK_root_degs$up, tmTK_root_degs$up) %>% length()
`traunsteineri vs majalis (Kitzbuhel) Root Mdown_Tdown` <- intersect(tmMK_root_degs$down, tmTK_root_degs$down) %>% length()
`traunsteineri vs majalis (Kitzbuhel) Root Mup_Tdown` <- intersect(tmMK_root_degs$up, tmTK_root_degs$down) %>% length()
`traunsteineri vs majalis (Kitzbuhel) Root Mdown_Tup` <- intersect(tmMK_root_degs$down, tmTK_root_degs$up) %>% length()

`traunsteineri vs majalis (St. Ulrich) Root Mup_Tup` <- intersect(tmMS_root_degs$up, tmTS_root_degs$up) %>% length()
`traunsteineri vs majalis (St. Ulrich) Root Mdown_Tdown` <- intersect(tmMS_root_degs$down, tmTS_root_degs$down) %>% length()
`traunsteineri vs majalis (St. Ulrich) Root Mup_Tdown` <- intersect(tmMS_root_degs$up, tmTS_root_degs$down) %>% length()
`traunsteineri vs majalis (St. Ulrich) Root Mdown_Tup` <- intersect(tmMS_root_degs$down, tmTS_root_degs$up) %>% length()

result_table %>% filter(comparison %in% c("traunsteineri vs majalis Mup_Tup", "traunsteineri vs majalis Mdown_Tdown"))


result_table <- data.frame(comparison=c("traunsteineri vs majalis (M) up", 
                                      "traunsteineri vs majalis (M) down",
                                      "traunsteineri vs majalis (T) up", 
                                      "traunsteineri vs majalis (T) down",
                                      "traunsteineri vs majalis Mup_Tup", 
                                      "traunsteineri vs majalis Mdown_Tdown", 
                                      "traunsteineri vs majalis Mup_Tdown", 
                                      "traunsteineri vs majalis Mdown_Tup",
                                      "majalis (M-T) up",
                                      "majalis (M-T) down",
                                      "traunsteineri (T-M) up",
                                      "traunsteineri (T-M) down"),
                           `Leaf Kitzbuhel`=c(`traunsteineri vs majalis (M Kitzbuhel) leaf up`, 
                                              `traunsteineri vs majalis (M Kitzbuhel) leaf down`,
                                              `traunsteineri vs majalis (T Kitzbuhel) leaf up`,
                                              `traunsteineri vs majalis (T Kitzbuhel) leaf down`,
                                              `traunsteineri vs majalis (Kitzbuhel) Leaf Mup Tup`, 
                                              `traunsteineri vs majalis (Kitzbuhel) Leaf Mdown Tdown`, 
                                              `traunsteineri vs majalis (Kitzbuhel) Leaf Mup Tdown`, 
                                              `traunsteineri vs majalis (Kitzbuhel) Leaf Mdown Tup`,
                                              `majalis (M-T Kitzbuhel) leaf up`, 
                                              `majalis (M-T Kitzbuhel) leaf down`, 
                                              `traunsteineri (T-M Kitzbuhel) leaf up`, 
                                              `traunsteineri (T-M Kitzbuhel) leaf down`),
                           `Leaf St. Ulrich`=c(`traunsteineri vs majalis (M St. Ulrich) leaf up`,
                                               `traunsteineri vs majalis (M St. Ulrich) leaf down`,
                                               `traunsteineri vs majalis (T St. Ulrich) leaf up`,
                                               `traunsteineri vs majalis (T St. Ulrich) leaf down`,
                                               `traunsteineri vs majalis (St. Ulrich) Leaf Mup Tup`, 
                                               `traunsteineri vs majalis (St. Ulrich) Leaf Mdown Tdown`, 
                                               `traunsteineri vs majalis (St. Ulrich) Leaf Mup Tdown`, 
                                               `traunsteineri vs majalis (St. Ulrich) Leaf Mdown Tup`,
                                               `majalis (M-T St. Ulrich) leaf up`, 
                                               `majalis (M-T St. Ulrich) leaf down`,
                                               `traunsteineri (T-M St. Ulrich) leaf up`,
                                               `traunsteineri (T-M St. Ulrich) leaf down`),
                           `Root Kitzbuhel`=c(`traunsteineri vs majalis (M Kitzbuhel) root up`,
                                              `traunsteineri vs majalis (M Kitzbuhel) root down`,
                                              `traunsteineri vs majalis (T Kitzbuhel) root up`,
                                              `traunsteineri vs majalis (T Kitzbuhel) root down`,
                                              `traunsteineri vs majalis (Kitzbuhel) Root Mup_Tup`, 
                                              `traunsteineri vs majalis (Kitzbuhel) Root Mdown_Tdown`, 
                                              `traunsteineri vs majalis (Kitzbuhel) Root Mup_Tdown`, 
                                              `traunsteineri vs majalis (Kitzbuhel) Root Mdown_Tup`,
                                              `majalis (M-T Kitzbuhel) root up`,
                                              `majalis (M-T Kitzbuhel) root down`,
                                              `traunsteineri (T-M Kitzbuhel) root up`,
                                              `traunsteineri (T-M Kitzbuhel) root down`),
                           `Root St. Ulrich`=c(`traunsteineri vs majalis (M St. Ulrich) root up`, 
                                               `traunsteineri vs majalis (M St. Ulrich) root down`,
                                               `traunsteineri vs majalis (T St. Ulrich) root up`, 
                                               `traunsteineri vs majalis (T St. Ulrich) root down`,
                                               `traunsteineri vs majalis (St. Ulrich) Root Mup_Tup`, 
                                               `traunsteineri vs majalis (St. Ulrich) Root Mdown_Tdown`, 
                                               `traunsteineri vs majalis (St. Ulrich) Root Mup_Tdown`, 
                                               `traunsteineri vs majalis (St. Ulrich) Root Mdown_Tup`,
                                               `majalis (M-T St. Ulrich) root up`,
                                               `majalis (M-T St. Ulrich) root down`,
                                               `traunsteineri (T-M St. Ulrich) root up`,
                                               `traunsteineri (T-M St. Ulrich) root down`))


write.xlsx(result_table, file = "SupplementaryTable4_differentially_expressed_genes.xlsx")

##################################################################################################
#    find number of DEGs in each environment per locality, regardless of tissue or direction     #
################################################################################################## 

# Majalis environment Kitzbuhel
result_table %>% filter(comparison %in% c("traunsteineri vs majalis (M) up", "traunsteineri vs majalis (M) down")) %>% dplyr::select("Leaf.Kitzbuhel", "Root.Kitzbuhel") %>% colSums() %>% sum()
# Majalis environment St. Ulrich
result_table %>% filter(comparison %in% c("traunsteineri vs majalis (M) up", "traunsteineri vs majalis (M) down")) %>% dplyr::select("Leaf.St..Ulrich", "Root.St..Ulrich") %>% colSums() %>% sum()

# Traunsteineri environment Kitzbuhel
result_table %>% filter(comparison %in% c("traunsteineri vs majalis (T) up", "traunsteineri vs majalis (T) down")) %>% dplyr::select("Leaf.Kitzbuhel", "Root.Kitzbuhel") %>% colSums() %>% sum()
# Traunsteineri environment St. Ulrich
result_table %>% filter(comparison %in% c("traunsteineri vs majalis (T) up", "traunsteineri vs majalis (T) down")) %>% dplyr::select("Leaf.St..Ulrich", "Root.St..Ulrich") %>% colSums() %>% sum()






#############################################
#            GO term enrichment             #
#############################################

#######################################
#        majalis leaf kitzbuhl        #
#######################################
transplant_majalis_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE)$up, mp) %>% mutate(species="D. majalis", locality="Kitzbuhel", direction="up", tissue="leaf")
transplant_majalis_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_leaf, directional = TRUE)$down, mp) %>% mutate(species="D. majalis", locality="Kitzbuhel", direction="down", tissue="leaf")

#######################################
#        majalis leaf st ulrich        #
#######################################
transplant_majalis_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE)$up, mp) %>% mutate(species="D. majalis", locality="St. Ulrich", direction="up", tissue="leaf")
transplant_majalis_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_leaf, directional = TRUE)$down, mp) %>% mutate(species="D. majalis", locality="St. Ulrich", direction="down", tissue="leaf") 



#######################################
#        majalis root kitzbuhl        #
#######################################
transplant_majalis_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE)$up, mp) %>% mutate(species="D. majalis", locality="Kitzbuhel", direction="up", tissue="root")
transplant_majalis_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_kitzbuhl_root, directional = TRUE)$down, mp) %>% mutate(species="D. majalis", locality="Kitzbuhel", direction="down", tissue="root") 

#######################################
#        majalis root st ulrich       #
#######################################
transplant_majalis_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE)$up, mp) %>% mutate(species="D. majalis", locality="St. Ulrich", direction="up", tissue="root")
transplant_majalis_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_majalis_stulrich_root, directional = TRUE)$down, mp) %>% mutate(species="D. majalis", locality="St. Ulrich", direction="down", tissue="root") 

#######################################
#     traunsteineri leaf kitzbuhl     #
#######################################
transplant_traunsteineri_kitzbuhl_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE)$up, mp) %>% mutate(species="D. traunsteineri", locality="Kitzbuhel", direction="up", tissue="leaf")
transplant_traunsteineri_kitzbuhl_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_leaf, directional = TRUE)$down, mp) %>% mutate(species="D. traunsteineri", locality="Kitzbuhel", direction="down", tissue="leaf")


#######################################
#     traunsteineri leaf st ulrich    #
#######################################
transplant_traunsteineri_stulrich_leaf_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE)$up, mp) %>% mutate(species="D. traunsteineri", locality="St. Ulrich", direction="up", tissue="leaf")
transplant_traunsteineri_stulrich_leaf_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_leaf, directional = TRUE)$down, mp) %>% mutate(species="D. traunsteineri", locality="St. Ulrich", direction="down", tissue="leaf")

#######################################
#     traunsteineri root kitzbuhl     #
#######################################
transplant_traunsteineri_kitzbuhl_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE)$up, mp) %>% mutate(species="D. traunsteineri", locality="Kitzbuhel", direction="up", tissue="root")
transplant_traunsteineri_kitzbuhl_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_kitzbuhl_root, directional = TRUE)$down, mp) %>% mutate(species="D. traunsteineri", locality="Kitzbuhel", direction="down", tissue="root")

#######################################
#     traunsteineri root st ulrich    #
#######################################
transplant_traunsteineri_stulrich_root_up<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE)$up, mp) %>% mutate(species="D. traunsteineri", locality="St. Ulrich", direction="up", tissue="root")
transplant_traunsteineri_stulrich_root_down<-get_enriched_terms(get_significant_genes(transplant_traunsteineri_stulrich_root, directional = TRUE)$down, mp) %>% mutate(species="D. traunsteineri", locality="St. Ulrich", direction="down", tissue="root")


#########################################################
#     traunsteineri root st ulrich testing extraction   #
#########################################################



# bind all annotated go enrichment results togethert, annotated
all_enrichment_results<-rbind(transplant_majalis_kitzbuhl_root_up,
      transplant_majalis_kitzbuhl_leaf_down,
      transplant_majalis_stulrich_root_up,
      transplant_majalis_stulrich_root_down,
      transplant_traunsteineri_kitzbuhl_leaf_up,
      transplant_traunsteineri_kitzbuhl_leaf_down,
      transplant_traunsteineri_stulrich_leaf_up,
      transplant_traunsteineri_stulrich_leaf_down,
      transplant_traunsteineri_kitzbuhl_root_up,
      transplant_traunsteineri_kitzbuhl_root_down,
      transplant_traunsteineri_stulrich_root_up,
      transplant_traunsteineri_stulrich_root_down) %>% filter(as.numeric(classicFisher) < 0.05)



# make a function to get the GOslim terms for an input of GO terms
get_goslim<-function(mylist){
  fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
  slim <- getOBOCollection(fl)
  myIDs<-mylist
  myCollection <- GOCollection(myIDs)
  slimbp<-goSlim(myCollection, slim, "BP") %>% filter(Count > 0)
  gomap <- as.list(GOBPOFFSPRING[rownames(slimbp)])
  mapped <- lapply(gomap, intersect, ids(myCollection))
  mapped_dataframe<-data.frame(GO.ID=mapped %>% unlist() %>% as.character(),
                               GO_slim=mapped %>% unlist2() %>% names())
  return(list(go_table=slimbp, go_dataframe=mapped_dataframe))
}

majalis_in<-all_enrichment_results %>% filter(species == "D. majalis") %>% dplyr::select("GO.ID") %>% pull()
trauns_in<-all_enrichment_results %>% filter(species == "D. traunsteineri") %>% dplyr::select("GO.ID") %>% pull()

majalis_out<-get_goslim(majalis_in)
trauns_out<-get_goslim(trauns_in)



toplot<-inner_join(majalis_out$go_table, trauns_out$go_table, by="Term") %>% 
  set_colnames(c("count_majalis", "percent_majalis", "term", "count_traunsteineri", "percent_traunsteineri")) %>%
  dplyr::select(term, percent_majalis, percent_traunsteineri) %>%
  mutate(delta=abs(percent_majalis-percent_traunsteineri)/(percent_majalis+percent_traunsteineri), 
         percent_total=percent_majalis+percent_traunsteineri,
         to_asterisk=ifelse(delta > 0.5, "*", ""))

toplot$term <- reorder(toplot$term, toplot$percent_total)
nonspecific_terms<-c("biological_process", "cellular process", "metabolic process", "catabolic process", "multicellular organism development")
toplot %<>% filter(!(term %in% nonspecific_terms))


test_toplot<-toplot %>% dplyr::select(-c("delta", "percent_total")) %>% 
  melt() %>% 
  mutate(to_asterisk=ifelse((variable == "percent_majalis" & to_asterisk == "*"), "*", "")) %>%
  ggplot(aes(x=term, y=value, fill=variable)) +
  geom_bar(stat="identity", color = "gray34") + 
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1, size=65),
        axis.text.y = element_text(size=50),
        axis.title.y = element_text(size=50),
        legend.position = c(0.18, 0.8),
        legend.text=element_text(size=80),
        legend.background=element_blank(),
        legend.title=element_blank(),
        #plot.margin = margin("top", "right, "bottom", "left, "cm")) +
        plot.margin = margin(5, 2, 2, 2, "cm")) +
  geom_text(aes(label = to_asterisk), group="term",
            vjust = 0, size = 50, position = "stack") + 
  scale_fill_manual(values = c("gold", "dodgerblue"),
                    labels = c("D. majalis", "D. traunsteineri")) +
  ylab("% terms in GOSlim Category") +
  xlab("") + 
  guides(colour=guide_legend(override.aes=list(shape=17)))



interesting <- toplot %>% filter(to_asterisk == "*") %>% dplyr::select(term) %>% pull() %>% as.character()
interesting_go<-unique(majalis_out$go_table %>% filter(Term %in% interesting) %>% rownames(), 
       trauns_out$go_table %>% filter(Term %in% interesting) %>% rownames())

majalis_leaf_interesting <- all_enrichment_results %>% filter(species == "D. majalis" & tissue == "leaf" & GO.ID %in% (majalis_out$go_dataframe %>% filter(GO_slim %in% interesting_go) %>% dplyr::select(GO.ID) %>% pull()))
majalis_root_interesting <- all_enrichment_results %>% filter(species == "D. majalis" & tissue == "root" & GO.ID %in% (majalis_out$go_dataframe %>% filter(GO_slim %in% interesting_go) %>% dplyr::select(GO.ID) %>% pull()))
traunsteineri_leaf_interesting <- all_enrichment_results %>% filter(species == "D. traunsteineri" & tissue == "leaf" & GO.ID %in% (trauns_out$go_dataframe %>% filter(GO_slim %in% interesting_go) %>% dplyr::select(GO.ID) %>% pull()))
traunsteineri_root_interesting <- all_enrichment_results %>% filter(species == "D. traunsteineri" & tissue == "root" & GO.ID %in% (trauns_out$go_dataframe %>% filter(GO_slim %in% interesting_go) %>% dplyr::select(GO.ID) %>% pull()))


leaf_go_bound <- rbind(majalis_leaf_interesting, traunsteineri_leaf_interesting) %>% dplyr::select("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "direction", "locality", "species")
root_go_bound <- rbind(majalis_root_interesting, traunsteineri_root_interesting) %>% dplyr::select("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "direction", "locality", "species")




colnames(leaf_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Direction", "Locality", "Species")
colnames(root_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Direction", "Locality", "Species")

# I am adding a new column with a placeholder the same across all rows
# I will use this column as my X axis for the GO term plots
# if I use the comparison column, it staggers the plots, and it looks better if the points are stacked
leaf_go_bound_newcol <- leaf_go_bound %>% mutate(newcol="placeholder", tissue="Leaf")
root_go_bound_newcol <- root_go_bound %>% mutate(newcol="placeholder", tissue="Root")

leaf_go_bound_newcol$Term <- leaf_go_bound_newcol$Term %>% as.character() %>% make.unique(sep = " ")
root_go_bound_newcol$Term <- root_go_bound_newcol$Term %>% as.character() %>% make.unique(sep = " ")

leaf_go_bound_newcol$Term<-str_replace(leaf_go_bound_newcol$Term, "1", " ")
root_go_bound_newcol$Term<-str_replace(root_go_bound_newcol$Term, "1", " ")
leaf_go_bound_newcol$Term<-str_replace(leaf_go_bound_newcol$Term, "2", "")
root_go_bound_newcol$Term<-str_replace(root_go_bound_newcol$Term, "2", "")

root_go_bound_newcol[root_go_bound_newcol$Term == "positive regulation of transcription from RNA polymerase II promoter in response to heat stress",]$Term <- "+ve reg. transcription from RNApolII prmtr in response to heat stress"
root_go_bound_newcol[root_go_bound_newcol$Term == "positive regulation of transcription from RNA polymerase II promoter in response to heat stress  ",]$Term <- "+ve reg. transcription from RNApolI prmtr in response to heat stress"


leaf_go_bound_newcol_order <- leaf_go_bound_newcol %>%
  arrange(Species, Locality) %>%
  mutate(Direction = case_when(Direction == "up" ~ "Plastic Upregulated",
                               Direction == "down" ~"Plastic Downregulated"),
         Direction = factor(Direction, c("Plastic Upregulated", "Plastic Downregulated")),
         Term = factor(Term, levels = unique(Term)))

root_go_bound_newcol_order <- root_go_bound_newcol %>%
  arrange(Species, Locality) %>%
  mutate(Direction = case_when(Direction == "up" ~ "Plastic Upregulated",
                               Direction == "down" ~"Plastic Downregulated"),
         Direction = factor(Direction, c("Plastic Upregulated", "Plastic Downregulated")),
         Term = factor(Term, levels = unique(Term)))
 


a<-ggplot(leaf_go_bound_newcol_order, aes(x=newcol, y=Term, fill = Species, shape=Locality), colour="black") + 
  geom_point(size=30) + facet_grid(rows=vars(Direction), scales="free_y", space= "free", switch = "y", cols=vars(tissue)) + 
  theme(text = element_text(size = 77), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 61),
        #strip.text.x = element_text(size = 70),
        strip.text.x = element_blank(),
        panel.spacing=unit(1, "lines"),
        #plot.margin = margin("top", "right, "bottom", "left, "cm")) +
        plot.margin = margin(2, 2, 3, 2, "cm")) + 
  #scale_size_continuous(range = c(13, 33)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_fill_manual(values=c("gold", "dodgerblue")) +
  scale_shape_manual(values=c(21, 24)) +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(size = 29)),
         shape = guide_legend(override.aes = list(size = 29)))



b<-ggplot(root_go_bound_newcol_order, aes(x=newcol, y=Term, fill = Species, shape=Locality), colour="black") + 
  geom_point(size=30) + facet_grid(rows=vars(Direction), scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 77), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 61),
        #strip.text.x = element_text(size = 70),
        strip.text.x = element_blank(),
        legend.text=element_text(size=65),
        legend.title=element_text(size=65),
        panel.spacing=unit(1, "lines"),
        plot.margin = margin(2, 2, 3, 2, "cm")) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  theme(legend.position = "none") + 
  scale_fill_manual(values=c("gold", "dodgerblue")) +
  scale_colour_manual(values=c("gold", "dodgerblue")) +
  scale_shape_manual(values=c(21, 24)) +
  scale_y_discrete(position = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 29)),
         shape = guide_legend(override.aes = list(size = 29)))

c<-ggplot(root_go_bound_newcol_order, aes(x=newcol, y=Term, colour = Species, shape=Locality)) + 
  geom_point(size=30) + facet_grid(rows=vars(Direction), scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 40), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 63),
        #strip.text.x = element_text(size = 70),
        strip.text.x = element_blank(),
        legend.text=element_text(size=50),
        legend.title=element_text(size=47),
        panel.spacing=unit(1, "lines"),
        legend.position = c(0.45, 0.57),) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_fill_manual(values=c("gold", "dodgerblue")) +
  scale_colour_manual(values=c("gold", "dodgerblue")) +
  #scale_shape_manual(values=c(21, 24)) +
  scale_y_discrete(position = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 29)),
         shape = guide_legend(override.aes = list(size = 29)))



leg <- cowplot::get_legend(c)
leg<-as_ggplot(leg)





# the null grid here is a neat trick to make the differetn laots closer together, especially when using negative rel width and heights for it
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure6_v2.pdf", height=66, width=62, onefile=FALSE)
fst<-plot_grid(NULL, test_toplot, NULL, ncol = 3, rel_heights = c(0.01, 2.3, 0.01), rel_widths=c(0.3, 2.2, 0.3))
scd<-(plot_grid(a, b, NULL, leg,
          ncol = 4,
          rel_widths = c(2, 2, -1, 1),
          rel_heights = c(2, 2, -1, 1)))
plot_grid(fst, scd, ncol = 1, rel_heights = c(1,0.9), labels=c("A", "B"), label_size=150)
dev.off()





















#######################################
#        Plot GO term enrichment      #
#######################################


#mKLeafUp<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_up) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. majalis")
#mSLeafUp<-prepare_go_df(transplant_majalis_stulrich_leaf_up) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. majalis")
#tKLeafUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_up) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
#tSLeafUp<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_up) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. traunsteineri")
#mKLeafDown<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_down) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. majalis")
#mSLeafDown<-prepare_go_df(transplant_majalis_stulrich_leaf_down) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. majalis")
#tKLeafDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_down) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
#tSLeafDown<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_down) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. traunsteineri")
#
#mKRootUp<-prepare_go_df(transplant_majalis_kitzbuhl_root_up) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. majalis")
#mSRootUp<-prepare_go_df(transplant_majalis_stulrich_root_up) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. majalis")
##tKRootUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up) %>% mutate(comparison="D.t.Kitz", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
#tKRootUp<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
#
#tSRootUp<-prepare_go_df(transplant_traunsteineri_stulrich_root_up) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. traunsteineri")
#mKRootDown<-prepare_go_df(transplant_majalis_kitzbuhl_root_down) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. majalis")
#mSRootDown<-prepare_go_df(transplant_majalis_stulrich_root_down) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. majalis")
##tKRootDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down) %>% mutate(comparison="D.t.Kitz", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
#tKRootDown<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
#tSRootDown<-prepare_go_df(transplant_traunsteineri_stulrich_root_down) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. traunsteineri")



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




colnames(leaf_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Rich factor", "comparison", "Direction", "Locality", "Species")
colnames(root_go_bound) <-c("Term", "GO.ID", "Significant", "Annotated", "classicFisher", "Rich factor", "comparison", "Direction", "Locality", "Species")

# I am adding a new column with a placeholder the same across all rows
# I will use this column as my X axis for the GO term plots
# if I use the comparison column, it staggers the plots, and it looks better if the points are stacked
leaf_go_bound_newcol <- leaf_go_bound %>% mutate(newcol="placeholder", tissue="Leaf")
root_go_bound_newcol <- root_go_bound %>% mutate(newcol="placeholder", tissue="Root")


# I had to shorten the name of one of the facet grids (D.t. Kitz due to soace
# this alphabetically rearranges the orderr so locking in order of localitoies this way
#root_go_bound_newcol$comparison <- factor(root_go_bound$comparison, levels = unique(root_go_bound$comparison))
root_go_bound_newcol[root_go_bound_newcol$Term == "positive regulation of transcription from RNA polymerase II promoter in response to heat stress",]$Term <- "+ve reg. transcription from RNApolII prmtr in response to heat stress"
root_go_bound_newcol[root_go_bound_newcol$Term == "regulation of transcription from RNA polymerase II promoter in response to stress",]$Term <- "reg. transcription from RNApolII promoter in response to stress"
#leaf_go_bound_newcol[leaf_go_bound_newcol$Term == "peptidyl-diphthamide biosynthetic process from peptidyl-histidine",]$Term <- "peptidyl-diphthamide biosyn. proc. from peptidyl-histidine"

leaf_go_bound_newcol[leaf_go_bound_newcol$Term == "negative regulation of phosphoprotein phosphatase activity",]$Term <- "-ve regulation of phosphoprotein phosphatase activity"


leaf_go_bound_newcol$Term <- leaf_go_bound_newcol$Term %>% as.character() %>% make.unique(sep = " ")
root_go_bound_newcol$Term <- root_go_bound_newcol$Term %>% as.character() %>% make.unique(sep = " ")

leaf_go_bound_newcol$Term<-str_replace(leaf_go_bound_newcol$Term, "1", " ")
root_go_bound_newcol$Term<-str_replace(root_go_bound_newcol$Term, "1", " ")

leaf_go_bound_newcol_order <- leaf_go_bound_newcol %>%
  arrange(comparison) %>%   
  mutate(Term = factor(Term, unique(Term)),
         Direction = factor(Direction, c("Plastic Upregulated", "Plastic Downregulated")))

root_go_bound_newcol_order <- root_go_bound_newcol %>%
  arrange(comparison) %>%   
  mutate(Term = factor(Term, unique(Term)),
         Direction = factor(Direction, c("Plastic Upregulated", "Plastic Downregulated")))

a<-ggplot(leaf_go_bound_newcol_order, aes(x=newcol, y=Term, color = Species, shape=Locality, size=`Rich factor`)) + 
  geom_point() + facet_grid(rows=vars(Direction), scales="free_y", space= "free", switch = "y", cols=vars(tissue)) + 
  theme(text = element_text(size = 73), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 65),
        #strip.text.x = element_text(size = 70),
        strip.text.x = element_blank(),
        legend.text=element_text(size=40),
        legend.title=element_text(size=40),
        panel.spacing=unit(1, "lines")) + 
  scale_size_continuous(range = c(13, 33)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("gold", "dodgerblue")) +
  theme(legend.position = "none")



# I had to shorten the name of one of the facet grids (D.t. Kitz due to soace
# this alphabetically rearranges the orderr so locking in order of localitoies this way
#root_go_bound_newcol$comparison <- factor(root_go_bound$comparison, levels = unique(root_go_bound$comparison))


#root_go_bound_newcol[root_go_bound_newcol$Term == "positive regulation of transcription from RNA polymerase II promoter in response to heat stress",]$Term = "+ve reg. of transcription from RNApolII promoter in response to heat stress"
#root_go_bound_newcol[root_go_bound_newcol$Term == "regulation of transcription from RNA polymerase II promoter in response to stress" ,]$Term = "reg. of transcription from RNApolII promoter in response to stress" 


b<-ggplot(root_go_bound_newcol_order, aes(x=newcol, y=Term, color = Species, shape=Locality, size=`Rich factor`)) + 
  geom_point() + facet_grid(rows=vars(Direction), scales="free", space= "free", cols=vars(tissue)) + 
  theme(text = element_text(size = 73), 
        axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        strip.text.y = element_text(size = 70),
        #strip.text.x = element_text(size = 70),
        strip.text.x = element_blank(),
        legend.text=element_text(size=45),
        legend.title=element_text(size=45),
        panel.spacing=unit(1, "lines"),
        legend.key.size = unit(5, 'cm'), 
        legend.key.height = unit(3, 'cm'), 
        legend.key.width = unit(3, 'cm')) + 
  scale_size_continuous(range = c(13, 33)) + 
  guides(colour = guide_legend(override.aes = list(size=17))) + 
  scale_color_manual(values=c("gold", "dodgerblue")) +
  scale_y_discrete(position = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 29)),
         shape = guide_legend(override.aes = list(size = 29)))

#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure6.png", height=2800, width=4200)
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure6.pdf", height=35, width=60, onefile=FALSE)
egg::ggarrange(a, b, ncol=2)
dev.off()




###############################################################################
#     Make a table of plastic gene GO terms in full supplementary table 6     #
###############################################################################



mKLeafUp_top<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_up, top_terms = 100000) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. majalis")
mSLeafUp_top<-prepare_go_df(transplant_majalis_stulrich_leaf_up, top_terms = 100000) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. majalis")
tKLeafUp_top<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_up, top_terms = 100000) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSLeafUp_top<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_up, top_terms = 100000) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. traunsteineri")
mKLeafDown_top<-prepare_go_df(transplant_majalis_kitzbuhl_leaf_down, top_terms = 100000) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. majalis")
mSLeafDown_top<-prepare_go_df(transplant_majalis_stulrich_leaf_down, top_terms = 100000) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. majalis")
tKLeafDown_top<-prepare_go_df(transplant_traunsteineri_kitzbuhl_leaf_down, top_terms = 100000) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSLeafDown_top<-prepare_go_df(transplant_traunsteineri_stulrich_leaf_down, top_terms = 100000) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. traunsteineri")

mKRootUp_top<-prepare_go_df(transplant_majalis_kitzbuhl_root_up, top_terms = 100000) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. majalis")
mSRootUp_top<-prepare_go_df(transplant_majalis_stulrich_root_up, top_terms = 100000) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. majalis")
tKRootUp_top<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_up, top_terms = 100000) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Upregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSRootUp_top<-prepare_go_df(transplant_traunsteineri_stulrich_root_up, top_terms = 100000) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Upregulated", locality="St. Ulrich", species="D. traunsteineri")
mKRootDown_top<-prepare_go_df(transplant_majalis_kitzbuhl_root_down, top_terms = 100000) %>% mutate(comparison="D. majalis Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. majalis")
mSRootDown_top<-prepare_go_df(transplant_majalis_stulrich_root_down, top_terms = 100000) %>% mutate(comparison="D. majalis St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. majalis")
tKRootDown_top<-prepare_go_df(transplant_traunsteineri_kitzbuhl_root_down, top_terms = 100000) %>% mutate(comparison="D. traunsteineri Kitzbuhel", Direction="Plastic Downregulated", locality="Kitzbuhel", species="D. traunsteineri")
tSRootDown_top<-prepare_go_df(transplant_traunsteineri_stulrich_root_down, top_terms = 100000) %>% mutate(comparison="D. traunsteineri St. Ulrich", Direction="Plastic Downregulated", locality="St. Ulrich", species="D. traunsteineri")




rbind(mKLeafUp_top,
                     mSLeafUp_top,
                     tKLeafUp_top,
                     tSLeafUp_top,
                     mKLeafDown_top,
                     mSLeafDown_top,
                     tKLeafDown_top,
                     tSLeafDown_top,
                    mKRootUp_top,
                    mSRootUp_top,
                    tKRootUp_top,
                    tSRootUp_top,
                    mKRootDown_top,
                    mSRootDown_top,
                    tKRootDown_top,
                    tSRootDown_top) %>% write.xlsx(file = "SupplementaryTable6_plastic_enriched_GOterms.xlsx")













packageDescription("DeSeq2")$Version













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



#library(GSEABase)
#fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
#äslim <- getOBOCollection(fl)
#myIDs<-transplant_majalis_kitzbuhl_root_up$GO.ID
#myCollection <- GOCollection(myIDs)
#goSlim(myCollection, slim, "MF")





#ensembl = useMart("plants_mart", host="plants.ensembl.org")
#ensembl = useDataset("athaliana_eg_gene", mart=ensembl)
#query = getBM(attributes=c("go_id", "name_1006"), values=test$Representative, mart=ensembl) %>% set_colnames(c("Representative", "name_1006"))

