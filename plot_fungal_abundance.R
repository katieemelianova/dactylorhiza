

library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)


read_tax<-function(file, sampleid, rank, full_taxonomy=FALSE){
  exam<-read_tsv(file) %>% mutate(sample=sampleid) %>% filter(RANK == rank)
  if(full_taxonomy == FALSE){ 
    exam %<>% separate(TAXPATHSN, c(NA, NA, NA, NA, "B", NA))
  }
  return(exam)
}



mTS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS2r_abundances.txt", "mTS2r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tMK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK2r_abundances.txt", "tMK2r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK5r_abundances.txt", "mMK5r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
tMS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS1r_abundances.txt", "tMS1r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mTK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK1r_abundances.txt", "mTK1r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tMK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK3r_abundances.txt", "tMK3r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mTS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS3r_abundances.txt", "mTS3r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tTS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS4r_abundances.txt", "tTS4r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK4r_abundances.txt", "mMK4r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
tMS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS3r_abundances.txt", "tMS3r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mTK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK3r_abundances.txt", "mTK3r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tTK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK4r_abundances.txt", "tTK4r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS4r_abundances.txt", "mMS4r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
mTK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK2r_abundances.txt", "mTK2r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tMS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS2r_abundances.txt", "tMS2r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
tTK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK5r_abundances.txt", "tTK5r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
tMK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK1r_abundances.txt", "tMK1r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mTS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS1r_abundances.txt", "mTS1r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tMS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS4r_abundances.txt", "tMS4r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mTK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK4r_abundances.txt", "mTK4r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tTK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK3r_abundances.txt", "tTK3r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS3r_abundances.txt", "mMS3r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
tTS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS1r_abundances.txt", "tTS1r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK1r_abundances.txt", "mMK1r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
mTK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK5r_abundances.txt", "mTK5r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
mMS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS2r_abundances.txt", "mMS2r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
tTK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK2r_abundances.txt", "tTK2r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
tTK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK1r_abundances.txt", "tTK1r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS1r_abundances.txt", "mMS1r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
tMK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK5r_abundances.txt", "tMK5r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mMK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK2r_abundances.txt", "mMK2r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")
tTS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS2r_abundances.txt", "tTS2r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
tMK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK4r_abundances.txt", "tMK4r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. majalis environment")
mTS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS4r_abundances.txt", "mTS4r", "family") %>% mutate(Species="D. majalis", Environment="D. traunsteineri environment")
tTS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS3r_abundances.txt", "tTS3r", "family") %>% mutate(Species="D. traunsteineri", Environment="D. traunsteineri environment")
mMK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK3r_abundances.txt", "mMK3r", "family") %>% mutate(Species="D. majalis", Environment="D. majalis environment")






all<-rbind(mTS2r_data,
           tMK2r_data,
           mMK5r_data,
           tMS1r_data,
           mTK1r_data,
           tMK3r_data,
           mTS3r_data,
           tTS4r_data,
           mMK4r_data,
           tMS3r_data,
           mTK3r_data,
           tTK4r_data,
           mMS4r_data,
           mTK2r_data,
           tMS2r_data,
           tTK5r_data,
           tMK1r_data,
           mTS1r_data,
           tMS4r_data,
           mTK4r_data,
           tTK3r_data,
           mMS3r_data,
           tTS1r_data,
           mMK1r_data,
           mTK5r_data,
           mMS2r_data,
           tTK2r_data,
           tTK1r_data,
           mMS1r_data,
           tMK5r_data,
           mMK2r_data,
           tTS2r_data,
           tMK4r_data,
           mTS4r_data,
           tTS3r_data,
           mMK3r_data
           ) %>% 
  filter(B != "Unclassified") %>%
  dplyr::select(B, PERCENTAGE, sample, Species, Environment) %>%
  filter(PERCENTAGE > 5) %>%
  melt()

##################################################
#            plot abundances barplot             #
##################################################

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure9.png", width = 1100, height = 500)
ggplot(data=all, aes(x=sample, y=value, fill=B)) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(
          size = 13,
          face="bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.text=element_text(size=13),
        legend.title=element_text(size=15, face="bold"),
        axis.title.y = element_text(size = 15)) +
  ylab("% Abundance") +
  labs(fill = "Fungal order") +
  facet_wrap(~ Species + Environment, scales = "free", ncol=4)+
  scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, 0)))
dev.off()








#all_temp<-all_5pct %>% mutate(V5 = case_when(species == "majalis" ~ value*(-1),
#                              species == "traunsteineri" ~ value
#                      ))
#brks <- seq(-15000000, 15000000, 5000000)
#lbls = paste0(as.character(c(seq(15, 0, -5), seq(5, 15, 5))), "m")
#
#ggplot(all_temp, aes(x = B, y = V5, fill = species)) +
#  geom_bar(stat = "identity", width = .6) +   # draw the bars
#  scale_y_continuous(breaks = brks,
#                     labels = lbls) +
#  coord_flip() +  # Flip axes
#  scale_fill_brewer(palette = "Accent") +  # Color palette
#  facet_wrap(~environment)


all$species<-case_when(substr(all$sample,1,1) == "m" ~ "majalis",
                       substr(all$sample,1,1) == "t" ~ "traunsteineri")
all$locality<-case_when(substr(all$sample,3,3) == "S" ~ "St Ulrich",
                        substr(all$sample,3,3) == "K" ~ "Kitzbuhl")
all$environment<-case_when(substr(all$sample,2,2) == "M" ~ "majalis",
                           substr(all$sample,2,2) == "T" ~ "traunsteineri")



env_labs <- c("D. traunsteineri environment", "D. majalis environment")
names(env_labs) <- c("traunsteineri", "majalis")
sp_labs <- c("D. traunsteineri", "D. majalis")
names(sp_labs) <- c("traunsteineri", "majalis")

#pie_all<-aggregate(all$value, by=list(Category=all$B), FUN=sum)
pie_all<-aggregate(all$value ~all$B+species+locality+environment, all, FUN=sum)
colnames(pie_all)<-c("B", "species", "locality", "environment", "value")


qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) %in% c("Dark2", "Spectral", "Paired", "PuOr"),]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure10.png", width = 900, height = 600)
ggplot(pie_all, aes(x="", y=value, fill=B)) +
  geom_bar(stat="identity", width=1, position = "fill") +
  coord_polar("y", start=0) +
  facet_wrap(~environment + species + locality, labeller = ggplot2::labeller(environment = env_labs, species = sp_labs), ncol=4) +
  theme_void() +
  scale_fill_manual(values=col_vector[1:38], name = "Fungal family") +
  theme(strip.text.x = element_text(size = 13, face="bold"), 
        axis.text.x=element_blank(),
        legend.text=element_text(size=13),
        legend.title=element_text(size=15, face="bold"))
dev.off()

#############################################################################
#     statistical tests (repetitive code but needs to be a bit formatted)   #
#############################################################################

mTS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS2r_abundances.txt", "mTS2r", "familyy", full_taxonomy=TRUE)
tMK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK2r_abundances.txt", "tMK2r", "family", full_taxonomy=TRUE)
mMK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK5r_abundances.txt", "mMK5r", "family", full_taxonomy=TRUE)
tMS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS1r_abundances.txt", "tMS1r", "family", full_taxonomy=TRUE)
mTK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK1r_abundances.txt", "mTK1r", "family", full_taxonomy=TRUE)
tMK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK3r_abundances.txt", "tMK3r", "family", full_taxonomy=TRUE)
mTS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS3r_abundances.txt", "mTS3r", "family", full_taxonomy=TRUE)
tTS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS4r_abundances.txt", "tTS4r", "family", full_taxonomy=TRUE)
mMK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK4r_abundances.txt", "mMK4r", "family", full_taxonomy=TRUE)
tMS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS3r_abundances.txt", "tMS3r", "family", full_taxonomy=TRUE)
mTK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK3r_abundances.txt", "mTK3r", "family", full_taxonomy=TRUE)
tTK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK4r_abundances.txt", "tTK4r", "family", full_taxonomy=TRUE)
mMS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS4r_abundances.txt", "mMS4r", "family", full_taxonomy=TRUE)
mTK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK2r_abundances.txt", "mTK2r", "family", full_taxonomy=TRUE)
tMS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS2r_abundances.txt", "tMS2r", "family", full_taxonomy=TRUE)
tTK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK5r_abundances.txt", "tTK5r", "family", full_taxonomy=TRUE)
tMK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK1r_abundances.txt", "tMK1r", "family", full_taxonomy=TRUE)
mTS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS1r_abundances.txt", "mTS1r", "family", full_taxonomy=TRUE)
tMS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS4r_abundances.txt", "tMS4r", "family", full_taxonomy=TRUE)
mTK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK4r_abundances.txt", "mTK4r", "family", full_taxonomy=TRUE)
tTK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK3r_abundances.txt", "tTK3r", "family", full_taxonomy=TRUE)
mMS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS3r_abundances.txt", "mMS3r", "family", full_taxonomy=TRUE)
tTS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS1r_abundances.txt", "tTS1r", "family", full_taxonomy=TRUE)
mMK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK1r_abundances.txt", "mMK1r", "family", full_taxonomy=TRUE)
mTK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK5r_abundances.txt", "mTK5r", "family", full_taxonomy=TRUE)
mMS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS2r_abundances.txt", "mMS2r", "family", full_taxonomy=TRUE)
tTK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK2r_abundances.txt", "tTK2r", "family", full_taxonomy=TRUE)
tTK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK1r_abundances.txt", "tTK1r", "family", full_taxonomy=TRUE)
mMS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS1r_abundances.txt", "mMS1r", "family", full_taxonomy=TRUE)
tMK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK5r_abundances.txt", "tMK5r", "family", full_taxonomy=TRUE)
mMK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK2r_abundances.txt", "mMK2r", "family", full_taxonomy=TRUE)
tTS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS2r_abundances.txt", "tTS2r", "family", full_taxonomy=TRUE)
tMK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK4r_abundances.txt", "tMK4r", "family", full_taxonomy=TRUE)
mTS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS4r_abundances.txt", "mTS4r", "family", full_taxonomy=TRUE)
tTS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS3r_abundances.txt", "tTS3r", "family", full_taxonomy=TRUE)
mMK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK3r_abundances.txt", "mMK3r", "family", full_taxonomy=TRUE)

all_stat<-rbind(mTS2r_data,
                tMK2r_data,
                mMK5r_data,
                tMS1r_data,
                mTK1r_data,
                tMK3r_data,
                mTS3r_data,
                tTS4r_data,
                mMK4r_data,
                tMS3r_data,
                mTK3r_data,
                tTK4r_data,
                mMS4r_data,
                mTK2r_data,
                tMS2r_data,
                tTK5r_data,
                tMK1r_data,
                mTS1r_data,
                tMS4r_data,
                mTK4r_data,
                tTK3r_data,
                mMS3r_data,
                tTS1r_data,
                mMK1r_data,
                mTK5r_data,
                mMS2r_data,
                tTK2r_data,
                tTK1r_data,
                mMS1r_data,
                tMK5r_data,
                mMK2r_data,
                tTS2r_data,
                tMK4r_data,
                mTS4r_data,
                tTS3r_data,
                mMK3r_data
                )

all_stat$species<-dplyr::case_when(substr(all_stat$sample,1,1) == "m" ~ "majalis",
                       substr(all_stat$sample,1,1) == "t" ~ "traunsteineri")
all_stat$locality<-case_when(substr(all_stat$sample,3,3) == "S" ~ "St Ulrich",
                        substr(all_stat$sample,3,3) == "K" ~ "Kitzbuhl")
all_stat$environment<-case_when(substr(all_stat$sample,2,2) == "M" ~ "majalis",
                           substr(all_stat$sample,2,2) == "T" ~ "traunsteineri")


unclassified<-all_stat$TAXPATHSN[endsWith(all_stat$TAXPATHSN, "Unclassified")]

all_stat %<>% filter(RANK == "family" & sample != "tMK4r" & !(TAXPATHSN %in% unclassified)) %>% dplyr::select(TAXPATHSN, 
                                                       PERCENTAGE,
                                                       sample
                                                       )

all_stat_spread <- spread(all_stat, sample, PERCENTAGE) %>% filter %>% data.frame()

species_line<-case_when(substr(colnames(all_stat_spread),1,1) == "m" ~ "species_majalis",
          substr(colnames(all_stat_spread),1,1) == "t" ~ "speciestraunsteineri")
env_line<-case_when(substr(colnames(all_stat_spread),2,2) == "M" ~ "environment_majalis",
          substr(colnames(all_stat_spread),2,2) == "T" ~ "environment_traunsteineri")

all_stat_spread

species_line[1]<-"Species"
env_line[1]<-"Environemt"
colnames(all_stat_spread)[1]<-"Sample ID"

all_stat_spread_annotated<-rbind(species_line, env_line, all_stat_spread) %>% na.omit()



write_tsv(all_stat_spread_annotated, "/Users/katieemelianova/Desktop/Dactylorhiza/dactylorhiza/all_stat_spread_annotated.tsv")

#######################################################################
#    After this I uploaded the dataset to the LEFSE analysis here:    #
#            http://huttenhower.sph.harvard.edu/galaxy/               #
#######################################################################


