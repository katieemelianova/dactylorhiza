

library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(ggh4x)
library(Polychrome)
library(ggpubr)
library(vegan)


read_tax<-function(file, sampleid, rank, full_taxonomy=FALSE){
  exam<-read_tsv(file) %>% mutate(sample=sampleid) %>% filter(RANK == rank)
  if(full_taxonomy == FALSE){ 
    exam %<>% separate(TAXPATHSN, c(NA, NA, NA, "B", NA, NA))
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




samples_to_exclude<-c("mTS2r", "tTS4r", "tTK1r", "tMK4r")



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
  #filter(B != "Unclassified" & !(sample %in% samples_to_exclude)) %>%
  filter(!(sample %in% samples_to_exclude)) %>%
  dplyr::select(B, PERCENTAGE, sample, Species, Environment) %>%
  mutate(B2 = case_when(B == "Unclassified" ~ "Other",
                        PERCENTAGE >= 1 ~ B,
                        PERCENTAGE <= 1 ~ "Other",
                        )) %>%
  melt()

# rename enviorinment to be shorter
all %<>% mutate(Environment = ifelse(Environment == "D. majalis environment", "M", "T"))
# order orders, put other at the end
order_of_orders=(c(unique(all$B2)[which(unique(all$B2) != "Other")], "Other"))
all %<>% 
  arrange(B2) %>% 
  mutate(B2 = factor(B2, levels=order_of_orders))





#####################################################################
#     define matching colour scheme for bar and pie charts          #
#####################################################################


# I used these lines to pick seed colours for 
#qual_col_pals = brewer.pal.info[rownames(brewer.pal.info) %in% c("Spectral"),]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

P36 = createPalette(36,  c("#5E4FA2", "#3288BD", "#66C2A5")) %>% as.character()

# get as many colours as you need, minus one to put grey at the end for Other
col_vector<-c(P36[1:length(levels(all$B2)) -1], "gray63")




##################################################
#            plot abundances barplot             #
##################################################

strip <- strip_themed(background_x = elem_list_rect(fill = c("gold", "gold", "dodgerblue", "dodgerblue", "gold", "dodgerblue", "gold", "dodgerblue")))

#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure7.png", width = 1600, height = 700)
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure7.pdf", height=20, width=45)
ggplot(data=all, aes(x=sample, y=value, fill=B2), colour=B2) +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(
          size = 50),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.text=element_text(size=50),
        legend.title=element_text(size=60, face="bold"),
        axis.title.y = element_text(size = 65)) +
  ylab("% Abundance") +
  labs(fill = "Fungal order") +
  facet_wrap2(~ Species + Environment, scales = "free", strip = strip, ncol=4) +
  #facet_wrap(~ Species + Environment, scales = "free", ncol=4)+
  #scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values=col_vector[1:20], name = "Fungal order")
dev.off()





# calculation of average non-Wallemiales taxa
all %>% filter(Environment == "D. majalis environment" & Species == "D. majalis" & B != "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()
all %>% filter(Environment == "D. traunsteineri environment" & Species == "D. majalis" & B != "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()
all %>% filter(Environment == "D. majalis environment" & Species == "D. traunsteineri" & B != "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()
all %>% filter(Environment == "D. traunsteineri environment" & Species == "D. traunsteineri" & B != "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()


# calculation of how much wallemiales in each comparison
all %>% filter(Environment == "D. majalis environment" & Species == "D. majalis" & value > 1 & B == "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()
all %>% filter(Environment == "D. traunsteineri environment" & Species == "D. majalis" & value > 1 & B == "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()
all %>% filter(Environment == "D. majalis environment" & Species == "D. traunsteineri" & value > 1 & B == "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()
all %>% filter(Environment == "D. traunsteineri environment" & Species == "D. traunsteineri" & value > 1 & B == "Wallemiales") %>% dplyr::select(c("sample", "value")) %>% aggregate(.$value ~.$sample, ., FUN=sum) %>% set_colnames(c("sample", "value")) %>% dplyr::select("value") %>% pull() %>% mean()






##################################################################
#                    Shannon Diversity index                     #
##################################################################



BCI

diversity(mM_K)

##################
#     majalis    #
##################

mM_K<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "majalis" & locality == "Kitzbuhl" & environment == "majalis") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")


mT_K<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "majalis" & locality == "Kitzbuhl" & environment == "traunsteineri") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")


mM_S<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "majalis" & locality == "St Ulrich" & environment == "majalis") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")

mT_S<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "majalis" & locality == "St Ulrich" & environment == "traunsteineri") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")


##################
#  traunsteineri #
##################

tT_K<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "traunsteineri" & locality == "Kitzbuhl" & environment == "traunsteineri") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")

tM_K<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "traunsteineri" & locality == "Kitzbuhl" & environment == "majalis") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")

tT_S<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "traunsteineri" & locality == "St Ulrich" & environment == "traunsteineri") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")

tM_S<-all %>% 
  dplyr::select(sample, B2, value, species, locality, environment) %>% 
  filter(species == "traunsteineri" & locality == "St Ulrich" & environment == "majalis") %>% 
  group_by(B2, sample) %>% 
  summarise(mean(value)) %>% 
  set_colnames(c("Order", "Sample", "Mean")) %>% 
  pivot_wider(names_from = Order, values_from = Mean) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("Sample")


shannon_boxplot<-rbind(diversity(mM_K) %>% data.frame() %>% mutate(Species="D. majalis", Environment="M", Locality="Kitzbuhel") %>% rename("Shannon Index"="."),
  diversity(mT_K) %>% data.frame() %>% mutate(Species="D. majalis", Environment="T", Locality="Kitzbuhel") %>% rename("Shannon Index"="."),
  diversity(mM_S) %>% data.frame() %>% mutate(Species="D. majalis", Environment="M", Locality="St. Ulrich") %>% rename("Shannon Index"="."),
  diversity(mT_S) %>% data.frame() %>% mutate(Species="D. majalis", Environment="T", Locality="St. Ulrich") %>% rename("Shannon Index"="."),
  diversity(tT_K) %>% data.frame() %>% mutate(Species="D. traunsteineri", Environment="T", Locality="Kitzbuhel") %>% rename("Shannon Index"="."),
  diversity(tM_K) %>% data.frame() %>% mutate(Species="D. traunsteineri", Environment="M", Locality="Kitzbuhel") %>% rename("Shannon Index"="."),
  diversity(tT_S) %>% data.frame() %>% mutate(Species="D. traunsteineri", Environment="T", Locality="St. Ulrich") %>% rename("Shannon Index"="."),
  diversity(tM_S) %>% data.frame() %>% mutate(Species="D. traunsteineri", Environment="M", Locality="St. Ulrich") %>% rename("Shannon Index"="."))

jitter_colour<-ifelse(shannon_boxplot$Species == "D. majalis", "gold", "dodgerblue")

pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure8_v3.pdf", height=22, width=20)
ggplot(shannon_boxplot, aes(x=Environment, y=`Shannon Index`, fill=Environment)) + 
  geom_boxplot(outlier.shape =NA, lwd=3, key_glyph = "rect") + 
  geom_jitter(shape=21, size=13, position=position_jitter(0.2), fill=jitter_colour, colour="black", stroke=5) +
  facet_wrap(~Locality) +
  scale_fill_manual(values=c("gold", "dodgerblue")) +
  theme(legend.position = c(0.8, 0.9),
        #legend.position = "none",
        text = element_text(size = 70),
        legend.background=element_blank(),
        legend.title=element_blank()) +
  scale_fill_discrete(labels=c('D. majalis', 'D. traunsteineri'), type=c("gold", "dodgerblue"))
dev.off()

# check if there is any difference in standard deviation for species in home vs away fungal diversity
shannon_boxplot %>% group_by(Species, Environment, Locality) %>% summarise_at(c("Shannon Index"), mean) 

###########################################################################################
#       t-test to show statistical difference between M and T in both K and S             #
###########################################################################################

T_environment_K<-shannon_boxplot %>% filter(Environment == "T" & Locality == "Kitzbuhel") %>% pull("Shannon Index")
T_environment_S<-shannon_boxplot %>% filter(Environment == "T" & Locality == "St. Ulrich") %>% pull("Shannon Index")
M_environment_K<-shannon_boxplot %>% filter(Environment == "M" & Locality == "Kitzbuhel") %>% pull("Shannon Index")
M_environment_S<-shannon_boxplot %>% filter(Environment == "M" & Locality == "St. Ulrich") %>% pull("Shannon Index")

res_K <- t.test(T_environment_K, M_environment_K)
res_S <- t.test(M_environment_S, T_environment_S)

T_environment<-shannon_boxplot %>% filter(Environment == "T") %>% pull("Shannon Index")
M_environment<-shannon_boxplot %>% filter(Environment == "M") %>% pull("Shannon Index")

res<-t.test(T_environment, M_environment)



##################################################################
#               Pie charts (might not include these)             #
##################################################################



all$species<-case_when(substr(all$sample,1,1) == "m" ~ "majalis",
                       substr(all$sample,1,1) == "t" ~ "traunsteineri")
all$locality<-case_when(substr(all$sample,3,3) == "S" ~ "St Ulrich",
                        substr(all$sample,3,3) == "K" ~ "Kitzbuhl")
all$environment<-case_when(substr(all$sample,2,2) == "M" ~ "majalis",
                           substr(all$sample,2,2) == "T" ~ "traunsteineri")



env_labs <- c("D. traunsteineri environment", "D. majalis environment")
names(env_labs) <- c("traunsteineri", "majalis")

env_labs_blank <- c("", "")
names(env_labs_blank) <- c("traunsteineri", "majalis")

sp_labs <- c("D. traunsteineri", "D. majalis")
names(sp_labs) <- c("traunsteineri", "majalis")

pie_all<-aggregate(all$value ~all$B2+species+locality+environment, all, FUN=mean)
colnames(pie_all)<-c("B2", "species", "locality", "environment", "value")

maj_pie<-pie_all %>% filter(environment == "majalis") %>% 
  ggplot(aes(x="", y=value, fill=B2)) +
  geom_bar(stat="identity", width=1, position = "fill") +
  coord_polar("y", start=0) +
  facet_wrap(~environment + species + locality, labeller = ggplot2::labeller(environment = env_labs_blank, species = sp_labs), ncol=2) +
  theme_void() +
  scale_fill_manual(values=col_vector[1:20], name = "Fungal order") +
  theme(#strip.text.x = element_blank(),
        strip.text.x = element_text(size = 13), 
        legend.position = "none",
        axis.text.x=element_blank(),
        legend.text=element_text(size=13),
        legend.title=element_text(size=15, face="bold"),
        plot.background = element_rect(fill=alpha('gold', 0.2), colour="white"))



tra_pie<-pie_all %>% filter(environment == "traunsteineri") %>% 
  ggplot(aes(x="", y=value, fill=B2)) +
  geom_bar(stat="identity", width=1, position = "fill") +
  coord_polar("y", start=0) +
  facet_wrap(~environment + species + locality, labeller = ggplot2::labeller(environment = env_labs_blank, species = sp_labs), ncol=2) +
  theme_void() +
  scale_fill_manual(values=col_vector[1:20], name = "Fungal order") +
  theme(strip.text.x = element_text(size = 13), 
        legend.position = "none",
        axis.text.x=element_blank(),
        legend.text=element_text(size=13),
        legend.title=element_text(size=15, face="bold"),
        plot.background = element_rect(fill=alpha('dodgerblue', 0.2), colour="white"))


legend_plot<-ggplot(pie_all, aes(x="", y=value, fill=B2)) +
  geom_bar(stat="identity", width=1, position = "fill") +
  coord_polar("y", start=0) +
  facet_wrap(~environment + species + locality, labeller = ggplot2::labeller(environment = env_labs, species = sp_labs), ncol=4) +
  theme_void() +
  scale_fill_manual(values=col_vector[1:20], name = "Fungal order") +
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15, face="bold"))

leg <- cowplot::get_legend(legend_plot)
leg<-as_ggplot(leg)


grid_a<-cowplot::plot_grid(maj_pie, tra_pie, ncol = 1)
#png(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure8_v2.png", width = 800, height = 2000)
pdf(file="~/Desktop/Dactylorhiza/dactylorhiza/Figure8_v2.pdf", height=30, width=30)
cowplot::plot_grid(grid_a, NULL, 
                   leg, ncol=3, 
                   rel_widths = c(3, -1, 1),
                   rel_heights = c(3, -1, 1))
dev.off()




















pie_all %>% filter(environment == "traunsteineri" & B2 == "Spizellomycetales")
pie_all %>% filter(environment == "majalis" & B2 == "Spizellomycetales")


pie_all %>% filter(environment == "traunsteineri" & B2 == "Wallemiales")
pie_all %>% filter(environment == "majalis" & B2 == "Wallemiales")

pie_all %>% filter(environment == "traunsteineri" & B2 == "Mortierellales")
pie_all %>% filter(environment == "majalis" & B2 == "Mortierellales")


Metschnikowiaceae<-c("4.95051140598397", "4.14190720823853", "3.75429195682506", "4.0666008024767", "7.79768874779775", "3.13711986917807", "4.25065818326258", "3.19785504237587", "2.15101319638948", "1.79415277851309", "11.4381103015332", "4.08167398824527", "5.71679467233381", "5.77298416061022", "5.84661267855805", "3.62544403139846", "3.40470192870115", "3.0021293822947", "5.65171167602384", "13.3713331029679", "4.60768721930311", "10.6179750634019", "5.11088233317787", "7.16394204051486", "2.0412422486734", "5.83200730326674", "7.97951064993315", "9.45357383821955", "2.11726229298103", "89.4459467113006", "5.03049809433911", "4.20752980303311", "6.91381717835442")
Wallemiacetes<-c("59.1668055600235", "36.215710507056", "91.6381422089486", "37.7644107085191", "60.2817218353655", "29.8170541231346", "52.1817920427708", "89.6808540756771", "5.09064911206817", "4.80126153808208", "40.1337025619341", "0.647357012426202", "90.157455908127", "87.9530483735043", "63.1897303489549", "17.4025675282367", "89.8215713061064", "17.5676206314639", "64.3362659956943", "35.3479484510778", "46.4053243706007", "22.459840273908 13.4983058215499", "68.5944663413997", "13.762026347436 87.5202440883756", "77.3058008512956", "54.7379409580447", "22.0912617026523", "0.313134018027669", "83.0301613666154", "55.8514507087588", "81.4127642328277")

Wallemiacetes<-c("59.1668055600235", "36.215710507056 91.6381422089486", "37.7644107085191", "60.2817218353655", "29.8170541231346", "52.1817920427708", "89.6808540756771", "5.09064911206817", "4.80126153808208", "40.1337025619341", "0.647357012426202", "90.157455908127 87.9530483735043", "63.1897303489549", "17.4025675282367", "89.8215713061064", "17.5676206314639", "64.3362659956943", "35.3479484510778", "46.4053243706007", "22.459840273908 13.4983058215499", "68.5944663413997", "13.762026347436 87.5202440883756", "77.3058008512956", "54.7379409580447", "22.0912617026523", "0.313134018027669", "83.0301613666154", "55.8514507087588", "81.4127642328277")

Metschnikowiaceae %>% length()

en<-c("M", "M", "M", "T", "M", "T", "M", "M", "T", "T", "M", "T", "M", "T", "M", "T", "M", "T", "T", "M", "T", "M", "T", "M", "T", "M", "M", "M", "T", "M", "T", "T", "M")
test<-data.frame(m=Metschnikowiaceae, w=wall, e=en)

test%>% filter(e=="T") %>% dplyr::select(m) %>% pull() %>% as.numeric() %>% mean()
test%>% filter(e=="T") %>% dplyr::select(w) %>% pull() %>% as.numeric() %>% mean()


#########################################################
#                      LEFSE input                      #
#########################################################

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


