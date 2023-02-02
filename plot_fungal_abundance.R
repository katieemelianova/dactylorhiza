

library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(viridis)


read_tax<-function(file, sampleid, rank){
  exam<-read_tsv(file) %>% mutate(sample=sampleid) %>% filter(RANK == rank)
  exam %<>% separate(TAXPATHSN, c(NA, NA, NA, "B", NA, NA))
  return(exam)
}



mTS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS2r_abundances.txt", "mTS2r", "order")
tMK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK2r_abundances.txt", "tMK2r", "order")
mMK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK5r_abundances.txt", "mMK5r", "order")
tMS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS1r_abundances.txt", "tMS1r", "order")
mTK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK1r_abundances.txt", "mTK1r", "order")
tMK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK3r_abundances.txt", "tMK3r", "order")
mTS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS3r_abundances.txt", "mTS3r", "order")
tTS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS4r_abundances.txt", "tTS4r", "order")
mMK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK4r_abundances.txt", "mMK4r", "order")
tMS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS3r_abundances.txt", "tMS3r", "order")
mTK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK3r_abundances.txt", "mTK3r", "order")
tTK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK4r_abundances.txt", "tTK4r", "order")
mMS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS4r_abundances.txt", "mMS4r", "order")
mTK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK2r_abundances.txt", "mTK2r", "order")
tMS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS2r_abundances.txt", "tMS2r", "order")
tTK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK5r_abundances.txt", "tTK5r", "order")
tMK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK1r_abundances.txt", "tMK1r", "order")
mTS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS1r_abundances.txt", "mTS1r", "order")
tMS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMS4r_abundances.txt", "tMS4r", "order")
mTK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK4r_abundances.txt", "mTK4r", "order")
tTK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK3r_abundances.txt", "tTK3r", "order")
mMS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS3r_abundances.txt", "mMS3r", "order")
tTS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS1r_abundances.txt", "tTS1r", "order")
mMK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK1r_abundances.txt", "mMK1r", "order")
mTK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTK5r_abundances.txt", "mTK5r", "order")
mMS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS2r_abundances.txt", "mMS2r", "order")
tTK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK2r_abundances.txt", "tTK2r", "order")
tTK1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTK1r_abundances.txt", "tTK1r", "order")
mMS1r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMS1r_abundances.txt", "mMS1r", "order")
tMK5r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK5r_abundances.txt", "tMK5r", "order")
mMK2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK2r_abundances.txt", "mMK2r", "order")
tTS2r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS2r_abundances.txt", "tTS2r", "order")
tMK4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tMK4r_abundances.txt", "tMK4r", "order")
mTS4r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mTS4r_abundances.txt", "mTS4r", "order")
tTS3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/tTS3r_abundances.txt", "tTS3r", "order")
mMK3r_data<-read_tax("/Users/katieemelianova/Desktop/Dactylorhiza/abundances/mMK3r_abundances.txt", "mMK3r", "order")






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
  dplyr::select(B, PERCENTAGE, sample) %>%
  filter(PERCENTAGE > 5) %>%
  melt()



ggplot(data=all, aes(x=sample, y=value, fill=B)) +
  geom_bar(stat="identity")


