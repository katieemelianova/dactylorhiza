# plotting sargasso loci mapped after partitioning
setwd("/Users/katieemelianova/Desktop/Dactylorhiza/dactylorhiza/loci_counts")

traunsteineri_fuchsii_min1<-read_tsv("traunsteineri_fuchsii_min1", col_names = "mapped_reads") %>% mutate(species="traunsteineri", reference="fuchsii", minimum_mapped="1")
traunsteineri_fuchsii_min10<-read_tsv("traunsteineri_fuchsii_min10", col_names = "mapped_reads") %>% mutate(species="traunsteineri", reference="fuchsii", minimum_mapped="10")
traunsteineri_fuchsii_min100<-read_tsv("traunsteineri_fuchsii_min100", col_names = "mapped_reads") %>% mutate(species="traunsteineri", reference="fuchsii", minimum_mapped="100")

traunsteineri_incarnata_min1<-read_tsv("traunsteineri_incarnata_min1", col_names = "mapped_reads") %>% mutate(species="traunsteineri", reference="incarnata", minimum_mapped="1")
traunsteineri_incarnata_min10<-read_tsv("traunsteineri_incarnata_min10", col_names = "mapped_reads") %>% mutate(species="traunsteineri", reference="incarnata", minimum_mapped="10")
traunsteineri_incarnata_min100<-read_tsv("traunsteineri_incarnata_min100", col_names = "mapped_reads") %>% mutate(species="traunsteineri", reference="incarnata", minimum_mapped="100")

majalis_fuchsii_min1<-read_tsv("majalis_fuchsii_min1", col_names = "mapped_reads") %>% mutate(species="majalis", reference="fuchsii", minimum_mapped="1")
majalis_fuchsii_min10<-read_tsv("majalis_fuchsii_min10", col_names = "mapped_reads") %>% mutate(species="majalis", reference="fuchsii", minimum_mapped="10")
majalis_fuchsii_min100<-read_tsv("majalis_fuchsii_min100", col_names = "mapped_reads") %>% mutate(species="majalis", reference="fuchsii", minimum_mapped="100")

majalis_incarnata_min1<-read_tsv("majalis_incarnata_min1", col_names = "mapped_reads") %>% mutate(species="majalis", reference="incarnata", minimum_mapped="1")
majalis_incarnata_min10<-read_tsv("majalis_incarnata_min10", col_names = "mapped_reads") %>% mutate(species="majalis", reference="incarnata", minimum_mapped="10")
majalis_incarnata_min100<-read_tsv("majalis_incarnata_min100", col_names = "mapped_reads") %>% mutate(species="majalis", reference="incarnata", minimum_mapped="100")



bound<-rbind(traunsteineri_fuchsii_min1,
             traunsteineri_fuchsii_min10,
             traunsteineri_fuchsii_min100,
             traunsteineri_incarnata_min1,
             traunsteineri_incarnata_min10,
             traunsteineri_incarnata_min100,
             majalis_fuchsii_min1,
             majalis_fuchsii_min10,
             majalis_fuchsii_min100,
             majalis_incarnata_min1,
             majalis_incarnata_min10,
             majalis_incarnata_min100) %>% melt()                                                                                               



png("root_sargasso_mapped_loci_boxplot.png", height = 900, width=1300)
ggplot(bound, aes(x=minimum_mapped, y=value, fill=reference)) + 
  geom_boxplot() +
  xlab("Minimum mapped reads (root)") +
  ylab("Number of loci (root)") +
  scale_fill_brewer(palette = "Accent") +
  facet_wrap(~species) +
  theme(text = element_text(size = 28),
        strip.text.x = element_text(size = 30),
        legend.text=element_text(size=30))
dev.off()
