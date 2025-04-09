#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#read in metadata and asv table
asv_table <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt", header=T)

#remove contams
contam<-which(asv_table$PCRCtrl>0)
length(contam)
#105

asv_table2<-asv_table[-contam,]
asv_table2<-asv_table2[,-54]
dim(asv_table2)
#54 53

#look at sequenceing depth
min(colSums(asv_table2))
#2769 is the lowest good depth

#rarefy data 
set.seed(515)
cort_rare<-rrarefy(t(asv_table2), sample=2769)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(cort_rare  ~ 1, distance='jaccard')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#19.2
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#14.1

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color = Antibiotic, shape=Fusarium))+
  geom_point(size=3)+
  #geom_text()+
  theme_bw()+
  xlab("PC1- 19.2%")+
  ylab("PC2- 14.1%")

###Calculate alpha diversity
#load in asv table and metadata
#read in metadata and asv table
asv_table <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt", header=T)

#remove contams
contam<-which(asv_table$PCRCtrl>0)
length(contam)
#105

asv_table2<-asv_table[-contam,]
asv_table2<-asv_table2[,-54]
dim(asv_table2)
#54 53

#look at sequenceing depth
min(colSums(asv_table2))
#2769 is the lowest good depth

#rarefy data 
set.seed(515)
cort_rare<-rrarefy(t(asv_table2), sample=2769)

cort.alph<-as.data.frame(specnumber(cort_rare))
cort.alph$SampleID<-row.names(cort.alph)
cort.alph<-merge(cort.alph, meta, by='SampleID')

#plot richness
ggplot(cort.alph, aes(Antibiotic, `specnumber(cort_rare)`, fill=Fusarium))+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("sOTU Richness")

alphy<-cort.alph[which(cort.alph$Fusarium == 'Yes'),]

ggplot(alphy, aes(as.numeric(FusariumDose), `specnumber(cort_rare)`))+
  geom_point()+
  scale_y_log10()+
  stat_cor(method='spearman')+
  geom_smooth(method='lm')+
  theme_bw()

#percent inhibitory analysis
#calculate in vsearch first
vsearch -usearch_global rep_seqs/dna-sequences.fasta -db antiBd_db/AmphibBac_Inhibitory_2023.2r.fasta --strand plus --id 0.99 --blast6out synth_com_inhib_otus.txt
#Matching unique query sequences: 30 of 66 (45.45%)

#read in results from vsearch clustering against AmphiBac database
inhibitory<-read.delim("~/GitHub/Newt-cort-microbiome/Synth_com_data/synth_com_inhib_otus.txt", header=F)
#ASVs are V1 in df

#load in asv table and metadata
asv_table <- read.delim("~/GitHub/Newt-cort-microbiome/Synth_com_data/asv_table_synth_com.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Newt-cort-microbiome//Synth_com_data/cort_synth_com_map.txt", header=T)

#rarefy data 
set.seed(515)
library(vegan)
asv_table<-rrarefy(t(asv_table), sample=4591)
asv_table<-as.data.frame(t(asv_table))

#subset ASV table to only include ASVs that match database
inhib_tb<-asv_table[row.names(asv_table) %in% inhibitory$V1,]

#calculate inhibitory richness
cort.alph<-as.data.frame(specnumber(t(inhib_tb)))
cort.alph$SampleID<-row.names(cort.alph)
cort.alph<-merge(cort.alph, meta, by='SampleID')

#plot antifungal richness
library(ggplot2)
synth_com_antifungal_richness<-ggplot(cort.alph, aes(Type, `specnumber(t(inhib_tb))`, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  scale_fill_manual(values=c('blue', 'grey'))+
  theme_bw()+
  ggtitle("(D)")+
  xlab("")+
  coord_flip()+
  ylab("Antifungal richness")

#math it
t.test(cort.alph$`specnumber(t(inhib_tb))` ~ cort.alph$Type)
# = -1.7821, df = 24.724, p-value = 0.08702

#calculate colSums for total/inhibitory communities
total_sum<-as.data.frame(colSums(asv_table))
inhib_sum<-as.data.frame(colSums(inhib_tb))

#bind data together
inhib_tb2<-cbind(total_sum, inhib_sum)

#calculate percent inhibitory
inhib_tb2$per_inhib<-inhib_tb2$`colSums(inhib_tb)`/inhib_tb2$`colSums(asv_table)`

#add column for SampleID and merge metadata
inhib_tb2$SampleID<-row.names(inhib_tb2)
inhib_tb2<-merge(inhib_tb2, meta, by='SampleID')
inhib_tb2$per_inhib2<-inhib_tb2$per_inhib*100

#plot per species
library(ggplot2)
synthcom_bd_per_inhib<-
  ggplot(inhib_tb2, aes(Type, per_inhib2, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  guides(fill="none")+
  theme_bw()+
  ggtitle("(C)")+
  xlab("")+
  coord_flip()+
  scale_fill_manual(values = c('blue', 'grey'))+
  ylab("Percent Inhibitory towards Bd")

#test it with maths
t.test(inhib_tb2$per_inhib ~ inhib_tb2$Type)
#t = -1.4482, df = 27.617, p-value = 0.1588

##Make publucation figure
library(ggpubr)
ggarrange(synth_com_pcoa, synth_com_richness, synthcom_bd_per_inhib, synth_com_antifungal_richness, legend = c('none','none', 'none', 'right'))


###Calculate ASVs that respond (+/-) to cort
#read in metadata and asv table
asv_table <- read.delim("~/GitHub/Newt-cort-microbiome/asv_table_synth_com.txt", row.names=1, header=T)
meta<-read.delim("~/GitHub/Newt-cort-microbiome/cort_synth_com_map.txt", header=T)

#rarefy table
set.seed(515)
asv_table<-rrarefy(t(asv_table), sample=4591)
asv_table<-as.data.frame(t(asv_table))

#add ASV names as a column
asv_table$OTU<-row.names(asv_table)

#reshape table for calculation
library(reshape2)
asv_m<-melt(asv_table)

#add metadata
asv_m<-merge(asv_m, meta, by.x='variable', by.y='SampleID')
asv_m<-asv_m[,c(2,3,7)]

#calculate kruskal-wallis to ID differential abundant taxa with Benjamini-hochberg corrected p-value
library(ggpubr)
kruskal_results<-as.data.frame(compare_means(value ~ Type, group.by = 'OTU', p.adjust.method='fdr', method = 'kruskal.test', data=asv_m))
length(which(kruskal_results$p.adj<0.01))
#12 differentially abundant OTUs because of cort

#read in taxonomy data for each OTU & append to Kruskal-wallis results
tax<-read.delim("~/GitHub/Newt-cort-microbiome/Synth_com_data/taxonomy.tsv", header=T)
kruskal_results<-merge(kruskal_results, tax, by.x='OTU', by.y='Feature.ID', all.y=F)

#create a table of only significant OTUs
sig_krusk<-kruskal_results[which(kruskal_results$p.adj<0.01),]

library(plyr)
library(stringi)
library(tidyr)
#add and split taxonomy for Kruskal-Wallis test
sig_krusk<-separate(sig_krusk, Taxon, sep=';', , into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#calculate average abundance of each  for n-enrichment/reference sites
asv_means<-ddply(asv_m, c("OTU", 'Type'), summarize, mean=mean(value), sd=sd(value), n=length(value), se=sd/n)

#add averages to significant kruskal wallis results
sig_krusk<-merge(sig_krusk, asv_means, by.x='OTU', by.y='OTU', all.y=F)

#change to relative abundance
sig_krusk$rel_abun<-10*(sig_krusk$mean/4591)
sig_krusk$se2<-sig_krusk$se/10

#fix genus names to file and fix genus names
#write.table(sig_krusk, '~/GitHub/Newt-cort-microbiome/sig_krusk_synth_com.txt', row.names=F, sep='\t', quote=F)
sig_krusk <- read.delim("~/GitHub/Newt-cort-microbiome/Synth_com_data/sig_krusk_synth_com.txt")

#plot mean abundance
ggplot(sig_krusk, aes(OTU, rel_abun, color=Type))+
  geom_point()+
  #facet_wrap(~Phylum)+
  #scale_color_manual(values = c('#f58231', '#4363d8'))+
  theme_bw()+
  coord_flip()+
  ylab("Relative Abundance")+
  #facet_wrap(~Phylum)+
  xlab("")+
  geom_errorbar(aes(ymin=rel_abun-se2, ymax=rel_abun+se2, colour=Type), width=.2)

#heat map
ggplot(sig_krusk, aes(OTU, Type, fill=rel_abun))+
  geom_tile()+
  coord_flip()+
  facet_wrap(~Genus, scales='free_y')+
  xlab("")+
  ylab("")+
  theme_bw()+
  scale_fill_gradient2(low = "#075AFF", mid = "#FFFFCC", high = "#FF0000")
