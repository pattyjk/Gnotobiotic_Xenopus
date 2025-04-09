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

#run homogeneity of variances
bartlett.test(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Antibiotic)
#Bartlett's K-squared = 8.8968, df = 7, p-value = 0.2602

bartlett.test(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Fusarium)
#Bartlett's K-squared = 0.54816, df = 1, p-value = 0.4591

#run AOV
summary(aov(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Antibiotic))
TukeyHSD(aov(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Antibiotic))


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
