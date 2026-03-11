#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(phyloseq)
library(decontam)

#read in metadata and asv table
asv_table <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/xenopus_asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt", header=T)

#remove contams
# Convert to phyloseq objects, makes things easier for decontam package
OTU <- otu_table(as.matrix(asv_table), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
sample_names(SAM) <- meta$SampleID  # ensure they match
sample_names(OTU) <- colnames(asv_table)
ps <- phyloseq(OTU, SAM)

# Identify negative controls, in a column called 'Control', who would've thought?
sample_data(ps)$is.neg <- sample_data(ps)$Control %in% c("Yes")

# Run decontam based on prevalance method
#check the other methods in the help file, but the prevalance method is most commonly used
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")

# Output contaminants
contaminants <- as.data.frame(rownames(contamdf.prev[contamdf.prev$contaminant, ]))
names(contaminants)<-'OTUS'
write.table(contaminants, "~/Documents/GitHub/Gnotobiotic_Xenopus/contaminants.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#how many? 
dim(contaminants)
#[1] 87  1

#remove contaminants
dim(asv_table)
#[1] 169  54

otu_no_contam<-asv_table[-which(row.names(asv_table) %in% contaminants$OTUS), ]
dim(otu_no_contam)
#[1] 82 54

#the dim isn't necessary, I'm just paranoid about not removing what I say I want removed

#write to file
write.table(otu_no_contam, '~/Documents/GitHub/Gnotobiotic_Xenopus/no_contam_asv_table.txt', row.names = T, sep='\t', quote=F)

#remove the PCR control from OTU table
otu_no_contam<-otu_no_contam[,-54]

#look at sequenceing depth
min(colSums(otu_no_contam))
#2832 is the lowest good depth

#rarefy data 
set.seed(515)
cort_rare<-rrarefy(t(otu_no_contam), sample=2832)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(cort_rare  ~ 1, distance='bray')
#ko_pcoa<-capscale(cort_rare  ~ 1, distance='jaccard')

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
#24.9
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#17.2

#make color pallet
cb_palette <- c(
  "#332288", # dark blue
  "#88CCEE", # light blue
  "#44AA99", # teal
  "#117733", # green
  "#DDCC77", # sand/yellow
  "#CC6677", # muted red
  "#AA4499", # purple
  "#882255"  # dark rose
)

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, color = Antibiotic, shape=Fusarium))+
  geom_point(size=3)+
  #geom_text()+
  theme_bw()+
  scale_color_manual(values=cb_palette)+
  xlab("PC1- 24.9%")+
  ylab("PC2- 17.2%")

###Calculate alpha diversity
#load in asv table and metadata
#read in metadata and asv table
asv_table <- read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/no_contam_asv_table.txt", row.names=1, header=T)
meta<-read.delim("~/Documents/GitHub/Gnotobiotic_Xenopus/molly_xenopus_map.txt", header=T)

#remove PCR control
asv_table <- asv_table[,-54]

#look at sequenceing depth
min(colSums(asv_table))
#2832 is the lowest good depth

#rarefy data 
set.seed(515)
cort_rare<-rrarefy(t(asv_table), sample=2832)

cort.alph<-as.data.frame(specnumber(cort_rare))
cort.alph$SampleID<-row.names(cort.alph)
cort.alph<-merge(cort.alph, meta, by='SampleID')

#plot richness
cort.alph2<-cort.alph[which(cort.alph$Fusarium == 'No'),]

ggplot(cort.alph2, aes(Antibiotic, `specnumber(cort_rare)`))+
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
  theme_bw()+
  ylab('sOTUS Richness')+
  xlab('Fusarium Dose (spores/mL')

ggplot(alphy, aes(FusariumDose, `specnumber(cort_rare)`))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ylab('sOTUS Richness')+
  xlab('Fusarium Dose (spores/mL')

#run homogeneity of variances
bartlett.test(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Antibiotic)
#Bartlett's K-squared = 8.8968, df = 7, p-value = 0.2602

bartlett.test(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Fusarium)
#Bartlett's K-squared = 0.54816, df = 1, p-value = 0.4591

#run AOV
summary(aov(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Antibiotic))
#                     Df  Sum Sq Mean Sq F value   Pr(>F)    
#cort.alph$Antibiotic  7   1748  249.73    5.63    0.000104 ***
#Residuals            45   1996   44.35

TukeyHSD(aov(cort.alph$`specnumber(cort_rare)` ~ cort.alph$Antibiotic))
#sig hits: Trimethoprim-AllAB, Trimethoprim-Ampicillin, Trimethoprim-FuskNoAB, Control-AllAB 
