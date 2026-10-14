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
