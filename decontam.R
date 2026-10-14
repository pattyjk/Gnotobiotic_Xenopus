# ---- paths: edit `dir` only -------------------------------------------------
dir <- "~/Documents/GitHub/Gnotobiotic_Xenopus/"

#load packages (only phyloseq + decontam are used in this script)
library(phyloseq)
library(decontam)

#read in metadata and asv table
asv_table <- read.delim(file.path(dir, "xenopus_asv_table.txt"), row.names=1, header=T)
meta<-read.delim(file.path(dir, "molly_xenopus_map.txt"), header=T)

#remove contams
# Convert to phyloseq objects, makes things easier for decontam package
OTU <- otu_table(as.matrix(asv_table), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
stopifnot(setequal(meta$SampleID, colnames(asv_table)))   # same samples in both files
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
write.table(contaminants, file.path(dir, "contaminants.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#how many? 
dim(contaminants)
#[1] 87  1

#remove contaminants
dim(asv_table)
#[1] 169  54

# logical index (a negative which() would drop EVERY row if no contaminant were found)
otu_no_contam<-asv_table[!(row.names(asv_table) %in% contaminants$OTUS), ]
dim(otu_no_contam)
#[1] 82 54

#the dim isn't necessary, I'm just paranoid about not removing what I say I want removed

#write to file
write.table(otu_no_contam, file.path(dir, "no_contam_asv_table.txt"), row.names = T, sep='\t', quote=F)
