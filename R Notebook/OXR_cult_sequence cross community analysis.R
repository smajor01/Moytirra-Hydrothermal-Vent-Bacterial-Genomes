
# This script will separate the abundant and the Representative sequences. Re-name them, and then combine them and make a fasta file
library(dada2)
library(phyloseq)
library(dplyr)
library(tidyr)
library(reshape2)

oxr.table <- read.csv("c:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-221129_bact_genus_pure.csv", 
                      row.names = "CULT_NAMES")
oxr.table$X <- NULL
oxr.table <- oxr.table[! row.names(oxr.table) %in% c("SP12-124","SP12-125B", "SP12-126"),]

##### Representative ASVs
rep.asv <- oxr.table[c("TAXA_REP_SEQ","TAXA_ASV")]
names(rep.asv) <- c("sequence","ASV")
# rep.asv$TAXA_ASV <- gsub("ASV","REP-", rep.asv$TAXA_ASV) # rename the asv's

##### Abundant ASV
abun.asv <- oxr.table[c("ABUN_ASV_SEQ","ABUN_ASV")]
names(abun.asv) <- c("sequence","ASV")
# abun.asv$ABUN_ASV <- gsub("ASV","ABUN-", abun.asv$ABUN_ASV)# rename the asv's

seqs <- rbind(rep.asv, abun.asv)
seqs <- seqs[! duplicated(seqs$sequence),]
row.names(seqs) <- seqs$ASV
seqs$ASV <- NULL 

library(seqinr)
write.fasta(as.list(seqs$sequence),rownames(seqs), "OXR_cult.fa")

# in Terminal, run the following code

# makeblastdb -in "OXR_ASV_search\OXR-16S_metabarcodes\OXR-16S_metabarcodes\OXR-16S_zotus.fa" -dbtype nucl -out OXR_ASV_search\OXR-16S_metabarcodes\OXR_metabarcode_db -title "OXR Database"

# blastn -query OXR_ASV_search\OXR_cult.fa -db OXR_ASV_search\OXR-16S_metabarcodes\OXR_metabarcode_db -outfmt 6 -out OXR_ASV_search\OXR_blastn_search.txt -max_target_seqs 5


# I will take the hit ID's from the blast search and find the % representation in the Ocean X metabarcodes
# read in the hits
hits <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR_blastn_search-5.txt", sep = "\t", header = F)
# new
hits <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/barrnap_trycycler/OXR_16S_trycycler_metabarcode_blastn-5.txt", sep = "\t", header = F)
names(hits) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
hits$qseqid <- gsub("_consensus_polypolish","",hits$qseqid)

# Extract the Hit IDs
sseq <- hits$sseqid

# read in the metabarcodes hit table
meta.tab <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR-16S_metabarcodes/OXR-16S_zotutab.txt", sep = "\t", row.names = 1)

# Transform values into percentatge and investigate the by Sample coverage
meta.tab.sel <- meta.tab[rownames(meta.tab) %in% sseq,]

meta.perc <- round(prop.table(as.matrix(meta.tab), margin = 2),4)*100

meta.sel <- meta.perc[rownames(meta.perc) %in% sseq,]

perc.5.tot <- sort(colSums(meta.sel)) # the total percentage the top 5 hits in the metabarcoding data yield.
  
# Transform values into percentage and investigate the total community coverage
meta.perc.tots <- round(prop.table(as.matrix(meta.tab)),6)*100

meta.tots.sel <- meta.perc.tots[rownames(meta.perc.tots) %in% sseq,]

head(sort(rowSums(meta.perc.tots)))

# Next, Take the top single hits of each ASV, 
# Find and download the OXR-16S_zotus.nr_v138_1.wang.taxonomy
# Manually replace the "tabs" and the ";" with commas. Add headers "ASV","Kingdom","Phylum","Class","Order","Family","Genus" 
# Replace "zotu" with "ASV" and remove the appended parenthetic numbers "(123)"

# extract the ASV's from the Identification table 
rep.gen <- oxr.table[c("PHYLUM","CLASS","ORDER","FAMILY","GENUS","TAXA_ASV")]
names(rep.gen) <- c("cultured_Phylum","cultured_Class","cultured_Order","cultured_Family","cultured_Genus","ASV")

##### Abundant ASV
abun.gen <- oxr.table[c("PHYLUM","CLASS","ORDER","FAMILY","GENUS","ABUN_ASV")]
names(abun.gen) <- c("cultured_Phylum","cultured_Class","cultured_Order","cultured_Family","cultured_Genus","ASV")
# abun.asv$ABUN_ASV <- gsub("ASV","ABUN-", abun.asv$ABUN_ASV)# rename the asv's

gen <- rbind(rep.gen, abun.gen)
gen <- gen[! duplicated(gen$ASV),]
# names(gen) <- c("Cultured Genus","Cultured ASV")
row.names(gen) <- gen$ASV
# gen$ASV <- NULL 

# read in the single hit file
singles <- read.csv("c:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR_blastn_search.txt"
                    , sep = "\t", header = F)
names(singles) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
row.names(singles) <- singles$qseqid

# Merge the single hits with the ASV and genus
single.cult <- merge(singles, gen, by = "row.names")

# read in the taxonomic information of the metabarcoding data
meta.tax <- read.csv("c:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR-16S_metabarcodes/OXR-16S_zotus.nr_v138_1.csv")
names(meta.tax) <- c("sseqid","meta_Kingdom","meta_Phylum","meta_Class","meta_Order","meta_Family","meta_Genus")

# merge the single.cult data frame with the metabarcoding taxonomy
single.cult.2 <- merge(single.cult, meta.tax, by = "sseqid")

# Add color
single.cult.2$color <- ifelse(single.cult.2$qseqid %in% c("ASV50","ASV18","ASV71",
                                                          "ASV32","ASV5","ASV38",
                                                          "ASV46","ASV41","ASV40",
                                                          "ASV81","ASV12"),
                              "blue",
                              ifelse(single.cult.2$qseqid %in% c("ASV45","ASV64"),
                                     "darkblue","darkgray"))


# jpeg(filename = "c:/Users/smajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR ASV hit identities.jpeg"
#      , width = 1500
#      , height = 800)
png(filename = "c:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR ASV hit identities.png",
    width = 1000,
    height = 1000/(4/3))

ggplot(single.cult.2, aes(x = reorder(qseqid,pident), y= pident, fill = color)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90
                                   , hjust = 1
                                   , vjust = 0.4
                                   , size = 15)
        , axis.title.x = element_text(size = 20)
        , axis.title.y = element_text(size = 20)
        , plot.title = element_text(size = 25)) +
  # scale_x_discrete(labels = single.cult.2$cultured_Genus)+
  # single.cult.2$meta_Genus) +
  xlab("Cultured Genus") +
  ylab("% Identity") +
  ggtitle("Cultured Bacterial 16S rRNA gene Identity in Moytirra Metabarcoding Data") +
  geom_hline(yintercept = c(90,97,100), linetype = c("solid","dashed","solid"), color = c("red","blue","red")) +
  # scale_y_continuous(limits = c(75,100), breaks = seq(80,100, by =10))
  coord_cartesian(ylim = c(70,100)) +
  geom_text(aes(label= cultured_Genus)
  , color = "white"
  , size = 5
  , position = position_stack(vjust = 0.75)
  , angle = 90
  , hjust = "left")+
  scale_fill_identity()

dev.off()

# jpeg(filename = "c:/Users/smajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR ASV hit ASV.jpeg"
#      , width = 1500
#      , height = 800)


png(filename = "c:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-16S_metabarcode/OXR_ASV_search/OXR ASV hit ASV.png"
     , width = 1000
     , height = 1000/(4/3))

# reorder(qseqid, pident)

ggplot(single.cult.2, aes(x = reorder(qseqid,pident) ,y= pident, fill = color)) +
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90
                                   , hjust = 1
                                   , vjust = 0.4
                                   , size = 15)
        , axis.title.x = element_text(size = 20)
        , axis.title.y = element_text(size = 20)
        , plot.title = element_text(size = 25)) +
  # scale_x_discrete(labels = single.cult.2$qseqid) +
  xlab("Cultured ASV") +
  ylab("% Identity") +
  ggtitle("Cultured Bacterial 16S rRNA gene Identity in Moytirra Metabarcoding Data") +
  geom_hline(yintercept = c(90,97,100), linetype = c("solid","dashed","solid"), color = c("red","blue","red")) +
  # scale_y_continuous(limits = c(75,100), breaks = seq(80,100, by =10))
  coord_cartesian(ylim = c(70,100)) +
  geom_text(aes(label=sseqid)
            , color = "white"
            , size = 5
            , position = position_stack(vjust = 0.75)
            , angle = 90
            , hjust = "left")+
  scale_fill_identity()
dev.off()

