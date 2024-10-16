
kegg.pathway.queries <- function(genome){
  packages <- c("dplyr")
  sapply(c(packages), require, character.only = TRUE)

# Read in the 
path.to.ko <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                              genome,"/trycycler/eggNOG_annotations/kegg_query_all_transcripts_ALL/ALL_Finaloutput.txt"),
                       header=F, sep = "")

names(path.to.ko) <- c("pathway","ko")

# read in pathway_name.txt
path.name <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                             genome,
                             "/trycycler/eggNOG_annotations/kegg_query_all_transcripts_ALL/pathway_name.txt"),
                      header=F, sep = "\t")

names(path.name) <- c("pathway","name")

# read in pathway_name.txt
path.class <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                              genome,
                              "/trycycler/eggNOG_annotations/kegg_query_all_transcripts_ALL/pathway-to-class.txt"), header=F, sep = "\t")

names(path.class) <- c("pathway","class","description")
# something was funky, so I fixed it
path.class[path.class == "Folding"] <- "Folding sorting and degradation"
path.class <- path.class[- grep("sorting and degradation", path.class$pathway),]

# Merge the pathway information with the KO's
pathway.all <- merge(path.to.ko, path.name, by = "pathway", all = T)
pathway.all <- merge(pathway.all, path.class, by = "pathway", all = T)

# read in the KO's for the organism
trans.ko <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                            genome,
                            "/trycycler/eggNOG_annotations/",
                            genome,
                            "_ALL_ko.txt"), header = T, sep = "\t")

names(trans.ko)[names(trans.ko) == "value"] <- "ko"

# Merge the pathway information with the ko list
all <- merge(trans.ko, pathway.all, by = "ko", all = T)

# read in the annos with the DESeq info
deseq <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         genome,
                         "/trycycler/eggNOG_annotations/",
                         genome,
                         "_ALL_annos_w-DESeq.txt"),
                  header = T, sep = "\t")

deseq <- deseq[,c(1:8,117,118,130)]

# merge with all the pathway information
WOW <- merge(deseq, all, by = "eggNOG_name", all = T)
# remove NA's in the eggNOG_name field
WOW <- WOW[! is.na(WOW$eggNOG_name),]

######
regions <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                           genome,
                           "/trycycler/eggNOG_annotations/cds.regions.txt")
                    , sep = " ", header = F)

test <- splitstackshape::cSplit(regions, 1, sep="\t") %>% .[,c("V1_1","V1_2","V2","V3")]
names(test) <- c("eggNOG_name","start","stop","sense")

WOW <- merge(test, WOW, by = "eggNOG_name", all = T)
###


output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                           genome,
                           "/trycycler/eggNOG_annotations/",
                           genome,
                           "_ALL_all_data.txt"), "wb")
write.table(WOW, file = output.file, row.names = F, col.names = T, sep = "\t")
close(output.file)
}
