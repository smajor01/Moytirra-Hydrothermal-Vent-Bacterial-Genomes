#
#
# ~~~ In this script I will try to take the gff file from the prodigal out put 
# ~~~ and take the non-CDS annotations

gff.test <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-11/trycycler/prokka_annotation/OXR-11.gff",
                       sep = "\t", skip = 10, header = F)
gff.test <- gff.test[c(1:grep("FASTA",gff.test$V1)-1),] # Remove the excess fasta file below the gff
non.cds <- gff[! gff$V3 %in% "CDS",]
non.cds$isolate <- "TEST"

# ~~~ Function to load in prokka gff and remove the CDS ~~~ #
#####
nonCDS.annot.prodigal <- function(OXR = "OXR-9"){
  library(plyr)
  library(dplyr)
  gff <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         OXR,
                         "/trycycler/prokka_annotation/",
                         OXR,
                         ".gff"),
                  sep = "\t", skip = 10, header = F)

  gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
  
  non.cds <- gff[! gff$V3 %in% "CDS",]
  non.cds$isolate <- OXR
  names(non.cds) <- c("Contig","Prg","Classification","Start","Stop","unk1","Sense","unk2","Information","Isolate")
  non.cds <- non.cds[c(10,1,2,3,4,5,6,7,8,9)]
  
  
  assign(OXR, non.cds, envir = globalenv())
  
}
#####
# ~~

# ~~~ Setting up the non-CDS table from the prokka gff
#####
OXR.table.list <- list()

files <- c("OXR-9","OXR-11","OXR-76",
           "OXR-85","OXR-96","OXR-134",
           "OXR-137","OXR-159","OXR-189",
           "OXR-199","OXR-203","OXR-209")

for (i in files){
 OXR.table.list[[i]] <- nonCDS.annot.prodigal(OXR = i)
}

OXR.nonCDS <- bind_rows(OXR.table.list)

OXR.nonCDS$Information <- stringr::str_remove(OXR.nonCDS$Information, ".*product=")

##~~~ Rename the plasmids
plasmid.names <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/OceanX Deep Water Culturing/Notes and Genome Assemblies/plasmid_names.csv", header = T)
plasmid.names$plasmid <- paste(plasmid.names$Genome, plasmid.names$CLUSTER_NAME, sep = " ")

OXR.nonCDS$Contig <- sub("_consensus_polypolish", "", OXR.nonCDS$Contig)
OXR.nonCDS$plasmid <- paste(OXR.nonCDS$Isolate,OXR.nonCDS$Contig, sep = " ")
  
new.names <- merge(OXR.nonCDS, plasmid.names, by = "plasmid")

names(new.names)

new.2 <- new.names[c(2,14,4,5,6,7,9,11)]

colnames(new.2)[colnames(new.2) == "PLASMID_NAME"] <- "Contig"

names()

# write.csv(new.2, file = "C:/Users/SamuelMajor/Box/Science/Microbe Repository/OceanX Deep Water Culturing/Manuscript/Supplementary/Table SX - nonCDS RNA.csv")
#####

# ~~~ Taking the barrnap RNA's
#####
barrnap.rrna <- function(OXR = "OXR-9"){
  library(plyr)
  library(dplyr)
  gff <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         OXR,
                         "/trycycler/",
                         OXR,
                         "_barrnap.gff"),
                  sep = "\t", skip = 1, header = F)
  
  gff$V9 <- stringr::str_remove(gff$V9, ".*product=")
  
  names(gff) <- c("Contig","Prg","Classification","Start","Stop","unk1","Sense","unk2","Information")
  gff <- gff[c(1,4,5,9)] # This is for circos. Remove for making tables with more information and uncomment other.
  # gff$Isolate <- OXR
  # gff <- gff[c("Isolate","Contig","Prg","Classification","Start","Stop","Sense","Information")]
  
  assign(OXR, gff, envir = globalenv())
  
}

barrnap.rrna("OXR-209")

OXR.table.list <- list()

files <- c("OXR-9","OXR-11","OXR-76",
           "OXR-85","OXR-96","OXR-134",
           "OXR-137","OXR-159","OXR-189",
           "OXR-199","OXR-203","OXR-209")

for (i in files){
  OXR.table.list[[i]] <- barrnap.rrna(OXR = i)
  write.table(OXR.table.list[[i]], 
            file = paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/Circos/OXR Genomes/",
                          i,
                          "/rRNA.txt"),
            col.names = FALSE, 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")
}

oxr.new <- lapply(OXR.table.list, function(df){
  if("Contig" %in% colnames(df)) {
    df %>%
    mutate(Contig = stringr::str_remove(Contig, "_consensus_polypolish")) 
  } else{
    df
    }
  })

oxr.new <- bind_rows(oxr.new)

oxr.new$plasmid <- paste(oxr.new$Isolate, oxr.new$Contig, sep = " ") 

oxr.new.2 <- merge(oxr.new, plasmid.names, by = "plasmid")

oxr.new.2 <- oxr.new.2[c(2,13,4,5,6,7,8,9)]

colnames(oxr.new.2)[colnames(oxr.new.2) == "PLASMID_NAME"] <- "Contig"

all <- rbind(new.2, oxr.new.2)

write.csv(all, file = "C:/Users/SamuelMajor/Box/Science/Microbe Repository/OceanX Deep Water Culturing/Manuscript/Supplementary/Table SX - nonCDS RNA.csv")
