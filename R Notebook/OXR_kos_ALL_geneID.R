

kos.geneID <- function(genome){
  packages <- c("dplyr","data.table")
  sapply(c(packages), require, character.only = TRUE)
annos <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         genome,
                         "/trycycler/eggNOG_annotations/",
                         genome,
                         "_ALL_annos_w-DESeq.txt"), 
                  sep = "\t",  header = T)

#####
# KEGG_ko

ko <- annos[,c(1,2,9:110,121)]

ko[is.na(ko)] <- "-"

ko$KEGG_ko <- gsub("ko:","",ko$KEGG_ko)

test <- splitstackshape::cSplit(ko, "KEGG_ko", sep=",")

test2 <- melt(setDT(test), id.vars = c(1:104), variable.name = "KEGG_ko")

test2$KEGG_ko <- NULL

test3 <- na.omit(test2)

output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                           genome,
                           "/trycycler/eggNOG_annotations/",
                           genome,
                           "_ALL_ko.txt"), "wb")
write.table(test3, file = output.file, row.names = F, col.names = T, sep = "\t")
close(output.file)

assign("ko",test3, envir = globalenv())
#####
# KEGG_Pathway
ko <- annos[,c(1,2,9:110,122)]

ko[is.na(ko)] <- "-"

test <- splitstackshape::cSplit(ko, "KEGG_Pathway", sep=",")

test2 <- melt(setDT(test), id.vars = c(1:104),variable.name = "KEGG_Pathway")
# c("eggNOG_name","eggNOG_alt_name","TPM.C1B11F1","TPM.C1B11F2","TPM.C1B11F3","COUNT.C1B11F1","COUNT.C1B11F2","COUNT.C1B11F3"), variable.name = "KEGG_Pathway")

test2$KEGG_Pathway <- NULL

test3 <- na.omit(test2)

output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                           genome,
                           "/trycycler/eggNOG_annotations/",
                           genome,
                           "_ALL_pathway_DE.txt"), "wb")
write.table(test3, file = output.file, row.names = F, col.names = T, sep = "\t")
close(output.file)

assign("pathway",test3, envir = globalenv())
}
