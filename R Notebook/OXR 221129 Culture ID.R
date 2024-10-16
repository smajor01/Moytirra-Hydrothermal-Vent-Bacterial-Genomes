#
#
BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("msa")
install.packages("phangorn")

library(dada2)

source("C:/Users/smajor/Box/Science/SOPs and Protocols/Microbe/CODES for MiSeq Data Analysis/DADA-2 Pipeline v2.R")

setwd("C:/Users/smajor/Box/Science/Microbe Repository/16s-data/OXR cultures/221129_OXR/")

path <- ("C:/Users/smajor/Box/Science/Microbe Repository/16s-data/OXR cultures/221129_OXR/221129_OXR sequences/")

db <- "C:/Databases/silva_nr99_v138.1_train_set.fa.gz"

DBspec <- "C:/Databases/silva_species_assignment_v138.1.fa.gz"

DADA2qualplots(dataPath = path, PlotsofPairs = 3)

OXR.221129.SeqTab <- DADA2QfiltandTaxa(dataPath = path,
                                  FwrdStrtTrim = 10, RevStrtTrim = 10,
                                  FwrdEndTrim = 200, RevEndTrim = 230,
                                  forEE = 2,
                                  revEE = 2,
                                  PickTaxa = T,
                                  database = db)

OXR.221129.tree <- phylotree_UPGMA(OXR.221129.SeqTab)

OXR.221129.PhySeq <- buildPhyloseq(OXR.221129.SeqTab, taxa, OXR.221129.tree)

source("C:/Users/smajor/Box/Science/SOPs and Protocols/Microbe/CODES for MiSeq Data Analysis/writePhySeqs.R")

writePhySeqs(OXR.221129.PhySeq, path, Sample.Name = "OXR-221129")

source("C:/Users/smajor/Box/Experiments/R functions/SampleID-6.R")

OXR.221129.ID <- cultureID.2(OXR.221129.PhySeq, "Genus")

write.csv(file = "C:/Users/smajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-221129_bact_genus.csv", OXR.221129.ID)

pure.OXR.221129 <- OXR.221129.ID[! OXR.221129.ID$CULTURE_QUALITY %in% c("Mixed","Unknown"),]

write.csv(file = "C:/Users/smajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-221129_bact_genus_pure.csv", pure.OXR.221129)
