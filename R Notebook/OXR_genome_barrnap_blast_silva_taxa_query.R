library(dplyr)

res <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/barrnap_taxonomy/OXR_barrnap_silva_blast_formatted.csv")
tax <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/taxmap_slv_ssu_ref_nr_138.1.txt", sep = "\t")

names(tax)[names(tax) == "primaryAccession"] <- "sseqid"

tax.short <- tax[tax$sseqid %in% res$sseqid,]

tax.short <- tax.short %>% select(sseqid, organism_name)

res2 <- merge(res, tax.short, by = "sseqid")

write.csv(res2, file = "C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/barrnap_taxonomy/OXR_barrnap_silva_blast_formatted-names.csv")
