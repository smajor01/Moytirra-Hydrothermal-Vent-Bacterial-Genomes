library(dplyr)

bact.pure <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/16s-data/OXR cultures/OXR-221129_bact_genus_pure.csv", row.names = 2)

meta <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/Culture Tracking Records/OXR/OXR Cultures-3.csv", row.names=2)

meta <- meta[row.names(meta) %in% row.names(bact.pure),]
bact.pure <- bact.pure[row.names(bact.pure) %in% row.names(meta),]

all <- cbind(bact.pure, meta)

(table(bact.pure$CLASS) / 176)*100

sort((table(bact.pure$GENUS)/ 176)*100)

unique(bact.pure$PHYLUM)

rhodo.5 <- all[all$INC_TEMP %in% "5" & all$ORDER %in% "Rhodobacterales",]
rhodo.5.agar <- all[all$INC_TEMP %in% "5" & all$ORDER %in% "Rhodobacterales" & all$SOLIDIFYING_AGENT %in% "Agar" ,]
rhodo.5.gum <- all[all$INC_TEMP %in% "5" & all$ORDER %in% "Rhodobacterales" & all$SOLIDIFYING_AGENT %in% "Gellan Gum" ,]
rhodo.20 <- all[all$INC_TEMP %in% "20" & all$ORDER %in% "Rhodobacterales",]
rhodo.20.agar <- all[all$INC_TEMP %in% "20" & all$ORDER %in% "Rhodobacterales" & all$SOLIDIFYING_AGENT %in% "Agar" ,]
rhodo.20.gum <- all[all$INC_TEMP %in% "20" & all$ORDER %in% "Rhodobacterales" & all$SOLIDIFYING_AGENT %in% "Gellan Gum" ,]

table(rhodo.5.agar$GENUS)
table(rhodo.5.gum$GENUS)
table(rhodo.20.agar$GENUS)
table(rhodo.20.gum$GENUS)

actino <- all[all$CLASS %in% "Actinobacteria",]

actino.5 <- all[all$INC_TEMP %in% "5" & all$CLASS %in% "Actinobacteria",]
actino.5.agar <- all[all$INC_TEMP %in% "5" & all$CLASS %in% "Actinobacteria" & all$SOLIDIFYING_AGENT %in% "Agar" ,]
actino.5.gum <- all[all$INC_TEMP %in% "5" & all$CLASS %in% "Actinobacteria" & all$SOLIDIFYING_AGENT %in% "Gellan Gum" ,]
actino.20 <- all[all$INC_TEMP %in% "20" & all$CLASS %in% "Actinobacteria",]
actino.20.agar <- all[all$INC_TEMP %in% "20" & all$CLASS %in% "Actinobacteria" & all$SOLIDIFYING_AGENT %in% "Agar" ,]
actino.20.gum <- all[all$INC_TEMP %in% "20" & all$CLASS %in% "Actinobacteria" & all$SOLIDIFYING_AGENT %in% "Gellan Gum" ,]

table(actino.5.agar$GENUS)
table(actino.5.gum$GENUS)
table(actino.20.agar$GENUS)
table(actino.20.gum$GENUS)

table(rhodo.20$GENUS)
str(all)

oxr.genomes <- all[row.names(all) %in% c("OXR-9","OXR-11","OXR-76",
                                         "OXR-85","OXR-96","OXR-134",
                                         "OXR-137","OXR-159","OXR-189",
                                         "OXR-199", "OXR-203","OXR-209"),]
