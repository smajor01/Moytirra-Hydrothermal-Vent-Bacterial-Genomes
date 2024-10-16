---
title: "OXR-cds_tRNA count tables"
output: word_document
date: "2023-12-04"
---

This script is collecting the counts of each isolates CDS and tRNA

OXR-9
```{r}
# OXR-9
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-9/trycycler/prokka_annotation/OXR-9.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff

hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

clust2 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_002_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust2$Var1 <- NULL
clust2 <- as.data.frame(t(clust2))

clust3 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_003_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust3$Var1 <- NULL
clust3 <- as.data.frame(t(clust3))

clust4 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_004_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust4$Var1 <- NULL
clust4 <- as.data.frame(t(clust4))

clust5 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_005_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust5$Var1 <- NULL
clust5 <- as.data.frame(t(clust5))

clust6 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_006_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust6$Var1 <- NULL
clust6 <- as.data.frame(t(clust6))

clust11 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_011_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust11$Var1 <- NULL
clust11 <- as.data.frame(t(clust11))

clust12 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_012_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust12$Var1 <- NULL
clust12 <- as.data.frame(t(clust12))


p <- rbind.fill(clust1,clust2, clust3, clust4,clust5,clust6,clust11,clust12)
p$isolate <- rep("OXR-9", nrow(p))
p$chrom <- c("cluster_001_consensus_polypolish","cluster_002_consensus_polypolish",
             "cluster_003_consensus_polypolish","cluster_004_consensus_polypolish",
             "cluster_005_consensus_polypolish","cluster_006_consensus_polypolish",
             "cluster_011_consensus_polypolish","cluster_012_consensus_polypolish")
p <- p[,c(4,5,1,3,2)]
```


OXR-11
```{r}
# OXR-11
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-11/trycycler/prokka_annotation/OXR-11.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff

hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

oxr.11 <- rbind.fill(clust1)
oxr.11$isolate <- rep("OXR-11", nrow(q))
oxr.11$chrom <- c("cluster_001_consensus_polypolish")
oxr.11 <- oxr.11[,c(5,6,1,4,3,2)]
```

OXR-76
```{r}
# OXR-76
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-76/trycycler/prokka_annotation/OXR-76.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

oxr.76 <- rbind.fill(clust1)
oxr.76$isolate <- rep("OXR-76", nrow(q))
oxr.76$chrom <- c("cluster_001_consensus_polypolish")
oxr.76 <- oxr.76[,c(5,6,1,4,3,2)]
```

OXR-85
```{r}
# OXR-85
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-85/trycycler/prokka_annotation/OXR-85.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

clust2 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_002_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust2$Var1 <- NULL
clust2 <- as.data.frame(t(clust2))

clust3 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_003_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust3$Var1 <- NULL
clust3 <- as.data.frame(t(clust3))

oxr.85 <- rbind.fill(clust1,clust2,clust3)
oxr.85$isolate <- rep("OXR-85", nrow(oxr.85))
oxr.85$chrom <- c("cluster_001_consensus_polypolish","cluster_002_consensus_polypolish","cluster_003_consensus_polypolish")
oxr.85 <- oxr.85[,c(4,5,1,3,2)]
```

OXR-96
```{r}
# OXR-96
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-96/trycycler/prokka_annotation/OXR-96.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

oxr.96 <- rbind.fill(clust1)
oxr.96$isolate <- rep("OXR-96", nrow(oxr.96))
oxr.96$chrom <- c("cluster_001_consensus_polypolish")
oxr.96 <- oxr.96[,c(4,5,1,3,2)]
```

OXR-134
```{r}
# OXR-134
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-134/trycycler/prokka_annotation/OXR-134.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

oxr.134 <- rbind.fill(clust1)
oxr.134$isolate <- rep("OXR-134", nrow(oxr.134))
oxr.134$chrom <- c("cluster_001_consensus_polypolish")
oxr.134 <- oxr.134[,c(4,5,1,3,2)]
```

OXR-137
```{r}
# OXR-137
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-137/trycycler/prokka_annotation/OXR-137.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))

clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

clust4 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_004_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust4$Var1 <- NULL
clust4 <- as.data.frame(t(clust4))

clust9 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_009_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust9$Var1 <- NULL
clust9 <- as.data.frame(t(clust9))

oxr.137 <- rbind.fill(clust1,clust4,clust9)
oxr.137$isolate <- rep("OXR-137", nrow(oxr.137))
oxr.137$chrom <- c("cluster_001_consensus_polypolish","cluster_002_consensus_polypolish","cluster_003_consensus_polypolish")
oxr.137 <- oxr.137[,c(4,5,1,3,2)]
```
OXR-159
```{r}
# OXR-159
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-159/trycycler/prokka_annotation/OXR-159.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

clust3 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_003_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust3$Var1 <- NULL
clust3 <- as.data.frame(t(clust3))

clust4 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_004_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust4$Var1 <- NULL
clust4 <- as.data.frame(t(clust4))

clust7 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_007_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust7$Var1 <- NULL
clust7 <- as.data.frame(t(clust7))

clust8 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_008_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust8$Var1 <- NULL
clust8 <- as.data.frame(t(clust8))

clust9 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_009_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust9$Var1 <- NULL
clust9 <- as.data.frame(t(clust9))

clust10 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_010_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust10$Var1 <- NULL
clust10 <- as.data.frame(t(clust10))

clust13 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_013_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust13$Var1 <- NULL
clust13 <- as.data.frame(t(clust13))


oxr.159 <- rbind.fill(clust1,clust3, clust4, clust7,clust8,clust9,clust10,clust13)
oxr.159$isolate <- rep("OXR-159", nrow(oxr.159))
oxr.159$chrom <- c("cluster_001_consensus_polypolish","cluster_004_consensus_polypolish",
             "cluster_004_consensus_polypolish","cluster_007_consensus_polypolish",
             "cluster_008_consensus_polypolish","cluster_009_consensus_polypolish",
             "cluster_010_consensus_polypolish","cluster_013_consensus_polypolish")
oxr.159 <- oxr.159[,c(4,5,1,3,2)]
```

OXR-189
```{r}
# OXR-189
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-189/trycycler/prokka_annotation/OXR-189.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

oxr.189 <- rbind.fill(clust1)
oxr.189$isolate <- rep("OXR-189", nrow(oxr.189))
oxr.189$chrom <- c("cluster_001_consensus_polypolish")
oxr.189 <- oxr.189[,c(4,5,1,3,2)]
```

OXR-199
```{r}
# OXR-199
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-199/trycycler/prokka_annotation/OXR-199.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

clust2 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_002_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust2$Var1 <- NULL
clust2 <- as.data.frame(t(clust2))

clust3 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_003_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust3$Var1 <- NULL
clust3 <- as.data.frame(t(clust3))

clust4 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_004_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust4$Var1 <- NULL
clust4 <- as.data.frame(t(clust4))

clust5 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_005_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust5$Var1 <- NULL
clust5 <- as.data.frame(t(clust5))

clust6 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_006_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust6$Var1 <- NULL
clust6 <- as.data.frame(t(clust6))

clust7 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_007_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust7$Var1 <- NULL
clust7 <- as.data.frame(t(clust7))


oxr.199 <- rbind.fill(clust1,clust2,clust3,clust4,clust5,clust6,clust7)
oxr.199$isolate <- rep("OXR-199", nrow(oxr.199))
oxr.199$chrom <- c("cluster_001_consensus_polypolish","cluster_002_consensus_polypolish",
                   "cluster_003_consensus_polypolish","cluster_004_consensus_polypolish",
                   "cluster_005_consensus_polypolish","cluster_006_consensus_polypolish",
                   "cluster_007_consensus_polypolish")
oxr.199 <- oxr.199[,c(4,5,1,3,2)]
```

OXR-203
```{r}
# OXR-203
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-203/trycycler/prokka_annotation/OXR-203.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

oxr.203 <- rbind.fill(clust1)
oxr.203$isolate <- rep("OXR-203", nrow(oxr.203))
oxr.203$chrom <- c("cluster_001_consensus_polypolish")
oxr.203 <- oxr.203[,c(4,5,1,3,2)]
```

OXR-209
```{r}
# OXR-209
library(plyr)
library(dplyr)
gff <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-209/trycycler/prokka_annotation/OXR-209.gff",
                sep = "\t", skip = 10, header = F)

gff <- gff[c(1:grep("FASTA",gff$V1)-1),] # Remove the excess fasta file below the gff
# Look at the hypothetical proteins
hyp <- gff[grep("product=hypothetical", gff$V9, ignore.case = T),]
table(hyp$V1)

# Show the chromosomes (contigs)
unique(sort(gff$V1))
# Make a table of unique strings in column 3 which shows "CDS","tRNA","tmRNA" (in this sample)
clust1 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_001_consensus_polypolish"]))  %>% data.frame(., row.names = .$Var1)
# Just table manipulation to get it in to the form I'd like. 
clust1$Var1 <- NULL
clust1 <- as.data.frame(t(clust1))

clust3 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_003_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust3$Var1 <- NULL
clust3 <- as.data.frame(t(clust3))

clust4 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_004_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust4$Var1 <- NULL
clust4 <- as.data.frame(t(clust4))

clust5 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_005_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust5$Var1 <- NULL
clust5 <- as.data.frame(t(clust5))

clust6 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_006_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust6$Var1 <- NULL
clust6 <- as.data.frame(t(clust6))

clust7 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_007_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust7$Var1 <- NULL
clust7 <- as.data.frame(t(clust7))

clust10 <- as.data.frame(table(gff$V3[gff$V1 %in% "cluster_010_consensus_polypolish"])) %>% data.frame(., row.names = .$Var1)
clust10$Var1 <- NULL
clust10 <- as.data.frame(t(clust10))

oxr.209 <- rbind.fill(clust1,clust3,clust4,clust5,clust6,clust7,clust10)
oxr.209$isolate <- rep("OXR-209", nrow(oxr.209))
oxr.209$chrom <- c("cluster_001_consensus_polypolish","cluster_003_consensus_polypolish",
                   "cluster_004_consensus_polypolish","cluster_005_consensus_polypolish",
                   "cluster_006_consensus_polypolish","cluster_007_consensus_polypolish",
                   "cluster_010_consensus_polypolish")
oxr.209 <- oxr.209[,c(4,5,1,3,2)]
```
Combine all the tables
```{r}
full <- rbind.fill(p,oxr.11,oxr.76,oxr.85,oxr.96,oxr.134,oxr.137,oxr.159,oxr.189,oxr.199,oxr.203,oxr.209)

output.file <- file("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-cds_tRNA.txt", "wb")
write.table(full, file = output.file, row.names = F, col.names = F)
close(output.file)
```