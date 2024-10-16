total.diff.analysis <- function(genome){
  packages <- c("dplyr","DESeq2")
  sapply(c(packages), require, character.only = TRUE)
  dta <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         genome,
                         "/",
                         genome,
                         "_expt.ct.txt"), sep = "\t",
                  header=TRUE, row.names = "gene_id")
  ###-------- define data structure --------###
  samples <- data.frame(row.names=c("C1B10F1","C1B10F2","C1B10F3",
                                    "C1B11F1","C1B11F2","C1B11F3",
                                    "C1B12F1","C1B12F2","C1B12F3",
                                    "C1B1F1","C1B1F2","C1B1F3",
                                    "C1B2F1","C1B2F2","C1B2F3",
                                    "C1B3F1","C1B3F2","C1B3F3",
                                    "C1B4F1","C1B4F2","C1B4F3",
                                    "C1B5F1","C1B5F2","C1B5F3",
                                    "C1B6F1","C1B6F2","C1B6F3",
                                    "C1B7F1","C1B7F2","C1B7F3",
                                    "C1B8F1","C1B8F2","C1B8F3",
                                    "C1B9F1","C1B9F2","C1B9F3",
                                    "C2B1F1","C2B1F2","C2B1F3",
                                    "C2B2F1","C2B2F2","C2B2F3",
                                    "C2B3F1","C2B3F2","C2B3F3",
                                    "C2B4F1","C2B4F2","C2B4F3",
                                    "C2B5F1","C2B5F2","C2B5F3"), 
                        condition=as.factor(c(rep("In",45),
                                              # rep("C1B12",3),rep("C2B1",3))))
                                              rep("Out",6))))
  ### Get all samples in the same order
  all(rownames(samples) %in% colnames(dta))
  all(rownames(samples) == colnames(dta))
  dta <- dta[, rownames(samples)]
  all(rownames(samples) == colnames(dta))
  ### Make data integers
  dta2 <- as.data.frame(lapply(dta, as.integer))
  rownames(dta2) <- rownames(dta)
  colSums(dta)
  ###-------- create deseq object and process --------###
  dta2["fakey",] <- rep(1, ncol(dta2))
  
  cds <- DESeqDataSetFromMatrix(countData = dta2, colData=samples, design=~condition)
  cds1 <- DESeq(cds)
  
  ref <- results(cds1)
  ref[order(ref$padj),]
  
  c1b11vOut <- results(cds1, contrast=c("condition","In","Out"))
  data <- as.data.frame(c1b11vOut[order(c1b11vOut$padj),])
  # View(as.data.frame(c1b11vOut[order(c1b11vOut$padj),]))
  
  output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                             genome,
                             "/trycycler/eggNOG_annotations/",
                             genome,
                             "_ALL_DESeq.txt"), "wb")
  write.table(data, file = output.file, row.names = T, col.names = T)
  close(output.file)
  
  
  ids <- row.names(data) %>% as.factor()
  ids <- as.factor(ids)
  
  
  output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                             genome,
                             "/trycycler/eggNOG_annotations/",
                             genome,
                             "_ALL_trans.txt"), "wb")
  write.table(ids, file = output.file, row.names = F, col.names = F)
  close(output.file)
  
  #####
  # Add TPM data to the "data" data frame
  tpm <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         genome,
                         "/",
                         genome,
                         "_TPM.txt"), 
                  sep = "\t",
                  header=TRUE, row.names = "gene_id")
  tpm2 <- tpm[,c("C1B10F1","C1B10F2","C1B10F3",
                 "C1B11F1","C1B11F2","C1B11F3",
                 "C1B12F1","C1B12F2","C1B12F3",
                 "C1B1F1","C1B1F2","C1B1F3",
                 "C1B2F1","C1B2F2","C1B2F3",
                 "C1B3F1","C1B3F2","C1B3F3",
                 "C1B4F1","C1B4F2","C1B4F3",
                 "C1B5F1","C1B5F2","C1B5F3",
                 "C1B6F1","C1B6F2","C1B6F3",
                 "C1B7F1","C1B7F2","C1B7F3",
                 "C1B8F1","C1B8F2","C1B8F3",
                 "C1B9F1","C1B9F2","C1B9F3",
                 "C2B1F1","C2B1F2","C2B1F3",
                 "C2B2F1","C2B2F2","C2B2F3",
                 "C2B3F1","C2B3F2","C2B3F3",
                 "C2B4F1","C2B4F2","C2B4F3",
                 "C2B5F1","C2B5F2","C2B5F3")] 
  tpm3 <- tpm2[row.names(tpm2) %in% row.names(data),]
  names(tpm3) <- c("TPM-C1B10F1","TPM-C1B10F2","TPM-C1B10F3",
                   "TPM-C1B11F1","TPM-C1B11F2","TPM-C1B11F3",
                   "TPM-C1B12F1","TPM-C1B12F2","TPM-C1B12F3",
                   "TPM-C1B1F1","TPM-C1B1F2","TPM-C1B1F3",
                   "TPM-C1B2F1","TPM-C1B2F2","TPM-C1B2F3",
                   "TPM-C1B3F1","TPM-C1B3F2","TPM-C1B3F3",
                   "TPM-C1B4F1","TPM-C1B4F2","TPM-C1B4F3",
                   "TPM-C1B5F1","TPM-C1B5F2","TPM-C1B5F3",
                   "TPM-C1B6F1","TPM-C1B6F2","TPM-C1B6F3",
                   "TPM-C1B7F1","TPM-C1B7F2","TPM-C1B7F3",
                   "TPM-C1B8F1","TPM-C1B8F2","TPM-C1B8F3",
                   "TPM-C1B9F1","TPM-C1B9F2","TPM-C1B9F3",
                   "TPM-C2B1F1","TPM-C2B1F2","TPM-C2B1F3",
                   "TPM-C2B2F1","TPM-C2B2F2","TPM-C2B2F3",
                   "TPM-C2B3F1","TPM-C2B3F2","TPM-C2B3F3",
                   "TPM-C2B4F1","TPM-C2B4F2","TPM-C2B4F3",
                   "TPM-C2B5F1","TPM-C2B5F2","TPM-C2B5F3")
  
  # Add the count Data 
  dta2 <- dta[,c("C1B10F1","C1B10F2","C1B10F3",
                 "C1B11F1","C1B11F2","C1B11F3",
                 "C1B12F1","C1B12F2","C1B12F3",
                 "C1B1F1","C1B1F2","C1B1F3",
                 "C1B2F1","C1B2F2","C1B2F3",
                 "C1B3F1","C1B3F2","C1B3F3",
                 "C1B4F1","C1B4F2","C1B4F3",
                 "C1B5F1","C1B5F2","C1B5F3",
                 "C1B6F1","C1B6F2","C1B6F3",
                 "C1B7F1","C1B7F2","C1B7F3",
                 "C1B8F1","C1B8F2","C1B8F3",
                 "C1B9F1","C1B9F2","C1B9F3",
                 "C2B1F1","C2B1F2","C2B1F3",
                 "C2B2F1","C2B2F2","C2B2F3",
                 "C2B3F1","C2B3F2","C2B3F3",
                 "C2B4F1","C2B4F2","C2B4F3",
                 "C2B5F1","C2B5F2","C2B5F3")] 
  dta3 <- dta2[row.names(dta2) %in% row.names(data),]
  names(dta3) <- c("COUNT-1B10F1","COUNT-1B10F2","COUNT-1B10F3",
                   "COUNT-1B11F1","COUNT-1B11F2","COUNT-1B11F3",
                   "COUNT-1B12F1","COUNT-1B12F2","COUNT-1B12F3",
                   "COUNT-1B1F1","COUNT-1B1F2","COUNT-1B1F3",
                   "COUNT-1B2F1","COUNT-1B2F2","COUNT-1B2F3",
                   "COUNT-1B3F1","COUNT-1B3F2","COUNT-1B3F3",
                   "COUNT-1B4F1","COUNT-1B4F2","COUNT-1B4F3",
                   "COUNT-1B5F1","COUNT-1B5F2","COUNT-1B5F3",
                   "COUNT-1B6F1","COUNT-1B6F2","COUNT-1B6F3",
                   "COUNT-1B7F1","COUNT-1B7F2","COUNT-1B7F3",
                   "COUNT-1B8F1","COUNT-1B8F2","COUNT-1B8F3",
                   "COUNT-1B9F1","COUNT-1B9F2","COUNT-1B9F3",
                   "COUNT-2B1F1","COUNT-2B1F2","COUNT-2B1F3",
                   "COUNT-2B2F1","COUNT-2B2F2","COUNT-2B2F3",
                   "COUNT-2B3F1","COUNT-2B3F2","COUNT-2B3F3",
                   "COUNT-2B4F1","COUNT-2B4F2","COUNT-2B4F3",
                   "COUNT-2B5F1","COUNT-2B5F2","COUNT-2B5F3")
  # Remove the "Fakey" data
  data <- data[! (row.names(data) %in% "fakey"),]
  
  assign("data", data, envir = globalenv())
  assign("tpm3", tpm3, envir = globalenv())
  assign("dta3", dta3, envir = globalenv())
}

# new <- cbind(data[sort(row.names(data)),],tpm3[sort(row.names(tpm3)),],dta3[sort(row.names(dta3)),])

# Can't do this until the "OXR-X_trans.txt" is run through the bash script "DESeq_diff_expr_from_eggNOG-2/3.sh" which collects the expressed genes with the eggNOG annotations

part.2 <- function(genome){ 
  codex <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                           genome,
                           "/trycycler/eggNOG_annotations/",
                           genome,
                           "_eggNOG_codex-3.txt"),
                    sep = " ", header = F, row.names = 1)
  
  new2 <- cbind(data[sort(row.names(data)),],
                tpm3[sort(row.names(tpm3)),],
                dta3[sort(row.names(dta3)),],
                codex[sort(row.names(codex)),])
  assign("new2",new2, envir = globalenv())
  
  colnames(new2)[109] <- "eggNOG_alt_name"
  
  new2$eggNOG_name <- row.names(new2)
  
  new2 <- new2[,c("eggNOG_name","eggNOG_alt_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",
                  "TPM-C1B10F1","TPM-C1B10F2","TPM-C1B10F3",
                  "TPM-C1B11F1","TPM-C1B11F2","TPM-C1B11F3",
                  "TPM-C1B12F1","TPM-C1B12F2","TPM-C1B12F3",
                  "TPM-C1B1F1","TPM-C1B1F2","TPM-C1B1F3",
                  "TPM-C1B2F1","TPM-C1B2F2","TPM-C1B2F3",
                  "TPM-C1B3F1","TPM-C1B3F2","TPM-C1B3F3",
                  "TPM-C1B4F1","TPM-C1B4F2","TPM-C1B4F3",
                  "TPM-C1B5F1","TPM-C1B5F2","TPM-C1B5F3",
                  "TPM-C1B6F1","TPM-C1B6F2","TPM-C1B6F3",
                  "TPM-C1B7F1","TPM-C1B7F2","TPM-C1B7F3",
                  "TPM-C1B8F1","TPM-C1B8F2","TPM-C1B8F3",
                  "TPM-C1B9F1","TPM-C1B9F2","TPM-C1B9F3",
                  "TPM-C2B1F1","TPM-C2B1F2","TPM-C2B1F3",
                  "TPM-C2B2F1","TPM-C2B2F2","TPM-C2B2F3",
                  "TPM-C2B3F1","TPM-C2B3F2","TPM-C2B3F3",
                  "TPM-C2B4F1","TPM-C2B4F2","TPM-C2B4F3",
                  "TPM-C2B5F1","TPM-C2B5F2","TPM-C2B5F3",
                  "COUNT-1B10F1","COUNT-1B10F2","COUNT-1B10F3",
                  "COUNT-1B11F1","COUNT-1B11F2","COUNT-1B11F3",
                  "COUNT-1B12F1","COUNT-1B12F2","COUNT-1B12F3",
                  "COUNT-1B1F1","COUNT-1B1F2","COUNT-1B1F3",
                  "COUNT-1B2F1","COUNT-1B2F2","COUNT-1B2F3",
                  "COUNT-1B3F1","COUNT-1B3F2","COUNT-1B3F3",
                  "COUNT-1B4F1","COUNT-1B4F2","COUNT-1B4F3",
                  "COUNT-1B5F1","COUNT-1B5F2","COUNT-1B5F3",
                  "COUNT-1B6F1","COUNT-1B6F2","COUNT-1B6F3",
                  "COUNT-1B7F1","COUNT-1B7F2","COUNT-1B7F3",
                  "COUNT-1B8F1","COUNT-1B8F2","COUNT-1B8F3",
                  "COUNT-1B9F1","COUNT-1B9F2","COUNT-1B9F3",
                  "COUNT-2B1F1","COUNT-2B1F2","COUNT-2B1F3",
                  "COUNT-2B2F1","COUNT-2B2F2","COUNT-2B2F3",
                  "COUNT-2B3F1","COUNT-2B3F2","COUNT-2B3F3",
                  "COUNT-2B4F1","COUNT-2B4F2","COUNT-2B4F3",
                  "COUNT-2B5F1","COUNT-2B5F2","COUNT-2B5F3")]
  # Load the eggNOG annotation file that was parsed from the "DESeq_diff_expr_from_eggNOG.sh" script
  anno <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                          genome,
                          "/trycycler/eggNOG_annotations/",
                          genome,
                          "_trans_annotations_3.txt"), sep = "\t", header = T)
  
  all <- left_join(new2, anno, by = join_by(eggNOG_alt_name == X.query))
  
  output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                             genome,"/trycycler/eggNOG_annotations/",
                             genome,
                             "_ALL_annos_w-DESeq.txt"), "wb")
  write.table(all, file = output.file, row.names = F, col.names = T, sep = "\t")
  close(output.file)
}


total.diff.analysis.76 <- function(genome){
  packages <- c("dplyr","DESeq2")
  sapply(c(packages), require, character.only = TRUE)
  dta <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         genome,
                         "/",
                         genome,
                         "_expt.ct.txt"), sep = "\t",
                  header=TRUE, row.names = "gene_id")
  ###-------- define data structure --------###
  samples <- data.frame(row.names=c("C1B10F1","C1B10F2","C1B10F3",
                                    "C1B11F1","C1B11F2","C1B11F3",
                                    "C1B12F1","C1B12F2","C1B12F3",
                                    "C1B1F1","C1B1F2","C1B1F3",
                                    "C1B2F1","C1B2F2","C1B2F3",
                                    "C1B3F1","C1B3F2","C1B3F3",
                                    "C1B4F1","C1B4F2","C1B4F3",
                                    "C1B5F1","C1B5F2","C1B5F3",
                                    "C1B6F1","C1B6F2","C1B6F3",
                                    "C1B7F1","C1B7F2","C1B7F3",
                                    "C1B8F1","C1B8F2","C1B8F3",
                                    "C1B9F1","C1B9F2","C1B9F3",
                                    "C2B1F1","C2B1F2","C2B1F3",
                                    "C2B2F1","C2B2F2","C2B2F3",
                                    "C2B3F1","C2B3F2","C2B3F3",
                                    "C2B4F1","C2B4F2","C2B4F3",
                                    "C2B5F1","C2B5F2","C2B5F3"), 
                        condition=as.factor(c(rep("In",45),
                                              # rep("C1B12",3),rep("C2B1",3))))
                                              rep("Out",6))))
  ### Get all samples in the same order
  all(rownames(samples) %in% colnames(dta))
  all(rownames(samples) == colnames(dta))
  dta <- dta[, rownames(samples)]
  all(rownames(samples) == colnames(dta))
  ### Make data integers
  dta2 <- as.data.frame(lapply(dta, as.integer))
  rownames(dta2) <- rownames(dta)
  colSums(dta)
  ###-------- create deseq object and process --------###
  # dta2["fakey",] <- rep(1, ncol(dta2))
  # 
  # cds <- DESeqDataSetFromMatrix(countData = dta2, colData=samples, design=~condition)
  # cds1 <- DESeq(cds)
  # 
  # ref <- results(cds1)
  # ref[order(ref$padj),]
  # 
  # c1b11vOut <- results(cds1, contrast=c("condition","In","Out"))
  # data <- as.data.frame(c1b11vOut[order(c1b11vOut$padj),])
  
  data <- data.frame(row.names = row.names(dta2),
                     baseMean = rep(0,nrow(dta2)),
                     log2FoldChange = rep(0,nrow(dta2)),
                     lfcSE = rep(0,nrow(dta2)),
                     stat = rep(0,nrow(dta2)),
                     pvalue =rep(0,nrow(dta2)),
                     padj=rep(0,nrow(dta2)))
  # View(as.data.frame(c1b11vOut[order(c1b11vOut$padj),]))
  
  output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                             genome,
                             "/trycycler/eggNOG_annotations/",
                             genome,
                             "_ALL_DESeq.txt"), "wb")
  write.table(data, file = output.file, row.names = T, col.names = T)
  close(output.file)
  
  
  ids <- row.names(data) %>% as.factor()
  ids <- as.factor(ids)
  
  
  output.file <- file(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                             genome,
                             "/trycycler/eggNOG_annotations/",
                             genome,
                             "_ALL_trans.txt"), "wb")
  write.table(ids, file = output.file, row.names = F, col.names = F)
  close(output.file)
  
  #####
  # Add TPM data to the "data" data frame
  tpm <- read.csv(paste0("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/",
                         genome,
                         "/",
                         genome,
                         "_TPM.txt"), 
                  sep = "\t",
                  header=TRUE, row.names = "gene_id")
  tpm2 <- tpm[,c("C1B10F1","C1B10F2","C1B10F3",
                 "C1B11F1","C1B11F2","C1B11F3",
                 "C1B12F1","C1B12F2","C1B12F3",
                 "C1B1F1","C1B1F2","C1B1F3",
                 "C1B2F1","C1B2F2","C1B2F3",
                 "C1B3F1","C1B3F2","C1B3F3",
                 "C1B4F1","C1B4F2","C1B4F3",
                 "C1B5F1","C1B5F2","C1B5F3",
                 "C1B6F1","C1B6F2","C1B6F3",
                 "C1B7F1","C1B7F2","C1B7F3",
                 "C1B8F1","C1B8F2","C1B8F3",
                 "C1B9F1","C1B9F2","C1B9F3",
                 "C2B1F1","C2B1F2","C2B1F3",
                 "C2B2F1","C2B2F2","C2B2F3",
                 "C2B3F1","C2B3F2","C2B3F3",
                 "C2B4F1","C2B4F2","C2B4F3",
                 "C2B5F1","C2B5F2","C2B5F3")] 
  tpm3 <- tpm2[row.names(tpm2) %in% row.names(data),]
  names(tpm3) <- c("TPM-C1B10F1","TPM-C1B10F2","TPM-C1B10F3",
                   "TPM-C1B11F1","TPM-C1B11F2","TPM-C1B11F3",
                   "TPM-C1B12F1","TPM-C1B12F2","TPM-C1B12F3",
                   "TPM-C1B1F1","TPM-C1B1F2","TPM-C1B1F3",
                   "TPM-C1B2F1","TPM-C1B2F2","TPM-C1B2F3",
                   "TPM-C1B3F1","TPM-C1B3F2","TPM-C1B3F3",
                   "TPM-C1B4F1","TPM-C1B4F2","TPM-C1B4F3",
                   "TPM-C1B5F1","TPM-C1B5F2","TPM-C1B5F3",
                   "TPM-C1B6F1","TPM-C1B6F2","TPM-C1B6F3",
                   "TPM-C1B7F1","TPM-C1B7F2","TPM-C1B7F3",
                   "TPM-C1B8F1","TPM-C1B8F2","TPM-C1B8F3",
                   "TPM-C1B9F1","TPM-C1B9F2","TPM-C1B9F3",
                   "TPM-C2B1F1","TPM-C2B1F2","TPM-C2B1F3",
                   "TPM-C2B2F1","TPM-C2B2F2","TPM-C2B2F3",
                   "TPM-C2B3F1","TPM-C2B3F2","TPM-C2B3F3",
                   "TPM-C2B4F1","TPM-C2B4F2","TPM-C2B4F3",
                   "TPM-C2B5F1","TPM-C2B5F2","TPM-C2B5F3")
  
  # Add the count Data 
  dta2 <- dta[,c("C1B10F1","C1B10F2","C1B10F3",
                 "C1B11F1","C1B11F2","C1B11F3",
                 "C1B12F1","C1B12F2","C1B12F3",
                 "C1B1F1","C1B1F2","C1B1F3",
                 "C1B2F1","C1B2F2","C1B2F3",
                 "C1B3F1","C1B3F2","C1B3F3",
                 "C1B4F1","C1B4F2","C1B4F3",
                 "C1B5F1","C1B5F2","C1B5F3",
                 "C1B6F1","C1B6F2","C1B6F3",
                 "C1B7F1","C1B7F2","C1B7F3",
                 "C1B8F1","C1B8F2","C1B8F3",
                 "C1B9F1","C1B9F2","C1B9F3",
                 "C2B1F1","C2B1F2","C2B1F3",
                 "C2B2F1","C2B2F2","C2B2F3",
                 "C2B3F1","C2B3F2","C2B3F3",
                 "C2B4F1","C2B4F2","C2B4F3",
                 "C2B5F1","C2B5F2","C2B5F3")] 
  dta3 <- dta2[row.names(dta2) %in% row.names(data),]
  names(dta3) <- c("COUNT-1B10F1","COUNT-1B10F2","COUNT-1B10F3",
                   "COUNT-1B11F1","COUNT-1B11F2","COUNT-1B11F3",
                   "COUNT-1B12F1","COUNT-1B12F2","COUNT-1B12F3",
                   "COUNT-1B1F1","COUNT-1B1F2","COUNT-1B1F3",
                   "COUNT-1B2F1","COUNT-1B2F2","COUNT-1B2F3",
                   "COUNT-1B3F1","COUNT-1B3F2","COUNT-1B3F3",
                   "COUNT-1B4F1","COUNT-1B4F2","COUNT-1B4F3",
                   "COUNT-1B5F1","COUNT-1B5F2","COUNT-1B5F3",
                   "COUNT-1B6F1","COUNT-1B6F2","COUNT-1B6F3",
                   "COUNT-1B7F1","COUNT-1B7F2","COUNT-1B7F3",
                   "COUNT-1B8F1","COUNT-1B8F2","COUNT-1B8F3",
                   "COUNT-1B9F1","COUNT-1B9F2","COUNT-1B9F3",
                   "COUNT-2B1F1","COUNT-2B1F2","COUNT-2B1F3",
                   "COUNT-2B2F1","COUNT-2B2F2","COUNT-2B2F3",
                   "COUNT-2B3F1","COUNT-2B3F2","COUNT-2B3F3",
                   "COUNT-2B4F1","COUNT-2B4F2","COUNT-2B4F3",
                   "COUNT-2B5F1","COUNT-2B5F2","COUNT-2B5F3")
  # Remove the "Fakey" data
  data <- data[! (row.names(data) %in% "fakey"),]
  
  assign("data", data, envir = globalenv())
  assign("tpm3", tpm3, envir = globalenv())
  assign("dta3", dta3, envir = globalenv())
}
