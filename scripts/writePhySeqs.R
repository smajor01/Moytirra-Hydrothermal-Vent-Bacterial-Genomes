writePhySeqs <- function(PhySeq, path, Sample.Name = "Sample") {
  taxa <- as.data.frame((tax_table(PhySeq)))
  otu <- as.data.frame(otu_table(PhySeq))
  tree <- phy_tree(PhySeq)
  seqs <- as.data.frame(refseq(PhySeq))
  write.csv(taxa, paste0(path, Sample.Name, "-taxa.csv"))
  write.csv(otu, paste0(path, Sample.Name,"-otu.csv"))
  write.tree(tree, paste0(path, Sample.Name, "-tree.tre"))
  Biostrings::writeXStringSet(refseq(PhySeq), 
                              paste0(path, Sample.Name, "-seqs.fasta", 
                              format = "fasta"))
}
