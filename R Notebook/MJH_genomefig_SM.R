###-------- load packages you need --------###
library(pheatmap)

###-------- load z-score function --------###
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

###-------- load data --------###
dta <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/OceanX Deep Water Culturing/Notes and Genome Assemblies/genomePathOut.csv", header=TRUE, row.names = 1)
annot <- read.csv("C:/Users/SamuelMajor/Box/Science/Microbe Repository/OceanX Deep Water Culturing/Notes and Genome Assemblies/GenomePathClass.csv", header=TRUE, row.names = 1)
names(annot) <- c("Secondary Classification","Classification")
###-------- remove rows with all zeros --------###
dta <- dta[rowSums(dta[])>0,]

j <- names(dta)

j <- gsub("\\.","-",j)
j <- gsub("_"," ",j)
names(dta) <- j
# pick top 25 rows with highest values
top <- dta[order(rowSums(dta), decreasing = TRUE)[1:30], ]

#scale the rows
dta_norm <- t(apply(top, 1, cal_z_score))

#custom annotation colors
###-------- set colors and other formatting --------###
ann_colors = list(`Secondary Classification` = c("Amino acid metabolism" = "pink",
                                                                     # "Biosynthesis of other secondary metabolites" = "#00CDCD",
                                                                     "Carbohydrate metabolism" = "#7FFFD4",
                                                                     # "Cell motility" = "#FF7F24",
                                                                     "Cellular community - prokaryotes" = "#0000FF",
                                                                     "Energy metabolism" = "#C0FF3E",
                                                                     # "Glycan biosynthesis and metabolism" = "#8B3A62",
                                                                     "Lipid metabolism" = "#FA8072",
                                                                     "Membrane transport" = "#BF3EFF",
                                                                     "Metabolism of cofactors and vitamins" = "#2E8B57",
                                                                     "Metabolism of other amino acids" = "#FFE1FF",
                                                                     "Metabolism of terpenoids and polyketides" = "royalblue4",
                                                                     "Nucleotide metabolism" = "springgreen",
                                                                     "Replication and repair" = "plum4",
                                                                     "Signal transduction" = "#EE1289",
                                                                     "Translation" = "turquoise"),
                                                                     # "Xenobiotics biodegradation and metabolism" = "orange"),
                                      # "Folding sorting and degradation", 
                                      # "Environmental adaptation",
                                      # "Transport and catabolism", 
                                      # "Transcription",
                                      # "Information processing in viruses",
                                      # "Signaling molecules and interaction" = "black"),
                  Classification = c("Cellular Processes"="#d81b60ff",
                                     "Environmental Information Processing"="#1e88e5ff",
                                     "Genetic Information Processing"="#ffc107ff",
                                     "Metabolism"="#00d133ff"))

name.order <- c('OXR-9 in','OXR-11 in','OXR-85 in','OXR-96 in',
                'OXR-134 in','OXR-137 in','OXR-159 in','OXR-189 in',
                'OXR-199 in','OXR-203 in','OXR-209 in',
                'OXR-9 out','OXR-11 out','OXR-85 out','OXR-96 out',
                'OXR-134 out','OXR-137 out','OXR-159 out','OXR-189 out',
                'OXR-199 out','OXR-203 out','OXR-209 out')

dta_norm.2 <- dta_norm[, match(name.order, colnames(dta_norm))]

k <- colnames(dta_norm.2)
k <- gsub(" in","",k)
k <- gsub(" out","",k)
colnames(dta_norm.2) <- k
###-------- generate pheatmap plot --------###
fig <- pheatmap(dta_norm.2,
         cluster_cols=F, 
         cluster_rows = T,
         gaps_col = c(11),
         annotation_col = primary.size, # This can be found in the R file "Finding_KEGG_functions_20240310.RMD"
         annotation_row = annot,
         annotation_colors=ann_colors,
         fontsize_row = 10, 
         show_rownames = T,
         show_colnames = T,
         cellwidth = 10,
         cellheight = 10,
         angle_col = 90,
         main = "In                             Out"
         )

#add title to heat color legend?
# g <- fig$figtable
# grid.newpage()
# grid.draw(fig)
# grid.text("Z-score", x=0.6, y=0.88, rot=0, gp=gpar(fontsize=10, fontface="bold"))

save_pheatmap_png(plot = fig,
                  filename = "C:/Users/SamuelMajor/Box/Science/Microbe Repository/OceanX Deep Water Culturing/Genome_MJH_trans_heatmap_SM.png",
                  width = 1700, height = 1000)

library(ggplot2)
###-------- load data --------###
df <- read.csv("input/genometotalcount.csv", header=TRUE)

b <- ggplot(data=df, aes(x=isolate, y=In)) +
  geom_bar(stat="identity")

