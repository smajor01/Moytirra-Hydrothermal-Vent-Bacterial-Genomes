oxr.9.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-9/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.11.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-11/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.76.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-76/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.85.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-85/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.96.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-96/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.134.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-134/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.137.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-137/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.159.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-159/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.189.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-189/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.199.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-199/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.203.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-203/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")
oxr.209.ddh <- read.csv("C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/OXR-209/TYGS_taxonomy/pairwise comparison genome vs type strain.csv")

ddh <- rbind(oxr.9.ddh,oxr.11.ddh,oxr.76.ddh,oxr.85.ddh,
      oxr.96.ddh,oxr.134.ddh,oxr.137.ddh,oxr.159.ddh,
      oxr.189.ddh,oxr.199.ddh,oxr.203.ddh,oxr.209.ddh)

write.csv(file = "C:/Users/SamuelMajor/OneDrive - Gloucester Marine Genomics Institute/Desktop/OXR genomes/TYGS_ddh_results.csv",
          ddh)
