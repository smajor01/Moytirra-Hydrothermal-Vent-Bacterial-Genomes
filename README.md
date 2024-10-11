# Moytirra-Hydrothermal-Vent-Bacterial-Genomes

Genomes for this project can be found at NCBI BioProject accession number PRJNA1046213.

## Sample collection
Sample collection is described in Polinski et al. 2023 (DOI 10.3389/fmars.2023.1219784)

## Microbe Culturing
Filters were submerged in 25% glycerol in sterile filtered seawater for cryopreservation. The filters  were thawed on ice and vortexed 
Bacteria were cultured on a variety of media composed of either 2% NaCl or artifical seawater (ASW).
Water was supplemented with 50% or 10% of yeast extract and peptone in DIFCO Marine Broth
The media was solidified with either gellan gum (Phytagel) or agar.
All media contained 200 µg/mL of Cycloheximide.
Isolated cultures were grown in broth media for 3-18 days before being cryopreserved in 25% glycerol. Separate aliquots were collected and frozen prior to 16S identification.

## DNA Extraction and Sequencing
### DNA Extraction and 16S rRNA gene identification: 
Frozen aliquots were thawed and lysed in ultrapure water at 95°C for 15 minutes then immediately placed on ice. The 16S rRNA gene was amplified with forward primer v4_515F and reverse primer v4_806R, both containing Illumina linking sequences then sequenced according to Polinski et al. 2023 (DOI 10.3389/fmars.2023.1219784) on a 2x250 MiSeq run with v2 reagent kit (Illumina, Inc., San Diego, CA, USA).

Read trimming and ASV calling was performed with the DADA2 v1.30 R package, removing the first 10 bases from each read and removing the last 50 and 20 bases of the forward and reverse reads, respectively, allowing a maximum expected error rate of 2. Taxonomy was assigned with the Silva v138.1 database. To determine the presence of the cultured isolates within the Moytirra HVP plume community, ASVs were aligned against metabarcoding data of the Moytirra HVP (DOI 10.3389/fmars.2023.1219784) using BLASTn with default parameters.

### DNA Extraction and Whole Genome Sequencing of 12 Isolates:
Frozen 20 µL aliquots of each of the 12 isolates was thawed at room temperature and extracted using a modified CTAB protocol.

**Illumina Short Read Sequencing**: 1 µg of DNA was sheared to 500 bp, size selected, end-repaired, A-tailed, ligated with Illumina indexes, size selected, then sent to the University of Connectericut Center for Genomic Innovation for paired-end sequencing on an Illumina NovaSeq.

**Oxford Nanopore Long Read Sequencing**: DNA was size selected for fragments >1.5 kbp and 1 µg of high molecular weight was preparred following the Oxford Nanopore protocol for Native Barcoding Genomic DNA with the Ligation Sequencing kit (SQK-LSK109) and Native Barcoding Expansion kit (EXP-NBD104). Libraries were sequenced on a MinION R9.4.1 flow cell, basecalling was performed with Guppy v6.4.8 using the fast basecalling algorithm (community.nanoportech.com/download). The resulting sequences were quality filtered using the NanoPack programs’ nanofilt (25), retaining sequences greater than 150 bp with a FAST5 q score of 10.

## Genome Assembly, Annotation, and Taxonomy
Hybrid de novo assembly of the bacterial genomes was achieved using Trycycler v0.5.4 (https://github.com/rrwick/Trycycler) by clustering assemblies from 6 different pipelines using their default parameters: MaSuRCA v4.1.0 hybrid assembler (https://github.com/alekseyzimin/masurca); SPAdes v3.15.5 hybrid assembler (https://github.com/ablab/spades); Canu v2.1.1 long-read assembler (https://github.com/marbl/canu) combined with short-read polishing with POLCA (https://github.com/alekseyzimin/masurca); Flye v2.9.2 long-read assembler (https://github.com/mikolmogorov/Flye); Miniasm v0.3 long-read assembler (https://github.com/lh3/miniasm) followed by short-read polishing with Minipolish v0.1.3 (https://github.com/rrwick/Minipolish); and Raven v1.8.1 (https://github.com/lbcb-sci/raven). 

Clusters of contigs were manually curated by discarding clusters with 3 or fewer contigs and then reconciled into circular elements. When reconciliation failed, contigs required manual trimming, particularly in instances where contigs within a cluster were approximately double the size of contigs from the same cluster. Curation was performed by searching for sequences that appear in the beginning or end that are repeated within the contig and deleting all bases after or before the repeated sequence, respectively. Trycycler used Muscle (https://github.com/rcedgar/muscle) to generate a multiple sequence alignment of each cluster for which the entire read set can be aligned to generate a consensus sequence for each replicon. Illumina short reads were aligned to the Trycycler consensus sequence with BWA v0.7.17 (37) and polished with PolyPolish v0.5.0 (https://github.com/rrwick/Polypolish) and POLCA (https://github.com/alekseyzimin/masurca). The resulting consensus genomes for each isolate were assessed for completion and contamination with BUSCO v5.4.5 with the bacteria_odb10 lineage and CheckM v1.2.2 (https://github.com/Ecogenomics/CheckM). 

Annotations were performed with Prokka v1.14.5 (https://github.com/tseemann/prokka) and eggNOG-mapper v2.1.9 (http://eggnog-mapper.embl.de/). The resulting KEGG orthology terms were aggregated by their functions and counts were visualized with the ‘pheatmap’ v1.0.12 package in R. Genome-based taxonomy was performed with the largest replicon (containing the 16S rRNA genes) using the Type (Strain) Genome Server (TYGS) to predict the closest type-strain by providing a digital DNA-DNA hybridization (dDDH), where <70% dDDH predicts a novel species. Information on nomenclature and associated taxonomic literature was provided by the List of Prokaryotic names with Standing in Nomenclature (LPSN) (44). 

The putatively novel species, Thalassobaculum sp. OXR-137, Sulfitobacter sp. OXR-159, Idiomarina sp. OXR-189, and Christiangramia sp. OXR-203, were further compared with their closest predicted type-strains and a complete genome of a different species from the type-strain with high 16S rRNA sequence identity from GenBank.  The dDDH between the curated genomes and the novel isolates was calculated with the Genome-to-Genome Distance Calculator v3.0 and the average nucleotide identity (ANI) was calculated with FastANI v1.3.3 (https://github.com/ParBLiSS/FastANI). BLASTx aligned proteins with an E-value cutoff of 1.0e-20. Synteny analysis was performed using MCscan v1.3.8 (https://github.com/tanghaibao/mcscan), and non-syntenic genomic islands (GIs) were identified as 10 or more non-syntenic genes separated by fewer than 10 syntenic genes. GI annotations and functions were supplemented with NCBI’s PGAP pipeline and GhostKoala (https://www.kegg.jp/ghostkoala/) against the “genus_prokaryotes+virus” database, respectively . To putatively explain the occurrence of GIs, horizontal gene transfer events were predicted with Alien Hunter v1.7, mobile genetic elements were predicted with mobileOG-db v1.6, and prophage elements were predicted with Phigaro v2.3.0 and VirSorter2 v2.2.4 through the Proksee webtool (https://proksee.ca/). 

## Biosynthetic Gene Cluster Prediction
Biosynthetic gene clusters (BGCs) were first predicted using DeepBGC v0.1.23 (https://github.com/Merck/deepbgc) and the resulting JSON file was uploaded with its respective genome FASTA file to the antiSMASH v7.0.1 webserver (https://antismash.secondarymetabolites.org/#!/start). The resulting analysis was then compiled into a table using antismash_coverter.py (https://github.com/sandragodinhosilva/bgc-analysis/tree/main) and plotted in R with ggplot2.

## Transcriptional Activity within the Vent Plume
The activity of the isolates within the plume was assessed by aligning metatranscriptomic reads from Polinski et al. (DOI 10.3389/fmars.2023.1219784) to each genome with Bowtie2 (https://github.com/BenLangmead/bowtie2), and count tables generated with RSEM v1.3.3 (https://github.com/deweylab/RSEM). Resulting count tables were evaluated for differential expression with DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) by comparing samples taken within the plume against samples taken outside of the plume.
To assess overall metabolic activity, normalized read counts (TPM) for all genes identified with KEGG orthology terms were first binned by the KEGG orthology term and then KEGG pathway to evaluate pathway-level responses within the plume (summed normalized read counts of samples within the plume) vs outside the plume (summed normalized read counts of samples out of the plume). The resulting table for each genome was then combined and ordered by abundance, and the top 30 pathways receiving the most counts were visualized with the ‘pheatmap’ v1.0.12 package in R. 

