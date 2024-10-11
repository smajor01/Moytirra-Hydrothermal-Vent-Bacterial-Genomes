#!/bin/bash

# $1 = R1 illumina
# $2 = R2 illumina
# $3 = OXR-X_polca_polypolish_trycycler.fasta
# $4 = OXR-X
# must be in trycycler folder

bwa index trycycler_consensus.fasta && bwa mem -t 16 -a trycycler_consensus.fasta $1 > alignments_1.sam && bwa mem -t 16 -a trycycler_consensus.fasta $2 > alignments_2.sam && /data/app/polypolish-v0.5.0/polypolish trycycler_consensus.fasta alignments_1.sam alignments_2.sam > polypolish_trycycler.fasta

/data/app/MaSuRCA-4.1.0/bin/polca.sh -a polypolish_trycycler.fasta -r "$1 $2" -t 16 -m 1G

mv *.PolcaCorrected.fa $3
rm *sam
cp $3 ../assemblies/

cd ../assemblies/
quast --glimmer *fasta
checkm lineage_wf -t 8 -x fasta ./ checkm > CheckM_table.txt


cd ../trycycler
busco -i $3 -o busco -l ../../OXR-11/busco_downloads/lineages/bacteria_odb10/ --mode genome --long -c 12
prokka --outdir prokka_annotation $3 --prefix $4
