#!/bin/bash

# $1 = ref/the final genome assembly
# $2 = OXR-X_ref
#
# Usage...
# ~/RSEM_OXR.sh ref/OXR-9_final_assembly.fsa OXR-9_ref

/data/app/gff-utilities/gffread-0.11.6/gffread -T ref/out.emapper.genepred.gff -o ref/out.emapper.genepred.gtf

sed -n '/CDS/!p' ref/out.emapper.genepred.gtf > ref/transcript.gtf
sed "s/transcript/exon/g" ref/transcript.gtf > ref/exon.gtf
sed "s/exon_id/transcript_id/g" ref/exon.gtf > ref/exon_2.gtf
cat ref/transcript.gtf ref/exon_2.gtf > ref/out.emapper.genepred_2.gtf

/data/app/RSEM-1.3.3/rsem-prepare-reference --gtf ref/out.emapper.genepred_2.gtf --bowtie2 --bowtie2-path /data/app/bowtie2-2.5.1/ $1 ref/$2

for file in /data/prj/ecosystem-diversity/OceanX/YEP1_2021/metatranscriptome/trimmed-reads/*.gz
do

    prefix="${file%R1*}"
    suffix="${file#*R1}"

        folder="${prefix##/data/prj/ecosystem-diversity/OceanX/YEP1_2021/metatranscriptome/trimmed-reads/}"
        folder2="${folder%%_*}"
        mkdir "${folder2}"
    if [ -f "${prefix}R2${suffix}" ]
    then
        /data/app/RSEM-1.3.3/rsem-calculate-expression -p 8 --paired-end --bowtie2 --bowtie2-path /data/app/bowtie2-2.5.1/ --estimate-rspd --append-names --output-genome-bam "${prefix}R1${suffix}" "${prefix}R2${suffix}" ref/$2 "${folder2}"/"${2%_*}_${folder2}"
    fi
done
