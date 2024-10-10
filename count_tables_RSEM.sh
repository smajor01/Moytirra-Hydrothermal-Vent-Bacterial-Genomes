#!/bin/bash

# $1 = OXR-X
# Agglomerate all of the expected_counts from the .gene.results files into a single table

cut -f 1 C1B1F1/*genes.results > gene_id.tmp

mkdir results.tmp/

cp C*/*genes.results results.tmp/

for files in results.tmp/*genes.results
do

        sample="${files##*_}"
        sample2="${sample%%.*}"

        cut -f 5 $files > "${files}_expt.ct.tmp"

        sed -e "1s/expected_count/${sample2}/" "${files}_expt.ct.tmp" > "${files}_list.tmp"

        paste gene_id.tmp results.tmp/*list.tmp > "${1}_expt.ct.txt"


        cut -f 6 $files > "${files}_tpm.tmp"

        sed -e "1s/TPM/${sample2}/" "${files}_tpm.tmp" > "${files}_list2.tmp"

        paste gene_id.tmp results.tmp/*list2.tmp > "${1}_TPM.txt"


        cut -f 7 $files > "${files}_fpkm.tmp"

        sed -e "1s/FPKM/${sample2}/" "${files}_fpkm.tmp" > "${files}_list3.tmp"

        paste gene_id.tmp results.tmp/*list3.tmp > "${1}_FPKM.txt"

done

rm -r results.tmp/
rm *tmp

