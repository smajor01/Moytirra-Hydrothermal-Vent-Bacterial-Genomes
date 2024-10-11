#!/bin/bash
# Kegg querying
# Be in the trycycler/eggNOG_annotations folder
# Using the files generated from the R scripts "OXR_kos and gene ID all transcripts.Rmd"

# $1 = OXR-X_ko.txt

mkdir kegg_query_all_transcripts_ALL

# First remove the erroneous ""
sed -i 's/"//g' $1

# Then subset only the gene ID and the Kegg ortholog
awk '{print $1, $105}' $1 > "${1}_ALL_simple.txt"

# Remove lines with dashes
sed -i '/-/d' "${1}_ALL_simple.txt"

# Remove rows that are empty in column 2 
awk '$2 !=""' "${1}_ALL_simple.txt" > text.txt

# make spaces into tabs
awk -v OFS="\t" '$1=$1' text.txt > ALL_ko_list.txt

# remove gene_id column
awk '{print $2}' ALL_ko_list.txt > tmp; mv tmp ALL_ko_list.txt

# remove first row with headers
sed -i '1d' ALL_ko_list.txt

cd kegg_query_all_transcripts_ALL

python3 ../../../../../scripts/kegg_query.py ../ALL_ko_list.txt

# Make a copy so not to mess with the original.
cp pathway-to-ko.txt pathway-to-ko-2.txt
# add ";" to the end of each line
sed -i 's/$/;/' pathway-to-ko-2.txt
# split file into columns listing the KEGG pathway map ID to each of the possible KO ids
awk '{for(a=2;a<=NF;){printf"%s%c\n",$1" "$a,a++-NF?";":""}}' pathway-to-ko-2.txt > ALL_Finaloutput.txt
# remove the ";"
sed -i 's/;//g' ALL_Finaloutput.txt
rm pathway-to-ko-2.txt

cd ../

rm text.txt

