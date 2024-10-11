#!/bin/bash
# Kegg querying
# Using the files generated from the R scripts "OXR_kos and gene ID.Rmd"

# $1 = OXR-X_KEGG_egg.txt

mkdir kegg_query

# First remove the erroneous ""
sed -i 's/"//g' $1

# Then subset only the gene ID and the Kegg ortholog
awk '{print $1, $3}' $1 > "${1}_simple.txt"

# Remove lines with dashes
sed -i '/-/d' "${1}_simple.txt"

# Remove rows that are empty in column 2 
awk '$2 !=""' "${1}_simple.txt" > text.txt

# make spaces into tabs
awk -v OFS="\t" '$1=$1' text.txt > ko_list.txt

# remove gene_id column
awk '{print $2}' ko_list.txt > tmp; mv tmp ko_list.txt

# remove first row with headers
sed -i '1d' ko_list.txt

cd kegg_query

python3 ../../../../../scripts/kegg_query.py ../ko_list.txt

cd ../

rm text.txt
