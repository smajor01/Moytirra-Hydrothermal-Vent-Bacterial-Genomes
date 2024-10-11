#!/bin/bash

# $1 = OXR-X

# Take the OXR-9_diff_expr.txt file and keep the 9th column to query for em_KEGG_Pathway
awk 'BEGIN {FS="\t"} {print $9}' trycycler/eggNOG_annotations/out.emapper.decorated.gff > trycycler/eggNOG_annotations/tmp.txt

# take the em_ID, which corresponds to another alternative identification in the out.emapper.annotation table
awk 'BEGIN {FS=";"} {print $14}' trycycler/eggNOG_annotations/tmp.txt > trycycler/eggNOG_annotations/em_ID.txt

# remove blank lines ; remove th em_ID= from each lines
sed -i '/^\s*$/d' trycycler/eggNOG_annotations/em_ID.txt ; sed -i 's/em_ID=//' trycycler/eggNOG_annotations/em_ID.txt

# make a file with em_ID and the original id
awk 'BEGIN {FS=";"} {print $1, $14}' trycycler/eggNOG_annotations/tmp.txt > trycycler/eggNOG_annotations/"${1}_eggNOG_codex.txt"

# Remove erroneous stuff
sed -i 's/em_ID=//' trycycler/eggNOG_annotations/"${1}_eggNOG_codex.txt" ; sed -i 's/ID=//' trycycler/eggNOG_annotations/"${1}_eggNOG_codex.txt" ; sed -i '/^\s*$/d' trycycler/eggNOG_annotations/"${1}_eggNOG_codex.txt"


rm trycycler/eggNOG_annotations/tmp.txt
rm trycycler/eggNOG_annotations/em_ID.txt

