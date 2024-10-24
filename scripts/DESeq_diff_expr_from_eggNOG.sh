#!/bin/bash

# $1 = OXR-X

# Remove the quotations from the file written from R.
sed -ie 's/"//g' "${1}_padj_trans.txt"

# copy paste from the padj_trans file to new file in Notepad++ to avoid windows v linux issues.
grep -wf "${1}_padj_trans.txt" ../trycycler/eggNOG_annotations/out.emapper.decorated.gff > "${1}_diff_expr.txt"

# Take the OXR-9_diff_expr.txt file and keep the 9th column to query for em_KEGG_Pathway
awk 'BEGIN {FS="\t"} {print $9}' "${1}_diff_expr.txt" > "${1}_diff_expr_2.txt"

# take the em_ID, which corresponds to another alternative identification in the out.emapper.annotation table
awk 'BEGIN {FS=";"} {print $14}' "${1}_diff_expr_2.txt" > em_ID.txt

# remove blank lines
sed -ie '/^\s*$/d' em_ID.txt

# remove teh em_ID= from each lines
sed -i 's/em_ID=//' em_ID.txt

# make a file with em_ID and the original id
awk 'BEGIN {FS=";"} {print $1, $14}' "${1}_diff_expr_2.txt" > "${1}_eggNOG_codex.txt"
sed -i 's/em_ID=//' "${1}_eggNOG_codex.txt"
sed -i 's/ID=//' "${1}_eggNOG_codex.txt"


# take the names in the em_ID.txt file and match them to the annotations file
grep -wf em_ID.txt ../trycycler/eggNOG_annotations/out.emapper.annotations > diff_expr_annotations.txt

# take the 5th row that from the out.emapper.annotations file
awk 'NR == 5' ../trycycler/eggNOG_annotations/out.emapper.annotations > headers.txt

# concatenate the headers with the diff_expr_annotations.txt
cat headers.txt diff_expr_annotations.txt > "${1}_diff_expr_annotations.txt"

rm "${1}_padj_trans.txt"
rm "${1}_padj_trans.txte"
rm "${1}_diff_expr.txt"
rm "${1}_diff_expr_2.txt"
rm em_ID.txt
rm em_ID.txte
rm diff_expr_annotations.txt
rm headers.txt
