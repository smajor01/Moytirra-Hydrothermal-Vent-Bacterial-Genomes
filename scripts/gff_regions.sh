#!/bin/bash

# Use this file to get the start and stop of each gene from a gff file


# Take the 4th and 5th and 7th column which corresponds to the start and stop and orientation
awk '{print $4, $5, $7}' out.emapper.genepred.gff > regions.txt
# remove blank lines
sed -i '/^\s*$/d' regions.txt
sed -i '/seqnum/d' regions.txt
sed -i '/version/d' regions.txt

# take the 9th column which corresponds to the special information in the gff
awk '{print $9}' out.emapper.genepred.gff > tmp.id.txt
# take the em_ID, which corresponds to another alternative identification in the out.emapper.annotation table
awk 'BEGIN {FS=";"} {print $1}' tmp.id.txt > ID.txt
# Removingnthe ID=
sed -i 's/ID=//g' ID.txt
# Make the file tab delimited
sed -i 's/\\s+$/\t/g' ID.txt

# remove blank lines
sed -i '/^\s*$/d' ID.txt

paste ID.txt regions.txt > cds.regions.txt

