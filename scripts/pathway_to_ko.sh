#!/bin/bash

# run this code in the "kegg_query_all_transcripts" directory after running kegg_querying-2.bash script

cp pathway-to-ko.txt pathway-to-ko-2.txt
# add ";" to the end of each line
sed -i 's/$/;/' pathway-to-ko-2.txt
# split file into columns listing the KEGG pathway map ID to each of the possible KO ids
awk '{for(a=2;a<=NF;){printf"%s%c\n",$1" "$a,a++-NF?";":""}}' pathway-to-ko-2.txt > Finaloutput.txt
# remove the ";"
sed -i 's/;//g' Finaloutput.txt
rm pathway-to-ko-2.txt
