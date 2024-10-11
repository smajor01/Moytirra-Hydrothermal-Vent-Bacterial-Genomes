#!/bin/bash

# Get the non-syntenic regions from novel genome along with annotations

# $1 = OXR-x
# $2 = Syntenous genome

python3 -m jcvi.compara.synteny mcscan "${1}.bed" "${2}.${1}.lifted.anchors" --iter=1 -o "${2}.${1}.i1.blocks"

# take sequence names that don't have a syntenous genes
awk '$2 == "."' "${2}.${1}.i1.blocks" > "${2}.dot.test"

# take the annotations of the proteins from the .ffn file
grep ">" "${1}.ffn" > "${1}.annos.txt"
sed -i 's/>//g' "${1}.annos.txt"

# take the annotations that match the non-syntenic regions.
awk '{print $1}' "${2}.dot.test" | grep -Ff - "${1}.annos.txt" > "${2}.${1}.no-synteny.annos.txt"
# add a new first column with some identifying information
awk -v new_column="'${1}_${2}'" '{print new_column, $0}' "${2}.${1}.no-synteny.annos.txt" > tmp && mv tmp "${2}.${1}.no-synteny.annos.txt"
