#!/bin/bash

/data/app/Flye-2.9.2/bin/flye --nano-hq $1 --threads 16 --out-dir flye_assembly_01 && cp flye_assembly_01/assembly.fasta ../assemblies/flye_assembly_01.fasta && cp flye_assembly_01/flye_assembly_graph.gfa ../assemblies/flye_assembly_01.gfa && rm -r flye_assembly_01

~/miniasm_and_minipolish.sh $1 16 > ../assemblies/miniasm_assembly_02.gfa 

~/any2fasta/any2fasta ../assemblies/miniasm_assembly_02.gfa > ../assemblies/miniasm_assembly_02.fasta

raven --threads 16 --disable-checkpoints --graphical-fragment-assembly ../assemblies/raven_assembly_03.gfa $1 > ../assemblies/raven_assembly_03.fasta

# Cluster all of the assembly reads together
trycycler cluster --assemblies ../assemblies/*.fasta --reads $1 --out_dir trycycler_clusters

echo "This is how many contigs are in your Flye assembly"
grep -Fo ">" ../assemblies/flye_assembly_01.fasta | wc -l 

echo "This is how many contigs are in your Miniasm assembly"
grep -Fo ">" ../assemblies/miniasm_assembly_02.fasta | wc -l 

echo "This is how many contigs are in your Raven assembly"
grep -Fo ">" ../assemblies/raven_assembly_03.fasta | wc -l 
