threads=16  # change as appropriate for your system
mkdir assemblies

flye --nano-hq read_subsets/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01
~/miniasm_and_minipolish.sh read_subsets/sample_02.fastq "$threads" > assemblies/assembly_02.gfa && ~/any2fasta/any2fasta assemblies/assembly_02.gfa > assemblies/assembly_02.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa read_subsets/sample_03.fastq > assemblies/assembly_03.fasta

flye --nano-hq read_subsets/sample_04.fastq --threads "$threads" --out-dir assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && cp assembly_04/assembly_graph.gfa assemblies/assembly_04.gfa && rm -r assembly_04
~/miniasm_and_minipolish.sh read_subsets/sample_05.fastq "$threads" > assemblies/assembly_05.gfa && ~/any2fasta/any2fasta assemblies/assembly_05.gfa > assemblies/assembly_05.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_06.gfa read_subsets/sample_06.fastq > assemblies/assembly_06.fasta

flye --nano-hq read_subsets/sample_07.fastq --threads "$threads" --out-dir assembly_07 && cp assembly_07/assembly.fasta assemblies/assembly_07.fasta && cp assembly_07/assembly_graph.gfa assemblies/assembly_07.gfa && rm -r assembly_07
~/miniasm_and_minipolish.sh read_subsets/sample_08.fastq "$threads" > assemblies/assembly_08.gfa && ~/any2fasta/any2fasta assemblies/assembly_08.gfa > assemblies/assembly_08.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_09.gfa read_subsets/sample_09.fastq > assemblies/assembly_09.fasta

flye --nano-hq read_subsets/sample_10.fastq --threads "$threads" --out-dir assembly_10 && cp assembly_10/assembly.fasta assemblies/assembly_10.fasta && cp assembly_10/assembly_graph.gfa assemblies/assembly_10.gfa && rm -r assembly_10
~/miniasm_and_minipolish.sh read_subsets/sample_11.fastq "$threads" > assemblies/assembly_11.gfa && ~/any2fasta/any2fasta assemblies/assembly_11.gfa > assemblies/assembly_11.fasta
raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_12.gfa read_subsets/sample_12.fastq > assemblies/assembly_12.fasta