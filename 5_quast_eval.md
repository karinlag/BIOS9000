# Comparing assemblies to the reference

The Quast program can be used to generate similar metrics as the
assemblathon-stats script, plus some more and some visualisations.
QUAST works with or without a genome to compare it to. We will use it to
compare our assemblies to a reference. Thus, we need to give it a
reference genome and a reference annotation. We will also give it the
short reads, which the program will map to the assemblies using bwa.
This will show coverage, like IGV.

When running QUAST, we can look at one or more assemblies at the time.

`Quast` options:

* `-o`: output directory
* `--glimmer`: tell the program to predict genes
* `--rna-finding`: tell the program to find rRNAs
* `-1/2`: path to the reads.
* `-G`: File with positions of genes in the reference (see manual)
* `-t`: number of threads (cpu's) to use
* `sequence.fasta`: one or more files with assembled sequences
* `-l`:  comma-separated list of names for the assemblies, e.g. `"assembly 1,
assembly 2"` (in the same order as the sequence files)
* `--scaffolds`: input sequences are scaffolds, not contigs. They will be split
at 10 N's or more to analyse contigs ('broken' assembly)

[See the manual for information on the output of Quast]:
(http://quast.sourceforge.net/docs/manual.html#sec3)

### Getting the reference genomes

First, go to the `assembly` folder, and then into the `rawdata` folder.
Then link in the reference data:

```
ln -s /work/IN-BIOSx/data/assembly/GCA* .
```

### Class activity: comparing velvet and spades. 

Go to the `assembly` folder, make a folder called `quast` and move into it. Run this command:

```
quast -o spades_velvet \
-t 8 \
--glimmer \
-1 ../rawdata/SRR10015223_1_1mill_mut.fastq \
-2 ../rawdata/SRR10015223_2_1mill_mut.fastq \
-r ../rawdata/GCA_012098105.1_PDT000380306.1_genomic.fna \
-g ../rawdata/GCA_012098105.1_PDT000380306.1_genomic.gff \
PATH/TO/SPADES_ASM \
PATH/TO/VELVET_ASM \
-l "spades, velvet"
```
Note that the `--scaffold` option is not used here for simplification. Also,
if there are multiple assemblies being compared, make sure you name the
assemblies (`-l`) in the same order as you give them to quast!

### Quast output
Quast will produce several files in the output directory. Transfer the directory
to your local computer, and open the `report.html` in a browser.

Questions:
  * What is the N50 for these assemblies?
  * What is the NA50 for these genomes?
  * How many whole genes have been found? How many partial genes?
  * Have you found significantly fewer or more genes than in the reference?
  * What about rRNA genes?
  * How many misassmblies are there?
  * How much of the genomes are unaligned?
  * How many mismatches and indels do we have?

Click on the "Icarus contig browser"
  * Can you figure out what each of the panels show?
  * Can you see regions of your reference that are not in the assemblies?
  * Can you figure out what the various contig colors mean?

### Group activity: comparing all three assemblies

Go to the `assembly/quast` folder. Run this command.

```
quast -o allthree \
-t 8 \
--glimmer \
-1 ../rawdata/SRR10015223_1_1mill_mut.fastq.gz \
-2 ../rawdata/SRR10015223_2_1mill_mut.fastq.gz \
-r ../rawdata/GCA_012098105.1_PDT000380306.1_genomic.fna \
-g ../rawdata/GCA_012098105.1_PDT000380306.1_genomic.gff \
PATH/TO/SPADES_ASM \
PATH/TO/FLYE_ASM \
PATH/TO/UNICYCLER_ASM \
-l "spades, flye, unicycler"
```

Answer the same kinds of questions here as you did for the velvet vs spades comparison.
