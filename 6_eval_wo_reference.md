# Mapping reads to an assembly and visualising the results


Mapping is the process of aligning short reads to a reference genome, with the
goal of figuring where they could come from. We will use `bwa` for mapping.

## Mapping short and long reads

We will first map both the short and the long reads to the Unicycler assembly.
This assembly was created using both of those read sets. We can use this
to see how good our assembly is.


### Indexing the assembly

We will use `bwa` for mapping. This program needs an assembly to
map against, this is our unicycler assembly. `bwa` needs to create an
index of the assembly to make mapping
go faster. For large genomes such as the human genome, this takes a long
time. For the small bacterial genome we work with here this is very fast.

Move (using `cd`) to the folder with your hybrid assembly in it.

Index the fasta file with:

```
bwa index -a bwtsw ASSEMBLY.FASTA
```

Replace `ASSEMBLY.FASTA` with the name of your fasta file. Run `ls` to check
the results, you should see a couple of new files.

For the visualization, you will also need to create an index file for your.
You can do that with this command:

```
samtools faidx assembly.fasta
```



### Mapping paired end reads

Mapping the reads using `bwa mem` yields SAM output. Instead of saving this
output to disk, we will immediately convert it to a sorted (binary) BAM
file by piping into the `samtools`program. 'Sorted' here means that the
alignments of the mapped reads are in the order of the reference sequences,
rather than random. Finally, we will generate an index of the sorted BAM
file for faster searching later on.

First, create a new folder *in the same folder as the `ASSEMBLY.FASTA` file*  
and `cd` into it:

```
mkdir bwa
cd bwa
```
Then do the mapping:

```
bwa mem -t 8 ../ASSEMBLY.FASTA \
../../../rawdata/SRR10015223_1_1mill_mut.fastq.gz \
../../../rawdata/SRR10015223_2_1mill_mut.fastq.gz \
| samtools view -buS - | samtools sort - -o hybrid_ilm.sorted.bam
```

Explanation of some of the parameters:

* `../` means 'look in the folder one level up', i.e. where the fasta file and
its index files are
* `-t n` tells `bwa mem` to use n threads (cpus)
* `-buS` tells `samtools view` that the input is in SAM format (`S`) and to
output uncompressed (`u`) BAM format (`b`).
* the `-` for both `samtools` commands indicate that instead of using a file
as input, the input comes from a pipe (technically, from 'standard in', or
'STDIN').
* ` -o hybrid_ilm.sorted.bam` tells `samtools view` the name of the outputfile

Note: the difference between the mut files and the regular files is just that
the read names in the "normal" look like this: `@SRR10015223.1066259.2`, while
this in the mut files look like this: `@SRR10015223.1066259 2`. The first part
is then the name of the SRA project, then the read name, then whether the read
is number 1 or 2 in the pair.  

Next, generate an index of the BAM file:

```
samtools index hybrid_ilm.sorted.bam
```

If you would like to have a look at the alignments in the BAM file (which is in
binary format), use `samtools view`:

```
samtools view hybrid_ilm.sorted.bam |less -S
```

### Mapping long reads

Mapping long reads works pretty much the same as the short reads, but we
need some extra options.

```
bwa mem -t 8 -x ont2d ../ASSEMBLY.FASTA \
../../../rawdata/SRR10015224.fastq.gz  | \
samtools view -buS - | \
samtools sort - -o hybrid_nano.sorted.bam
```

Then do the indexing again:
```
samtools index hybrid_nano.sorted.bam

```

### Downloading results

We will visualize these results in a genome viewer. But, first we need to
download the results to your computer. Download the entire bwa folder, and
also the assembly fasta file.


### Visualising the assembly in a genome browser

For this part, we will use Integrative Genomics Viewer (IGV), a genome browser
developed by the Broad Institute. This program needs to be installed on your
computer. Go [to this location, and download and install the program](https://software.broadinstitute.org/software/igv/download).

### Loading the data

* Start the IGV program
* Choose `Genomes --> Load Genome from File…` (**NB** not File --> Load from
File...)
* Select the `fasta` file with your assembly (**NB** the same file as you used
for mapping the reads against!), and also the `fai` file.

**Adding the mapped reads**  
Adding tracks to the browser is as simple as uploading a new file:

* Choose `File --> Load from File…`
* Choose the sorted `bam` file of the paired end mapping
* Repeat this for the `bam` file of the nanopore mapping
* You can choose different sequences (contigs/scaffolds) from the drop-down
menu at the top. Start by selecting (one of) the longest scaffold(s)
* Start browsing!
* Zoom in to see the alignments

You can find more information about interpreting what you see
[on this website](http://software.broadinstitute.org/software/igv/PopupMenus#AlignmentTrack).


**Question:**

* Do you see differences between some of the reads relative to the reference?
What are these?
* Is coverage even? Are there gaps in the coverage, or peaks? Where?
* Can you find some regions where the short reads are "white"? What does
that mean?
