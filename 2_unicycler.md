# Hybrid assemblies

Unicycler is a hybrid assembler, which can also be used just on long or short
reads. It uses both long (PacBio and Nanopore) and short read data to create
assemblies. It first performs a SPAdes assembly with the short reads, and then
uses the long reads to sort out the resulting scaffolds (this is a
simplification). It can sort out plasmids, and also circularize genomes. Note
that Unicycler can take a long time to run.

## Documentation

For more about [the Unicycler program see here](https://github.com/rrwick/Unicycler).


## How to use Unicycer

This program does not have too many options, however, it does have some options
for tuning the output that you get. This can give you more information on what
is happening.  

**Input:**
* -1 SHORT1 - FASTQ file of first short reads in each pair (required)
* -2 SHORT2 - FASTQ file of second short reads in each pair (required)
* -s UNPAIRED - FASTQ file of unpaired short reads (optional)
* -l LONG - FASTQ or FASTA file of long reads (optional)

**Options for running:**

By adjusting the value to `--verbosity` from 0 - 3 you will adjust how much
info it will give while it is running.

In addition, you need to regulate how many threads (cpus) the program uses, by
adjusting the `-t` option.


### Setting up the assembly

Ensure you are in the `assembly` directory. Create a subdirectory named
`unicycler`. Go into that directory.

We will save the output from the command using `> unicycler.log` in a file to
be able to follow progress. `&>` makes sure any error-messages are written to
the same file.

### Running the assembly

* First, start `screen`
* Then load the software  

We need to adjust how much memory the Pilon program can use. We do that with
this command:

```
export _JAVA_OPTIONS='-Xmx12G'
```

Then, we run the assembly using:

```
unicycler \
-1 ../rawdata/SRR10015223_1_1mill.fastq.gz \
-2 ../rawdata/SRR10015223_2_1mill.fastq.gz \
-l ../rawdata/SRR10015224.fastq.gz \
-o hybrid --verbosity 2 -t 8 &> hybrid.log
```

### Output files

The most important output files are
* assembly.gfa - the assembly in an assembly graph file format,
* assembly.fasta - the assembly
* unicycler.log  - this is the log file
