# Assembly using SPAdes

SPAdes was written as an assembly program for bacterial genomes, both for
data from single cells, as well as from whole-genome amplified samples.
SPAdes also works well, sometimes even best, when given high-coverage datasets.
It is generally considered to do better on assemblies than velvet.

Before assembly, SPADES will error-correct the reads.

## Using SPADES

Spades can be used for many different kinds of data. In our case, we will
in this round use it for single isolate data. Thus, we will want to use
with this option:

* The `--isolate` flag is used for high coverage single isolate data.

In addition, we want to ensure that our assemblies are as good as possible:

* The `--cov-cutoff auto` flag is used to tell SPAdes to leave out low coverage
contigs.

For each read file, a flag is used to indicate what kind of dataset it is.


**Input data:**
* -1 and -2: files with forward and reverse reads
* -pe-1/2 #n : in case of multiple datasets with forward and reverse reads,
this denotes the first/second file of reads in the #nth read set

In addition, there are options for Nanopore and Pacbio data, SPAdes can also
be used for hybrid assemblies.      

**Other parameters:**
* `-t` number of threads (CPUs) to use for calculations
* `-m` maximum memory usage in Gb
* `-k` k-mers to use (this gives room for experimenting!)
* `-o` name of the output folder


### Setting up the assembly

First, go to your `assembly` folder, and create a new folder called
`spades`. Then `cd` into it. We will save the output from the command using
`>spades_miseq.log` in a file to be able to follow progress. `&>` makes sure any
error-messages are written to the same file.

**NOTE** the assembly will take several hours, so use the `screen` command!
`screen` is a program which will start a terminal kind of "within" a terminal.
You can start one, and exit and enter it at will. Whatever you run in it will
continue running in it even if you close the terminal or lose your internet
connection.

See [this quick guide on screen usage](0_tech.md)

### Running the assembly

* First, start `screen`
* Then load the software
* Then, enter the following command:

```
spades.py -t 4 --isolate --cov-cutoff auto -m 10 \
-1 ../rawdata/SRR10015223_1_1mill.fastq \
-2 ../rawdata/SRR10015223_1_1mill.fastq \
-o spades_miseq &> spades_miseq.log
```
