# Velvet practical

Velvet is one of the first assembly programs that was created for short read data.
It is not frequently used today and is generally known for producing fragmented
assemblies. However, one benefit with velvet is that it is a very fast program.
We will thus first experiment a bit with velvet to see how various factors
can affect an assembly.


## Class exercise

In this first exercise, we will:
* log in to the server
* create a file structure that will serve the assembly exercise
* link in the data we will work on
* do a first velvet assembly
* see how we can evaluate an assembly
* figure out what is the best k-mer value for this data set

### Data structure setup

First, log on to the server, and go "home".

```
cd ~
```

or simply type

```
cd
```

Then: create a directory for the assembly exercise.

```
mkdir assembly
cd assembly
mkdir rawdata
```

We will now create a "shortcut" to the data that we will be working on.

```
cd rawdata
ln -s /work/IN-BIOSx/data/assembly/* .
```

When you now do `ls`, you will see nine files there. These files are all from
the same _E. coli_ isolate.

MiSeq files:
* SRR10015223 - Full MiSeq files
* SRR10015223_?_1mill - a 1 million reads subsample of the full MiSeq files
* SRR10015223_?_1mill_mut - we need a slight fix in the header for one program,
this has been fixed for these files

Nanopore:
* SRR10015224 - the full nanopore sequence set
* SRR10015224_400k - a downsampled version of the one above

In addition there is a `perl` program that we will use during the exercises.

### Use of velvet

Velvet consists of two different programs:

* `velveth` - this program creates the index of the reads based on the value of **k**
* `velvetg` - this program builds the graph from the k-mers produced by `velveth`
and outputs the contigs.



### A first velvet assembly - testing out the values of **k**

Velvet is a `de Bruijn` assembler. Thus, it is quite sensitive to the value
of **k**. We will first experiment and figure out what value is best for this
dataset.

Note, the reads in question are 150 bp long. The value of **k** must be shorter
than that, and it has to be an odd number.

First, go to where we will do the assemblies:

```
cd
cd assembly
mkdir velvet
cd velvet
```

Look [in this google spreadsheet](https://docs.google.com/spreadsheets/d/14fuw7fXjBp3A8xCUXmmwkJBb9E-r2B-jqi_Q3Aq6RrQ/edit?usp=sharing),
find a value of *k*, and record your choice in the spreadsheet. Run
`velveth` to build the hash index (see below).

`velveth` options:

* `ASM_NAME`: select your own name for the folder where the results should go
* `value_of_k`: use k-mers of this size
* `-short`: short reads (as opposed to long, Sanger-like reads)
* `-separate`: read1 and read2 are in separate files
* `-fastq`: read type is fastq

Build the index as follows:

```
velveth ASM_NAME VALUE_OF_K \  
-short -separate -fastq \  
../rawdata/SRR10015223_1_1mill.fastq \
../rawdata/SRR10015223_2_1mill.fastq
```
**NOTES**

* Change `ASM_NAME` to a directory name of your choosing
* Change `VALUE_OF_K` to the value you have picked
* The command is split over several lines by adding a space, and a `\`
(backslash) to each line. This trick makes long commands more readable.
Sometimes copy pasting will not work, it is best to write in the command
yourself.

After `velveth` is finished, use `ls` to look in the new folder that has the
name you chose. You should see the following files:

```
Log
Roadmaps
Sequences
```


The '`Log`' file has a useful reminder of what commands you typed to get this
assembly result, for reproducing results later on. '`Sequences`' contains the
sequences we put in, and '`Roadmaps`' contains the index you just created.

Now we will run the assembly with default parameters:

```
velvetg ASM_NAME
```

Velvet will end with a text like this:

`Final graph has ... nodes and n50 of ..., max ..., total ..., using .../... reads`

The number of nodes represents the number of nodes in the graph.
Nodes in the graph get converted to contigs, however, nodes with a length
below that of the k-mer size do not become contigs.
Velvet reports its N50 (as well as everything
else) in 'kmer' space. The conversion to 'basespace' is as simple as adding k-1
to the reported length.

Look again at the folder `ASM_NAME`, you should see the following extra files:

`contigs.fa`  
`Graph`  
`LastGraph`  
`PreGraph`  
`stats.txt`

The important files are:

`contigs.fa` - the assembly itself  
`Graph` - a textual representation of the contig graph  
`stats.txt` - a file containing statistics on each contig

**Questions**

* What k-mer did you use?
* What is the N50 of the assembly? (in basespace, not in k-mer space)
* What is the size of the largest contig? (in basespace, not in k-mer space)
* How many contigs are there in the `contigs.fa` file? Use `grep -c NODE contigs.fa`. Is this the same number as velvet reported?

Log your results in the google spreadsheet.

**We will discuss the results together and determine *the optimal* k-mer for
this dataset.**

FROM NOW ON: keep track of the values that velvet reports on your own!



### Assembly 1: Optimal k-mer assembly

You will now create your own assembly with using this optimal k-mer. Use the
commands shown above. For all assemblies for here on out, use this k-mer size.

### Assembly 2: Estimating and setting `exp_cov`

Much better assemblies are produced if Velvet understands the expected coverage
for unique regions of your genome. This allows it to try and resolve repeats.
The data to determine this is in the `stats.txt` file. The full description of
this file is [in the Velvet Manual](http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf).

First, we will have a look at the file and figure out how it works.

```
less -S stats.txt
```

Here you see wh the cntent of the file is. Have a look at the manual to figure
out what the fields are.

Press `q` to finish looking at the file.

Next, we will fetch the `short1_cov` column and extract that. We will also
sort it, so that we see how high the coverage goes.

```
cat stats.txt |awk '{print $6}' |sort -n
```

You will notice that there is one node (contig) with very high coverage. We will plot this as a histogram. Having one very high number when making a
histogram will just squish things. Thus we need to filter that out. We wil
do that by sending the results to the `grep` command.

```
 cat stats.txt |awk '{print $6}' |grep -v NUMBER |sort -n
```

Put the coverage number you found in instead of NUMBER in the command above.

Then we will plot the histogram.

```
cat stats.txt |awk '{print $6}' |grep -v NUMBER | ../../rawdata/hist.pl
```


**Question:**

* What do you think is the approximate expected k-mer coverage for your assembly?

Stop here, we will discuss things together here.

Now run velvet again, supplying the value for `exp_cov` (k-mer coverage)
corresponding to your answer:

```
velvetg ASM_NAME -exp_cov PEAK_K_MER_COVERAGE
```

**Question:**

* What improvements do you see in the assembly by setting a value for `exp_cov`?

## Group exercise

From now on, you will work in teams of three. Each of you will run one of the
subsequent assemblies, and keep track of your results. The teacher will drop
in on your sessions, and you can ask questions via Mattermost.

In order to keep track of who is doing what, please fill in
[this Google spreadsheet](https://docs.google.com/spreadsheets/d/1Av_T9Tm3b_o_PIcdSmSAUarfzci7IYpXI9nrn0r4daA/edit?usp=sharing).

### Assembly 3: Setting `cov_cutoff`

You can also clean up the graph by removing low-frequency nodes from the
*de Bruijn* graph using the `cov_cutoff` parameter. Low-frequency nodes
can result from sequencing errors, or from parts of the genome with very
little sequencing coverage. Removing them will often result in better
assemblies, but setting the cut-off too high will also result in losing
useful parts of the assembly. Using the histogram from previously, estimate
a good value for `cov_cutoff`.

```
velvetg ASM_NAME -exp_cov YOUR_VALUE -cov_cutoff YOUR_VALUE  
```

Try some different values for `cov_cutoff`, keeping `exp_cov` the same and
record your assembly results.

**Question:**

* What improvements do you see in the assembly by setting a value for
  `cov_cutoff`?


### Assembly 4: Asking velvet to determine the parameters

You can also ask Velvet to predict the values for you:

```
velvetg ASM_NAME -exp_cov auto -cov_cutoff auto
```

You will see the estimated best values reported just above `Final graph` in
the output.


**Questions:**

* What values of *exp_cov* and *cov_cutoff* did Velvet choose?
* Check the output to the screen. Is this assembly better than your best one?

### Assembly 5: Incorporating paired end information

Paired end information contributes additional information to the assembly,
allowing contigs to be scaffolded. We will first re-index your reads telling
Velvet to use paired-end information, by using `-shortPaired` instead
of `-short` for `velveth`. Then, re-run velvetg with the two coverage
related options set to auto.

**!!! IMPORTANT Pick a new name for your assembly !!!**

```
velveth ASM_NAME2 VALUE_OF_K \  
-shortPaired -fastq -separate \  
../rawdata/SRR7896249_1.fastq \
../rawdata/SRR7896249_2.fastq

velvetg ASM_NAME2 -exp_cov auto -cov_cutoff auto  
```
Velvet will towards the bottom of the output report the estimated insert size
for the paired library.

**Questions:**

* How does doing this affect the assembly?
* What does velvet say about the insert size of the paired end library?

### Scaffold and contig metrics

The sequences in the contigs.fa file are actually scaffolds. We will use the
program `assembly-stats` to generate metrics for this, and all following assemblies.

```
assembly-stats -t contigs.fa
```

Run this command on all of the assemblies you have created so far, and report
them in the spreadsheet.


### Looking for repeats

Have a look for contigs which are long and have a much higher coverage than
the average for your genome. One tedious way to do this is to look into
the `contigs.fa` file (with `less`). You will see the name of the contig
('NODE'), it's length and the kmer coverage. However, trying to find
contigs with high coverage this way is not very efficient.  

A faster was is to again use the `stats.txt` file.

Relevant columns are:

1) ID --> sequence ID, same as 'NODE' number in the `contigs.fa` file   
2) lgth --> sequence 'length'  
6) short1_cov --> kmer coverage (column 6)  


Knowing this, we can use the `awk` command to select lines for contigs at
least 1kb, with k-mer coverage greater than 100:

```
awk '($2>=1000 && $6>=100)' stats.txt
```

`awk` is an amazing program for tabular data. In this case, we ask it to check
that column 2 ($2, the length) is at least 1000 and column 6 ($6, coverage) at
least 100. If this is the case, awk will print the entire line. See
[http://bit.ly/QjbWr7](http://bit.ly/QjbWr7) for more information on awk.

Find the contigs with the highest coverage in the `contigs.fa` file. Perform a
BLAST search using NCBI. Look at the `Graphics` results to see what is in the
region of the results that you got.

**Question:**

* What is it?
* Is this surprising? Why, or why not?
