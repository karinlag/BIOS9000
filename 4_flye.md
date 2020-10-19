# Long read assembly

Flye is a program that does long read assembly of both Nanopore and Pacbio
data. Most notable it is able to skip some steps that are very expensive in OLC
assemblers, that is error correction and the all-against-all comparison steps.
It does that by creating `disjointigs`, which are contigs which have been
created by taking reads at random and trying to join them together by using
overlaps. This means that it is quite fast. It also has options for fishing
out haplotypes and plasmids, and can also be used on metagenomes.  


## Using Flye

Flye can use both PacBio and Nanopore data, both raw and error corrected. It
can also take in a high quality assembly to assist.

**Input data:**
* --pacbio-raw - PacBio raw reads
* --pacbio-corr - PacBio corrected reads
* --pacbio-hifi - PacBio HiFi reads
* --nano-raw - ONT raw reads
* --nano-corr - ONT corrected reads
* --subassemblies - high-quality contigs input

**Run options**
* -t #n - number of threads
* -i #n - number of rounds of polishing


### Setting up the assembly

First, make sure you are in the assembly directory. Then, create a directory
named `flye`. Go into that directory.

### Running Flye

* Start `screen`
* Then, run the following command:

```
flye --nano-raw ../rawdata/SRR10015224_400 -o flye -t 8 &> flye.log
```
