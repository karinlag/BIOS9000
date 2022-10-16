# Polishing an assembly

As you noticed, there were many mismatches in the `flye` assembly. We will
now try to error correct this assembly with the short reads. We will
do this with a program that is called `pilon`. What we do here, is to map
the short reads to the flye assembly with `bwa`. The output is then given to
the `pilon` program, which uses it to correct small errors. Note, this program
is not able to do something with misassemblies.

[The manual of pilon can be found here]( https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage).

### Mapping to the assembly

The first thing we have to do is to map the reads to the `flye` assembly.

First, go to the `flye` assembly directory. Then, create an index for your
assembly.

```
bwa index -a bwtsw assembly.fasta
```

Then create a new directory for doing the mapping, and go into it.

```
mkdir bwa
cd bwa
```

Then we will do the actual mapping.

```
bwa mem -t 8 ../assembly.fasta \
../../../rawdata/SRR10015223_1_1mill_mut.fastq.gz \
../../../rawdata/SRR10015223_2_1mill_mut.fastq.gz \
| samtools view -buS - | samtools sort - -o flye.bam
```

After that is done, we need to index the resulting bam file:

```
samtools index flye.bam
```

### Pilon correction

Now that we have the bam file, we can do the correction.

Go one level up. First we need to set the memory for pilon, like we did earlier:

```
export _JAVA_OPTIONS='-Xmx12G'
```

Then we start `pilon`:

```
java -Xmx8G -jar $EBROOTPILON/pilon.jar \
--genome assembly.fasta \
--bam bwa/flye.bam \
--changes --vcfqe &> flye_pilon.log
```

### QUAST comparison

We can now do a new quast comparison to see what happens. This is the
abbreviated quast command:

```
quast -o flye_flyepilon \
-t 5 \
-r ../../rawdata/GCA_12098104.0_PDT380305.0_genomic.fna \
-g ../../rawdata/GCA_12098104.0_PDT380305.0_genomic.gff \
assembly.fasta pilon.fasta \
-l "Original, Pilon"
```

The results can then be transferred and compared again as before.
