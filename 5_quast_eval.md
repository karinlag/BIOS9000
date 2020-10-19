# Assembly evaluation

quast -o spades_velvet -k -f --rna-finding -b -1 ../rawdata/1 -2 ../rawdata  



quast -o spades_velvet3 -k -b --rna-finding --glimmer -1 ../rawdata/SRR10015223_1_1mill_mut.fastq -2 ../rawdata/SRR10015223_2_1mill_mut.fastq ../spades/spades_miseq/scaffolds.fasta ../velvet/optimal_k_expvals_paired/contigs.fa


quast -o spades_velvet4 -k -b --rna-finding --glimmer --nanopore ../rawdata/SRR10015224_400k.fastq ../flye/downsampled/assembly.fasta

Notice difference in rRNA genes
Notice difference in number of genes
