#!/bin/bash

#references sont deja indexées
#dezipage des fastq

gzip -d fastq/EchF*.gz
mkdir output
mkdir output/sam/

all_read=(fastq/*R1.fastq fastq/*R2.fastq)

./soft/bowtie2 --fast --end-to-end -x databases/all_genome.fasta -p 6 -1 ${all_read[0]} -2 ${all_read[1]} -S output/sam/EchF.sam

#2)

mkdir output/bam/

samtools view -Sb output/sam/*.sam > output/bam/aln.bam

samtools sort output/bam/aln.bam > output/bam/aln_sorted.bam

samtools index output/bam/aln_sorted.bam

samtools idxstats output/bam/aln_sorted.bam > output/bam/stat.txt
mkdir output/draw/
grep ">" databases/all_genome.fasta|cut -f 2 -d ">" >output/draw/association.tsv

#3)

./soft/megahit -1 ${all_read[0]} -2 ${all_read[1]} --k-list 21 --mem-flag 0 -o output/assembly

#4)

mkdir output/gene_prediction/

./soft/prodigal -i output/assembly/final.contigs.fa -d output/gene_prediction/gene.fasta

#5)

sed "s:>:*\n>:g" output/gene_prediction/gene.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > output/gene_prediction/genes_full.fasta

#6)

mkdir output/blast/

./soft/blastn -db databases/resfinder.fna -query output/gene_prediction/genes_full.fasta -evalue 0.001 -perc_identity 80 -qcov_hsp_perc 80 -out output/blast/blast_res.txt






