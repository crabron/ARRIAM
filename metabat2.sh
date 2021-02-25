#!/bin/bash

# consensus.fasta - fasta после сборки flye / после сборки и полировки racon + medaka
bwa index consensus.fasta
bwa mem  -t 84 consensus.fasta all.fq | samtools sort -o all.bam
jgi_summarize_bam_contig_depths --outputDepth true_depth.txt --referenceFasta consensus.fasta  all.bam
metabat2 -i consensus.fasta -a true_depth.txt -o metabat2_bins_3
# Выдача metabat2:
# MetaBAT 2 (2.15 (Bioconda)) using minContig 2500, minCV 1.0, minCVSum 1.0, maxP 95%, minS 60, maxEdges 200 and minClsSize 200000. with random seed=1603977749

maCMD -p 'Nanopore' -x mc_out/consensus.json -i ../all.fq -t 64 -o mc_out/consensus.sam

for filename in *.faa; do
    hmmscan --noali --notextw --acc -E 0.000000000000000001 --cpu 50 -o  hm_${filename}.txt hmm/glyco.hmm $filename
done

hmmscan --noali --notextw --acc -E 0.000000000000000001 --cpu 50 -o  hm_02.txt glyco.hmm pr_02.faa

for filename in *.faa; do
    echo $filename
done

awk '{s++}END{print s/4}' file.fastq


 cat ../c1_ann.tsv | cut -f2 | grep  pr_01.gff | awk '{print $1}' | taxonkit --data-dir ~/storage/temp/  -j 50 lineage | taxonkit --data-dir ~/storage/temp/ -j 50 reformat


 cat ../c1_ann.tsv | cut -f2 |grep -f /dev/stdin pr_01.gff | awk '{print $1}' | 


 cat ../c1_ann.tsv | cut -f2 |grep -f /dev/stdin pr_01.gff | awk '{print $1}' | grep -f /dev/stdin ka_01.tsv | cut -f3 | taxonkit --data-dir ~/storage/temp/  -j 50 lineage | taxonkit --data-dir ~/storage/temp/ -j 50 reformat > kid_01.txt


checkm ssu_finder consensus.fasta bins/ out/ -x fa -t 20

checkm qa ./out/lineage.ms ./out -o 2 -f res --tab_table -t 20


wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml
conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2019.10-py36-linux-conda.yml

qiime picrust2 full-pipeline  --p-threads 30 --p-max-nsti --i-table feature-table.qza --i-seq rep_seq.qza

qiime feature-table summarize \
   --i-table picrust/pathway_abundance.qza \
   --o-visualization picrust/pathway_abundance.qzv

   unzip -p filtered_table.qza > pathway_abundance.biom
   biom convert -i  pathway_abundance.biom -o pathway_abundance.txt --to-tsv


#    try to use MiSeq/Pilon

 bbduk.sh in1=qtrim=r trimq=8 46_S10_L001_R1_001.fastq.gz 46_S10_L001_R2_001.fastq.gz
 
 
 kraken2 --db $KRAKEN2_DB_PATH --threads 50 gbk_concensus/Ng110_gene.faa > Ng110_kraken.tsv