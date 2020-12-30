#!/bin/bash

while getopts u:d:p:f: option
do
case "${option}" #switch statement:start
in
u) USER=${OPTARG};;
d) DATE=${OPTARG};;
p) PRODUCT=${OPTARG};;
f) FORMAT=${OPTARG};;
esac #switch statement: end
done



conda activate qiime2-2019.1
sed -i '1s/^/#OTU_ID/' otu_table.txt 
biom convert -i otu_table.txt -o otu_table.biom --to-hdf5
qiime tools import \
--input-path rep_seq.fasta \
--output-path rep_seq.qza \
--type 'FeatureData[Sequence]'
qiime fragment-insertion sepp --i-representative-sequences rep_seq.qza  --o-tree data/insertion-tree.qza  --o-placements data/insertion-placements.qza --p-threads 20
qiime tools import \
--input-path otu_table.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path feature-table.qza
qiime fragment-insertion filter-features \
--i-table feature-table.qza \
--i-tree insertion-tree.qza \
--o-filtered-table filtered_table.qza \
--o-removed-table removed_table.qza
unzip -p filtered_table.qza */data/* > filtered_table.biom
biom convert -i  filtered_table.biom -o filtered_table.txt --to-tsv
sed -i '1d' filtered_table.txt
unzip -p insertion-tree.qza */data/* > tree.nwk

