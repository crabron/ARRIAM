#!/usr/bin/env bash
set -e

map=''
workdirectory=''
input=''
databasefna=''
databasetax=''
databasetree=''
trim_sw=''
trim_minlen=''


# main module. Starts with arguments, check the map, pathways and files. Call other modules one by one

# Parameters module
while [ -n "$1" ]; do
case $1 in
    -h) echo ""
        echo "Some help"
        echo ""
        echo "QIIME, vsearch, Trimmomatic and fastq-join REQUIRED"
        echo ""
        echo "script parameters:"
        echo ""
        echo "-i - input directory with raw .fastq.gz files. You can catch them from all sequence pool"
        echo "-m - map file in QIIME format. Second column should be <SampleID>.fna, third - current filename prefix (without _L001..)"
        echo "-o - output folder for processing"
        echo "-p - parameter file. It must content full path for reference database and trimming parameters:"
        echo "          export databasefna='/home/alexey/tax_n_refs/silva_132_97_16S.fna'"
        echo "          export databasetax='/home/alexey/tax_n_refs/taxonomy_7_levels.txt'"
        echo "          export databasetree='/home/alexey/tax_n_refs/97_otus.tre'"
        echo "          export trim_sw='4:12'"
        echo "          export trim_minlen='180'"
        exit 0;;
    -i) input=$2
     shift;;
    -m) map=$2
     shift;;
    -o) workdirectory=$2
     shift;;
    -p) param=$2
     shift;;
    *) echo "$1 is not a valid parameter. Please, read -h annotation"
      exit 1
esac
shift
done

#export other variables
source $param

#...and check our values

#check a map as a file

map=`pwd $map`/$map
if [[ ${map} == "" ]]; then
    echo "Please, enter a map name, use -m parameter"
    exit 1
elif [[ ! -f "${map}" ]]; then
    echo "Error: wrong map file - ${map}";
    exit 2
fi

#check an input
if [[ ${input} == "" ]]; then
    echo "Please, enter an input directory, use -i parameter"
    exit 1
elif [[ ! -d "${input}" ]]; then
    echo "Error: wrong input directory - ${input}";
    exit 2
fi


#check a workdirectory
if [[ ${workdirectory} == "" ]]; then
    echo "Please, specify the work directory"
    exit 1
else
    mkdir ${workdirectory}
fi



#copy files to our directory and rename them

echo -n 'Copy raw files...'
mkdir ${workdirectory}/raw_data

lenMap=`wc -l $map | cut -f 1 -d " "`
for (( i=2; i <= ${lenMap}; i++ )); do
    line=`sed -n -e ${i}p $map`;
    sampleID=`echo ${line} | awk '{print $1}'`
    prefix=`echo ${line} | awk '{print $3}'`
    echo $sampleID >> $workdirectory/samples.txt


    if [[ -f "${input}/${prefix}_L001_R1_001.fastq.gz" ]] && [[ -f "${input}/${prefix}_L001_R2_001.fastq.gz" ]]; then
       cp ${input}/${prefix}_L001_R1_001.fastq.gz ${workdirectory}/raw_data/${sampleID}_L001_R1_001.fastq.gz
       cp ${input}/${prefix}_L001_R2_001.fastq.gz ${workdirectory}/raw_data/${sampleID}_L001_R2_001.fastq.gz
       else
       echo
       echo "Error in ${prefix} pair in ${input}: please, check a table or files in folder"
       exit 2
    fi
    done
echo 'done.'


cd $workdirectory




#trimming and joining

echo 'Trimming...'

mkdir tmp_trimmed
mkdir trimmed

for sample in $(cat samples.txt); do
echo ${sample};
dirname="tmp_trimmed/${sample}";
forward="raw_data/${sample}_L001_R1_001.fastq.gz";
reverse="raw_data/${sample}_L001_R2_001.fastq.gz";
mkdir ${dirname};
trimmomatic PE -phred33 2>>trimmed/trimming.log \
                ${forward} \
                ${reverse} \
                ${dirname}/trimed_paired_forward.fastq.gz \
                ${dirname}/trimed_upaired_forward.fastq.gz \
                ${dirname}/trimed_paired_reverse.fastq.gz \
                ${dirname}/trimed_unpaired_reverse.fastq.gz \
                SLIDINGWINDOW:${trim_sw} MINLEN:${trim_minlen};
fastq-join 1>>trimmed/trimming.log \
                ${dirname}/trimed_paired_forward.fastq.gz \
                ${dirname}/trimed_paired_reverse.fastq.gz \
                -o ${dirname}/trimmed_%.fastq.gz;
cp ${dirname}/trimmed_join.fastq.gz \
    trimmed/${sample}.fastq.gz;
done

rm tmp_trimmed -r
echo 'Trimming done'
echo


#chimera removing

echo 'Chimera removing...'
mkdir dechimered
mkdir tmp_dechimered

for sample in $(cat samples.txt); do
echo ${sample};
dirname="tmp_dechimered/${sample}";
mkdir ${dirname}
cp trimmed/${sample}.fastq.gz ${dirname}/${sample}.fastq.gz
gzip -d ${dirname}/${sample}.fastq.gz;
convert_fastaqual_fastq.py \
                -c fastq_to_fastaqual\
                -f ${dirname}/${sample}.fastq\
                -o ${dirname}/fastaqual;
vsearch 2>>dechimered/chimera.log --uchime_ref ${dirname}/fastaqual/${sample}.fna\
        --nonchimeras ${dirname}/dechim_${sample}.fna\
        --db ${databasefna} --threads 20;
mv ${dirname}/dechim_${sample}.fna\
    dechimered/${sample}.fna;

done

rm -r tmp_dechimered
echo 'Chimera removing done'
echo


#make OTU picking

echo -n 'adding QIIME labels...'
add_qiime_labels.py -i dechimered -m ${map} -c Filename -o qiime_labeled
cat qiime_labeled/*.fna >> all.fna
echo 'done'

echo -n 'Closed reference OTU picking...'
pick_closed_reference_otus.py -i all.fna -o otus -r ${databasefna} -t ${databasetax} -a -O 8
echo 'done'

echo -n 'Filtering, sorting and converting of BIOM table...'
rm all.fna
rm -r qiime_labeled
filter_taxa_from_otu_table.py -i otus/otu_table.biom -o otus/f_otu_table.biom -n D_3__Chloroplast,D_4__Mitochondria
sort_otu_table.py -i otus/f_otu_table.biom -o otus/sf_otu_table.biom
biom convert -i otus/sf_otu_table.biom -o otus/sf_otu_table.txt --to-tsv --header-key taxonomy

lenTable=`wc -l otus/sf_otu_table.txt | cut -f 1 -d " "`
for (( i=3; i <= ${lenTable}; i++ )); do
    line=`sed -n -e ${i}p otus/sf_otu_table.txt`;
    echo ${line} | awk '{print $1}' >>otus/otu_list.txt;
done
echo 'done'

echo -n 'Filtering tree...'
filter_tree.py -i ${databasetree} -t otus/otu_list.txt -o otus/filtered.tre
echo 'done'


echo -n 'Barplots...'
summarize_taxa_through_plots.py -o taxa_summary -i otus/sf_otu_table.biom -m ${map}
echo 'done'

echo -n 'Alpha diversity...'
alpha_rarefaction.py -i otus/sf_otu_table.biom -o alpha_div -t otus/filtered.tre -m ${map}
echo 'done'

echo -n 'Beta-diversity...'
beta_diversity_through_plots.py -i otus/sf_otu_table.biom -o beta_div -t otus/filtered.tre -m ${map}
echo 'done'

exit 0


