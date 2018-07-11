#$-N Rock_weathering_generateOTU
#$-cwd
#$-l vf=5G
#$-m beas
#$-M roey.angel@bc.cas.cz
#$-p -1
#$-pe parallel 8

#!/bin/bash
#title          :Rock_weathering_analysis.sh
#description    :Analysis of rock_weathering MiSeq sequences
#author         :Roey Angel
#date           :20160315
#version        :1.0    
#usage          :qsub Rock_weathering_analysis.sh ${DATASET} ${PLATFORM} ${DATAFILETYPE} ${AMPLENGTH}
#notes          :       
#bash_version   :4.3.30(1)-release
#============================================================================

. /etc/profile

# Set parameters
#DATASET=$1
DATASET="Rock_weathering"
#PLATFORM=$2
PLATFORM="MiSeq"
#DATAFILETYPE=$3
DATAFILETYPE="FASTAQ"
# AMPLENGTH=$4
# AMPLENGTH="357"

declare -a SAMPLES
SAMPLES=("GrfCr1S46" "GrfCr2S47" "GrfCr3S48" "ShivCr2S49" "ShivCr3S50" "ShivCr4S51" "GrfCr2S70" "GrfCr3S71" "GrfCr4S72" "TamSlF1S82" "TamSlF2S83" "TamSlF3S84" "TamSlF4S85" "ShivSlF1S86" "ShivSlF2S87" "ShivSlF3S88" "ShivSlF4S89" "GrfSlF1S90" "GrfSlF2S91" "GrfSlF3S92" "GrfSlF4S93")
COUNTMIN="2"
OTUdef="97"
OTU_id=`echo -e "(100-${OTUdef})/100" | bc`
OTU_radius=`echo -e "100-${OTUdef}" | bc`
CLASSIFY="Y"
CHIMCHECK="Y"

# Set paths
MOTHUR=/home/newapps/mothur/1.36.1/
USEARCH=/home/newapps/usearch/8.1.1861/
DBS=~/DBs/ # needs to contains: uchime - gold.fa, Silva.seed_v119.tgz
SCRIPTS=~/Scripts/MiSeq/
RUNDIR=./Rock_weathering/rock_weathering_front_samples/
# RUNDIR=`pwd`"/"
WORKDIR=/tmp/${DATASET}_work/
RESULTSDIR=/tmp/${DATASET}_results/

# Create dirs
if [[ -d ${WORKDIR} ]]
then
 rm -r ${WORKDIR}
fi
if [[ -d ${RESULTSDIR} ]]
then
 rm -r ${RESULTSDIR}
fi

mkdir ${WORKDIR}
mkdir ${RESULTSDIR}

# Process raw files according to platform and data type
echo -e "Analysis log of the ${DATASET} dataset.\nPlatform=${PLATFORM}, filetype=${DATAFILETYPE}" > ${RUNDIR}${DATASET}_Log.txt
echo -e "===========================================\n" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "Using: " >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch >> ${RUNDIR}${DATASET}_Log.txt
echo -e "Cite: Edgar, 2013; Edgar & Flyvbjerg, 2015" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

for fastqfile in *.fastq
do
    gzip ${fastqfile}
done

cp *.gz ${WORKDIR}
cd ${WORKDIR}
gunzip ${WORKDIR}*.gz 

echo -e "The following files were found in the folder and will be analysed:" >> ${RUNDIR}${DATASET}_Log.txt
FASTQFILES=`ls *.fastq`
echo -e $FASTQFILES >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt
  
# Merge fastq files according to run
echo -e "Merge fastq pairs by samples into one file" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "------------------------------------------" >> ${RUNDIR}${DATASET}_Log.txt

echo -e "Prepending sample names to file names:" >> ${RUNDIR}${DATASET}_Log.txt
FASTQFILESR1=($(ls *R1*.fastq)) # this is bad programming but what the hell
FASTQFILESR2=($(ls *R2*.fastq))
for (( i=0; i<${#SAMPLES[@]}; i++ )); do
    mv ${FASTQFILESR1[i]} ${SAMPLES[i]}${FASTQFILESR1[i]}
    mv ${FASTQFILESR2[i]} ${SAMPLES[i]}${FASTQFILESR2[i]}
done
echo -e "The following samples names have been prepended: ${SAMPLES[@]}\n" >> ${RUNDIR}${DATASET}_Log.txt

echo -e "Merge fastq pairs:" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -fastq_mergepairs *_R1_*.fastq -relabel @  -report mergereport.txt -fastqout ${DATASET}_merged.fq" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -fastq_mergepairs *_R1_*.fastq -relabel @ -report mergereport.txt -fastqout ${DATASET}_merged.fq
cat mergereport.txt >> ${RUNDIR}${DATASET}_Log.txt

echo -e "\nCharacters distribuion:" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -fastq_chars  ${DATASET}_merged.fq -log chars.log"
${USEARCH}usearch -fastq_chars  ${DATASET}_merged.fq -log chars.log
cat chars.log >> ${RUNDIR}${DATASET}_Log.txt

echo -e "\nQuality statistics:" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -fastq_eestats2 ${DATASET}_merged.fq -output eestats2.txt -length_cutoffs 200,300,10" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -fastq_eestats2 ${DATASET}_merged.fq -output eestats2.txt -length_cutoffs 200,300,10
cat eestats2.txt >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# Quality filter
echo -e "Quality filter fastq pairs" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "--------------------------" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -fastq_filter ${DATASET}_merged.fq -fastq_maxee 1.0 -fastqout ${DATASET}_filtered.fq" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -fastq_filter ${DATASET}_merged.fq -fastq_minlen 400 -fastq_maxee 1.0 -relabel Filt -fastaout ${DATASET}_filtered.fa
echo -e "Quality statistics" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -fastq_eestats2 ${DATASET}_filtered.fa -output eestats2.txt -length_cutoffs 200,300,10
echo -e "usearch -fastq_eestats2 ${DATASET}_filtered.fa -output eestats2.txt -length_cutoffs 200,300,10" >> ${RUNDIR}${DATASET}_Log.txt
cat eestats2.txt >> ${RUNDIR}${DATASET}_Log.txt

${MOTHUR}mothur "#summary.seqs(fasta=${DATASET}_filtered.fa, processors=10)" > tmp.log
grep -m 1 'mothur >' tmp.log >> ${RUNDIR}${DATASET}_Log.txt 
grep -m 1 -A1000 'Start' tmp.log | head -n -3  >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# Dereplicate sequences for OTU table
echo -e "Dereplicate sequences for OTU table" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-----------------------------------" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -derep_fulllength ${DATASET}_filtered.fa -minuniquesize ${COUNTMIN} -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -derep_fulllength ${DATASET}_filtered.fa -relabel Uniq -minuniquesize ${COUNTMIN} -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log
if [[ ${COUNTMIN} -gt 1 ]]; then
  echo -e "Sequences smaller than ${COUNTMIN} will not be allowed to form OTUs" >> ${RUNDIR}${DATASET}_Log.txt
fi
cat tmp.log >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# Sort sequences by size
echo -e "Sort sequences by size" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "----------------------" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -sortbysize ${DATASET}_unique.fa -fastaout ${DATASET}_seqs_sorted.fa --log tmp.log" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -sortbysize ${DATASET}_unique.fa -fastaout ${DATASET}_seqs_sorted.fa --log tmp.log
cat tmp.log >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# Find OTU representatives
echo -e "Find OTU representatives" >> ${RUNDIR}${DATASET}_Log.txt
# if chimeras are being filtered:
if [[ $CHIMCHECK == "Y" ]]
then
 echo -e "-----------------------" >> ${RUNDIR}${DATASET}_Log.txt
 echo -e "Chimeras are being checked" >> ${RUNDIR}${DATASET}_Log.txt
 echo -e " usearch -cluster_otus ${DATASET}_seqs_sorted.fa -otu_radius_pct $OTU_radius -sizeout -otus ${DATASET}_otus.fa -relabel OTU -maxrejects 0 -maxaccepts 0 --log tmp.log" >> ${RUNDIR}${DATASET}_Log.txt
 ${USEARCH}usearch -cluster_otus ${DATASET}_seqs_sorted.fa -otu_radius_pct $OTU_radius -sizeout -otus ${DATASET}_otus.fa -relabel OTU -maxrejects 0 -maxaccepts 0 --log tmp.log
 cat tmp.log >> ${RUNDIR}${DATASET}_Log.txt
 echo -e "\n" >> ${RUNDIR}${DATASET}_Log.txt
 # make a ref-based chimera check
 echo -e "Make a ref-based chimera check:" >> ${RUNDIR}${DATASET}_Log.txt
 echo e- "usearch -uchime_ref ${DATASET}_otus.fa -db ${DBS}gold.fa -nonchimeras nochimeras.fa -strand plus --log tmp.log" >> ${RUNDIR}${DATASET}_Log.txt
 ${USEARCH}usearch -uchime_ref ${DATASET}_otus.fa -db ${DBS}gold.fa -nonchimeras ${DATASET}_otuReps.fa -strand plus --log tmp.log
 cat tmp.log >> ${RUNDIR}${DATASET}_Log.txt
 echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt
 
else
 echo -e "Chimeras are not being checked" >> ${RUNDIR}${DATASET}_Log.txt
 echo -e "------------------------------" >> ${RUNDIR}${DATASET}_Log.txt
 echo -e "usearch -cluster_smallmem ${DATASET}_seqs_sorted.fa -id $OTU_id -sizeout -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -centroids ${DATASET}_otus.fa -usersort -maxrejects 0 -maxaccepts 0 --log tmp.log" >> ${RUNDIR}${DATASET}_Log.txt
 ${USEARCH}usearch -cluster_smallmem ${DATASET}_seqs_sorted.fa -id $OTU_id -sizeout -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -centroids ${DATASET}_otuReps.fa -usersort -maxrejects 0 -maxaccepts 0 --log tmp.log
 cat tmp.log >> ${RUNDIR}${DATASET}_Log.txt
 echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt
fi
cp ${DATASET}_otuReps.fa ${RESULTSDIR}${DATASET}_otuReps.fa

# Determine abundance of OTUs in samples  
echo -e "Determine abundance of OTUs in samples" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "--------------------------------------" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -usearch_global ${DATASET}_merged.fq -db ${RESULTSDIR}${DATASET}_otuReps.fa -strand plus -id $OTU_id -maxrejects 0 -maxaccepts 0 -maxhits 1 -uc ${DATASET}_readmap.uc -matched ${RESULTSDIR}${DATASET}_matched_seqs.fa -mothur_shared_out ${DATASET}_otuTab.shared -otutabout ${DATASET}_otuTab.txt" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -usearch_global ${DATASET}_merged.fq -db ${RESULTSDIR}${DATASET}_otuReps.fa -strand plus -id $OTU_id -maxrejects 0 -maxaccepts 0 -maxhits 1 -uc ${DATASET}_readmap.uc -matched ${RESULTSDIR}${DATASET}_matched_seqs.fa -mothur_shared_out ${DATASET}_otuTab.shared -otutabout ${DATASET}_otuTab.txt
cp ${DATASET}_otuTab.* ${RESULTSDIR}
echo -e `grep -v -c "*" ${DATASET}_readmap.uc`" merged read pairs map to clusters" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# How many reads did not match any OTU?
echo -e "How many sequences did not match any OTU?" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-----------------------------------------" >> ${RUNDIR}${DATASET}_Log.txt
grep -c "^N" ${DATASET}_readmap.uc >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# Classify sequences 
echo -e "Classify sequences" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "------------------" >> ${RUNDIR}${DATASET}_Log.txt
if [[ $CLASSIFY == "Y" ]]
then
 echo -e "Classifying OTUs:" >> ${RUNDIR}${DATASET}_Log.txt
 cp ${DBS}Silva.seed_v119.tgz ./
 tar -xzf Silva.seed_v119.tgz
 cp ${RESULTSDIR}/${DATASET}_otuReps.fa ${DATASET}_otuReps.fa
 ${MOTHUR}mothur "#classify.seqs(fasta=${DATASET}_otuReps.fa,template=silva.seed_v119.align, taxonomy=silva.seed_v119.tax, processors=10)" > tmp.log
 grep -m 1 'mothur >' tmp.log >> ${RUNDIR}${DATASET}_Log.txt 
 grep -m 1 -A1000 'Output File Names:' tmp.log | head -n -3  >> ${RUNDIR}${DATASET}_Log.txt   
 echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt
 cp ${DATASET}_otuReps.nr_v119.wang.taxonomy ${RESULTSDIR}${DATASET}_silva.v119.taxonomy
fi

# Creat singleton table
echo -e "Creating an OTU table of rejected singletons" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "----------------------------------------------" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "This is used for evaluation only. It should worry us if the singletons show some pattern." >> ${RUNDIR}${DATASET}_Log.txt
# Dereplicate sequences for singleton table
echo -e "Dereplicate sequences for singleton table:" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "{USEARCH}usearch -derep_fulllength ${RESULTSDIR}${DATASET}_matched_seqs.fa -minuniquesize 1 -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -derep_fulllength ${RESULTSDIR}${DATASET}_matched_seqs.fa -minuniquesize 1 -fastaout ${DATASET}_unique.fa -sizeout --log tmp.log

cat tmp.log >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

cp ${DATASET}_unique.fa ${RESULTSDIR}${DATASET}_unique.fa

# Determine abundance of singletons in samples  
echo -e "Determine abundance of singletons in samples:" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "--------------------------------------------" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "usearch -usearch_global ${RESULTSDIR}${DATASET}_matched_seqs.fa -db ${RESULTSDIR}${DATASET}_unique.fa -strand plus -id 1.0 -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -maxrejects 0 -maxhits 1 -uc ${DATASET}_unique.uc -mothur_shared_out ${DATASET}_uniqueOtuTab.shared -otutabout ${DATASET}_uniqueOtuTab.txt" >> ${RUNDIR}${DATASET}_Log.txt
${USEARCH}usearch -usearch_global ${RESULTSDIR}${DATASET}_matched_seqs.fa -db ${RESULTSDIR}${DATASET}_unique.fa -strand plus -id 1.0 -gapopen 10.0I/10.0E -gapext 1.0I/1.0E -maxrejects 0 -maxhits 1 -uc ${DATASET}_unique.uc -mothur_shared_out ${DATASET}_uniqueOtuTab.shared -otutabout ${DATASET}_uniqueOtuTab.txt

cp ${DATASET}_uniqueOtuTab.* ${RESULTSDIR}
echo -e `grep -v -c "*" ${DATASET}_unique.uc`" merged unique read pairs map to clusters" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt
# echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt

# Classify sequences 
echo -e "Classify sequences:" >> ${RUNDIR}${DATASET}_Log.txt
echo -e "------------------" >> ${RUNDIR}${DATASET}_Log.txt
if [[ $CLASSIFY == "Y" ]]
then
 echo -e "Classifying OTUs:" >> ${RUNDIR}${DATASET}_Log.txt
 cp ${RESULTSDIR}${DATASET}_unique.fa ${DATASET}_unique.fa
 ${MOTHUR}mothur "#classify.seqs(fasta=${DATASET}_unique.fa,template=silva.nr_v119.align, taxonomy=silva.nr_v119.tax, processors=10)" > tmp.log
 grep -m 1 'mothur >' tmp.log >> ${RUNDIR}${DATASET}_Log.txt 
 grep -m 1 -A1000 'Output File Names:' tmp.log | head -n -3  >> ${RUNDIR}${DATASET}_Log.txt   
 echo -e "-- \n" >> ${RUNDIR}${DATASET}_Log.txt
 cp ${DATASET}_unique.nr_v119.wang.taxonomy ${RESULTSDIR}${DATASET}_unique_silva.nrv119.taxonomy
fi
cp 

# Return results 
cp ${RUNDIR}${DATASET}_Log.txt ${RESULTSDIR}
tar -zcvf ${RUNDIR}${DATASET}.OTU.tar.gz ${RESULTSDIR}/*
#rm *.logfile silva.nr* tmp.log
cd ..
#tar -zcvf$ ${RUNDIR}${DATASET}.OTU.tar.gz ${WORKDIR}/*

# Cleanup
rm -r ${WORKDIR}
rm -r ${RESULTSDIR}
