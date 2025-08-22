#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=QC
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8

# activate conda env "RNASeq"
# conda activate RNASeq

#-------------------- 1.fastp ------------------------------------
############################
###                      ###
### for Single-end Reads ###
###                      ###
############################
WD=~/PanPathogen_TE

mkdir -p $WD/02.fastp_out &&
cat $WD/test_sample_name_SE.txt | while read i;
do
fastp -i $WD/01.raw_data/${i}.fq.gz \
      -o $WD/02.fastp_out/${i}.clean.fq.gz \
      -q 20 \
      -l 50 \
      --thread 16 \
      -h $WD/02.fastp_out/${i}.html \
      -j $WD/02.fastp_out/${i}.json \

done &&

############################
###                      ###
### for Paired-end Reads ###
###                      ###
###########################

mkdir -p $WD/02.fastp_out &&
cat $WD/test_sample_name_PE.txt | while read i;
do
 fastp -i $WD/01.raw_data/${i}_1.fq.gz \
       -I $WD/01.raw_data/${i}_2.fq.gz \
       -o $WD/02.fastp_out/${i}_1.clean.fq.gz \
       -O $WD/02.fastp_out/${i}_2.clean.fq.gz \
       -q 20 \
       -l 50 \
       --thread 16 \
       -h $WD/02.fastp_out/${i}.html \
       -j $WD/02.fastp_out/${i}.json
done &&


#-------------------- 2.multiqc -----------------------------------
mkdir -p $WD/02.multiqc_out &&
multiqc $WD/02.fastp_out -o $WD/02.multiqc_out/

