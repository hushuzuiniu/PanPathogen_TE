#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=rename_sample_name_20240620
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10


cd ~/Pathogen_TE/01.raw_data/pathogen/
#--- modify fastq as fq ---------------
ls *fastq.gz | cut -d"." -f 1 | while read id;
do
        mv ${id}.fastq.gz ${id}.fq.gz
done &&

#--- change sample name ------------
while read a b;
do
    mv ${a}.fq.gz ${b}.fq.gz
    mv ${a}_1.fq.gz ${b}_1.fq.gz
    mv ${a}_2.fq.gz ${b}_2.fq.gz
done < pathogen_SampleName_old_new_list.txt

