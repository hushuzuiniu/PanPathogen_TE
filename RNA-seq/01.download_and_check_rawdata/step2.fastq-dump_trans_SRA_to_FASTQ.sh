#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=fastq-dump
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=6

PRJ=~/pan_pathogen/bacteria/MTB

cd $PRJ
#for sra format -- SRRXXXXX.sra
ls SRR* -dt -r | while read i;
do
        fastq-dump --split-3 $PRJ/${i}/${i}.sra --gzip -O $PRJ/
done

