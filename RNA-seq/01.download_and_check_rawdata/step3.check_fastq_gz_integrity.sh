#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=check_integrity
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=6

PRJ=~/pan_pathogen/bacteria/MTB

cd $PRJ
#for sra format -- SRRXXXXX.fastq.gz
ls SRR*fastq.gz -dt -r | while read i;
do
     gzip -t ${i} 2>> check_out.txt
done &&

if [ -s "check_out.txt" ]; then
     sed 's/gzip: //g' check_out.txt | sed 's/:.*//g' | sed '/^[[:space:]]*$/d' | sort | uniq > \
		integrity_failed_samples.txt
else 
     echo "check_out.txt file is 0, exit!"
     exit 0
fi


