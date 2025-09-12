#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=Rickettsia_STAR_mapping
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=39

WD=/share/home/cuijie/PanPathogen_TE/
export WD=/share/home/cuijie/PanPathogen_TE/
#----------- step3: STAR to map reads to reference genome  --------------------------------------

##-------------3.1 generate indexs for STAR ----------------------
# mkdir -p $WD/03.mapping/index &&
# STAR  --runThreadN 20 --runMode genomeGenerate --genomeDir $WD/03.mapping/index  --genomeFastaFiles $WD/00.ref/hg38.p13.fa --sjdbGTFfile $WD/00.ref/hg38.p13.gene.anno.gtf &&

##-------------3.2 reads mapping by STAR--------------------------
############################
###                      ###
### for Single-end Reads ###
###                      ###
############################
mkdir -p $WD/03.mapping/mapping_out

# Create a function to run STAR for a given sample
run_star_se() {
    i=$1
    STAR --runThreadN 8 \
         --genomeDir "$WD/03.mapping/index" \
         --readFilesIn "$WD/02.fastp_out/${i}.clean.fq.gz" \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix "$WD/03.mapping/mapping_out/${i}." \
         --quantMode GeneCounts \
         --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
         --outSAMtype BAM SortedByCoordinate
}

export -f run_star_se

cat $WD/test_sample_name_SE.txt | parallel -j 4 run_star_se {}



############################
###                      ###
### for Paired-end Reads ###
###                      ###
###########################

mkdir -p $WD/03.mapping/mapping_out

# Create a function to run STAR for a given sample
run_star_pe() {
    i=$1
    STAR --runThreadN 8 -p \
         --genomeDir "$WD/03.mapping/index" \
         --readFilesIn "$WD/02.fastp_out/${i}_1.clean.fq.gz" "$WD/02.fastp_out/${i}_2.clean.fq.gz" \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix "$WD/03.mapping/mapping_out/${i}." \
         --quantMode GeneCounts \
         --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
         --outSAMtype BAM SortedByCoordinate
}

export -f run_star_pe

cat $WD/test_sample_name_PE.txt | parallel -j 4 run_star_pe {}

