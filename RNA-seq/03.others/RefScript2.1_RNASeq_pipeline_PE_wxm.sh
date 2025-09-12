#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=RNASeq_for_Gene&TE_Pipeline_2024614
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20


############################
###                      ###
### for Paired-end Reads ###
###                      ###
############################

#set the parameters:
PRJ=~/IFN_TE/GSE165428_THP1_6h_RNASeq/
sjdbOverhang=149
Thread=40

:<<EOF
#--------  step1:  fastqc to qc; multiqc to merge all qc results of all samples ---------------------

# qc
mkdir -p $PRJ/02.qc_before_filter &&
cd $PRJ/01.raw_data
ls *gz | cut -d"." -f 1 | while read i;
do
  fastqc -t 40 $PRJ/01.raw_data/${i}.fq.gz -o $PRJ/02.qc_before_filter/
done &&

multiqc $PRJ/02.qc_before_filter/* &&



# -------  step2: trim-golare to trim low quantity of reads and adapter, and to obtain the paired mapped reads, then perform fastqc again ----------------------
mkdir -p $PRJ/03.clean_data &&
cd $PRJ/01.raw_data/
ls *gz | cut -d"." -f 1 | sed 's/_[1-2]//g' | uniq |  while read i;
do
     trim_galore -j 4 -q 25 --phred33 --fastqc --paired --length 25 -e 0.1 --stringency 4 \
                 -o $PRJ/03.clean_data/ \
                 $PRJ/01.raw_data/${i}_1.fq.gz  $PRJ/01.raw_data/${i}_2.fq.gz \
                 >> $PRJ/03.clean_data/filter_trim-galore.paired.out.log
done &&

#----------- step3: STAR to map reads to reference genome  --------------------------------------

##-------------3.1 generate indexs for STAR ----------------------
mkdir -p $PRJ/04.mapping_quantify/index &&
STAR  --runThreadN $Thread --runMode genomeGenerate --genomeDir $PRJ/04.mapping_quantify/index  --genomeFastaFiles $PRJ/00.ref/hg38.p13.fa --sjdbGTFfile $PRJ/00.ref/hg38.p13.gene.anno.gtf --sjdbOverhang $sjdbOverhang  &&

EOF

##-------------3.2 reads mapping by STAR--------------------------
mkdir -p $PRJ/04.mapping_quantify/mapping &&
cd  $PRJ/01.raw_data/ &&
ls *gz | cut -d"." -f 1 | sed 's/_[1-2]//g' | uniq  |  while read i;
do
    STAR --runThreadN $Thread -p \
         --genomeDir $PRJ/04.mapping_quantify/index \
         --readFilesIn $PRJ/03.clean_data/${i}_1_val_1.fq.gz  $PRJ/03.clean_data/${i}_2_val_2.fq.gz \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix $PRJ/04.mapping_quantify/mapping/${i}. \
         --quantMode TranscriptomeSAM GeneCounts
done &&


###############
## for Gene  ##
###############
##-----------step4: featureCounts to quantify the raw readscounts --------------------------
mkdir -p $PRJ/04.mapping_quantify/quantify_Gene &&
cd  $PRJ/01.raw_data/ &&
ls *gz | cut -d"." -f 1 | sed 's/_[1-2]//g' | uniq  |  while read i;
do
     featureCounts -T $Thread \
                   -a $PRJ/00.ref/hg38.p13.gene.anno.gtf \
                   -t exon \
                   -g gene_name \
                   -p -B -C \
                   -o $PRJ/04.mapping_quantify/quantify_Gene/${i}_readscounts.txt $PRJ/04.mapping_quantify/mapping/${i}.Aligned.out.sam
done &&

##-----------step5: combine readscounts of all samples as one matrix-----------------------
mkdir -p $PRJ/04.mapping_quantify/readscounts_matrix_Gene &&
cd  $PRJ/01.raw_data/ &&
ls *gz | cut -d"." -f 1 | sed 's/_[1-2]//g' | uniq  |  while read i;
do
      sed '1d' $PRJ/04.mapping_quantify/quantify_Gene/${i}_readscounts.txt | awk '{print $7}' | sed 's/\/.*\///g' | sed 's/\.Aligned.out.sam//g'  \
             > $PRJ/04.mapping_quantify/readscounts_matrix_Gene/${i}_readscounts_column_temp.txt &&
      sed '1d' $PRJ/04.mapping_quantify/quantify_Gene/${i}_readscounts.txt | cut -f1-6 \
             > $PRJ/04.mapping_quantify/readscounts_matrix_Gene/gene_id_column_temp.txt &&
      paste $PRJ/04.mapping_quantify/readscounts_matrix_Gene/gene_id_column_temp.txt $PRJ/04.mapping_quantify/readscounts_matrix_Gene/*_readscounts_column_temp.txt \
             > $PRJ/04.mapping_quantify/readscounts_matrix_Gene/THP1_6h_IFNa_vs_Mock_readscounts_matrix_Gene.txt
done &&
      rm -rf $PRJ/04.mapping_quantify/readscounts_matrix_Gene/*temp.txt

###############
## for TE    ##
###############

##-----------step4: featureCounts to quantify the raw readscounts --------------------------
mkdir -p $PRJ/04.mapping_quantify/quantify_TE &&
cd  $PRJ/01.raw_data/ &&
ls *gz | cut -d"." -f 1 | sed 's/_[1-2]//g' | uniq  |  while read i;
do
     featureCounts -T $Thread \
                   -a $PRJ/00.ref/hg38_te_anno_custom_v20230418.gtf \
                   -t exon \
                   -g transcript_id \
                   -p -B -C \
                   -o $PRJ/04.mapping_quantify/quantify_TE/${i}_readscounts.txt $PRJ/04.mapping_quantify/mapping/${i}.Aligned.out.sam
done &&

##-----------step5: combine readscounts of all samples as one matrix-----------------------
mkdir -p $PRJ/04.mapping_quantify/readscounts_matrix_TE &&
cd  $PRJ/01.raw_data/ &&
ls *gz | cut -d"." -f 1 | sed 's/_[1-2]//g' | uniq  |  while read i;
do
      sed '1d' $PRJ/04.mapping_quantify/quantify_TE/${i}_readscounts.txt | awk '{print $7}' | sed 's/\/.*\///g' | sed 's/\.Aligned.out.sam//g'  \
             > $PRJ/04.mapping_quantify/readscounts_matrix_TE/${i}_readscounts_column_temp.txt &&
      sed '1d' $PRJ/04.mapping_quantify/quantify_TE/${i}_readscounts.txt | cut -f1-6 \
             > $PRJ/04.mapping_quantify/readscounts_matrix_TE/gene_id_column_temp.txt &&
      paste $PRJ/04.mapping_quantify/readscounts_matrix_TE/gene_id_column_temp.txt $PRJ/04.mapping_quantify/readscounts_matrix_TE/*_readscounts_column_temp.txt \
             > $PRJ/04.mapping_quantify/readscounts_matrix_TE/THP1_6h_IFNa_vs_Mock_readscounts_matrix_TE.txt
done &&
      rm -rf $PRJ/04.mapping_quantify/readscounts_matrix_TE/*temp.txt

