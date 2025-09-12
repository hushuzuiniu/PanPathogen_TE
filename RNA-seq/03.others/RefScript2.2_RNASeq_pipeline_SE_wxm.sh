#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=RNASeq_for_Gene&TE_Pipeline_20240103
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40


############################
###                      ###
### for Singled-end Reads ###
###                      ###
############################

#set the parameters:
PRJ=~/IFN_TE/GSE159462_K562_IFNa_bulkRNA-Seq
sjdbOverhang=74
Thread=40


#--------  step1:  fastqc to qc; multiqc to merge all qc results of all samples ---------------------
:<<EOF
# qc
mkdir -p $PRJ/02.qc_before_filter &&
cd $PRJ/01.raw_data
ls *gz | cut -d"." -f 1 | while read i;
do
  fastqc -t $Thread $PRJ/01.raw_data/${i}.fq.gz -o $PRJ/02.qc_before_filter/
done &&

cd $PRJ/02.qc_before_filter &&
multiqc * &&

# -------  step2: trim-golare to trim low quantity of reads and adapter, and to obtain the paired mapped reads, then perform fastqc again ----------------------
mkdir -p $PRJ/03.clean_data &&
cd $PRJ/01.raw_data/
for i in Mock_rep1 Mock_rep2 Mock_rep3 IFNa_rep1 IFNa_rep2 IFNa_rep3
do
     trim_galore -j 4 -q 25 --phred33 --fastqc --length 25 -e 0.1 --stringency 4 \
                 -o $PRJ/03.clean_data/ \
                 $PRJ/01.raw_data/${i}.fq.gz  \
                 >> $PRJ/03.clean_data/filter_trim-galore.singled.out.log
done &&

cd $PRJ/03.clean_data &&
multiqc * &&

#----------- step3: STAR to map reads to reference genome  --------------------------------------

##-------------3.1 generate indexs for STAR ----------------------
mkdir -p $PRJ/04.mapping_quantify/index &&
STAR  --runThreadN $Thread --runMode genomeGenerate --genomeDir $PRJ/04.mapping_quantify/index  --genomeFastaFiles $PRJ/00.ref/hg38.p13.fa --sjdbGTFfile $PRJ/00.ref/hg38.p13.gene.anno.gtf --sjdbOverhang $sjdbOverhang  &&

##-------------3.2 reads mapping by STAR--------------------------
mkdir -p $PRJ/04.mapping_quantify/mapping &&
for i in Mock_rep1 Mock_rep2 Mock_rep3 IFNa_rep1 IFNa_rep2 IFNa_rep3
do
    STAR --runThreadN $Thread  \
         --genomeDir $PRJ/04.mapping_quantify/index \
         --readFilesIn $PRJ/03.clean_data/${i}_trimmed.fq.gz \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix $PRJ/04.mapping_quantify/mapping/${i}. \
         --quantMode TranscriptomeSAM GeneCounts
done &&

###############
## for Gene  ##
###############
##-----------step4: featureCounts to quantify the raw readscounts --------------------------
mkdir -p $PRJ/04.mapping_quantify/quantify &&
for i in Mock_rep1 Mock_rep2 Mock_rep3 IFNa_rep1 IFNa_rep2 IFNa_rep3
do
     featureCounts -T $Thread \
                   -a $PRJ/00.ref/hg38.p13.gene.anno.gtf \
                   -t exon \
                   -g gene_name \
                   -o $PRJ/04.mapping_quantify/quantify/${i}_readscounts.txt $PRJ/04.mapping_quantify/mapping/${i}.Aligned.out.sam
done &&

##-----------step5: combine readscounts of all samples as one matrix-----------------------
mkdir -p $PRJ/04.mapping_quantify/readscounts_matrix &&
for i in Mock_rep1 Mock_rep2 Mock_rep3 IFNa_rep1 IFNa_rep2 IFNa_rep3
do
      sed '1d' $PRJ/04.mapping_quantify/quantify/${i}_readscounts.txt | awk '{print $7}' | sed 's/\/.*\///g' | sed 's/\.Aligned.out.sam//g'  \
             > $PRJ/04.mapping_quantify/readscounts_matrix/${i}_readscounts_column_temp.txt &&
      sed '1d' $PRJ/04.mapping_quantify/quantify/${i}_readscounts.txt | cut -f1-6 \
             > $PRJ/04.mapping_quantify/readscounts_matrix/gene_id_column_temp.txt &&
      paste $PRJ/04.mapping_quantify/readscounts_matrix/gene_id_column_temp.txt $PRJ/04.mapping_quantify/readscounts_matrix/*_readscounts_column_temp.txt \
             > $PRJ/04.mapping_quantify/readscounts_matrix/K562_IFNa_vs_Mock_readscounts_matrix.txt
done &&
      rm -rf $PRJ/04.mapping_quantify/readscounts_matrix/*temp.txt

EOF

###############
## for TE    ##
###############
##-----------step4: featureCounts to quantify the raw readscounts --------------------------
mkdir -p $PRJ/04.mapping_quantify/quantify_TE &&
for i in Mock_rep1 Mock_rep2 Mock_rep3 IFNa_rep1 IFNa_rep2 IFNa_rep3
do
     featureCounts -T $Thread \
                   -a $PRJ/00.ref/hg38_te_anno_custom_v20230418.gtf \
                   -t exon \
                   -g transcript_id \
                   -o $PRJ/04.mapping_quantify/quantify_TE/${i}_readscounts.txt $PRJ/04.mapping_quantify/mapping/${i}.Aligned.out.sam
done &&

##-----------step5: combine readscounts of all samples as one matrix-----------------------
mkdir -p $PRJ/04.mapping_quantify/readscounts_matrix_TE &&
for i in Mock_rep1 Mock_rep2 Mock_rep3 IFNa_rep1 IFNa_rep2 IFNa_rep3
do
      sed '1d' $PRJ/04.mapping_quantify/quantify_TE/${i}_readscounts.txt | awk '{print $7}' | sed 's/\/.*\///g' | sed 's/\.Aligned.out.sam//g'  \
             > $PRJ/04.mapping_quantify/readscounts_matrix_TE/${i}_readscounts_column_temp.txt &&
      sed '1d' $PRJ/04.mapping_quantify/quantify_TE/${i}_readscounts.txt | cut -f1-6 \
             > $PRJ/04.mapping_quantify/readscounts_matrix_TE/gene_id_column_temp.txt &&
      paste $PRJ/04.mapping_quantify/readscounts_matrix_TE/gene_id_column_temp.txt $PRJ/04.mapping_quantify/readscounts_matrix_TE/*_readscounts_column_temp.txt \
             > $PRJ/04.mapping_quantify/readscounts_matrix_TE/K562_IFNa_vs_Mock_readscounts_matrix_TE.txt
done &&
      rm -rf $PRJ/04.mapping_quantify/readscounts_matrix_TE/*temp.txt





