#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=featureCounts_multi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=39

WD=/share/home/cuijie/PanPathogen_TE/


##-----------1. featureCounts to quantify the readscounts of each sample -------------------------

############################
###                      ###
### for Single-end Reads ###
###                      ###
############################

mkdir -p $WD/05.quantify_TE_loci/featureCounts_out &&

cat $WD/test_sample_name_SE.txt | while read i;
do
     featureCounts -T 10 \
                   -a $WD/00.ref/hg38_TE_anno_custom_v20240110.gtf \
                   -t exon -M \
                   -g transcript_id \
                   -o $WD/05.quantify_TE_loci/featureCounts_out/${i}_readscounts.txt $WD/03.mapping/mapping_out/${i}.Aligned.sortedByCoord.out.bam
done &&


############################
###                      ###
### for Paired-end Reads ###
###                      ###
###########################

mkdir -p $WD/05.quantify_TE_loci/featureCounts_out &&

cat $WD/test_sample_name_PE.txt | while read i;
do
     featureCounts -T 10 \
                   -a $WD/00.ref/hg38_TE_anno_custom_v20240110.gtf \
                   -t exon -M \
		   -p -B -C \
                   -g transcript_id \
                   -o $WD/05.quantify_TE_loci/featureCounts_out/${i}_readscounts.txt $WD/03.mapping/mapping_out/${i}.Aligned.sortedByCoord.out.bam
done &&

	
#----------------------2. combine featureCounts readscounts of all samples for each dataset -----------------------------------------------------

mkdir -p $WD/05.quantify_TE_loci/readscounts_combined &&

cd $WD/05.quantify_TE_loci/featureCounts_out/

ls *readscounts.txt | sed 's/_Pre_SRR.*//g' | sed 's/_Post_SRR.*//g' | sort | uniq | while read dataset;
do
   ls ${dataset}_*readscounts.txt | sed 's/.*_Pre/Pre/g' | sed 's/.*_Post/Post/g' | sed 's/_readscounts.txt//g' | while read id;
   do

          sed '1d' $WD/05.quantify_TE_loci/featureCounts_out/${dataset}_${id}_readscounts.txt | cut -f1-6  > $WD/05.quantify_TE_loci/readscounts_combined/geneID_column_tmp.txt
          sed '1d' $WD/05.quantify_TE_loci/featureCounts_out/${dataset}_${id}_readscounts.txt | awk '{print $7}' | sed 's/\/.*\///g' | sed 's/\.Aligned.out.sam//g' > $WD/05.quantify_TE_loci/readscounts_combined/${dataset}_${id}_readscounts_column_tmp.txt
   done
   paste -d "\t" $WD/05.quantify_TE_loci/readscounts_combined/geneID_column_tmp.txt $WD/05.quantify_TE_loci/readscounts_combined/${dataset}*_readscounts_column_tmp.txt > \
                 $WD/05.quantify_TE_loci/readscounts_combined/${dataset}_readscounts_matrix_TE_Loci.txt
done &&

rm -rf $WD/05.quantify_TE_loci/readscounts_combined/*tmp.txt
