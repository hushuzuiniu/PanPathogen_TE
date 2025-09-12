#!/bin/bash
#SBATCH --output=job.%j.out
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --job-name=tecount_multi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=39

WD=/share/home/cuijie/PanPathogen_TE/


#------------------------1. quantify Gene&TE readscounts per sample by using TEcount (TEtranscrits software)------------------------
mkdir -p $WD/04.quantify_Gene_and_TE_subfamily/TEcount_out

cd $WD/03.mapping/mapping_out

# Export necessary variables
export WD

# Generate a list of BAM files and run TEcount in parallel
ls *.bam | cut -d "." -f1 | parallel -j 39 '/share/home/cuijie/software/TEtranscripts-master/bin/TEcount --sortByPos -b $WD/03.mapping/mapping_out/{}.Aligned.sortedByCoord.out.bam --GTF $WD/00.ref/hg38.p13.gene.anno.gtf --TE $WD/00.ref/hg38_TE_anno_custom_v20240110.gtf --format BAM --mode multi --project {} --outdir $WD/04.quantify_Gene_and_TE_subfamily/TEcount_out/'



#----------------------2. combine  Gene&TE readscounts of all samples for each dataset -----------------------------------------------------
mkdir -p $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined &&

cd $WD/04.quantify_Gene_and_TE_subfamily/TEcount_out/

ls *cntTable | sed 's/_Pre_SRR.*//g' | sed 's/_Post_SRR.*//g' | sort | uniq | while read dataset;
do
   ls ${dataset}_*cntTable | sed 's/.*_Pre/Pre/g' | sed 's/.*_Post/Post/g' | sed 's/.cntTable//g' | while read id;
   do
#---------- 2.1.combine TEcount tables of all samples for each dataset---------------------------------------------------

          cut -f1  $WD/04.quantify_Gene_and_TE_subfamily/TEcount_out/${dataset}_${id}.cntTable | sed 's/"//g' | sed 's/gene\/TE/GeneID/g' > $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/geneID_column_tmp.txt
          awk '{print $2}' $WD/04.quantify_Gene_and_TE_subfamily/TEcount_out/${dataset}_${id}.cntTable | sed 's/\/.*\///g' | sed 's/Aligned.out.sam//g' > $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}_${id}_readscounts_column_tmp.txt
   done
   paste -d "\t" $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/geneID_column_tmp.txt $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}*_readscounts_column_tmp.txt > \
                 $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}_readscounts_matrix_GeneTE.txt

#----------2.2.get Gene and TE readscounts respectively for each dataset -----------------------------------------------

   grep -e "GeneID" -e "ENSG" $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}_readscounts_matrix_GeneTE.txt > $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}_readscounts_matrix_Gene.txt
   grep -v "ENSG" $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}_readscounts_matrix_GeneTE.txt >  $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/${dataset}_readscounts_matrix_TE.txt
done

rm -rf $WD/04.quantify_Gene_and_TE_subfamily/readscounts_combined/*tmp.txt

