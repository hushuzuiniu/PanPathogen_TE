#!/bin/bash


WD=~/PanPathogen_TE/

mkdir -p $WD/04.quantify/readscounts_combined &&

cd $WD/04.quantify/TEcount_out/

ls *cntTable | sed 's/_Pre_SRR.*//g' | sed 's/_Post_SRR.*//g' | sort | uniq | while read dataset;
do 
   ls ${dataset}_*cntTable | sed 's/.*_Pre/Pre/g' | sed 's/.*_Post/Post/g' | sed 's/.cntTable//g' | while read id;
   do

   #---------- step1.combine TEcount tables of all samples for each dataset---------------------------------------------------

	  cut -f1  $WD/04.quantify/TEcount_out/${dataset}_${id}.cntTable | sed 's/"//g' | sed 's/gene\/TE/GeneID/g' > $WD/04.quantify/readscounts_combined/geneID_column_tmp.txt
	  awk '{print $2}' $WD/04.quantify/TEcount_out/${dataset}_${id}.cntTable | sed 's/\/.*\///g' | sed 's/Aligned.out.sam//g' > $WD/04.quantify/readscounts_combined/${dataset}_${id}_readscounts_column_tmp.txt 
   done
   paste -d "\t" $WD/04.quantify/readscounts_combined/geneID_column_tmp.txt $WD/04.quantify/readscounts_combined/${dataset}*_readscounts_column_tmp.txt > \
		 $WD/04.quantify/readscounts_combined/${dataset}_readscounts_matrix_GeneTE.txt

   #------------step2.get Gene and TE readscounts respectively for each dataset -----------------------------------------------
   
   grep -e "GeneID" -e "ENSG" $WD/04.quantify/readscounts_combined/${dataset}_readscounts_matrix_GeneTE.txt > $WD/04.quantify/readscounts_combined/${dataset}_readscounts_matrix_Gene.txt
   grep -v "ENSG" $WD/04.quantify/readscounts_combined/${dataset}_readscounts_matrix_GeneTE.txt >  $WD/04.quantify/readscounts_combined/${dataset}_readscounts_matrix_TE.txt
done 

rm -rf $WD/04.quantify/readscounts_combined/*tmp.txt


	

