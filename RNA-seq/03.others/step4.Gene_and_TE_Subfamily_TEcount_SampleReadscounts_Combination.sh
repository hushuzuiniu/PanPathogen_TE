#/bin/bash!

PRJ=~/wxm/04.STAR_default_TEcount_uniq_Virus_n1846

cd $PRJ/03.readscounts_all/

ls *cntTable | sed 's/_Pre_SRR.*//g' | sed 's/_Post_SRR.*//g' | sort | uniq | while read dataset;
do 
   ls ${dataset}_*cntTable | sed 's/.*_Pre/Pre/g' | sed 's/.*_Post/Post/g' | sed 's/.cntTable//g' | while read id;
   do

   #---------- step1.combine TEcount tables of all samples for each dataset---------------------------------------------------

	  cut -f1  $PRJ/03.readscounts_all/${dataset}_${id}.cntTable | sed 's/"//g' | sed 's/gene\/TE/GeneID/g' > $PRJ/04.combine_readscounts_all/geneID_column_tmp.txt
	  awk '{print $2}' $PRJ/03.readscounts_all/${dataset}_${id}.cntTable | sed 's/\/.*\///g' | sed 's/Aligned.out.sam//g' > $PRJ/04.combine_readscounts_all/${dataset}_${id}_readscounts_column_tmp.txt 
   done
   paste -d "\t" $PRJ/04.combine_readscounts_all/geneID_column_tmp.txt $PRJ/04.combine_readscounts_all/${dataset}*_readscounts_column_tmp.txt > \
		  $PRJ/04.combine_readscounts_all/${dataset}_Post_vs_Pre_readscounts_matrix_GeneTE.txt

   #------------step2.get Gene and TE readscounts respectively for each dataset -----------------------------------------------
   
   grep -e "GeneID" -e "ENSG" $PRJ/04.combine_readscounts_all/${dataset}_Post_vs_Pre_readscounts_matrix_GeneTE.txt > $PRJ/04.combine_readscounts_all/${dataset}_Post_vs_Pre_readscounts_matrix_Gene.txt
   grep -v "ENSG" $PRJ/04.combine_readscounts_all/${dataset}_Post_vs_Pre_readscounts_matrix_GeneTE.txt >  $PRJ/04.combine_readscounts_all/${dataset}_Post_vs_Pre_readscounts_matrix_TE.txt
done 

rm -rf $PRJ/04.combine_readscounts_all/*tmp.txt


	

