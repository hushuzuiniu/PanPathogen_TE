#!/bin/bash


#---------------- ruTE (recurrent up-TEs intergenic) vs RE ------------------------

# v1: overlap > 0bp
bedtools intersect -a ruTE_n29646_0based.bed -b RE_n339693_0based.bed -wo > ruTE_vs_RE_gt0bp_wo_results.txt
# v2: overlap > 50% TE length
bedtools intersect -a ruTE_n29646_0based.bed -b RE_n339693_0based.bed -f 0.5  -wo > ruTE_vs_RE_gt0.5TE_wo_results.txt
# V3: overlap > 100% TE length
bedtools intersect -a ruTE_n29646_0based.bed -b RE_n339693_0based.bed -f 1 -wo > ruTE_vs_RE_gt1TE_wo_results.txt



#-------------- other intergenic TE vs RE ----------------------------------------
# v1: overlap > 0bp
bedtools intersect -a otherTE_n1881375_0based.bed -b RE_n339693_0based.bed -wo > otherTE_vs_RE_gt0bp_wo_results.txt
# v2: overlap > 50% TE length
bedtools intersect -a otherTE_n1881375_0based.bed -b RE_n339693_0based.bed -f 0.5 -wo > otherTE_vs_RE_gt0.5TE_wo_results.txt
# v3: overlap > 100% TE length
bedtools intersect -a otherTE_n1881375_0based.bed -b RE_n339693_0based.bed -f 1 -wo > otherTE_vs_RE_gt1TE_wo_results.txt


# v1: overlap > 0bp
cat ruTE_vs_RE_gt0bp_wo_results.txt | awk -F '\t' '{print $4}' | sort | uniq | wc -l
# 4415
cat  otherTE_vs_RE_gt0bp_wo_results.txt | awk -F '\t' '{print $10}' | sort | uniq | wc -l
# 221218
# v2: overlap > 50% TE length
cat ruTE_vs_RE_gt0.5TE_wo_results.txt | awk -F '\t' '{print $4}' | sort | uniq | wc -l
# 3688
cat  otherTE_vs_RE_gt0.5TE_wo_results.txt | awk -F '\t' '{print $10}' | sort | uniq | wc -l
# 172379
# v3: overlap > 100% TE length
cat ruTE_vs_RE_gt1TE_wo_results.txt | awk -F '\t' '{print $4}' | sort | uniq | wc -l
# 3133
cat  otherTE_vs_RE_gt1TE_wo_results.txt | awk -F '\t' '{print $10}' | sort | uniq | wc -l
# 138814


# v1: overlap > 0bp
cat ruTE_vs_RE_gt0bp_wo_results.txt |awk -F '\t' '$11 == 1 {print $4}'  | sort | uniq | wc -l
# 3099
cat  otherTE_vs_RE_gt0bp_wo_results.txt | awk -F '\t' '$17 == 1 {print $10}' | sort | uniq | wc -l
# 155171
# v2: overlap > 50% TE length
cat ruTE_vs_RE_gt0.5TE_wo_results.txt | awk -F '\t' '$11 == 1 {print $4}' | sort | uniq | wc -l
# 2785
cat  otherTE_vs_RE_gt0.5TE_wo_results.txt |awk -F '\t' '$17 == 1 {print $10}'  | sort | uniq | wc -l
# 130160
# v3: overlap > 100% TE length
cat ruTE_vs_RE_gt1TE_wo_results.txt | awk -F '\t' '$11 == 1 {print $4}' | sort | uniq | wc -l
# 2481
cat  otherTE_vs_RE_gt1TE_wo_results.txt | awk -F '\t' '$17 == 1 {print $10}' | sort | uniq | wc -l
# 110590

