#!/bin/bash

# MTB 
for i in GSE227851 GSE194017 GSE126614 GSE164287 GSE173560 GSE132283 GSE143731
do
   prefetch --option-file ${i}_SRR_Acc_List.txt --force all 1> ${i}_prefetch_download_srr.log
done

