#!/bin/bash
set -e 

mkdir -p 01.fastqc
mkdir -p 02.trim_galore
mkdir -p 03.mapping/index
mkdir -p 04.filter_after_mapping
mkdir -p 05.peaks
mkdir -p 05.bigwig
# replace Bacteria_1.txt to other dataset files to run other datasets
cat  Bacteria_1.txt | while read id; do
    fastqc -t 20 ./00.raw_data/${id}_1.fq.gz ./00.raw_data/${id}_2.fq.gz -o 01.fastqc/
done

cd /data2t_2/hushu/ATAC-seq/01.multiqc
multiqc ../01.fastqc/


THREADS=8  
cat  Bacteria_1.txt | while read id; do
    trim_galore -j 4 -q 25 --phred33 --fastqc --paired --length 25 -e 0.1 --stringency 4 \
        --cores 8 \
        -o ./02.trim_galore/ \
        ./00.raw_data/${id}_1.fq.gz ./00.raw_data/${id}_2.fq.gz
done


# ln -s $(realpath 00.ref/hg38.p13.fa) 03.mapping/index/
# cd 03.mapping/index/
# bwa index -a bwtsw hg38.p13.fa  

cat  Bacteria_1.txt | while read id; do
    bwa mem -M -t 20 03.mapping/index/hg38.p13.fa \
        ./02.trim_galore/${id}_1_val_1.fq.gz \
        ./02.trim_galore/${id}_2_val_2.fq.gz \
        -o ./03.mapping/${id}.sam
done


cat  Bacteria_1.txt | while read id; do
    samtools sort -@ 40 \
        ./03.mapping/${id}.sam \
        -o ./03.mapping/${id}.bam
    

    samtools index -@ 40 ./03.mapping/${id}.bam
    # rm ./03.mapping/${id}.sam
done



cat  Bacteria_1.txt | while read id; do
    samtools view -@ 10 -bh -q 20 -F 260 \
        03.mapping/${id}.bam \
        > ./04.filter_after_mapping/${id}.filtered.bam
done


cd ./04.filter_after_mapping/
cat ../Bacteria_1.txt | while read id; do
    if [ -f "${id}.filtered.bam" ]; then
        samtools addreplacerg -r "@RG\tID:${id}\tSM:${id}\tPL:ILLUMINA\tLB:lib1" -o "${id}.filtered.tmp.bam" "${id}.filtered.bam"
    else
        echo " error ${id}.filtered.bam"
    fi
done

cat ../Bacteria_1.txt | while read id; do
    picard MarkDuplicates \
        REMOVE_DUPLICATES=true \
        I=./${id}.filtered.tmp.bam \
        O=./${id}.filtered.rmdup.bam \
        M=./${id}.picard.MarkDuplicates.out.txt
done
cd /data2t_2/hushu/ATAC-seq/

cd ./04.filter_after_mapping/
cat ../Bacteria_1.txt | while read id; do
    bedtools intersect -v -abam ./${id}.filtered.rmdup.bam \
        -b ../00.ref/hg38_unified_blacklist.bed \
        > ./${id}.filtered.rmdup.rmbl.bam
done
cd /data2t_2/hushu/ATAC-seq/

cat  ./Bacteria_1.txt | parallel -j 4 'macs2 callpeak \
    -t ./04.filter_after_mapping/{}.filtered.rmdup.rmbl.bam \
    -n {} \
    --outdir ./05.peaks \
    --nomodel \
    --extsize 200 \
    --SPMR \
    -f BAMPE -B -g hs \
    --keep-dup all --shift -100'

cat ./Bacteria_1.txt | while read id; do
    samtools index -M ./04.filter_after_mapping/${id}.filtered.rmdup.rmbl.bam
done

cat  ./Bacteria_1.txt | while read id; do
    bamCoverage \
        -b ./04.filter_after_mapping/${id}.filtered.rmdup.rmbl.bam \
        -o ./05.bigwig/${id}.bw \
        --normalizeUsing RPKM \
        -p 40
done
