# 创建输出目录
mkdir -p TE_loci_bed
mkdir -p TE_sequences

# 读取subfamily列表，为每个subfamily提取loci
while read subfamily; do
    echo "Processing $subfamily..."
    
    # 从GTF文件提取该subfamily的所有loci
    grep -w "gene_id \"$subfamily\"" /data2t_2/hushu/00.ref/hg38_TE_anno_custom_v20240110.gtf| \
    awk -v OFS='\t' '{
        # 计算中心点
        center = int(($4 + $5) / 2)
        # 输出: chr, start, end, name, score, strand
        print $1, $4, $5, $12, ".", $7
    }' | \
    # 去掉引号
    sed 's/"//g' | sed 's/;//g' > TE_loci_bed/${subfamily}_loci.bed
    
    echo "Found $(wc -l < TE_loci_bed/${subfamily}_loci.bed) loci for $subfamily"
    
done < Recurrent_TE_loci_and_subfamily_uniq_n55.txt

GENOME_FA="/data2t_2/hushu/00.ref/hg38.fa" 

# 为每个subfamily提取序列
while read subfamily; do
    echo "Extracting sequences for $subfamily..."
    # 使用bedtools getfasta提取序列，保持链方向
    bedtools getfasta -fi $GENOME_FA \
                      -bed TE_loci_bed/${subfamily}_loci.bed \
                      -s \
                      -name \
                      -fo TE_sequences/${subfamily}_loci.fa
    
    echo "Extracted $(grep -c ">" TE_sequences/${subfamily}_loci.fa) sequences for $subfamily"
done < Recurrent_TE_loci_and_subfamily_uniq_n55.txt
