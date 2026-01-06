#!/bin/bash
# 为每个TE subfamily单独运行HOMER findMotifsGenome.pl
# Target: 每个subfamily的所有loci
# Background: 统一背景 (所有intergenic TE排除34个recurrent subfamily)

WORKDIR="/data2t_2/pathogen_TE_2025_New/08.TF_prediction/TF_prediction_add_bg_v3"
HOMER_PATH="/data2t_2/pathogen_TE_2025_New/HOMER/bin"
export PATH="$HOMER_PATH:$PATH"
GENOME="hg38"
THREADS=6
PARALLEL_JOBS=10  # 10个任务 × 6线程 = 60线程

# 自定义去重后的motif数据库 (去除了与HOMER内置重复的TF)
# HOMER会自动加载内置数据库，所以只需要提供去重后的自定义数据库
# 统一的背景文件
BG_BED="$WORKDIR/background_bed/background_intergenic_TEs.bed"

cd $WORKDIR

# 创建输出目录
mkdir -p HOMER_results_per_subfamily

echo "=== Running HOMER for each TE subfamily ==="
echo "Using unified background: $BG_BED"
echo ""

# 定义处理函数
run_homer_subfamily() {
    local subfamily=$1
    local target_bed="target_bed_per_subfamily/${subfamily}_target.bed"
    local output_dir="HOMER_results_per_subfamily/${subfamily}"

    # 检查target文件是否存在且非空
    if [ ! -s "$target_bed" ]; then
        echo "Skipping $subfamily: empty target file"
        return
    fi

    # 检查loci数量
    target_count=$(wc -l < "$target_bed")
    if [ "$target_count" -lt 5 ]; then
        echo "Skipping $subfamily: only $target_count loci (need >= 5)"
        return
    fi

    echo "Processing: $subfamily ($target_count target loci)"

    mkdir -p "$output_dir"

    # 运行HOMER findMotifsGenome.pl
    ${HOMER_PATH}/findMotifsGenome.pl \
        "$target_bed" \
        $GENOME \
        "$output_dir" \
        -bg "$BG_BED" \
        -size given \
        -p $THREADS \
        -nomotif \
        2>&1 | tee "${output_dir}/homer.log"

    if [ $? -eq 0 ]; then
        echo "✓ Completed: $subfamily"
    else
        echo "✗ Failed: $subfamily"
    fi
}

export -f run_homer_subfamily
export HOMER_PATH GENOME THREADS BG_BED PATH

# 使用parallel运行
cat target_subfamily_list.txt | parallel -j $PARALLEL_JOBS run_homer_subfamily {}

echo ""
echo "=== HOMER analysis completed ==="
echo "Results saved in: HOMER_results_per_subfamily/"
