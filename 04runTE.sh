#!/bin/bash
set -euo pipefail

# -------------------- 参数配置 -------------------- 
BAM_DIR="/home/yuetong/YANGPeida/bank/ERV/GSE227647/En150bamSTAR"
OUTPUT_DIR="/home/yuetong/YANGPeida/bank/ERV/GSE227647/TEtranscripts_results"
GENOME_GTF="/home/yuetong/YANGPeida/bank/index/Ensembl/GRCm39/GRCm39_Ensembl.gtf"
TE_GTF="/home/yuetong/YANGPeida/bank/index/Ensembl/GRCm39/TE/GRCm39_Ensembl_rmsk_TE.gtf.ind"
CONTROL_PREFIXES=("SRR23907619" "SRR23907620" "SRR23907621")

# -------------------- 初始化检查 -------------------- 
mkdir -p "${OUTPUT_DIR}"
declare -a treatment_bams control_bams

# 智能样本分类（支持任意BAM后缀）
while IFS= read -r -d '' bam; do
    sample=$(basename "$bam" | grep -oE '^[^_]+')
    if [[ " ${CONTROL_PREFIXES[@]} " =~ " $sample " ]]; then
        control_bams+=("$bam")
    else
        treatment_bams+=("$bam")
    fi
done < <(find "${BAM_DIR}" -maxdepth 1 -name "*.bam" -print0)

# -------------------- 执行TEtranscripts（兼容旧版参数）------------------- 
TEtranscripts --mode multi \
    -t "${treatment_bams[@]}" \
    -c "${control_bams[@]}" \
    --GTF "${GENOME_GTF}" \
    --TE "${TE_GTF}" \
    --project "${OUTPUT_DIR}/TEtranscripts_out" \
    --sortByPos

# -------------------- 结果合并 -------------------- 
# 使用首行标题+合并计数列
awk 'BEGIN {OFS="\t"; print "TE_ID"} 
    FNR==1 {next}
    {count[$1]=$2} 
    END {for (te in count) print te, count[te]}' \
    "${OUTPUT_DIR}"/TEtranscripts_out/*.cntTable \
    > "${OUTPUT_DIR}/combined_ERV_counts.tsv"

echo "分析完成！结果文件: ${OUTPUT_DIR}/combined_ERV_counts.tsv"
