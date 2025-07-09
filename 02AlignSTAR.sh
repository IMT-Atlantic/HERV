#!/bin/bash
set -euo pipefail

# ------------ 参数设置 ------------
# MultiThread is forbidden
FASTQ_DIR="/home/yuetong/YANGPeida/bank/ERV/GSE227647/fastq_output"
STAR_INDEX_DIR="/home/yuetong/YANGPeida/bank/index/Ensembl/GRCm39/STAR150"
ALIGN_OUTPUT_DIR="/home/yuetong/YANGPeida/bank/ERV/GSE227647/En150bamSTAR"
THREADS_STAR=1
THREADS_SORT=1
ERROR_LOG="${ALIGN_OUTPUT_DIR}/alignment_errors.log"  # 错误日志文件

# ------------ 目录检查 ------------
mkdir -p "${ALIGN_OUTPUT_DIR}"
touch "${ERROR_LOG}"  # 确保错误日志存在

# 进入FASTQ目录
cd "${FASTQ_DIR}" || { echo "无法进入FASTQ目录: ${FASTQ_DIR}"; exit 1; }

# 提取所有样本名前缀（如SRR23907607）
samples=$(ls *_1.fastq | sed 's/_1.fastq//')

# ------------ 循环比对 ------------
for sample in ${samples}; do
    # 检查是否已存在排序后的BAM文件（标志已完成）
    output_bam="${ALIGN_OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
    if [[ -s "${output_bam}" ]]; then
        echo "样本 ${sample} 已处理，跳过..."
        continue
    fi

    echo "正在处理样本: ${sample}..."
    
    # 使用if判断捕获错误，|| true防止set -e导致脚本终止
    if STAR \
        --genomeDir "${STAR_INDEX_DIR}" \
        --readFilesIn "${sample}_1.fastq" "${sample}_2.fastq" \
        --runThreadN "${THREADS_STAR}" \
        --outFileNamePrefix "${ALIGN_OUTPUT_DIR}/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN "${THREADS_SORT}" \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 100 \
        --outSAMattributes All \
        --twopassMode Basic \
        2>> "${ERROR_LOG}"; then  # 将错误信息追加到日志
        echo "样本 ${sample} 比对完成。输出文件: ${ALIGN_OUTPUT_DIR}/${sample}_*"
    else
        echo "样本 ${sample} 处理失败！详情见日志: ${ERROR_LOG}"
        # 可选：删除可能生成的不完整文件
        rm -f "${ALIGN_OUTPUT_DIR}/${sample}_*" 2>/dev/null || true
    fi
done

echo "全部样本处理完成！检查错误日志: ${ERROR_LOG}"
