#!/bin/bash
set -euo pipefail

# ------------ 参数设置 ------------
FASTQ_DIR="/home/yangpeida/BanqueDonee/HERV/GSE227647/fastq_output"
HISAT_INDEX="/home/yangpeida/BanqueDonee/index/mm10/mm10"  
ALIGN_OUTPUT_DIR="/home/yangpeida/BanqueDonee/HERV/GSE227647/ERVbamHisat"
THREADS_HISAT=15
THREADS_SORT=15
ERROR_LOG="${ALIGN_OUTPUT_DIR}/alignment_errors.log"

# ------------ 目录检查 ------------
mkdir -p "${ALIGN_OUTPUT_DIR}"
touch "${ERROR_LOG}"

cd "${FASTQ_DIR}" || { echo "无法进入FASTQ目录: ${FASTQ_DIR}"; exit 1; }

samples=$(ls *_1.fastq | sed 's/_1.fastq//')

# ------------ 循环比对 ------------
for sample in ${samples}; do
    # 最终输出文件
    sorted_bam="${ALIGN_OUTPUT_DIR}/${sample}_sorted.bam"
    if [[ -s "${sorted_bam}" ]]; then
        echo "样本 ${sample} 已处理，跳过..."
        continue
    fi

    echo "正在处理样本: ${sample}..."
    
    # 临时文件
    sam_output="${ALIGN_OUTPUT_DIR}/${sample}.sam"
    unsorted_bam="${ALIGN_OUTPUT_DIR}/${sample}_unsorted.bam"

    # 使用if判断捕获错误
    if (
        # Hisat2比对步骤
        hisat2 -x "${HISAT_INDEX}" \
            -1 "${sample}_1.fastq" \
            -2 "${sample}_2.fastq" \
            -p ${THREADS_HISAT} \
            -S "${sam_output}" \
            2>> "${ERROR_LOG}"
        
        # 转换SAM到BAM
        samtools view -@ ${THREADS_SORT} -Sb "${sam_output}" > "${unsorted_bam}"
        
        # 排序BAM
        samtools sort -@ ${THREADS_SORT} -o "${sorted_bam}" "${unsorted_bam}"
        
        # 创建索引
        samtools index -@ ${THREADS_SORT} "${sorted_bam}"
        
        # 清理中间文件
        rm "${sam_output}" "${unsorted_bam}"
    ) ; then
        echo "样本 ${sample} 比对完成。输出文件: ${sorted_bam}"
    else
        echo "样本 ${sample} 处理失败！详情见日志: ${ERROR_LOG}"
        # 删除可能生成的不完整文件
        rm -f "${sam_output}" "${unsorted_bam}" "${sorted_bam}"* 2>/dev/null || true
    fi
done

echo "全部样本处理完成！检查错误日志: ${ERROR_LOG}"
