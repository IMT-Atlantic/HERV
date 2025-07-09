#!/bin/bash
# SRA智能处理脚本（完全兼容现有目录结构）
# SRA提取并改写为fastq文件

# 全局配置
THREADS=15
SRR_LIST=(SRR239076{07..21}) 
BASE_DIR="/home/yuetong/YANGPeida/bank/ERV/GSE227647/"  # 固定基础路径
LOG_FILE="${BASE_DIR}/sra_processing.log"

# 日志初始化
echo "==== 任务启动: $(date) ====" > $LOG_FILE

# 函数：智能下载器（适配嵌套目录）
download_sra() {
    local srr=$1
    local sra_dir="${BASE_DIR}/${srr}"
    
    # 预生成的目标路径（适配prefetch的嵌套结构）
    local prefetch_dir="${sra_dir}/${srr}"  # 新增层级
    local sra_path="${prefetch_dir}/${srr}.sra"  # 实际存储路径

    mkdir -p "$sra_dir"

    # 存在性检查（适配新路径）
    if [ -f "$sra_path" ]; then
        echo "[$(date)] 发现现有SRA文件: ${srr}" | tee -a $LOG_FILE
        return 0
    fi

    # 执行下载（允许prefetch生成嵌套目录）
    echo "[$(date)] 开始下载: ${srr}" | tee -a $LOG_FILE
    if prefetch -O "$sra_dir" "$srr" >> $LOG_FILE 2>&1; then
        # 确认嵌套目录结构
        if [ -f "$sra_path" ]; then
            echo "[$(date)] 下载成功: ${srr}" | tee -a $LOG_FILE
            return 0
        else
            # 处理可能的路径异常
            echo "[ERROR] 文件路径异常，尝试修复..." | tee -a $LOG_FILE
            local found_path=$(find "$sra_dir" -name "${srr}.sra" -print -quit)
            if [ -f "$found_path" ]; then
                echo "[WARNING] 检测到非常规路径，自动迁移文件: $found_path -> $sra_path" | tee -a $LOG_FILE
                mkdir -p "$prefetch_dir"
                mv "$found_path" "$prefetch_dir/"
                return 0
            else
                echo "[ERROR] 无法定位下载文件: ${srr}" | tee -a $LOG_FILE
                return 1
            fi
        fi
    else
        echo "[ERROR] 下载失败: ${srr}" | tee -a $LOG_FILE
        return 1
    fi
}

# 函数：智能转换器
convert_to_fastq() {
    local srr=$1
    local sra_path="${BASE_DIR}/${srr}/${srr}/${srr}.sra"
    local output_dir="${BASE_DIR}/fastq_output"
    
    # 创建输出目录
    mkdir -p "$output_dir"

    # 结果文件检查（兼容单端/双端）
    local output_files=(
        "${output_dir}/${srr}.fastq"           # 单端情况
        "${output_dir}/${srr}_1.fastq"         # 双端情况
        "${output_dir}/${srr}_2.fastq"
    )

    # 已存在文件检测
    if [ -f "${output_files[1]}" ] && [ -f "${output_files[2]}" ]; then
        echo "[$(date)] 检测到完整FASTQ文件: ${srr}" | tee -a $LOG_FILE
        return 0
    elif [ -f "${output_files[0]}" ]; then
        echo "[$(date)] 检测到单端FASTQ文件: ${srr}" | tee -a $LOG_FILE
        return 0
    fi

    # 执行转换（自动检测布局）
    echo "[$(date)] 开始转换: ${srr}" | tee -a $LOG_FILE
    fasterq-dump --split-3 --threads $THREADS --outdir "$output_dir" "$sra_path" >> $LOG_FILE 2>&1

    # 结果验证
    if [ $? -eq 0 ]; then
        # 检测实际生成的文件
        local generated_files=($(ls "${output_dir}/${srr}"*fastq 2>/dev/null))
        if [ ${#generated_files[@]} -ge 1 ]; then
            echo "[$(date)] 转换成功: 生成 ${#generated_files[@]} 个文件" | tee -a $LOG_FILE
            return 0
        else
            echo "[ERROR] 无输出文件生成: ${srr}" | tee -a $LOG_FILE
            return 1
        fi
    else
        echo "[ERROR] 转换失败: ${srr}" | tee -a $LOG_FILE
        return 1
    fi
}

# 主流程控制
export BASE_DIR LOG_FILE THREADS
export -f download_sra convert_to_fastq

# 阶段1：并行下载（限制同时下载数）
printf "%s\n" "${SRR_LIST[@]}" | xargs -P 2 -I {} bash -c 'download_sra "$@"' _ {}

# 阶段2：顺序转换（避免IO瓶颈）
for srr in "${SRR_LIST[@]}"; do
    convert_to_fastq "$srr"
done

echo "==== 任务完成: $(date) ====" | tee -a $LOG_FILE
