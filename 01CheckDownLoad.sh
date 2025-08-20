#!/bin/bash

# 全局配置
THREADS=15
SRR_LIST=(SRR239076{07..21}) 
BASE_DIR="/home/yuetong/YANGPeida/bank/ERV/GSE227647/" 
LOG_FILE="${BASE_DIR}/sra_processing.log"

# 日志初始化
echo "==== 任务启动: $(date) ====" > $LOG_FILE

# download data
download_sra() {
    local srr=$1
    local sra_dir="${BASE_DIR}/${srr}"
    
    local prefetch_dir="${sra_dir}/${srr}"  
    local sra_path="${prefetch_dir}/${srr}.sra"  

    mkdir -p "$sra_dir"

    if [ -f "$sra_path" ]; then
        echo "[$(date)] 发现现有SRA文件: ${srr}" | tee -a $LOG_FILE
        return 0
    fi

    # 执行下载
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

# 转换
convert_to_fastq() {
    local srr=$1
    local sra_path="${BASE_DIR}/${srr}/${srr}/${srr}.sra"
    local output_dir="${BASE_DIR}/fastq_output"
    
    mkdir -p "$output_dir"

    local output_files=(
        "${output_dir}/${srr}.fastq"           
        "${output_dir}/${srr}_1.fastq"       
        "${output_dir}/${srr}_2.fastq"
    )

    # avoid repeat
    if [ -f "${output_files[1]}" ] && [ -f "${output_files[2]}" ]; then
        echo "[$(date)] 检测到完整FASTQ文件: ${srr}" | tee -a $LOG_FILE
        return 0
    elif [ -f "${output_files[0]}" ]; then
        echo "[$(date)] 检测到单端FASTQ文件: ${srr}" | tee -a $LOG_FILE
        return 0
    fi

    # 转换
    echo "[$(date)] 开始转换: ${srr}" | tee -a $LOG_FILE
    fasterq-dump --split-3 --threads $THREADS --outdir "$output_dir" "$sra_path" >> $LOG_FILE 2>&1

    # 结果验证
    if [ $? -eq 0 ]; then
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

export BASE_DIR LOG_FILE THREADS
export -f download_sra convert_to_fastq

# 并行下载
printf "%s\n" "${SRR_LIST[@]}" | xargs -P 2 -I {} bash -c 'download_sra "$@"' _ {}

# 顺序转换
for srr in "${SRR_LIST[@]}"; do
    convert_to_fastq "$srr"
done

echo "==== 任务完成: $(date) ====" | tee -a $LOG_FILE
