import sys
import re
import gzip
import collections


class GFFReader:
    """解析GTF文件的核心读取器"""

    def __init__(self, filename, id_attribute):
        self.filename = filename
        self.id_attribute = id_attribute
        self._attr_pattern = re.compile(r'\s*([^\s=]+)[\s=]+"(.*?)"')

    def __iter__(self):
        # 自动处理gzip文件
        opener = gzip.open if self.filename.endswith(('.gz', '.gzip')) else open
        with opener(self.filename, 'rt') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.strip().split('\t')
                if len(parts) != 9:
                    continue

                # 解析字段
                seqname = parts[0]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = self._parse_attributes(parts[8])

                yield (
                    attributes.get(self.id_attribute),  # ID
                    seqname,
                    strand,
                    start,
                    end,
                    feature_type
                )

    def _parse_attributes(self, attr_str):
        """解析属性字段"""
        attrs = {}
        for match in self._attr_pattern.finditer(attr_str):
            key, value = match.groups()
            attrs[key.strip()] = value.strip()
        return attrs


class GTFStructureChecker:
    def __init__(self, gtf_file, id_attribute='gene_id', feature_type='exon', stranded=False):
        self.gtf_file = gtf_file
        self.id_attribute = id_attribute
        self.feature_type = feature_type
        self.stranded = stranded

        # 原始数据结构
        self.temp_plus = collections.defaultdict(lambda: collections.defaultdict(list))
        self.temp_minus = collections.defaultdict(lambda: collections.defaultdict(list))
        self.temp_nostrand = collections.defaultdict(lambda: collections.defaultdict(list))
        self.feature_ids = set()

    def analyze(self):
        """执行分析"""
        reader = GFFReader(self.gtf_file, self.id_attribute)

        processed = 0
        for record in reader:
            feature_id, seqname, strand, start, end, ftype = record

            # 过滤特征类型
            if ftype != self.feature_type:
                continue

            # 原始处理逻辑
            self._process_record(feature_id, seqname, strand, start, end)

            # 进度显示
            processed += 1
            if processed % 100000 == 0:
                print(f"Processed {processed} features...", file=sys.stderr)

        self._show_summary()

    def _process_record(self, feature_id, seqname, strand, start, end):
        """处理单条记录"""
        if not feature_id:
            return

        # 保存feature ID
        self.feature_ids.add(feature_id)

        # 原始分类逻辑
        try:
            if strand == '+':
                self.temp_plus[seqname][feature_id].append((start, end))
            elif strand == '-':
                self.temp_minus[seqname][feature_id].append((start, end))
            else:
                self.temp_nostrand[seqname][feature_id].append((start, end))

            # 严格链模式检查
            if self.stranded and strand == '.':
                print(f"WARNING: Feature {feature_id} missing strand information", file=sys.stderr)

        except Exception as e:
            print(f"Error processing {feature_id}: {str(e)}", file=sys.stderr)
            raise

    def _show_summary(self):
        """显示统计摘要"""
        print("\nGTF Structure Analysis Report")
        print("=" * 40)
        print(f"File: {self.gtf_file}")
        print(f"Feature type analyzed: {self.feature_type}")
        print(f"ID attribute used: {self.id_attribute}")
        print(f"\nBasic Statistics:")
        print(f"Total features processed: {len(self.feature_ids)}")
        print(f"Chromosomes (+ strand): {len(self.temp_plus)}")
        print(f"Chromosomes (- strand): {len(self.temp_minus)}")
        print(f"Chromosomes (no strand): {len(self.temp_nostrand)}")

        # 示例数据展示
        sample_chrom = next(iter(self.temp_plus), 3)  # 显示前3个示例
        if sample_chrom in self.temp_plus:
            print(f"\nSample + strand features on {sample_chrom}:")
            for fid, intervals in list(self.temp_plus[sample_chrom].items())[:3]:
                print(f"  {fid}: {len(intervals)} intervals")


if __name__ == '__main__':
    # 直接在PyCharm中运行时配置参数
    INPUT_FILE = "F:/BioBank/Index/GRCm39/Ensembl/GRCm39_Ensembl.gtf"  # 修改为实际文件路径

    checker = GTFStructureChecker(
        gtf_file=INPUT_FILE,
        id_attribute='gene_id',
        feature_type='exon',
        stranded=True
    )

    try:
        checker.analyze()
    except Exception as e:
        print(f"Analysis failed: {str(e)}", file=sys.stderr)
        sys.exit(1)