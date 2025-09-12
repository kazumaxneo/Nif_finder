
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import glob
import csv
import os
from collections import defaultdict

# ターゲット遺伝子
NIF_GENES = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"]
VNF_GENES = ["VnfH", "VnfD", "VnfK"]
STATUSES = ["Full", "Full_operon", "Fragment"]

# 出力ファイル
output_file = "nif_vnf_detailed_summary.tsv"

# 出力ヘッダー（nifHDKENBごとに Full→Full_operon→Fragment の順）
header = []
for status in STATUSES:
    for g in NIF_GENES:
        header.append(f"{g}_{status}")
header.append("nif_total")

for status in STATUSES:
    for g in VNF_GENES:
        header.append(f"{g}_{status}")
header.append("vnf_total")

header.append("Category")  # potential_diazotroph or no_hit or blank

# 集計
with open(output_file, "w") as out_f:
    writer = csv.writer(out_f, delimiter="\t")
    writer.writerow(["Genome"] + header)

    for file in sorted(glob.glob("*result*")):
        filename = os.path.basename(file)
        genome_name = os.path.splitext(filename)[0]
        counts = defaultdict(int)

        with open(file) as f:
            for line in f:
                if line.startswith("Query") or line.strip() == "":
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                prediction = parts[4]
                completeness = parts[5]

                if prediction in NIF_GENES or prediction in VNF_GENES:
                    key = f"{prediction}_{completeness}"
                    counts[key] += 1

        # 合計を計算
        nif_total = sum(counts[f"{g}_{s}"] for g in NIF_GENES for s in STATUSES)
        vnf_total = sum(counts[f"{g}_{s}"] for g in VNF_GENES for s in STATUSES)

        # カテゴリ分類
        all_full = all(counts[f"{g}_Full"] >= 1 for g in NIF_GENES)
        all_no_hit = all(sum(counts[f"{g}_{s}"] for s in STATUSES) == 0 for g in NIF_GENES)
        if all_full:
            category = "potential_diazotroph"
        elif all_no_hit:
            category = "no_hit"
        else:
            category = ""

        # 出力行作成
        row = [genome_name]
        for s in STATUSES:
            for g in NIF_GENES:
                row.append(counts[f"{g}_{s}"])
        row.append(nif_total)

        for s in STATUSES:
            for g in VNF_GENES:
                row.append(counts[f"{g}_{s}"])
        row.append(vnf_total)

        row.append(category)

        writer.writerow(row)

print(f"Detailed summary written to {output_file}")
