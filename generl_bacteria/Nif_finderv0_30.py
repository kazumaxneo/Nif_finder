#!/usr/bin/env python3

from Bio import SeqIO
import subprocess
import math
import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import sqrt
import argparse
import glob
import tempfile
import numpy as np
from xml.sax.saxutils import escape
# v0.21 改善内容:
#   1. 1-NNのスケーリング: z-scoreによる正規化（参照データロード時に一度だけ計算）
#   2. Full_operon判定の再設計:
#      - 同一クエリに2種類以上の異なるnif遺伝子がヒット
#      - かつ、それらのヒットのalignment_lengthが遺伝子ごとの閾値以上
#      両方を満たす場合にFull_operonと判定（ポスト処理として実装）
#   3. --plot オプション追加: alignment_length vs -log10(Evalue) の散布図を出力
#      参照データ（灰色）とクエリヒット（遺伝子ごとに色分け）を6パネルで表示
# v0.22 追加:
#   4. -g / --genome オプション追加: ゲノムDNA FASTAを直接クエリに指定可能
#      内部で6フレーム翻訳を行い、一時FASTAを生成してhmmscanに渡す
#      翻訳後ORFの最小長は --min_orf_len で指定（デフォルト: 10 aa）
#      -q / -d / -g は排他オプション
# v0.23 修正:
#   5. [バグ修正] Full_operonクエリのFragment重複プロット問題を解消:
#      plot_scatter にoperon_queriesを渡し、Full_operonクエリは自パネル（★）のみ表示、
#      他パネルへの△重複表示をスキップ
#   6. 散布図のX軸をアライメント長からタンパク質長（query_full_length）に変更:
#      クエリ点のX座標をタンパク質長で描画し、どのようなタンパク質がヒットしたか推定可能に
#      参照データはalignment_lengthのまま維持
#   7. 1-NN距離閾値による「判定不能」(unclassifiable) 追加:
#      z-score正規化後の最近傍距離が NN_DISTANCE_THRESHOLD を超えた場合、
#      nif/otherではなく "unclassifiable" と判定。散布図では◇マーカーで表示
# v0.24 追加:
#   8. NIF_FINDER_DB 環境変数でDBルートを指定すると、
#      -t / -r 未指定時に標準nifモデルと分類参照ファイルを自動検出
# v0.25 追加:
#   9. GenBankからnif-encoding regionを抽出するCLI機能を追加
#      (--genbank / --genbank_dir / --context_size_kb)
#      注釈付きGBK、whole-genome overview SVG、local context SVGを出力
# v0.30 追加:
#  10. vnf/vup関連モデルを拡張し、vnfHDGKEN/vupABC配列出力オプションを追加
#  11. --save_vnf_hdgken_vupabc_fasta ではnifモデルを除外する排他的vnf/vup専用モードで解析

NIF_GENES = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"]
VNF_GENES = ["vnfH/nifH", "vnfD", "vnfK", "vnfE/nifE", "vnfN/nifN", "vnfG", "vnfDG"]
ACCESSORY_GENES = [
    "nifZ", "nifX", "nifP/cysE", "nifT", "nifV", "nifS", "nifU", "nifU_like",
    "modB/vupB", "modC/vupC", "modAlike", "vupA/modA", "vupB/modB",
    "vupC/modC", "cnfR/patB", "cnfR/patB_like",
]
TARGET_GENES = NIF_GENES + VNF_GENES + ACCESSORY_GENES
VNF_HDGKEN_VUPABC_GENES = {
    "vnfH/nifH", "vnfD", "vnfDG", "vnfG", "vnfK", "vnfE/nifE", "vnfN/nifN",
    "vupA/modA", "vupB/modB", "vupC/modC",
}
VNF_CONTEXT_ANCHOR_GENES = {"vnfG", "vnfDG"}
VNF_CONTEXT_FILTER_GENES = {"vnfH/nifH", "vnfD", "vnfK", "vnfE/nifE", "vnfN/nifN"}
DEFAULT_VNF_CONTEXT_DISTANCE_BP = 20000
VNF_REGION_GBK_GENES = VNF_CONTEXT_ANCHOR_GENES | VNF_CONTEXT_FILTER_GENES
VNF_EXCLUSIVE_MODEL_GENES = ["vnfH", "vnfD", "vnfK", "vnfE", "vnfN", "vnfG", "vupA", "vupB", "vupC"]
MODEL_GENES = [
    "nifH", "nifD", "nifK",
    "nifE", "nifN", "nifB",
    "nifZ", "nifX", "nifP",
    "nifT", "nifV",
    "nifS", "nifU",
    "vnfH", "vnfD", "vnfK",
    "vnfE", "vnfN", "vnfG",
    "modA", "modB", "modC",
    "vupA", "vupB", "vupC",
    "cnfR-patB",
]
PANEL_GENES = MODEL_GENES
FUSION_GENES = {"vnfDG"}
VNF_DG_COMPONENTS = {"nifD", "vnfD"}
OPERON_CANDIDATE_GENES = set(NIF_GENES + ["vnfD", "vnfK"])
GENE_PANEL_MAP = {
    "vnfDG": "vnfG",
    "vnfH/nifH": "vnfH",
    "vnfE/nifE": "vnfE",
    "vnfN/nifN": "vnfN",
    "nifP/cysE": "nifP",
    "nifU_like": "nifU",
    "modB/vupB": "modB",
    "modC/vupC": "modC",
    "modAlike": "modA",
    "vupA/modA": "vupA",
    "vupB/modB": "vupB",
    "vupC/modC": "vupC",
    "cnfR/patB": "cnfR-patB",
    "cnfR/patB_like": "cnfR-patB",
}
TARGET_GENE_BY_LOWER = {gene.lower(): gene for gene in TARGET_GENES}
NIF_FINDER_DB_ENV = "NIF_FINDER_DB"
DEFAULT_PROFILE_NAME = "proteins_hmm"
DEFAULT_REFERENCE_SUFFIX = "classification"

DEFAULT_E_VALUE = 1e-10
DEFAULT_CPU = 8
DEFAULT_JOBS = 1
DEFAULT_LOG_E_VALUE_MIN = 250
DEFAULT_MIN_ORF_LEN = 10  # 6フレーム翻訳後に保持するORFの最小アミノ酸長
DEFAULT_CONTEXT_SIZE_KB = 10
MAX_CONTEXT_SIZE_KB = 30
LOCAL_CONTEXT_MERGE_DISTANCE = 20000
OPERON_LENGTH_MULTIPLIER = 1.6

NIF_GENE_THRESHOLDS = {
    "nifH": 240,
    "nifD": 370,
    "nifK": 410,
    "nifE": 400,
    "nifN": 400,
    "nifB": 370,
    "nifZ": 50,
    "nifX": 100,
    "vnfD": 200,
    "vnfK": 400,
    "vnfG": 100,
    "vnfDG": 200,
    "vnfH": 240,
    "vnfH/nifH": 240,
    "vnfE": 400,
    "vnfE/nifE": 400,
    "vnfN": 400,
    "vnfN/nifN": 400,
    "nifP": 200,
    "nifP/cysE": 200,
    "nifT": 50,
    "nifV": 300,
    "nifS": 350,
    "nifU": 250,
    "nifU_like": 250,
    "modB": 200,
    "modB/vupB": 200,
    "modC": 200,
    "modC/vupC": 200,
    "modA": 200,
    "modAlike": 200,
    "vupA": 200,
    "vupA/modA": 200,
    "vupB": 200,
    "vupB/modB": 200,
    "vupC": 200,
    "vupC/modC": 200,
    "cnfR-patB": 400,
    "cnfR/patB": 400,
    "cnfR/patB_like": 400,
}

PANEL_X_LIMITS = {
    "nifT": 500,
    "nifV": 1000,
    "nifS": 1000,
    "nifU": 1000,
}

NIF_GENE_COLORS = {
    "nifH": "#E74C3C",
    "nifD": "#3498DB",
    "nifK": "#2ECC71",
    "nifE": "#F39C12",
    "nifN": "#9B59B6",
    "nifB": "#1ABC9C",
    "nifZ": "#0B85D4",
    "nifX": "#D81B60",
    "vnfD": "#1F77B4",
    "vnfK": "#1E8449",
    "vnfG": "#130F73",
    "vnfDG": "#0F0B5C",
    "vnfH": "#B03A2E",
    "vnfH/nifH": "#B03A2E",
    "vnfE": "#B9770E",
    "vnfE/nifE": "#B9770E",
    "vnfN": "#7D3C98",
    "vnfN/nifN": "#7D3C98",
    "nifP": "#FB8CFF",
    "nifP/cysE": "#FB8CFF",
    "nifT": "#8E44AD",
    "nifV": "#6C8E23",
    "nifS": "#CFEA35",
    "nifU": "#599E87",
    "nifU_like": "#7DBBA7",
    "modB": "#F19AA0",
    "modB/vupB": "#F19AA0",
    "modC": "#E86F7A",
    "modC/vupC": "#E86F7A",
    "modA": "#E072BF",
    "modAlike": "#E072BF",
    "vupA": "#B65FCF",
    "vupA/modA": "#B65FCF",
    "vupB": "#C75B73",
    "vupB/modB": "#C75B73",
    "vupC": "#D46F5C",
    "vupC/modC": "#D46F5C",
    "cnfR-patB": "#4D4D4D",
    "cnfR/patB": "#4D4D4D",
    "cnfR/patB_like": "#7A7A7A",
}

PROFILE_ALLOWED_PREDICTIONS = {
    "nifP": {"nifP/cysE", "other"},
    "nifT": {"nifT", "other"},
    "nifV": {"nifV", "other"},
    "nifZ": {"nifZ", "other"},
    "nifX": {"nifX", "other"},
    "nifS": {"nifS", "other"},
    "vnfK": {"vnfK", "other"},
    "vnfD": {"vnfD", "other"},
    "vnfH": {"vnfH/nifH", "other"},
    "vnfE": {"vnfE/nifE", "other"},
    "vnfN": {"vnfN/nifN", "other"},
    "nifU": {"nifU", "nifU_like", "other"},
    "modB": {"modB/vupB", "other"},
    "modC": {"modC/vupC", "other"},
    "modA": {"modAlike", "other"},
    "vupA": {"vupA/modA", "other"},
    "vupB": {"vupB/modB", "other"},
    "vupC": {"vupC/modC", "other"},
    "cnfR-patB": {"cnfR/patB", "cnfR/patB_like", "other"},
}

# z-score正規化後の1-NN距離がこの閾値を超えた場合、"unclassifiable"と判定
# 参照クラスタから著しく離れた点（nifENオペロンがnifD/nifKパネルに誤ヒット等）を除外する。
# 値はz-score空間でのユークリッド距離。参照データの広がりに応じて調整可能。
# デフォルト 2.0: nifENとnifDK参照クラスタ間の距離が概ね3〜5 z-score程度を想定。
NN_DISTANCE_THRESHOLD = 2.0


def normalize_prediction_attribute(attribute):
    value = (attribute or "").strip()
    return TARGET_GENE_BY_LOWER.get(value.lower(), value)


def normalize_profile_prediction(profile_gene, attribute):
    if profile_gene == "vnfG" and attribute == "vnfD":
        return "vnfDG"
    allowed = PROFILE_ALLOWED_PREDICTIONS.get(profile_gene)
    if allowed and attribute not in allowed and attribute != "unclassifiable":
        return "other"
    return attribute


# ============================================================
# 標準DBパス解決
# ============================================================

def build_default_model_paths(db_root):
    """
    NIF_FINDER_DB で指定されたDBルートから標準モデル一式を構築する。

    期待する構成:
      $NIF_FINDER_DB/nifH/proteins_hmm
      $NIF_FINDER_DB/nifH/nifHclassification
      $NIF_FINDER_DB/vnfG/proteins_hmm
      $NIF_FINDER_DB/vnfG/vnfGclassification
      ...
    """
    profile_files = []
    reference_files = []
    for gene in MODEL_GENES:
        gene_dir = os.path.join(db_root, gene)
        reference_gene = "cnfR_patB" if gene == "cnfR-patB" else gene
        profile_files.append(os.path.join(gene_dir, DEFAULT_PROFILE_NAME))
        reference_files.append(os.path.join(gene_dir, f"{reference_gene}{DEFAULT_REFERENCE_SUFFIX}"))
    return profile_files, reference_files


def profile_gene_from_path(profile_file):
    return os.path.basename(os.path.dirname(profile_file))


def filter_model_paths(profile_files, reference_files, model_genes):
    wanted = set(model_genes)
    order = {gene: idx for idx, gene in enumerate(model_genes)}
    pairs = [
        (profile_gene_from_path(profile_file), profile_file, reference_file)
        for profile_file, reference_file in zip(profile_files, reference_files)
        if profile_gene_from_path(profile_file) in wanted
    ]
    pairs.sort(key=lambda item: order[item[0]])
    return [item[1] for item in pairs], [item[2] for item in pairs]


def apply_vnf_exclusive_model_filter(profile_files, reference_files, parser):
    filtered_profiles, filtered_references = filter_model_paths(
        profile_files, reference_files, VNF_EXCLUSIVE_MODEL_GENES
    )
    if not filtered_profiles:
        parser.error(
            "--save_vnf_hdgken_vupabc_fasta requires standard vnf/vup model folders "
            "under NIF_FINDER_DB."
        )
    return filtered_profiles, filtered_references


def panel_genes_for_profiles(profile_files):
    return [profile_gene_from_path(profile_file) for profile_file in profile_files]


def validate_model_files(profile_files, reference_files):
    missing = []
    for path in profile_files + reference_files:
        if not os.path.isfile(path):
            missing.append(path)
    return missing


def resolve_model_paths(profile_files, reference_files, parser):
    """
    明示指定された -t / -r を優先し、両方未指定の場合は NIF_FINDER_DB から補完する。
    """
    if profile_files or reference_files:
        if not profile_files or not reference_files:
            parser.error("Please specify both --profile and --reference, or set NIF_FINDER_DB and omit both.")
        if len(profile_files) != len(reference_files):
            parser.error("The number of --profile files must match the number of --reference files.")
        missing = validate_model_files(profile_files, reference_files)
        if missing:
            parser.error("Model file(s) not found:\n  " + "\n  ".join(missing))
        return profile_files, reference_files

    db_root = os.environ.get(NIF_FINDER_DB_ENV)
    if not db_root:
        parser.error(
            "Please specify --profile/--reference, or set NIF_FINDER_DB to the directory "
            "containing standard target model folders."
        )

    db_root = os.path.abspath(os.path.expanduser(db_root))
    profile_files, reference_files = build_default_model_paths(db_root)
    missing = validate_model_files(profile_files, reference_files)
    if missing:
        parser.error(
            f"{NIF_FINDER_DB_ENV} is set to {db_root}, but required model file(s) are missing:\n  "
            + "\n  ".join(missing)
        )
    return profile_files, reference_files


# ============================================================
# HMMscan実行
# ============================================================

def run_hmmscan(query_file, target_file, output_file, e_value=DEFAULT_E_VALUE, cpu=DEFAULT_CPU):
    command = [
        "hmmscan",
        "--domtblout", output_file,
        "--cpu", str(cpu),
        "-E", str(e_value),
        target_file,
        query_file
    ]
    try:
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmscan: {e}")
        print(f"Command: {' '.join(command)}")
        raise


# ============================================================
# 参照データのロード & z-scoreパラメータ計算
# ============================================================

def load_reference_data(reference_file):
    """
    参照TSVを読み込み、z-scoreスケーリング用の統計量を同時に計算して返す。
    戻り値: (reference_data, scale_params)
      reference_data : list of dict {"x", "y", "attribute"}
      scale_params   : dict {"mean_x", "std_x", "mean_y", "std_y"}
    """
    reference_data = []
    try:
        with open(reference_file, 'r') as ref_file:
            reader = csv.DictReader(ref_file, delimiter='\t')
            for row in reader:
                try:
                    alignment_length = int(row["length"])
                    e_value = float(row["Evalue"])
                    log_e_value = -math.log10(e_value) if e_value > 0 else DEFAULT_LOG_E_VALUE_MIN
                    attribute = normalize_prediction_attribute(row["attribute"])
                    reference_data.append({"x": alignment_length, "y": log_e_value, "attribute": attribute})
                except (ValueError, KeyError) as e:
                    print(f"Skipping malformed row in {reference_file}: {row}. Error: {e}")
    except FileNotFoundError:
        print(f"Reference file not found: {reference_file}")
        raise

    if reference_data:
        xs = np.array([r["x"] for r in reference_data], dtype=float)
        ys = np.array([r["y"] for r in reference_data], dtype=float)
        scale_params = {
            "mean_x": float(np.mean(xs)),
            "std_x":  float(np.std(xs)) or 1.0,
            "mean_y": float(np.mean(ys)),
            "std_y":  float(np.std(ys)) or 1.0,
        }
    else:
        scale_params = {"mean_x": 0.0, "std_x": 1.0, "mean_y": 0.0, "std_y": 1.0}

    return reference_data, scale_params


# ============================================================
# 1-NN分類（z-score正規化済み距離）
# ============================================================

def find_nearest_attribute(alignment_length, log_e_value, reference_data, scale_params):
    mean_x = scale_params["mean_x"]
    std_x  = scale_params["std_x"]
    mean_y = scale_params["mean_y"]
    std_y  = scale_params["std_y"]

    norm_qx = (alignment_length - mean_x) / std_x
    norm_qy = (log_e_value      - mean_y) / std_y

    min_distance = float('inf')
    nearest_attribute = "unknown"
    for ref in reference_data:
        rx = (ref["x"] - mean_x) / std_x
        ry = (ref["y"] - mean_y) / std_y
        distance = sqrt((norm_qx - rx) ** 2 + (norm_qy - ry) ** 2)
        if distance < min_distance:
            min_distance = distance
            nearest_attribute = ref["attribute"]

    # 最近傍参照点からの距離が閾値を超える場合は判定不能とする
    # （パラログ間の誤判定、例: nifENオペロンがnifD/nifKパネルに誤ヒット 等を防ぐ）
    if min_distance > NN_DISTANCE_THRESHOLD:
        nearest_attribute = "unclassifiable"

    return nearest_attribute, min_distance


# ============================================================
# FASTAからクエリ長取得
# ============================================================

def get_query_lengths(fasta_file):
    lengths = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            lengths[record.id] = len(record.seq)
    except FileNotFoundError:
        print(f"FASTA file not found: {fasta_file}")
        raise
    except Exception as e:
        print(f"Error parsing FASTA file {fasta_file}: {e}")
        raise
    return lengths


# ============================================================
# domtblout解析 → レコード生成
# Full/Fragmentのみ判定（Full_operonはポスト処理）
# ============================================================

def convert_to_single_tab(input_file, reference_data, scale_params, query_lengths, profile_gene=None):
    records = []
    try:
        with open(input_file, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                columns = line.strip().split()
                if len(columns) < 18:
                    continue
                try:
                    query_name = columns[3]
                    e_value = float(columns[6])
                    aln_from = int(columns[15])
                    aln_to = int(columns[16])
                    alignment_length = abs(aln_to - aln_from) + 1
                    log_e_value = -math.log10(e_value) if e_value > 0 else DEFAULT_LOG_E_VALUE_MIN
                    query_full_length = query_lengths.get(query_name, "N/A")
                    nearest_attribute, nn_distance = find_nearest_attribute(
                        alignment_length, log_e_value, reference_data, scale_params
                    )
                    nearest_attribute = normalize_profile_prediction(profile_gene, nearest_attribute)
                    gene_status = "N/A"
                    if nearest_attribute == "unclassifiable":
                        gene_status = "unclassifiable"
                    elif nearest_attribute in NIF_GENE_THRESHOLDS:
                        threshold = NIF_GENE_THRESHOLDS[nearest_attribute]
                        gene_status = "Full" if alignment_length >= threshold else "Fragment"

                    record = {
                        "query_name":        query_name,
                        "log_e_value":       log_e_value,
                        "alignment_length":  alignment_length,
                        "query_full_length": query_full_length,
                        "prediction":        nearest_attribute,
                        "gene_status":       gene_status,
                        "nn_distance":       nn_distance,
                        "profile_gene":      profile_gene,
                    }
                    records.append(record)
                except (ValueError, IndexError) as e:
                    print(f"Skipping malformed HMMscan output line: {line.strip()}. Error: {e}")
    except FileNotFoundError:
        print(f"HMMscan output file not found: {input_file}")
        raise
    return records


# ============================================================
# profileごとのHMMscan実行
# ============================================================

def _scan_single_profile(query_file, target_file, reference_data, scale_params, query_lengths, cpu, e_value):
    profile_gene = os.path.basename(os.path.dirname(target_file))
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as tmp:
        tmp_out = tmp.name
    try:
        run_hmmscan(query_file, target_file, tmp_out, e_value=e_value, cpu=cpu)
        return convert_to_single_tab(tmp_out, reference_data, scale_params, query_lengths, profile_gene=profile_gene)
    finally:
        if os.path.exists(tmp_out):
            os.remove(tmp_out)


def run_profile_scans(query_file, profile_files, references_loaded, query_lengths, cpu,
                      jobs=DEFAULT_JOBS, e_value=DEFAULT_E_VALUE):
    """
    profileごとのhmmscanを実行し、既定では従来どおり逐次処理する。
    jobs > 1 の場合のみ、独立したprofile scanを並列実行する。

    --cpu は総CPU予算として扱い、並列scan数で割って各hmmscanに渡す。
    HMM profile、E-value、分類ロジックは変更しない。
    """
    jobs = max(1, int(jobs))
    if jobs == 1 or len(profile_files) <= 1:
        all_records = []
        for target_file, (reference_data, scale_params) in zip(profile_files, references_loaded):
            try:
                all_records.extend(
                    _scan_single_profile(query_file, target_file, reference_data,
                                         scale_params, query_lengths, cpu, e_value)
                )
            except Exception as e:
                print(f"Processing of {query_file} with {target_file} failed: {e}")
        return all_records

    workers = min(jobs, len(profile_files))
    cpu_per_scan = max(1, cpu // workers)
    records_by_profile = [[] for _ in profile_files]

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for idx, (target_file, (reference_data, scale_params)) in enumerate(
                zip(profile_files, references_loaded)):
            future = executor.submit(
                _scan_single_profile, query_file, target_file, reference_data,
                scale_params, query_lengths, cpu_per_scan, e_value
            )
            futures[future] = (idx, target_file)

        for future in as_completed(futures):
            idx, target_file = futures[future]
            try:
                records_by_profile[idx] = future.result()
            except Exception as e:
                print(f"Processing of {query_file} with {target_file} failed: {e}")

    all_records = []
    for records in records_by_profile:
        all_records.extend(records)
    return all_records


# ============================================================
# ベストレコード選択
# ============================================================

def _is_full_hit(record):
    prediction = record.get("prediction")
    threshold = NIF_GENE_THRESHOLDS.get(prediction)
    return threshold is not None and record.get("alignment_length", 0) >= threshold


def build_vnfdg_fusion_record(records):
    d_hits = [record for record in records if record.get("prediction") in VNF_DG_COMPONENTS and _is_full_hit(record)]
    g_hits = [record for record in records if record.get("prediction") == "vnfG" and _is_full_hit(record)]
    if not d_hits or not g_hits:
        return None

    best_d = max(d_hits, key=lambda record: record["log_e_value"])
    best_g = max(g_hits, key=lambda record: record["log_e_value"])
    representative = best_d if best_d["log_e_value"] >= best_g["log_e_value"] else best_g
    query_full_length = representative.get("query_full_length", "N/A")
    try:
        fusion_length = int(query_full_length)
    except (TypeError, ValueError):
        fusion_length = best_d["alignment_length"] + best_g["alignment_length"]

    fusion = representative.copy()
    fusion["prediction"] = "vnfDG"
    fusion["alignment_length"] = fusion_length
    fusion["gene_status"] = "Full" if fusion_length >= NIF_GENE_THRESHOLDS["vnfDG"] else "Fragment"
    fusion["log_e_value"] = max(best_d["log_e_value"], best_g["log_e_value"])
    return fusion


def select_best_records(all_records):
    records_by_query = {}
    for record in all_records:
        records_by_query.setdefault(record["query_name"], []).append(record)

    results = {}
    for query, query_records in records_by_query.items():
        fusion_record = build_vnfdg_fusion_record(query_records)
        if fusion_record is not None:
            results[query] = fusion_record
            continue
        for record in query_records:
            current = results.get(query)
            if current is None:
                results[query] = record
                continue
            record_is_fusion = record.get("prediction") in FUSION_GENES
            current_is_fusion = current.get("prediction") in FUSION_GENES
            if record_is_fusion and not current_is_fusion:
                results[query] = record
            elif record_is_fusion == current_is_fusion:
                record_target = record.get("prediction") in TARGET_GENES
                current_target = current.get("prediction") in TARGET_GENES
                if record_target and not current_target:
                    results[query] = record
                elif current_target and not record_target:
                    continue
                elif record_target and current_target and record.get("prediction") != current.get("prediction"):
                    record_distance = record.get("nn_distance")
                    current_distance = current.get("nn_distance")
                    if (
                        isinstance(record_distance, (int, float))
                        and isinstance(current_distance, (int, float))
                    ):
                        if record_distance + 0.05 < current_distance:
                            results[query] = record
                        elif current_distance + 0.05 < record_distance:
                            continue
                        elif record["log_e_value"] > current["log_e_value"]:
                            results[query] = record
                    elif record["log_e_value"] > current["log_e_value"]:
                        results[query] = record
                elif record["log_e_value"] > current["log_e_value"]:
                    results[query] = record
    return results


# ============================================================
# Full_operon検出（ポスト処理）
# 条件1: 同一クエリに2種類以上の異なるnif遺伝子がヒット
# 条件2: それらのヒットのalignment_length >= NIF_GENE_THRESHOLDS[gene]
# 条件3: クエリ長が構成遺伝子の通常full lengthの1.6倍以上
# ============================================================

def build_operon_label(genes):
    ordered_genes = [gene for gene in TARGET_GENES if gene in genes]
    if not ordered_genes:
        return ""
    prefixes = {gene[:3] for gene in ordered_genes}
    if len(prefixes) == 1:
        prefix = ordered_genes[0][:3]
        return prefix + "".join(gene[3:] for gene in ordered_genes)
    return "+".join(ordered_genes)


def is_operon_length(query_full_length, full_genes):
    if query_full_length == "N/A" or not full_genes:
        return False
    try:
        query_full_length = int(query_full_length)
    except (TypeError, ValueError):
        return False
    normal_full_length = max(NIF_GENE_THRESHOLDS[gene] for gene in full_genes)
    return query_full_length >= normal_full_length * OPERON_LENGTH_MULTIPLIER


def detect_operon_queries(all_records):
    """
    all_records（select_best_records前の全レコード）を使って
    Full_operon候補のクエリ名セットと、クエリ→構成nif遺伝子セットのマップを返す。

    戻り値:
      operon_queries   : set of query names
      operon_gene_map  : dict { query_name -> set of nif gene names that qualify as Full }
    """
    query_gene_aln = {}
    query_lengths = {}
    for r in all_records:
        pred = r["prediction"]
        if pred in FUSION_GENES:
            continue
        if pred not in OPERON_CANDIDATE_GENES:
            continue
        qname = r["query_name"]
        query_lengths.setdefault(qname, r.get("query_full_length", "N/A"))
        aln = r["alignment_length"]
        if qname not in query_gene_aln:
            query_gene_aln[qname] = {}
        # 同一遺伝子への複数ヒットは最長アライメントを保持
        if pred not in query_gene_aln[qname] or aln > query_gene_aln[qname][pred]:
            query_gene_aln[qname][pred] = aln

    operon_queries = set()
    operon_gene_map = {}   # { qname -> set of gene names that pass threshold }
    for qname, gene_aln in query_gene_aln.items():
        full_genes = {
            gene for gene, aln in gene_aln.items()
            if aln >= NIF_GENE_THRESHOLDS[gene]
        }
        if len(full_genes) >= 2 and is_operon_length(query_lengths.get(qname, "N/A"), full_genes):
            operon_queries.add(qname)
            operon_gene_map[qname] = full_genes

    return operon_queries, operon_gene_map


def apply_operon_status(best_records, operon_queries, operon_gene_map=None):
    """best_recordsのうちoperon_queriesに含まれるFullをFull_operonに更新"""
    operon_gene_map = operon_gene_map or {}
    for qname, rec in best_records.items():
        rec["operon_label"] = ""
        if qname in operon_queries and rec["gene_status"] == "Full":
            rec["gene_status"] = "Full_operon"
            rec["operon_label"] = build_operon_label(operon_gene_map.get(qname, set()))
    return best_records


# ============================================================
# マトリクス生成・FASTA書き出し
# ============================================================

def get_gene_status_matrix(records):
    status_sets = {gene: set() for gene in TARGET_GENES}
    status_order = ["Full_operon", "Full", "Fragment", "unclassifiable"]
    for r in records.values():
        pred = r["prediction"]
        gene_status = r["gene_status"]
        if pred in status_sets and gene_status != "N/A":
            status_sets[pred].add(gene_status)
    return {
        gene: "+".join(s for s in status_order if s in status_sets[gene]) or "N/A"
        for gene in TARGET_GENES
    }


def write_selected_fasta(fasta_file, records, output_fasta, selected_genes=None, empty_label="target genes"):
    selected_genes = set(selected_genes or TARGET_GENES)
    selected_map = {r["query_name"]: r["prediction"] for r in records.values() if r["prediction"] in selected_genes}
    if not selected_map:
        print(f"No {empty_label} found in {os.path.basename(fasta_file)}. Skipping FASTA output.")
        return
    try:
        with open(output_fasta, "w") as out:
            for rec in SeqIO.parse(fasta_file, "fasta"):
                if rec.id in selected_map:
                    new_id = f"{rec.id}|{selected_map[rec.id]}"
                    rec.id = new_id
                    rec.description = ""
                    SeqIO.write(rec, out, "fasta")
    except FileNotFoundError:
        print(f"Input FASTA file not found: {fasta_file}")
    except Exception as e:
        print(f"Error writing selected FASTA to {output_fasta}: {e}")


# ============================================================
# GenBank local nif-cluster extraction
# ============================================================

def normalize_identifier(value):
    parts = (value or "").strip().lstrip(">").split()
    return parts[0] if parts else ""


def query_identifier_candidates(query):
    first = normalize_identifier(query)
    pieces = [first]
    for separator in ("|", ";", ","):
        if separator in first:
            pieces.extend(part.strip() for part in first.split(separator))
    if "_" in first:
        pieces.append(first.rsplit("_", 1)[0])
    seen = set()
    candidates = []
    for piece in pieces:
        if piece and piece not in seen and piece not in TARGET_GENES:
            candidates.append(piece)
            seen.add(piece)
    return candidates


def query_cds_index_candidate(query):
    first = normalize_identifier(query)
    suffix = first.rsplit("_", 1)[-1]
    if not suffix.isdigit():
        return None
    index = int(suffix)
    return index if index > 0 else None


def qualifier_value(feature, key):
    values = feature.qualifiers.get(key, [])
    if not values:
        return None
    return str(values[0])


def parse_fasta_sequences(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[normalize_identifier(record.id)] = str(record.seq).replace("*", "").upper()
    return sequences


def parse_genbank_cds_features(genbank_file):
    features = []
    records = list(SeqIO.parse(genbank_file, "genbank"))
    if not records:
        raise ValueError("No GenBank records were found.")

    for record in records:
        contig = record.id or record.name or "record"
        contig_length = len(record.seq)
        if contig_length <= 0:
            continue
        cds_index = 0
        for feature in record.features:
            if feature.type != "CDS" or feature.location is None:
                continue
            cds_index += 1
            start = int(feature.location.start)
            end = int(feature.location.end)
            if end <= start:
                continue
            strand = feature.location.strand if feature.location.strand in (-1, 1) else 1
            features.append({
                "contig": contig,
                "contig_length": contig_length,
                "cds_index": cds_index,
                "start": max(0, start),
                "end": min(contig_length, end),
                "strand": strand,
                "protein_id": qualifier_value(feature, "protein_id"),
                "locus_tag": qualifier_value(feature, "locus_tag"),
                "old_locus_tag": qualifier_value(feature, "old_locus_tag"),
                "gene": qualifier_value(feature, "gene"),
                "translation": (qualifier_value(feature, "translation") or "").replace("*", "").upper() or None,
            })
    return features


def match_nif_hits_to_genbank(records, features, fasta_sequences=None):
    by_protein_id = {}
    by_locus_tag = {}
    by_old_locus_tag = {}
    by_gene = {}
    by_cds_index = {}
    by_translation = {}
    for feature in features:
        by_cds_index.setdefault(feature["cds_index"], feature)
        if feature["translation"]:
            by_translation.setdefault(feature["translation"], []).append(feature)
        for mapping, value in (
            (by_protein_id, feature["protein_id"]),
            (by_locus_tag, feature["locus_tag"]),
            (by_old_locus_tag, feature["old_locus_tag"]),
            (by_gene, feature["gene"]),
        ):
            key = normalize_identifier(value)
            if key:
                mapping.setdefault(key, []).append(feature)

    record_list = list(records.values()) if isinstance(records, dict) else list(records)
    matches = []
    used = set()
    prediction_counts = {}
    for record in record_list:
        prediction = str(record.get("prediction") or "")
        if prediction in TARGET_GENES:
            prediction_counts[prediction] = prediction_counts.get(prediction, 0) + 1

    for record in record_list:
        prediction = str(record.get("prediction") or "")
        if prediction not in TARGET_GENES:
            continue
        query = str(record.get("query_name") or record.get("query") or "")
        candidates = query_identifier_candidates(query)
        found = []
        for mapping in (by_protein_id, by_locus_tag, by_old_locus_tag):
            for candidate in candidates:
                found = mapping.get(candidate, [])
                if found:
                    break
            if found:
                break
        if not found:
            exact_query = normalize_identifier(query)
            gene_matches = by_gene.get(exact_query, [])
            if len(gene_matches) == 1:
                found = gene_matches
            elif exact_query == prediction and len(by_gene.get(prediction, [])) == 1:
                found = by_gene[prediction]
            elif prediction_counts.get(prediction) == 1 and len(by_gene.get(prediction, [])) == 1:
                found = by_gene[prediction]
        if not found:
            cds_index = query_cds_index_candidate(query)
            feature = by_cds_index.get(cds_index) if cds_index else None
            if feature and normalize_identifier(feature["gene"]) == prediction:
                found = [feature]
        if not found and fasta_sequences:
            for candidate in candidates:
                sequence = fasta_sequences.get(candidate)
                if not sequence:
                    continue
                translated_matches = [
                    feature
                    for feature in by_translation.get(sequence, [])
                    if normalize_identifier(feature["gene"]) in ("", prediction)
                ]
                if translated_matches:
                    found = translated_matches
                    break

        for feature in found:
            key = (feature["contig"], feature["start"], feature["end"], prediction)
            if key in used:
                continue
            used.add(key)
            matches.append({
                "feature": feature,
                "prediction": prediction,
                "query": query,
                "completeness": str(record.get("gene_status") or record.get("completeness") or ""),
                "operon_label": record.get("operon_label") or record.get("operonLabel") or None,
            })
    return matches


def feature_distance_bp(feature_a, feature_b):
    if feature_a["contig"] != feature_b["contig"]:
        return None
    if feature_a["end"] < feature_b["start"]:
        return feature_b["start"] - feature_a["end"]
    if feature_b["end"] < feature_a["start"]:
        return feature_a["start"] - feature_b["end"]
    return 0


def apply_vnf_context_filter(records, genbank_file, query_file, max_distance_bp=DEFAULT_VNF_CONTEXT_DISTANCE_BP,
                             enabled=True):
    if not enabled or not genbank_file:
        return records
    if not os.path.isfile(genbank_file):
        print(f"GenBank file not found: {genbank_file}. Skipping vnf context filter.")
        return records

    try:
        features = parse_genbank_cds_features(genbank_file)
        fasta_sequences = parse_fasta_sequences(query_file) if query_file and os.path.isfile(query_file) else None
        matches = match_nif_hits_to_genbank(records, features, fasta_sequences)
    except Exception as e:
        print(f"Error applying vnf context filter with {genbank_file}: {e}")
        return records

    if not matches:
        print("No GenBank CDS matches found for vnf context filter. Keeping unfiltered predictions.")
        return records

    anchor_features = [
        match["feature"]
        for match in matches
        if match["prediction"] in VNF_CONTEXT_ANCHOR_GENES
    ]
    feature_by_query_prediction = {
        (match["query"], match["prediction"]): match["feature"]
        for match in matches
    }
    feature_by_query = {}
    for match in matches:
        feature_by_query.setdefault(match["query"], match["feature"])

    filtered = []
    demoted = 0
    checked_vnf_hits = 0
    for record in records:
        prediction = record.get("prediction")
        if prediction not in VNF_CONTEXT_FILTER_GENES:
            filtered.append(record)
            continue
        checked_vnf_hits += 1

        query = str(record.get("query_name") or "")
        feature = feature_by_query_prediction.get((query, prediction)) or feature_by_query.get(query)
        keep_as_vnf = False
        if feature and anchor_features:
            keep_as_vnf = any(
                distance is not None and distance <= max_distance_bp
                for distance in (feature_distance_bp(feature, anchor) for anchor in anchor_features)
            )

        if keep_as_vnf:
            filtered.append(record)
            continue

        demoted_record = record.copy()
        demoted_record["vnf_context_filtered_prediction"] = prediction
        demoted_record["prediction"] = "other"
        demoted_record["gene_status"] = "N/A"
        demoted_record["operon_label"] = ""
        filtered.append(demoted_record)
        demoted += 1

    if demoted:
        print(
            f"vnf context filter: demoted {demoted} vnf hit(s) outside "
            f"{max_distance_bp} bp of vnfG/vnfDG anchors."
        )
    elif anchor_features:
        print(f"vnf context filter: all vnf hits are within {max_distance_bp} bp of vnfG/vnfDG anchors.")
    elif checked_vnf_hits:
        print("vnf context filter: no vnfG/vnfDG anchor found; non-anchor vnf hits were not reported as vnf.")
    return filtered


def group_local_regions(matches, context_padding):
    regions = []
    by_contig = {}
    for match in matches:
        by_contig.setdefault(match["feature"]["contig"], []).append(match)

    for contig, contig_matches in sorted(by_contig.items()):
        ordered = sorted(contig_matches, key=lambda match: match["feature"]["start"])
        current = []
        current_end = -1
        for match in ordered:
            if not current or match["feature"]["start"] - current_end <= LOCAL_CONTEXT_MERGE_DISTANCE:
                current.append(match)
                current_end = max(current_end, match["feature"]["end"])
            else:
                contig_length = current[0]["feature"]["contig_length"]
                start = max(0, min(m["feature"]["start"] for m in current) - context_padding)
                end = min(contig_length, max(m["feature"]["end"] for m in current) + context_padding)
                regions.append((contig, start, end, current))
                current = [match]
                current_end = match["feature"]["end"]
        if current:
            contig_length = current[0]["feature"]["contig_length"]
            start = max(0, min(m["feature"]["start"] for m in current) - context_padding)
            end = min(contig_length, max(m["feature"]["end"] for m in current) + context_padding)
            regions.append((contig, start, end, current))
    return regions


def local_region_label(contig, start, end):
    return f"{contig}: {start + 1:,}-{end:,}"


def estimated_svg_text_width(label, font_size):
    return len(label) * font_size * 0.58


def choose_label_y(x1, x2, y, label, occupied_ranges, font_size=7):
    center = (x1 + x2) / 2
    half_width = max((x2 - x1) / 2, estimated_svg_text_width(label, font_size) / 2)
    candidate_range = (center - half_width - 3, center + half_width + 3)
    for level in range(4):
        ranges = occupied_ranges.setdefault(level, [])
        if all(candidate_range[1] < start or candidate_range[0] > end for start, end in ranges):
            ranges.append(candidate_range)
            return y - 15 - level * 12
    occupied_ranges.setdefault(3, []).append(candidate_range)
    return y - 15 - 3 * 12


def svg_arrow(x1, x2, y, strand, color, label="", label_font_size=7, min_label_width=10, label_y=None):
    x1, x2 = sorted((x1, x2))
    if x2 - x1 < 2:
        x2 = x1 + 2
    height = 18
    head = min(10, max(4, (x2 - x1) * 0.35))
    if strand == -1:
        points = [
            (x1, y),
            (x1 + head, y - height / 2),
            (x2, y - height / 2),
            (x2, y + height / 2),
            (x1 + head, y + height / 2),
        ]
    else:
        points = [
            (x1, y - height / 2),
            (x2 - head, y - height / 2),
            (x2, y),
            (x2 - head, y + height / 2),
            (x1, y + height / 2),
        ]
    point_text = " ".join(f"{x:.1f},{py:.1f}" for x, py in points)
    parts = [f'<polygon points="{point_text}" fill="{color}" stroke="{color}" stroke-width="0.6" />']
    force_label = label in {"T"}
    if label and (force_label or x2 - x1 >= min_label_width):
        text_y = y - 15 if label_y is None else label_y
        parts.append(
            f'<text x="{(x1 + x2) / 2:.1f}" y="{text_y:.1f}" text-anchor="middle" '
            f'font-size="{label_font_size}" font-family="Arial">{escape(label)}</text>'
        )
    return "\n".join(parts)


def write_svg(path, width, height, body):
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<rect width="100%" height="100%" fill="white"/>\n'
        f"{body}\n"
        "</svg>\n"
    )
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(svg)
    print(f"Genomic context SVG saved to: {path}")


def unique_target_legend_entries(matches):
    predictions = {match["prediction"] for match in matches if match["prediction"] in TARGET_GENES}
    return [
        (gene, NIF_GENE_COLORS.get(gene, "#0f766e"))
        for gene in TARGET_GENES
        if gene in predictions
    ]


def display_gene_label(label):
    return {
        "nifT": "T",
        "cnfR/patB_like": "cnfR_like",
    }.get(label, label)


def svg_legend(entries, x, y, width, columns=4):
    if not entries:
        return "", 0
    item_width = max(170, (width - x - 40) / max(1, columns))
    row_height = 24
    rows = math.ceil(len(entries) / columns)
    parts = [
        f'<text x="{x}" y="{y}" font-size="13" font-family="Arial" '
        f'font-weight="bold">Detected target genes</text>'
    ]
    for idx, (label, color) in enumerate(entries):
        col = idx % columns
        row = idx // columns
        item_x = x + col * item_width
        item_y = y + 20 + row * row_height
        parts.append(
            f'<rect x="{item_x:.1f}" y="{item_y - 10:.1f}" width="18" height="10" '
            f'fill="{color}" stroke="{color}" stroke-width="0.6" />'
        )
        parts.append(
            f'<text x="{item_x + 25:.1f}" y="{item_y:.1f}" font-size="12" '
            f'font-family="Arial">{escape(display_gene_label(label))}</text>'
        )
    return "\n".join(parts), 25 + rows * row_height


def build_overview_svg(matches, output_svg):
    contigs = sorted({m["feature"]["contig"]: m["feature"]["contig_length"] for m in matches}.items())
    if not contigs:
        return False
    width = 1600
    margin_left = 160
    margin_right = 70
    row_height = 90
    height = max(130, 45 + row_height * len(contigs))
    body = []
    for idx, (contig, contig_length) in enumerate(contigs):
        y = 60 + idx * row_height
        x_start = margin_left
        x_end = width - margin_right
        body.append(f'<text x="20" y="{y + 5}" font-size="14" font-family="Arial">{escape(contig)}</text>')
        body.append(f'<line x1="{x_start}" x2="{x_end}" y1="{y}" y2="{y}" stroke="#20242a" stroke-width="1"/>')
        body.append(f'<text x="{x_start}" y="{y + 34}" font-size="11" font-family="Arial">1</text>')
        body.append(
            f'<text x="{x_end}" y="{y + 34}" text-anchor="end" font-size="11" '
            f'font-family="Arial">{contig_length:,}</text>'
        )
        scale = (x_end - x_start) / max(1, contig_length)
        for match in sorted((m for m in matches if m["feature"]["contig"] == contig), key=lambda m: m["feature"]["start"]):
            feature = match["feature"]
            x1 = x_start + feature["start"] * scale
            x2 = x_start + feature["end"] * scale
            label = display_gene_label(match["prediction"][3:] if match["prediction"].startswith(("nif", "vnf")) else match["prediction"])
            color = NIF_GENE_COLORS.get(match["prediction"], "#0f766e")
            body.append(svg_arrow(x1, x2, y, feature["strand"], color, label))
    write_svg(output_svg, width, height, "\n".join(body))
    return True


def build_local_context_svg(matches, features, context_padding, output_svg):
    regions = group_local_regions(matches, context_padding)
    if not regions:
        return False
    hit_by_feature = {
        (match["feature"]["contig"], match["feature"]["start"], match["feature"]["end"]): match
        for match in matches
    }
    width = 1800
    margin_left = 210
    margin_right = 45
    row_height = 110
    legend_entries = unique_target_legend_entries(matches)
    legend_height = 0
    if legend_entries:
        legend_height = 55 + math.ceil(len(legend_entries) / 4) * 24
    height = max(150, 50 + row_height * len(regions) + legend_height)
    body = []
    for idx, (contig, start, end, _region_matches) in enumerate(regions):
        y = 70 + idx * row_height
        x_start = margin_left
        x_end = width - margin_right
        span = max(1, end - start)
        body.append(
            f'<text x="20" y="{y + 5}" font-size="13" font-family="Arial">'
            f'{escape(local_region_label(contig, start, end))}</text>'
        )
        body.append(f'<line x1="{x_start}" x2="{x_end}" y1="{y}" y2="{y}" stroke="#20242a" stroke-width="0.8"/>')
        overlapping = [
            feature for feature in features
            if feature["contig"] == contig and feature["end"] >= start and feature["start"] <= end
        ]
        occupied_label_ranges = {}
        for feature in sorted(overlapping, key=lambda item: item["start"]):
            clipped_start = max(feature["start"], start)
            clipped_end = min(feature["end"], end)
            x1 = x_start + (clipped_start - start) / span * (x_end - x_start)
            x2 = x_start + (clipped_end - start) / span * (x_end - x_start)
            hit = hit_by_feature.get((feature["contig"], feature["start"], feature["end"]))
            if hit:
                label = hit["operon_label"] if hit["completeness"] == "Full_operon" and hit["operon_label"] else display_gene_label(hit["prediction"])
                color = NIF_GENE_COLORS.get(hit["prediction"], "#0f766e")
            else:
                label = ""
                color = "#c7cdd4"
            label_y = choose_label_y(x1, x2, y, label, occupied_label_ranges) if label else None
            body.append(svg_arrow(x1, x2, y, feature["strand"], color, label, label_y=label_y))
            if hit and hit["completeness"]:
                suffix = {"Full": "(full)", "Fragment": "(frag.)", "Full_operon": "(operon)"}.get(hit["completeness"], "")
                if suffix:
                    body.append(
                        f'<text x="{(x1 + x2) / 2:.1f}" y="{y + 27}" text-anchor="middle" '
                        f'font-size="10" font-family="Arial" fill="{color}">{suffix}</text>'
                    )
    legend, _legend_height = svg_legend(legend_entries, margin_left, 45 + row_height * len(regions), width)
    if legend:
        body.append(legend)
    write_svg(output_svg, width, height, "\n".join(body))
    return True


def build_local_context_genbank(genbank_file, matches, context_padding, output_gbk,
                                output_label="Nif-encoding"):
    records = list(SeqIO.parse(genbank_file, "genbank"))
    records_by_contig = {
        (record.id or record.name or "record"): record
        for record in records
    }
    hit_by_feature = {
        (match["feature"]["contig"], match["feature"]["start"], match["feature"]["end"]): match
        for match in matches
    }
    region_records = []
    regions = group_local_regions(matches, context_padding)

    for index, (contig, start, end, region_matches) in enumerate(regions, start=1):
        record = records_by_contig.get(contig)
        if record is None:
            continue
        sub_record = record[start:end]
        region_id = f"{contig}_{start + 1}_{end}"
        sub_record.id = region_id
        sub_record.name = region_id[:16]
        sub_record.description = f"Nif-Finder local context region {index}: {contig}:{start + 1}-{end}"
        sub_record.annotations["molecule_type"] = record.annotations.get("molecule_type", "DNA")
        matched_labels = ", ".join(
            f"{match['prediction']}:{match['query']}:{match['completeness']}"
            for match in sorted(region_matches, key=lambda item: item["feature"]["start"])
        )
        comment = f"Generated by Nif-Finder from {contig}:{start + 1}-{end}. Matched hits: {matched_labels}."
        if sub_record.annotations.get("comment"):
            comment = f"{sub_record.annotations['comment']}\n{comment}"
        sub_record.annotations["comment"] = comment

        for feature in sub_record.features:
            if feature.type != "CDS" or feature.location is None:
                continue
            original_start = int(feature.location.start) + start
            original_end = int(feature.location.end) + start
            hit = hit_by_feature.get((contig, original_start, original_end))
            if not hit:
                continue
            feature.qualifiers["nif_finder_gene"] = [hit["prediction"]]
            feature.qualifiers["nif_finder_status"] = [hit["completeness"]]
            feature.qualifiers["nif_finder_query"] = [hit["query"]]
            if hit["operon_label"]:
                feature.qualifiers["nif_finder_operon"] = [hit["operon_label"]]
        region_records.append(sub_record)

    if not region_records:
        return False
    SeqIO.write(region_records, output_gbk, "genbank")
    print(f"{output_label} GenBank region saved to: {output_gbk}")
    return True


def select_vnf_anchor_region_matches(matches, max_distance_bp):
    anchor_matches = [
        match for match in matches
        if match["prediction"] in VNF_CONTEXT_ANCHOR_GENES
    ]
    if not anchor_matches:
        return []

    selected = []
    used = set()
    for match in matches:
        if match["prediction"] not in VNF_REGION_GBK_GENES:
            continue
        if match["prediction"] in VNF_CONTEXT_ANCHOR_GENES:
            keep = True
        else:
            keep = any(
                distance is not None and distance <= max_distance_bp
                for distance in (
                    feature_distance_bp(match["feature"], anchor["feature"])
                    for anchor in anchor_matches
                )
            )
        if not keep:
            continue
        key = (match["feature"]["contig"], match["feature"]["start"], match["feature"]["end"], match["prediction"])
        if key in used:
            continue
        used.add(key)
        selected.append(match)
    return selected


def write_nif_cluster_genbank(genbank_file, records, query_file, output_prefix, context_size_kb):
    if not genbank_file:
        return
    if not os.path.isfile(genbank_file):
        print(f"GenBank file not found: {genbank_file}. Skipping nif-cluster GenBank output.")
        return
    try:
        features = parse_genbank_cds_features(genbank_file)
        fasta_sequences = parse_fasta_sequences(query_file) if query_file and os.path.isfile(query_file) else None
        matches = match_nif_hits_to_genbank(records, features, fasta_sequences)
        if not matches:
            print(
                "A GenBank file was provided, but no Nif-Finder hits could be matched to CDS "
                "protein_id, locus_tag, old_locus_tag, a unique gene qualifier, or CDS translation."
            )
            return
        context_padding = context_size_kb * 1000
        output_gbk = f"{output_prefix}_nif_encoding_region.gbk"
        build_local_context_genbank(genbank_file, matches, context_padding, output_gbk)
        build_overview_svg(matches, f"{output_prefix}_genome_overview.svg")
        build_local_context_svg(matches, features, context_padding, f"{output_prefix}_local_context.svg")
    except Exception as e:
        print(f"Error extracting nif-cluster GenBank region from {genbank_file}: {e}")


def write_vnf_cluster_genbank(genbank_file, records, query_file, output_prefix, context_size_kb,
                              vnf_context_distance_bp=DEFAULT_VNF_CONTEXT_DISTANCE_BP):
    if not genbank_file:
        return
    if not os.path.isfile(genbank_file):
        print(f"GenBank file not found: {genbank_file}. Skipping vnf-cluster GenBank output.")
        return
    try:
        features = parse_genbank_cds_features(genbank_file)
        fasta_sequences = parse_fasta_sequences(query_file) if query_file and os.path.isfile(query_file) else None
        matches = match_nif_hits_to_genbank(records, features, fasta_sequences)
        if not matches:
            print(
                "A GenBank file was provided, but no Nif-Finder vnf hits could be matched to CDS "
                "protein_id, locus_tag, old_locus_tag, a unique gene qualifier, or CDS translation."
            )
            return
        vnf_matches = select_vnf_anchor_region_matches(matches, vnf_context_distance_bp)
        if not vnf_matches:
            print("No vnfG or vnfDG anchor-containing vnf region found. Skipping vnf GenBank output.")
            return
        context_padding = context_size_kb * 1000
        output_gbk = f"{output_prefix}_vnf_encoding_region.gbk"
        build_local_context_genbank(
            genbank_file, vnf_matches, context_padding, output_gbk,
            output_label="Vnf-encoding",
        )
    except Exception as e:
        print(f"Error extracting vnf-cluster GenBank region from {genbank_file}: {e}")


def find_matching_genbank(genbank_dir, base_name):
    if not genbank_dir:
        return None
    for ext in (".gbk", ".gb", ".gbff", ".genbank"):
        candidate = os.path.join(genbank_dir, f"{base_name}{ext}")
        if os.path.isfile(candidate):
            return candidate
    return None


# ============================================================
# 散布図出力
# target遺伝子ごとのパネル
# 参照データ: nif該当（遺伝子色・薄め）/ other（灰色）
# クエリヒット: 遺伝子色（濃い）・ステータス別マーカー
# ============================================================

def _lighten_hex(hex_color, factor=0.55):
    """HEXカラーをfactorの割合で白に近づけて返す（参照nif色用）"""
    hex_color = hex_color.lstrip("#")
    r, g, b = [int(hex_color[i:i+2], 16) for i in (0, 2, 4)]
    r = int(r + (255 - r) * factor)
    g = int(g + (255 - g) * factor)
    b = int(b + (255 - b) * factor)
    return f"#{r:02X}{g:02X}{b:02X}"


def plot_scatter(all_records, references_loaded, output_png,
                 operon_queries=None, operon_gene_map=None, panel_genes=None):
    """
    target遺伝子ごとの散布図を出力する。
    X軸: クエリのタンパク質長 (query_full_length)
    Y軸: -log10(E-value)

    Full_operonクエリの表示ルール（問題1修正）:
      operon_gene_map { qname -> set of nif genes } を使い、
      Full_operonクエリはそのオペロンを構成する各nif遺伝子のパネルに★でプロット。
      all_records内の対応レコード（prediction == 構成遺伝子）を使って座標を決定。
      他の（構成遺伝子ではない）パネルへのFragment/Full重複表示はスキップ。

    unclassifiableの表示ルール（問題2・3修正）:
      prediction == "unclassifiable" または gene_status == "unclassifiable" のレコードは
      全パネルに◇でまとめて表示（パラログ誤ヒット等を視覚化）。
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.lines as mlines
    except ImportError:
        print("matplotlib is not installed. Skipping scatter plot.")
        print("Install with: pip install matplotlib")
        return

    if operon_queries is None:
        operon_queries = set()
    if operon_gene_map is None:
        operon_gene_map = {}
    panel_genes = list(panel_genes or PANEL_GENES)

    ncols = 3
    nrows = math.ceil(len(panel_genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows))
    axes = axes.flatten()

    # ------------------------------------------------------------------
    # all_records をパネルごとに振り分け
    #
    # Full_operonクエリの処理:
    #   operon_gene_map[qname] に含まれる遺伝子パネルにのみ★でプロット。
    #   all_records から prediction == 構成遺伝子 のレコードを探し座標に使う。
    #   構成遺伝子以外のパネルへの重複ヒットはスキップ。
    #
    # non-operon通常クエリ: predictionパネルにそのままプロット。
    #
    # unclassifiable: 全パネル共通の◇リストに収集。
    # ------------------------------------------------------------------

    # gene_hits[gene] = list of (record, display_status)
    # display_status: "Full" / "Fragment" / "Full_operon" / "N/A"
    gene_hits = {gene: [] for gene in panel_genes}
    unclassifiable_hits = []

    # Full_operonクエリの代表レコードを遺伝子ごとに収集
    # operon_records[qname][gene] = best record for that gene (highest log_e_value)
    operon_records = {}
    for r in all_records:
        qname = r["query_name"]
        pred  = r["prediction"]
        if qname not in operon_queries:
            continue
        if pred not in operon_gene_map.get(qname, set()):
            continue
        if qname not in operon_records:
            operon_records[qname] = {}
        prev = operon_records[qname].get(pred)
        if prev is None or r["log_e_value"] > prev["log_e_value"]:
            operon_records[qname][pred] = r

    # Full_operon★を構成遺伝子パネルに追加
    for qname, gene_rec_map in operon_records.items():
        for gene, rec in gene_rec_map.items():
            panel_gene = GENE_PANEL_MAP.get(gene, gene)
            if panel_gene in gene_hits:
                gene_hits[panel_gene].append((rec, "Full_operon"))

    # 通常クエリ（非operon）の振り分け
    for r in all_records:
        pred   = r["prediction"]
        qname  = r["query_name"]
        status = r["gene_status"]

        # unclassifiable → 全パネル共通リスト
        if pred == "unclassifiable" or status == "unclassifiable":
            unclassifiable_hits.append(r)
            continue

        # Full_operonクエリのレコードは上で処理済みなのでスキップ
        if qname in operon_queries:
            continue

        panel_gene = GENE_PANEL_MAP.get(pred, pred)
        if panel_gene not in gene_hits:
            continue

        gene_hits[panel_gene].append((r, status))

    # ------------------------------------------------------------------
    # 参照データをパネルごとに振り分け
    # ------------------------------------------------------------------
    ref_by_gene = {gene: {"target": [], "other": []} for gene in panel_genes}
    for i, (ref_data, _scale) in enumerate(references_loaded):
        if i < len(panel_genes):
            panel_gene = panel_genes[i]
        else:
            continue
        for ref in ref_data:
            attr = ref["attribute"]
            if panel_gene == "vnfG" and attr in ("vnfG", "vnfD", "vnfDG"):
                ref_by_gene[panel_gene]["target"].append(ref)
            elif GENE_PANEL_MAP.get(attr, attr) == panel_gene and attr in NIF_GENE_COLORS:
                ref_by_gene[panel_gene]["target"].append(ref)
            else:
                ref_by_gene[panel_gene]["other"].append(ref)

    # ステータス別マーカー設定: (marker, size, edgecolor, linewidth)
    status_markers = {
        "Full":            ("o",  70, "black", 1.2),
        "Fragment":        ("^",  70, "black", 1.2),
        "Full_operon":     ("*", 130, "black", 1.2),
        "N/A":             ("x",  40, "black", 1.0),
        "unclassifiable":  ("D",  60, "gray",  1.0),
    }

    for idx, gene in enumerate(panel_genes):
        ax = axes[idx]
        color       = NIF_GENE_COLORS[gene]
        color_light = _lighten_hex(color, factor=0.35 if gene == "nifX" else 0.60)
        color_other = "#CCCCCC"

        refs = ref_by_gene[gene]

        # 参照 other（灰色・最背面）
        if refs["other"]:
            colored_other = {}
            gray_other = []
            for ref in refs["other"]:
                attr = ref.get("attribute")
                panel_attr = GENE_PANEL_MAP.get(attr, attr)
                if panel_attr in panel_genes and panel_attr != gene and panel_attr in NIF_GENE_COLORS:
                    colored_other.setdefault(panel_attr, []).append(ref)
                else:
                    gray_other.append(ref)
            if gray_other:
                ox = [r["x"] for r in gray_other]
                oy = [r["y"] for r in gray_other]
                ax.scatter(ox, oy, c=color_other, s=14, alpha=0.40,
                           zorder=1, label="Ref: other")
            for attr_gene, attr_refs in colored_other.items():
                ox = [r["x"] for r in attr_refs]
                oy = [r["y"] for r in attr_refs]
                attr_color = _lighten_hex(NIF_GENE_COLORS.get(attr_gene, color_other), factor=0.45)
                ax.scatter(ox, oy, c=attr_color, s=16, alpha=0.55,
                           zorder=2, label=f"Ref: {attr_gene}")

        # 参照 target（同系色・薄め）
        if refs["target"]:
            nx = [r["x"] for r in refs["target"]]
            ny = [r["y"] for r in refs["target"]]
            ax.scatter(nx, ny, c=color_light, s=20, alpha=0.75 if gene == "nifX" else 0.60,
                       zorder=2, label=f"Ref: {gene}/target")

        # クエリヒット（X軸: タンパク質長）
        hits_by_status = {}
        for (rec, disp_status) in gene_hits[gene]:
            hits_by_status.setdefault((disp_status, rec.get("prediction", gene)), []).append(rec)

        for (status, prediction), hits in hits_by_status.items():
            marker, size, ecol, elw = status_markers.get(status, ("o", 70, "black", 1.2))
            hx = []
            for h in hits:
                qfl = h["query_full_length"]
                hx.append(qfl if isinstance(qfl, (int, float)) else 0)
            hy = [h["log_e_value"] for h in hits]
            hit_color = NIF_GENE_COLORS.get(prediction, color)
            hit_label = f"Query: {status}" if prediction == gene else f"Query: {prediction} {status}"
            ax.scatter(hx, hy, c=hit_color, s=size, marker=marker,
                       alpha=1.0, zorder=4, edgecolors=ecol, linewidths=elw,
                       label=hit_label)

        # 判定不能ヒット（◇・灰色）を全パネルに表示
        if unclassifiable_hits:
            ux = []
            for h in unclassifiable_hits:
                qfl = h["query_full_length"]
                ux.append(qfl if isinstance(qfl, (int, float)) else 0)
            uy = [h["log_e_value"] for h in unclassifiable_hits]
            ax.scatter(ux, uy, c="#AAAAAA", s=60, marker="D",
                       alpha=0.6, zorder=3, edgecolors="gray", linewidths=0.8,
                       label="Query: unclassifiable")

        # Full/Fragment閾値の垂直破線（参照データスケール基準の目安）
        if gene in NIF_GENE_THRESHOLDS:
            thr = NIF_GENE_THRESHOLDS[gene]
            ax.axvline(x=thr, color=color, linestyle="--",
                       linewidth=1.2, alpha=0.7,
                       label=f"Aln threshold ({thr} aa)")
        if gene == "vnfG" and "vnfDG" in NIF_GENE_THRESHOLDS:
            thr = NIF_GENE_THRESHOLDS["vnfDG"]
            ax.axvline(x=thr, color=NIF_GENE_COLORS["vnfDG"], linestyle=":",
                       linewidth=1.2, alpha=0.7,
                       label=f"vnfDG threshold ({thr} aa)")

        ax.set_title(gene, fontsize=13, fontweight="bold", color=color)
        ax.set_xlabel("Protein Length (aa)", fontsize=9)
        ax.set_ylabel("-log10(E-value)", fontsize=9)
        if gene in PANEL_X_LIMITS:
            ax.set_xlim(right=PANEL_X_LIMITS[gene])
        ax.tick_params(labelsize=8)

        # パネルごと凡例（重複除去）
        handles, labels = ax.get_legend_handles_labels()
        seen = {}
        for h, l in zip(handles, labels):
            if l not in seen:
                seen[l] = h
        ax.legend(seen.values(), seen.keys(), fontsize=7, loc="lower right",
                  framealpha=0.35, markerscale=0.9)

    for ax in axes[len(panel_genes):]:
        ax.axis("off")

    # 図全体の共通凡例（下部）
    rep_color       = NIF_GENE_COLORS["nifH"]
    rep_color_light = _lighten_hex(rep_color, factor=0.55)
    legend_elements = [
        mlines.Line2D([0], [0], marker='o', color='w', markerfacecolor='#CCCCCC',
                      markeredgecolor='#AAAAAA', markersize=7, label='Ref: other'),
        mlines.Line2D([0], [0], marker='o', color='w', markerfacecolor=rep_color_light,
                      markeredgecolor=rep_color_light, markersize=7, label='Ref: target'),
        mlines.Line2D([0], [0], marker='o', color='w', markerfacecolor=rep_color,
                      markersize=7, label='Query: Full'),
        mlines.Line2D([0], [0], marker='^', color='w', markerfacecolor=rep_color,
                      markersize=7, label='Query: Fragment'),
        mlines.Line2D([0], [0], marker='*', color='w', markerfacecolor=rep_color,
                      markersize=11, label='Query: Full_operon'),
        mlines.Line2D([0], [0], marker='D', color='w', markerfacecolor='#AAAAAA',
                      markeredgecolor='gray', markersize=7, label='Query: unclassifiable'),
        mlines.Line2D([0], [0], color='gray', linestyle='--',
                      linewidth=1.2, label='Aln threshold (ref scale)'),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=7,
               fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, 0.01))

    plt.suptitle("Nif/Vnf Gene HMMscan Results", fontsize=15, fontweight="bold", y=1.01)
    plt.tight_layout(rect=[0, 0.06, 1, 1])

    try:
        plt.savefig(output_png, dpi=150, bbox_inches="tight")
        print(f"Scatter plot saved to: {output_png}")
    except Exception as e:
        print(f"Error saving scatter plot: {e}")
    plt.close(fig)


# ============================================================
# 単一クエリ処理
# ============================================================

def process_single_query(query_file, profile_files, reference_files, output_prefix,
                         save_fasta, save_vnf_hdgken_vupabc_fasta, cpu, plot, jobs=DEFAULT_JOBS, e_value=DEFAULT_E_VALUE,
                         genbank_file=None, context_size_kb=DEFAULT_CONTEXT_SIZE_KB,
                         vnf_context_filter=True, vnf_context_distance_bp=DEFAULT_VNF_CONTEXT_DISTANCE_BP,
                         save_vnf_region_gbk=False):
    base_name = os.path.splitext(os.path.basename(query_file))[0]
    output_prefix = output_prefix or f"{base_name}_results"
    query_lengths = get_query_lengths(query_file)
    references_loaded = [load_reference_data(r_file) for r_file in reference_files]
    active_panel_genes = panel_genes_for_profiles(profile_files)

    all_records = run_profile_scans(query_file, profile_files, references_loaded,
                                    query_lengths, cpu, jobs=jobs, e_value=e_value)
    all_records = apply_vnf_context_filter(
        all_records, genbank_file, query_file,
        max_distance_bp=vnf_context_distance_bp,
        enabled=vnf_context_filter,
    )

    operon_queries, operon_gene_map = detect_operon_queries(all_records)
    unique_records = select_best_records(all_records)
    unique_records = apply_operon_status(unique_records, operon_queries, operon_gene_map)

    output_summary_file = f"{output_prefix}.txt"
    try:
        with open(output_summary_file, "w") as f:
            f.write("Query\t-log_Evalue\tAlign_Len\tQuery_Length\tPrediction\tCompleteness\tOperon\n")
            for rec in unique_records.values():
                f.write(f"{rec['query_name']}\t{rec['log_e_value']:.2f}\t"
                        f"{rec['alignment_length']}\t{rec['query_full_length']}\t"
                        f"{rec['prediction']}\t{rec['gene_status']}\t"
                        f"{rec.get('operon_label', '')}\n")
    except IOError as e:
        print(f"Error writing summary output to {output_summary_file}: {e}")

    if save_fasta:
        fasta_output_file = f"{output_prefix}_nif_vnf_related.faa"
        write_selected_fasta(query_file, unique_records, fasta_output_file)
    if save_vnf_hdgken_vupabc_fasta:
        fasta_output_file = f"{output_prefix}_vnfHDGKEN_vupABC.faa"
        write_selected_fasta(query_file, unique_records, fasta_output_file,
                             selected_genes=VNF_HDGKEN_VUPABC_GENES, empty_label="vnfHDGKEN/vupABC genes")

    if genbank_file:
        write_nif_cluster_genbank(genbank_file, unique_records, query_file, output_prefix, context_size_kb)
        if save_vnf_region_gbk:
            write_vnf_cluster_genbank(
                genbank_file, unique_records, query_file, output_prefix, context_size_kb,
                vnf_context_distance_bp=vnf_context_distance_bp,
            )

    if plot:
        plot_png = f"{output_prefix}_scatter.png"
        plot_scatter(all_records, references_loaded, plot_png,
                     operon_queries=operon_queries, operon_gene_map=operon_gene_map,
                     panel_genes=active_panel_genes)


# ============================================================
# ディレクトリ処理
# ============================================================

def process_query_directory(query_dir, profile_files, reference_files, matrix_output_file,
                            save_fasta, save_vnf_hdgken_vupabc_fasta, cpu, plot, jobs=DEFAULT_JOBS, e_value=DEFAULT_E_VALUE,
                            genbank_dir=None, context_size_kb=DEFAULT_CONTEXT_SIZE_KB,
                            vnf_context_filter=True, vnf_context_distance_bp=DEFAULT_VNF_CONTEXT_DISTANCE_BP,
                            save_vnf_region_gbk=False):
    faa_files = sorted(glob.glob(os.path.join(query_dir, "*.faa")))
    if not faa_files:
        print(f"No .faa files found in directory: {query_dir}")
        return

    references_loaded = [load_reference_data(r_file) for r_file in reference_files]
    active_panel_genes = panel_genes_for_profiles(profile_files)

    try:
        with open(matrix_output_file, "w") as out:
            out.write("Genome\t" + "\t".join(TARGET_GENES) + "\n")
            for faa in faa_files:
                genome_name = os.path.basename(faa)
                base = os.path.splitext(genome_name)[0]
                query_lengths = get_query_lengths(faa)

                all_records = run_profile_scans(faa, profile_files, references_loaded,
                                                query_lengths, cpu, jobs=jobs, e_value=e_value)
                genbank_file = find_matching_genbank(genbank_dir, base) if genbank_dir else None
                all_records = apply_vnf_context_filter(
                    all_records, genbank_file, faa,
                    max_distance_bp=vnf_context_distance_bp,
                    enabled=vnf_context_filter,
                )

                operon_queries, operon_gene_map = detect_operon_queries(all_records)
                best = select_best_records(all_records)
                best = apply_operon_status(best, operon_queries, operon_gene_map)

                row = get_gene_status_matrix(best)
                out.write(f"{genome_name}\t" + "\t".join(row[g] for g in TARGET_GENES) + "\n")

                if save_fasta:
                    fasta_output_file = os.path.join(query_dir, f"{base}_nif_vnf_related.faa")
                    write_selected_fasta(faa, best, fasta_output_file)
                if save_vnf_hdgken_vupabc_fasta:
                    fasta_output_file = os.path.join(query_dir, f"{base}_vnfHDGKEN_vupABC.faa")
                    write_selected_fasta(faa, best, fasta_output_file,
                                         selected_genes=VNF_HDGKEN_VUPABC_GENES, empty_label="vnfHDGKEN/vupABC genes")

                if genbank_dir and genbank_file:
                    gbk_prefix = os.path.join(query_dir, f"{base}_results")
                    write_nif_cluster_genbank(genbank_file, best, faa, gbk_prefix, context_size_kb)
                    if save_vnf_region_gbk:
                        write_vnf_cluster_genbank(
                            genbank_file, best, faa, gbk_prefix, context_size_kb,
                            vnf_context_distance_bp=vnf_context_distance_bp,
                        )
                elif genbank_dir:
                    print(f"No matching GenBank file found for {base} in {genbank_dir}.")

                if plot:
                    plot_png = os.path.join(query_dir, f"{base}_scatter.png")
                    plot_scatter(all_records, references_loaded, plot_png,
                                 operon_queries=operon_queries, operon_gene_map=operon_gene_map,
                                 panel_genes=active_panel_genes)

    except IOError as e:
        print(f"Error writing matrix output to {matrix_output_file}: {e}")


# ============================================================
# 6フレーム翻訳（-g モード用）
# ============================================================

def translate_six_frames(genome_fasta, output_faa, min_orf_len=DEFAULT_MIN_ORF_LEN):
    """
    ゲノムDNA FASTAを6フレーム翻訳してタンパク質FASTAに変換する。
    停止コドン（*）でORFを分割し、min_orf_len(aa)以上のものを保持。

    query_nameフォーマット:
        {contig_id}_frame{+1/+2/+3/-1/-2/-3}_start{nt_pos}
    nt_posはフォワード鎖上の0-based塩基座標。
    """
    written = 0
    try:
        with open(output_faa, "w") as out:
            for rec in SeqIO.parse(genome_fasta, "fasta"):
                seq_len = len(rec.seq)
                for strand in (+1, -1):
                    nuc = rec.seq if strand == 1 else rec.seq.reverse_complement()
                    for frame in range(3):
                        trans = str(nuc[frame:].translate())
                        parts = trans.split("*")
                        # 各ORFの開始位置（フォワード鎖0-based）を追跡
                        orf_start_in_frame = frame  # フレーム内のnt位置（正鎖基準）
                        for part in parts:
                            if len(part) >= min_orf_len:
                                if strand == 1:
                                    nt_pos = orf_start_in_frame
                                else:
                                    # リバース鎖の場合は正鎖座標に変換
                                    nt_pos = seq_len - orf_start_in_frame - len(part) * 3
                                frame_label = f"+{frame + 1}" if strand == 1 else f"-{frame + 1}"
                                orf_id = f"{rec.id}_frame{frame_label}_start{nt_pos}"
                                out.write(f">{orf_id}\n{part}\n")
                                written += 1
                            orf_start_in_frame += (len(part) + 1) * 3
    except FileNotFoundError:
        print(f"Genome FASTA file not found: {genome_fasta}")
        raise
    except Exception as e:
        print(f"Error during 6-frame translation of {genome_fasta}: {e}")
        raise
    print(f"6-frame translation: {written} ORFs (>= {min_orf_len} aa) written to {output_faa}")
    return written


# ============================================================
# ゲノムDNAクエリ処理（-g モード）
# ============================================================

def process_genome_query(genome_file, profile_files, reference_files, output_prefix,
                         save_fasta, save_vnf_hdgken_vupabc_fasta, cpu, plot, min_orf_len, jobs=DEFAULT_JOBS,
                         e_value=DEFAULT_E_VALUE, genbank_file=None,
                         context_size_kb=DEFAULT_CONTEXT_SIZE_KB,
                         vnf_context_filter=True, vnf_context_distance_bp=DEFAULT_VNF_CONTEXT_DISTANCE_BP,
                         save_vnf_region_gbk=False):
    """
    ゲノムDNA FASTAを6フレーム翻訳してからhmmscanに渡す。
    翻訳済み一時FASTAを生成後、process_single_queryと同じフローで処理する。
    """
    base_name = os.path.splitext(os.path.basename(genome_file))[0]
    output_prefix = output_prefix or f"{base_name}_genome_results"

    # 6フレーム翻訳 → 一時FASTAファイル
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_6frame.faa') as tmp_faa:
        tmp_faa_path = tmp_faa.name

    try:
        n_orfs = translate_six_frames(genome_file, tmp_faa_path, min_orf_len=min_orf_len)
        if n_orfs == 0:
            print(f"No ORFs found in {genome_file} with min_orf_len={min_orf_len}. Aborting.")
            return

        # 翻訳済みFASTAを通常のタンパク質クエリとして処理
        query_lengths = get_query_lengths(tmp_faa_path)
        references_loaded = [load_reference_data(r_file) for r_file in reference_files]
        active_panel_genes = panel_genes_for_profiles(profile_files)

        all_records = run_profile_scans(tmp_faa_path, profile_files, references_loaded,
                                        query_lengths, cpu, jobs=jobs, e_value=e_value)
        all_records = apply_vnf_context_filter(
            all_records, genbank_file, tmp_faa_path,
            max_distance_bp=vnf_context_distance_bp,
            enabled=vnf_context_filter,
        )

        operon_queries, operon_gene_map = detect_operon_queries(all_records)
        unique_records = select_best_records(all_records)
        unique_records = apply_operon_status(unique_records, operon_queries, operon_gene_map)

        output_summary_file = f"{output_prefix}.txt"
        try:
            with open(output_summary_file, "w") as f:
                f.write("Query\t-log_Evalue\tAlign_Len\tQuery_Length\tPrediction\tCompleteness\tOperon\n")
                for rec in unique_records.values():
                    f.write(f"{rec['query_name']}\t{rec['log_e_value']:.2f}\t"
                            f"{rec['alignment_length']}\t{rec['query_full_length']}\t"
                            f"{rec['prediction']}\t{rec['gene_status']}\t"
                            f"{rec.get('operon_label', '')}\n")
        except IOError as e:
            print(f"Error writing summary output to {output_summary_file}: {e}")

        if save_fasta:
            # -s の場合は翻訳済み配列（ORF）を保存
            fasta_output_file = f"{output_prefix}_nif_vnf_related.faa"
            write_selected_fasta(tmp_faa_path, unique_records, fasta_output_file)
        if save_vnf_hdgken_vupabc_fasta:
            fasta_output_file = f"{output_prefix}_vnfHDGKEN_vupABC.faa"
            write_selected_fasta(tmp_faa_path, unique_records, fasta_output_file,
                                 selected_genes=VNF_HDGKEN_VUPABC_GENES, empty_label="vnfHDGKEN/vupABC genes")

        if genbank_file:
            write_nif_cluster_genbank(genbank_file, unique_records, tmp_faa_path, output_prefix, context_size_kb)
            if save_vnf_region_gbk:
                write_vnf_cluster_genbank(
                    genbank_file, unique_records, tmp_faa_path, output_prefix, context_size_kb,
                    vnf_context_distance_bp=vnf_context_distance_bp,
                )

        if plot:
            plot_png = f"{output_prefix}_scatter.png"
            plot_scatter(all_records, references_loaded, plot_png,
                         operon_queries=operon_queries, operon_gene_map=operon_gene_map,
                         panel_genes=active_panel_genes)

    finally:
        if os.path.exists(tmp_faa_path):
            os.remove(tmp_faa_path)


# ============================================================
# エントリーポイント
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Nif Finder: Identifies nif, vnf, and related genes in protein sequences using HMMscan."
    )
    parser.add_argument("-q", "--query",
                        help="Path to a single protein FASTA file.")
    parser.add_argument("-d", "--query_dir",
                        help="Path to a directory containing multiple .faa files.")
    parser.add_argument("-g", "--genome",
                        help="Path to a genome DNA FASTA file. "
                             "Internally performs 6-frame translation before running hmmscan. "
                             "Useful for detecting nif genes interrupted by intervening sequences.")
    parser.add_argument("-t", "--profile", nargs='+',
                        help="Paths to HMM profile files (e.g., nifH.hmm nifD.hmm). "
                             "If omitted, standard nif models are loaded from NIF_FINDER_DB.")
    parser.add_argument("-r", "--reference", nargs='+',
                        help="Paths to reference files (e.g., nifH_ref.tsv nifD_ref.tsv). "
                             "Must correspond to --profile files. "
                             "If omitted, standard references are loaded from NIF_FINDER_DB.")
    parser.add_argument("-o", "--outprefix",
                        help="Prefix for output files. Default: input file name.")
    parser.add_argument("-m", "--matrix_output", default="nif_matrix.tsv",
                        help="Output file for gene status matrix (directory mode). "
                             "Default: nif_matrix.tsv")
    parser.add_argument("--genbank",
                        help="Path to a GenBank file for extracting the nif/vnf-related local region "
                             "as an annotated .gbk file (single -q or -g mode).")
    parser.add_argument("--genbank_dir",
                        help="Directory containing GenBank files matching .faa basenames "
                             "for directory mode. Supported extensions: .gbk, .gb, .gbff, .genbank.")
    parser.add_argument("--context_size_kb", type=int, default=DEFAULT_CONTEXT_SIZE_KB,
                        help="Size of the flanking region added to detected nif hits when extracting "
                             f"GenBank local context. Range: 1-{MAX_CONTEXT_SIZE_KB} kb. "
                             f"Default: {DEFAULT_CONTEXT_SIZE_KB} kb.")
    parser.add_argument("--vnf_context_filter", dest="vnf_context_filter", action="store_true", default=True,
                        help="Filter non-anchor vnf calls using GenBank coordinates. When vnfG or vnfDG "
                             "is detected, only vnf hits within the configured distance are reported as vnf. "
                             "Default: on.")
    parser.add_argument("--no_vnf_context_filter", dest="vnf_context_filter", action="store_false",
                        help="Disable the GenBank-coordinate vnf context filter.")
    parser.add_argument("--vnf_context_distance_bp", type=int, default=DEFAULT_VNF_CONTEXT_DISTANCE_BP,
                        help="Maximum distance from vnfG/vnfDG anchors for reporting non-anchor vnf hits. "
                             f"Default: {DEFAULT_VNF_CONTEXT_DISTANCE_BP} bp.")
    parser.add_argument("-s", "--save_fasta", action="store_true",
                        help="Save predicted target nif, vnf, and related sequences to FASTA.")
    parser.add_argument("--save_vnf_hdgken_vupabc_fasta", action="store_true",
                        help="Run exclusive vnf/vup mode and save predicted vnfHDGKEN/vupABC sequences "
                             "to FASTA (vnfH/nifH, vnfD, vnfDG, vnfG, vnfK, vnfE/nifE, "
                             "vnfN/nifN, vupA/modA, vupB/modB, vupC/modC).")
    parser.add_argument("--save_vnf_region_gbk", action="store_true",
                        help="Save a separate GenBank file for vnf anchor regions only. "
                             "A region is written only when vnfG or vnfDG is detected; "
                             "vup/mod-only regions are not exported.")
    parser.add_argument("--save_vnf_hdk_en_fasta", action="store_true",
                        help=argparse.SUPPRESS)
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Save scatter plot PNG (alignment_length vs -log10(Evalue)). "
                             "Target panels by gene; reference data in gray, "
                             "query hits color-coded; Full/Fragment threshold shown.")
    parser.add_argument("-c", "--cpu", type=int, default=DEFAULT_CPU,
                        help=f"Number of threads for HMMscan. Default: {DEFAULT_CPU}")
    parser.add_argument("-j", "--jobs", type=int, default=DEFAULT_JOBS,
                        help="Number of HMM profile scans to run in parallel. "
                             f"Default: {DEFAULT_JOBS} (sequential, previous behavior).")
    parser.add_argument("-e", "--evalue", type=float, default=DEFAULT_E_VALUE,
                        help=f"HMMscan E-value threshold. Default: {DEFAULT_E_VALUE:g}. "
                             "Changing this can alter sensitivity and specificity.")
    parser.add_argument("--min_orf_len", type=int, default=DEFAULT_MIN_ORF_LEN,
                        help=f"Minimum ORF length (aa) retained after 6-frame translation "
                             f"(-g mode only). Default: {DEFAULT_MIN_ORF_LEN}")
    args = parser.parse_args()

    if args.cpu < 1:
        parser.error("--cpu must be 1 or greater.")
    if args.jobs < 1:
        parser.error("--jobs must be 1 or greater.")
    if args.min_orf_len < 1:
        parser.error("--min_orf_len must be 1 or greater.")
    if args.evalue <= 0:
        parser.error("--evalue must be greater than 0.")
    if args.context_size_kb < 1 or args.context_size_kb > MAX_CONTEXT_SIZE_KB:
        parser.error(f"--context_size_kb must be between 1 and {MAX_CONTEXT_SIZE_KB}.")
    if args.vnf_context_distance_bp < 0:
        parser.error("--vnf_context_distance_bp must be 0 or greater.")
    if args.genbank and args.query_dir:
        parser.error("--genbank is for single -q or -g mode. Use --genbank_dir with -d.")
    if args.genbank_dir and not args.query_dir:
        parser.error("--genbank_dir can only be used with -d/--query_dir.")
    if args.save_vnf_region_gbk and not (args.genbank or args.genbank_dir):
        parser.error("--save_vnf_region_gbk requires --genbank for -q/-g or --genbank_dir for -d.")
    save_vnf_hdgken_vupabc_fasta = args.save_vnf_hdgken_vupabc_fasta or args.save_vnf_hdk_en_fasta

    modes = [args.query, args.query_dir, args.genome]
    if sum(bool(m) for m in modes) > 1:
        parser.error("Please specify only one of -q, -d, or -g.")

    if args.query:
        profile_files, reference_files = resolve_model_paths(args.profile, args.reference, parser)
        if save_vnf_hdgken_vupabc_fasta:
            profile_files, reference_files = apply_vnf_exclusive_model_filter(profile_files, reference_files, parser)
        process_single_query(args.query, profile_files, reference_files,
                             args.outprefix, args.save_fasta, save_vnf_hdgken_vupabc_fasta, args.cpu, args.plot,
                             jobs=args.jobs, e_value=args.evalue,
                             genbank_file=args.genbank, context_size_kb=args.context_size_kb,
                             vnf_context_filter=args.vnf_context_filter,
                             vnf_context_distance_bp=args.vnf_context_distance_bp,
                             save_vnf_region_gbk=args.save_vnf_region_gbk)
    elif args.query_dir:
        profile_files, reference_files = resolve_model_paths(args.profile, args.reference, parser)
        if save_vnf_hdgken_vupabc_fasta:
            profile_files, reference_files = apply_vnf_exclusive_model_filter(profile_files, reference_files, parser)
        process_query_directory(args.query_dir, profile_files, reference_files,
                                args.matrix_output, args.save_fasta, save_vnf_hdgken_vupabc_fasta, args.cpu, args.plot,
                                jobs=args.jobs, e_value=args.evalue,
                                genbank_dir=args.genbank_dir, context_size_kb=args.context_size_kb,
                                vnf_context_filter=args.vnf_context_filter,
                                vnf_context_distance_bp=args.vnf_context_distance_bp,
                                save_vnf_region_gbk=args.save_vnf_region_gbk)
    elif args.genome:
        profile_files, reference_files = resolve_model_paths(args.profile, args.reference, parser)
        if save_vnf_hdgken_vupabc_fasta:
            profile_files, reference_files = apply_vnf_exclusive_model_filter(profile_files, reference_files, parser)
        process_genome_query(args.genome, profile_files, reference_files,
                             args.outprefix, args.save_fasta, save_vnf_hdgken_vupabc_fasta, args.cpu, args.plot,
                             args.min_orf_len, jobs=args.jobs, e_value=args.evalue,
                             genbank_file=args.genbank, context_size_kb=args.context_size_kb,
                             vnf_context_filter=args.vnf_context_filter,
                             vnf_context_distance_bp=args.vnf_context_distance_bp,
                             save_vnf_region_gbk=args.save_vnf_region_gbk)
    else:
        parser.print_help()
        print("\nError: Please specify one of -q, -d, or -g.")


if __name__ == "__main__":
    main()
