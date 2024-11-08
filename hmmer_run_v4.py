import pandas as pd
import numpy as np
import argparse
import subprocess
import os
import matplotlib.pyplot as plt

# コマンドライン引数の設定
def parse_arguments():
    parser = argparse.ArgumentParser(description="Run HMMER, format output, find closest point in TSV data, and plot scatter plot.")
    parser.add_argument("--input", required=True, help="Input .faa file for HMMER scan")
    parser.add_argument("--hmm", required=True, help="HMM profile file")
    parser.add_argument("--output", required=True, help="Output file for hmmscan")
    parser.add_argument("--cpu", type=int, default=8, help="Number of CPU cores to use (default: 8)")
    parser.add_argument("--evalue", type=str, default="1e-10", help="E-value threshold (default: 1e-10)")
    parser.add_argument("--tsv", required=True, help="TSV file with length, log_Evalue, and attribute columns")
    parser.add_argument("--savefile", required=True, help="File to save the closest attributes")
    parser.add_argument("--plotfile", required=True, help="PDF file to save the scatter plot")  # 保存するPDFファイルを指定
    return parser.parse_args()

# HMMER実行とフォーマット処理
def run_hmmscan_and_format(args):
    # Step 1: HMMERの実行
    hmmscan_cmd = [
        "hmmscan", 
        "--domtblout", args.output, 
        "--cpu", str(args.cpu), 
        "-E", args.evalue, 
        args.hmm, 
        args.input
    ]
    
    print("Running HMMER scan...")
    subprocess.run(hmmscan_cmd, check=True)
    
    # Step 2: フォーマット処理（コメント除去とTSVフォーマット）
    formatted_output = "formatted_" + os.path.basename(args.output)
    grep_sed_cmd = f"grep -v '#' {args.output} | tr -s ' ' '\\t' > {formatted_output}"
    
    print("Filtering spaces and formatting output...")
    subprocess.run(grep_sed_cmd, shell=True, check=True)

    print(f"Formatted output saved to {formatted_output}")
    return formatted_output

# log_EvalueとAlignment_Lengthの計算
def calculate_log_evalue_and_alignment_length(formatted_output):
    df = pd.read_csv(formatted_output, sep="\t", header=None)
    
    # 5列目のEvalueの-log10を計算し、新しい列を追加
    df['log_Evalue'] = -np.log10(df[5].astype(float))
    
    # アライメント長を計算し、新しい列を追加 (to - from)
    df['Alignment_Length'] = df[18] - df[17]
    
    return df

# 最も近い点を見つけて視覚化
def find_closest_points_and_plot(df, tsv_file, savefile, plotfile):
    # TSVデータの読み込み
    tsv_df = pd.read_csv(tsv_file, sep="\t")

    # 結果をファイルに保存（追記モード）
    with open(savefile, 'w') as f:
        for index, row in df.iterrows():
            target_evalue = row['log_Evalue']
            target_length = row['Alignment_Length']

            # ユークリッド距離を計算
            tsv_df['distance'] = np.sqrt((tsv_df['Evalue'] - target_evalue)**2 + (tsv_df['length'] - target_length)**2)

            # 最も近い点のインデックスを取得
            min_index = tsv_df['distance'].idxmin()

            # 最も近い点の属性を取得
            closest_point = tsv_df.loc[min_index]

            # 結果をテキストファイルに追記
            f.write(f"Row {index + 1}: Closest point -> Attribute={closest_point['attribute']}\n")
            print(f"Row {index + 1}: Closest point saved to {savefile}")

            # 散布図のプロット (元の点と最も近い点を色分け)
            plt.scatter(tsv_df['Evalue'], tsv_df['length'], label="Data Points", color='blue')
            plt.scatter(target_evalue, target_length, color='red', label="Target Point")
            plt.scatter(closest_point['Evalue'], closest_point['length'], color='green', label="Closest Point")
            plt.xlabel("log_Evalue")
            plt.ylabel("Alignment_Length")
            plt.legend()

            # プロットをPDF形式で保存
            pdf_filename = f"{plotfile}_row_{index + 1}.pdf"
            plt.savefig(pdf_filename)
            print(f"Scatter plot saved to {pdf_filename}")
            plt.clf()  # 次のプロットのためにクリア

# メイン処理
def main():
    # コマンドライン引数の取得
    args = parse_arguments()

    # HMMERの実行とフォーマット処理
    formatted_output = run_hmmscan_and_format(args)

    # log_EvalueとAlignment_Lengthの計算
    df = calculate_log_evalue_and_alignment_length(formatted_output)

    # 最も近い点を見つけ、属性をテキストファイルに保存し、散布図をPDFとして保存
    find_closest_points_and_plot(df, args.tsv, args.savefile, args.plotfile)

if __name__ == "__main__":
    main()
