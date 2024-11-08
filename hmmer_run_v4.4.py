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
    parser.add_argument("--hmm", required=True, nargs='+', help="HMM profile files (multiple profiles can be specified)")
    parser.add_argument("--output", required=True, help="Output file for hmmscan (temporary file, will be overwritten)")
    parser.add_argument("--cpu", type=int, default=8, help="Number of CPU cores to use (default: 8)")
    parser.add_argument("--evalue", type=str, default="1e-10", help="E-value threshold (default: 1e-10)")
    parser.add_argument("--tsv", required=True, nargs='+', help="TSV files with length, log_Evalue, and attribute columns (multiple TSVs)")
    parser.add_argument("--savefile", required=True, help="File to save the closest attributes (one file for all results)")  # まとめる出力ファイル
    parser.add_argument("--plotdir", required=True, help="Directory to save the scatter plot images")  # サブディレクトリ指定
    return parser.parse_args()

# HMMER実行とフォーマット処理
def run_hmmscan_and_format(hmm_profile, input_file, output_file, cpu, evalue):
    # Step 1: HMMERの実行（標準出力非表示）
    hmmscan_cmd = [
        "hmmscan", 
        "--domtblout", output_file, 
        "--cpu", str(cpu), 
        "-E", evalue, 
        hmm_profile, 
        input_file
    ]
    
    print(f"Running HMMER scan for {hmm_profile}...")

    # stdoutとstderrをDEVNULLに設定してログを抑制
    with open(os.devnull, 'w') as devnull:
        subprocess.run(hmmscan_cmd, check=True, stdout=devnull, stderr=devnull)
    
    # Step 2: フォーマット処理（コメント除去とTSVフォーマット）
    formatted_output = "formatted_" + os.path.basename(output_file)
    grep_sed_cmd = f"grep -v '#' {output_file} | tr -s ' ' '\\t' > {formatted_output}"
    
    print("Filtering spaces and formatting output...")
    subprocess.run(grep_sed_cmd, shell=True, check=True)

    print(f"Formatted output saved to {formatted_output}")
    return formatted_output

# log_EvalueとAlignment_Lengthの計算
def calculate_log_evalue_and_alignment_length(formatted_output):
    df = pd.read_csv(formatted_output, sep="\t", header=None)
    
    # 7列目のEvalueの-log10を計算し、新しい列を追加
    df['log_Evalue'] = -np.log10(df[6].astype(float))  # 6番目のインデックスが正しい7列目
    
    # アライメント長を計算し、新しい列を追加 (to - from)
    df['Alignment_Length'] = df[18] - df[17]
    
    return df

# 最も近い点を見つけて視覚化
def find_closest_points_and_plot(df, tsv_file, savefile, plotdir, profile_name):
    # TSVデータの読み込み
    tsv_df = pd.read_csv(tsv_file, sep="\t")

    # サブディレクトリが存在しない場合、作成する
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    # 結果をファイルに追記モードで保存
    with open(savefile, 'a') as f:  # 'a'モードで追記
        for index, row in df.iterrows():
            # target point の Evalue を 7列目から取得して -log10 変換
            target_evalue_raw = float(row[6])  # 7列目から正しくEvalueを取得
            if target_evalue_raw <= 1e-306:
                target_evalue = 306
            else:
                target_evalue = -np.log10(target_evalue_raw)

            target_length = row['Alignment_Length']

            # 変換した target_evalue をプリント
            print(f"Row {index + 1}: Raw Evalue = {target_evalue_raw}, -log10(Evalue) = {target_evalue}")

            # ユークリッド距離を計算
            tsv_df['distance'] = np.sqrt((tsv_df['Evalue'] - target_evalue)**2 + (tsv_df['length'] - target_length)**2)

            # 最も近い点のインデックスを取得
            min_index = tsv_df['distance'].idxmin()

            # 最も近い点の属性を取得
            closest_point = tsv_df.loc[min_index]

            # 結果をテキストファイルに追記
            f.write(f"{profile_name} - Row {index + 1}: Closest point -> Attribute={closest_point['attribute']}\n")
            print(f"Row {index + 1}: Closest point saved to {savefile}")

            # 散布図のプロット (元の点と最も近い点を色分け)
            plt.figure(figsize=(12, 9))  # サイズを3倍に設定
            plt.scatter(tsv_df['Evalue'], tsv_df['length'], label="Data Points", color='blue')
            plt.scatter(target_evalue, target_length, color='red', label="Target Point")  # target point の -log変換後の値
            plt.scatter(closest_point['Evalue'], closest_point['length'], color='green', label="Closest Point")
            plt.xlabel("log_Evalue")
            plt.ylabel("Alignment_Length")
            plt.legend()

            # プロットをPNG形式で保存
            png_filename = os.path.join(plotdir, f"{profile_name}_row_{index + 1}.png")
            plt.savefig(png_filename)
            print(f"Scatter plot saved to {png_filename}")
            plt.clf()  # 次のプロットのためにクリア

# メイン処理
def main():
    # コマンドライン引数の取得
    args = parse_arguments()

    # 出力ファイルが存在していれば削除（新しく上書きするため）
    if os.path.exists(args.savefile):
        os.remove(args.savefile)

    # 複数のHMMプロファイルに対して処理をループで実行
    for hmm_profile, tsv_file in zip(args.hmm, args.tsv):
        profile_name = os.path.basename(hmm_profile).split("/")[0]  # プロファイル名の取得

        # HMMERの実行とフォーマット処理
        formatted_output = run_hmmscan_and_format(hmm_profile, args.input, args.output, args.cpu, args.evalue)

        # log_EvalueとAlignment_Lengthの計算
        df = calculate_log_evalue_and_alignment_length(formatted_output)

        # 最も近い点を見つけ、属性をテキストファイルに保存し、散布図をPNGとして保存
        find_closest_points_and_plot(df, tsv_file, args.savefile, args.plotdir, profile_name)

if __name__ == "__main__":
    main()


