import collections
import argparse

# コマンドライン引数の設定
def parse_arguments():
    parser = argparse.ArgumentParser(description="Count attributes from combined hits file.")
    parser.add_argument("--input", required=True, help="Input file containing combined hits (e.g., combined_hits.txt)")
    parser.add_argument("--countfile", required=True, help="Output file for attribute counts (e.g., attribute_counts.txt)")
    return parser.parse_args()

# Attributeを集計する関数
def count_attributes(input_file, count_file):
    # 入力ファイルの内容を読み込む
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Attribute（3列目）の出現回数をカウント
    attribute_counts = collections.Counter()
    for line in lines:
        # Attributeが "-> Attribute=" 以降に記載されているので、それを取得
        attribute = line.split('-> Attribute=')[-1].strip()
        attribute_counts[attribute] += 1

    # 集計結果をファイルに保存
    with open(count_file, 'w') as f:
        for attribute, count in attribute_counts.items():
            f.write(f"{attribute}: {count}\n")

    print(f"Attribute counts saved to {count_file}")

# メイン処理
if __name__ == "__main__":
    # コマンドライン引数を解析
    args = parse_arguments()

    # 属性の出現回数を集計して保存
    count_attributes(args.input, args.countfile)
