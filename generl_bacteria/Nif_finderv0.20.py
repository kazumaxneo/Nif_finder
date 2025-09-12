
from Bio import SeqIO
import subprocess
import math
import csv
import os
from math import sqrt
import argparse
import glob
import tempfile
#v0.20 修正: Full_operon判定におけるquery_full_lengthの初期化位置を修正
NIF_GENES = ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB"]

# 定数の定義
DEFAULT_E_VALUE = 1e-10
DEFAULT_CPU = 8
DEFAULT_LOG_E_VALUE_MIN = 250 # E_valueが0の場合の-log10(Evalue)のデフォルト値

# Nif遺伝子ごとの閾値（Full/Fragment判定用）
NIF_GENE_THRESHOLDS = {
    "nifH": 240,
    "nifD": 370,
    "nifK": 410,
    "nifE": 400,
    "nifN": 400,
    "nifB": 370
}

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

def load_reference_data(reference_file):
    reference_data = []
    try:
        with open(reference_file, 'r') as ref_file:
            reader = csv.DictReader(ref_file, delimiter='\t')
            for row in reader:
                try:
                    alignment_length = int(row["length"])
                    e_value = float(row["Evalue"])
                    log_e_value = -math.log10(e_value) if e_value > 0 else DEFAULT_LOG_E_VALUE_MIN
                    attribute = row["attribute"]
                    reference_data.append({"x": alignment_length, "y": log_e_value, "attribute": attribute})
                except (ValueError, KeyError) as e:
                    print(f"Skipping malformed row in {reference_file}: {row}. Error: {e}")
    except FileNotFoundError:
        print(f"Reference file not found: {reference_file}")
        raise
    return reference_data

def find_nearest_attribute(alignment_length, log_e_value, reference_data):
    min_distance = float('inf')
    nearest_attribute = "unknown"
    for ref in reference_data:
        ref_x = ref["x"]
        ref_y = ref["y"]
        distance = sqrt((alignment_length - ref_x) ** 2 + (log_e_value - ref_y) ** 2)
        if distance < min_distance:
            min_distance = distance
            nearest_attribute = ref["attribute"]
    return nearest_attribute

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

def convert_to_single_tab(input_file, reference_data, query_lengths):
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
                    nearest_attribute = find_nearest_attribute(alignment_length, log_e_value, reference_data)
                    gene_status = "N/A"
                    if nearest_attribute in NIF_GENE_THRESHOLDS:
                        threshold = NIF_GENE_THRESHOLDS[nearest_attribute]
                        if alignment_length >= threshold:
                            if query_full_length != "N/A" and query_full_length >= 1.8 * alignment_length:
                                gene_status = "Full_operon"
                            else:
                                gene_status = "Full"
                        else:
                            gene_status = "Fragment"
                    record = {
                        "query_name": query_name,
                        "log_e_value": log_e_value,
                        "alignment_length": alignment_length,
                        "query_full_length": query_full_length,
                        "prediction": nearest_attribute,
                        "gene_status": gene_status
                    }
                    records.append(record)
                except (ValueError, IndexError) as e:
                    print(f"Skipping malformed HMMscan output line: {line.strip()}. Error: {e}")
    except FileNotFoundError:
        print(f"HMMscan output file not found: {input_file}")
        raise
    return records

def select_best_records(all_records):
    results = {}
    for record in all_records:
        query = record["query_name"]
        if query not in results or record["log_e_value"] > results[query]["log_e_value"]:
            results[query] = record
    return results

def get_gene_status_matrix(records):
    status = {gene: "N/A" for gene in NIF_GENES}
    for r in records.values():
        pred = r["prediction"]
        if pred in status:
            status[pred] = r["gene_status"]
    return status

def write_selected_fasta(fasta_file, records, output_fasta):
    selected_map = {r["query_name"]: r["prediction"] for r in records.values() if r["prediction"] in NIF_GENES}
    if not selected_map:
        print(f"No Nif genes found in {os.path.basename(fasta_file)}. Skipping FASTA output.")
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

def process_single_query(query_file, profile_files, reference_files, output_prefix, save_fasta, cpu):
    base_name = os.path.splitext(os.path.basename(query_file))[0]
    output_prefix = output_prefix or f"{base_name}_results"
    query_lengths = get_query_lengths(query_file)
    all_records = []
    references_loaded = [load_reference_data(r_file) for r_file in reference_files]
    for i, (target_file, reference_data) in enumerate(zip(profile_files, references_loaded)):
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as temp_output_file:
            temp_output = temp_output_file.name
        try:
            run_hmmscan(query_file, target_file, temp_output, cpu=cpu)
            all_records.extend(convert_to_single_tab(temp_output, reference_data, query_lengths))
        except Exception as e:
            print(f"Processing of {query_file} with {target_file} failed: {e}")
        finally:
            if os.path.exists(temp_output):
                os.remove(temp_output)
    unique_records = select_best_records(all_records)
    output_summary_file = f"{output_prefix}.txt"
    try:
        with open(output_summary_file, "w") as f:
            f.write("Query	-log_Evalue	Align_Len	Query_Length	Prediction	Completeness\n")
            for rec in unique_records.values():
                f.write(f"{rec['query_name']}\t{rec['log_e_value']:.2f}\t{rec['alignment_length']}\t{rec['query_full_length']}\t{rec['prediction']}\t{rec['gene_status']}\n")
    except IOError as e:
        print(f"Error writing summary output to {output_summary_file}: {e}")
    if save_fasta:
        fasta_output_file = f"{base_name}_nifHDKENB.faa"
        write_selected_fasta(query_file, unique_records, fasta_output_file)

def process_query_directory(query_dir, profile_files, reference_files, matrix_output_file, save_fasta, cpu):
    faa_files = sorted(glob.glob(os.path.join(query_dir, "*.faa")))
    if not faa_files:
        print(f"No .faa files found in directory: {query_dir}")
        return
    references_loaded = [load_reference_data(r_file) for r_file in reference_files]
    try:
        with open(matrix_output_file, "w") as out:
            out.write("Genome\t" + "\t".join(NIF_GENES) + "\n")
            for faa in faa_files:
                genome_name = os.path.basename(faa)
                base = os.path.splitext(genome_name)[0]
                query_lengths = get_query_lengths(faa)
                all_records = []
                for i, (tfile, reference_data) in enumerate(zip(profile_files, references_loaded)):
                    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.txt') as tmp_out_file:
                        tmp_out = tmp_out_file.name
                    try:
                        run_hmmscan(faa, tfile, tmp_out, cpu=cpu)
                        all_records.extend(convert_to_single_tab(tmp_out, reference_data, query_lengths))
                    except Exception as e:
                        print(f"Processing of {faa} with {tfile} failed: {e}")
                    finally:
                        if os.path.exists(tmp_out):
                            os.remove(tmp_out)
                best = select_best_records(all_records)
                row = get_gene_status_matrix(best)
                out.write(f"{genome_name}\t" + "\t".join(row[g] for g in NIF_GENES) + "\n")
                if save_fasta:
                    fasta_output_file = os.path.join(query_dir, f"{base}_nifHDKENB.faa")
                    write_selected_fasta(faa, best, fasta_output_file)
    except IOError as e:
        print(f"Error writing matrix output to {matrix_output_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Nif Finder: Identifies Nif genes in protein sequences using HMMscan.")
    parser.add_argument("-q", "--query", help="Path to a single protein FASTA file.")
    parser.add_argument("-d", "--query_dir", help="Path to a directory containing multiple .faa files.")
    parser.add_argument("-t", "--profile", required=True, nargs='+', help="Paths to HMM profile files (e.g., nifH.hmm nifD.hmm).")
    parser.add_argument("-r", "--reference", required=True, nargs='+', help="Paths to reference files (e.g., nifH_ref.tsv nifD_ref.tsv). Must correspond to --profile files.")
    parser.add_argument("-o", "--outprefix", help="Prefix for single query output files. Default is based on query file name.")
    parser.add_argument("-m", "--matrix_output", default="nif_matrix.tsv", help="Output file name for the gene status matrix when processing a directory. Default: nif_matrix.tsv")
    parser.add_argument("-s", "--save_fasta", action="store_true", help="Save predicted NifHDKENB sequences to a new FASTA file for each input file.")
    parser.add_argument("-c", "--cpu", type=int, default=DEFAULT_CPU, help=f"Number of threads for HMMscan. Default: {DEFAULT_CPU}")
    args = parser.parse_args()
    if len(args.profile) != len(args.reference):
        parser.error("The number of --profile files must match the number of --reference files.")
    if args.query and args.query_dir:
        parser.error("Please specify either -q (single query file) or -d (query directory), not both.")
    elif args.query:
        process_single_query(args.query, args.profile, args.reference, args.outprefix, args.save_fasta, args.cpu)
    elif args.query_dir:
        process_query_directory(args.query_dir, args.profile, args.reference, args.matrix_output, args.save_fasta, args.cpu)
    else:
        parser.print_help()
        print("\nError: Please specify either -q (single query file) or -d (query directory).")

if __name__ == "__main__":
    main()
