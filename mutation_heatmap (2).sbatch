#!/usr/bin/env bash
#SBATCH --job-name=mutation_profile
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/%u/slurm_logs/%x_%A.out
#SBATCH --error=/scratch/%u/slurm_logs/%x_%A.err

set -euo pipefail
shopt -s nullglob

SCRATCH_BASE="/scratch/$USER"
mkdir -p "$SCRATCH_BASE/slurm_logs"

FINAL_BASE="${FINAL_BASE:-$HOME/bucket_lists}"
BUCKET_NAME="${BUCKET_NAME:-}"
PDB_STEM="${PDB_STEM:-}"
OUT_BASE="${OUT_BASE:-$FINAL_BASE/mutation_profiles}"
MODEL_TAGS="${MODEL_TAGS:-protein ligand soluble}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --final-base) FINAL_BASE="${2:-}"; shift 2 ;;
    --bucket) BUCKET_NAME="${2:-}"; shift 2 ;;
    --pdb-stem) PDB_STEM="${2:-}"; shift 2 ;;
    --out-base) OUT_BASE="${2:-}"; shift 2 ;;
    --models) MODEL_TAGS="${2:-}"; shift 2 ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

[[ -n "$BUCKET_NAME" ]] || { echo "ERROR: missing --bucket or BUCKET_NAME" >&2; exit 2; }
[[ -n "$PDB_STEM" ]] || { echo "ERROR: missing --pdb-stem or PDB_STEM" >&2; exit 2; }

OUT_DIR="${OUT_BASE}/${BUCKET_NAME}/${PDB_STEM}"
mkdir -p "$OUT_DIR"

read -r -a MODELS <<< "$MODEL_TAGS"
MODEL_DIRS=()
for model in "${MODELS[@]}"; do
  model_dir="${FINAL_BASE}/${model}/${BUCKET_NAME}/${PDB_STEM}"
  wt_fasta="${model_dir}/seqs/seq_0000_wt.fasta"
  metrics_csv="${model_dir}/mpnn_metrics.csv"

  [[ -d "$model_dir" ]] || { echo "ERROR: missing model dir: $model_dir" >&2; exit 3; }
  [[ -f "$wt_fasta" ]] || { echo "ERROR: missing WT fasta: $wt_fasta" >&2; exit 3; }
  [[ -f "$metrics_csv" ]] || { echo "ERROR: missing metrics CSV: $metrics_csv" >&2; exit 3; }

  MODEL_DIRS+=("$model_dir")
done

python3 - "$OUT_DIR" "${MODEL_DIRS[@]}" <<'PY'
import csv
import glob
import os
import sys
from collections import OrderedDict

out_dir = sys.argv[1]
model_dirs = sys.argv[2:]

classes = OrderedDict([
    ("hydrophobic", set("AVLIM")),
    ("hydrophilic", set("STNQCG")),
    ("aromatic", set("FWY")),
    ("charged", set("KRHDE")),
    ("special", set("P")),
])
class_names = list(classes.keys())


def aa_class(aa):
    aa = aa.upper()
    for name, letters in classes.items():
        if aa in letters:
            return name
    return None


def read_single_seq(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    s = "".join(seq).strip().upper()
    if not s:
        raise RuntimeError(f"Empty FASTA sequence: {path}")
    return s


def write_matrix_csv(path, matrix, row_labels, col_labels):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + col_labels)
        for r in row_labels:
            w.writerow([r] + [matrix[r][c] for c in col_labels])


def write_float_matrix_csv(path, matrix, row_labels, col_labels):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([""] + col_labels)
        for r in row_labels:
            w.writerow([r] + [f"{matrix[r][c]:.6f}" for c in col_labels])


def write_text_summary(path, model_tag, counts, matrix, row_norm):
    total_positions = sum(counts["source_total"].values())
    no_change_total = sum(counts["no_change"].values())
    no_change_pct = (100.0 * no_change_total / total_positions) if total_positions else 0.0
    with open(path, "w") as f:
        f.write(f"model={model_tag}\n")
        f.write(f"sequences_processed={counts['sequences_processed']}\n")
        f.write(f"total_positions={total_positions}\n")
        f.write(f"no_change_total={no_change_total}\n")
        f.write(f"no_change_percent={no_change_pct:.2f}\n")
        f.write(f"skipped_non_target={counts['skipped_non_target']}\n")
        f.write(f"length_mismatch_positions={counts['length_mismatch_positions']}\n")
        f.write("\nraw_matrix\n")
        for r in class_names:
            f.write(r + "," + ",".join(str(matrix[r][c]) for c in class_names) + "\n")
        f.write("\nrow_normalized\n")
        for r in class_names:
            f.write(r + "," + ",".join(f"{row_norm[r][c]:.6f}" for c in class_names) + "\n")


def init_counts():
    return {
        "matrix": {r: {c: 0 for c in class_names} for r in class_names},
        "source_total": {r: 0 for r in class_names},
        "no_change": {r: 0 for r in class_names},
        "skipped_non_target": 0,
        "length_mismatch_positions": 0,
        "sequences_processed": 0,
    }


summary_rows = []

for model_dir in model_dirs:
    # model_dir looks like .../protein/helix_rich/107m_assembly1
    # we want model_tag = protein, ligand, soluble
    parts = model_dir.rstrip("/").split(os.sep)
    model_tag = parts[-3] if len(parts) >= 3 else os.path.basename(parts[0])

    wt_fasta = os.path.join(model_dir, "seqs", "seq_0000_wt.fasta")
    wt_seq = read_single_seq(wt_fasta)

    seq_dir = os.path.join(model_dir, "seqs")
    design_fastas = sorted(
        p for p in glob.glob(os.path.join(seq_dir, "seq_*.fasta"))
        if not os.path.basename(p).startswith("seq_0000_wt")
    )

    counts = init_counts()

    for fasta_path in design_fastas:
        mut_seq = read_single_seq(fasta_path)
        n = min(len(wt_seq), len(mut_seq))
        if len(wt_seq) != len(mut_seq):
            counts["length_mismatch_positions"] += abs(len(wt_seq) - len(mut_seq))

        counts["sequences_processed"] += 1

        for i in range(n):
            wt = wt_seq[i]
            mut = mut_seq[i]
            src = aa_class(wt)
            dst = aa_class(mut)

            if src is None or dst is None:
                counts["skipped_non_target"] += 1
                continue

            counts["matrix"][src][dst] += 1
            counts["source_total"][src] += 1
            if wt == mut:
                counts["no_change"][src] += 1

    row_norm = {r: {c: 0.0 for c in class_names} for r in class_names}
    for r in class_names:
        total = counts["source_total"][r]
        for c in class_names:
            row_norm[r][c] = (counts["matrix"][r][c] / total) if total else 0.0

    model_out = os.path.join(out_dir, model_tag)
    os.makedirs(model_out, exist_ok=True)

    write_matrix_csv(
        os.path.join(model_out, "transition_matrix_raw.csv"),
        counts["matrix"],
        class_names,
        class_names,
    )
    write_float_matrix_csv(
        os.path.join(model_out, "transition_matrix_row_norm.csv"),
        row_norm,
        class_names,
        class_names,
    )
    write_matrix_csv(
        os.path.join(model_out, "no_change_counts.csv"),
        {r: {"no_change": counts["no_change"][r]} for r in class_names},
        class_names,
        ["no_change"],
    )
    write_float_matrix_csv(
        os.path.join(model_out, "no_change_rate.csv"),
        {
            r: {
                "no_change_rate": (
                    counts["no_change"][r] / counts["source_total"][r]
                    if counts["source_total"][r]
                    else 0.0
                )
            }
            for r in class_names
        },
        class_names,
        ["no_change_rate"],
    )
    write_text_summary(
        os.path.join(model_out, "summary.txt"),
        model_tag,
        counts,
        counts["matrix"],
        row_norm,
    )

    total_positions = sum(counts["source_total"].values())
    no_change_total = sum(counts["no_change"].values())
    no_change_pct = (100.0 * no_change_total / total_positions) if total_positions else 0.0

    summary_rows.append({
        "model": model_tag,
        "model_dir": model_dir,
        "sequences_processed": counts["sequences_processed"],
        "total_positions": total_positions,
        "no_change_total": no_change_total,
        "no_change_percent": f"{no_change_pct:.2f}",
        "skipped_non_target": counts["skipped_non_target"],
        "length_mismatch_positions": counts["length_mismatch_positions"],
    })

    print(f"[OK] {model_tag}: wrote CSV outputs to {model_out}")

summary_csv = os.path.join(out_dir, "model_summary.csv")
with open(summary_csv, "w", newline="") as f:
    fieldnames = list(summary_rows[0].keys()) if summary_rows else []
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(summary_rows)

print(f"[DONE] CSV outputs written under {out_dir}")
PY

echo "DONE_PROFILE $OUT_DIR"
