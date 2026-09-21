#!/usr/bin/env bash
#SBATCH --job-name=mpnn_only
#SBATCH --partition=htc
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/%u/slurm_logs/%x_%A_%a.out
#SBATCH --error=/scratch/%u/slurm_logs/%x_%A_%a.err

set -euo pipefail
shopt -s nullglob

SCRATCH_BASE="/scratch/$USER"
mkdir -p "$SCRATCH_BASE/slurm_logs"

INDEX="${SLURM_ARRAY_TASK_ID:-}"
PDB_PATH=""
NUM_SEQS="${NUM_SEQS:-200}"
TEMPERATURE="${TEMPERATURE:-0.1}"
MPNN_MODEL="${MPNN_MODEL:-protein_mpnn}"
FINAL_BASE="${FINAL_BASE:-$HOME/protein_start}"
RUN_TAG="${RUN_TAG:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) INDEX="${2:-}"; shift 2 ;;
    -p|--pdb) PDB_PATH="${2:-}"; shift 2 ;;
    -n|--num_seqs) NUM_SEQS="${2:-}"; shift 2 ;;
    -t|--temp) TEMPERATURE="${2:-}"; shift 2 ;;
    -m|--model)
      case "${2:-}" in
        ligand)  MPNN_MODEL="ligand_mpnn" ;;
        protein) MPNN_MODEL="protein_mpnn" ;;
        soluble) MPNN_MODEL="soluble_mpnn" ;;
        *)       MPNN_MODEL="${2:-}" ;;
      esac
      shift 2
      ;;
    --final-base) FINAL_BASE="${2:-}"; shift 2 ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

case "$MPNN_MODEL" in
  protein_mpnn) MODEL_TAG="protein" ;;
  ligand_mpnn)  MODEL_TAG="ligand" ;;
  soluble_mpnn) MODEL_TAG="soluble" ;;
  *)            MODEL_TAG="$MPNN_MODEL" ;;
esac

[[ -n "${INDEX}" ]] || { echo "ERROR: missing array index" >&2; exit 2; }
[[ -n "${LIST_FILE:-}" ]] || { echo "ERROR: LIST_FILE not set" >&2; exit 2; }
[[ -f "${LIST_FILE}" ]] || { echo "ERROR: LIST_FILE not found: ${LIST_FILE}" >&2; exit 2; }

PDB_PATH="$(
  grep -vE '^[[:space:]]*($|#)' "${LIST_FILE}" | sed -n "${INDEX}p"
)"
PDB_PATH="$(echo "${PDB_PATH}" | tr -d '\r' | xargs || true)"
[[ -n "$PDB_PATH" ]] || { echo "ERROR: no usable line for INDEX=${INDEX}" >&2; exit 2; }

# --- NEW ASSEMBLY BLOCK HERE ---
if [[ "$PDB_PATH" =~ ^[0-9A-Za-z]{4}$ ]]; then
  PDB_ID="$(echo "$PDB_PATH" | tr '[:lower:]' '[:upper:]')"
  ID_LOWER="$(echo "$PDB_ID" | tr '[:upper:]' '[:lower:]')"
  mkdir -p pdb_cache

  RAW_ASM_GZ="pdb_cache/${ID_LOWER}_assembly1.cif.gz"
  RAW_CIF="pdb_cache/${ID_LOWER}_assembly1.cif"
  RAW_PDB="pdb_cache/${ID_LOWER}_assembly1.pdb"

  curl -fsSL --max-time 60 -o "${RAW_ASM_GZ}" \
    "https://files.wwpdb.org/download/${PDB_ID}-assembly1.cif.gz" \
    || { echo "ERROR: download failed"; exit 3; }

  gunzip -c "${RAW_ASM_GZ}" > "${RAW_CIF}"

  python3 - <<PY
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser(QUIET=True)
s = parser.get_structure("${PDB_ID}", "${RAW_CIF}")
io = PDBIO()
io.set_structure(s)
io.save("${RAW_PDB}")
PY

  [[ -s "${RAW_PDB}" ]] || { echo "ERROR: conversion failed"; exit 3; }


# --- CHAIN CHECK ---
CHAIN_COUNT=$(awk '/^ATOM/ {print substr($0,22,1)}' "${RAW_PDB}" | sort -u | wc -l)

if [[ "${CHAIN_COUNT}" -gt 1 ]]; then
  echo "Multimer detected ($CHAIN_COUNT chains) — selecting largest"

  BEST_CHAIN=$(
    awk '
      BEGIN { max=0 }
      /^ATOM/ {
        c=substr($0,22,1)
        count[c]++
      }
      END {
        for (c in count) {
          if (count[c] > max) {
            max = count[c]
            best = c
          }
        }
        print best
      }
    ' "${RAW_PDB}"
  )

  awk -v c="$BEST_CHAIN" '
    /^ATOM/ {
      if (substr($0,22,1) != c) next
      line = $0
      line = substr(line,1,21) "A" substr(line,23)
      print line
    }
    END { print "TER" }
  ' "${RAW_PDB}" > "${RAW_PDB}.single"

  mv -f "${RAW_PDB}.single" "${RAW_PDB}"
fi


PDB_PATH="${RAW_PDB}"

fi  
module load mamba/latest

MPNN_ENV="${MPNN_ENV:-ligandmpnn}"
LIGANDMPNN_PATH="${LIGANDMPNN_PATH:-$HOME/LigandMPNN}"
LIGANDMPNN_PY="${LIGANDMPNN_PY:-python}"

run_in_env() {
  local env="$1"; shift
  set +u
  source activate "$env"
  set -u
  "$@"
}

PDB_ABS="$(realpath "$PDB_PATH")"
PDB_BASE="$(basename "$PDB_ABS")"
PDB_STEM="${PDB_BASE%.*}"

FINAL_DIR="${FINAL_BASE}/${MODEL_TAG}/${PDB_STEM}"
mkdir -p "$FINAL_DIR/seqs"

cp -f "$PDB_ABS" "$FINAL_DIR/input.pdb"

awk '
function is_std(res) {
  return (res ~ /^(ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL)$/)
}
BEGIN { prev_resid=""; newi=0 }
(/^ATOM/) {
  res=substr($0,18,3); gsub(/^ +| +$/,"",res)
  if (!is_std(res)) next
  chain="A"
  resid=substr($0,23,4)
  reskey=resid res
  if (reskey != prev_resid) { newi++; prev_resid=reskey }
  line=$0
  line = substr(line,1,21) chain substr(line,23)
  printf("%s%4d%s\n", substr(line,1,22), newi, substr(line,27))
}
END { print "TER" }
' "$FINAL_DIR/input.pdb" > "$FINAL_DIR/ref_clean.pdb"

run_in_env "$MPNN_ENV" python - "$FINAL_DIR/ref_clean.pdb" "$FINAL_DIR/seqs/seq_0000_wt.fasta" <<'PY'
import sys
pdb_path, fasta_path = sys.argv[1], sys.argv[2]
aa3_to_1 = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}
seq=[]
seen=set()
with open(pdb_path) as f:
  for line in f:
    if not line.startswith("ATOM"): continue
    resn=line[17:20].strip()
    if resn not in aa3_to_1: continue
    resid=line[22:26].strip()
    key=(resid,resn)
    if key in seen: continue
    seen.add(key)
    seq.append(aa3_to_1[resn])
s="".join(seq)
with open(fasta_path,"w") as o:
  o.write(">WT\n")
  for i in range(0,len(s),80):
    o.write(s[i:i+80]+"\n")
PY

WORKDIR="$SCRATCH_BASE/mpnn_${MODEL_TAG}_${PDB_STEM}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$WORKDIR/mpnn_out"

START_TS=$(date +%s)

run_in_env "$MPNN_ENV" bash -lc "
  set -e
  cd '$LIGANDMPNN_PATH' || exit 1
  $LIGANDMPNN_PY run.py \
    --pdb_path '$(realpath "$FINAL_DIR/ref_clean.pdb")' \
    --out_folder '$(realpath "$WORKDIR/mpnn_out")' \
    --batch_size '$NUM_SEQS' \
    --temperature '$TEMPERATURE' \
    --model_type '$MPNN_MODEL'
"

END_TS=$(date +%s)
ELAPSED_SECONDS=$((END_TS - START_TS))

FOUND_FA="$(find "$WORKDIR/mpnn_out" -maxdepth 2 -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fas' \) -print -quit)"
[[ -n "${FOUND_FA}" && -s "${FOUND_FA}" ]] || { echo "ERROR: no MPNN fasta found" >&2; exit 4; }

run_in_env "$MPNN_ENV" python - "$FOUND_FA" "$FINAL_DIR/seqs" "$FINAL_DIR/mpnn_metrics.csv" "$FINAL_DIR/seqs/seq_0000_wt.fasta" "$ELAPSED_SECONDS" <<'PY'
import sys, os, csv, re

fasta, outdir, metrics_csv, wt_fasta, elapsed_seconds = sys.argv[1:6]
elapsed_seconds = int(elapsed_seconds)
os.makedirs(outdir, exist_ok=True)

KV_RE = re.compile(r"(\w+)=([^,>]+)")

def parse_fasta(path):
  h = None
  s = []
  with open(path) as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      if line.startswith('>'):
        if h is not None:
          yield h, "".join(s)
        h = line[1:].strip()
        s = []
      else:
        s.append(line)
  if h is not None:
    yield h, "".join(s)

def parse_header(header):
  items = {}
  for k, v in KV_RE.findall(header):
    items[k] = v
  items['name'] = header.split(',', 1)[0].strip()
  return items

def read_single_seq(path):
  seq = []
  with open(path) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('>'):
        continue
      seq.append(line)
  return "".join(seq)

def seq_identity(a, b):
  n = min(len(a), len(b))
  if n == 0:
    return ""
  return f"{sum(x == y for x, y in zip(a[:n], b[:n])) / n:.4f}"

wt_seq = read_single_seq(wt_fasta)
rows = []
i = 1
for h, seq in parse_fasta(fasta):
  if not seq:
    continue
  fn = os.path.join(outdir, f"seq_{i:04d}.fasta")
  with open(fn, 'w') as o:
    o.write(f">seq_{i:04d}\n{seq}\n")

  meta = parse_header(h)
  rows.append({
    'sequence_file': f"seq_{i:04d}.fasta",
    'header_name': meta.get('name', ''),
    'score': meta.get('score', ''),
    'global_score': meta.get('global_score', ''),
    'overall_confidence': meta.get('overall_confidence', ''),
    'ligand_confidence': meta.get('ligand_confidence', ''),
    'seq_rec': meta.get('seq_rec', ''),
    'seq_to_wt_similarity': seq_identity(seq, wt_seq),
    'seq_len': len(seq),
    'seconds': elapsed_seconds,
    'sequence': seq,
  })
  i += 1

with open(metrics_csv, 'w', newline='') as f:
  fieldnames = [
    'sequence_file', 'header_name', 'score', 'global_score',
    'overall_confidence', 'ligand_confidence', 'seq_rec',
    'seq_to_wt_similarity', 'seq_len', 'seconds', 'sequence'
  ]
  writer = csv.DictWriter(f, fieldnames=fieldnames)
  writer.writeheader()
  writer.writerows(rows)

print(f'Wrote {i-1} MPNN seq files')
print(f'Wrote metrics CSV: {metrics_csv}')
PY

echo "$PDB_STEM" > "$FINAL_DIR/pdb_stem.txt"
echo "DONE_MPNN_ONLY $FINAL_DIR"
