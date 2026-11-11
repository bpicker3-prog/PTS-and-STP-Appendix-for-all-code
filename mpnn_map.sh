#!/usr/bin/env bash
set -euo pipefail

# Terminal colors (disabled for redirected output or when NO_COLOR is set)
if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
  BOLD='\033[1m'; GREEN='\033[0;32m'; CYAN='\033[0;36m'; YELLOW='\033[0;33m'; RED='\033[0;31m'; RESET='\033[0m'
else
  BOLD=''; GREEN=''; CYAN=''; YELLOW=''; RED=''; RESET=''
fi
info() { printf "%b[SETUP]%b %s\n" "$CYAN" "$RESET" "$*"; }
warn() { printf "%b[WARNING]%b %s\n" "$YELLOW" "$RESET" "$*" >&2; }
fail() { printf "%b[ERROR]%b %s\n" "$RED" "$RESET" "$*" >&2; }

RUNS=3
BASE_SEED=42
OUT_BASE="mpnn_benchmark_results"
NUM_SEQS=200
TEMP=0.1
MODELS="protein ligand soluble"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNTIME_DIR="${MPNN_BENCHMARK_RUNTIME:-$SCRIPT_DIR/.mpnn_benchmark_runtime}"
PDB_CACHE_DIR="${MPNN_PDB_CACHE:-$RUNTIME_DIR/pdb_cache}"
export MPNN_PDB_CACHE="$PDB_CACHE_DIR"
LIGANDMPNN_PATH="${LIGANDMPNN_PATH:-$RUNTIME_DIR/LigandMPNN}"
ANALYSIS_PYTHON_BIN="${ANALYSIS_PYTHON_BIN:-$RUNTIME_DIR/benchmark-analysis-python}"
MPNN_PYTHON_BIN="${MPNN_PYTHON_BIN:-$RUNTIME_DIR/benchmark-mpnn-python}"
INPUT_ARGS=()

usage(){ cat <<'EOF'
MPNN all-in-one replicate-aware mutation benchmark

Input may be one protein, many proteins, or a text list containing local PDB/mmCIF paths or 4-character PDB IDs. Biological assemblies are reduced to the largest standard-amino-acid protein chain before MPNN design.

Usage:
  ./mpnn_mutation_benchmark.sh --pdb protein.pdb [options]
  ./mpnn_mutation_benchmark.sh --pdb 1CRN --pdb 1UBQ [options]
  ./mpnn_mutation_benchmark.sh --list proteins.txt [options]

Options:
  -p, --pdb FILE_OR_ID    Repeatable protein input
  -l, --list FILE         One protein path/PDB ID per line
  -r, --runs N            Independent MPNN runs per model (default: 3)
      --seed N            Base seed (default: 42); run r uses seed+r-1
  -n, --num-seqs N        Designs per model per run (default: 200)
  -t, --temp FLOAT        Sampling temperature (default: 0.1)
  -m, --models STRING     protein ligand soluble (default: all three)
  -o, --out DIR           Benchmark output directory
      --runtime DIR      Pixi runtime/cache directory (three isolated workspaces)
      --ligandmpnn PATH   Existing LigandMPNN checkout (optional)
  -h, --help              Show this help

Example:
  ./mpnn_mutation_benchmark.sh --list proteins.txt --runs 5 --seed 12345 --num-seqs 200 --temp 0.1
EOF
}

while [[ $# -gt 0 ]]; do
 case "$1" in
  -p|--pdb) INPUT_ARGS+=(--pdb "${2:?}"); shift 2;;
  -l|--list) INPUT_ARGS+=(--list "${2:?}"); shift 2;;
  -r|--runs) RUNS="${2:?}"; shift 2;;
  --seed) BASE_SEED="${2:?}"; shift 2;;
  -n|--num-seqs) NUM_SEQS="${2:?}"; shift 2;;
  -t|--temp) TEMP="${2:?}"; shift 2;;
  -m|--models) MODELS="${2:?}"; shift 2;;
  -o|--out) OUT_BASE="${2:?}"; shift 2;;
  --runtime) RUNTIME_DIR="${2:?}"; LIGANDMPNN_PATH="$RUNTIME_DIR/LigandMPNN"; ANALYSIS_PYTHON_BIN="$RUNTIME_DIR/benchmark-analysis-python"; MPNN_PYTHON_BIN="$RUNTIME_DIR/benchmark-mpnn-python"; shift 2;;
  --ligandmpnn) LIGANDMPNN_PATH="${2:?}"; shift 2;;
  -h|--help) usage; exit 0;;
  *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 2;;
 esac
done
[[ ${#INPUT_ARGS[@]} -gt 0 ]] || { echo "ERROR: provide --pdb or --list" >&2; exit 2; }
[[ "$RUNS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: --runs must be positive integer" >&2; exit 2; }
[[ "$BASE_SEED" =~ ^[0-9]+$ ]] || { echo "ERROR: --seed must be non-negative integer" >&2; exit 2; }

# -----------------------------------------------------------------------------
# Reproducible Pixi runtime bootstrap
# -----------------------------------------------------------------------------
bootstrap_pixi() {
  if command -v pixi >/dev/null 2>&1; then return 0; fi
  if [[ -x "$HOME/.pixi/bin/pixi" ]]; then export PATH="$HOME/.pixi/bin:$PATH"; return 0; fi
  echo "Pixi is not installed; installing Pixi with the official installer..."
  command -v curl >/dev/null 2>&1 || { echo "ERROR: curl is required to install Pixi." >&2; exit 2; }
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="$HOME/.pixi/bin:$PATH"
  command -v pixi >/dev/null 2>&1 || { echo "ERROR: Pixi installation completed but pixi is not on PATH." >&2; exit 2; }
}

setup_runtime() {
  bootstrap_pixi
  mkdir -p "$RUNTIME_DIR" "$RUNTIME_DIR/inference" "$RUNTIME_DIR/analysis"

  local inference_manifest="$RUNTIME_DIR/inference/pixi.toml"
  local analysis_manifest="$RUNTIME_DIR/analysis/pixi.toml"

  cat > "$inference_manifest" <<'__PIXI_INFERENCE__'
[workspace]
name = "mpnn-benchmark-inference"
channels = ["conda-forge"]
platforms = ["linux-64"]

[dependencies]
python = "3.11.*"
numpy = ">=1.24.4,<2"
scipy = "1.12.*"
biopython = "1.79.*"
prody = "2.4.1.*"
pytorch-cpu = "*"
setuptools = "<81"
pip = "*"
git = "*"
curl = "*"
gzip = "*"
wget = "*"

[pypi-dependencies]
ml-collections = "==0.1.1"
dm-tree = "==0.1.8"
__PIXI_INFERENCE__

  cat > "$analysis_manifest" <<'__PIXI_ANALYSIS__'
[workspace]
name = "mpnn-benchmark-analysis"
channels = ["conda-forge"]
platforms = ["linux-64"]

[dependencies]
python = "3.11.*"
numpy = ">=1.24,<2"
scipy = ">=1.10,<2"
biopython = ">=1.79,<2"
matplotlib = ">=3.7,<4"
__PIXI_ANALYSIS__

  info "Preparing isolated Pixi LigandMPNN inference environment..."
  pixi install --manifest-path "$inference_manifest"
  info "Preparing isolated Pixi analysis environment..."
  pixi install --manifest-path "$analysis_manifest"

  cat > "$MPNN_PYTHON_BIN" <<EOF
#!/usr/bin/env bash
exec pixi run --manifest-path "$inference_manifest" python "\$@"
EOF
  chmod +x "$MPNN_PYTHON_BIN"

  cat > "$ANALYSIS_PYTHON_BIN" <<EOF
#!/usr/bin/env bash
MPLBACKEND=Agg exec pixi run --manifest-path "$analysis_manifest" python "\$@"
EOF
  chmod +x "$ANALYSIS_PYTHON_BIN"

  if [[ ! -f "$LIGANDMPNN_PATH/run.py" ]]; then
    info "Cloning official LigandMPNN..."
    pixi run --manifest-path "$inference_manifest" git clone --depth 1 https://github.com/dauparas/LigandMPNN.git "$LIGANDMPNN_PATH"
  fi
  if [[ ! -f "$LIGANDMPNN_PATH/model_params/proteinmpnn_v_48_020.pt" || ! -f "$LIGANDMPNN_PATH/model_params/ligandmpnn_v_32_010_25.pt" || ! -f "$LIGANDMPNN_PATH/model_params/solublempnn_v_48_020.pt" ]]; then
    info "Downloading LigandMPNN model parameters..."
    (cd "$LIGANDMPNN_PATH" && pixi run --manifest-path "$inference_manifest" bash get_model_params.sh "./model_params")
  fi

  # Compatibility patch for LigandMPNN bundled OpenFold with NumPy >=1.24.
  # Replace only the removed np.int alias; preserve np.int32 / np.int64.
  if [[ -f "$LIGANDMPNN_PATH/openfold/np/residue_constants.py" ]]; then
    sed -i -E 's/np\.int([^0-9A-Za-z_]|$)/int\1/g' "$LIGANDMPNN_PATH/openfold/np/residue_constants.py"
  fi

  info "Validating isolated runtimes..."
  "$MPNN_PYTHON_BIN" - <<'PYVALIDATE'
import numpy, scipy, Bio, prody, torch
print("  inference stack: OK")
print("  torch:", torch.__version__)
print("  ProDy:", prody.__version__)
PYVALIDATE
  "$ANALYSIS_PYTHON_BIN" - <<'PYVALIDATE'
import numpy, scipy, matplotlib, Bio
print("  analysis stack: OK")
print("  matplotlib:", matplotlib.__version__)
PYVALIDATE
}

setup_runtime

ENGINE="$(mktemp "${TMPDIR:-/tmp}/mpnn_benchmark_engine.XXXXXX.sh")"
cleanup_engine() { rm -f "$ENGINE"; }
trap cleanup_engine EXIT
cat > "$ENGINE" <<'__MPNN_ENGINE__'
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# MPNN Mutation Benchmark - local workstation version
# Core benchmark: sequence mutation analysis; no PyMOL/DSS/structural bucket dependency
# Requires a local LigandMPNN checkout containing run.py and its model weights.

NUM_SEQS=200
TEMPERATURE=0.1
RANDOM_SEED="${MPNN_SEED:-42}"
MODELS="protein ligand soluble"
OUT_BASE="mpnn_benchmark_results"
LIGANDMPNN_PATH="${LIGANDMPNN_PATH:-$HOME/LigandMPNN}"
ANALYSIS_PYTHON_BIN="${ANALYSIS_PYTHON_BIN:-python3}"
MPNN_PYTHON_BIN="${MPNN_PYTHON_BIN:-python3}"
INPUTS=()
LIST_FILE=""

usage() {
cat <<'EOF'
Usage:
  ./mpnn_mutation_benchmark.sh --pdb protein.pdb [--pdb protein2.pdb ...] [options]
  ./mpnn_mutation_benchmark.sh --pdb 1ABC [--pdb 2XYZ ...] [options]
  ./mpnn_mutation_benchmark.sh --list proteins.txt [options]

Options:
  -p, --pdb FILE_OR_ID       PDB file/mmCIF file or 4-character PDB ID; repeatable
  -l, --list FILE            Text file with one PDB path/ID per line
  -n, --num-seqs N           Designs per model (default: 200)
  -t, --temp FLOAT           Sampling temperature (default: 0.1)
      --seed N               MPNN random seed (default: 42)
  -m, --models STRING        Space/comma separated: protein ligand soluble
  -o, --out DIR              Output directory (default: mpnn_benchmark_results)
      --ligandmpnn PATH      Local LigandMPNN checkout (default: $LIGANDMPNN_PATH or ~/LigandMPNN)
      --analysis-python PATH Analysis Python executable (internal/advanced)
      --mpnn-python PATH     LigandMPNN Python executable (internal/advanced)
  -h, --help                 Show help

Examples:
  ./mpnn_mutation_benchmark.sh -p 1CRN
  ./mpnn_mutation_benchmark.sh -p myprotein.pdb -n 100 -m "protein ligand soluble"
  ./mpnn_mutation_benchmark.sh --list pdbs.txt --out benchmark_run
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--pdb) INPUTS+=("${2:?missing value for $1}"); shift 2 ;;
    -l|--list) LIST_FILE="${2:?missing value for $1}"; shift 2 ;;
    -n|--num-seqs) NUM_SEQS="${2:?missing value for $1}"; shift 2 ;;
    -t|--temp) TEMPERATURE="${2:?missing value for $1}"; shift 2 ;;
    --seed) RANDOM_SEED="${2:?missing value for $1}"; shift 2 ;;
    -m|--models) MODELS="${2:?missing value for $1}"; shift 2 ;;
    -o|--out) OUT_BASE="${2:?missing value for $1}"; shift 2 ;;
    --ligandmpnn) LIGANDMPNN_PATH="${2:?missing value for $1}"; shift 2 ;;
    --analysis-python) ANALYSIS_PYTHON_BIN="${2:?missing value for $1}"; shift 2 ;;
    --mpnn-python) MPNN_PYTHON_BIN="${2:?missing value for $1}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

MODELS="${MODELS//,/ }"
if [[ -n "$LIST_FILE" ]]; then
  [[ -f "$LIST_FILE" ]] || { echo "ERROR: list not found: $LIST_FILE" >&2; exit 2; }
  while IFS= read -r x; do
    x="${x%%#*}"; x="$(printf '%s' "$x" | xargs)"
    [[ -n "$x" ]] && INPUTS+=("$x")
  done < "$LIST_FILE"
fi
[[ ${#INPUTS[@]} -gt 0 ]] || { echo "ERROR: provide --pdb or --list" >&2; exit 2; }
[[ -f "$LIGANDMPNN_PATH/run.py" ]] || { echo "ERROR: $LIGANDMPNN_PATH/run.py not found. Use --ligandmpnn PATH." >&2; exit 2; }
[[ -x "$ANALYSIS_PYTHON_BIN" ]] || command -v "$ANALYSIS_PYTHON_BIN" >/dev/null || { echo "ERROR: analysis Python not found: $ANALYSIS_PYTHON_BIN" >&2; exit 2; }
[[ -x "$MPNN_PYTHON_BIN" ]] || command -v "$MPNN_PYTHON_BIN" >/dev/null || { echo "ERROR: MPNN Python not found: $MPNN_PYTHON_BIN" >&2; exit 2; }
command -v curl >/dev/null || { echo "ERROR: curl is required" >&2; exit 2; }

"$ANALYSIS_PYTHON_BIN" - <<'PY'
mods=[]
for m in ('numpy','scipy','matplotlib','Bio'):
    try: __import__(m)
    except Exception: mods.append(m)
if mods:
    raise SystemExit('ERROR: missing Python packages: '+', '.join(mods)+'\nInstall with: pip install numpy scipy matplotlib biopython')
PY

mkdir -p "$OUT_BASE/inputs" "$OUT_BASE/runs" "$OUT_BASE/profiles" "$OUT_BASE/summary"
OUT_BASE="$(cd "$OUT_BASE" && pwd)"
FAIL_LOG="$OUT_BASE/summary/failed_inputs.tsv"
printf 'input\tprotein\tstage\tmodel\treason\tdetail\n' > "$FAIL_LOG"

normalize_model() {
  case "$1" in
    protein|protein_mpnn) echo "protein_mpnn protein" ;;
    ligand|ligand_mpnn) echo "ligand_mpnn ligand" ;;
    soluble|soluble_mpnn) echo "soluble_mpnn soluble" ;;
    *) echo "ERROR unknown_model" ;;
  esac
}

prepare_pdb() {
  local input="$1" idx="$2"
  local raw stem
  if [[ "$input" =~ ^[0-9A-Za-z]{4}$ ]]; then
    local id="${input^^}" lower="${input,,}"
    local cache_dir="${MPNN_PDB_CACHE:-$OUT_BASE/pdb_cache}"
    local cache_gz="$cache_dir/${lower}_assembly1.cif.gz"
    local cache_cif="$cache_dir/${lower}_assembly1.cif"
    mkdir -p "$cache_dir"
    if [[ -s "$cache_cif" ]]; then
      echo "[cache] $id biological assembly 1" >&2
    else
      echo "[download] $id biological assembly 1" >&2
      rm -f "$cache_gz.tmp" "$cache_gz"
      if ! timeout --signal=TERM --kill-after=10s 150s \
          curl -fsSL --connect-timeout 20 --max-time 120 \
          --retry 2 --retry-delay 2 --retry-all-errors \
          -o "$cache_gz.tmp" "https://files.wwpdb.org/download/${id}-assembly1.cif.gz"; then
        rm -f "$cache_gz.tmp"
        echo "ERROR: download failed or timed out for $id" >&2
        return 6
      fi
      mv "$cache_gz.tmp" "$cache_gz"
      if ! gzip -t "$cache_gz"; then
        rm -f "$cache_gz"
        echo "ERROR: invalid gzip downloaded for $id" >&2
        return 7
      fi
      gzip -dc "$cache_gz" > "$cache_cif.tmp"
      mv "$cache_cif.tmp" "$cache_cif"
    fi
    raw="$cache_cif"
    stem="${lower}_assembly1"
  else
    [[ -f "$input" ]] || { echo "ERROR: input not found: $input" >&2; return 3; }
    raw="$(cd "$(dirname "$input")" && pwd)/$(basename "$input")"
    stem="$(basename "$input")"; stem="${stem%.*}"
  fi
  local out="$OUT_BASE/inputs/${idx}_${stem}.pdb"
  "$ANALYSIS_PYTHON_BIN" - "$raw" "$out" <<'PY'
import os,sys
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
src,dst=sys.argv[1],sys.argv[2]
parser=MMCIFParser(QUIET=True) if src.lower().endswith(('.cif','.mmcif')) else PDBParser(QUIET=True)
s=parser.get_structure('input',src)
model=next(s.get_models())
std={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
counts={c.id:sum(1 for r in c if r.resname.strip() in std) for c in model}
if not counts: raise SystemExit('ERROR: no standard protein residues found in '+src)
best=max(counts,key=counts.get)
class ProteinSelect(Select):
    def accept_model(self,m): return int(m.id==model.id)
    def accept_chain(self,c): return int(c.id==best)
    def accept_residue(self,r): return int(r.resname.strip() in std)
io=PDBIO(); io.set_structure(s); io.save(dst,ProteinSelect())
print(f'[prepare] {os.path.basename(src)} -> {os.path.basename(dst)}; protein_chains={len(counts)}; largest_chain={best}; residues={counts[best]}')
PY
  printf '%s\n' "$out"
}

validate_backbone() {
  local pdb="$1"
  "$ANALYSIS_PYTHON_BIN" - "$pdb" <<'PYVALID'
import sys
from Bio.PDB import PDBParser
pdb=sys.argv[1]
aa={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
required={'N','CA','C','O'}
s=PDBParser(QUIET=True).get_structure('x',pdb)
model=next(s.get_models())
problems=[]
for chain in model:
    for r in chain:
        if r.resname.strip() not in aa:
            continue
        names={a.get_name().strip() for a in r.get_atoms()}
        missing=sorted(required-names)
        if missing:
            problems.append((chain.id,r.id[1],r.resname.strip(),','.join(missing)))
if problems:
    for chain,resnum,resname,missing in problems:
        print(f'{chain}{resnum}:{resname}:missing={missing}')
    raise SystemExit(10)
print('OK')
PYVALID
}

make_reference() {
  local pdb="$1" outdir="$2"
  mkdir -p "$outdir/seqs"
  cp -f "$pdb" "$outdir/input.pdb"
  "$ANALYSIS_PYTHON_BIN" - "$pdb" "$outdir/ref_clean.pdb" "$outdir/seqs/seq_0000_wt.fasta" <<'PY'
import sys
from Bio.PDB import PDBParser, PDBIO
src,dst,fa=sys.argv[1:]
aa={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
s=PDBParser(QUIET=True).get_structure('x',src); model=next(s.get_models())
chains=list(model.get_chains())
if not chains: raise SystemExit('ERROR: cleaned structure has no protein chain')
chain=chains[0]; residues=[r for r in chain if r.resname.strip() in aa]
if not residues: raise SystemExit('ERROR: cleaned structure has no standard amino acids')
# Reproduce the original preprocessing convention: one chain named A and consecutive residue numbering.
chain.id='A'
for i,r in enumerate(residues,1): r.id=(' ',10000+i,' ')
for i,r in enumerate(residues,1): r.id=(' ',i,' ')
# Remove any non-standard residues left in the selected chain.
for r in list(chain):
 if r.resname.strip() not in aa: chain.detach_child(r.id)
io=PDBIO(); io.set_structure(s); io.save(dst)
seq=''.join(aa[r.resname.strip()] for r in residues)
with open(fa,'w') as f: f.write('>WT\n'+seq+'\n')
PY
}

run_model() {
  local pdb="$1" protein="$2" model_type="$3" tag="$4"
  local final="$OUT_BASE/runs/$tag/$protein" tmp="$OUT_BASE/.tmp_${tag}_${protein}_$$"
  mkdir -p "$final" "$tmp"
  make_reference "$pdb" "$final"
  local start end elapsed
  start=$(date +%s)
  echo "[MPNN] $protein :: $tag ($NUM_SEQS designs, T=$TEMPERATURE)"
( cd "$LIGANDMPNN_PATH" && "$MPNN_PYTHON_BIN" run.py \
    --pdb_path "$final/ref_clean.pdb" \
    --out_folder "$tmp" \
    --batch_size "$NUM_SEQS" \
    --temperature "$TEMPERATURE" \
    --seed "$RANDOM_SEED" \
    --model_type "$model_type" \
    --parse_atoms_with_zero_occupancy 1 )
  end=$(date +%s); elapsed=$((end-start))
  local found
  found="$(find "$tmp" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fas' \) -print -quit)"
  [[ -n "$found" && -s "$found" ]] || { echo "ERROR: no FASTA produced for $protein/$tag" >&2; return 4; }
  "$ANALYSIS_PYTHON_BIN" - "$found" "$final" "$elapsed" "$tag" "$protein" <<'PY'
import csv,os,re,sys
fasta,out,elapsed,model,protein=sys.argv[1:]; seqdir=os.path.join(out,'seqs'); os.makedirs(seqdir,exist_ok=True)
def parse(path):
 h=None;s=[]
 for line in open(path):
  line=line.strip()
  if not line: continue
  if line.startswith('>'):
   if h is not None: yield h,''.join(s)
   h=line[1:];s=[]
  else:s.append(line)
 if h is not None: yield h,''.join(s)
def one(path): return ''.join(x.strip() for x in open(path) if x.strip() and not x.startswith('>'))
wt=one(os.path.join(seqdir,'seq_0000_wt.fasta')); rows=[]; kv=re.compile(r'(\w+)=([^,>]+)')
records=list(parse(fasta))
# Some MPNN FASTAs include the native/reference record first. Drop it if sequence equals WT.
if records and records[0][1].replace('/','')==wt: records=records[1:]
for i,(h,s) in enumerate(records,1):
 s=s.replace('/','')
 with open(os.path.join(seqdir,f'seq_{i:04d}.fasta'),'w') as f:f.write(f'>seq_{i:04d}\n{s}\n')
 meta=dict(kv.findall(h)); n=min(len(wt),len(s)); ident=sum(a==b for a,b in zip(wt[:n],s[:n]))/n if n else float('nan')
 rows.append(dict(protein=protein,model=model,sequence_file=f'seq_{i:04d}.fasta',score=meta.get('score',''),global_score=meta.get('global_score',''),overall_confidence=meta.get('overall_confidence',''),ligand_confidence=meta.get('ligand_confidence',''),seq_rec=meta.get('seq_rec',''),seq_to_wt_similarity=f'{ident:.6f}',seq_len=len(s),seconds=elapsed,sequence=s))
fields=['protein','model','sequence_file','score','global_score','overall_confidence','ligand_confidence','seq_rec','seq_to_wt_similarity','seq_len','seconds','sequence']
with open(os.path.join(out,'mpnn_metrics.csv'),'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
print(f'[metrics] {protein}/{model}: {len(rows)} designs')
PY
  rm -rf "$tmp"
}

PROFILE_PY="$OUT_BASE/.profile.py"
cat > "$PROFILE_PY" <<'PY'
import csv,glob,os,sys
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
root,protein,out=sys.argv[1:]
models=[m for m in ('protein','ligand','soluble') if os.path.isdir(os.path.join(root,m,protein))]
aas=list('ACDEFGHIKLMNPQRSTVWY')
classes=OrderedDict([('hydrophobic',set('AVLIM')),('hydrophilic',set('STNQCG')),('aromatic',set('FWY')),('charged',set('KRHDE')),('special',set('P'))])
def readfa(p):return ''.join(x.strip() for x in open(p) if x.strip() and not x.startswith('>')).upper()
def cls(a):
 for k,v in classes.items():
  if a in v:return k
 return None
def csvmat(path,mat,rows,cols,floaty=False):
 with open(path,'w',newline='') as f:
  w=csv.writer(f);w.writerow(['WT\\Design']+cols)
  for r in rows:w.writerow([r]+[(f'{mat[r][c]:.6f}' if floaty else mat[r][c]) for c in cols])
def heat(path,mat,labels,title):
 a=np.array([[mat[r][c] for c in labels] for r in labels],float);fig,ax=plt.subplots(figsize=(10,8));im=ax.imshow(a,aspect='auto',interpolation='nearest',vmin=0,vmax=1);ax.set_xticks(range(len(labels)),labels=labels,rotation=45,ha='right');ax.set_yticks(range(len(labels)),labels=labels);ax.set_xlabel('Designed residue');ax.set_ylabel('WT residue');ax.set_title(title);fig.colorbar(im,ax=ax,label='Within-WT transition frequency');fig.tight_layout();fig.savefig(path,dpi=200);plt.close(fig)
os.makedirs(out,exist_ok=True); summaries=[]
model_names={'protein':'ProteinMPNN','ligand':'LigandMPNN','soluble':'SolubleMPNN'}
for model in models:
 display_model=model_names.get(model,model)
 d=os.path.join(root,model,protein); wt=readfa(os.path.join(d,'seqs','seq_0000_wt.fasta')); seqs=[readfa(p) for p in sorted(glob.glob(os.path.join(d,'seqs','seq_*.fasta'))) if '0000_wt' not in p]
 aa={r:{c:0 for c in aas} for r in aas}; cm={r:{c:0 for c in classes} for r in classes}; pos=np.zeros((len(wt),len(aas)),int); changed=np.zeros(len(wt),int); valid=np.zeros(len(wt),int)
 for s in seqs:
  for i,(w,x) in enumerate(zip(wt,s)):
   if w not in aa or x not in aa:continue
   aa[w][x]+=1;pos[i,aas.index(x)]+=1;valid[i]+=1;changed[i]+=int(w!=x)
   cw,cx=cls(w),cls(x)
   if cw and cx:cm[cw][cx]+=1
 aan={r:{c:(aa[r][c]/sum(aa[r].values()) if sum(aa[r].values()) else 0) for c in aas} for r in aas}; cn={r:{c:(cm[r][c]/sum(cm[r].values()) if sum(cm[r].values()) else 0) for c in classes} for r in classes}
 mo=os.path.join(out,model);os.makedirs(mo,exist_ok=True);csvmat(os.path.join(mo,'aa_transition_matrix_raw.csv'),aa,aas,aas);csvmat(os.path.join(mo,'aa_transition_matrix_row_norm.csv'),aan,aas,aas,True);csvmat(os.path.join(mo,'class_transition_matrix_raw.csv'),cm,list(classes),list(classes));csvmat(os.path.join(mo,'class_transition_matrix_row_norm.csv'),cn,list(classes),list(classes),True);heat(os.path.join(mo,'aa_transition_heatmap.png'),aan,aas,f'{protein} — {display_model} amino-acid transition profile');heat(os.path.join(mo,'class_transition_heatmap.png'),cn,list(classes),f'{protein} — {display_model} residue-class transition profile')
 with open(os.path.join(mo,'position_mutation_map.csv'),'w',newline='') as f:
  w=csv.writer(f);w.writerow(['position','wt_residue','n_designs','mutation_count','mutation_rate']+['design_'+a for a in aas])
  for i,wtaa in enumerate(wt):w.writerow([i+1,wtaa,int(valid[i]),int(changed[i]),f'{changed[i]/valid[i]:.6f}' if valid[i] else '']+list(pos[i]))
 # position mutation-rate map
 rates=np.divide(changed,valid,out=np.zeros_like(changed,dtype=float),where=valid>0)[None,:];fig,ax=plt.subplots(figsize=(max(10,len(wt)/12),2.8));im=ax.imshow(rates,aspect='auto',vmin=0,vmax=1);ax.set_yticks([0],labels=[display_model]);ax.set_xlabel('Residue position');ax.set_title(f'{protein} — {display_model} per-position mutation rate');fig.colorbar(im,ax=ax,label='Mutation rate');fig.tight_layout();fig.savefig(os.path.join(mo,'position_mutation_map.png'),dpi=200);plt.close(fig)
 summaries.append([protein,model,len(seqs),len(wt),float(changed.sum()/valid.sum()) if valid.sum() else np.nan])
with open(os.path.join(out,'profile_summary.csv'),'w',newline='') as f:w=csv.writer(f);w.writerow(['protein','model','n_designs','protein_length','overall_mutation_rate']);w.writerows(summaries)
PY

SUMMARY_PY="$OUT_BASE/.summarize.py"
cat > "$SUMMARY_PY" <<'PY'
import csv,glob,os,sys
import numpy as np
from scipy.stats import friedmanchisquare,wilcoxon,kruskal,mannwhitneyu
root,out=sys.argv[1:]; models=['protein','ligand','soluble']; labels={'protein':'ProteinMPNN','ligand':'LigandMPNN','soluble':'SolubleMPNN'}
def fnum(x):
 try:return float(x)
 except:return np.nan
rows=[]
for m in models:
 for p in glob.glob(os.path.join(root,m,'*','mpnn_metrics.csv')):
  protein=os.path.basename(os.path.dirname(p)); rr=list(csv.DictReader(open(p)))
  for metric in ('seq_to_wt_similarity','overall_confidence'):
   a=np.array([fnum(r.get(metric,'')) for r in rr]);a=a[np.isfinite(a)]
   if len(a): rows.append(dict(protein=protein,model=labels[m],metric=metric,n=len(a),mean=float(a.mean()),median=float(np.median(a)),sd=float(a.std(ddof=1)) if len(a)>1 else 0.0))
os.makedirs(out,exist_ok=True)
with open(os.path.join(out,'per_protein_summary.csv'),'w',newline='') as f:
 fields=['protein','model','metric','n','mean','median','sd'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
# Across-protein inference: paired by protein, matching the original benchmark logic.
results=[]
for metric in ('seq_to_wt_similarity','overall_confidence'):
 d={}
 for r in rows:
  if r['metric']==metric:d.setdefault(r['protein'],{})[r['model']]=r['mean']
 complete={p:v for p,v in d.items() if all(labels[m] in v for m in models)}
 if len(complete)>=2:
  ps=sorted(complete); arr=[np.array([complete[p][labels[m]] for p in ps]) for m in models]
  if len(complete)>=3:
   st,pv=friedmanchisquare(*arr);results.append([metric,'across proteins','Friedman','all 3 models',len(ps),st,pv,''])
  pairs=[(0,1),(0,2),(1,2)]; raw=[]
  for i,j in pairs:
   try: st,pv=wilcoxon(arr[i],arr[j],zero_method='wilcox',alternative='two-sided')
   except ValueError: st,pv=np.nan,1.0
   raw.append([i,j,st,pv])
  # Holm correction for three planned pairwise tests
  order=sorted(range(len(raw)),key=lambda k:raw[k][3]);adj=[None]*len(raw);run=0
  for rank,k in enumerate(order):run=max(run,min(1,(len(raw)-rank)*raw[k][3]));adj[k]=run
  for k,(i,j,st,pv) in enumerate(raw):results.append([metric,'across proteins','Wilcoxon signed-rank',f'{labels[models[i]]} vs {labels[models[j]]}',len(ps),st,pv,adj[k]])
# Design-level descriptive/exploratory comparisons for each single protein (unpaired samples).
proteins=sorted({r['protein'] for r in rows})
for protein in proteins:
 for metric in ('seq_to_wt_similarity','overall_confidence'):
  groups=[]; names=[]
  for m in models:
   p=os.path.join(root,m,protein,'mpnn_metrics.csv')
   if not os.path.exists(p):continue
   a=np.array([fnum(r.get(metric,'')) for r in csv.DictReader(open(p))]);a=a[np.isfinite(a)]
   if len(a):groups.append(a);names.append(labels[m])
  if len(groups)>=2:
   st,pv=kruskal(*groups);results.append([metric,'within one protein; designs are generated samples, not independent proteins','Kruskal-Wallis',' vs '.join(names),sum(map(len,groups)),st,pv,''])
with open(os.path.join(out,'statistical_tests.csv'),'w',newline='') as f:
 w=csv.writer(f);w.writerow(['metric','scope','test','comparison','n','statistic','p_raw','p_adjusted']);w.writerows(results)
with open(os.path.join(out,'README_stats.txt'),'w') as f:
 f.write('Across-protein Friedman/Wilcoxon tests are the primary inferential benchmark when the same proteins are evaluated by all three models.\n')
 f.write('Single-protein Kruskal-Wallis results compare generated design distributions and should be treated as exploratory because generated sequences are not independent protein replicates.\n')
print(f'[summary] {len({r["protein"] for r in rows})} proteins summarized; {len(results)} statistical-test rows')
PY

printf 'input\tprepared_pdb\tprotein\n' > "$OUT_BASE/summary/input_manifest.tsv"
idx=0
for input in "${INPUTS[@]}"; do
  idx=$((idx+1))
  prep_tmp="$OUT_BASE/summary/.prepare_${idx}.log"
  if ! prepare_pdb "$input" "$idx" >"$prep_tmp" 2>&1; then
    detail="$(tr '\n\t' '  ' < "$prep_tmp" | sed 's/[[:space:]]\+/ /g')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$input" "" "prepare" "" "preparation_failed" "$detail" >> "$FAIL_LOG"
    echo "[WARNING] skipping $input: preparation failed" >&2
    rm -f "$prep_tmp"
    continue
  fi
  prepared="$(tail -n 1 "$prep_tmp")"
  rm -f "$prep_tmp"
  protein="$(basename "$prepared" .pdb)"
  printf '%s\t%s\t%s\n' "$input" "$prepared" "$protein" >> "$OUT_BASE/summary/input_manifest.tsv"

  validation="$OUT_BASE/summary/.validate_${idx}.log"
  if ! validate_backbone "$prepared" >"$validation" 2>&1; then
    detail="$(tr '\n\t' '; ' < "$validation" | sed 's/[[:space:]]\+/ /g')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$input" "$protein" "validation" "" "incomplete_backbone" "$detail" >> "$FAIL_LOG"
    echo "[WARNING] skipping $protein: incomplete N/CA/C/O backbone ($detail)" >&2
    rm -f "$validation"
    continue
  fi
  rm -f "$validation"

  completed=0
  for m in $MODELS; do
    read -r model_type tag <<< "$(normalize_model "$m")"
    [[ "$model_type" != ERROR ]] || { echo "ERROR: model must be protein, ligand, or soluble: $m" >&2; exit 2; }
    model_log="$OUT_BASE/summary/.model_${idx}_${tag}.log"
    set +e
    run_model "$prepared" "$protein" "$model_type" "$tag" > >(tee "$model_log") 2>&1
    rc=$?
    set -e
    if [[ $rc -eq 0 ]]; then
      completed=$((completed+1))
      rm -f "$model_log"
    else
      detail="$(tail -n 20 "$model_log" | tr '\n\t' '  ' | sed 's/[[:space:]]\+/ /g')"
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$input" "$protein" "model" "$tag" "model_failed_exit_${rc}" "$detail" >> "$FAIL_LOG"
      echo "[WARNING] $protein/$tag failed (exit $rc); continuing benchmark" >&2
      rm -f "$model_log"
    fi
  done

  if (( completed > 0 )); then
    set +e
    "$ANALYSIS_PYTHON_BIN" "$PROFILE_PY" "$OUT_BASE/runs" "$protein" "$OUT_BASE/profiles/$protein"
    prc=$?
    set -e
    if [[ $prc -ne 0 ]]; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$input" "$protein" "profile" "" "profile_failed" "profile generation exit $prc" >> "$FAIL_LOG"
      echo "[WARNING] profile generation failed for $protein; continuing" >&2
    fi
  else
    echo "[WARNING] no models completed for $protein; no profile generated" >&2
  fi
done

"$ANALYSIS_PYTHON_BIN" "$SUMMARY_PY" "$OUT_BASE/runs" "$OUT_BASE/summary"
rm -f "$PROFILE_PY" "$SUMMARY_PY"

echo
echo "DONE"
echo "Results: $OUT_BASE"
echo "  Mutation maps: $OUT_BASE/profiles/<protein>/<model>/"
echo "  Per-protein metrics: $OUT_BASE/summary/per_protein_summary.csv"
echo "  Statistical tests: $OUT_BASE/summary/statistical_tests.csv"
echo "  Failed/skipped inputs: $FAIL_LOG"

__MPNN_ENGINE__
chmod +x "$ENGINE"
mkdir -p "$OUT_BASE/replicates" "$OUT_BASE/summary"
OUT_BASE="$(cd "$OUT_BASE" && pwd)"
printf 'run\tseed\ttemperature\tn_sequences\tmodels\tresult_directory\n' > "$OUT_BASE/summary/run_metadata.tsv"

for ((r=1;r<=RUNS;r++)); do
 printf -v run_name 'run_%03d' "$r"
 seed=$((BASE_SEED+r-1))
 run_out="$OUT_BASE/replicates/$run_name"
 echo "=== $run_name/$RUNS | seed=$seed | T=$TEMP ==="
 cmd=("$ENGINE" "${INPUT_ARGS[@]}" --num-seqs "$NUM_SEQS" --temp "$TEMP" --models "$MODELS" --out "$run_out" --ligandmpnn "$LIGANDMPNN_PATH" --analysis-python "$ANALYSIS_PYTHON_BIN" --mpnn-python "$MPNN_PYTHON_BIN" --seed "$seed")
  "${cmd[@]}"
 printf '%s\t%d\t%s\t%s\t%s\t%s\n' "$run_name" "$seed" "$TEMP" "$NUM_SEQS" "$MODELS" "$run_out" >> "$OUT_BASE/summary/run_metadata.tsv"
done

# Core mutation benchmark intentionally performs no secondary-structure bucketing.
# Mutation maps and statistics are computed directly from WT and designed sequences.

"$ANALYSIS_PYTHON_BIN" - "$OUT_BASE" <<'PY'
import csv,glob,os,sys
import numpy as np
from scipy.stats import friedmanchisquare,wilcoxon
root=sys.argv[1]; summary=os.path.join(root,'summary')
labels={'protein':'ProteinMPNN','ligand':'LigandMPNN','soluble':'SolubleMPNN'}
# Gather design metrics from every replicate.
records=[]
for p in glob.glob(os.path.join(root,'replicates','run_*','runs','*','*','mpnn_metrics.csv')):
 run=os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(p)))))
 model=os.path.basename(os.path.dirname(os.path.dirname(p))); protein=os.path.basename(os.path.dirname(p))
 rows=list(csv.DictReader(open(p,newline='')))
 for metric in ('seq_to_wt_similarity','overall_confidence'):
  vals=[]
  for row in rows:
   try: vals.append(float(row.get(metric,'')))
   except: pass
  if vals:
   a=np.asarray(vals,float);records.append(dict(run=run,protein=protein,model=labels.get(model,model),metric=metric,n=len(a),mean=float(a.mean()),median=float(np.median(a)),sd=float(a.std(ddof=1)) if len(a)>1 else 0.0))
fields=['run','protein','model','metric','n','mean','median','sd']
with open(os.path.join(summary,'per_run_summary.csv'),'w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
# Across-run descriptive summary for every protein/model.
agg=[]
for protein in sorted({x['protein'] for x in records}):
 for model in sorted({x['model'] for x in records}):
  for metric in ('seq_to_wt_similarity','overall_confidence'):
   x=np.array([z['mean'] for z in records if z['protein']==protein and z['model']==model and z['metric']==metric])
   if len(x):agg.append(dict(protein=protein,model=model,metric=metric,n_runs=len(x),mean_across_runs=float(x.mean()),sd_across_runs=float(x.std(ddof=1)) if len(x)>1 else 0.0,min_run_mean=float(x.min()),max_run_mean=float(x.max())))
af=['protein','model','metric','n_runs','mean_across_runs','sd_across_runs','min_run_mean','max_run_mean']
with open(os.path.join(summary,'across_run_summary.csv'),'w',newline='') as f:w=csv.DictWriter(f,fieldnames=af);w.writeheader();w.writerows(agg)
# Paired model comparisons using run as the paired computational unit within each protein.
def holm(ps):
 order=sorted(range(len(ps)),key=ps.__getitem__);out=[None]*len(ps);running=0
 for rank,k in enumerate(order):running=max(running,min(1.0,(len(ps)-rank)*ps[k]));out[k]=running
 return out
tests=[]; model_order=['ProteinMPNN','LigandMPNN','SolubleMPNN']
for protein in sorted({x['protein'] for x in records}):
 for metric in ('seq_to_wt_similarity','overall_confidence'):
  d={}
  for z in records:
   if z['protein']==protein and z['metric']==metric:d.setdefault(z['run'],{})[z['model']]=z['mean']
  complete={r:v for r,v in d.items() if all(m in v for m in model_order)}
  if len(complete)>=2:
   rs=sorted(complete); arr=[np.array([complete[r][m] for r in rs]) for m in model_order]
   if len(rs)>=3:
    st,p=friedmanchisquare(*arr);tests.append([protein,metric,'run-level model effect','Friedman','all three models',len(rs),st,p,''])
   raw=[]
   for i,j in ((0,1),(0,2),(1,2)):
    try:st,p=wilcoxon(arr[i],arr[j],zero_method='wilcox',alternative='two-sided')
    except ValueError:st,p=np.nan,1.0
    raw.append((i,j,st,p))
   adj=holm([x[3] for x in raw])
   for k,(i,j,st,p) in enumerate(raw):tests.append([protein,metric,'run-level model effect','Wilcoxon signed-rank',model_order[i]+' vs '+model_order[j],len(rs),st,p,adj[k]])
with open(os.path.join(summary,'run_level_statistical_tests.csv'),'w',newline='') as f:
 w=csv.writer(f);w.writerow(['protein','metric','scope','test','comparison','n_paired_runs','statistic','p_raw','p_holm']);w.writerows(tests)
# Across-protein inference after averaging replicate runs within protein/model.
primary=[]
for metric in ('seq_to_wt_similarity','overall_confidence'):
 d={}
 for z in agg:
  if z['metric']==metric:d.setdefault(z['protein'],{})[z['model']]=z['mean_across_runs']
 complete={p:v for p,v in d.items() if all(m in v for m in model_order)}
 if len(complete)>=2:
  ps=sorted(complete);arr=[np.array([complete[p][m] for p in ps]) for m in model_order]
  if len(ps)>=3:
   st,p= friedmanchisquare(*arr);primary.append([metric,'Friedman','all three models',len(ps),st,p,''])
  raw=[]
  for i,j in ((0,1),(0,2),(1,2)):
   try:st,p=wilcoxon(arr[i],arr[j],zero_method='wilcox',alternative='two-sided')
   except ValueError:st,p=np.nan,1.0
   raw.append((i,j,st,p))
  adj=holm([x[3] for x in raw])
  for k,(i,j,st,p) in enumerate(raw):primary.append([metric,'Wilcoxon signed-rank',model_order[i]+' vs '+model_order[j],len(ps),st,p,adj[k]])
with open(os.path.join(summary,'primary_across_protein_tests.csv'),'w',newline='') as f:
 w=csv.writer(f);w.writerow(['metric','test','comparison','n_matched_proteins','statistic','p_raw','p_holm']);w.writerows(primary)
with open(os.path.join(summary,'README_statistics.txt'),'w') as f:
 f.write('Analysis hierarchy: protein > MPNN model > independent seeded run > generated sequence.\n\n')
 f.write('per_run_summary.csv: mean/median/SD of generated designs within each run.\n')
 f.write('across_run_summary.csv: reproducibility/variation of run means for each protein and model.\n')
 f.write('run_level_statistical_tests.csv: paired model comparisons within a protein, pairing by seeded run.\n')
 f.write('primary_across_protein_tests.csv: primary benchmark inference when multiple proteins are supplied; replicate runs are averaged first so proteins remain the independent units.\n')
 f.write('Generated sequences are not treated as independent proteins.\n')
print('[replicate summary] wrote run-level, across-run, and across-protein statistics')
PY

echo
printf "%b%b✓ BENCHMARK COMPLETE%b\n\n" "$GREEN" "$BOLD" "$RESET"
printf "%bResults%b\n" "$BOLD" "$RESET"
printf "  %-24s %b%s%b\n" "Benchmark root:" "$CYAN" "$OUT_BASE" "$RESET"
printf "  %-24s %b%s%b\n" "Replicate runs:" "$CYAN" "$OUT_BASE/replicates/run_XXX/" "$RESET"
printf "  %-24s %b%s%b\n" "Mutation maps:" "$CYAN" "$OUT_BASE/replicates/run_XXX/profiles/<protein>/<model>/" "$RESET"
printf "  %-24s %b%s%b\n" "Run metadata:" "$CYAN" "$OUT_BASE/summary/run_metadata.tsv" "$RESET"
printf "  %-24s %b%s%b\n" "Run reproducibility:" "$CYAN" "$OUT_BASE/summary/across_run_summary.csv" "$RESET"
printf "  %-24s %b%s%b\n" "Run-level tests:" "$CYAN" "$OUT_BASE/summary/run_level_statistical_tests.csv" "$RESET"
printf "  %-24s %b%s%b\n" "Across-protein tests:" "$CYAN" "$OUT_BASE/summary/primary_across_protein_tests.csv" "$RESET"
printf "  %-24s %b%s%b\n" "Failures by replicate:" "$CYAN" "$OUT_BASE/replicates/run_XXX/summary/failed_inputs.tsv" "$RESET"
printf "  %-24s %b%s%b\n" "Shared PDB cache:" "$CYAN" "$PDB_CACHE_DIR" "$RESET"
printf "\n%bTip:%b set NO_COLOR=1 to disable terminal colors.\n" "$YELLOW" "$RESET"
