# phylo_core.py
## THIS SCRIPT IS TO TEST THE MAIN LOGIC FOR DISTANCE TREE MATRIX BEFORE CONNECTING TO API

from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO
from Bio import Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.Phylo.BaseTree import Tree

def load_fasta(fp: str, fmt: str = "fasta") -> List[Tuple[str, str]]:
    p = Path(fp)
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(f"FASTA file missing or empty: {fp}")
    # If your input may contain leading comments, use "fasta-blast" here.
    records = list(SeqIO.parse(fp, fmt))

    if len(records) < 2:
        raise ValueError("Need at least two sequences to build a tree.")
    return [(r.id or r.name or f"seq_{i}", str(r.seq)) for i, r in enumerate(records)]

def compute_distance_matrix(ids: List[str], seqs: List[str]) -> DistanceMatrix:
    """
    Build a Biopython DistanceMatrix in LOWER-TRIANGULAR form:
      For n taxa, pass a list of n rows:
        row 0: []
        row 1: [d(1,0)]
        row 2: [d(2,0), d(2,1)]
        ...
        row i: [d(i,0), d(i,1), ..., d(i,i-1)]
    """
    print("[DEBUG] compute_distance_matrix() active file is being used")

    n = len(seqs)
    if n != len(ids):
        raise ValueError("ids and seqs must have the same length")
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicate IDs not allowed: {ids}")

    aligner = PairwiseAligner()
    aligner.mode = "global"

    dm_rows: List[List[float]] = []
    for i in range(n):
        row: List[float] = []
        for j in range(i):  # only distances to previous taxa
            score = aligner.score(seqs[i], seqs[j])
            max_len = max(len(seqs[i]), len(seqs[j]))
            sim = (score / max_len) if max_len else 0.0
            dist = max(0.0, 1.0 - sim)
            row.append(float(dist))
        row.append(0.0)  # diagonal element d(i,i)
        dm_rows.append(row)

    # (Optional) sanity print
    for k, r in enumerate(dm_rows):
        print(f"row {k} length={len(r)} -> {r}")

    return DistanceMatrix(ids, dm_rows)

def build_tree(dm: DistanceMatrix, method: str = "nj") -> Tree:
    constructor = DistanceTreeConstructor()
    if method.lower() == "upgma":
        return constructor.upgma(dm)
    return constructor.nj(dm)  # default

def tree_to_newick(tree: Tree) -> str:
    return tree.format("newick")

def draw(tree: Tree) -> None:
    Phylo.draw(tree)