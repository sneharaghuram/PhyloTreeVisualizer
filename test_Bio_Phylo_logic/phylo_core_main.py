# main.py
## ## THIS SCRIPT IS TO TEST THE MAIN LOGIC FOR DISTANCE TREE MATRIX BEFORE CONNECTING TO API

import argparse
from phylo_core import load_fasta, compute_distance_matrix, build_tree, tree_to_newick, draw

def main():
    ap = argparse.ArgumentParser(description="Build a phylogenetic tree from a FASTA file")
    ap.add_argument("fasta", help="Path to FASTA file")
    ap.add_argument("--format", default="fasta-pearson", choices=["fasta", "fasta-blast", "fasta-pearson"],
                    help="FASTA parser (use fasta-blast if file has leading comments)")
    ap.add_argument("--method", default="nj", choices=["nj", "upgma"], help="Tree method")
    args = ap.parse_args()

    ids_seqs = load_fasta(args.fasta, fmt=args.format)
    ids = [x[0] for x in ids_seqs]
    seqs = [x[1] for x in ids_seqs]

    print(f"[DEBUG] Loaded {len(ids)} sequences: {ids}")
    dm = compute_distance_matrix(ids, seqs)
    print(f"[DEBUG] Distance matrix taxa: {dm.names}")

    tree = build_tree(dm, method=args.method)
    newick = tree_to_newick(tree)
    print("\n=== NEWICK TREE ===")
    print(newick)
    # draw(tree)

if __name__ == "__main__":
    main()
