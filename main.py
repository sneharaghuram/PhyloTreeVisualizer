from fastapi import FastAPI, File, UploadFile
import aiofiles
import os
from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import PairwiseAligner

app = FastAPI()
UPLOAD_DIR = "uploads"

# Ensure upload directory exists
os.makedirs(UPLOAD_DIR, exist_ok=True)


@app.post("/analyze/")
async def analyze_fasta(file: UploadFile = File(...)):
    """Processes a FASTA file, aligns sequences, and generates a phylogenetic tree."""
    
    file_path = os.path.join(UPLOAD_DIR, file.filename)

    # Step 1: Save uploaded FASTA file
    async with aiofiles.open(file_path, "wb") as f:
        content = await file.read()
        await f.write(content)

    # Step 2: Read sequences from FASTA file
    sequences = list(SeqIO.parse(file_path, "fasta"))
    
    if len(sequences) < 2:
        return {"error": "At least two sequences are required to build a phylogenetic tree."}

    # Step 3: Align sequences using Biopython's PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"

    # Step 4: Compute distance matrix using identity score
    num_sequences = len(sequences)
    matrix = [[0] * num_sequences for _ in range(num_sequences)]

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            score = aligner.score(sequences[i].seq, sequences[j].seq)
            max_length = max(len(sequences[i].seq), len(sequences[j].seq))
            distance = 1 - (score / max_length)  # Convert similarity to evolutionary distance
            matrix[i][j] = matrix[j][i] = distance

    # Convert matrix to Biopython's DistanceCalculator format
    taxa = [seq.id for seq in sequences]
    distance_matrix = DistanceCalculator("identity")._format_matrix(matrix, taxa)

    # Step 5: Build phylogenetic tree using Neighbor-Joining
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Step 6: Convert tree to Newick format
    newick_str = tree.format("newick")

    return {"tree": newick_str}
