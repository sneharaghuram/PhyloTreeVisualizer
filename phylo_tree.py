from fastapi import FastAPI, File, UploadFile
import aiofiles
import os
from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator

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

    for i in sequences:
        print(i)

    # Step 3: Align sequences using Biopython's PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    # Create a list to store the alignments in the correct format
    alignments = []

    # Align each pair of sequences and store the alignments
    for seq1 in sequences:
        for seq2 in sequences:
            if seq1.id != seq2.id:  # Avoid aligning the same sequences
                alignment = aligner.align(seq1.seq, seq2.seq)
                alignments.append(alignment[0])  # Store the first alignment (best match)

    # Write alignments to a temporary file in FASTA format
    with open("temp.aln", "w") as aln:
        for alignment in alignments:
            # Write each alignment in FASTA format (adjust the format as needed)
            aln.write(f"Str\n")
            aln.write(str(alignment[0]) + "\n")

    
    # Step 4: Read the alignment file
    with open("temp.aln", "r") as aln:
        alignment = AlignIO.read(aln, "fasta")  # Or use 'clustal' if that's the format

    print(type(alignment))  # Should print: <class 'Bio.Align.MultipleSeqAlignment'>

    # Step 5: Compute the distance matrix
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    print(distance_matrix)

    return distance_matrix

    # Step X: Delete temp.aln

    # # Step 4: Compute distance matrix using identity score
    # num_sequences = len(sequences)
    # matrix = [[0] * num_sequences for _ in range(num_sequences)]

    # for i in range(num_sequences):
    #     for j in range(i + 1, num_sequences):
    #         score = aligner.score(sequences[i].seq, sequences[j].seq)
    #         max_length = max(len(sequences[i].seq), len(sequences[j].seq))
    #         distance = 1 - (score / max_length)  # Convert similarity to evolutionary distance
    #         matrix[i][j] = matrix[j][i] = distance

    # # Step 5: Build phylogenetic tree using Neighbor-Joining
    # # Create a list of taxa (sequence identifiers)
    # taxa = [seq.id for seq in sequences]
    
    # # Convert the distance matrix to a DistanceMatrix object
    # dist_matrix = DistanceMatrix(taxa, matrix)

    # # Now use the DistanceTreeConstructor to build the tree
    # constructor = DistanceTreeConstructor()
    # tree = constructor.nj(dist_matrix)  # Pass the DistanceMatrix object

    # # Step 6: Convert tree to Newick format
    # newick_str = tree.format("newick")

    # return {"tree": newick_str}
