from fastapi import FastAPI, File, UploadFile
import aiofiles
import os
import subprocess
from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

app = FastAPI()

UPLOAD_DIR = "uploads"
CLUSTALW_EXEC = "clustalw2"  # Adjust this if needed

# Ensure the upload directory exists
os.makedirs(UPLOAD_DIR, exist_ok=True)


@app.post("/analyze/")
async def analyze_fasta(file: UploadFile = File(...)):
    """Processes a FASTA file, aligns sequences, and generates a phylogenetic tree."""
    
    file_path = os.path.join(UPLOAD_DIR, file.filename)
    
    # Save uploaded file
    async with aiofiles.open(file_path, "wb") as f:
        content = await file.read()
        await f.write(content)
    
    # Step 1: Align sequences using ClustalW
    clustal_output = file_path.replace(".fasta", ".aln")
    dnd_output = file_path.replace(".fasta", ".dnd")

    try:
        subprocess.run([CLUSTALW_EXEC, "-INFILE=" + file_path], check=True)
    except Exception as e:
        return {"error": f"Alignment failed: {str(e)}"}

    # Step 2: Parse aligned sequences
    alignment = AlignIO.read(clustal_output, "clustal")

    # Step 3: Compute evolutionary distances
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    # Step 4: Construct the phylogenetic tree using Neighbor-Joining
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Step 5: Convert tree to Newick format
    newick_str = tree.format("newick")

    return {"tree": newick_str}
