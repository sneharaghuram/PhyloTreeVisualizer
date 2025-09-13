from fastapi import FastAPI, File, UploadFile
from fastapi.middleware.cors import CORSMiddleware
import aiofiles, os
from typing import List
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    # allow all Vercel previews + your prod domain
    allow_origin_regex=r"https://.*\.vercel\.app",
    allow_origins=[  # explicit prod & local
        "https://phylo-tree-generator.vercel.app/",
        "http://localhost:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=False, #don't need cookies/auth
    allow_methods=["*"],
    allow_headers=["*"],
)


UPLOAD_DIR = "uploads"
os.makedirs(UPLOAD_DIR, exist_ok=True)

@app.post("/analyze")
@app.post("/analyze/")
async def analyze_fasta(file: UploadFile = File(...)):
    """Processes a FASTA file and returns a phylogenetic tree (Newick)."""

    # --- save upload ---
    file_path = os.path.join(UPLOAD_DIR, file.filename or "upload.fasta")
    async with aiofiles.open(file_path, "wb") as f:
        await f.write(await file.read())

    # --- read sequences ---
    records = list(SeqIO.parse(file_path, "fasta"))  # tolerant of comment lines
    if len(records) < 2:
        return {"error": "At least two sequences are required."}

    ids: List[str] = [r.id.strip() for r in records]
    seqs: List[str] = [str(r.seq) for r in records]

    # --- check for duplicate IDs ---
    if len(ids) != len(set(ids)):
        return {
            "error": "Duplicate sequence names found in FASTA. "
                     "Please ensure each sequence header is unique.",
            "duplicates": [i for i in ids if ids.count(i) > 1],
        }

    # --- require equal lengths for this simple path ---
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        return {
            "error": (
                "Sequences must be equal length for this endpoint (no external MSA used). "
                "Provide pre-aligned sequences or align them first (e.g., MAFFT/MUSCLE)."
            ),
            "lengths_found": sorted(lengths),
            "taxa": ids,
        }

    # --- build alignment in-memory ---
    alignment = MultipleSeqAlignment(
        [SeqRecord(Seq(s), id=i) for i, s in zip(ids, seqs)]
    )

    # --- distances & tree ---
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    return {
        "tree_newick": tree.format("newick"),
        "taxa": ids,
        "num_taxa": len(ids),
        "distance_names": dm.names,
        "method": "nj",
        "distance_model": "identity",
    }
