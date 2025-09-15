# phylo-gen-visualizer 🧬🌳

A lightweight web application for **phylogenetic tree generation and visualization**.  
This project ingests raw sequencing data, processes it through an automated backend pipeline, and produces **interactive phylogenetic trees** for exploration and analysis.

---

## 🌍 Why This Project

Evolutionary biology seeks to answer fundamental questions about how species arise, adapt, and diverge over time.  
One of its most powerful tools is the **phylogenetic tree**, which reveals evolutionary relationships through genetic data.

However, most existing tools for phylogenetic analysis are either:

- **Command-line only** (requiring deep bioinformatics expertise), or
- **Not easily reproducible/deployable** across platforms.

This project was born out of my interest in **evolutionary biology** and my passion for **computer science**.  
By combining the two, I wanted to create a tool that:

- Automates the tedious parts of **data parsing, alignment, and visualization**,
- Makes phylogenetic analysis **accessible and interactive** via a modern web app,
- Demonstrates how **software engineering practices** (Docker, serverless deployment, CI/CD) can strengthen reproducibility in biology research.

In short: it’s a step toward making evolutionary analysis more usable, scalable, and available to both researchers and students.

## 🚀 Features

- **Sequencing Data Ingestion**: Upload FASTA/sequence data for analysis.
- **Automated Workflow**:
  - Parsing and preprocessing
  - Phylogenetic tree construction
- **Interactive Visualization**: Explore tree structures directly in the browser.
- **Serverless Hosting**: Deployed on **Vercel** and **Railway** for easy access without managing infrastructure.

---

## 🛠️ Tech Stack

- **Backend**: Python (BioPython, NumPy, Matplotlib)
- **Workflow Automation**: Custom scripts for parsing, alignment, and visualization
- **Containerization**: Docker
- **Frontend / Deployment**: Next.js (React) + Vercel
- **Visualization**: D3.js & Matplotlib

# Usage

## Backend

Start your server using uvicorn phylo_tree:app --reload
Terminate Ctrl + C

To test locally, run this command in the terminal:
curl -i -X POST "http://127.0.0.1:8000/analyze/" \
 -F "file=@example.fasta"

To test the deployed server on railway in the command line:
curl -i -X POST "https://phylotreevisualizer-production.up.railway.app/analyze/" \
 -F "file=@example.fasta"

## Frontend

npm run dev

The file needs to contain at least 2 genomic sequences and they have to be of equal length. Example,

```
>Human
ATGCGTACGTTAGCGTACGTAGCTAGCTAGTACG
>Chimpanzee
ATGCGTACGTTAGCGTACGTAGCTAGCTGGTACG
>Mouse
ATGCGTACGTTAGCGTACGTGGCTAGCTGGTTCG
```

## Methods for Improvement

- Generate Cleaner trees for smaller datasets
- Incorporate MultipleSeqAlignment so that you don't have to manually create the matrix
- Use a model-corrected distance (e.g., JC69/K80) before tree building (scikit-bio provides these easily), or run a real MSA (MAFFT/MUSCLE) if sequences aren’t aligned.
- Better UI
  - d3.zoom
  - export SVG, PNG

## Other Bugs and Fixes

- Distance matrix should have zeroes along the diagonal
- FASTA Blast or FASTA Pearson format does not allow files that start with #, you need to update to include the FATSA-BLAST or FASTA-PEARSON format
- FASTA Blast only for parsing FASTA inside BLAST output files, so it will collapse header names to the same name
- React import vs export default
- Adding CORS Middleware to prevent browser from blocking requests to different origins, in our case we add the FastAPI url to make sure that request is allowed
- The tree uses Bezier curves and rather than rectangular cladograms
- Refactored custom UI components (`Input`, `Button`, `Card`) to use proper TypeScript types with optional `className` and uncontrolled file inputs, ensuring clean Next.js builds and better type safety.
- Add Vercel deployed URL in CORS and railway generated URL in Vercel environment variables
- Watch out for endpoint analyze versus analyze/
- In Network, Response Headers should see Access-Control-Allow-Origin: your-url, add a regex URL for Costco

## Learnings

### CORS

CORS stands for **Cross-Origin Resource Sharing**. It’s a security mechanism built into browsers that decides **which websites are allowed to make requests to which servers**.

Imagine you’re logged into your bank’s website. Without CORS, a malicious site you visit in another tab could secretly make a request to your bank’s API (with your cookies) and transfer money. That’s **cross-site request forgery (CSRF)**.

CORS prevents this by default: browsers **block requests to a different origin** (different domain, port, or protocol) unless the server explicitly says “yes, I allow this.”

When your frontend (running at `http://localhost:3000`) tries to call your FastAPI backend (`http://127.0.0.1:8000`), the browser notices:

- These are different **origins** (`localhost:3000` vs `127.0.0.1:8000`).
- Before sending the real request, the browser first sends an **OPTIONS preflight request** asking:

  > “Hey server, do you allow requests from [http://localhost:3000](http://localhost:3000) with method POST and header Content-Type=…?”

The FastAPI server must respond with headers like:

```
Access-Control-Allow-Origin: http://localhost:3000
Access-Control-Allow-Methods: POST
Access-Control-Allow-Headers: Content-Type
```

Only then will the browser proceed with the actual POST.

If those headers aren’t present, you’ll see a “Network Error” in your frontend (even if curl works fine — because curl doesn’t enforce CORS).

- curl works (no CORS in command-line clients).
- Browser blocks the request (CORS protection).
- Fix: add the `CORSMiddleware` in FastAPI so it replies with the right headers for `http://localhost:3000`.

### Alternate Method - ClustalW (Did not use due dependency issues)

March 2025
Note: This project does not use clustalW because of dependency issues but it is a viable option

Clustal-W Usage and Installation
Clustal-W is a tool for multiple sequence alignment, commonly used in bioinformatics for aligning DNA or protein sequences. To use it, simply install the tool and run it on your sequence files.

Installation Issues:
Dependency Conflicts:

libcxx version conflicts may occur when installing Clustal-W via Conda, especially with Python 3.8. Clustal-W requires libcxx >=18, while Python 3.8 often requires an older version.
Solution: Create a new Conda environment with Python 3.7:
bash
Copy
conda create -n clustal-env python=3.7 -c conda-forge
conda activate clustal-env
conda install -c bioconda clustalw
Homebrew on macOS:

Clustal-W may not be available directly via Homebrew for macOS (Apple Silicon).
Solution: Download precompiled binaries from Clustal-W Download and manually install them by placing the binary in a directory included in your PATH.
Usage:
After installation, run Clustal-W from the terminal:
bash
Copy
clustalw
