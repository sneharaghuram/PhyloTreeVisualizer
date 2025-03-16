This is a [Next.js](https://nextjs.org) project bootstrapped with [`create-next-app`](https://nextjs.org/docs/app/api-reference/cli/create-next-app).

## Getting Started

First, run the development server:

```bash
npm run dev
# or
yarn dev
# or
pnpm dev
# or
bun dev
```

Open [http://localhost:3000](http://localhost:3000) with your browser to see the result.

You can start editing the page by modifying `app/page.tsx`. The page auto-updates as you edit the file.

This project uses [`next/font`](https://nextjs.org/docs/app/building-your-application/optimizing/fonts) to automatically optimize and load [Geist](https://vercel.com/font), a new font family for Vercel.

## Learn More

To learn more about Next.js, take a look at the following resources:

- [Next.js Documentation](https://nextjs.org/docs) - learn about Next.js features and API.
- [Learn Next.js](https://nextjs.org/learn) - an interactive Next.js tutorial.

You can check out [the Next.js GitHub repository](https://github.com/vercel/next.js) - your feedback and contributions are welcome!

## Deploy on Vercel

The easiest way to deploy your Next.js app is to use the [Vercel Platform](https://vercel.com/new?utm_medium=default-template&filter=next.js&utm_source=create-next-app&utm_campaign=create-next-app-readme) from the creators of Next.js.

Check out our [Next.js deployment documentation](https://nextjs.org/docs/app/building-your-application/deploying) for more details.

### Phylogenetic Analysis

## Tech Stack

[BioPython Module - Phylo](https://biopython.org/wiki/Phylo)
[FastAPI Endpoint](https://fastapi.tiangolo.com/)
Visualization using D3.js or Plotly

# Alternate Method - ClustalW

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
