export default function AboutPage() {
  return (
    <div className="max-w-3xl mx-auto p-8">
      <h2 className="text-2xl font-bold mb-4">About PhyloTree Visualizer</h2>
      <p className="mb-4">
        <strong>PhyloTree Visualizer</strong> is a web application that allows you to
        upload DNA or protein sequences in FASTA format and generate phylogenetic
        trees using Neighbor-Joining (NJ) or UPGMA algorithms. The trees are exported
        in Newick format and visualized interactively with D3.js.
      </p>
      <h3 className="text-xl font-semibold mb-2">Features</h3>
      <ul className="list-disc pl-6 mb-4">
        <li>Upload nucleotide or protein FASTA files</li>
        <p> In this version, the FASTA Files have to be equal in length</p>
        <li>Automatic tree generation with Biopython</li>
        <li>Interactive D3.js visualization (rectangular cladograms)</li>
        <li>Drop Down Menu of Algorithms that create the Phylogenetic Tree: Neighbor Joining or UPGMA</li>
      </ul>
      <h3 className="text-xl font-semibold mb-2">Future Improvements</h3>
      <ul className="list-disc pl-6 mb-4">
        <li>Export Results in Newick Format</li>
        <li>Download tree in pdf/png format</li>
        <li>Use MultSeqAlignment to allow multiple sequences of varying length and perform genome sequencing alignment in the backend</li>
        <li>Improve UI/UX</li>
        <li>Incorporate Character Based Methods such as maximum likelihood or maximum parsimony</li>
      </ul>
      <h3 className="text-xl font-semibold mb-2">Tech Stack</h3>
      <ul className="list-disc pl-6">
        <li>Frontend: Next.js (TypeScript), D3.js, Tailwind</li>
        <li>Backend: FastAPI (Python), Biopython</li>
        <li>Deployment: Vercel (frontend), Railway (backend)</li>
      </ul>
    
      <section className="mt-8">
        <h3 className="text-xl font-semibold mb-2">Some Terms and Definitions</h3>
        <p className="mb-4">
            A <strong>phylogenetic tree</strong> is like a family tree for species or
            genes. It shows how different organisms (or sequences) are related through
            evolution, with branches representing shared ancestry.
        </p>
        <p className="mb-4">
            A <strong>FASTA file</strong> is a simple text format for storing DNA,
            RNA, or protein sequences. Each entry begins with a label (e.g., “Human”)
            on a line starting with <code>&gt;</code>, followed by the sequence itself.
        </p>
        <p className="mb-4">
            <strong>Newick format</strong> is a compact way of writing trees using
            parentheses and commas. It’s widely used for exchanging phylogenetic trees
            between programs.
        </p>
        <p className="mb-4">
            The <strong>Neighbor-Joining algorithm</strong> is a method for building
            trees by clustering sequences based on their pairwise distances. It’s fast
            and commonly used to infer evolutionary relationships.
        </p>
        </section>
    </div>
  );
}
