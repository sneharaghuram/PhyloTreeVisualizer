export default function Footer() {
  return (
    <footer className="fixed bottom-0 left-0 w-full bg-lime-800 text-white py-4 px-6 text-sm shadow-inner">
      <div className="max-w-4xl mx-auto flex flex-col sm:flex-row justify-between items-center space-y-2 sm:space-y-0">
        <p>&copy; {new Date().getFullYear()} PhyloTree Visualizer</p>
        <nav className="space-x-4">
          <a
            href="https://github.com/sneharaghuram/PhyloTreeVisualizer"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:text-lime-400 transition-colors"
          >
            GitHub
          </a>
          <a href="/aboutme" className="hover:text-lime-400 transition-colors">
            About Me
          </a>
          <a
            href="mailto:sneha.raghuram3@gmail.com"
            className="hover:text-lime-400 transition-colors"
          >
            Contact
          </a>
        </nav>
      </div>
    </footer>
  );
}
