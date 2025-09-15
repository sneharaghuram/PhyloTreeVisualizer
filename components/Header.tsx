"use client";
import Link from "next/link";

export default function Header() {
  return (
    <header className="fixed top-0 left-0 w-full bg-lime-800 text-white px-6 py-3 flex justify-between items-center shadow-md z-50">
      <h1 className="text-lg font-bold">PhyloTree Visualizer</h1>
      <nav className="space-x-6">
        <Link href="/" className="hover:text-lime-400 transition-colors">
          Generate Tree
        </Link>
        <Link href="/about" className="hover:text-lime-400 transition-colors">
          About
        </Link>
      </nav>
    </header>
  );
}
