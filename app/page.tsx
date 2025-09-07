"use client";
import { useState } from "react";
import axios from "axios";
import Button from "@/components/ui/button";
import Input from "@/components/ui/input";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import dynamic from "next/dynamic";

// dynamic to avoid SSR issues with window/D3
const PhyloTree = dynamic(() => import("@/components/PhyloTree"), { ssr: false });

type PhyloMethod = "nj" | "upgma"; // ML not implemented on backend yet

interface ApiResponse {
  tree_newick: string;
  taxa: string[];
  num_taxa: number;
  distance_names: string[];
  method: string;
  distance_model: string;
  error?: string;
}

export default function Home() {
  const [file, setFile] = useState<File | null>(null);
  const [method, setMethod] = useState<PhyloMethod>("nj");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<ApiResponse | null>(null);

  const API_BASE =
    process.env.NEXT_PUBLIC_API_BASE?.replace(/\/+$/, "") || "http://127.0.0.1:8000";

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    setResult(null);
    if (e.target.files && e.target.files.length > 0) {
      setFile(e.target.files[0]);
    }
  };

  const handleSubmit = async () => {
    if (!file) return alert("Please upload a FASTA file.");

    setLoading(true);
    setResult(null);

    const formData = new FormData();
    formData.append("file", file);
    formData.append("method", method); // backend expects "method"

    try {
      const { data } = await axios.post<ApiResponse>(`${API_BASE}/analyze/`, formData);
      setResult(data);
    } catch (err: any) {
      const msg =
        err?.response?.data?.error ||
        err?.message ||
        "Error processing the file.";
      setResult({
        tree_newick: "",
        taxa: [],
        num_taxa: 0,
        distance_names: [],
        method,
        distance_model: "identity",
        error: msg,
      });
      console.error("Analyze error:", err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="flex flex-col items-center min-h-screen p-10 bg-gray-100">
      <Card className="w-full max-w-4xl p-4"> {/* widened for SVG */}
        <CardHeader>
          <CardTitle>Phylogenetic Tree Generator</CardTitle>
        </CardHeader>
        <CardContent>
          {/* File Upload */}
          <Input type="file" accept=".fasta,.fa,.txt" onChange={handleFileChange} />

          {/* Method Selection */}
          <select
            className="w-full p-2 mt-4 border rounded"
            value={method}
            onChange={(e) => setMethod(e.target.value as PhyloMethod)}
          >
            <option value="nj">Neighbor-Joining</option>
            <option value="upgma">UPGMA</option>
          </select>

          {/* Submit Button */}
          <Button className="w-full mt-4" onClick={handleSubmit} disabled={loading}>
            {loading ? "Processing..." : "Generate Tree"}
          </Button>

          {/* Response */}
          {result && (
            <div className="mt-4 space-y-4">
              {result.error ? (
                <p className="text-red-600">{result.error}</p>
              ) : (
                <>
                  {result.taxa?.length > 0 && (
                    <p className="text-sm">
                      <span className="font-semibold">Taxa:</span>{" "}
                      {result.taxa.join(", ")}
                    </p>
                  )}

                  {/* Newick (keep for debugging) */}
                  {result.tree_newick && (
                    <>
                      <p className="text-sm font-semibold">Newick:</p>
                      <pre className="bg-white rounded p-3 text-xs overflow-x-auto">
                        {result.tree_newick}
                      </pre>
                    </>
                  )}

                  {/* 🔽 D3 Tree render (this is the block you asked to add) */}
                  {result.tree_newick && (
                    <>
                      <h3 className="mt-6 font-semibold">Tree</h3>
                      <div className="border rounded bg-white p-2 overflow-auto">
                        <PhyloTree newick={result.tree_newick} width={900} height={500} />
                      </div>
                    </>
                  )}
                </>
              )}
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  );
}
