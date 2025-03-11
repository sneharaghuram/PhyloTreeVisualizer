"use client";
import { useState } from "react";
import axios from "axios";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";

// Define valid phylogenetic model types
type PhyloModel = "nj" | "upgma" | "ml";

// Define API response type
interface ApiResponse {
  result: string;
}

export default function Home() {
  const [file, setFile] = useState<File | null>(null);
  const [model, setModel] = useState<PhyloModel>("nj");
  const [loading, setLoading] = useState<boolean>(false);
  const [response, setResponse] = useState<string>("");

  // Handle file selection
  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>): void => {
    if (e.target.files && e.target.files.length > 0) {
      setFile(e.target.files[0]);
    }
  };

  // // Handle form submission
  // const handleSubmit = async (): Promise<void> => {
  //   if (!file) {
  //     alert("Please upload a FASTA file.");
  //     return;
  //   }

  //   setLoading(true);
  //   const formData = new FormData();
  //   formData.append("file", file);
  //   formData.append("model", model);

  //   try {
  //     const res = await axios.post<ApiResponse>("/api/analyze", formData);
  //     setResponse(res.data.result);
  //   } catch (error) {
  //     console.error("Error processing file:", error);
  //     setResponse("Error processing the file.");
  //   } finally {
  //     setLoading(false);
  //   }
  // };

  return (
    <div className="flex flex-col items-center min-h-screen p-10 bg-gray-100">
      <Card className="w-full max-w-md p-4">
        <CardHeader>
          <CardTitle>Phylogenetic Tree Generator</CardTitle>
        </CardHeader>
        <CardContent>
          {/* File Upload */}
          <Input type="file" accept=".fasta" />

          {/* Model Selection Dropdown */}
          <select
            className="w-full p-2 mt-4 border rounded"
            value={model}
            onChange={(e) => setModel(e.target.value as PhyloModel)}
          >
            <option value="nj">Neighbor-Joining</option>
            <option value="upgma">UPGMA</option>
            <option value="ml">Maximum Likelihood</option>
          </select>

          {/* Submit Button */}
          <Button className="w-full mt-4">
            {loading ? "Processing..." : "Generate Tree"}
          </Button>

          {/* Response Output */}
          {response && <p className="mt-4 text-green-600">{response}</p>}
        </CardContent>
      </Card>
    </div>
  );
}

