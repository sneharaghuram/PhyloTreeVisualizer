"use client";

import { useEffect, useMemo, useRef } from "react";
import * as d3 from "d3";
import { NewickNode, parseNewick } from "@/lib/newick";

type Props = {
  newick: string;
  width?: number;
  height?: number;
  margin?: { top: number; right: number; bottom: number; left: number };
  fontSize?: number;
};

type HNode = d3.HierarchyNode<NewickNode> & {
  x: number;
  y: number;
  data: NewickNode;
};

// Compute cumulative branch lengths for x-scaling
function assignCumulativeLength(node: HNode, parentLen = 0) {
  const myLen = (node.data.length ?? 0) + parentLen;
  // store in y temporarily; we'll re-map later
  // but we’ll keep D3 cluster’s x for vertical spacing
  (node as any).cum = myLen;
  node.children?.forEach((c) => assignCumulativeLength(c as HNode, myLen));
}

export default function PhyloTree({
  newick,
  width = 800,
  height = 400,
  margin = { top: 20, right: 180, bottom: 20, left: 120 },
  fontSize = 12,
}: Props) {
  const ref = useRef<SVGSVGElement | null>(null);

  const root = useMemo(() => {
    const parsed = parseNewick(newick);
    // Build a hierarchy (D3 expects children[])
    const h = d3.hierarchy<NewickNode>(parsed, (d) => d.children);
    return h;
  }, [newick]);

  useEffect(() => {
    if (!ref.current || !root) return;

    const innerW = width - margin.left - margin.right;
    const innerH = height - margin.top - margin.bottom;

    // Layout: use cluster for tidy leaf spacing (y is horizontal by default for dendrograms,
    // but we'll handle x/y manually to support branch-length scaling).
    const cluster = d3.cluster<NewickNode>().size([innerH, innerW]);

    // D3’s cluster uses node.y as layout depth. We'll compute our own 'y' using branch lengths.
    const rootCopy = root.copy();
    cluster(rootCopy as any);

    // Compute cumulative branch lengths (distance from root)
    assignCumulativeLength(rootCopy as HNode);

    // Find max cumulative distance for scaling x-axis
    let maxCum = 0;
    rootCopy.each((n: any) => {
      if (n.cum != null && n.cum > maxCum) maxCum = n.cum;
    });
    const xScale = d3.scaleLinear().domain([0, maxCum || 1]).range([0, innerW]);

    // Create svg/g
    const svg = d3.select(ref.current);
    svg.selectAll("*").remove();

    const g = svg
      .attr("width", width)
      .attr("height", height)
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    // Convert cluster x (vertical position) + scaled distance to coordinates:
    // We'll render as a left-to-right tree: x = scaled distance, y = cluster x
    type P = { x: number; y: number };
    function nodePoint(n: any): P {
      return { x: xScale(n.cum ?? 0), y: n.x };
    }

    // Classic elbow link generator
    const links = rootCopy.links();
    // Classic rectangular cladogram elbow: vertical first, then horizontal
    function elbowVH(d: any) {
    const s = nodePoint(d.source); // parent
    const t = nodePoint(d.target); // child
    return `M${s.x},${s.y} V${t.y} H${t.x}`;
    }

    g.append("g")
    .selectAll("path")
    .data(links)
    .join("path")
    .attr("fill", "none")
    .attr("stroke", "#555")
    .attr("stroke-width", 1.2)
    .attr("stroke-linecap", "square")
    .attr("d", elbowVH);


    // Nodes
    const nodes = rootCopy.descendants();
    const node = g
      .append("g")
      .selectAll("g")
      .data(nodes)
      .join("g")
      .attr("transform", (d: any) => {
        const p = nodePoint(d);
        return `translate(${p.x},${p.y})`;
      });

    node
      .append("circle")
      .attr("r", 2.5)
      .attr("fill", (d) => (d.children ? "#555" : "#111"));

    // Labels for leaves and named internals
    node
    .append("text")
    .attr("dy", "0.32em")
    .attr("x", (d) => (d.children ? -8 : 8))
    .attr("text-anchor", (d) => (d.children ? "end" : "start"))
    .style("font-size", `${fontSize}px`)
    .text((d) => {
        // Show only leaf names, or skip "Inner#"
        if (!d.children && d.data.name) return d.data.name;
        return ""; // no label for internal nodes
    });


    // Axis for branch length scale (optional, below the tree)
    const axis = d3.axisBottom(xScale).ticks(6);
    g
      .append("g")
      .attr("transform", `translate(0,${innerH})`)
      .call(axis)
      .call((gAxis) =>
        gAxis
          .append("text")
          .attr("x", innerW)
          .attr("y", 30)
          .attr("fill", "#000")
          .attr("text-anchor", "end")
          .style("font-size", `${fontSize}px`)
          .text("branch length")
      );

  }, [root, width, height, margin, fontSize]);

  return <svg ref={ref} role="img" aria-label="Phylogenetic tree" />;
}
