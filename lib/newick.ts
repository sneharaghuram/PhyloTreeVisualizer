// lib/newick.ts
export type NewickNode = {
  name?: string;
  length?: number;
  children?: NewickNode[];
};

export function parseNewick(newick: string): NewickNode {
  // Based on a small, classic recursive descent parser
  let i = 0;

  function readName(): string {
    let s = "";
    while (i < newick.length) {
      const c = newick[i];
      if (",():;".includes(c)) break;
      s += c;
      i++;
    }
    return s.trim();
  }

  function readLength(): number | undefined {
    if (newick[i] !== ":") return undefined;
    i++; // skip :
    const start = i;
    while (i < newick.length && !",():;".includes(newick[i])) i++;
    const num = newick.slice(start, i).trim();
    const v = Number(num);
    return isFinite(v) ? v : undefined;
  }

  function readSubtree(): NewickNode {
    const node: NewickNode = {};
    if (newick[i] === "(") {
      i++; // skip '('
      node.children = [];
      while (true) {
        node.children.push(readSubtree());
        if (newick[i] === ",") {
          i++;
          continue;
        }
        if (newick[i] === ")") {
          i++;
          break;
        }
        throw new Error("Unexpected character in Newick near: " + newick.slice(i, i + 20));
      }
      // optional name after children
      const nm = readName();
      if (nm) node.name = nm;
      node.length = readLength();
      return node;
    } else {
      // leaf
      const nm = readName();
      if (nm) node.name = nm;
      node.length = readLength();
      return node;
    }
  }

  const root = readSubtree();
  if (newick[i] === ";") i++;
  return root;
}
