import matplotlib.pyplot as plt
import networkx as nx
from Bio.Data import CodonTable

# === Codon → Amino Acid mapping ===
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
codon_to_aa = {c: a for c, a in standard_table.forward_table.items()}
for c in standard_table.stop_codons:
    codon_to_aa[c] = "Stop"

# === Build tree structure ===
G = nx.DiGraph()
root = "Genetic Code"
G.add_node(root)

# Group codons by amino acid
aa_to_codons = {}
for codon, aa in codon_to_aa.items():
    aa_to_codons.setdefault(aa, []).append(codon)

# Add amino acids and codons to graph
for aa, codons in aa_to_codons.items():
    G.add_edge(root, aa)
    for codon in codons:
        G.add_edge(aa, codon)

# === Draw tree ===
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, seed=42)  # layout for visualization
nx.draw(G, pos, with_labels=True, node_size=1500, node_color="lightblue", font_size=8, arrows=False)
plt.title("Amino Acid → Codon Tree (Genetic Code)")
plt.tight_layout()
plt.savefig("amino_acid_codon_tree.png")
plt.close()
