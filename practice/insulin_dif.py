from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

Entrez.email = "email"


def extract_insulin_chains(id_seq):
    with Entrez.efetch(db="protein", id=id_seq, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    chaines = []
    for feature in record.features:
        if feature.type == "mat_peptide":
            sequence_peptide = feature.location.extract(record.seq)
            nom_peptide = feature.qualifiers.get("product", ["Inconnu"])[0]
            chaines.append((nom_peptide, sequence_peptide))

    return chaines


id_humain = "NP_000198.1"
id_porc = "NP_001103242.1"

chaines_h = extract_insulin_chains(id_humain)
chaines_p = extract_insulin_chains(id_porc)

print(chaines_h)
print(chaines_p)

aligner = PairwiseAligner()

aligner.substitution_matrices = substitution_matrices.load("BLOSUM62")

for i in range(len(chaines_h)):
    nom_h, seq_h = chaines_h[i]
    nom_p, seq_p = chaines_p[i]

    print(f"\n--- Alignement : {nom_h} ---")
    alignement = aligner.align(seq_h, seq_p)[0]
    print(alignement)
