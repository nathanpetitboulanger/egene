import os
from Bio import Entrez, SeqIO, Align

# --- CONFIGURATION ---
Entrez.email = "nath39.np@gmail.com"

print("=== ÉTAPE 1 : Récupération des données via NCBI (Entrez) ===")

ids = {"Humain": "P01308", "Souris": "P01325"}
sequences = {}

for espece, accession in ids.items():
    print(f"Téléchargement de la séquence FASTA pour : {espece} ({accession})...")
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    sequences[espece] = record
    handle.close()

print("\n=== ÉTAPE 2 : Traitement et Exploration de base ===")

for espece, seq_record in sequences.items():
    print(f"[{espece}] Longueur : {len(seq_record.seq)} AA")

print("\n=== ÉTAPE 3 : Alignement Global ===")

aligner = Align.PairwiseAligner()
aligner.mode = "global"
human_seq = sequences["Humain"].seq
mouse_seq = sequences["Souris"].seq

alignments = aligner.align(human_seq, mouse_seq)
best_alignment = alignments[0]
print(best_alignment)

matches = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a == b)
identity = (matches / len(human_seq)) * 100
print(f"Identité globale : {identity:.2f}%")

print("\n=== ÉTAPE 4 : Extraction des régions (GenBank) ===")

regions_sequences = {"Humain": {}, "Souris": {}}

for espece, accession in ids.items():
    print(f"Extraction des annotations pour {espece}...")
    handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in record.features:
        # On accepte 'Chain', 'mat_peptide' ou 'Region'
        if feature.type in ["Chain", "mat_peptide", "Region"]:
            name = feature.qualifiers.get(
                "product", feature.qualifiers.get("note", ["Unknown"])
            )[0]
            regions_sequences[espece][name] = feature.extract(record.seq)


# Vérification du contenu
print(f"Régions Humain : {list(regions_sequences['Humain'].keys())}")
print(f"Régions Souris : {list(regions_sequences['Souris'].keys())}")


print("\n=== ÉTAPE 5 : Comparaison ciblée (Chaîne A vs B) ===")

for region_nom in ["Insulin A chain", "Insulin B chain"]:
    h_key = next(
        (
            k
            for k in regions_sequences["Humain"].keys()
            if region_nom.lower() in k.lower()
        ),
        None,
    )
    m_key = next(
        (
            k
            for k in regions_sequences["Souris"].keys()
            if region_nom.lower() in k.lower()
        ),
        None,
    )

    if h_key and m_key:
        h_seq = regions_sequences["Humain"][h_key]
        m_seq = regions_sequences["Souris"][m_key]
        reg_alignments = aligner.align(h_seq, m_seq)
        reg_id = (
            sum(1 for a, b in zip(reg_alignments[0][0], reg_alignments[0][1]) if a == b)
            / len(h_seq)
        ) * 100
        print(f"Identité [{region_nom}] : {reg_id:.1f}%")

print("\n=== ÉTAPE 6 : Le cas du Peptide C (moins conservé) ===")

# Le Peptide C est la partie centrale qui est découpée pour libérer l'insuline active.
# On va voir s'il a autant "résisté" à l'évolution que les chaînes A et B.
region_c = "peptide"  # Souvent nommé 'Insulin C-peptide' ou juste 'Peptide'
h_key_c = next(
    (
        k
        for k in regions_sequences["Humain"].keys()
        if "c-peptide" in k.lower() or "peptide" in k.lower()
    ),
    None,
)
m_key_c = next(
    (
        k
        for k in regions_sequences["Souris"].keys()
        if "c-peptide" in k.lower() or "peptide" in k.lower()
    ),
    None,
)

if h_key_c and m_key_c:
    h_seq_c = regions_sequences["Humain"][h_key_c]
    m_seq_c = regions_sequences["Souris"][m_key_c]
    reg_alignments_c = aligner.align(h_seq_c, m_seq_c)
    reg_id_c = (
        sum(1 for a, b in zip(reg_alignments_c[0][0], reg_alignments_c[0][1]) if a == b)
        / len(h_seq_c)
    ) * 100
    print(f"Identité [Peptide C] : {reg_id_c:.1f}%")
    print(
        "Observation : Le peptide C est souvent moins conservé que A et B car sa séquence précise"
    )
    print("est moins cruciale pour la survie de l'individu.")

print("\n=== ÉTAPE 7 : Passage à l'ADN (mRNA) ===")

# Pour comprendre les mutations "silencieuses", il faut regarder l'ADN (ou l'ARNm).
# NM_000618 : ARNm Insuline humaine
# NM_008386 : ARNm Insuline-1 souris
ids_adn = {"Humain": "NM_000618", "Souris": "NM_008386"}
sequences_adn = {}

for espece, accession in ids_adn.items():
    print(f"Téléchargement de l'ARNm pour : {espece} ({accession})...")
    handle = Entrez.efetch(
        db="nucleotide", id=accession, rettype="fasta", retmode="text"
    )
    record = SeqIO.read(handle, "fasta")
    sequences_adn[espece] = record
    handle.close()

# Alignement des séquences d'ADN
adn_human = sequences_adn["Humain"].seq
adn_mouse = sequences_adn["Souris"].seq

# Pour l'ADN, on peut utiliser des scores différents (Match=2, Mismatch=-3, etc.)
# mais on reste simple ici avec le PairwiseAligner par défaut.
alignments_adn = aligner.align(adn_human, adn_mouse)
best_adn = alignments_adn[0]

# Calcul de l'identité nucléotidique
matches_adn = sum(1 for a, b in zip(best_adn[0], best_adn[1]) if a == b)
identity_adn = (matches_adn / len(adn_human)) * 100
print(f"Identité au niveau ADN (ARNm) : {identity_adn:.2f}%")

print("\n--- Conclusion finale du Notebook ---")
print("Nous avons vu que :")
print("1. La protéine est très conservée (pression de sélection).")
print("2. Certaines régions (A et B) le sont plus que d'autres (C).")
print("3. L'identité au niveau ADN est souvent plus faible qu'au niveau protéine")
print("   à cause de la redondance du code génétique (mutations synonymes).")
