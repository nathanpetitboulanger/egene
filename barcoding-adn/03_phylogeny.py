import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import math


def calculate_k2p_distance(seq1, seq2):
    """
    Calcule la distance de Kimura 2-Paramètres (K2P) entre deux séquences.
    K2P prend en compte le fait que les transitions (A<->G, C<->T)
    sont plus fréquentes que les transversions (A<->C, A<->T, G<->C, G<->T).
    """
    n = 0
    transitions = 0
    transversions = 0

    purines = ["A", "G"]
    pyrimidines = ["C", "T"]

    for b1, b2 in zip(seq1, seq2):
        if b1 == "-" or b2 == "-" or b1 == "N" or b2 == "N":
            continue
        n += 1
        if b1 != b2:
            if (b1 in purines and b2 in purines) or (
                b1 in pyrimidines and b2 in pyrimidines
            ):
                transitions += 1
            else:
                transversions += 1

    if n == 0:
        return 0

    P = transitions / n
    Q = transversions / n

    try:
        val1 = 1 - 2 * P - Q
        val2 = 1 - 2 * Q
        if val1 <= 0 or val2 <= 0:
            return 1.0
        dist = -0.5 * math.log(val1 * math.sqrt(val2))
        return dist
    except ValueError:
        return 1.0


def build_tree(file_path):
    print("Chargement des données pour la phylogénie...")
    df = pd.read_csv(file_path, sep="\t", low_memory=False)
    df = df.dropna(subset=["nucleotides"])

    # On utilise Species_name + ProcessID pour identifier les branches
    # On nettoie un peu le nom d'espèce pour l'affichage
    df["label"] = df["species_name"].fillna("Unknown") + "_" + df["processid"]
    labels = df["label"].tolist()
    sequences = df["nucleotides"].tolist()

    n = len(labels)
    print(f"Calcul de la matrice de distance pour {n} séquences...")

    matrix = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                dist = calculate_k2p_distance(sequences[i], sequences[j])
                row.append(dist)
        matrix.append(row)

    dm = DistanceMatrix(names=labels, matrix=matrix)

    print("Construction de l'arbre Neighbor-Joining...")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    Phylo.write(tree, "bumblebee_tree.nwk", "newick")
    print("Arbre sauvegardé dans bumblebee_tree.nwk")

    print("\\n--- Aperçu de l'Arbre (NJ) ---")
    Phylo.draw_ascii(tree)


if __name__ == "__main__":
    build_tree("bold_data.tsv")

    print("\\nConcepts clés :")
    print("1. K2P (Kimura 2-Paramètres) : Un modèle d'évolution de l'ADN.")
    print("2. Neighbor-Joining (NJ) : Méthode de clustering pour arbres.")
    print("3. Topologie : Permet d'identifier les erreurs d'identification.")

