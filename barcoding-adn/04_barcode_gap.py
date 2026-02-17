import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math


def calculate_k2p_distance(seq1, seq2):
    """Calcule la distance K2P entre deux séquences."""
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
            return 0.2
        return -0.5 * math.log(val1 * math.sqrt(val2))
    except ValueError:
        return 0.2


def barcode_gap_analysis(file_path):
    print("Analyse du Barcode Gap...")
    df = pd.read_csv(file_path, sep="\t", low_memory=False).dropna(
        subset=["nucleotides", "species_name"]
    )

    species = df["species_name"].tolist()
    sequences = df["nucleotides"].tolist()
    ids = df["processid"].tolist()
    n = len(species)

    intra_dist = []
    inter_dist = []
    problems = []

    for i in range(n):
        for j in range(i + 1, n):
            dist = calculate_k2p_distance(sequences[i], sequences[j])
            if species[i] == species[j]:
                intra_dist.append(dist)
                if dist > 0.03:
                    problems.append(
                        f"Alerte Intra: {species[i]} ({ids[i]} vs {ids[j]}) dist={dist:.4f}"
                    )
            else:
                inter_dist.append(dist)
                if dist < 0.01:
                    problems.append(
                        f"Alerte Inter: {species[i]} vs {species[j]} ({ids[i]}, {ids[j]}) dist={dist:.4f}"
                    )

    plt.figure(figsize=(10, 6))
    plt.hist(
        inter_dist, bins=50, alpha=0.5, label="Distances Inter-spécifiques", color="red"
    )
    plt.hist(
        intra_dist,
        bins=50,
        alpha=0.5,
        label="Distances Intra-spécifiques",
        color="blue",
    )
    plt.axvline(x=0.02, color="black", linestyle="--", label="Seuil suggéré (2%)")
    plt.title("Barcode Gap Analysis - Bombus dataset")
    plt.xlabel("Distance génétique (K2P)")
    plt.ylabel("Fréquence")
    plt.legend()
    plt.grid(axis="y", alpha=0.3)
    plt.savefig("barcode_gap_plot.png")
    print("Graphique sauvegardé dans barcode_gap_plot.png")

    print("""
--- Alertes de potentiels erreurs d'identification ---""")
    for p in problems[:10]:
        print(p)

    print("""
Concepts clés :
1. Barcode Gap : C'est l'écart entre la variation maximale au sein d'une espèce et la variation minimale entre espèces.
2. Red line : Souvent, une distance > 2% suggère une nouvelle espèce ou une erreur d'identification.
3. Si les distributions se chevauchent, le barcoding seul ne suffit pas.""")


if __name__ == "__main__":
    barcode_gap_analysis("bold_data.tsv")
