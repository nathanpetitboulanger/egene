import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os

def calculate_simple_distance(seq1, seq2):
    """Calculer la distance de Hamming simple (approximation rapide du modèle Kimura 2P)"""
    if pd.isna(seq1) or pd.isna(seq2): return np.nan
    # Nettoyage minimal (caractères non-bases)
    s1 = str(seq1).strip().upper().replace("-", "")
    s2 = str(seq2).strip().upper().replace("-", "")
    min_len = min(len(s1), len(s2))
    if min_len == 0: return np.nan
    mismatches = sum(1 for a, b in zip(s1[:min_len], s2[:min_len]) if a != b)
    return (mismatches / min_len) * 100

def run_barcode_gap_analysis(file_path):
    print(f"Chargement des données réelles depuis : {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    
    # Filtrage des séquences nulles
    df = df[df['nucleotides'].notnull()]
    
    # Filtrer les espèces avec au moins 2 spécimens pour le calcul intra
    species_counts = df['species_name'].value_counts()
    valid_species = species_counts[species_counts >= 2].index.tolist()
    df_filtered = df[df['species_name'].isin(valid_species)]
    
    print(f"Analyse de {len(df_filtered)} spécimens pour {len(valid_species)} espèces.")
    
    intra_distances = []
    inter_distances = []
    
    # 1. Distances intra-spécifiques
    for species in valid_species:
        group = df_filtered[df_filtered['species_name'] == species]['nucleotides'].tolist()
        for s1, s2 in combinations(group, 2):
            d = calculate_simple_distance(s1, s2)
            if not np.isnan(d): intra_distances.append(d)
            
    # 2. Distances inter-spécifiques (échantillonnage des paires d'espèces)
    species_list = df_filtered['species_name'].unique()
    for sp1, sp2 in combinations(species_list, 2):
        s1 = df_filtered[df_filtered['species_name'] == sp1]['nucleotides'].iloc[0]
        s2 = df_filtered[df_filtered['species_name'] == sp2]['nucleotides'].iloc[0]
        d = calculate_simple_distance(s1, s2)
        if not np.isnan(d): inter_distances.append(d)

    # 3. Tracé du graphique
    if intra_distances or inter_distances:
        plt.figure(figsize=(10, 6))
        plt.hist(intra_distances, bins=30, alpha=0.5, label='Intra-spécifique', color='blue')
        plt.hist(inter_distances, bins=30, alpha=0.5, label='Inter-spécifique', color='red')
        plt.axvline(x=2.0, color='green', linestyle='--', label='Seuil typique (2%)')
        plt.title('Distribution des distances génétiques (Données réelles)')
        plt.xlabel('Distance génétique (%)')
        plt.ylabel('Fréquence')
        plt.legend()
        plt.grid(alpha=0.3)
        
        output_img = "barcoding_2/td_python/barcode_gap_plot.png"
        plt.savefig(output_img)
        print(f"\nAnalyse terminée. Graphique sauvegardé dans : {output_img}")
        
        # Statistiques
        print("\n--- RÉSULTATS BARCODE GAP ---")
        if intra_distances:
            print(f"Distance intra moyenne : {np.mean(intra_distances):.2f}% (max: {np.max(intra_distances):.2f}%)")
        if inter_distances:
            print(f"Distance inter moyenne : {np.mean(inter_distances):.2f}% (min: {np.min(inter_distances):.2f}%)")
    else:
        print("Erreur : Aucune distance n'a pu être calculée. Vérifiez le format du fichier TSV.")

if __name__ == "__main__":
    file_path = "barcoding_2/td_python/bold_dataset_bbb.tsv"
    if os.path.exists(file_path):
        run_barcode_gap_analysis(file_path)
    else:
        print(f"ERREUR : Fichier {file_path} introuvable.")
        print("Veuillez télécharger le dataset DS-BBBABSV manuellement depuis le Workbench BOLD et l'enregistrer à cet emplacement.")
