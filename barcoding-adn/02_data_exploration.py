import pandas as pd


def explore_bold_data(file_path):
    print(f"Analyse du fichier : {file_path}\\n")

    # Chargement des données. On utilise le séparateur tabulation (\t) car c'est un TSV.
    # On spécifie low_memory=False pour éviter les warnings sur les types de colonnes.
    df = pd.read_csv(file_path, sep="\t", low_memory=False)

    # --- Question 1 : Description du dataset ---
    print("--- Question 1 : Description du Dataset ---")

    num_specimens = len(df)
    # On compte les séquences présentes (colonne 'nucleotides' non vide)
    num_sequences = df["nucleotides"].dropna().count()
    # Nombre d'espèces uniques
    num_species = df["species_name"].nunique()

    print(f"Nombre de spécimens : {num_specimens}")
    print(f"Nombre de séquences : {num_sequences}")
    print(f"Nombre d'espèces : {num_species}")

    # Origine géographique (Pays)
    countries = df["country"].unique()
    print(
        f"Origines géographiques : {', '.join([str(c) for c in countries if pd.notna(c)])}"
    )

    # Marqueur génétique (souvent COI-5P pour le barcoding animal)
    marker = df["markercode"].unique()
    print(
        f"Marqueur génétique utilisé : {', '.join([str(m) for m in marker if pd.notna(m)])}"
    )

    # Longueur moyenne des séquences
    seq_lengths = df["nucleotides"].dropna().str.len()
    print(f"Longueur moyenne des fragments : {seq_lengths.mean():.1f} bp")
    print(f"Longueur min/max : {seq_lengths.min()} / {seq_lengths.max()} bp")

    # --- Question 2 : Champs et colonnes ---
    print("\\n--- Question 2 : Champs et Signification ---")
    print("Principales colonnes présentes :")
    important_cols = [
        "processid",
        "sampleid",
        "recordset",
        "species_name",
        "bin_uri",
        "country",
        "nucleotides",
    ]
    for col in important_cols:
        if col in df.columns:
            print(f"- {col}")

    # Pourquoi certaines entrées n'ont pas de séquence ?
    missing_seq = num_specimens - num_sequences
    print(f"\\nSpécimens sans séquence : {missing_seq}")

    print("\\nConcepts clés :")
    print("1. processid : Identifiant unique du processus de barcoding chez BOLD.")
    print("2. bin_uri : Code BIN (Barcode Index Number). C'est un système qui")
    print("   regroupe les séquences similaires pour suggérer des espèces,")
    print("   indépendamment des noms taxonomiques classiques.")
    print("3. Pourquoi pas de séquence ? Cela peut être dû à un échec de l'extraction")
    print("   de l'ADN, une dégradation de l'échantillon (trop vieux), ou un échec")
    print("   de l'amplification PCR.")


if __name__ == "__main__":
    explore_bold_data("bold_data.tsv")
