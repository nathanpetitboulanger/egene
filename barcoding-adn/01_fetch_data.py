import requests
import os


def download_bold_dataset(dataset_code, output_file):
    """
    Télécharge un dataset complet depuis BOLD Systems au format TSV.
    Le format TSV de BOLD contient les métadonnées et les séquences.
    """
    print(f"Début du téléchargement du dataset {dataset_code}...")

    # URL de l'API BOLD pour récupérer les données combinées (spécimens + séquences)
    url = f"http://v4.boldsystems.org/index.php/API_Public/combined?container={dataset_code}&format=tsv"

    try:
        response = requests.get(url)
        # On vérifie si la requête a réussi (code 200)
        response.raise_for_status()

        # On enregistre le contenu dans un fichier
        with open(output_file, "wb") as f:
            f.write(response.content)

        print(f"Succès ! Fichier enregistré sous : {output_file}")

        # Petite vérification de la taille du fichier
        file_size = os.path.getsize(output_file) / 1024
        print(f"Taille du fichier : {file_size:.2f} KB")

    except requests.exceptions.RequestException as e:
        print(f"Erreur lors du téléchargement : {e}")


if __name__ == "__main__":
    # Le code du dataset fourni dans les consignes
    DATASET_ID = "DS-BBBABSV"
    OUTPUT_PATH = "bold_data.tsv"

    download_bold_dataset(DATASET_ID, OUTPUT_PATH)

    print("\nConcept : BOLD (Barcode of Life Data Systems) est une plateforme cloud")
    print("qui centralise les séquences d'ADN 'barcodes' (généralement le gène COI)")
    print("associées à des spécimens identifiés par des experts. L'API nous permet")
    print("d'extraire ces données sans passer par l'interface web.")

