import pandas as pd
import requests
import io

# Dataset ID specified in the TD
DATASET_ID = "DS-BBBABSV"

def fetch_and_explore_dataset(dataset_id):
    # Important: Use HTTPS and specify the format
    url = f"https://www.boldsystems.org/index.php/API_Public/combined?dataset={dataset_id}&format=tsv"
    
    # Adding a User-Agent header is often necessary for BOLD
    headers = {
        'User-Agent': 'Mozilla/5.0 (compatible; BOLD-Python-Client; +https://github.com/nathan/egene)'
    }
    
    print(f"Tentative de téléchargement du dataset {dataset_id} depuis BOLD (HTTPS)...")
    
    try:
        # allow_redirects=True is default in requests, but let's be explicit
        response = requests.get(url, headers=headers, timeout=60, allow_redirects=True)
        
        if response.status_code == 200:
            if not response.text.strip():
                print("L'API a répondu avec succès mais le contenu est vide. Le dataset est peut-être privé.")
                return None
                
            # Parse TSV data
            data = pd.read_csv(io.StringIO(response.text), sep='\t')
            
            print("\n--- STATISTIQUES DU JEU DE DONNÉES RÉEL ---")
            print(f"Nombre de spécimens : {len(data)}")
            
            # Check for necessary columns
            if 'species_name' in data.columns:
                print(f"Nombre d'espèces : {data['species_name'].nunique()}")
            if 'nucleotides' in data.columns:
                print(f"Nombre de séquences (COI-5P) : {data['nucleotides'].notnull().sum()}")
            
            return data
        else:
            print(f"Erreur HTTP {response.status_code} sur {url}")
            return None
    except Exception as e:
        print(f"Erreur lors de la requête : {e}")
        return None

if __name__ == "__main__":
    df = fetch_and_explore_dataset(DATASET_ID)
    if df is not None:
        # Save to local file
        output_path = "barcoding_2/td_python/bold_dataset_bbb.tsv"
        df.to_csv(output_path, sep='\t', index=False)
        print(f"\nSuccès ! Données réelles sauvegardées dans : {output_path}")
    else:
        print("\nImpossible de récupérer les données réelles. Veuillez vérifier si le dataset est public ou si les serveurs BOLD sont accessibles.")
