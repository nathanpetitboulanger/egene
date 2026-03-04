import requests

UNKNOWN_SEQ = (
    "TATACTATATTTTATATTTGCTATATGATCGGGAATAGTTGGAGCTTCTTTAAGTATAATTATTCGTTTAGAATTA"
    "AGTGCTCCAGGAAGATGAATTAATAATGATCAAATTTATAATACTATTGTAACTTCTCATGCTTTTGTTATAATT"
    "TTTTTTATAGTTATACCATTTATAATTGGAGGATTCGGAAATTGACTAGTACCATTAATAATTGGAGCCCCTGAT"
    "ATAGCATTTCCACGAATAAATAATATAAGATTTTGACTATTAATACCATCTTTATTTATACTTATTATAAGAGTAA"
    "CTTTATCAACAGGTTCAGGAACAGGTTGAACTATATATCCTCCTTTATCATCAATTATATATCACTCATCCTATGC"
    "TATAGATTTTACTATTTTTTCTCTACATATTGCAGGAATTTCTTCAATTATAGGAGCTATTAATTTTATTGTTTCAA"
    "TTTTATTAATGAAAAATATTTCACTTAATTTAAATCAAATTCCTTTATTCCCATGATCAGTTAAAATTACTGCAAT"
    "TTTATTACTATTATCTTTACCTGTCTTAGCAGGAGCTATCACTATATTATTAACTGACCGAAATTTAAATACTTCC"
    "TTTTTTGACCCCTCAGGAGGAGGAGATCCTATTTTATATCAACATTTATTT"
)


def identify_sequence(sequence):
    print("Envoi de la séquence à l'API BOLD (Ids_Service)...")
    url = "http://www.boldsystems.org/index.php/Ids_Service"
    params = {"db": "COX1_L640bp", "sequence": sequence}

    try:
        response = requests.get(url, params=params, timeout=30)
        if response.status_code == 200:
            print("\n--- RÉSULTATS D'IDENTIFICATION ---")
            lines = response.text.split("\n")
            for line in lines[:5]:
                if line.strip():
                    print(line)
        else:
            print(
                f"Erreur HTTP : {response.status_code}. L'API est peut-être temporairement indisponible."
            )
    except Exception as e:
        print(f"Erreur : {e}")


if __name__ == "__main__":
    identify_sequence(UNKNOWN_SEQ)
