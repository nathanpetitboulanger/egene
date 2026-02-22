# Test des méthodes d'importation
from Bio.Seq import Seq

adn_brut = "ATGCGTGCCTCGATC"
ma_sequence = Seq(adn_brut)

print(f"Voici ma séquence :{ma_sequence} ")

print(f"Voici sa longueur {len(ma_sequence)}")

ma_seq_rev = ma_sequence.reverse_complement()
print(f"Et voici son complément inverse {ma_seq_rev}")

ma_seq_transcrite = ma_sequence.transcribe()
print(f"Et voici son ARNm {ma_seq_transcrite}")

ma_seq_prot = ma_sequence.translate()
print(f"Et voici sa protéine {ma_seq_prot}")

######################################################

from Bio import Entrez, SeqIO

Entrez.email = "email"

handle = Entrez.efetch(db="nucleotide", id="NM_000207", rettype="gb", retmode="text")

record = SeqIO.read(handle, "genbank")

handle.close()

print(f"Id du gène : {record.id}")
print(f"Description : {record.description}")
print(f"Longueur du gène : {len(record.seq)}")

print(f"Séquence du gène : {record.seq}")
print(f"Et voici les 10 premiers an : {record.seq[:10]}")

insulin_transcript = record.seq.transcribe()
print(insulin_transcript)

insulin_prot = record.seq.translate()
print(insulin_prot)


######################################################

print(record.features)

for feature in record.features:
    if feature.type == "CDS":
        print(f"Voici le type : {feature.type}")
        print(f"La position : {feature.location}")

        gene_name = feature.qualifiers.get("gene", ["Inconnu"])[0]
        print(f"Nom : {gene_name}")

        traduction = feature.qualifiers.get("translation", [""])[0]
        print(f"traduction : {traduction[:60]}")


######################################################

nom_sortie = "fasta_test"

count = SeqIO.write([record], nom_sortie, "fasta")

SeqIO.write([record], "insulin_fasta", "fasta")
