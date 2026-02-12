from os import access
from Bio import Entrez, SeqIO


Entrez.email = "nath39.np@gmail.com"


# On récupère deux versions d'une même protéine (exemples d'IDs NCBI)
# Séquence A (Référence) et Séquence B (Variant)
from Bio import Entrez, SeqIO, Align
import matplotlib.pyplot as plt
import numpy as np

ids = ["P01308", "P01325"]
sequences = []

for access_id in ids:
    handle = Entrez.efetch(
        db="protein",
        id=access_id,
        rettype="fasta",
        retmode="text",
    )
    record = SeqIO.read(handle, "fasta")
    sequences.append(record)
    handle.close()


handle = Entrez.efetch(
    db="protein",
    id=ids[0],
    rettype="fasta",
    retmode="text",
)
record = SeqIO.read(handle, "fasta")
sequences.append(record)
handle.close()

seq = sequences
