import os
from Bio import Entrez, SeqIO

# --- CONFIGURATION ---
# Utilisation de l'email configuré dans les autres scripts du projet
Entrez.email = "nath39.np@gmail.com"

# %% [markdown]
# # Interaction avec le Génome du Riz (Oryza sativa) via NCBI
# Ce notebook (script style notebook) explore les données génomiques du riz.
# Le riz a l'un des plus petits génomes parmi les céréales, ce qui en fait un modèle idéal.

# %% [markdown]
# ## ÉTAPE 1 : Recherche de l'assemblage du génome
# Nous allons chercher les informations sur l'assemblage de référence pour *Oryza sativa Japonica Group*.

# %%
print("Recherche de l'assemblage de référence pour le riz...")
search_handle = Entrez.esearch(
    db="assembly",
    term="Oryza sativa Japonica Group[Organism] AND latest[filter]",
    retmax=5,
)
search_results = Entrez.read(search_handle)
search_handle.close()

ids = search_results["IdList"]
print(f"IDs d'assemblage trouvés : {ids}")

# %% [markdown]
# ## ÉTAPE 2 : Récupération des détails de l'assemblage
# On récupère les résumés pour comprendre quelle version du génome est la plus récente.

# %%
if ids:
    summary_handle = Entrez.esummary(db="assembly", id=ids[0])
    summary_record = Entrez.read(summary_handle)
    summary_handle.close()

    assembly_info = summary_record["DocumentSummarySet"]["DocumentSummary"][0]
    print(f"Nom de l'assemblage : {assembly_info['AssemblyName']}")
    print(f"Statut : {assembly_info['AssemblyStatus']}")
    print(f"Espèce : {assembly_info['SpeciesName']}")
    print(f"Date de dépôt : {assembly_info['AsmReleaseDate_GenBank']}")

# %% [markdown]
# ## ÉTAPE 3 : Recherche d'un gène spécifique (ex: Hexokinase)
# On cherche un gène spécifique pour illustrer la récupération de séquence.

# %%
gene_name = "Oryza sativa hexokinase"
print(f"\nRecherche du gène : {gene_name}...")
gene_search_handle = Entrez.esearch(
    db="nucleotide", term=f"{gene_name}[Title] AND Oryza sativa[Organism]", retmax=1
)
gene_search_results = Entrez.read(gene_search_handle)
gene_search_handle.close()

if gene_search_results["IdList"]:
    gene_id = gene_search_results["IdList"][0]
    print(f"ID du gène trouvé : {gene_id}")

    # Récupération de la séquence au format FASTA
    fetch_handle = Entrez.efetch(
        db="nucleotide", id=gene_id, rettype="fasta", retmode="text"
    )
    gene_record = SeqIO.read(fetch_handle, "fasta")
    fetch_handle.close()

    print(f"\nSéquence récupérée : {gene_record.description}")
    print(f"Longueur de la séquence : {len(gene_record.seq)} bp")
    print(f"Début de la séquence : {gene_record.seq[:100]}...")

# %% [markdown]
# ## ÉTAPE 4 : Analyse de l'annotation (GenBank)
# On récupère le format GenBank pour voir les annotations (exons, introns, etc.)

# %%
if gene_search_results["IdList"]:
    print(f"\nRécupération des annotations GenBank pour l'ID : {gene_id}...")
    gb_handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
    gb_record = SeqIO.read(gb_handle, "genbank")
    gb_handle.close()

    print(f"Nombre de features (annotations) : {len(gb_record.features)}")
    for feature in gb_record.features[:5]:  # Affiche les 5 premières
        print(f"- Type: {feature.type}, Location: {feature.location}")
