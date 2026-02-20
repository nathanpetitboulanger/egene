from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

Entrez.email = "tauzin_pierre@orange.fr"


def telecharger_sequence(id_ncbi):
    """Petite fonction pour éviter de copier-coller le code Entrez"""
    with Entrez.efetch(
        db="nucleotide", id=id_ncbi, rettype="fasta", retmode="text"
    ) as handle:
        return SeqIO.read(handle, "fasta")


seq_human = telecharger_sequence("NM_000207")
seq_pig = telecharger_sequence("NM_001109772.2")

print(f"Séquence humaine : {len(seq_human)}")
print(f"Séquence du porc : {len(seq_pig)}")

aligner = PairwiseAligner()

aligner.mode = "global"

aligner.match_score = 1  # Identique (+1 point)
aligner.mismatch_score = -1  # Différent (-1 point)
aligner.open_gap_score = -2  # Ouvrir un trou (-2 points)
aligner.extend_gap_score = -1  # Agrandir un trou (-1 point)

alignement = aligner.align(seq_human, seq_pig)
alignement

best_align = alignement[0]
print(best_align.score)

# A noter : ici, on obtient un score négatif, donc les séquences sont fortement différentes. Ce n'est pas le résultat attendu car ces deux séquences sont assez proches l'une de l'autre.
# Ces différences viennent de la différence entre les deux fstas, l'un ne prend que les séquences codantes et l'autres prends les utr aussi

aligner_local = PairwiseAligner()

aligner_local.mode = "local"

aligner_local.match_score = 1
aligner_local.mismatch_score = -1
aligner_local.open_gap_score = -2
aligner_local.extend_gap_score = -1

align_local = aligner_local.align(seq_human.seq, seq_pig.seq)
align_local

best_align_local = align_local[0]
print(best_align_local)
print(best_align_local.score)

# Comme le score reste très bas, je vais essayer de traduire pour comparer la séquence protéique.
# La différence pourrait s'expliquer par la différence du code génétique ? me parait difficile d'expliquer autant de différences mais c'est à tester


prot_human = seq_human.seq.translate(to_stop=True)
prot_pig = seq_pig.seq.translate(to_stop=True)

len(prot_human)
len(prot_pig)

# Il y a encore un problème important de longueur entre les deux, il est donc probable que j'ai pris un gène de humain mRNA (donc post épissage) et une CDS de porc (donc plus longue), l'idée est de transformer cette séquence CDS en mRNA
# Avant de passer à une telle étape, on va d'abord juste voir si ce n'est pas un problème d'ATG (si ça se trouve la protéine du porc commence par le 5'UTR et je l'ai traduite depuis cette zone, expliquant sa longueur)


start_human = seq_human.seq.find("ATG")
print(start_human)

human_marn = seq_human.seq[start_human:]
prot_human_processed = human_marn.translate(to_stop=True)
len(prot_human_processed)

start_pig = seq_pig.seq.find("ATG")
print(start_pig)

pig_marn = seq_pig.seq[start_pig:]
prot_pig_processed = pig_marn.translate(to_stop=True)
len(prot_pig_processed)


aligner.mode = "global"

aligner.match_score = 1  # Identique (+1 point)
aligner.mismatch_score = -1  # Différent (-1 point)
aligner.open_gap_score = -2  # Ouvrir un trou (-2 points)
aligner.extend_gap_score = -1  # Agrandir un trou (-1 point)

prot_align = aligner.align(prot_pig_processed, prot_human_processed)
prot_align

best_prot_align = prot_align[0]
print(best_prot_align.score)

print(best_prot_align)


# Au final, le problème venait de la mauvaise acession pour l'insuline du porc au début.
# Note à moi même : apprendre à utiliser le ncbi


# Pour les comparaisons de protéines, on n'utilise pas des matrices +1 -1, mais des matrices spécifiques qui donnent un score spécifique à chaque remplacement, car un remplacement de leucine/isoleucine entraine beaucoup moins de problèmes qu'un remplacement leucine / glycine
from Bio.Align import substitution_matrices

matrice = substitution_matrices.load("BLOSUM62")

aligner = PairwiseAligner()
aligner.substitution_matrix = matrice

aligner.open_gap_score = -10
aligner.extend_gap_score = -0, 5
aligner.mode = "global"

prot_align_matrix = aligner.align(prot_human_processed, prot_pig_processed)
best_prot_align_matrix = prot_align_matrix[0]
print(best_prot_align_matrix.score)
