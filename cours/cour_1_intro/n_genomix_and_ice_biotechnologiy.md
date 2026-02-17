# Cours : Génomique et Biotechnologies Végétales

## Focus : La Génomique du Riz (*Oryza sativa*)

Le riz est l'organisme modèle par excellence pour les céréales en raison de ses caractéristiques biologiques et génomiques uniques.

---

### 1. Pourquoi le Riz comme modèle ?

Le riz possède plusieurs atouts majeurs qui facilitent son étude génomique :

* **Taille du génome** : C'est l'un des plus petits génomes parmi les graminées (environ 430 Mb), ce qui le rend "simple" à manipuler et à séquencer par rapport au blé ou au maïs.
* **Génome de référence** : L'intégralité du génome est séquencée et annotée.
* **QTL Orthologues** : L'identification de zones d'intérêt (Quantitative Trait Loci) chez le riz permet souvent de retrouver des fonctions similaires chez d'autres espèces.
* **Marqueurs moléculaires** : L'abondance de marqueurs permet un positionnement précis sur la carte génétique.
* **Annotation par transcriptomique** : L'utilisation de bases de données d'ARN (RNA-seq) est cruciale pour identifier les régions codantes (exons) et valider l'expression des gènes.

### 2. Colinéarité et Synténie

L'un des concepts les plus puissants en génomique végétale est la **synténie** (conservation de l'ordre des gènes entre espèces).

* Grâce à la colinéarité, les informations découvertes sur le génome du riz peuvent être transférées à d'autres espèces plus complexes.
* **Comparaison Riz vs Arabidopsis** : Bien que très éloignées, ces deux plantes partagent de nombreux gènes communs impliqués dans les fonctions métaboliques de base.

### 3. Stratégies de Séquençage

On distingue deux approches majeures dans l'histoire de la génomique :

1. **Génome de Référence** : Génération d'une séquence "étalon" très précise.
2. **Séquençage Individuel** : On utilise la "carte" (le génome de référence) comme guide pour reconstruire le génome d'un individu spécifique par ré-alignement.

*Note historique* : Le Projet Génome Humain utilisait initialement une approche par "BAC par BAC" (hiérarchique), alors que les techniques modernes (Next Generation Sequencing) reposent sur le séquençage massif de fragments courts avec un fort recouvrement pour obtenir un consensus global.

---

### 4. Construction d'une banque BAC (Bacterial Artificial Chromosome)

La construction d'une banque BAC est une étape clé pour la cartographie physique. Voici le processus détaillé :

#### A. Préparation de l'ADN

* **Extraction** : L'ADN de haut poids moléculaire est extrait.
* **PFGE (Pulsed-Field Gel Electrophoresis)** : Une électrophorèse en champ pulsé est nécessaire pour manipuler de très grands fragments d'ADN. Contrairement à une électrophorèse classique, le champ électrique change de direction, permettant de "débobiner" les fragments et d'atteindre une résolution de l'ordre de 25 kb à plusieurs mégabases.

#### B. Clonage et Transformation

* **Insertion** : Les fragments d'ADN (inserts) sont insérés dans des vecteurs BAC.
* **Transformation** : On utilise des bactéries *E. coli* pour multiplier ces vecteurs.
* **Criblage Bleu/Blanc** : Les bactéries sont étalées sur des boîtes de Pétri carrées contenant du X-Gal.
  * **Points Bleus** : Colonies sans insert (le gène lacZ est fonctionnel).
  * **Points Blancs** : Colonies contenant un insert d'ADN (le gène lacZ est inactivé).

#### C. Automatisation et Stockage

* Un **automate de repiquage** (mechanical head) identifie et prélève les colonies blanches.
* Les clones sont transférés dans des plaques multi-puits (96 ou 384 puits) pour constituer la banque.

---

### 5. Statistiques de Couverture : Formule de Clarke et Carbon

Pour s'assurer que la banque BAC couvre l'intégralité du génome avec une probabilité donnée, on utilise la formule de **Clarke et Carbon** :

$$N = \frac{\ln(1 - P)}{\ln(1 - f)}$$

Où :

* **N** : Nombre de clones nécessaires.
* **P** : Probabilité souhaitée de trouver n'importe quelle séquence du génome (ex: 0.99).
* **f** : Fraction du génome représentée par un seul clone ($f = \frac{\text{taille insert}}{\text{taille génome}}$).

### 6. Le 3D Pooling

Pour retrouver un clone spécifique sans séquencer chaque puits individuellement, on utilise le **3D pooling**. On mélange l'ADN des clones selon trois axes (X, Y, Z : lignes, colonnes, plaques). Une simple PCR sur ces mélanges permet d'identifier par intersection les coordonnées exactes du clone d'intérêt.

---

### 7. Séquençage Classique (Méthode de Sanger)

La méthode de Sanger (ou méthode des didésoxyribonucléotides) repose sur l'arrêt de la synthèse d'ADN :

1. On utilise des **ddNTPs** (didésoxyribonucléotides) marqués par fluorescence.
2. Ces ddNTPs ne possèdent pas de groupement -OH en 3', empêchant l'ajout d'un nouveau nucléotide.
3. La synthèse s'arrête aléatoirement à chaque position, créant des fragments de toutes les tailles possibles.
4. L'électrophorèse capillaire sépare ces fragments par taille, et le laser lit la fluorescence pour déduire la séquence.
