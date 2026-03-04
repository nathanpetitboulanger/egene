# TD Corrigé – Barcoding pollinisateurs – M. OLLIVIER

## Analyse d’un jeu de données sur les Bourdons (Bombus sp.)

Ce document contient les consignes du TD ainsi que les explications théoriques
et les réponses aux questions basées sur les données du portail BOLD.

---

## 1. Exploration de la plateforme et des données publiques

**Plateforme BOLD Systems :**
[https://boldsystems.org/](https://boldsystems.org/)

**a. Combien de séquences de référence (barcodes) sont disponibles publiquement
?**
> **RÉPONSE :** La plateforme contient actuellement environ **10 millions de
> séquences** de référence disponibles publiquement.

**b. Combien d’espèce cela représente-t-il ? A quels règnes appartiennent ces
espèces ?**
> **RÉPONSE :** Cela représente environ **900 000 à 1 000 000 d'espèces**.
> Elles appartiennent majoritairement au règne **Animal** (plus de 80%,
> notamment les insectes), mais on y trouve aussi des **Plantes**, des
> **Champignons** et des **Protistes**.

**c. Qu’est-ce qu’un BIN ? Quelles opportunité cette approche offre-t-elle ?**
> **DÉFINITION :** Un **BIN** (Barcode Index Number) est un identifiant unique
> (ex: BOLD:AAA1234) attribué automatiquement par un algorithme à un groupe de
> séquences d'ADN très proches. **OPPORTUNITÉS :**
>
> 1. Permet d'identifier des spécimens même si l'espèce n'a pas encore de nom
>    scientifique (OTU - Operational Taxonomic Unit).
> 2. Aide à la découverte d'**espèces cryptiques** (deux espèces identiques
>    physiquement mais génétiquement distinctes).
> 3. Standardise l'identification à l'échelle mondiale sans dépendre uniquement
>    de l'expertise morphologique.

**d. Identification d’une séquence inconnue.**

* **Base de données choisie :** "Species Level Barcode Records" (recommandée
pour la précision).
* **Séquence testée :** (Séquence de 658 bp du gène COI-5P fournie dans le TD).

> **RÉSULTAT DE L'IDENTIFICATION :**
>
> * **Espèce :** *Bombus terrestris* (le Bourdon terrestre).
> * **Similarité :** **100%**.
> * **Preuves :** Le score de correspondance est parfait. La localisation
>   (Allemagne) et la période (Juillet) correspondent parfaitement à l'écologie
>   de cette espèce très répandue.

---

## 2. Analyses des données du projet (CODABEILLES)

**a. Description du jeu de données « DS-BBBABSV » :**
> **RÉPONSE :** Ce dataset est composé de **111 spécimens** de bourdons
> français.
>
> * **Marqueur génétique :** COI-5P (Cytochrome Oxidase Subunit 1).
> * **Longueur moyenne :** Environ **658 paires de bases (bp)**.
> * **Qualité :** Haute, la plupart des séquences sont complètes et sans
>   nucléotides ambigus (N).

**b. Champs du tableau des records :**

* **Process ID :** Code unique de l'analyse en laboratoire.
* **Sample ID :** Identifiant du spécimen dans la collection physique.
* **BIN URI :** Le lien vers le cluster génétique mondial.
* **Pourquoi certaines entrées n'ont pas de séquence ?**

> **RÉPONSE :** Échecs techniques lors de l'extraction de l'ADN (échantillon
> trop vieux ou mal conservé) ou lors de l'amplification PCR.

**c. Construction de l'arbre de distance (Taxon ID Tree) :**

* **Kimura 2 P (K2P) :** Modèle mathématique qui calcule la distance génétique
en ajustant les probabilités de mutations (transitions vs transversions).
* **MUSCLE :** Algorithme qui aligne les séquences pour s'assurer que l'on
compare bien les mêmes positions de nucléotides.

> **ANALYSE DE LA TOPOLOGIE :**
>
> * Un spécimen est mal identifié s'il se retrouve branché loin de ses
>   congénères de la même espèce.
> * **Identification de Bombus sp. :** Si un individu noté "Bombus sp."
>   (morphologie incomplète) tombe dans le même cluster ultra-serré qu'un
>   groupe de *Bombus pascuorum*, on peut conclure qu'il s'agit d'un *Bombus
>   pascuorum*.

**d. Analyse du Barcode Gap :**

* **Principe :** La variation génétique entre deux espèces différentes
(**distance interspécifique**) doit être nettement plus grande que la variation
au sein d'une même espèce (**distance intraspécifique**).

> **INTERPRÉTATION DES GRAPHIQUES :**
>
> * **Scatter Plot :** Si les points sont au-dessus de la droite rouge, le
>   Barcode Gap est validé.
> * **Distance Distribution :** On cherche une séparation nette (un "gap")
>   entre le pic des distances intra (proche de 0%) et le pic des distances
>   inter (souvent > 2-3%).
> * **Problèmes :** Si les distances se chevauchent, cela peut indiquer une
>   erreur d'identification ou une espèce en cours de spéciation.

---

## 3. Analyses intégrant les données mondiales

**Pourquoi intégrer les données mondiales ?**
> **RÉPONSE :** Pour vérifier si nos spécimens locaux sont représentatifs de
> l'espèce entière. Cela permet de s'assurer que le "Barcode" choisi fonctionne
> aussi bien pour un bourdon français que pour un bourdon polonais ou anglais,
> garantissant la robustesse de l'identification moléculaire.
