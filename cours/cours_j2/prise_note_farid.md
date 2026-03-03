# Cours : Accès aux Séquences et Introduction à la Protéomique

Ce module explore les méthodes modernes d'acquisition de données génomiques et les techniques fondamentales de l'analyse des protéines.

## 1. Accès aux Séquences Génomiques

### 1.1. Les Browsers de Génomes (Genome Browsers)

L'accès aux séquences ne se limite pas à la lecture de fichiers bruts. Les "Genome Browsers" (comme UCSC Genome Browser, Ensembl ou NCBI Genome Data Viewer) sont des outils de visualisation interactive permettant d'explorer le génome à différentes échelles. Ils permettent de superposer des annotations (gènes, promoteurs, variants) sur la séquence d'ADN.

### 1.2. Le Séquençage de Nouvelle Génération (NGS)

Le paysage technologique du séquençage a évolué vers des méthodes à haut débit (Next Generation Sequencing) :

* **Illumina (Solexa)** : Basé sur le séquençage par synthèse avec des terminateurs de fluorescence réversibles. C'est la technologie la plus utilisée pour sa précision élevée et son débit massif, bien que les lectures soient courtes (75-300 bp).
* **Oxford Nanopore Technologies (ONT)** : Séquençage de "troisième génération". Il mesure les perturbations du courant électrique lorsqu'une molécule d'ADN traverse un nanopore protéique. Permet des lectures ultra-longues (plusieurs centaines de kb).
* **Pacific Biosciences (PacBio) - SMRT Cell** : *Single Molecule Real Time*. Utilise des puits microscopiques (ZMW - Zero-Mode Waveguides) pour observer l'incorporation de nucléotides en temps réel sur une seule molécule d'ADN. Offre une grande précision (HiFi reads) et des lectures longues.

---

## 2. Introduction à la Protéomique

### 2.1. Définition

La protéomique est l'étude à grande échelle de l'ensemble des protéines d'un organisme, d'un tissu ou d'une cellule (le protéome). Contrairement au génome qui est statique, le protéome est dynamique et varie en fonction de l'environnement, du stade de développement et de l'état pathologique.

### 2.2. Séparation des Protéines : L'Électrophorèse

Avant l'analyse, les protéines doivent souvent être séparées selon leurs propriétés physico-chimiques.

* **Électrophorèse sur gel de polyacrylamide (PAGE)** : Technique utilisant un champ électrique pour faire migrer les protéines à travers une matrice de gel.
* **SDS-PAGE** : En ajoutant du SDS (détergent), les protéines sont dénaturées et chargées négativement de manière uniforme, permettant une séparation basée uniquement sur leur **masse moléculaire**.

### 2.3. Détermination de la Séquence par Clivage

Pour identifier une protéine, on procède souvent à sa digestion enzymatique ou chimique en peptides plus petits :

* **Digestion enzymatique** : Utilisation de la trypsine (coupe après Lysine et Arginine).
* **Clivage chimique** : Utilisation du **bromure de cyanogène (CNBr)**, qui clive spécifiquement les liaisons peptidiques au niveau des résidus méthionine.

---

## 3. Analyse par Spectrométrie de Masse (MS)

La spectrométrie de masse est l'outil central de la protéomique moderne.

### 3.1. MALDI-TOF (Matrix-Assisted Laser Desorption/Ionization - Time of Flight)

Cette technique permet d'ioniser des biomolécules fragiles sans les casser :

1. **Désorption/Ionisation** : L'échantillon est mélangé à une matrice cristalline et frappé par un faisceau laser, ce qui éjecte les molécules sous forme d'ions.
2. **Temps de Vol (TOF)** : Les ions sont accélérés dans un tube sous vide par un champ électrique. Les ions les plus légers voyagent plus vite que les plus lourds. Le temps mis pour atteindre le détecteur permet de calculer précisément la masse (rapport m/z).

### 3.2. Spectrométrie de Masse en Tandem (MS/MS)
Le MS-MS permet d'aller au-delà de la simple mesure de masse pour obtenir la séquence d'un peptide :
1.  **MS1** : Sélection d'un ion peptidique spécifique (ion précurseur) parmi un mélange.
2.  **Cellule de collision** : Le peptide sélectionné est fragmenté par collision avec un gaz inerte.
3.  **MS2** : Analyse de la masse des fragments produits. En mesurant la différence de masse entre les fragments successifs (qui correspond à la masse d'un acide aminé), on peut déduire la séquence primaire du peptide.

---

## 4. Récapitulatif des Techniques

| Technique | Domaine | Caractéristique principale |
| :--- | :--- | :--- |
| **Illumina (Solexa)** | Séquençage (NGS) | Lectures courtes, très haute précision, débit massif. |
| **Oxford Nanopore (ONT)** | Séquençage (3G) | Lectures ultra-longues, temps réel, appareil portable. |
| **PacBio (SMRT)** | Séquençage (3G) | Molécule unique, lectures longues et haute fidélité (HiFi). |
| **SDS-PAGE** | Protéomique | Séparation des protéines selon leur masse moléculaire. |
| **MALDI-TOF** | Spectrométrie de masse | Ionisation douce par laser, mesure du temps de vol. |
| **MS/MS** | Spectrométrie de masse | Fragmentation sélective pour le séquençage peptidique. |
| **CNBr (Bromure de cyanogène)** | Biochimie | Clivage chimique spécifique au niveau des Méthionines. |

