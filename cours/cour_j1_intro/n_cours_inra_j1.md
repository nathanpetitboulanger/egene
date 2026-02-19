# Cours de Génomique : Fondamentaux, Stratégies et Applications

## I. Introduction à la Génomique

### 1. Distinction entre Génétique et Génomique

Bien que liées, ces deux disciplines diffèrent par leur échelle et leur approche :

* **Génétique** : Étude de l'hérédité et de la fonction des gènes individuels. Elle se concentre sur la manière dont les traits sont transmis de génération en génération.
* **Génomique** : Étude globale de l'ensemble du matériel génétique (le génome) d'un organisme. Elle analyse l'organisation, la structure et les interactions entre tous les gènes, ainsi que les régions non codantes et l'influence de l'environnement sur le génome.

### 2. Le Domaine des Eucaryotes

Chez les eucaryotes, la génomique est complexifiée par :

* La présence d'**introns** (séquences non codantes au sein des gènes) et d'**exons**.
* Une grande proportion de **séquences répétées** et d'éléments transposables.
* Des génomes souvent très vastes par rapport aux procaryotes.
* La compartimentation cellulaire (génomes nucléaire, mitochondrial et chloroplastique).

### 3. Histoire et Évolutions du Séquençage

L'histoire du séquençage est marquée par trois ruptures technologiques majeures :

* **1ère Génération (Sanger, 1977)** : Utilisation de didésoxyribonucléotides (ddNTP) pour stopper la synthèse d'ADN. Très précis mais lent et limité à de petits fragments.
* **2ème Génération (NGS - Illumina)** : Séquençage massif en parallèle de fragments courts. A permis une chute drastique des coûts et une augmentation phénoménale du débit.
* **3ème Génération (TGS - Long Reads)** : Technologies comme **PacBio** (HiFi reads) ou **Oxford Nanopore**. Elles permettent de lire de très longs fragments d'ADN (plusieurs dizaines de kb), facilitant l'assemblage de régions répétitives complexes.

---

## II. Stratégies de Séquençage et Assemblage

### 1. Stratégie "BAC-to-BAC" (Hiérarchique)

Cette approche a été utilisée pour le premier génome humain.

* **Principe** : Le génome est découpé en grands fragments (100 à 200 kb) insérés dans des vecteurs **BAC** (Bacterial Artificial Chromosomes). Ces BAC sont cartographiés physiquement avant d'être séquencés individuellement.
* **Avantage** : Très robuste pour les génomes complexes et riches en répétitions, car on sait d'où vient chaque fragment.
* **Inconvénient** : Coûteux, lent et nécessite une infrastructure lourde (automates de repiquage, banques BAC).

### 2. Stratégie "Whole Genome Shotgun" (WGS)

* **Principe** : Tout le génome est cassé aléatoirement en millions de petits fragments qui sont séquencés sans cartographie préalable. L'ordinateur réassemble ensuite le puzzle par bioinformatique en cherchant les chevauchements.
* **Avantage** : Extrêmement rapide et beaucoup moins cher. C'est la norme actuelle, surtout couplée aux *long reads*.

### 3. De l'ADN Brut au Chromosome : Contigs et Scaffolds

L'assemblage se fait par étapes :

* **Contig** : Une séquence continue obtenue par l'alignement de fragments se chevauchant sans aucune lacune.
* **Scaffold** : Un assemblage de contigs orientés et ordonnés. Les contigs sont séparés par des "gaps" (lacunes de taille connue, représentées par des 'N') dont la continuité a été établie par d'autres preuves (ex: paires de lectures, cartes optiques).

### 4. Métriques de Qualité des Génomes

Pour juger de la qualité d'un assemblage, on utilise deux types de mesures :

* **N50** : C'est une métrique de continuité. Si le N50 est de 1 Mb, cela signifie que 50% de l'assemblage total est contenu dans des pièces (contigs ou scaffolds) d'au moins 1 Mb. Plus le N50 est élevé, plus l'assemblage est continu.
  * **BUSCO** (*Benchmarking Universal Single-Copy Orthologs*) : C'est une métrique de complétude. On cherche dans l'assemblage des gènes "universels" censés être présents en une seule copie chez l'espèce étudiée. Un score BUSCO de 98% indique un génome très complet.

### 5. Annotation des Séquences

Une fois le génome assemblé, il faut l'annoter pour lui donner un sens biologique :

* **Annotation Structurale** : Identification des gènes, des exons, des introns et des sites d'épissage. On utilise souvent des données d'**ARN (Transcriptomique)** pour valider les régions exprimées.
* **Annotation Fonctionnelle** : Attribution d'une fonction aux gènes identifiés par comparaison avec des bases de données de protéines connues (Homologie).

### 6. Bulk vs Single-cell Sequencing

Les approches modernes permettent d'analyser l'expression à différentes échelles :

* **Séquençage Bulk** : On analyse un mélange de milliers de cellules. Le résultat est une **moyenne** de l'expression génique du tissu.
* **Séquençage Single-cell (scRNA-seq)** : On analyse chaque cellule individuellement. Cela permet de découvrir l'**hétérogénéité cellulaire**, d'identifier des types cellulaires rares ou de suivre des trajectoires de développement.

---

## III. Focus : La Complexité du Vivant

### 1. Le Génome Humain

* **Taille** : ~3,2 Gigabases (milliards de paires de bases).
* **Gènes** : ~20 000 gènes codant pour des protéines (bien moins que prévu initialement !).
* **Projet HGP** : Lancé en 1990, achevé en 2003. La version "Telomere-to-Telomere" (T2T) complétant les zones ultra-répétitives n'a été publiée qu'en 2022.

### 2. Le Génome du Blé (*Triticum aestivum*)

Le blé est l'un des plus grands défis de la génomique végétale :

* **Ploïdie** : Hexaploïde (AABBDD), issu de deux hybridations successives.
* **Taille** : ~16 à 17 Gb (soit 5 à 6 fois le génome humain).
* **Complexité** : Plus de 85% de séquences répétitives, rendant l'assemblage par fragments courts presque impossible.

### 3. Le Cas de la Vanille (*Vanilla planifolia*)

Le séquençage du vanillier a révélé des défis uniques :

* **Taille** : ~4 Gb.
* **Difficultés** : Forte hétérozygotie et phénomène d'**endoréplication partielle** (certaines cellules ne répliquent qu'une partie de leur ADN), compliquant l'estimation de la taille et l'assemblage.
* **Enjeu** : Faible diversité génétique mondiale, rendant la culture vulnérable aux maladies comme la fusariose.

### 4. La Pangénomique

Aujourd'hui, on ne se contente plus d'un seul "génome de référence".

* **Pangénome** : L'ensemble des gènes trouvés au sein d'une espèce.
* **Core Genome** : Gènes présents chez tous les individus (fonctions vitales).
* **Accessory/Variable Genome** : Gènes présents seulement chez certains individus (adaptation, résistance aux maladies).

---

## IV. Objectifs de la Restitution (Projet de Groupe)

À l'issue de ce module, chaque groupe devra présenter une analyse de projet de séquençage :

* **Format** : 15 minutes de présentation + questions.
* **Contenu attendu** :
  * **Question scientifique** : Pourquoi séquencer cet organisme ?
  * **Caractéristiques biologiques** : Ploïdie, taille estimée du génome, hétérozygotie.
  * **Stratégie technique** : Type de données (Long/Short reads), profondeur de couverture.
  * **Résultats d'assemblage** : Métriques de qualité (N50, BUSCO).
