# Cours 1 : Introduction à la Génomique et Bioinformatique

---

**Note Administrative :** Rendez-vous au CNRGV à 14h30 avec William Maran (sonner à l'entrée 1 ou 2)
---

## 1. Introduction Générale

Ce cours pose les bases de la génomique moderne et explore le rôle crucial de l'informatique dans l'analyse biologique.

* **Rappels des bases** : Structure de l'ADN, transcription, traduction.
* **Annotation de séquences** : Identification des régions fonctionnelles du génome.
* **Apport de l'informatique** : Automatisation de l'analyse, stockage de données massives (Big Data) et modélisation biologique.

## 2. Annotation des Séquences

L'annotation consiste à ajouter des informations biologiques aux séquences d'ADN brutes.

* **Méthodologie** : Utilisation de données complémentaires pour affiner la précision :
  * **Données ARN (Transcriptomique)** : Permettent de valider les régions réellement exprimées.
  * **Données Protéiques** : Permettent de confirmer la fonction biologique par homologie.
* **Objectif principal** : Localiser précisément les **exons** et définir la structure des gènes (introns/exons, sites d'épissage).

## 3. Techniques de Séquençage : Bulk vs Single-cell

Il est crucial de distinguer les deux approches majeures de séquençage actuel :

| Caractéristique | Séquençage Bulk | Séquençage Single-cell (scRNA-seq) |
| :--- | :--- | :--- |
| **Principe** | Analyse d'un mélange de milliers de cellules. | Analyse individuelle de chaque cellule. |
| **Résultat** | Moyenne de l'expression génique de l'échantillon. | Hétérogénéité cellulaire et identification de sous-populations. |
| **Usage** | Comparaison globale (ex: tissu sain vs malade). | Étude du développement, lignages cellulaires, rareté. |
