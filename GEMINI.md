# Egene - Génétique et Bioinformatique

Ce dépôt est organisé par projets, chacun traitant d'un aspect spécifique de la génétique ou de la bioinformatique.

## Projets

### 1. Barcoding ADN (`barcoding-adn/`)
Le premier projet du dépôt se concentre sur l'étude du barcoding moléculaire.
- `01_fetch_data.py` : Récupération de données depuis BOLD Systems API.
- `02_data_exploration.py` : Exploration et nettoyage des données.
- `03_phylogeny.py` : Analyse phylogénétique.
- `04_barcode_gap.py` : Analyse du "barcode gap".

### 2. Exercices de Génétique Classique
- `exercice_genetique_insuline.py` : Exercice complet et pédagogique sur la comparaison de l'insuline (Humain vs Souris). Inclut la récupération NCBI, l'alignement et l'analyse.

## Scripts Utilitaires

- `pull_genome.py` : Script de base pour récupérer des séquences de protéines depuis le NCBI (Entrez).

## Rôles de l'IA

- **Correcteur / Rédacteur** : Transformation de notes brutes en supports de cours structurés et correction de contenu pédagogique.

## Technologies et Bibliothèques

- **Python** : Langage principal.
- **Biopython** : Manipulation de séquences biologiques et accès aux bases de données (Entrez).
- **Requests** : Interaction avec les APIs (BOLD Systems).
- **Pandas/Matplotlib/NumPy** : Analyse et visualisation de données.

## Utilisation

### Prérequis

Le projet utilise `uv` pour la gestion des dépendances.

```bash
uv sync
```

### Exécution des scripts

- Pour l'exercice d'insuline : `python exercice_genetique_insuline.py`
- Pour le projet barcoding :
  ```bash
  cd barcoding-adn
  python 01_fetch_data.py
  # Puis suivre l'ordre des scripts
  ```

## Sources de données

- **NCBI (Entrez)** : Séquences de protéines et génomes.
- **BOLD Systems** : Données de barcoding (Barcode of Life Data Systems).
