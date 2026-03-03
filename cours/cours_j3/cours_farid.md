# Cours : Immunoprécipitation de la Chromatine (ChIP) et Régulation Transcriptionnelle

Ce cours traite des mécanismes d'interaction entre les protéines et l'ADN, et des méthodes permettant de cartographier ces interactions à l'échelle du génome.

## 1. Mécanismes de la Régulation Transcriptionnelle

La transcription des gènes est un processus complexe orchestré par l'interaction physique entre l'ADN et diverses protéines.

### 1.1. Les Acteurs de la Transcription

* **Promoteur** : Région d'ADN située à proximité du site d'initiation de la transcription. C'est la plateforme d'atterrissage pour la machinerie basale de transcription.
* **Éléments de transcription canoniques** : Séquences d'ADN hautement conservées, comme la **boîte TATA**, qui permettent le recrutement des facteurs de transcription généraux.
* **Facteurs de Transcription (FT)** : Protéines régulatrices qui se lient à l'ADN de manière séquence-spécifique pour activer ou réprimer l'expression d'un gène.

### 1.2. Architecture 3D et Boucles de Chromatine

L'ADN ne doit pas être vu comme une ligne droite. La régulation implique souvent des **boucles de chromatine (looping)**. Des éléments régulateurs distants (appelés *enhancers* ou amplificateurs) peuvent être rapprochés physiquement du promoteur grâce à des complexes protéiques (comme le complexe Mediator ou les cohésines), formant une boucle. Cela permet à des signaux cellulaires distants de contrôler précisément l'expression d'un gène.

---

## 2. L'Immunoprécipitation de la Chromatine (ChIP)

La technique de **ChIP** (*Chromatin ImmunoPrecipitation*) permet d'identifier les séquences d'ADN spécifiques auxquelles une protéine donnée est liée *in vivo*.

### 2.1. Protocole Expérimental

1. **Réticulation (Cross-linking)** : On utilise du formaldéhyde pour fixer de manière covalente les protéines sur l'ADN au sein des cellules vivantes.
2. **Fragmentation (Sonication)** : L'ADN génomique est extrait puis cassé en petits fragments (200 à 500 paires de bases) par ultrasons.
3. **Immunoprécipitation** : On utilise un anticorps spécifique dirigé contre la protéine d'intérêt (ex: un facteur de transcription ou une histone modifiée). L'anticorps "pêche" le complexe protéine-ADN.
4. **Défixation (Reverse cross-linking)** : Les liaisons protéines-ADN sont rompues, les protéines sont éliminées, et on récupère l'ADN purifié qui était lié à la cible.

---

## 3. ChIP-seq : De l'Expérience à la Bioinformatique

Le **ChIP-seq** combine la ChIP avec le séquençage à haut débit (NGS) pour obtenir une cartographie globale.

### 3.1. Workflow Analytique

* **Séquençage** : Les fragments d'ADN issus de la ChIP sont séquencés (généralement sur technologie Illumina).
* **Alignement (Mapping)** : Les lectures (*reads*) sont alignées sur le génome de référence de l'espèce étudiée.
* **Peak Calling** : On identifie les zones du génome où la densité de reads est anormalement élevée par rapport à un témoin (input). Ces "pics" correspondent aux sites de fixation de la protéine.

### 3.2. Applications

* Identifier les gènes cibles d'un nouveau facteur de transcription.
* Étudier le "code des histones" (marques épigénétiques).
* Comprendre les dérèglements transcriptionnels dans des pathologies comme le cancer.

---

## 4. Récapitulatif des Techniques de ChIP

| Technique / Étape | Fonction |
| :--- | :--- |
| **ChIP** | Étude des interactions physiques entre protéines et ADN. |
| **ChIP-seq** | Cartographie génomique globale des sites de fixation via NGS. |
| **Cross-linking** | Fixation chimique (formaldéhyde) des complexes protéine-ADN *in vivo*. |
| **Sonication** | Fragmentation physique de l'ADN par ultrasons pour obtenir de petits fragments. |
| **Peak Calling** | Analyse bioinformatique pour repérer les enrichissements significatifs de reads. |
| **Looping** | Formation de boucles pour rapprocher des régulateurs distants (Enhancers). |
