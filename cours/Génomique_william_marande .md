# La génomique

## Introduction - Qu'est-ce que la génomique

Système intégratif : c'est l'étude de l'ensemble des gènes d'un organisme.\
C'est donc la génétique (suivi de l'hérédité) pour plusieurs gèènes différents.

Il existe une grande diversité génétique dans le monde, et pour la majeure partie, cela vient d'EV microscopiques (archées, bactéries, protistes,...). Ce que l'on étudie en général ne constitue qu'une partie de la diversité.

La taille du génome n'est pas en lien avec la complexité ou la taille de l'organisme. Le génome le plus grand actuellement appartient à une petite amibe (protistes)

Parfois, au sein d'une même espèce, entre plusieurs variétés, on peut accumuler un très très grande proportion. C'est assez spécifiques aux plantes, car chez les animaux, entre deux individus, on n'observe que très peu de variabilité.

## Séquençage

Début du projet humain, avec deux stratégies qui s'opposent. \

- approche bac à bac (vue ce matin)
- approche shotgun : on prend l'adn, on fragmente tout. Pas d'organisation structurelle au préalable, demande beaucoup de puissance de calcul par rapport à une approche bac à bac qui demande une bien meilleure organisation.

L'approche bac à bac se faisant globalement jusqu'en 2015, mais ne se fait plus aujourd'hui. \
Le génome humain compte 3 milliards de pdb, pour plus ou moins 20 000 gènes. C'est à dire qu'on a seulement 1,5% du génome qui est codant.

Impact du projet génome humain :

- premier exemple de big science international en Biologie (contrairement à la physique qui montrait déjà de la collaboration mondiale).
- Se rendre compte de la quantité d'ADN poubelle

Aujourd'hui, le gap est énorme, et on ne parle pas de refaire le séquençage pour un individu, mais pour 100 000.

### Séquençage du génome du blé

C'est un énorme génome, bien plus grand que celui de l'humain. Il est constitué de plein de séquences répétées, il est donc très difficile de le reconstituer avec certitude. Ajoute à cela le fait qu'il est hexaploïde, ce qui rajoute encore une couche de complexité.

Pour le déchiffrer, on a du combiner toutes les technologies à disposition. \
Après travail en collaboration avec des informaticiens, ils ont créer un protocole permettant d'aboutir à ces informations.

### Séquençage shotgun

Différence entre contigs et scafold :

- contigs : assemblage brut des lectures du séquenceur. Ce sont les données brutes de séquençage
- scafolds : Données traitées ??

#### Révolution technologique

Maintenant, on fait des approches de pyroséquençage puis illumina, on divise les réactions dans des petites goutelettes, permettant de diviser les réactions et de réduire les coûts.

Next generation squencing (NGS) : on lit dix fois les mêmes molécules dan un même puit, ce qui a vraiment changer les choses. Aujourd'hui, on sait séquences un chromose entier sans trou, en fusionnant plusieurs technologies.

Métriques importantes :

- N50 : En lien avec la taille des contigs, plus il est grand, plus on de la qualité
- BUSCO : C'est la qualité qu'on met sur la fiabilité de la séquence. Dans un individu qui appartient au groupe des plantes, on doit trouver un ensemble de 500 gènes, qui sont toujours présents peu importe l'individu. On regarde alors ces gènes sur notre nouveau séquençage pour attribuer une note de fiabilité au séquençage qu'on vient de faire.

La génomique des plantes avance lentement, car la plupart des espèces cultivées ont des très grands génomes (pour certains même plus grand que celui de l'humain).

#### Problème de la Vanille

Modèle complexe car très peu de diversité, se développe par bouture comme la banane.
Existence d'un phénomène d'endoréplication ...

PanGénomique : domaine le plus porteur en bioinformatique. C'est regarder toutes la variabilité et la diversité génétique qu'on peut retrouver à l'intérieur d'une même espèce. Aujourd'hui on travaille avec des pangénomes.

Il est toujours très important de se poser la question de : pourquoi ? je fais ce que je fais.

## Modification des génomes

Une fois qu'on a la séquence génomique, on ouvre la porte de nombreuses choses.

Par exemple, depuis l'arrivée de kripr cas9, on s'intéresse beaucoup plus à la fonction du gène et à leurs modifications. \
C'est l'exemple des vaches qu'on a modifié pour qu'elles n'aient plus de cornes.

Autre exemple d'une enfant ayant une leucémie, à laquelle on a injecter des cellules préalablement modifiée par crispr-kas9.

Ces utilisations posent par contre de nombreuses questions éthiques.

# Restitution

- Descirption : question scientifique imaginée ?
- Comment je produis le génome du projet (parler des caractéristiques des l'espèce, parler de la technologie HiC)
- Parler du séquençage et du type de séquenceur choisi.
