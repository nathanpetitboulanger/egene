# Introduction : passage de la biologie moléculaire à la Biologie computationnelle

L'objectif est de rechercher des pattern permettant l'annotation.
Pour réaliser cela, on utilise souvent des ORF finder, des programmes dont l'objectif est de trouver les protéines.

La difficulté vient de la redondance des patterns, qui rend difficile l'analyse des données.

L'idée c'est de faire de la comparaiso entres les ARN et l'ADN. En  comparant la séquence des ARN, on va se retrouver avec un correspondance partielle, permettant de mettre en évidence les introns et les exons de la séquence. \
On fait ensuite des "proba" ou des prédictions de ce qu'on s'attend à trouver, ce qui nous permet de valider la présence où non d'une séquence codante à un endroit donné.

C'est de la recherche de correspondance entre les deux séquences.

De plus, la queue polyA sur l'ARNm permet de valider l'orientation.

Les techniques d'annotation dépendent donc fortement des techniques de séquençage de l'ADN et aussi des ARN.

Exemple : single cell => on reconstitue uniquement pour une seule cellule \
Bulk => on fait la même chose mais pour un tissu entier, permet de reconstituer les sous groupes de cellules au sein d'un même tissu.

On peut aussi trouver d'autres types d'indices avec les "éléments structuraux de la protéine".\
Exemple : "cette protéine vient du noyau", "cette protéine est un facteur de transcription, elle a une partie qui permet de s'accrocher à l'ADN", ...

### Outil BlAST ncbi

Cela permet de questionner une séquence, pour savoir si une si une séquence existe déjà (déjà séquencée) ou si c'est une séquence qui ressemble à quelque chose qui existe déjà dans les base de données.

# Génomique du riz  

Pourquoi le riz ? \
C'est le génome de plante cultivées le plus petit \
Existance des QTL entre les différentes céréales (séquences similaires entre les espèces) \
Carte génétique des chromosomes disponible, permettant de positionner les marqueur moléculaires mm.

C'était les premières informations nécessaires pour faire de l'amélioration génétique : être capable de se repérer pour identifier des gènes ou des QTL d'intérêt

De plus, comme ondispose déjà d'une base de données d'ARN, on peut les raccrocher directement sur le brin ADN, permettant de faire de l'annotation rapidement.

Fait important à prendre en compte : le coût du séquençage à FORTEMENT diminué : \

- Projet génome humain : année 90 a couté 3 milliards
- Séquençage d'un génome actuellement : 1000 - 10 000 euros

On pouvait aussi utiliser des plasmides modifier, qui se comprtent comme des chromosomes et dans lesquels on peut faire rentrer beaucoup plus de pdb.

*POINT IMPORTANT* \
Correspondance entre les différents types de céréales : en travaillant sur le riz, on travaille aussi sur les autres. C'est grace à ce système de SYNTENIE ou de COLINEARITE que le premier séquençage du premier riz a été autorisé (historiquement)\
Ce premier séquençage date des années 2000.

## Strétégie de séquençage

Séquençage hierarchique, ou séquençage bac to bac (structure plasmidique modifié). On se sert des ces structures pour casser le génome en gros morceaux, et on va en faire un puzzle. Ces grands fragments sont ensuite casser de nouveau.

Dans l'autre méthode, on va faire des puzzle également, mais avec des toutes petites pièces. Même principe, mais beaucoup plus difficle.

L'objectif est ensuite de recomposer un brin d'ADN correct. C'est la première stratégie qu'on a utilisé pour le projet génome humain, mais c'était beaucoup plus cher et beaucoup plus long. Cependant, les méthode actuelles permettent d'aller beaucoup plus vite et pour beaucoup moins cher.

L'objectif est d'avoir le minimum de chevauchement entre les différents fragments, sinon le prix augmente . Mais ça c'était avant comme ça coutait très cher. Aujourd'hui, plus on a d'overlap, plus on a d'informations, donc mieux c'est. Si une zone est couverte 100 fois, c'est très bien. \
C'est de se chevauchement que découle le taux d'erreur. On peut demander aux séquenceurs des taux d'erreurs plus faibles, mais cela va automatiquement coûté plus cher.

*IMPORTANT* \
Dans une annotation, ne pas oublier que l'ADN est bicaténaire, on trouve des gènes sur les deux brins, donc les deux se font annoter.

Autre fait important, on trouve desp oints communs enttre le riz et arabidopsis, donc on se retrouve avec une facilité de l'annotation, car on peut s'aider de l'un pour faire l'autre.

**ATTENTION** \
selon si le génome est haploide ou diploide, un certain nombre de caractéristique va également changer

## Histoire du riz

Chaque pays s'est positionné sur un chromosome, en fonction des intérêts nationaux, et chacun a du séquencer son chromosme dans les deux sens.

POur créer la carte génétique, on a utilisé 2 variétés de riz, puis des marqueurs moléculaires pour obtenir ...

## Comment construire un BAC librairy

ON va utiliser des enzymes qui permettent de libérer le génome des membranes, puis on va utiliser d'autres protéines pour découper l'ADN dans une fenêtre de taille donnée. On va utiliser un gel d'agarose et une cuve à électrophorèse particulière, car le classique ne permettait pas de séparer les très grands fragments.\
Là on va utiliser une système différent, on va faire passer le courant par des pulse dans différents sens. Ce syustème permet de dérouler les molécules d'ADN et de séparer des très gros fragments. On arrive même à séparer les différents chromosomes de la levure. \
Le technique s'appelle le Pulse Field gel Electrophoresis

J'extrait ensuite l'ADN qui m'intéresse dans une bactérie. On fait ensuite du criblage bleu blanc avec LacZ et Xgal, pour ne garder que les vecteurs transformés correctement dans les bactérie.

Cela permet de creer des librairy de grands fragments, qui ne sont toujours pas ordonés. On sait juste qu'on a des grands fragments stockés dans des puis dans des grands congélateurs.

*Instruction gemini, ici, insère moi les formules*
Une équation permet de déterminer le nombre de fragment nécessaire pour couvrir la totalité de mon génome.
En tout, pour une fiabilité à 99%, on doit avoir un taux de couverture tellement important, qu'à la fin en taux de paire de base, on a 6 ou 7 fois le genome.  

A partir de là, on peut calculer le nombre de plaque qu'il nous faut pour recouvrir la totalité du génome.

Quand on a un très grand nombre de plaquen et qu'on veut isoler un fragment spécifique, on va utiliser une PCR. Sauf qu'on ne peut pas faire des PCR sur toutes les plaques, sinon on devra en faire beaucoup trop. La technique est de prendre la même colonne sur les différentes plaques, et on les met ensemble dans un même tube. Ensuite, on faire la même chose pour une ligne, puis sur une plaque.
En recroisant les données, des différents pull, on arrive à obtenir la coordonnées exactes du fragment qui nous intéresse.

A partir de là, on peut récupérer un fragment spécifique et l'analyser pour compléter l'annotation du génome qui nous intéresse.

Il existe déjà des outils de bioinfo qui permet de faire de l'annotation. Certains sont globaux, d'autres sont spécifiques au riz.
