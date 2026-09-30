# Chemical Orbital Visualiser (COV)

[English](README.md) · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · **Français**

Explorez les énergies et les occupations des orbitales, ainsi que les liens entre orbitales atomiques et moléculaires.

COV présente les énergies et les occupations des orbitales dans des diagrammes de niveaux d’énergie interactifs. Vous pouvez exporter les figures et les données. Il lit les fichiers Gaussian FCHK/FCH et Molden.

La préversion v0.4 ajoute l’analyse NBO. Avec la sortie NBO du calcul, vous pouvez voir comment les orbitales atomiques et localisées contribuent aux orbitales moléculaires, suivre leurs liens dans le diagramme et examiner les charges et les informations sur les liaisons.

Sélectionnez une orbitale pour voir sa forme en 3D. Le rendu utilise un GPU NVIDIA. L’interface est disponible en anglais, chinois simplifié, japonais et français.

## Fonctionnalités

- **Lire les diagrammes de niveaux d’énergie.** Consultez les énergies et les occupations électroniques, changez d’unité d’énergie et exportez le diagramme en PNG ou SVG.
- **Suivre les liens entre orbitales — préversion v0.4.** Voyez comment les orbitales atomiques et localisées contribuent aux orbitales moléculaires. Sélectionnez un lien pour examiner sa contribution.
- **Examiner les charges et les liaisons — préversion v0.4.** Consultez les charges NPA, les populations de spin, les indices de liaison de Wiberg et les interactions donneur–accepteur lorsque les fichiers NBO contiennent ces données.
- **Exporter les figures et les données.** Enregistrez les diagrammes de niveaux d’énergie en PNG ou SVG, et les données correspondantes en CSV ou JSON.
- **Parcourir les orbitales.** Accédez à la HOMO ou à la LUMO, recherchez dans la liste des orbitales et affichez les orbitales occupées, virtuelles, de cœur ou de valence.
- **Voir les orbitales en 3D.** Sélectionnez une orbitale, réglez son isosurface, puis faites pivoter la vue moléculaire ou zoomez.

## Télécharger

| Version | Contenu | Téléchargement Windows |
|---|---|---|
| [Version stable v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | Énergies et occupations des orbitales, diagrammes de niveaux d’énergie et vues 3D | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [Préversion v0.4.0-pre.1](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.1) | Ajoute l’analyse NBO, la composition des orbitales et leurs liens | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.1/CUDA-Orbital-Visualisation-v0.4.0-pre.1-Windows-sm120.zip) |

Les téléchargements actuels affichent encore l’ancien nom du produit.

Les paquets Windows actuels sont destinés à la série RTX 50 ; les paquets pour les autres architectures ne sont pas encore disponibles.

## Premiers pas

1. Téléchargez et extrayez une archive ZIP pour Windows. Pour v0.3.0, lancez `cov.exe`. Pour v0.4.0-pre.1, lancez `program/cov.exe`, ou double-cliquez sur un fichier `Open-*.cmd` pour ouvrir un exemple.
2. Ouvrez un fichier Gaussian FCHK/FCH ou un fichier Molden compatible, ou glissez-le dans la fenêtre.
3. Parcourez les énergies et les occupations dans la liste des orbitales et le diagramme. Sélectionnez un niveau pour voir son orbitale, ou exportez le diagramme.
4. Dans la préversion v0.4, ouvrez le dossier du calcul pour charger ensemble les fichiers de fonction d’onde et NBO. Si le dossier contient plusieurs calculs, choisissez celui à ouvrir.

Si `formchk` est installé, COV peut convertir les fichiers Gaussian CHK en FCHK. Les fichiers FCHK/FCH et Molden fournissent les énergies, les occupations et les formes des MO. Les formes des orbitales NBO et l’analyse de leur composition nécessitent aussi le rapport, l’archive `.47` et les matrices orbitalaires correspondants.

## Fichiers d’entrée et configuration requise

- **Fonctions d’onde :** fichiers Gaussian `.fchk` / `.fch`, ou fichiers compatibles `.molden` / `.mol` / `.input`. Les fichiers Gaussian `.chk` peuvent être ouverts avec un utilitaire `formchk` installé.
- **Fichiers NBO — préversion v0.4 :** [Préparer les fichiers de calcul](docs/NBO_ONE_JOB.fr.md) explique comment les produire ; [Utiliser les résultats NBO](docs/AOMO_NBO.fr.md) indique les fichiers nécessaires à chaque vue.
- **Taille des molécules :** jusqu’à 100 atomes par fichier d’entrée.
- **Affichage :** un GPU NVIDIA, un pilote compatible et OpenGL 2.1 ou une version ultérieure. Utilisez une version compilée pour l’architecture de votre GPU.

## Documentation

- [Utiliser COV (en anglais)](docs/UI.md)
- [Utiliser les résultats NBO](docs/AOMO_NBO.fr.md) — préversion v0.4
- [Préparer les fichiers de calcul](docs/NBO_ONE_JOB.fr.md) — inclut le modèle de travail dans l’arborescence des sources
- [Compiler depuis les sources (en anglais)](docs/BUILD.md)
- [Notes de la version stable](docs/releases/v0.3.0.md) · [Notes de la préversion](docs/releases/v0.4.0-pre.1.md) · [Anciennes préversions v0.3](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0-pre-archive)

[Signaler un problème ou proposer une fonction](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose).

## Licence

[Licence Apache 2.0](LICENSE). Les licences des bibliothèques fournies avec l’application figurent dans les [mentions relatives aux composants tiers](THIRD_PARTY_NOTICES.md).
