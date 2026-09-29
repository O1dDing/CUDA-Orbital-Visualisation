# Chemical Orbital Visualiser (COV)

[English](README.md) · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · **Français**

Explorez les énergies et les occupations des orbitales, ainsi que les liens entre orbitales atomiques et moléculaires.

COV présente les énergies et les occupations des orbitales dans des diagrammes de niveaux d’énergie interactifs. Vous pouvez exporter les figures et les données. Il lit les fichiers Gaussian FCHK/FCH et Molden.

La préversion v0.4 ajoute l’analyse NBO. Avec la sortie NBO du calcul, vous pouvez voir comment les orbitales atomiques et localisées contribuent aux orbitales moléculaires, suivre leurs liens dans le diagramme et examiner les charges et les informations sur les liaisons.

Sélectionnez une orbitale pour voir sa forme en 3D. Le rendu utilise NVIDIA CUDA. L’interface est disponible en anglais, chinois simplifié, japonais et français.

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

Les téléchargements Windows actuels sont compilés pour les GPU NVIDIA GeForce de la série RTX 50. Les autres architectures NVIDIA nécessitent une version compilée adaptée ; voir [Compiler depuis les sources](docs/BUILD.md). Installez un pilote NVIDIA compatible. Le CUDA Toolkit n’est nécessaire que pour compiler depuis les sources.

## Premiers pas

1. Téléchargez et extrayez une archive ZIP pour Windows, puis lancez `cov.exe`.
2. Ouvrez un fichier Gaussian FCHK/FCH ou un fichier Molden compatible, ou glissez-le dans la fenêtre.
3. Parcourez les énergies et les occupations dans la liste des orbitales et le diagramme. Sélectionnez un niveau pour voir son orbitale, ou exportez le diagramme.
4. Dans la préversion v0.4, ouvrez le dossier du calcul pour charger ensemble les fichiers de fonction d’onde et NBO. Si le dossier contient plusieurs calculs, choisissez celui à ouvrir.

Avant d’ouvrir un fichier Gaussian `.chk`, convertissez-le en `.fchk` avec l’utilitaire Gaussian `formchk`. Pour les fonctions NBO, utilisez des fichiers correspondant à la même géométrie et au même état électronique que la fonction d’onde. COV lit les résultats ; il n’exécute pas Gaussian ni NBO.

## Fichiers d’entrée et configuration requise

- **Fonctions d’onde :** fichiers Gaussian `.fchk` / `.fch`, ou fichiers compatibles `.molden` / `.mol` / `.input`. Les fichiers Gaussian `.chk` doivent d’abord être convertis.
- **Fichiers NBO — préversion v0.4 :** le rapport, l’archive et les matrices orbitalaires fournissent différentes parties de l’analyse. Pour savoir quels fichiers sont nécessaires à chaque vue, consultez [Utiliser les résultats NBO (en anglais)](docs/AOMO_NBO.md).
- **Taille des molécules :** jusqu’à 100 atomes par fichier d’entrée.
- **Affichage :** un GPU NVIDIA, un pilote compatible et OpenGL 2.1 ou une version ultérieure. Utilisez une version compilée pour l’architecture de votre GPU.

COV lit des résultats de calcul existants. Les étiquettes attribuées automatiquement aux orbitales et aux liaisons peuvent être erronées ; consultez la sortie d’origine pour les interpréter.

## Documentation

- [Utiliser COV](docs/UI.md)
- [Utiliser les résultats NBO (en anglais)](docs/AOMO_NBO.md) — préversion v0.4
- [Compiler depuis les sources](docs/BUILD.md)
- [Notes de la version stable](docs/releases/v0.3.0.md) · [Notes de la préversion](docs/releases/v0.4.0-pre.1.md)

## Licence

[Licence Apache 2.0](LICENSE). Les licences des bibliothèques fournies avec l’application figurent dans les [mentions relatives aux composants tiers](THIRD_PARTY_NOTICES.md).
