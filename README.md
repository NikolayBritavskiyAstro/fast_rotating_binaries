<p align="center">
<a href="https://github.com/showyourwork/showyourwork">
<img width = "450" src="https://raw.githubusercontent.com/showyourwork/.github/main/images/showyourwork.png" alt="showyourwork"/>
</a>
<br>
<br>
<a href="https://github.com/NikolayBritavskiyAstro/fast_rotating_binaries/actions/workflows/build.yml">
<img src="https://github.com/NikolayBritavskiyAstro/fast_rotating_binaries/actions/workflows/build.yml/badge.svg?branch=main" alt="Article status"/>
</a>
<a href="https://github.com/NikolayBritavskiyAstro/fast_rotating_binaries/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/NikolayBritavskiyAstro/fast_rotating_binaries/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
</p>


This repository contains all files for reproducing Britavskiy et al. (2024) work regarding tracing the evolution of short-period fast-rotating massive binaries.
The input and output MESA data are available at https://zenodo.org/records/10479754. 



##### How to build the article locally?

1. Install conda (`version 23.11.0` ) and [showyourwork](https://show-your.work/en/latest/install/) workflow. 
2. Clone this repository and run `showyourwork build`.
3. The first time this will download the `data` folder (~4 GB) from [Zenodo](https://zenodo.org/records/10479754) and will unzip all the required data to `../src/data` folder. It is also possible to download manually the `data` folder and put its content (~20 GB) in the `../src/data` folder.

All file dependencies for reproducing the plots are available in `showyourwork.yml` file.

An open source scientific article created using the [showyourwork](https://github.com/showyourwork/showyourwork) workflow.
