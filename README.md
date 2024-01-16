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
The input and output MESA data are available at https://zenodo.org/records/10479754. For the local article compilation, it is necessary to unzip the `data` folder from the Zenodo archive and put its content in the `../src/data` folder.

All file dependencies for reproducing the plots are available in `showyourwork.yml` file.
Use `conda 23.11.0` to build.

An open source scientific article created using the [showyourwork](https://github.com/showyourwork/showyourwork) workflow.
