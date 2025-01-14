# BERT: Batch-Effect Reduction Trees

[![Build Status](https://bioconductor.org/shields/build/release/bioc/BERT.svg)]()
[![Supported Platforms](https://bioconductor.org/shields/availability/release/BERT.svg)]()
[![Bioconductor Availability](https://bioconductor.org/shields/years-in-bioc/BERT.svg)]()
[![Last Update](https://bioconductor.org/shields/lastcommit/release/bioc/BERT.svg)]()


> Data from high-throughput technologies assessing global patterns of biomolecules (*omic* data), is often afflicted with missing values and with measurement-specific biases (batch-effects), that hinder the quantitative comparison of independently acquired datasets. This repository provides the BERT algorithm, a high-performance method for data integration of incomplete omic profiles.

> [!IMPORTANT]
> This repository is primarily intended for development purposes. For typical users, BERT is provided via [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/BERT.html). Note that repository badges refer to the release version of BERT, which may be multiple commits behind the source code provided here. The latest CI/CD results for BERT may be obtained [here](https://www.bioconductor.org/packages/devel/bioc/html/BERT.html).

> [!WARNING]
> The R package provided here is neither affiliated with nor related to Bidirectional Encoder Representations from Transformers as published by Devlin et al in 2019 (_arXiv:1810.04805_).

# Installation

> [!TIP]
> It is recommended to install BERT via Bioconductor as described [here](https://www.bioconductor.org/packages/release/bioc/html/BERT.html).

For development purposes, the BERT package can be installed directly from this repository using _devtools_.

```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('S4Vectors', 'S4Arrays', 'XVector', 'genefilter', 'SparseArray'))
devtools::install_github('HSU-HPC/BERT')
```

Please compare the installed version of R to the required version for Bioconductor and install all build dependencies if compilation from source is required for your target[^1].


# Usage

The BERT library is designed to offer high user friendliness whilst providing maximum flexibility. The following example demonstrates how to use the software on a simulated dataset with batch-effects and missing values:


```R
# import library
library(BERT)
# simulate dataset with 10% missing values
dataset_raw <- generate_dataset(features=60, batches=10, samplesperbatch=10, mvstmt=0.1, classes=2)
# apply BERT with default arguments
dataset_corrected <- BERT(dataset_raw)
```

> [!TIP]
> A detailed explanation of all available parameters, their default values and optimal configurations for typical scenarios can be found in the [Bioconductor vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/BERT/inst/doc/BERT-Vignette.html).

# Support

Users may ask for assistance via the [Bioconductor support site](https://support.bioconductor.org/tag/bert/). Bug reports may be filed via the [Issues](https://github.com/HSU-HPC/BERT/issues) tab of this repository. For confidential or security-related problems, please send an email to 

_yannis_ [dot] _schumann_ [at] _desy_ [dot] _de_ .

# License

This code is published under the GPLv3.0 License.

# References

Citations make research visible. If you use BERT for your research, please cite the following publication:

- Computational Methods for Data Integration and Imputation of Missing Values in Omics Datasets, Y. Schumann Gocke / A. Gocke / J. E. Neumann, 2024-12 PROTEOMICS, Wiley, [https://doi.org/10.1002/pmic.202400100](https://doi.org/10.1002/pmic.202400100)

[^1]: On Ubuntu 24.04, a complete list of depencies would be: _wget_, _curl _, _build-essential_, _libssl-dev_, _libcurl4-openssl-dev_, _pkg-config_, _git_, _ca-certificates_, _libxml2_, _libxml2-dev_, _gnupg_, _software-properties-common_, _libfontconfig1-dev_, _libharfbuzz-dev_, _libfribidi-dev_, _libfreetype6-dev_, _libpng-dev_, _libtiff5-dev_, _libjpeg-dev_
