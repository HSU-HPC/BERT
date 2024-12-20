BERT <img src="https://user-images.githubusercontent.com/81758255/236138668-c422b935-ed7f-4f2c-82a5-69503d8416f4.png" width="120px" align="right" />
===========

 BERT (Batch-Effect Reduction Trees) offers flexible and efficient batch-effect correction of *omics* data, while providing maximum tolerance to missing values. As such, BERT is a valuable preprocessing tool for data analysis workflows. By providing BERT via Bioconductor, we make this tool available to a wider research community. An accompanying research paper is currently under preparation and will be made public soon.

BERT addresses the same fundamental data integration challenges as [HarmonizR](https://github.com/HSU-HPC/HarmonizR) package, which has been released on Bioconductor in November 2023. However, various algorithmic modications and optimizations of BERT provide better execution time, better data coverage and enhanced flexibility compared to *HarmonizR*. Moreover, BERT offers a more user-friendly design and a less error-prone input format.
 
**Please note that our package _BERT_ is neither affiliated with nor related to _Bidirectional Encoder Representations from Transformers_ as published by Google.**

> This GitHub README provides only a brief introduction to BERT and we refer the reader to the [Bioconductor vignette](https://bioconductor.org/packages/release/bioc/html/BERT.html) for more details and more thorough explanations.

## System Requirements
BERT supports all major operating systems, i.e. Linux (e.g., Ubuntu 22.04 LTS), Microsoft Windows (e.g., Windows 10 and Windows 11) and macOS (e.g., Monterey and Ventura). Further, it has been tested to work on all major CPU architectures (x86_64, x64, arm64). The Bioconductor version requires R version 4.4. All other relevant software dependencies are specified in the source `DESCRIPTION` file along with their respective version numbers and will be installed automatically.

## Installation Guide
To install BERT, start R (version "4.4") and enter
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BERT")
```
The execution time for the installation may vary greatly depending on your bandwidth and latency of your internet connection, as well as pre-installed R packages. At maximum, we expect an installation time of 20 minutes.

## Example Usage
BERT provides functionality to generate simulated data with missing values and batch-effects. This data is correctly formatted for direct batch-effect correction using BERT.

```R
library(BERT)
dataset_raw <- generate_dataset(features=60, batches=10, samplesperbatch=10, mvstmt=0.1, classes=2)
dataset_corrected <- BERT(dataset_raw)
```

For this example, the average silhouette width (ASW) with respect to batch should decrease and vice versa for the ASW with respect to class label. At maximum, we expect a runtime of 20 seconds for the above example. 

## Usage
For details on how to use BERT, please refer to the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/BERT/inst/doc/BERT-Vignette.html).

## Issues
Please report any issues in the GitHub forum, the Bioconductor forum, or contact [the authors](mailto:yannis.schumann@desy.de,schlumbohm@hsu-hh.de) directly.

## License
This code is published under the GPLv3.0 License.

## Reference
Please cite our manuscript, if you use BERT for your research:
> High Performance Data Integration for Large-Scale Analyses of Incomplete Omic Profiles Using Batch-Effect Reduction Trees (BERT)