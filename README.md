BERT <img src="https://user-images.githubusercontent.com/81758255/236138668-c422b935-ed7f-4f2c-82a5-69503d8416f4.png" width="120px" align="right" />
===========

 > BERT (Batch-Effect Removal with Trees) offers flexible and efficient batch effect correction of omics data, while providing maximum tolerance to missing values. Tested on multiple datasets from proteomic analyses, BERT offered a typical 5-10x runtime improvement over existing methods, while retaining more numeric values and preserving batch effect reduction quality.
 
As such, BERT is a valuable preprocessing tool for data analysis workflows, in particular for proteomic data. By providing BERT via Bioconductor, we make this tool available to a wider research community. An accompanying research paper is currently under preparation and will be made public soon.

BERT addresses the same fundamental data integration challenges than the [HarmonizR](https://github.com/HSU-HPC/HarmonizR) package, which has been released on Bioconductor in November 2023. However, various algorithmic modications and optimizations of BERT provide better execution time and better data coverage than HarmonizR. Moreover, BERT offers a more user-friendly design and a less error-prone input format.
 
**Please note that our package _BERT_ is neither affiliated with nor related to _Bidirectional Encoder Representations from Transformers_ as published by Google.**

## Usage
For details on how to use BERT, please refer to the [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/BERT/inst/doc/BERT-Vignette.html).

## Issues
Please report any issues in the GitHub forum or contact [the authors](mailto:schumany@hsu-hh.de,schlumbohm@hsu-hh.de) directly.

## License

This code is published under the GPLv3.0 License.

## Reference
Please cite our manuscript, if you use BERT for your research:
> Yannis Schumann, Simon Schlumbohm et al., BERT - Batch Effect Reduction Trees for Large-Scale Data Integration with Tolerance for Missing Values and Imbalanced Data, 2024 (in preparation)