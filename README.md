BERT
===========
 > BERT (Batch-Effect Removal with Trees) offers flexible and efficient batch effect correction of omics data, while providing maximum tolerance to missing values. Tested on datasets from proteomic analyses, BERT offered a typical 5-10x runtime improvement over existing methods, while retaining more numeric values and preserving batch effect reduction quality.
 
## Installation
### Core Functionality
Please download and install a current version of R ([Windows binaries](https://cran.r-project.org/bin/windows/base/release.html)). You might want to consider installing a development environment as well, e.g. [RStudio](https://posit.co/downloads/).

After setting up R, first install the R package `devtools` by executing the following command in the R environment:
```R
install.packages("devtools")
```
Next, install the *Bioconductor* packages `sva` and `limma` by executing
```R
if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install("sva")
BiocManager::install("limma")
```
Then, download and unzip this GitHub repository. Finally, set the working directory of your R interpreter to the base directory[^1] of the package and install the package, e.g.
```R
setwd("/home/my_user/Downloads/BERT/BERT")
devtools::install()
```
Alternatively, BERT may also be installed directly from GitHub using 
```R
devtools::install_github("HSU-HPC/BERT/BERT")
```
### Additional Features
In order to use the MPI backend, the user should install the packages `Rmpi` and `doMPI` on their system. Note, that this requires a working MPI installation on your system.
```R
install.packages("Rmpi")
install.packages("doMPI")
```
## Data Preparation
As input, BERT requires a dataframe[^2] with samples in rows and features in columns. For each sample, the respective batch should be indicated by an integer in a corresponding column labelled *Batch*. Missing values should be labelled as `NA`. A valid example dataframe could look like this:
|  | Feature2 | Feature3| ... | Batch 
--- | --- | --- | --- | ---
 **Sample1** | 0.01 | 0.3  | ... | 1 
 **Sample2** | -0.01 | 0.1  | ... | 1 
 **Sample3** | 0.02 | 0.2  | ... | 2 
 **Sample3** | 0.01 | NA  | ... | 2 
 
 Note that each batch should contain at least two samples. Optional columns that can be passed are
 - `Label` A column with integers indicating the (known) class for each sample. `NA` is not allowed. BERT may use this columns and `Batch` to compute quality metrics after batch effect correction.
 - `Sample` A sample name. This column is ignored by BERT and can be used to provide meta-information for further processing.
 - `Cov_1`, `Cov_2`, ..., `Cov_x`: One or multiple columns with integers, indicating one or several covariate levels. `NA` is not allowed. If this(these) column(s) is present, BERT will pass them as covariates to the the underlying batch effect correction method. As an example, this functionality can be used to preserve differences between healthy/tumorous sample, if some of the batches exhibit strongly variable class distributions. Note that BERT requires at least two numeric values per batch and unique covariate level to adjust a feature. Features that don't satisfy this condition in a specific batch are set to `NA` for that batch. 
- `Reference` A column with integers from $\mathbb{N}_0$ that indicate, whether a sample should be used for "learning" the transformation for batch effect correction or whether the sample should be co-adjusted using the learned transformation from the other samples. `NA` is not allowed. This feature can be used, if some batches contain unique classes or samples with unknown classes -- prohibiting the usage of the covariate columns. If the column contains a `0` for a sample, this sample will be co-adjusted. Otherwise, the sample should contain the respective class (encoded as integer). Note that BERT requires at least two references of common class per adjustment step and that the `Reference` column is mutually exclusive with covariate columns.

## Basic Usage
BERT can be invoked by importing the `BERT` library and calling the `BERT` function. The batch effect corrected data is returned as a dataframe that mirrors the input dataframe[^3].
```R
library(BERT)
# generate test data with 10% missing values as provided by the BERT library
dataset_raw <- generateDataset(features=600, batches=10, samplesperbatch=10, mvstmt=0.1, classes=2)
# apply BERT
dataset_adjusted <- BERT(dataset_raw)
```

```
2023-05-03 15:54:56 INFO::Formatting Data.
2023-05-03 15:54:56 INFO::Replacing NaNs with NAs.
2023-05-03 15:54:56 INFO::Removing potential empty rows and columns
2023-05-03 15:54:57 INFO::Found  6000  missing values.
2023-05-03 15:54:57 INFO::Introduced  0  missing values due to singular proteins at batch/covariate level.
2023-05-03 15:54:57 INFO::Done
2023-05-03 15:54:57 INFO::Acquiring quality metrics before batch effect correction.
2023-05-03 15:54:57 INFO::Starting hierarchical adjustment
2023-05-03 15:54:57 INFO::Found  10  batches.
2023-05-03 15:54:57 INFO::Adjusting the last 10 batches sequentially
2023-05-03 15:54:57 INFO::Adjusting sequential tree level 1 with 10 batches
2023-05-03 15:54:58 INFO::Adjusting sequential tree level 2 with 5 batches
2023-05-03 15:54:59 INFO::Adjusting sequential tree level 3 with 3 batches
2023-05-03 15:54:59 INFO::Adjusting sequential tree level 4 with 2 batches
2023-05-03 15:54:59 INFO::Done
2023-05-03 15:54:59 INFO::Acquiring quality metrics after batch effect correction.
2023-05-03 15:54:59 INFO::ASW Batch was 0.498198518557101 prior to batch effect correction and is now -0.0810976305927644 .
2023-05-03 15:54:59 INFO::ASW Label was 0.347010559909875 prior to batch effect correction and is now 0.775039508521937 .
2023-05-03 15:54:59 INFO::Total function execution time is  3.00484490394592  s and adjustment time is  1.84376502037048 s ( 61.36 )
```

BERT uses the  `logging` library to convey live information to the user during the adjustment procedure. The algorithm first verifies the shape and suitability of the input dataframe (lines 1-6) before continuiing with the actual batch effect correction (lines 8-15). BERT measure batch effects before and after the correction step by means of the average silhouette score (ASW) with respect to batch and labels (lines 7 and 16). The ASW Label should increase in a successful batch effect correction, whereas low values ($\leq 0$) are desireable for the ASW Batch[^4]. Finally, BERT prints the total function execution time (including the evaluation of the quality metrics).

## Advanced Options


[^1]: The base directory contains the folders *man*,*R* and *tests*. 
[^2]: Matrices work as well, but will automatically be converted to dataframes.
[^3]: In particular, the row and column names are in the same order and the optional columns are preserved.
[^4]: The optimum of ASW Label is 1, which is typically however not achieved on real-world datasets. Also, the optimum of ASW Batch can vary, depending on the class distributions of the batches.
