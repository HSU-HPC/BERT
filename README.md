BERT <img src="https://user-images.githubusercontent.com/81758255/236138668-c422b935-ed7f-4f2c-82a5-69503d8416f4.png" width="120px" align="right" />
===========

 > BERT (Batch-Effect Removal with Trees) offers flexible and efficient batch effect correction of omics data, while providing maximum tolerance to missing values. Tested on multiple datasets from proteomic analyses, BERT offered a typical 5-10x runtime improvement over existing methods, while retaining more numeric values and preserving batch effect reduction quality.
 
As such, BERT is a valuable preprocessing tool for data analysis workflows, in particular for proteomic data. By providing BERT via Bioconductor, we make this tool available to a wider research community. An accompanying research paper is currently under preparation and will be made public soon.

BERT addresses the same fundamental data integration challenges than the [HarmonizR][https://github.com/HSU-HPC/HarmonizR] package, which is released on Bioconductor in November 2023. However, various algorithmic modications and optimizations of BERT provide better execution time and better data coverage than HarmonizR. Moreover, BERT offers a more user-friendly design and a less error-prone input format.
 
**Please note that our package _BERT_ is neither affiliated with nor related to _Bidirectional Encoder Representations from Transformers_ as published by Google.**

## Installation
## Core Functionality
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
BERT can then be installed directly from GitHub using 
```R
devtools::install_github("HSU-HPC/BERT")
```
Alternatively, download and unzip this GitHub repository. Finally, set the working directory of your R interpreter to the base directory[^1] of the package and install the package, e.g.
```R
setwd("/home/my_user/Downloads/BERT-main/BERT-main")
devtools::install()
```
BERT has also been submitted to Bioconductor. Once it has been accepted, this page will be updated and BERT can then be installed via
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BERT")
```
## Additional Features
In order to use the MPI backend, users should install the packages `Rmpi` and `doMPI` on their system. Note, that this requires a working MPI installation on your system.
```R
install.packages("Rmpi")
install.packages("doMPI")
```
# Data Preparation
As input, BERT requires a dataframe[^2] with samples in rows and features in columns. For each sample, the respective batch should be indicated by an integer or string in a corresponding column labelled *Batch*. Missing values should be labelled as `NA`. A valid example dataframe could look like this:
|  | Feature2 | Feature3| ... | Batch 
--- | --- | --- | --- | ---
 **Sample1** | 0.01 | 0.3  | ... | 1 
 **Sample2** | -0.01 | 0.1  | ... | 1 
 **Sample3** | 0.02 | 0.2  | ... | 2 
 **Sample3** | 0.01 | NA  | ... | 2 
 
 Note that each batch should contain at least two samples. Optional columns that can be passed are
 - `Label` A column with integers or strings indicating the (known) class for each sample. `NA` is not allowed. BERT may use this columns and `Batch` to compute quality metrics after batch effect correction.
 - `Sample` A sample name. This column is ignored by BERT and can be used to provide meta-information for further processing.
 - `Cov_1`, `Cov_2`, ..., `Cov_x`: One or multiple columns with integers, indicating one or several covariate levels. `NA` is not allowed. If this(these) column(s) is present, BERT will pass them as covariates to the the underlying batch effect correction method. As an example, this functionality can be used to preserve differences between healthy/tumorous samples, if some of the batches exhibit strongly variable class distributions. Note that BERT requires at least two numeric values per batch and unique covariate level to adjust a feature. Features that don't satisfy this condition in a specific batch are set to `NA` for that batch. 
- `Reference` A column with integers or strings from $\mathbb{N}_0$ that indicate, whether a sample should be used for "learning" the transformation for batch effect correction or whether the sample should be co-adjusted using the learned transformation from the other samples. `NA` is not allowed. This feature can be used, if some batches contain unique classes or samples with unknown classes which would prohibit the usage of covariate columns. If the column contains a `0` for a sample, this sample will be co-adjusted. Otherwise, the sample should contain the respective class (encoded as integer). Note that BERT requires at least two references of common class per adjustment step and that the `Reference` column is mutually exclusive with covariate columns.

Note that BERT only allows `SummarizedExperiment`s with only one assay. All metadata information, including the mandatory batch information, must be contained in the metadata, which BERT accesses using `colData`. For instance, a valid `SummarizedExperiment` might be defined as
```R
nrows <- 200
ncols <- 8
expr_values <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# colData also takes further metadata information, such as Label, Sample,
# Reference or Covariables
colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2), Reference=c(1,1,0,0,1,1,0,0))
dataset_raw = SummarizedExperiment::SummarizedExperiment(assays=list(expr=expr_values), colData=colData)
```


# Basic Usage
BERT can be invoked by importing the `BERT` library and calling the `BERT` function. The batch effect corrected data is returned as a dataframe that mirrors the input dataframe[^3].
```R
library(BERT)
# generate test data with 10% missing values as provided by the BERT library
dataset_raw <- generate_dataset(features=600, batches=10, samplesperbatch=10, mvstmt=0.1, classes=2)
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

BERT uses the  `logging` library to convey live information to the user during the adjustment procedure. The algorithm first verifies the shape and suitability of the input dataframe (lines 1-6) before continuing with the actual batch effect correction (lines 8-15). BERT measure batch effects before and after the correction step by means of the average silhouette score (ASW) with respect to batch and labels (lines 7 and 16). The ASW Label should increase in a successful batch effect correction, whereas low values ($\leq 0$) are desireable for the ASW Batch[^4]. Finally, BERT prints the total function execution time (including the computation time for the quality metrics).

# Advanced Options

## Parameters

BERT offers a large number of parameters to customize the batch effect adjustment. The full function call, including all defaults is
```R
BERT(data, cores = 1, combatmode = 1,method = "ComBat",qualitycontrol = TRUE,verify = TRUE,mpi = FALSE,stopParBatches = 4,corereduction = 2,backend = "default", labelname="Label", batchname="Batch", referencename="Reference", covariatename=NULL)
```
In the following, we list the respective meaning of each parameter:
- `data`: The input dataframe or matrix to adjust. See [Data Preparation](#data-preparation) for detailed formatting instructions.
- `method`: The method to use for the underlying batch effect correction steps. Should be either `ComBat`, `limma` for `limma::removeBatchEffects` or `ref` for adjustment using specified references (cf. [Data Preparation](#data-preparation)). The underlying batch effect adjustment method for `ref` is a modified version of the `limma` method.
- `combatmode`: An integer that encodes the parameters to use for ComBat.

  | Value | par.prior | mean.only 
  | --- | --- | ---
  | 1 | TRUE | FALSE
  | 2 | TRUE | TRUE
  | 3 | FALSE | FALSE
  | 4 | FALSE | TRUE

  The value of this parameter will be ignored, if `method!="ComBat"`.
- `qualitycontrol`: A boolean to (de)activate the ASW computation. Deactivating the ASW computations accelerates the computations.
- `verify`: A boolean to (de)activate the initial format check of the input data. Deactivating this verification step accelerates the computations.
- `cores`: The number of cores (processes) to use for parallel adjustment. Increasing this parameter can speed up the batch effect adjustment considerably, in particular for large datasets. If possible, the processes are spawned by forking -- otherwise, BERT uses `PSOCKCluster`. A value between $2$ and $4$ is a reasonable choice for typical commodity hardware.
- `stopParBatches`: Positive integer indicating the minimum number of batches required at a hierarchy level to proceed with parallelized adjustment. If the number of batches is smaller, adjustment will be performed sequentially to avoid communication overheads.
- `corereduction`: Positive integer indicating the factor by which the number of processes should be reduced, once no further adjustment is possible for the current number of batches.[^5]
- `mpi`: A boolean to (de)activate the MPI backend. If `TRUE`, this will replace the default `ForkCluster` or `PSOCKCluster`. *Cores must set to the total number of processes when using the MPI backend.*
- `backend`: The backend to use for inter-process communication. Possible choices are `default` and `file`, where the former refers to the default communication backend of the requested parallelization mode and the latter will create temporary `.rds` files for data communication. 'default' is usually faster for small to medium sized datasets.
- `labelname`: A string containing the name of the column to use as class labels. The default is "Label".
- `batchname`: A string containing the name of the column to use as batch labels. The default is "Batch".
- `referencename`: A string containing the name of the column to use as reference labels. The default is "Reference".
- `covariatename`: A vector containing the names of columns with categorical covariables.The default is NULL, in which case all column names are matched agains the pattern "Cov".

## Verbosity
BERT utilizes the ```logging``` package for output. The user can easily specify the verbosity of BERT by setting the global logging level in the script. For instance

```R
logging::setLevel("WARN") # set level to warn and upwards
result <- BERT(data,cores = 1) # BERT executes silently
```


## Choosing the Optimal Number of Cores
BERT exhibits a large number of parameters for parallelisation as to provide users with maximum flexibility. For typical scenarios, however, the default parameters are well suited. For very large experiments ($>15$ batches), we recommend to increase the number of cores (a reasonable value is $4$ but larger values may be possible on your hardware). Most users should leave all parameters to their respective default.

# Examples

In the following, we present simple cookbook examples for BERT usage. Note that ASWs (and runtime) will most likely differ on your machine, since the data generating process involves multiple random choices.

## Sequential Adjustment with limma
Here, BERT uses limma as underlying batch effect correction algorithm (```method='limma'```) and performs all computations on a single process (```cores``` parameter is left on default).
```R
# import BERT
library(BERT)
# generate data with 30 batches, 600 features, 15 samples per batch, 15% missing values and 2 classes
dataset_raw <- generate_dataset(features=600, batches=20, samplesperbatch=15, mvstmt=0.15, classes=2)
# BERT
dataset_adjusted <- BERT(dataset_raw, method="limma")
```

```
2023-05-04 08:47:30 INFO::Formatting Data.
2023-05-04 08:47:30 INFO::Replacing NaNs with NAs.
2023-05-04 08:47:30 INFO::Removing potential empty rows and columns
2023-05-04 08:47:30 INFO::Found  12000  missing values.
2023-05-04 08:47:30 INFO::Introduced  0  missing values due to singular proteins at batch/covariate level.
2023-05-04 08:47:30 INFO::Done
2023-05-04 08:47:30 INFO::Acquiring quality metrics before batch effect correction.
2023-05-04 08:47:30 INFO::Starting hierarchical adjustment
2023-05-04 08:47:30 INFO::Found  20  batches.
2023-05-04 08:47:30 INFO::Adjusting the last 20 batches sequentially
2023-05-04 08:47:30 INFO::Adjusting sequential tree level 1 with 20 batches
2023-05-04 08:47:31 INFO::Adjusting sequential tree level 2 with 10 batches
2023-05-04 08:47:31 INFO::Adjusting sequential tree level 3 with 5 batches
2023-05-04 08:47:32 INFO::Adjusting sequential tree level 4 with 3 batches
2023-05-04 08:47:32 INFO::Adjusting sequential tree level 5 with 2 batches
2023-05-04 08:47:32 INFO::Done
2023-05-04 08:47:32 INFO::Acquiring quality metrics after batch effect correction.
2023-05-04 08:47:32 INFO::ASW Batch was 0.505805793208734 prior to batch effect correction and is now -0.123165780027324 .
2023-05-04 08:47:32 INFO::ASW Label was 0.300158617645049 prior to batch effect correction and is now 0.813535343433545 .
2023-05-04 08:47:32 INFO::Total function execution time is  1.79658007621765  s and adjustment time is  1.4183521270752 s ( 78.95 )
```

## Parallel Batch Effect Correction with ComBat
Here, BERT uses ComBat as underlying batch effect correction algorithm (```method``` is left on default) and performs all computations on a 2 processes (```cores=2```).
```R
# import BERT
library(BERT)
# generate data with 30 batches, 600 features, 15 samples per batch, 15% missing values and 2 classes
dataset_raw <- generate_dataset(features=600, batches=20, samplesperbatch=15, mvstmt=0.15, classes=2)
# BERT
dataset_adjusted <- BERT(dataset_raw, cores=2)
```

```
2023-05-04 08:51:31 INFO::Formatting Data.
2023-05-04 08:51:31 INFO::Replacing NaNs with NAs.
2023-05-04 08:51:31 INFO::Removing potential empty rows and columns
2023-05-04 08:51:31 INFO::Found  27000  missing values.
2023-05-04 08:51:32 INFO::Introduced  0  missing values due to singular proteins at batch/covariate level.
2023-05-04 08:51:32 INFO::Done
2023-05-04 08:51:32 INFO::Acquiring quality metrics before batch effect correction.
2023-05-04 08:51:32 INFO::Setting up cluster with  2  cores.
2023-05-04 08:51:32 INFO::Identified OS as Windows. Using Parallel Socket Cluster (PSOCK).
2023-05-04 08:51:32 INFO::Done
2023-05-04 08:51:32 INFO::Starting hierarchical adjustment
2023-05-04 08:51:32 INFO::Found  20  batches.
2023-05-04 08:51:32 INFO::Processing subtree level 1 with 20 batches using 2 cores.
2023-05-04 08:51:37 INFO::Adjusting the last 2 batches sequentially
2023-05-04 08:51:37 INFO::Adjusting sequential tree level 1 with 2 batches
2023-05-04 08:51:37 INFO::Done
2023-05-04 08:51:37 INFO::Stopping cluster gracefully.
2023-05-04 08:51:37 INFO::Done
2023-05-04 08:51:37 INFO::Acquiring quality metrics after batch effect correction.
2023-05-04 08:51:38 INFO::ASW Batch was 0.515034992184603 prior to batch effect correction and is now -0.108202226399361 .
2023-05-04 08:51:38 INFO::ASW Label was 0.315180231450925 prior to batch effect correction and is now 0.806086240878908 .
2023-05-04 08:51:38 INFO::Total function execution time is  6.19199800491333  s and adjustment time is  4.94132399559021 s ( 79.8 )
```

## Batch Effect Correction Using SummarizedExperiment
Here, BERT takes the input data using a ```SummarizedExperiment``` instead. Batch effect correction is then performed using ComBat as underlying algorithm (```method``` is left on default) and all computations are performed on a single process (```cores``` parameter is left on default).
```R
nrows <- 200
ncols <- 8
expr_values <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
# colData also takes further metadata information, such as Label, Sample,
# Reference or Covariables
colData <- data.frame(Batch=c(1,1,1,1,2,2,2,2))
dataset_raw = SummarizedExperiment::SummarizedExperiment(assays=list(expr=expr_values), colData=colData)
dataset_adjusted = BERT(dataset_raw)
```

```
2023-06-01 12:49:54.803452 INFO::Formatting Data.
2023-06-01 12:49:54.816053 INFO::Replacing NaNs with NAs.
2023-06-01 12:49:54.826527 INFO::Removing potential empty rows and columns
2023-06-01 12:49:55.027343 INFO::Found  27000  missing values.
2023-06-01 12:49:55.029614 INFO::BERT requires at least 2 numeric values per batch/covariate level. This may reduce the number of adjustable features considerably, depending on the quantification technique.
2023-06-01 12:49:55.338783 INFO::Introduced  0  missing values due to singular proteins at batch/covariate level.
2023-06-01 12:49:55.340449 INFO::Done
2023-06-01 12:49:55.341769 INFO::Acquiring quality metrics before batch effect correction.
2023-06-01 12:49:55.513225 INFO::Starting hierarchical adjustment
2023-06-01 12:49:55.515381 INFO::Found  20  batches.
2023-06-01 12:49:55.51778 INFO::Adjusting the last 20 batches sequentially
2023-06-01 12:49:55.52013 INFO::Adjusting sequential tree level 1 with 20 batches
2023-06-01 12:49:58.438396 INFO::Adjusting sequential tree level 2 with 10 batches
2023-06-01 12:49:59.006294 INFO::Adjusting sequential tree level 3 with 5 batches
2023-06-01 12:49:59.290571 INFO::Adjusting sequential tree level 4 with 3 batches
2023-06-01 12:49:59.475011 INFO::Adjusting sequential tree level 5 with 2 batches
2023-06-01 12:49:59.660278 INFO::Done
2023-06-01 12:49:59.662174 INFO::Acquiring quality metrics after batch effect correction.
2023-06-01 12:49:59.826449 INFO::ASW Batch was 0.507973184713031 prior to batch effect correction and is now -0.133957677497005 .
2023-06-01 12:49:59.828034 INFO::ASW Label was 0.326591041133387 prior to batch effect correction and is now 0.795799704493019 .
2023-06-01 12:49:59.829399 INFO::Total function execution time is  5.29782199859619  s and adjustment time is  4.14515495300293 s ( 78.24 )
```
## BERT with Covariables
BERT can utilize categorical covariables that are specified in columns ```Cov_1, Cov_2, ...```. These columns are automatically detected and integrated into the batch effect correction process.
```R
# import BERT
library(BERT)
# generate data with 30 batches, 600 features, 15 samples per batch, 15% missing values and 2 classes
dataset_raw <- generate_dataset(features=600, batches=20, samplesperbatch=15, mvstmt=0.15, classes=2)
# create covariable column with 2 possible values, e.g. male/female condition
dataset_raw["Cov_1"] = sample(c(1,2), size=dim(dataset_raw)[1], replace=TRUE)
# BERT
dataset_adjusted <- BERT(dataset_raw)
```

```
2023-05-04 09:04:38 INFO::Formatting Data.
2023-05-04 09:04:38 INFO::Replacing NaNs with NAs.
2023-05-04 09:04:38 INFO::Removing potential empty rows and columns
2023-05-04 09:04:38 INFO::Found  27000  missing values.
2023-05-04 09:04:38 INFO::BERT requires at least 2 numeric values per batch/covariate level. This may reduce the number of adjustable features considerably, depending on the quantification technique.
2023-05-04 09:04:38 INFO::Introduced  0  missing values due to singular proteins at batch/covariate level.
2023-05-04 09:04:38 INFO::Done
2023-05-04 09:04:38 INFO::Acquiring quality metrics before batch effect correction.
2023-05-04 09:04:39 INFO::Starting hierarchical adjustment
2023-05-04 09:04:39 INFO::Found  20  batches.
2023-05-04 09:04:39 INFO::Adjusting the last 20 batches sequentially
2023-05-04 09:04:39 INFO::Adjusting sequential tree level 1 with 20 batches
2023-05-04 09:04:39 INFO::Adjusting sequential tree level 2 with 10 batches
2023-05-04 09:04:40 INFO::Adjusting sequential tree level 3 with 5 batches
2023-05-04 09:04:41 INFO::Adjusting sequential tree level 4 with 3 batches
2023-05-04 09:04:41 INFO::Adjusting sequential tree level 5 with 2 batches
2023-05-04 09:04:41 INFO::Done
2023-05-04 09:04:41 INFO::Acquiring quality metrics after batch effect correction.
2023-05-04 09:04:41 INFO::ASW Batch was 0.493321585409235 prior to batch effect correction and is now -0.12597947286234 .
2023-05-04 09:04:41 INFO::ASW Label was 0.333404392178747 prior to batch effect correction and is now 0.803468947362174 .
2023-05-04 09:04:41 INFO::Total function execution time is  3.29042482376099  s and adjustment time is  2.47309303283691 s ( 75.16 )
```

## BERT with references
In rare cases, class distributions across experiments may be severely skewed. In particular, a batch might contain classes that other batches don't contain. In these cases, samples of common conditions may serve as references (*bridges*) between the batches (```method="ref"```). BERT utilizes those samples as references that have a condition specified in the "Reference" column of the input. All other samples are co-adjusted. Please note, that this strategy implicitly uses limma as underlying batch effect correction algorithm.

```R
# import BERT
library(BERT)
# generate data with 4 batches, 600 features, 15 samples per batch, 15% missing values and 2 classes
dataset_raw <- generate_dataset(features=600, batches=4, samplesperbatch=15, mvstmt=0.15, classes=2)
# create reference column with default value 0.  The 0 indicates, that the respective sample should be co-adjusted only.
dataset_raw[, "Reference"] <- 0
# randomly select 2 references per batch and class - in practice, this choice will be determined by external requirements (e.g. class known for only these samples)
batches <- unique(dataset_raw$Batch) # all the batches
for(b in batches){ # iterate over all batches
    # references from class 1
    ref_idx = sample(which((dataset_raw$Batch==b)&(dataset_raw$Label==1)), size=2, replace=FALSE)
    dataset_raw[ref_idx, "Reference"] <- 1
    # references from class 2
    ref_idx = sample(which((dataset_raw$Batch==b)&(dataset_raw$Label==2)), size=2, replace=FALSE)
    dataset_raw[ref_idx, "Reference"] <- 2
}
# BERT
dataset_adjusted <- BERT(dataset_raw, method="ref")
```

```
2023-05-04 09:16:05 INFO::Formatting Data.
2023-05-04 09:16:05 INFO::Replacing NaNs with NAs.
2023-05-04 09:16:05 INFO::Removing potential empty rows and columns
2023-05-04 09:16:05 INFO::Found  5400  missing values.
2023-05-04 09:16:05 INFO::Introduced  0  missing values due to singular proteins at batch/covariate level.
2023-05-04 09:16:05 INFO::Done
2023-05-04 09:16:05 INFO::Acquiring quality metrics before batch effect correction.
2023-05-04 09:16:05 INFO::Starting hierarchical adjustment
2023-05-04 09:16:05 INFO::Found  4  batches.
2023-05-04 09:16:05 INFO::Adjusting the last 4 batches sequentially
2023-05-04 09:16:05 INFO::Adjusting sequential tree level 1 with 4 batches
2023-05-04 09:16:05 INFO::Adjusting sequential tree level 2 with 2 batches
2023-05-04 09:16:05 INFO::Done
2023-05-04 09:16:05 INFO::Acquiring quality metrics after batch effect correction.
2023-05-04 09:16:05 INFO::ASW Batch was 0.544575830258036 prior to batch effect correction and is now -0.1478167400377 .
2023-05-04 09:16:05 INFO::ASW Label was 0.380573240992057 prior to batch effect correction and is now 0.923705947392327 .
2023-05-04 09:16:05 INFO::Total function execution time is  0.394418001174927  s and adjustment time is  0.241964817047119 s ( 61.35 )
```

## Issues
Please report any issues in the GitHub forum or contact [the authors](mailto:schumany@hsu-hh.de,schlumbohm@hsu-hh.de) directly.

## License

This code is published under the GPLv3.0 License.

## Reference
Please cite our manuscript, if you use BERT for your research:
> Yannis Schumann, Simon Schlumbohm et al., BERT - Batch Effect Reduction Trees with Tolerance to Missing Values, 2023

[^1]: The base directory contains the folders *man*,*R* and *tests*. 
[^2]: Matrices and SummarizedExperiments work as well, but will automatically be converted to dataframes.
[^3]: In particular, the row and column names are in the same order and the optional columns are preserved.
[^4]: The optimum of ASW Label is 1, which is typically however not achieved on real-world datasets. Also, the optimum of ASW Batch can vary, depending on the class distributions of the batches.
[^5]: E.g. consider a BERT call with 8 batches and 8 processes. Further adjustment is not possible with this number of processes, since batches are always processed in pairs. With `corereduction=2`, the number of processes for the following adjustment steps would be set to $8/2=4$, which is the maximum number of usable processes for this example.
