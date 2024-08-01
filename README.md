
# rectalNRT_deconv

## Overview

This repository contains scripts and data for analyzing single-cell RNA sequencing data from rectal cancer studies. The primary focus is on deconvolution and cell-cell communication analysis using various R packages.

## Repository Structure

```
    ├── result/
    ├── plot/
    ├── dataset/
    ├── data-analysis/
    ├── .Rproj.user/
    ├── .test_data/
    ├── README.md
    ├── requirements.txt
    ├── ref.txt
    ├── .Rhistory
    └── .gitattributes
```

### Key Directories and Files

- `data_codes/20231001/`: Contains scripts, data, and results for the analysis.
  - `.Rproj.user/`: R project configuration files.
  - `dataset/`: Raw and processed datasets used for the analysis.
  - 'data-analysis/': For cellular interactions in the dataset isobolographic letter analysis
  - `output/`: Generated output files from the analysis.
  - `plot/`: Scripts and files for plotting results.
  - `result/`: Results of various analyses.

- `README.md`: This file.

## Prerequisites

Make sure you have the following R packages installed:
- CellChat
- patchwork
- stringr
- dplyr
- ggplot2

You can install these packages using the following commands in R:
```R
install.packages(c("stringr", "patchwork", "dplyr", "ggplot2"))
# Install CellChat from GitHub
devtools::install_github("sqjin/CellChat")
```

## Example running the Analysis

1. **Data Preprocessing**

   Load the datasets and preprocess them as necessary. Example script:
   ```R
   load("./datasets/cellchat/GSE132465_tumor_original.RData")
   data.input = data
   colnames_data = str_split(colnames(data.input), '\.', simplify = TRUE)
   colnames(data.input) = paste0(colnames_data[,2], '-', colnames_data[,3])
   rownames(data.input) = str_split(rownames(data.input), '["]', simplify = TRUE)[,2]
   ```

2. **CellChat Analysis**

   Perform cell-cell communication analysis using the CellChat package. Example:
   ```R
   cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Cell_type")
   cellchat <- addMeta(cellchat, meta = meta)
   cellchat <- setIdent(cellchat, ident.use = "Cell_type")
   ```

3. **Results and Plots**

   Generate and save plots for the results. Example:
   ```R
   plotCommunication(cellchat)
   ggsave(filename = "communication_plot.png")
   ```

## Results

The results of the analysis are saved in the `result/` directory within `data_codes/20231001/`. This includes:
- Deconvolution results
- Network analysis files
- Plots and visualizations

## Contact

For any questions or further information, please contact:
- Xueling Li: xlli@cmpt.ac.cn, xuelingli16@foxmail.com
