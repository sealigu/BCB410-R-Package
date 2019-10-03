
# Ubiquitin Analysis

The goal of Ubiquitin Analysis is to analysis the relationships between two different ubiquitin proteins.

## Installation

``` r
require("devtools")
install_github("sealigu/UbiquitinAnalysis")
library("UbiquitinAnalysis"")
```

## Overview

# Function 1: PlotResModification()
![](./pic/)

Input two selected protein uniprot id (e.g. "Q9UHB7" and "Q9UKV5"), the output is a plot to show the amino acid modification position in the selected protein sequences. The output will help to compare whether these two proteins have the same modified position in both sequences.

# Function 2: PlotProteinInteractions()
https://github.com/sealigu/UbiquitinAnalysis/blob/master/pic/Protein%20Interaction%20Example.png

Input two selected protein uniprot id (e.g. "Q9UHB7" and "Q9UKV5"), the output is a plot to show the interacting proteins of the selected proteins. The output will help to compare whether these two proteins have the same interacting protein.

# Function 3: GetSimilarPercentage()
https://github.com/sealigu/UbiquitinAnalysis/blob/master/pic/Similar%20Percentage%20Example.png

Input two selected protein uniprot id (e.g. "Q9UHB7" and "Q9UKV5"), the output is a percentage number that shows the similarity of the two selected protein sequences.


## Example
``` r
library(Ubiquitin Analysis);
PlotResModification("Q9UHB7", "Q9UKV5");
PlotProteinInteractions("Q9UHB7", "Q9UKV5");
GetSimilarPercentage("Q9UHB7", "Q9UKV5");
```
- PlotResModification("Q9UHB7", "Q9UKV5")
- PlotProteinInteractions("Q9UHB7", "Q9UKV5")

The functions PlotResModification and PlotProteinInteractions were authored by Shiyun. Part of the code is used from StackOverFlow. The link and the authors are given in the code.

- GetSimilarPercentage("Q9UHB7", "Q9UKV5")

The function GetSimilarPercentage was authored by Shiyun. The Levenshtein distance algorithm was used to measuring the differences between two sequences.

- library(UniprotR)

The library UniprotR was used to retrieve data from UniProt services.
