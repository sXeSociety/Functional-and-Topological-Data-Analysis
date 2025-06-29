# Functional and Topological Data Analysis

This repository contains the full R code, paper, presentation, plots, and saved objects for my experimental project on Functional and Topological Data Analysis (FTDA) using Continuous Glucose Monitoring (CGM) data from adults with Type 1 Diabetes.

## Contents
- `FTDA_code - Corradini Andrea - 46953A.R`: complete pipeline, from preprocessing to regression and topological analysis
- `Plots/`: PNG images of all visual outputs
- `Objects/`: saved RDS objects (e.g., fitted models, cross-validation results)
- `Corradini Andrea - Functional and Topological Data Analysis - 46953A.pdf`: the paper
- `Presentation - Corradini Andrea - 46953A.pdf`: the PowerPoint presentation
- `README.md`: this project description

## Main techniques
- **Functional Data Analysis (FDA)**: FPCA, regression (scalar-on-function, function-on-scalar), clustering on FPCA scores, depth measures (MBD, BD2), functional boxplots
- **Functional ANOVA (FANOVA)**: functional ANOVA by clusters and HbA1c groups
- **Topological Data Analysis (TDA)**: persistence diagrams (alpha complex), bottleneck distance computation

## Data source
The CGM datasets and related clinical tables used in this project are not included in this repository due to GitHub size limitations.  
They can be obtained from:

- [Awesome-CGM Aleppo (2017) dataset](https://github.com/IrinaStatsLab/Awesome-CGM/wiki/Aleppo-(2017))

## Author
Andrea Corradini  
Università degli Studi di Milano
