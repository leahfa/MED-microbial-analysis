# MED-microbial-analysis
R scripts used for processing microbial data for the paper "A microbiota targeted Mediterranean diet based lifestyle intervention- 
a feasibility pilot in healthy participants"

Notes:

Please start by installing and loading my package of custom R scripts for microbiome analysis:
library(devtools)
install_github("leahfa/lmic")
library(lmic)
Additionally, source the functions in functions_for_downstream_analysis.R

The file dada.R provides my workflow for processing of raw data, based on: https://benjjneb.github.io/dada2/tutorial.html
(and many thanks to Benjamin Callahan the dada2 developer)
Downstream  analysis scripts are subdivided into:
1. data_prepration_scripts - these should be run first, in this order:
 -data prepration and variable defibition.R - the strating point for all analysis, formats the metadata key, microbiome data, etc.
 -make_subsets_by_timepoint.R  - creates perfectly-matched subsets of metadata and microbiome data ,by timepoint;
 -Calculate deltas Vis1 Vis0.R - create delta matrices for microbiome , nutritional and clinical data
 2. data analysis scripts - these use the objects created in '1' 
 3. plotting_scripts - scripts to plot the results of '2'
 
 -


