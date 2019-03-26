[![Travis-CI Build Status](https://travis-ci.org/ahmohamed/lipidr.svg?branch=master)](https://travis-ci.org/ahmohamed/lipidr)

# Analysis workflow for targeted lipidomics.

`lipidr` implements a series of functions to facilitate inspection,
analysis and visualization of targeted lipidomics datasets. `lipidr`
takes exported Skyline CSV as input, allowing for multiple methods to
be analyzed together.

## Input
`lipidr` represents Skyline files as annotated data.frames, which can be easily manipulated by a wide variety of R packages. Sample annotations, such as sample group or other clinical information can be easily loaded.

## Quality control & plotting
`lipidr` generates various plots, such as PCA and box plots, for quality control of samples and measured lipids. Normalization methods with and without internal standards are also supported.

## Differential Analysis
Differential analysis can be performed using any of the loaded clinical variables, which can be readily visualized as volcano plots. A novel lipid set enrichment analysis is implemented to detect preferential regulation of certain lipid classes, chain lengths or saturation patterns. Plots for visualization of enrichment results are also implemented.
