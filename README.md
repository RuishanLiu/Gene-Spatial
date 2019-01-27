# Spatial Correlation of Gene Expressions

Codes for paper [Modeling Spatial Correlation of Transcripts With Application to Developing Pancreas](https://www.biorxiv.org/content/biorxiv/early/2018/08/14/391433.full.pdf).

## Process and Analyze Data

```
run train.py --path_file <path to data expression csv file> --path_nuc <Path to nuclei position csv file> --path_gene <Path to gene conversion csv file> --path_save <Path to save results>
```

By running train.py, the gene expression data are processed and the spatial correlation results are stored in a PKL file (as a SpatData class, defined in functions.model_spatdata).

## Plot Results

How to access and visualize results is illustrated in

```
read.ipynb
```
