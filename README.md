# Spatial Correlation of Gene Expressions


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
