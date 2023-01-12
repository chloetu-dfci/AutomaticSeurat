This is the README file for AutomaticSeurat.R

This script will take in several 10X samples, and automatically filter cells from each sample based on mitochondrial gene contamination, HBB gene contamination and standard deviation of features and count. Each sample is normalized using SCTransform, then all samples are integrated or merged together. The combined object is plotted using UMAP, where the given genes are visualized. The combined object is clustered using three resolutions. At each resolution, firstly, the proportion of cells in each cluster is visualized and saved as a CSV file. Secondly, the DEG in each cluster is saved as a CSV file.

Author/Point of contact: Chloe Tu, chloetu@ds.dfci.harvard.edu, TIGL


###---
Input: 
1. Multiple 10X's "filtered_feature_bc_matrix" folders containing barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

###---
Output:
1. PDF of figures generated during analysis
2. Rds object after analysis
3. Three csv files containing the differentially expressed genes in each cluster
4. Three csv files containing the number of cells in each cluster and sample

###---
Arguments:
--sampleNames: Names to give your samples. Samples must be separated by commas, no spaces in between.
--fileNames: File names of your samples. File names must be separated by commas, no spaces in between.
--outputDir: Directory where output PDF file and CSV files will go
--dataDir: Directory where input data files can be found
--pdfName: Name to give your pdf
--geneVec: A list of genes which expression you want to visualize in your samples. Genes must be separated by commas, no spaces in between.
--integrate: "true" if you want to integrate the samples together. "false" if you want to merge the samples. (default = false)
--species: Species of your samples. Valid arguments are "human" or "mouse"

###---
Example of how to run AutomaticSeurat.R on the Broad's server:
You need these R packages: Seurat, ggplot2, dplyr, tidyverse, gridExtra, grid, stringr, patchwork, Polychrome, Matrix, optparse. For reasons, Seurat is only found on the Broad's R-4.1 dotkit, so you must load R-4.1.

> ssh <yourLogin>@login.broadinstitute.org
> use UGER
> ish -l h_vmem=16G
> use R-4.1
> cd <thelocationOfAutomaticSeurat.R>
> Rscript AutomaticSeurat.R --sampleNames sample1,sample2,sample3 --fileNames sampleFolder1,sampleFolder2,sampleFolder3 --outputDir /home/unix/outputFolder/ --dataDir /home/unix/dataFolder/ --pdfName outputFigures.pdf --geneVec gene1,gene2,gene3 --integrate false --species mouse


###---
FAQ:
1) "Warning messages: Removed 2 rows containing missing values (geom_vline)" error while running script
Ignore it. 

2) "Warning messages: In is_quosure(x = expression) :
  restarting interrupted promise evaluation" error while running script
Ignore it.

3) "Killed" error while running script
Increase memory allocated to your interactive shell. Example:
> ish -l h_vmem=32G

4) Files aren't loading in properly? Ensure that:
a. The folders containing the "barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz" files for each sample must be named "filtered_feature_bc_matrix".
b. You have multiple samples, thus multiple "filtered_feature_bc_matrix" folders.
c. --dataDir must point to the directory containing the "filtered_feature_bc_matrix" folders.

5) Stuck on "Calculating cluster X"?
Don't worry. It can take a while.


