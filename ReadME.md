# Research Project 

#### This repository includes the files regarding the research project including the Supplementary Material. Please do not hesitate to contact me on david.mandia@yahoo.com

#### I have divided the repository into three sub-folders. One per each part of the overall reseach, The gene set enrichment analysis, pathway analysis, and gene network. The MirCLip+Data folder contains the result of direct targeting obtained through data manipulation and Webgestalt ( Gene Enrichment). Similarly, the PI data +Webgestalt result folder containd the result of gene expression analysis following miRNA knockout.
#### Lastly, gene network contains the data and the script for the gene network analysis using gene expression, GO term, and clustering.
#### See below for more details 

# Gene Network Analysis with R

##Please visit repository "miRNA-29-targetome-analysis-for-improved-wound-healing-and-keratinocytes-adhesion"  for the overall the project including data, figure, results 

This R script performs gene network analysis, focusing on identifying correlated genes and clustering them based on shared biological processes (GO terms). It utilizes various R packages for data manipulation, visualization, and analysis.

## Dependencies

Ensure you have the following R packages installed:

- igraph
- reactome.db
- biomaRt
- pRoloc
- AnnotationDbi

You can install these packages using the `install.packages()` function if they are not already installed. The script includes code to automatically install the packages if they are missing.

## Usage

1. **Run the Script**: Execute the R script in your R environment. The script assumes that you have the necessary data files and packages installed.

2. **Input Data**: The script reads input data from a CSV file named "fast_abc_nsa.csv" located in the "RP1 data" directory. Ensure that your data file follows the expected format.

3. **Analysis Steps**:
   - The script performs gene symbol conversion using Ensembl and maps genes to GO terms related to biological processes.
   - It calculates the correlation between genes based on numerical data columns and identifies gene pairs with significant correlations.
   - Gene pairs with shared biological processes are selected, and a graph is constructed to visualize the gene network.
   - Clustering is performed using the label propagation algorithm to identify communities of highly interconnected genes.
   - The most common GO terms within each cluster are identified and visualized on the graph.

4. **Output**: The script generates a JPEG image ("final_test.jpeg") containing the visual representation of the gene network with labeled clusters and annotated GO terms.

   #### Example of the output
<img src="https://raw.githubusercontent.com/davidmandia/miRNA-29-targetome-analysis-for-improved-wound-healing-and-keratinocytes-adhesion/main/gene%20network/final_test.jpeg" alt="gene network final test" width="400">
## Notes

- Adjustments may be needed based on the size of the input data and the complexity of the gene network.
- Ensure that the necessary data files are available and correctly referenced in the script.
