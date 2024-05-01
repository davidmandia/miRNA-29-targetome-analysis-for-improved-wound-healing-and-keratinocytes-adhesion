# Load required libraries
# Ensure all required packages are installed
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("reactome.db", quietly = TRUE)) install.packages("reactome.db")
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
if (!requireNamespace("pRoloc", quietly = TRUE)) install.packages("pRoloc")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) install.packages("AnnotationDbi")

# Load required libraries after installation
library(igraph)
library(reactome.db)
library(biomaRt)
library(pRoloc)
library(AnnotationDbi)


# Assuming 'data' is your dataframe

##The line below needs to be changed to fit the other data 
data <- read.csv("RP1 data/fast_abc_nsa.csv")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


# Perform gene symbol conversion
gene2entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                     filters = "external_gene_name",
                     values = data$gene_name,
                     mart = ensembl)

#Covert gene symobol to character
entrez_ids <- as.character(gene2entrez$entrezgene_id[!is.na(gene2entrez$entrezgene_id)])

# Map gene name with GO term biolobial process
GO_mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = entrez_ids, columns = c("GENENAME","SYMBOL", "GO"))
GO_mapping <- GO_mapping[GO_mapping$ONTOLOGY == "BP", ]

# Merge the mapping information back into your origentrez_idsinal dataframe
data_with_GO<- merge(data, GO_mapping, by.x = "gene_name", by.y = "SYMBOL", all.x = TRUE)
data_with_GO <- data_with_GO[!is.na(data_with_GO$GO), ]



## Take numerical columns and find correlation 

# Extract numeric columns
numeric_data <- data[, -1]

# Transpose the numeric data
transposed_data <- t(numeric_data)

# Create correlation matrix
cor_matrix <- cor(transposed_data, use = "pairwise.complete.obs")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA  # Set lower triangle (including diagonal) to NA

# Set correlation threshold  -- Changes depending on types of cells comparison  
cor_threshold <- 0.6

# In cor matrix only select genes with correlation between the thresholds  
cor_edges <- which(cor_matrix > cor_threshold & cor_matrix < 1, arr.ind = TRUE, useNames = FALSE)

cor_edges <- data.frame(cor_edges)
# Covert the matrix row and cols back to gene names 
convert_to_genes <- function(x){ 
  return (data$gene_name[x])} 

gene_pairs <- data.frame(lapply(cor_edges,convert_to_genes))

##Create emppty DF for BP comparison
gene_pairs[ , 'GO_gene1'] = NA
gene_pairs[ , 'GO_gene2'] = NA
gene_pairs[ , 'shared_GO'] = NA

# Columns to keep following interactions

cols_keep <- c("gene1","gene2", "GO_gene1", "GO_gene2", "shared_GO")
colnames(gene_pairs) <- cols_keep
#find common path and exclude NAs, do the graph 

#  assign  path of gene 1 
gene_pairs<- merge(gene_pairs, GO_mapping, by.x = "gene1", by.y = "SYMBOL")
gene_pairs$GO_gene1 <- gene_pairs$GO

gene_pairs <- gene_pairs[cols_keep]

# Same for path in gene  2 
gene_pairs<- merge(gene_pairs, GO_mapping, by.x = "gene2", by.y = "SYMBOL")
gene_pairs$GO_gene2 <- gene_pairs$GO
gene_pairs <- gene_pairs[cols_keep]

# Only select if both genes have path in common 
gene_pairs$shared_GO <- ifelse(gene_pairs$GO_gene1 == gene_pairs$GO_gene2, TRUE, FALSE)


gene_pairs <- na.omit(gene_pairs)
#print(gene_pairs_no_na[gene_pairs_no_na$shared_GO == TRUE, ])

gene_pairs <-  gene_pairs[gene_pairs$shared_GO == TRUE, ]

test_uniqueness_edge_no_dupli <- gene_pairs[!duplicated(gene_pairs),]


only_edges <-test_uniqueness_edge_no_dupli [c("gene1", "gene2")]
only_edges <- only_edges[only_edges$gene1 != 'MT1H', ]

# Create a graph
graph <- graph_from_data_frame(only_edges, directed = FALSE)
#Removes double edges (sometimes helps )
graph <-simplify(graph)  

#Clp algotithm clustering. See link on top for explantion 
communities <- cluster_label_prop(graph)


membership_vector <- as.integer(membership(communities))


gene_go_mapping <- GO_mapping[c("SYMBOL", "GO")]


## This whole part finds the most common BP within a clsuter.
# If GO BP already picked it moves to the next one (otherwise most terms were Gene expression)
unique_terms <- c("a")
class_labels <- data.frame(
  ClassLabel = sapply(unique(membership(communities)), function(cluster) {
    cluster_genes <- names(which(membership(communities) == cluster))
    cluster_go_terms <- gene_go_mapping$GO[gene_go_mapping$SYMBOL %in% cluster_genes]
    # Get the most common GO term in the cluster
    most_common_go <- names(sort(table(cluster_go_terms), decreasing = TRUE))[1:10]
    # print(most_common_go)
    for (go in most_common_go){
      if( go %in% unique_terms){
        next
      } else {
        # print("here")
        unique_terms[length(unique_terms) + 1] <- go
        #print(unique_terms)
        return(go)
        
      }
    }
    
  })
)

unique_terms <- c("a")
class_labels <- data.frame(
  ClassLabel = sapply(unique(membership(communities)), function(cluster) {
    cluster_genes <- names(which(membership(communities) == cluster))
    cluster_go_terms <- gene_go_mapping$GO[gene_go_mapping$SYMBOL %in% cluster_genes]
    # Get the most common GO term in the cluster
    most_common_go <- names(sort(table(cluster_go_terms), decreasing = TRUE))[1:10]
    for (go in most_common_go){
      #print(go)
      if( !go %in% unique_terms){
        unique_terms <<- c(unique_terms, go)  # Append the new term to unique_terms
        #print(unique_terms)
        return(go)
      }
      
    }
    return(most_common_go)
  })
)


#Convert GO terms ID to Pathways names
class_labels$term <- goIdToTerm(class_labels$ClassLabel, names = TRUE, keepNA = TRUE)


# Color nodes based on the most common GO term in each cluster


community_colors <- rainbow(length(unique(communities$membership)))
node_colors <- rainbow(length(unique(communities$membership)))

# Get unique GO terms in the data
unique_go_terms <- class_labels$term


# Amendments are needed depending on the size of graphs and number of genes 

jpeg(file="final_test.jpeg", width = 5000, height= 5000)
layout <-layout_nicely(graph)
plot(communities, graph, layout= layout,vertex.label.cex = 7, vertex.size = 15 )
legend("topleft", legend = unique_go_terms, fill = rainbow(length(unique_go_terms)), title = "GO Terms", cex = 7)

dev.off()

