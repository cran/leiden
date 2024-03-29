## ----setup, include = FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("reticulate")
module <- py_available() && reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  if (!requireNamespace("devtools"))
#      install.packages("devtools")
#  devtools::install_github("TomKellyGenetics/leiden")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  install.packages("leiden")

## ----eval=FALSE, include=FALSE----------------------------------------------------------------------------------------
#  install.packages("leiden",  quiet = TRUE, repos = 1)
#  devtools::install_github("TomKellyGenetics/leiden", ref = "dev")

## ---------------------------------------------------------------------------------------------------------------------
library("leiden")

## ---------------------------------------------------------------------------------------------------------------------
adjacency_matrix <- rbind(cbind(matrix(round(rbinom(400, 1, 0.8)), 20, 20),
                                matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.8)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.2)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.9)), 20, 20)))
str(adjacency_matrix)
dim(adjacency_matrix )

## ---------------------------------------------------------------------------------------------------------------------
library("igraph")
rownames(adjacency_matrix) <- 1:60
colnames(adjacency_matrix) <- 1:60
graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
graph_object

## ----warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
plot(graph_object, vertex.color = "grey75")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  library("igraph")
#  adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(rep(1, 20), rep(2, 20), rep(3, 20))

## ----eval=module------------------------------------------------------------------------------------------------------
partition <- leiden(adjacency_matrix)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ----eval=module------------------------------------------------------------------------------------------------------
#run with defaults
  partition <- leiden(adjacency_matrix)


#run with ModularityVertexPartition"
  partition <- leiden(adjacency_matrix, partition_type = "ModularityVertexPartition")


#run with resolution parameter
  partition <- leiden(adjacency_matrix, resolution_parameter = 0.95)

## ----warning=FALSE, message=FALSE, eval=module------------------------------------------------------------------------
partition <- leiden(adjacency_matrix, resolution_parameter = 0.5)
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ----warning=FALSE, message=FALSE, eval=module------------------------------------------------------------------------
partition <- leiden(adjacency_matrix, resolution_parameter = 1.8)
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ----warning=FALSE, message=FALSE, eval=module------------------------------------------------------------------------
partition <- leiden(adjacency_matrix, max_comm_size = 12)
node.cols <- brewer.pal(min(c(9, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ----warning=FALSE, message=FALSE, eval=module------------------------------------------------------------------------
#generate example weights
weights <- sample(1:10, sum(adjacency_matrix!=0), replace=TRUE)
partition <- leiden(adjacency_matrix, weights = weights)
table(partition)

## ----warning=FALSE, message=FALSE, eval=module------------------------------------------------------------------------
#generate example weighted matrix
adjacency_matrix[adjacency_matrix == 1] <- weights
partition <- leiden(adjacency_matrix)
table(partition)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  adjacency_matrix <- as.matrix(object@snn)
#  membership <- leiden(adjacency_matrix)
#  object@ident <- as.factor(membership)
#  names(object@ident) <- rownames(object@meta.data)
#  object@meta.data$ident <- as.factor(membership)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  library("Seurat")
#  FindClusters(pbmc_small)
#  membership <- leiden(pbmc_small@graphs$RNA_snn)
#  table(membership)

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  FindClusters(pbmc_small, algorithm = "leiden")
#  table(pbmc_small@active.ident)

