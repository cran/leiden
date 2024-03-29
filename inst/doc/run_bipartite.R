## ----setup, include = FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.cap = "", fig.path = "Plot")
knitr::opts_chunk$set(fig.align = "center")
options(width = 120, cli.unicode = FALSE, cli.width = 120)
library("reticulate")
module <- py_available() && reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

## ----eval=module------------------------------------------------------------------------------------------------------
library("leiden")

## ----eval=FALSE-------------------------------------------------------------------------------------------------------
#  remotes::install_github("schochastics/networkdata")
#  bipartite_graph <- networkdata::southern_women
#  bipartite_graph

## ----echo=FALSE-------------------------------------------------------------------------------------------------------
library("igraph")
suppressWarnings(suppressMessages({
  #imported from networkdata::southernwomen
  bipartite_graph <- structure(list(32, FALSE,
                                    c(18, 19, 20, 21, 22, 23, 25, 26, 18,
                                      19, 20, 22, 23, 24, 25, 19, 20, 21, 22, 23, 24, 25, 26, 18, 20,
                                      21, 22, 23, 24, 25, 20, 21, 22, 24, 20, 22, 23, 25, 22, 23, 24,
                                      25, 23, 25, 26, 22, 24, 25, 26, 24, 25, 26, 29, 25, 26, 27, 29,
                                      25, 26, 27, 29, 30, 31, 24, 25, 26, 27, 29, 30, 31, 23, 24, 26,
                                      27, 28, 29, 30, 31, 24, 25, 27, 28, 29, 30, 31, 25, 26, 27, 29,
                                      26, 28, 26, 28),
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                                      1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5,
                                      5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10,
                                      10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 13,
                                      13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15,
                                      15, 15, 16, 16, 17, 17),
                                    c(0, 8, 23, 1, 9, 15, 2, 10, 16, 24,
                                      30, 34, 3, 17, 25, 31, 4, 11, 18, 26, 32, 35, 38, 45, 5, 12,
                                      19, 27, 36, 39, 42, 70, 13, 20, 28, 33, 40, 46, 49, 63, 71, 78,
                                      6, 14, 21, 29, 37, 41, 43, 47, 50, 53, 57, 64, 79, 85, 7, 22,
                                      44, 48, 51, 54, 58, 65, 72, 86, 89, 91, 55, 59, 66, 73, 80, 87,
                                      74, 81, 90, 92, 52, 56, 60, 67, 75, 82, 88, 61, 68, 76, 83, 62,
                                      69, 77, 84),
                                    c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                      14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                      30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                                      46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
                                      62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
                                      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92),
                                    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                      3, 6, 12, 16, 24, 32, 42, 56, 68, 74, 78, 85, 89, 93),
                                    c(0, 8, 15, 23, 30, 34, 38, 42, 45, 49, 53, 57, 63, 70, 78, 85,
                                      89, 91, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93,
                                      93, 93),
                                    list(c(1, 0, 1), structure(list(), .Names = character(0)),
                                         list(type = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                       FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                       FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                                       TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                                              name = c("EVELYN", "LAURA", "THERESA", "BRENDA",
                                                       "CHARLOTTE", "FRANCES", "ELEANOR", "PEARL", "RUTH",
                                                       "VERNE", "MYRA", "KATHERINE", "SYLVIA", "NORA", "HELEN",
                                                       "DOROTHY", "OLIVIA", "FLORA", "E1", "E2", "E3", "E4",
                                                       "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
                                                       "E13", "E14")), list())), class = "igraph")
  bipartite_graph <- upgrade_graph(bipartite_graph)
}))
bipartite_graph

## ---------------------------------------------------------------------------------------------------------------------
table(as.numeric(V(bipartite_graph)$type))

## ----warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("graphsim")
node.cols <- c("palevioletred", "lightblue")[as.integer(V(bipartite_graph)$type)+1]
plot(bipartite_graph, vertex.color = node.cols, layout = layout.kamada.kawai)

## ---------------------------------------------------------------------------------------------------------------------
library("igraph")
bipartite_matrix <- igraph::as_adjacency_matrix(bipartite_graph)
bipartite_matrix 

## ----eval=module------------------------------------------------------------------------------------------------------
partition <- leiden(bipartite_matrix, partition_type = "CPMVertexPartition.Bipartite", resolution_parameter = 0.2, seed = 42)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
#                 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----eval=module------------------------------------------------------------------------------------------------------
partition <- leiden(bipartite_graph, partition_type = "CPMVertexPartition.Bipartite", resolution_parameter = 0.2, seed = 42)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
#                 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(bipartite_graph, vertex.color = node.cols, layout = layout.kamada.kawai)

## ----eval=module------------------------------------------------------------------------------------------------------
partition <- leiden(bipartite_graph, partition_type = "ModularityVertexPartition.Bipartite", resolution_parameter = 0.02, seed = 42)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(1, 1, 1, 1, 1, 6, 6, 3, 6, 7, 4, 2, 2, 2, 2, 4, 5, 5, 1, 1,
#                 1, 1, 6, 3, 7, 3, 5, 4, 5, 4, 2, 2)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(bipartite_graph, vertex.color = node.cols, layout = layout.kamada.kawai)

## ----eval=module, warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
partition <- leiden(bipartite_graph, partition_type = "ModularityVertexPartition.Bipartite", resolution_parameter = 0.005, seed = 42)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
#                 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----eval=module, warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(bipartite_graph, vertex.color = node.cols, layout = layout.kamada.kawai)

## ----eval=module, warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
partition <- leiden(bipartite_graph, partition_type = "ModularityVertexPartition.Bipartite", max_comm_size = 6, resolution_parameter = 0.025, seed = 42)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
#                 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----eval=module, warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(11, partition)),"Pastel1")[partition]
plot(bipartite_graph, vertex.color = node.cols, layout = layout.kamada.kawai)

## ----eval=module------------------------------------------------------------------------------------------------------
partition <- leiden(bipartite_graph, partition_type = "CPMVertexPartition.Bipartite", resolution_parameter = 0.02, seed = 42,
                    degree_as_node_size = TRUE)

## ----eval=!module,echo=FALSE, message=FALSE, warning=FALSE, results="hide"--------------------------------------------
#  partition <- c(1, 1, 1, 1, 1, 6, 6, 3, 6, 7, 4, 2, 2, 2, 2, 4, 5, 5, 1, 1,
#                 1, 1, 6, 3, 7, 3, 5, 4, 5, 4, 2, 2)

## ---------------------------------------------------------------------------------------------------------------------
table(partition)

## ----eval=module, warning=FALSE, message=FALSE, fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5----
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(bipartite_graph, vertex.color = node.cols, layout = layout.kamada.kawai)

