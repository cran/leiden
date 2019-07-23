## ----setup, include = FALSE----------------------------------------------
library("leiden")
library("reticulate")
py_available()
module <- py_available() && py_numpy_available() && py_module_available("leidenalg") && py_module_available("igraph")
if(module){
  py_config()$python
  py_install("numpy")
  py_install("scipy")
  py_install("python-igraph")
  py_install("leidenalg")
  module <- py_module_available("leidenalg") && py_module_available("igraph")
}
if(module){
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)
}

## ------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = module
)

## ------------------------------------------------------------------------
print(module)

## ------------------------------------------------------------------------
paste(Sys.info()[c(4, 2, 1)])

## ------------------------------------------------------------------------
R.version$version.string

## ---- eval=FALSE---------------------------------------------------------
#  library("reticulate")
#  py_install("python-igraph")
#  py_install("leidenalg")

## ------------------------------------------------------------------------
partition <- py$partition$membership + 1
table(partition)

## ------------------------------------------------------------------------
library("igraph")
library("reticulate")
library("RColorBrewer")
graph_object <- graph.famous("Zachary")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols, layout=layout_with_kk)

## ------------------------------------------------------------------------
partition <- py$partition$membership + 1
table(partition)

## ------------------------------------------------------------------------
graph_object <- graph.famous("Zachary")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols, layout=layout_with_kk)

## ------------------------------------------------------------------------
partition <- py$partition$membership + 1
table(partition)

## ------------------------------------------------------------------------
graph_object <- graph.famous("Zachary")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols, layout=layout_with_kk)

## ------------------------------------------------------------------------
bash_py_time <- as.numeric(readLines("bash_py_time"))

## ------------------------------------------------------------------------
leidenalg <- import("leidenalg", delay_load = TRUE)
ig <- import("igraph", delay_load = TRUE)
G = ig$Graph$Famous('Zachary')
G$summary()
partition = leidenalg$find_partition(G, leidenalg$ModularityVertexPartition)
partition$membership

## ------------------------------------------------------------------------
leidenalg <- import("leidenalg", delay_load = TRUE)
ig <- import("igraph", delay_load = TRUE)
G = ig$Graph$Famous('Zachary')
G$summary()
start <- Sys.time()
for(ii in 1:100){
  partition = leidenalg$find_partition(G, leidenalg$ModularityVertexPartition)
}
end <- Sys.time()
partition$membership
reticulate_time <- difftime(end, start)[[1]]
print(paste(c("leiden time:", reticulate_time, "seconds"), collapse = " "))

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("leiden")

## ------------------------------------------------------------------------
R.version.string

## ------------------------------------------------------------------------
library("igraph")
library("leiden")

## ------------------------------------------------------------------------
G <- graph.famous("Zachary")
summary(G)

## ------------------------------------------------------------------------
partition <- leiden(G, "ModularityVertexPartition")
partition

## ------------------------------------------------------------------------
table(partition)

## ------------------------------------------------------------------------
library("igraph")
library("reticulate")
library("RColorBrewer")
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(G, vertex.color = node.cols, layout=layout_with_kk)

## ------------------------------------------------------------------------
partition <- leiden(G, "CPMVertexPartition", resolution_parameter = 0.5)
partition

## ------------------------------------------------------------------------
table(partition)

## ------------------------------------------------------------------------
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(G, vertex.color = node.cols, layout=layout_with_kk)

## ------------------------------------------------------------------------
partition <- leiden(G, "RBConfigurationVertexPartition", resolution_parameter = 1.5)
partition

## ------------------------------------------------------------------------
table(partition)

## ------------------------------------------------------------------------
node.cols <- brewer.pal(max(c(3, partition)),"Pastel1")[partition]
plot(G, vertex.color = node.cols, layout=layout_with_kk)

## ---- cache=TRUE---------------------------------------------------------
G <- graph.famous('Zachary')
summary(G)
start <- Sys.time()
for(ii in 1:100){
  partition = leiden(G, "ModularityVertexPartition")
}
end <- Sys.time()
table(partition)
R_graph_time = difftime(end, start)[[1]]
print(paste(c("leiden time:", R_graph_time, "seconds"), collapse = " "))

## ---- cache=TRUE---------------------------------------------------------
G <- graph.famous('Zachary')
summary(G)

start <- Sys.time()
for(ii in 1:100){
  adj_mat <- as_adjacency_matrix(G, sparse = FALSE)
}
end <- Sys.time()
dim(adj_mat)
R_mat_cast_time = difftime(end, start)[[1]]
paste(print(c("cast time:", R_mat_cast_time, "seconds"), collapse = " "))

start <- Sys.time()
for(ii in 1:100){
  partition <- leiden(adj_mat, "ModularityVertexPartition")
}
end <- Sys.time()
table(partition)
R_mat_time = difftime(end, start)[[1]]
print(paste(c("leiden time:", R_mat_time, "seconds"), collapse = " "))

## ---- cache=TRUE---------------------------------------------------------
G <- graph.famous('Zachary')
summary(G)

start <- Sys.time()
for(ii in 1:100){
  adj_mat <- as_adjacency_matrix(G, sparse = TRUE)
}
end <- Sys.time()
class(adj_mat)
dim(adj_mat)
R_sparse_mat_cast_time = difftime(end, start)[[1]]
paste(print(c("cast time:", R_sparse_mat_cast_time, "seconds"), collapse = " "))

start <- Sys.time()
for(ii in 1:100){
  partition <- leiden(adj_mat, "ModularityVertexPartition")
}
end <- Sys.time()
table(partition)
R_sparse_mat_time = difftime(end, start)[[1]]
print(paste(c("leiden time:", R_mat_time, "seconds"), collapse = " "))

## ---- cache=TRUE---------------------------------------------------------
adjacency_matrix <- rbind(cbind(matrix(round(rbinom(1000000, 1, 0.008)), 1000, 1000),
                                matrix(round(rbinom(1000000, 1, 0.003)), 1000, 1000),
                                matrix(round(rbinom(1000000, 1, 0.001)), 1000, 1000)),
                          cbind(matrix(round(rbinom(1000000, 1, 0.003)), 1000, 1000),
                                matrix(round(rbinom(1000000, 1, 0.008)), 1000, 1000),
                                matrix(round(rbinom(0000000, 1, 0.002)), 1000, 1000)),
                          cbind(matrix(round(rbinom(1000000, 1, 0.003)), 1000, 1000),
                                matrix(round(rbinom(1000000, 1, 0.001)), 1000, 1000),
                                matrix(round(rbinom(1000000, 1, 0.009)), 1000, 1000)))
rownames(adjacency_matrix) <- 1:3000
colnames(adjacency_matrix) <- 1:3000
G <- graph_from_adjacency_matrix(adjacency_matrix)

start <- Sys.time()
for(ii in 1:10){
  adj_mat <- as_adjacency_matrix(G, sparse = FALSE)
}
end <- Sys.time()
class(adj_mat)
dim(adj_mat)
R_mat_large_cast_time = difftime(end, start)[[1]]
paste(print(c("cast time:", R_mat_large_cast_time, "seconds"), collapse = " "))

start <- Sys.time()
for(ii in 1:10){
  partition <- leiden(adj_mat, "ModularityVertexPartition")
}
end <- Sys.time()
table(partition)
R_mat_large_time = difftime(end, start)[[1]]
print(paste(c("leiden time:", R_mat_large_time, "seconds"), collapse = " "))

## ---- cache=TRUE---------------------------------------------------------
start <- Sys.time()
for(ii in 1:100){
  adj_mat <- as_adjacency_matrix(G, sparse = TRUE)
}
end <- Sys.time()
class(adj_mat)
dim(adj_mat)
R_mat_large_cast_time = difftime(end, start)[[1]]
paste(print(c("cast time:", R_mat_large_cast_time, "seconds"), collapse = " "))

start <- Sys.time()
for(ii in 1:10){
  partition <- leiden(adj_mat, "ModularityVertexPartition")
}
end <- Sys.time()
table(partition)
R_mat_large_time = difftime(end, start)[[1]]
print(paste(c("leiden time:", R_mat_large_time, "seconds"), collapse = " "))

## ------------------------------------------------------------------------
partition_type <- "RBConfigurationVertexPartition"
initial_membership <- NULL
weights <- NULL
node_sizes = NULL
resolution_parameter = 1

G <- graph.famous('Zachary')
summary(G)
time1 <- Sys.time()
object <- as.matrix(as_adjacency_matrix(G))
time2 <- Sys.time()
timing = difftime(time2, time1)[[1]]
print(paste(c("cast to adjacent:", timing, "seconds"), collapse = " "))

#run matrix method
leidenalg <- import("leidenalg", delay_load = TRUE)
ig <- import("igraph", delay_load = TRUE)

#convert matrix input (corrects for sparse matrix input)
if(is.matrix(object) || is(adj_mat_sparse, "Matrix")){
  adj_mat <- object
} else{
  adj_mat <- as.matrix(object)
}

#compute weights if non-binary adjacency matrix given
is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
if (is.null(weights) && !is_pure_adj) {
  #assign weights to edges (without dependancy on igraph)
  t_mat <- t(adj_mat)
  weights <- t_mat[t_mat!=0]
  #remove zeroes from rows of matrix and return vector of length edges
}

time3 <- Sys.time()
##convert to python numpy.ndarray, then a list
adj_mat_py <- r_to_py(adj_mat)
adj_mat_py <- adj_mat_py$tolist()
time4 <- Sys.time()
timing = difftime(time4, time3)[[1]]
print(paste(c("pass to python matrix:", timing, "seconds"), collapse = " "))


#convert graph structure to a Python compatible object
GraphClass <- if (!is.null(weights) && !is_pure_adj){
  ig$Graph$Weighted_Adjacency
} else {
  ig$Graph$Adjacency
}
time5 <- Sys.time()
snn_graph <- GraphClass(adj_mat_py)
time6 <- Sys.time()
timing = difftime(time6, time5)[[1]]
reticulate_create_time = difftime(time6, time5)[[1]]
print(paste(c("generate graph in python:", timing, "seconds"), collapse = " "))


# test performance for computing matrix to graph in R
# other option is to passing snn_graph to Python

time7 <- Sys.time()
#compute partitions
source("../R/find_partition.R")

partition <- find_partition(snn_graph, partition_type = partition_type,
                            initial_membership = initial_membership ,
                            weights = weights,
                            node_sizes = node_sizes,
                            resolution_parameter = resolution_parameter
)
time8 <- Sys.time()
timing = difftime(time8, time7)[[1]]
print(paste(c("compute partitions:", timing, "seconds"), collapse = " "))
timing = difftime(time8, time1)[[1]]
print(paste(c("total:", timing, "seconds"), collapse = " "))
partition

## ------------------------------------------------------------------------
partition_type <- "RBConfigurationVertexPartition"
initial_membership <- NULL
weights <- NULL
node_sizes = NULL
resolution_parameter = 1

G <- graph.famous('Zachary')
summary(G)
time1 <- Sys.time()
object <- as.matrix(as_adjacency_matrix(G))
time2 <- Sys.time()
timing = difftime(time2, time1)[[1]]
print(paste(c("cast to adjacent:", timing, "seconds"), collapse = " "))

#run matrix method
leidenalg <- import("leidenalg", delay_load = TRUE)
ig <- import("igraph", delay_load = TRUE)

time3 <- Sys.time()
##convert to python numpy.ndarray, then a list
object <- graph_from_adjacency_matrix(adj_mat)
time4 <- Sys.time()
timing = difftime(time4, time3)[[1]]
print(paste(c("generate graph in R:", timing, "seconds"), collapse = " "))

#convert graph structure to a Python compatible object
time5 <- Sys.time()
##convert to list for python input
    if(!is.named(object)){
        vertices <- as.list(as.character(V(object)))
    } else {
        vertices <- as.list(names(V(object)))
    }

    edges <- as_edgelist(object)
    dim(edges)
    edgelist <- list(rep(NA, nrow(edges)))
    for(ii in 1:nrow(edges)){
        edgelist[[ii]] <- as.character(edges[ii,])
    }

    snn_graph <- ig$Graph()
    snn_graph$add_vertices(r_to_py(vertices))
    snn_graph$add_edges(r_to_py(edgelist))
time6 <- Sys.time()
timing = difftime(time6, time5)[[1]]
print(paste(c("pass to python graph:", timing, "seconds"), collapse = " "))



# test performance for computing matrix to graph in R
# other option is to passing snn_graph to Python

time7 <- Sys.time()
#compute partitions
partition <- find_partition(snn_graph, partition_type = partition_type,
                            initial_membership = initial_membership ,
                            weights = weights,
                            node_sizes = node_sizes,
                            resolution_parameter = resolution_parameter
)
time8 <- Sys.time()
timing = difftime(time8, time7)[[1]]
print(paste(c("compute partitions:", timing, "seconds"), collapse = " "))
timing = difftime(time8, time1)[[1]]
print(paste(c("total:", timing, "seconds"), collapse = " "))
partition

## ------------------------------------------------------------------------
partition_type <- "RBConfigurationVertexPartition"
initial_membership <- NULL
weights <- NULL
node_sizes = NULL
resolution_parameter = 1

G <- graph.famous('Zachary')
summary(G)
time1 <- Sys.time()
object <- as.matrix(as_adjacency_matrix(G))
time2 <- Sys.time()
timing = difftime(time2, time1)[[1]]
print(paste(c("cast to adjacent:", timing, "seconds"), collapse = " "))

time3 <- Sys.time()
##convert to python numpy.ndarray, then a list
object <- graph_from_adjacency_matrix(adj_mat)
time4 <- Sys.time()
timing = difftime(time4, time3)[[1]]
R_graph_create_time = difftime(time4, time3)[[1]]
print(paste(c("generate graph in R:", timing, "seconds"), collapse = " "))


#convert graph structure to a Python compatible object
time5 <- Sys.time()
##convert to list for python input
   snn_graph <- object
time6 <- Sys.time()
timing = difftime(time6, time5)[[1]]
print(paste(c("pass to R graph:", timing, "seconds"), collapse = " "))



# test performance for computing matrix to graph in R
# other option is to passing snn_graph to Python

time7 <- Sys.time()
#compute partitions
partition <- leiden(snn_graph, partition_type = partition_type,
                            initial_membership = initial_membership ,
                            weights = weights,
                            node_sizes = node_sizes,
                            resolution_parameter = resolution_parameter
)
time8 <- Sys.time()
timing = difftime(time8, time7)[[1]]
print(paste(c("compute partitions:", timing, "seconds"), collapse = " "))
timing = difftime(time8, time1)[[1]]
print(paste(c("total:", timing, "seconds"), collapse = " "))
partition

## ---- fig.align = 'center', fig.height = 3, fig.width = 6, fig.keep = 'last'----
barplot(c(bash_py_time, py$py_time, reticulate_time, R_graph_time, 
          R_mat_time, R_sparse_mat_time), 
        names = c("Python (shell)", "Python (Rmd)", "Reticulate", "R igraph",
                  "R matrix","R dgCMatrix"), 
        col = brewer.pal(9,"Pastel1"), las = 2, srt = 45,
        ylab = "time (seconds)", main = "benchmarking 100 computations")
abline(h=0)

## ---- fig.align = 'center', fig.height = 3, fig.width = 6, fig.keep = 'last'----
barplot(c(bash_py_time, py$py_time, reticulate_time, R_graph_time, R_mat_time+R_mat_cast_time, 
          R_sparse_mat_time+R_sparse_mat_cast_time), 
        names = c("Python (shell)", "Python (Rmd)", "Reticulate", "R igraph",
                  "R matrix","R dgCMatrix"), 
        col = "grey80", las = 2, srt = 45,
        ylab = "time (seconds)", main = "benchmarking 100 computations")
barplot(c(bash_py_time, py$py_time, reticulate_time, R_graph_time,
          R_mat_time,  R_sparse_mat_time), 
        names = c("Python (shell)", "Python (Rmd)", "Reticulate", "R igraph",
                  "R matrix","R dgCMatrix"), 
        col = brewer.pal(9,"Pastel1"), las = 2, srt = 45,
        ylab = "time (seconds)", main = "benchmarking 100 computations", add = TRUE)
abline(h=0)

## ---- fig.align = 'center', fig.height = 3, fig.width = 6, fig.keep = 'last'----
R_graph_create_time = difftime(time4, time3)[[1]]
barplot(c(bash_py_time, py$py_time+reticulate_create_time*100, reticulate_time+reticulate_create_time*100, R_graph_time+R_graph_create_time*100, R_mat_time, 
          R_sparse_mat_time), 
        names = c("Python (shell)", "Python (Rmd)", "Reticulate", "R igraph",
                  "R matrix","R dgCMatrix"), 
        col = "grey80", las = 2, srt = 45,
        ylab = "time (seconds)", main = "benchmarking 100 computations")
barplot(c(bash_py_time, py$py_time, reticulate_time, R_graph_time,
          R_mat_time,  R_sparse_mat_time), 
        names = c("Python (shell)", "Python (Rmd)", "Reticulate", "R igraph",
                  "R matrix","R dgCMatrix"), 
        col = brewer.pal(9,"Pastel1"), las = 2, srt = 45,
        ylab = "time (seconds)", main = "benchmarking 100 computations", add = TRUE)
abline(h=0)

