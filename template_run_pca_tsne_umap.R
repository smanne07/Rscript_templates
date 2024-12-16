library(umap)
library(Rtsne)
library(tidyverse)
library(bench)
library(ggplot2)
library(ggthemes)
custom.settings = umap.defaults
custom.settings$n_neighbors = 2
# stuff to compare algorithms -------------------------------------------------
embed <- function(labels, d) {
  times <- mark(
    um <- umap(d,custom.settings)$layout,
    ts <- Rtsne(d,perplexity = 2)$Y,
    ts_no_pca <- Rtsne(d, pca = FALSE,perplexity = 2)$Y,
    check = FALSE)

  pca <- prcomp(d)$x[,1:2]

  times$expression <- c("UMAP", "PCA + t-SNE", "t-SNE")

  combo <- function(embedding, name) {
    colnames(embedding) <- c("V1", "V2")
    embedding %>%
      as.data.frame() %>%
      mutate(Algorithm = name, Class = labels)
  }

  list(times = times,
       results = bind_rows(
         combo(pca, "PCA"),
         combo(um, "UMAP"),
         combo(ts, "PCA + t-SNE"),
         combo(ts_no_pca, "t-SNE")))
}

plot_embeddings <- function(embeddings, dataset) {
  ggplot(embeddings, aes(V1, V2, color = Class)) +
    geom_point(size=3) + facet_wrap(~ Algorithm, scales = "free") +
    ggtitle(dataset)+theme_par(base_size = 12)
}

# iris -----------------------------------------------------------------------
d <- iris
d <- d[!duplicated(d), ]
with_labels <- d
d <- as.matrix(d[ , 1:4])

iris_result <- embed(labels=with_labels$Species, d)

plot_embeddings(iris_result$results, "iris")
ggsave("img/multiple_algorithms_iris.png", width = 6, height = 5, dpi = 300)


##Running on normalized RNAseq data
#mergdata object as input from limma analysis

d <- as.matrix(mergData[ ,-c(1:8)])
labels=mergData[,5]
embed_result <- embed(labels=labels, d)

plot_embeddings(embed_result$results, "LN")
ggsave("output/multiple_algorithms_dimred_LN.png", width = 10, height = 10, dpi = 300)


plot_embeddings(embed_result$results, "CNS")
ggsave("output/multiple_algorithms_dimred_CNS.png", width = 10, height = 10, dpi = 300)
