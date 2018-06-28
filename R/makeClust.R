
## makeClust

makeClust <- function(data, clust.method, dist.method, p) {
  assert_choice(clust.method, choices = c("complete", "single", "average",
    "median", "ward.D", "ward.D2", "mcquitty", "centroid"))
  assert_choice(dist.method, choices = c("euclidean", "maximum", "manhattan",
    "canberra", "binary", "minkowski"))
  assert_numeric(p, lower = 1, len = 1, any.missing = FALSE)
  inds <- hclust(dist(data, method = dist.method, p = p), method = clust.method)$order
  return(inds)
}
