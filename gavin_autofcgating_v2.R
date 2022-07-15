library(gmodels)
library(Ckmeans.1d.dp)
library(tidyverse)
library(mclust)
library(fpc)
library(plotly)
library(ggplot2)
library(grDevices)
library(RColorBrewer)
library(foreach)

z_normalization <- function(){
  # https://jmotif.github.io/sax-vsm_site/morea/algorithm/znorm.html
  return(function(data) {
    if (length(dim(data)) == 0){
      data <- as.data.frame(data)
    }
    apply(data, 2, function(ts){
      ts.mean <- mean(ts)
      ts.dev <- sd(ts)
      (ts - ts.mean)/ts.dev
    })
  })
}

outlier_resistant_minmax <- function(pos=c()){
  pos <- c(0, pos, 1, 1)
  return(function(data) {
    if (length(dim(data)) == 0){
      data <- as.data.frame(data)
    }
    apply(data, 2, function(ts){
      qtls <- quantile(ts, pos)
      for (row in seq_along(ts)){
        for (i in seq_len(length(pos)-1)){
          if (qtls[i+1] >= ts[row]){
            if (qtls[i] == qtls[i+1]) ts[row] <- pos[i]
            else ts[row] <- (ts[row]-qtls[i])/(qtls[i+1]-qtls[i])*(pos[i+1]-pos[i])+pos[i]
            break();
          }
        }
      }
      ts
    })
  })
}

bic <- function(data, range){
  range[which.max(mclustBIC(data, G=range, modelNames="V", verbose=F))]
}

silhouette <- function(data, range){
  pamk(data, krange=range)$nc
}

gmm_clustering <- function(g=1:9, f=bic){
  return(function(data) Mclust(data, G=f(data, g), modelNames = "V", verbose = F)$classification)
}

default_kmeans_clustering <- function(k=1:9){
  return(function(data) suppressWarnings(Ckmeans.1d.dp(data, k = k))$cluster)
}

kmeans_clustering <- function(k=1:9, f=silhouette){
  return(function(data) suppressWarnings(Ckmeans.1d.dp(data, k = f(data, k)))$cluster)
}

make_label <- function(){
  runif(1, 1000000000000, 9999999999999) %>% round %>% as.character
}

jitter_height <- function(dendro, ...){
  attr(dendro, "height") <- jitter(attr(dendro, "height"), ...)
  if (!is.null(attr(dendro, "leaf")) && attr(dendro, "leaf")) return(dendro)
  res <- lapply(dendro, function(dend) jitter_height(dend, ...))
  attributes(res) <- attributes(dendro)
  return(res)
}

better_hclust <- function(
  data,
  normalization=outlier_resistant_minmax(),
  clustering=default_kmeans_clustering(),
  max_depth=-1,
  do_parallel=T
){
  if (max_depth == 0) {
    res <- list(0)
    attributes(res) <- list(
      members = 1,
      midpoint = 0,
      height = 0,
      label = make_label(),
      leaf = TRUE,
      data = rownames(data)
    )
    return(res)
  }
  og <- data
  if (nrow(data) <= 2){
    res <- list(0)
    attributes(res) <- list(
      members = 1,
      midpoint = 0,
      height = 0,
      label = make_label(),
      leaf = TRUE,
      data = rownames(og)
    )
    return(res)
  }
  data <- data %>% dplyr::select(where(is.numeric)) %>% normalization()
  pca <- fast.prcomp(data)$x
  #config <- umap::umap.defaults
  #config$n_neighbors <- min(nrow(data), config$n_neighbors)
  #config$n_components <- 1
  #pca <- umap::umap(data, config)$layout
  for (pc in seq_len(ncol(pca))){
    clusters <- clustering(pca[,pc])
    if (max(clusters) > 1){
      res <- foreach(clust=seq_len(max(clusters))) %do% better_hclust(
        data = og[clusters == clust,],
        normalization = normalization,
        clustering = clustering,
        max_depth = max_depth-1,
        do_parallel = F
      )
      mems <- res %>% lapply(function(x) attr(x, "members")) %>% unlist() %>% sum()
      attributes(res) <- list(
        members = mems,
        midpoint = mems/ 2,
        height = res %>% lapply(function(x) attr(x, "height")) %>% unlist() %>% max() + 1,
        label = make_label(),
        data = res %>% lapply(function(x) attr(x, "data")) %>% unlist()
      )
      class(res) <- "dendrogram"
      return(res)
    }
  }
  res <- list(0)
  attributes(res) <- list(
    members = 1,
    midpoint = 0.5,
    height = 0,
    label = make_label(),
    leaf = TRUE,
    data = rownames(og)
  )
  return(res)
}

testClust <- function(clust, gh, parent, popmem, paths, compare="parent", significant_pval = 0.01){
  membershipClust <- as.matrix(popmem[as.numeric(attributes(clust)$data),])
  membershipParent <- as.matrix(popmem[as.numeric(attributes(parent)$data),])
  if (length(attributes(clust)$data) != 1){
    poscounts <- colCounts(membershipClust)
  } else {
    poscounts <- rowCounts(membershipClust)
  }
  poscountsParent <- colCounts(membershipParent)
  names(poscounts) <- paths
  names(poscountsParent) <- paths
  rval <- lapply(paths, function(path) {
    ppath <- ifelse(path == "root", "root", gh_pop_get_parent(gh, path))
    if (!poscounts[[ppath]] || !poscountsParent[[ppath]]){
      return("No data available for parent")
    }
    binom.test(poscounts[[path]], n = poscounts[[ppath]], p = poscountsParent[[path]]/poscountsParent[[ppath]], alternative = "greater")
  })
  names(rval) <- paths
  attributes(clust)$tests <- rval
  pval <- lapply(rval, function(test) ifelse(is.character(test), 1, test$p.value))
  names(pval) <- paths
  attributes(clust)$pval <- pval
  attributes(clust)$classification <- ifelse(length(which(pval < significant_pval))!=0, paths[[which.min(unlist(pval))]], "root")
  attributes(clust)$hovertext <- do.call(paste0, as.list(do.call(c, lapply(names(pval)[pval < significant_pval], function(path) list(path, " ", pval[[path]], "\n")))))
  if (attributes(clust)$hovertext == character(0)) {
    attributes(clust)$hovertext <- "No significant categories"
  }
  return(clust)
}

getTests <- function(result, gh, root=NULL, compare="parent", significant_pval=0.01){
  library(matrixStats)
  library(flowWorkspace)
  paths <- rev(gh_get_pop_paths(gh))
  popmem <- foreach(pop=paths, .combine='cbind') %do% gh_pop_get_indices(gh, pop)
  colnames(popmem) <- paths
  if (is.null(root)){
    result <- testClust(result, gh, result, popmem, paths, compare=compare)
    root <- result
  }
  for (clustno in seq_len(length(result))){
    if (compare == "root") { # "not better than randomly pulling items from the entire dataset" as null hypothesis
      clust <- testClust(result[[clustno]], gh, root, popmem, paths, compare=compare, significant_pval = significant_pval)
    } else if (compare == "parent") { # "not better than randomly pulling items from the parent cluster" as null hypothesis
      clust <- testClust(result[[clustno]], gh, result, popmem, paths, compare=compare, significant_pval = significant_pval)
    }
    
    if (is.null(attributes(clust)$leaf) || !attributes(clust)$leaf){
      result[[clustno]] <- getTests(clust, gh, root, compare = compare, significant_pval = significant_pval)
    } else {
      result[[clustno]] <- clust
    }
  }
  return(result)
}



unlist_dendro <- function(d){
  if (!is.null(attributes(d)$leaf) && attributes(d)$leaf) return(list(d))
  return(c(list(d), do.call(c, lapply(d, unlist_dendro))))
}

plot_dendro <- function(d, set = "A", xmin = -50, height = 500, width = 500, ...) {
  print("a1")
  # get x/y locations of every node in the tree
  allXY <- get_xy(d)
  allXY <- allXY[order(allXY[,2]),]
  print("a2")
  # get non-zero heights so we can split on them and find the relevant labels
  non0 <- allXY[["y"]][allXY[["y"]] > 0]
  print("a3")
  # splitting on the minimum height would generate all terminal nodes anyway
  split <- non0[min(non0) < non0]
  print("a4")
  # label is a list-column since non-zero heights have multiple labels
  # for now, we just have access to terminal node labels
  print("a5")
  allXY$label <- vector("list", nrow(allXY))
  allXY$classification <- vector("list", nrow(allXY))
  allXY$hovertext <- vector("list", nrow(allXY))
  print("a6")
  print("a7")
  print("a8")
  
  # collect all the *unique* non-trivial nodes
  nodes <- unlist_dendro(d)
  nodes <- nodes[order(sapply(nodes, function(node) attributes(node)$height))]
  print("a10")
  print("a11")
  print("a12")
  # NOTE: this won't support nodes that have the same height 
  # but that isn't possible, right?
  for (i in seq_along(nodes)) {
    node <- nodes[[i]]
    if (is.null(attributes(node)$classification)){
      allXY$classification[[i]] <- "root"
      allXY$hovertext[[i]] <- "leaf"
    } else {
      allXY$classification[[i]] <- attributes(node)$classification
      allXY$hovertext[[i]] <-  attributes(node)$hovertext
    }
  }
  print("a13")
  tidy_segments <- dendextend::as.ggdend(d)$segments
  print("a14")
  allTXT <- allXY[allXY$y == 0, ]
  print("a15")
  blank_axis <- list(
    title = "",
    showticklabels = FALSE,
    zeroline = FALSE
  )
  print("a16")
  allXY$members <- sapply(allXY$label, length)
  print("a17")
  allTXT$label <- as.character(allTXT$label)
  print("a18")
  allXY %>% 
    plot_ly(x = ~y, y = ~x, color=I("black"), hoverinfo = "none") %>%
    add_segments(
      data = tidy_segments, xend = ~yend, yend = ~xend, showlegend = FALSE
    ) %>%
    add_markers(
      data = allXY, x=~y, y=~x, key = ~label,
      text = ~paste(hovertext), hoverinfo = "text", color = ~unlist(classification), colors = c("black", brewer.pal(8, "Dark2"))
    ) %>%
    layout(
      dragmode = "select", 
      xaxis = c(blank_axis, list(autorange = TRUE)),
      yaxis = c(blank_axis, list(autorange = TRUE))
    )
}

get_xy <- function(node) {
  m <- dendextend::get_nodes_xy(node)
  colnames(m) <- c("x", "y")
  tibble::as_tibble(m)
}
# stop()
if (getOption('run.main', default=interactive())) {
  # library(doParallel)
  # library(parallel)
  # no_cores <- detectCores(logical = TRUE)
  # cl <- makeCluster(no_cores-1)
  # registerDoParallel(cl)
  
  # remotes::install_github("RGLab/RProtoBufLib", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/cytolib", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/flowCore", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/flowWorkspace", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/flowWorkspaceData", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/flowStats", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/ncdfFlow", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/CytoML", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/openCyto", upgrade="never", quiet=F)
  # remotes::install_github("RGLab/ggcyto", upgrade="never", quiet=F)
  library(flowCore)
  library(CytoML)
  library(flowWorkspace)
  library(openCyto)
  
  dataDir <- system.file("extdata",package="flowWorkspaceData")
  #load raw FCS
  fs <- load_cytoset_from_fcs(file.path(dataDir,"CytoTrol_CytoTrol_1.fcs"))
  
  #load the original template for tcell panel
  tbl <- data.table::fread("./tcell.csv")
  #write the new template to disc
  gtFile <- tempfile()
  write.csv(tbl, file = gtFile)
  ##reload the new template
  gt <- gatingTemplate(gtFile)
  
  # do 'manual' gating
  gs <- GatingSet(fs)
  comp <- spillover(fs[[1]])[["SPILL"]]
  chnls <- colnames(comp)
  comp <- compensation(comp)
  gs <- compensate(gs, comp)
  trans <- flowjo_biexp_trans()
  trans <- transformerList(chnls, trans)
  gs <- transform(gs, trans)
  gt_gating(gt, gs)
  gt_toggle_helpergates(gt, gs)
  # do unsupervised gating
  res <- better_hclust(as.data.frame(exprs(fs[[1]])), max_depth=-1)
  View(res)
  # perform statistical tests
  res <- getTests(res, gs[[1]])
  View(res)
  # generate ptree
  #p$data <- cbind(newdata, p$data)
  plot_dendro(res)
}
