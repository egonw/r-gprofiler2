find_all_anchestors = function(best, parents_df){
  anchs = c()
  queue = c(best)
  while (length(queue) > 0){
    term = queue[1]
    anchs = unique(append(anchs, parents_df$parent[parents_df$child == term]))
    queue = unique(append(queue, parents_df$parent[parents_df$child == term]))
    if (length(queue) == 1){
      queue = c()
    } else {
      queue = queue[2:length(queue)]
    }
  }
  return(anchs)
}

find_all_children = function(best, parents_df){
  children = c()
  queue = c(best)
  while (length(queue) > 0){
    term = queue[1]
    children = unique(append(children, as.character(parents_df$child)[parents_df$parent == term]))
    queue = unique(append(queue, as.character(parents_df$child)[parents_df$parent == term]))
    if (length(queue) == 1){
      queue = c()
    } else {
      queue = queue[2:length(queue)]
    }
  }
  return(children)
}

# plotting function
graphplot = function(subres, keep = "", filename = "test.png", dir = "/Users/Liis/Desktop/hier_filtering"){
  if(dim(subres)[1] > 0){
    url = "http://docker.cs.ut.ee:60041/gprofiler_beta/api/util/draw"
    body <- jsonlite::toJSON((
      list(
        nodes = as.list(structure(subres$p_value, names = as.character(subres$term_id))),
        colored = paste0(keep, collapse = " ")
      )
    ),
    auto_unbox = TRUE, digits = 10)

    oldw <- getOption("warn")
    options(warn = -1)
    h2 = RCurl::getCurlHandle() # Get info about the request

    # Request
    img = RCurl::getBinaryURL(
      url = url,
      postfields = body,
      #httpheader = "image/png",
      customrequest = 'POST',
      verbose = FALSE,
      ssl.verifypeer = FALSE,
      curl = h2
    )
    rescode = RCurl::getCurlInfo(h2)[["response.code"]]

    if (rescode != 200) {
      stop("Bad request, response code ", rescode)
    }

    extension = strsplit(basename(filename), split="\\.")[[1]][2]
    if (!grepl(extension, "png")){
      filename = paste0(filename, ".png")
    }
    save_filename = file.path(dir, filename)
    f = file(save_filename, "wb")
    writeBin(object = as.raw(as.hexmode(as.character(img))), con = f)
    close(f)
  }
  else{
    stop("No data to use")
  }
  return(save_filename)
}

hybrid_filtering = function(gostres, group_nr, plot = TRUE, dir = "/Users/Liis/Desktop/PSB_article"){
  res = gostres$result
  subres = res[res$group_id==group_nr,]
  subres$logpval = -log10(subres$p_value)

  # Order the terms in a group by the p_value
  subres = subres[order(subres$logpval, decreasing = T), ]

  # Find the parent-child relations
  parents <- lapply(subres$parents, function(x) unlist(strsplit(x, ", ")))
  names(parents) <- subres$term_id
  parents <- stack(parents)
  names(parents) <- c("parent", "child")
  parents$child <- as.character(parents$child)

  # Find the significance threshold
  pval_thr = gostres$meta$result_metadata[[subres$source[1]]]$threshold

  # Vector for keeping the important terms
  keep = c()
  # Vector for keeping the marked genes
  markedgenes = c()

  while (dim(subres)[1] > 0){
    # Take the first term
    best = subres$term_id[1]
    # Find the p_value
    new_intersection = setdiff(strsplit(subres[subres$term_id==best,]$intersection, ",")[[1]], markedgenes)
    pval = phyper(length(new_intersection) - 1, subres[subres$term_id==best,]$term_size, subres[subres$term_id==best,]$effective_domain_size-subres[subres$term_id==best,]$term_size, subres[subres$term_id==best,]$query_size, lower.tail=F)

    # If term is significant then keep
    if (pval <= pval_thr){
      keep = append(keep, best)
      markedgenes = unique(append(markedgenes, strsplit(subres[subres$term_id==best,]$intersection, ",")[[1]]))
    }

    # Find the ancestors and children
    best_parents = find_all_anchestors(best, parents)
    best_children = find_all_children(best, parents)
    # Exclude them from the search
    subres = subres[!(subres$term_id %in% c(best, best_parents, best_children)),]
  }

  # Save graphs
  if (plot){
    graphplot(res[res$group_id==group_nr,], keep, paste0("filt_group_", group_nr, ".png"), dir = dir)
  }

  return(keep)
}

# Test the method
gostres <- my_gost(c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), evcodes = T)
res = gostres$result
filtered_res = list()
for (gr in unique(res$group_id)){
  print(gr)
  if(!res[res$group_id == gr,][["source"]][1] %in% c("GO:BP", "GO:CC", "GO:MF")){
    next
  }
  filtered_res[[gr]] = hybrid_filtering(gostres, gr, plot = TRUE, dir = "/Users/Liis/Desktop/PSB_article/gp_example_query")
}

# Test 2
gostres <- my_gost(c("GO:0005005", "GO:0035184"), evcodes = T)
res = gostres$result
filtered_res = list()
for (gr in unique(res$group_id)){
  print(gr)
  if(!res[res$group_id == gr,][["source"]][1] %in% c("GO:BP", "GO:CC", "GO:MF")){
    next
  }
  filtered_res[[gr]] = hybrid_filtering(gostres, gr, plot = F, dir = "/Users/Liis/Desktop/PSB_article/gp_example_query2")
}

graphplot(res[res$source %in% c("GO:MF"),], unlist(filtered_res[1:3]), "full_MF.png", dir = "/Users/Liis/Desktop/PSB_article")


# Heatmap of remaining terms in the row and query genes
subres = res[res$term_id %in% unlist(filtered_res),]
intersections = lapply(subres$intersection, function(x) strsplit(x, ",")[[1]])
names(intersections) = subres$term_id
df = stack(intersections)
library(reshape2)
df_wide = data.frame(dcast(df, ind ~ values))
#row.names(df_wide) = df_wide$ind
terms = df_wide$ind
df_wide$ind = NULL
df_wide[is.na(df_wide)] = 0
df_wide[df_wide != 0] = 1
query_genes = gostres$meta$genes_metadata$query$query_1$ensgs
df_wide = df_wide[,query_genes]
df_wide = sapply(df_wide, as.numeric)
row.names(df_wide) = terms
library(pheatmap)
pheatmap(df_wide, cluster_rows = T, cluster_cols = F)

## Plot Revigo terms

terms = c("GO:0003824", "GO:0005003", "GO:0008144", "GO:0060089", "GO:0016740", "GO:0046875", "GO:0005524",
          "GO:0097367", "GO:0036094", "GO:0043167", "GO:1901363", "GO:0097159", "GO:1901265", "GO:0035184",
          "GO:0035173", "GO:0016773", "GO:0005102", "GO:0043168", "GO:0008046", "GO:0016772",
          "GO:0004697", "GO:0004713", "GO:0019199", "GO:0038023", "GO:0035639", "GO:0017076")
highlight = c("GO:0005003", "GO:0003824", "GO:0008144", "GO:0060089", "GO:0016740", "GO:0046875", "GO:0005524", "GO:0097367")
graphplot(res[res$term_id %in% terms,], highlight, "revigo.png", dir = "/Users/Liis/Desktop/PSB_article")
graphplot(res[res$source %in% c("GO:MF"),], terms, "revigo2.png", dir = "/Users/Liis/Desktop/PSB_article")

## Example gene query
