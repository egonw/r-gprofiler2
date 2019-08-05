gp_globals = new.env()

gp_globals$version =
  tryCatch(
    utils::packageVersion("gprofiler2"),
    error = function(e) { return("unknown_version") }
  );

# Set SSL version to TLSv1_2 with fallback to TLSv1_1
# CURL_SSLVERSION_SSLv3 is not used due to the SSLv3 vulnerability <https://access.redhat.com/articles/1232123>
# CURL_SSLVERSION_TLSv1_3 is not widespread enough to have a built-in LibreSSL support yet.
# (curl's authors may decide to change it at some point, so links to the source are provided.)
gp_globals$CURL_SSLVERSION_TLSv1_1 <- 5L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1925>
gp_globals$CURL_SSLVERSION_TLSv1_2 <- 6L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1926>

gp_globals$rcurl_opts =
  RCurl::curlOptions(useragent = paste("gprofiler2/", gp_globals$version, sep=""), sslversion = gp_globals$CURL_SSLVERSION_TLSv1_2)
gp_globals$base_url = "http://biit.cs.ut.ee/gprofiler_beta"


#' Gene list functional enrichment.
#'
#' Interface to the g:Profiler tool g:GOSt for functional enrichments analysis of gene lists.
#' In case the input 'query' is a list of gene vectors, results for multiple queries will be returned in the same data frame with column 'query' indicating the corresponding query name.
#' If 'multi_query' is selected, the result is a data frame for comparing multiple input lists,
#' just as in the web tool.
#'
#' @param query vector, or a (named) list of vectors for multiple queries, that can consist of mixed types of gene IDs (proteins, transcripts, microarray IDs, etc), SNP IDs, chromosomal intervals or term IDs.
#' @param organism organism name. Organism names are constructed by concatenating the first letter of the name and the
#' family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#' @param ordered_query in case input gene lists are ranked this option may be
#'  used to get GSEA style p-values.
#' @param multi_query in case of multiple gene lists, returns comparison table of these lists.
#' If enabled, the result data frame has columns named 'p_values', 'query_sizes', 'intersection_sizes' with vectors showing values in the order of input queries.
#' To get the results in a long format set 'multi_query' to FALSE and just input query list of multiple gene vectors.
#' @param significant whether all or only statistically significant results should
#'  be returned.
#' @param exclude_iea exclude GO electronic annotations (IEA).
#' @param measure_underrepresentation measure underrepresentation.
#' @param evcodes include evidence codes to the results. Note
#'  that this can decrease performance and make the query slower.
#'  In addition, a column 'intersection' is created that contains the gene id-s that intersect between the query and term.
#'  This parameter does not work if 'multi_query' is set to TRUE.
#' @param user_threshold custom p-value threshold, results with a larger p-value are
#'  excluded.
#' @param correction_method the algorithm used for multiple testing correction, one of "gSCS" (synonyms: "analytical", "g_SCS"), "fdr" (synonyms: "false_discovery_rate"), "bonferroni".
#' @param domain_scope how to define statistical domain, one of "annotated", "known" or "custom".
#' @param custom_bg vector of gene names to use as a statistical background. If given, the domain_scope is set to 'custom'.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param sources a vector of data sources to use. Currently, these include
#'  GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF,
#'  MIRNA, CORUM, HP, HPA, WP. Please see the g:GOSt web tool for the comprehensive
#'  list and details on incorporated data sources.
#' @return A named list where 'result' contains data.frame with the enrichment analysis results and 'meta' contains metadata needed for Manhattan plot. If the input
#'  consisted of several lists the corresponding list is indicated with a variable
#'  'query'.
#'  When requesting a 'multi_query', either TRUE or FALSE, the columns of the resulting data frame differ.
#'  If 'evcodes' is set, the return value includes columns 'evidence_codes' and 'intersection'.
#'  The latter conveys info about the intersecting genes between the corresponding query and term.
#' @author  Liis Kolberg <liis.kolberg@@ut.ee>, Uku Raudvere <uku.raudvere@@ut.ee>
#' @examples
#' gostres <- gost(c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"))
#'
#' @export
my_gost <- function(query,
                    organism = "hsapiens",
                    ordered_query = FALSE,
                    multi_query = FALSE,
                    significant = TRUE,
                    exclude_iea = TRUE,
                    measure_underrepresentation = FALSE,
                    evcodes = FALSE,
                    user_threshold = 0.05,
                    correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS", "analytical"),
                    domain_scope = c("annotated", "known", "custom"),
                    custom_bg = NULL,
                    numeric_ns  = "",
                    sources = NULL
) {

  url = paste0(file.path(gp_globals$base_url, "api", "gost", "profile"), "/")

  # Query

  if (is.null(query)) {
    stop("Missing query")
  } else if (is.list(query)) {
    if (is.data.frame(query)){
      stop("Query can't be a data.frame. Please use a vector or list of identifiers.")
    }
    # Multiple queries
    qnames = names(query)
    if (is.null(qnames)) {
      qnames = paste("query", seq(1, length(query)), sep = "_")
      names(query) = qnames
    }
    query = lapply(query, function(x) x[!is.na(x)])
  }
  else{
    query = query[!is.na(query)]
  }

  # Parameters

  ## evaluate choices
  correction_method <- match.arg(correction_method)

  if (startsWith(organism, "gp__")){
    message("Detected custom GMT source request")
    if (!is.null(sources)){
      message("Sources selection is not available for custom GMT requests. All sources in the GMT upload will be used.")
      sources <- NULL
    }
  }

  if (multi_query & evcodes){
    message("Evidence codes are not supported with multi_query option and will not be included in the results.\nYou can get evidence codes and intersection genes by setting multi_query = FALSE while keeping the input query as a list of multiple gene vectors.")
  }

  if (!is.null(custom_bg)){
    if (!is.vector(custom_bg)){
      stop("custom_bg must be a vector")
    }
    message("Detected custom background input, domain scope is set to 'custom'")
    domain_scope <- "custom"
    t <- ifelse(length(custom_bg) == 1, custom_bg <- jsonlite::unbox(custom_bg), custom_bg <- custom_bg)
  }

  domain_scope <- match.arg(domain_scope)

  body <- jsonlite::toJSON((
    list(
      organism = jsonlite::unbox(organism),
      query = query,
      sources = sources,
      user_threshold = jsonlite::unbox(user_threshold),
      all_results = jsonlite::unbox(!significant),
      ordered = jsonlite::unbox(ordered_query),
      no_evidences = jsonlite::unbox(!evcodes),
      combined = jsonlite::unbox(multi_query),
      measure_underrepresentation = jsonlite::unbox(measure_underrepresentation),
      no_iea = jsonlite::unbox(!exclude_iea),
      domain_scope = jsonlite::unbox(domain_scope),
      numeric_ns = jsonlite::unbox(numeric_ns),
      significance_threshold_method = jsonlite::unbox(correction_method),
      background = custom_bg,
      output = jsonlite::unbox("json")
    )
  ),
  auto_unbox = FALSE,
  null = "null")
  # Headers

  headers <-
    list("Accept" = "application/json",
         "Content-Type" = "application/json",
         "charset" = "UTF-8")

  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer(.mapUnicode = FALSE)
  h2 = RCurl::getCurlHandle() # Get info about the request

  # Request
  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()

  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }

  res <- jsonlite::fromJSON(txt)
  df = res$result
  meta = res$meta

  if (is.null(dim(df))) {
    message("No results to show\n", "Please make sure that the organism is correct or set significant = FALSE")
    return(NULL)
  }

  # Re-order the data frame columns
  if (multi_query) {
    # add column that shows if significant
    df$significant <- lapply(df$p_values, function(x) x <= user_threshold)

    col_names <- c(
      "term_id",
      "p_values",
      "significant",
      "term_size",
      "query_sizes",
      "intersection_sizes",
      "source",
      "term_name",
      "effective_domain_size",
      "source_order",
      "parents"
    )

  } else {
    col_names <- c(
      "query",
      "significant",
      "p_value",
      "term_size",
      "query_size",
      "intersection_size",
      "precision",
      "recall",
      "term_id",
      "source",
      "term_name",
      "effective_domain_size",
      "source_order",
      "parents"
    )

    if (evcodes) {
      col_names <- append(col_names, c("evidence_codes", "intersection"))
      # Add column that lists the intersecting genes separated by comma
      df$intersection <- mapply(
        function(evcodes, query){
          genemap <- meta$genes_metadata$query[[query]]$mapping
          genes <- meta$genes_metadata$query[[query]]$ensgs[which(lengths(evcodes) > 0)]
          genes2 <- lapply(genes, function(x) ifelse(x %in% genemap, names(which(genemap == x)) , x))
          return(paste0(genes2, collapse = ","))
        },
        df$intersections,
        df$query,
        SIMPLIFY = TRUE
      )
      df$evidence_codes <- sapply(df$intersections, function(x)
        paste0(sapply(x[which(lengths(x) > 0)], paste0, collapse = " "), collapse = ","), USE.NAMES = FALSE)
    }

    # Order by query, source and p_value
    df <- df[with(df, order(query, source, p_value)), ]

  }

  # Rename native to term_id
  colnames(df)[colnames(df) == "native"] <- "term_id"
  colnames(df)[colnames(df) == "name"] <- "term_name"

  # Fix the row indexing to start from 1
  row.names(df) <- NULL
  #df <- df[, col_names]
  return(list("result" = df, "meta" = meta))
}


###### Implement the hybrid filtering for gP results
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

hybrid_filtering = function(gostres, group_nr){
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

  return(keep)
}

# Sample query

gostres <- my_gost(c("GO:0005005", "GO:0035184"), evcodes = T, sources = "GO")
res = gostres$result
filtered_res = list()
for (gr in unique(res$group_id)){
  filtered_res[[gr]] = hybrid_filtering(gostres, gr)
}

## Glioblastoma example gene queries


data = read.xls("aaf2666_Table-S3.xlsx")
# data obtained from http://www.sciencemag.org/content/360/6389/660/suppl/DC1?_ga=2.80703895.187038585.1565015276-2016237288.1564921104
data$gene = as.character(data$gene)
query = data %>% group_by(ROI) %>% summarise(query = paste0(gene, collapse = " ")) %>% data.frame()

query2 = as.list(query$query)
names(query2) = query$ROI
query2 = lapply(query2, function(x) strsplit(x, " ")[[1]])

# CT
gostres1 <- my_gost(query2[[1]], evcodes = T, sources = "GO")
res = gostres1$result
filtered_res1 = list()
for (gr in unique(res$group_id)){
  filtered_res1[[gr]] = hybrid_filtering(gostres1, gr)
}
subres = res[res$term_id %in% unlist(filtered_res1), c("term_name", "group_id", "p_value", "term_size", "intersection_size", "query_size", "source")]
subres %>% arrange(group_id)

#CTmvp
gostres2 <- my_gost(query2[[2]], evcodes = T, sources = "GO")
res = gostres2$result
filtered_res2 = list()
for (gr in unique(res$group_id)){
  print(gr)
  if(gr == 1){
    filtered_res2[[gr]] = hybrid_filtering(gostres2, gr)
  }
  else{
    filtered_res2[[gr]] = hybrid_filtering(gostres2, gr)
  }

}

subres = res[res$term_id %in% unlist(filtered_res2), c("term_name", "group_id", "p_value", "term_size", "intersection_size", "query_size", "source")]
subres %>% arrange(group_id)

#CTpan
gostres3 <- my_gost(query2[[3]], evcodes = T, sources = "GO")
res = gostres3$result
filtered_res3 = list()
for (gr in unique(res$group_id)){
  filtered_res3[[gr]] = hybrid_filtering(gostres3, gr)
}

subres = res[res$term_id %in% unlist(filtered_res3), c("term_name", "group_id", "p_value", "term_size", "intersection_size", "query_size", "source")]
subres %>% arrange(group_id)

#LE
gostres4 <- my_gost(query2[[4]], evcodes = T, sources = "GO")
res = gostres4$result
filtered_res4 = list()
for (gr in unique(res$group_id)){
  filtered_res4[[gr]] = hybrid_filtering(gostres4, gr)
}
subres = res[res$term_id %in% unlist(filtered_res4), c("term_name", "group_id", "p_value", "term_size", "intersection_size", "query_size", "source")]
subres %>% arrange(group_id)



