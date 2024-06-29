# Viz functions from R/ipr2viz.R, R/networks_domarch.R

## Function to obtain element counts (DA, GC)
count_bycol <- function(prot=prot, column="DomArch", min.freq=1) {
  counts <- prot |>
    select(column) |>
    table() |>
    as_tibble() |>
    `colnames<-`(c(column, "freq")) |>
    filter(!grepl("^-$", column)) |>		# remove "-"
    filter(!is.na(column)) |>
    arrange(-freq) |> filter(freq>=min.freq)
  return(counts)
}

## Function to break up ELEMENTS to WORDS for DA and GC
elements2words <- function(prot, column= "DomArch",
                           conversion_type="da2doms") {
  z1 <- prot |>
    select(column) |>
    str_replace_all("\\,"," ") |>
    str_replace_all("\""," ")
  switch(conversion_type,
         da2doms = { z2 <- z1 |>
           str_replace_all("\\+"," ")},
         gc2da = { z2 <- z1 |>
           str_replace_all("\\<-"," ") |>
           str_replace_all("-\\>"," ") |>
           str_replace_all("\\|"," ")})
  # str_replace_all("^c\\($", " ") |>		# remove "c("
  # str_replace_all("\\)$", " ") |>			# remove ")"
  # str_replace_all("\\(s\\)"," ") |>		# Ignoring repeats
  # str_replace_all("-"," ") |>
  ## replace \n, \r, \t
  z3 <- z2 |>
    str_replace_all("\n"," ") |>
    str_replace_all("\r"," ") |>
    str_replace_all("\t"," ") |>
    ## replace multiple spaces ...
    str_replace_all("    "," ") |>
    str_replace_all("   "," ") |>
    str_replace_all("  "," ") |>
    str_replace_all("  "," ")
  return(z3)
}

## Function to get WORD COUNTS [DOMAINS (DA) or DOMAIN ARCHITECTURES (GC)]
## to be used after elements2words
words2wc <- function(x){ x |>
    str_replace_all("   "," ") |>
    str_replace_all("  "," ") |> str_replace_all("  "," ") |>
    paste(collapse=" ") |>
    strsplit(" ") |>
    # filter(grepl(query.list[j], Query)) |> # Create separate WCs for each Query
    # select(DA.wc) |>
    table() |> as_tibble() |>
    `colnames<-`(c("words", "freq")) |>
    ## filter out 'spurious-looking' domains
    filter(!grepl(" \\{n\\}", words)) |>
    filter(!grepl("^c\\($", words)) |>		# remove "c("
    filter(!grepl("^\\)$", words)) |>		# remove ")"
    filter(!grepl("^-$", words)) |>			# remove "-"
    filter(!grepl("^$", words)) |>				# remove empty rows
    filter(!grepl("^\\?$", words)) |>		# remove "?"
    filter(!grepl("^\\?\\*$", words)) |>	# remove "?*"
    filter(!grepl("^tRNA$", words)) |>		# remove "tRNA"
    filter(!grepl("^ncRNA$", words)) |>	# remove "ncRNA"
    filter(!grepl("^rRNA$", words)) |>		# remove "rRNA"
    filter(!grepl("^X$|^X\\(s\\)$", words)) |> # remove "X" and "X(s)"

    # filter(!grepl("\\*", words)) |>			# Remove/Keep only Query
    arrange(-freq)
}
## Function to filter based on frequencies
filter_freq <- function(x, min.freq){ x |>
    filter(freq>=min.freq)
}

total_counts <- function(prot, column = "DomArch", #lineage_col = "Lineage",
                         cutoff = 90, RowsCutoff = FALSE, digits = 2
                         #type = "GC"
){
  #'Total Counts
  #'
  #'Creates a data frame with a totalcount column
  #'
  #'This function is designed to sum the counts column by either Genomic Context or Domain Architecture and creates a totalcount column from those sums.
  #'
  #' @param prot A data frame that must contain columns:
  #' \itemize{\item Either 'GenContext' or 'DomArch.norep' \item count}
  #' @param cutoff Numeric. Cutoff for total count. Counts below cutoff value will not be shown. Default is 0.
  #' @param column Character. The column to summarize
  #' @examples total_counts(pspa-gc_lin_counts,0,"GC")
  #' @note Please refer to the source code if you have alternate file formats and/or
  #' column names.
  column <- sym(column)

  prot <- select(prot, {{column}}) |> #, {{lineage_col}}
    filter(!is.na({{column}})) |> # & !is.na({{lineage_col}})
    filter({{column}} != "")

  #prot <- summarize_bylin(prot, column, by =  lineage_col, query = "all")
  col_count <-  prot |>
    group_by({{column}}) |>
    summarise(totalcount = sum(count))

  total <- left_join(prot,col_count, by = as_string(column))

  sum_count <- sum(total$count)
  total <- total |>
    mutate("IndividualCountPercent" = totalcount/sum_count*100) |>
    arrange(-totalcount,-count)

  cumm_percent <- total |>
    select({{column}}, totalcount) |>
    distinct() |>
    mutate("CumulativePercent"=0)

  total_counter = 0

  for(x in length(cumm_percent$totalcount):1){
    total_counter = total_counter + cumm_percent$totalcount[x]
    cumm_percent$CumulativePercent[x] = total_counter/sum_count * 100
  }

  cumm_percent <- cumm_percent|> select( CumulativePercent, {{column}})

  total <- total |>
    left_join(cumm_percent, by = as_string(column))

  # Round the percentage columns
  total$CumulativePercent <- total$CumulativePercent |>
    round(digits = digits)
  total$IndividualCountPercent <- total$IndividualCountPercent |>
    round(digits = digits)

  if(RowsCutoff)
  {
    # If total counts is being used for plotting based on number of rows,
    # don't include other observations that fall below the cummulative percent cutoff
    #, but that have the same 'totalcount' number as the cutoff observation
    total <- total |>
      filter(CumulativePercent >= 100-cutoff-.0001)
    return(total)
  }

  # Include observations that fall below the cummulative percent cutoff,
  # but that have the same 'totalcount' as the cutoff observation
  t <- total |>
    filter(CumulativePercent >= 100-cutoff)
  if(length(t) == 0){
    cutoff_count = 0
  }
  else{
    cutoff_count = t$totalcount[nrow(t)]
  }

  total <- total |>
    filter(totalcount >= cutoff_count) |>
    ungroup()

  return(total)
}


BinaryDomainNetwork <- function(prot, column = "DomArch",
                                #domains_of_interest,
                                cutoff = 70,
                                layout = "nice",
                                query_color = adjustcolor("yellow", alpha.f = .5),
                                partner_color = adjustcolor("skyblue", alpha.f = .5),
                                border_color = adjustcolor("grey", alpha.f = .8),
                                IsDirected = T)
{
  #'Domain Network
  #'
  #'This function creates a domain network from the 'DomArch' column.
  #'
  #'Domains that are part of the 'domains_of_interest' are a different node color than the other domains.
  #'
  #'A network of domains is returned based on shared domain architectures.
  #'
  #'@param prot A data frame that contains the column 'DomArch'.
  #'@param column Name of column containing Domain architecture from which nodes and edges are generated.
  #'@param cutoff Integer. Only use domains that occur at or above the cutoff for total counts if cutoff_type is "Total Count".
  #'Only use domains that appear in cutoff or greater lineages if cutoff_type is Lineage.
  #'@param layout Character. Layout type to be used for the network. Options are:
  #'\itemize{\item "grid" \item "circle" \item "random" \item "auto"}
  #'@param query_color Color that the nodes of the domains in the domains_of_interest vector are colored
  #'@param partnercolor Color that the nodes that are not part of the domains_of_interest vector are colored
  #'@param IsDirected Is the network directed? Set to false to eliminate arrows
  #'@examples domain_network(pspa)
  # by domain networks or all, as required.
  #print(domains_of_interest)

  column_name <- sym(column)

  prot_tc <- prot |>
    total_counts(column =  column, cutoff = cutoff, RowsCutoff = F, digits = 5)

  within_list <- prot_tc |>
    select({{column_name}}) |>
    distinct()
  within_list <- pull(within_list, {{column_name}})

  prot <- prot |>
    filter({{column_name}} %in% within_list)

  ####### Below should be part of the standardized cleanup process
  prot$DomArch.ntwrk <- as_vector(prot |> select({{column}})) |> # testing with prot$DomArch.orig
    str_replace_all(coll(pattern="\\?", ignore_case=T), "X")

  domains_of_interest_regex = paste(domains_of_interest, collapse = "|")
  # domain.list <- prot |>
  #   dplyr::filter(grepl(pattern=domains_of_interest_regex,
  #                       x=DomArch.ntwrk,
  #                       ignore.case=T))
  domain.list = prot

  ##Separating column and converting to atomic vector prevents coercion
  domain.list <- domain.list$DomArch.ntwrk |>
    str_split(pattern="\\+")

  # Get domain counts before eliminating domarchs with no edges
  wc = elements2words(prot = prot, column =  column,
                      conversion_type = "da2doms") |>
    words2wc()

  nodes = data.frame(id = wc$words, label = wc$words, size = wc$freq) |>
    mutate(group = purrr::map(id,
                              function(x)
                                ifelse(x %in% domains_of_interest, "Query", "Partner")))

  max_size = max(nodes$size)
  min_size = min(nodes$size)
  nodes <- nodes |>
    mutate(size = (size-min_size)/((max_size-min_size)) *20+10)
  max_font_size = 43
  nodes <- nodes |>
    mutate(font.size = purrr::map(size, function(x) min(x*2 , max_font_size)))

  domain.list = domain.list[-which(lengths(domain.list)==1)]

  if(length(domain.list) != 0)
  {
    from <- unlist(lapply(domain.list,
                          function(x) sapply(1:(length(x)-1),
                                             function(y) x[y]))) #list elements 1 through n-1
    to <- unlist(lapply(domain.list,
                        function(x) sapply(1:(length(x)-1),
                                           function(y) x[y+1]))) #list elements 2 through n
    pwise <- cbind(from,to)
    if(any(pwise=="")) {
      pwise=pwise[-which( (pwise[,1]=="") | (pwise[,2]=="") ),]
    }
    if(any(pwise == "X"))
      pwise=pwise[-which( (pwise[,1]=="X") | (pwise[,2]=="X") ),]

    edges <- data.frame(from = pwise[,1], to = pwise[,2]) |>
      group_by(from, to) |>
      summarize(width = n())
    edges <- edges |>
      mutate(width = ifelse(width==1, .3, log(width)))
    ew <- c(2.7,4.5)

    ColorEdges <- function(x)
    {
      if(x>=ew[1] && x<=ew[2])
      {
        adjustcolor("cadetblue", alpha.f = .7)
      }
      else if(x>ew[2])
      {
        adjustcolor("maroon", alpha.f = .5)
      }
      else
      {
        "gray55"
      }
    }

    edges <- edges |>
      mutate(color = unlist(purrr::map(width, ColorEdges)) )

  }

  if(IsDirected){
    vg <- visNetwork(nodes, edges, width = "100%", height = '600px') |>
      visEdges(arrows = 'to', smooth =T)
  }
  else
  {
    vg <- visNetwork(nodes, edges, width = "100%", height = '600px') |>
      visEdges(smooth =T)
  }
  vg <- vg |>
    visGroups(groupname = "Query", color = query_color) |>
    visGroups(groupname = "Partner", color = partner_color) |>
    visOptions(highlightNearest = TRUE) |>
    visLegend(position = "right",width = .1)

  vg <- switch(layout,
               "nice" = visIgraphLayout(vg, "layout_nicely" ),
               "random" = visIgraphLayout(vg, "layout_randomly"),
               "grid" = visIgraphLayout(vg, "layout_on_grid"),
               "circle" = visIgraphLayout(vg, "layout.circle"),
               "auto" =  visIgraphLayout(vg, "layout.auto")
  )
  vg
}
