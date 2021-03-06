MINSUPPORT <- 4

create_index_page <- function(outdir="./", runtag="run", codedir="./",
                              indexpath="index.html",
                              glossarypath="glossary.html",
                              imgstr="imgs/", txtstr="txts/"){

  dir.create(file.path(outdir))
  htmldir <- paste0(outdir, runtag, "/")
  dir.create(file.path(htmldir))
  file.copy(from = paste0(codedir, "sorttable.js"), to = htmldir)
  dir.create(file.path(paste0(htmldir, imgstr)))
  dir.create(file.path(paste0(htmldir, txtstr)))

  glossary <- as.matrix(read.table(paste0(codedir, "glossary.txt"),
                                   header = TRUE, sep = "\t", quote = ""))
  file.copy(from = paste0(codedir, "glossary.txt"), to = htmldir)
  abspath <- paste0(htmldir, indexpath)

  write(paste0("<br>"), file = abspath)
  write_html_table_page(glossary, paste0(htmldir, glossarypath), glossarypath)

  return(list(htmldir = htmldir, indexpath = indexpath, imgstr = imgstr,
              txtstr = txtstr, glossarypath = glossarypath, abspath = abspath))
}



showfirstlast <- function(mynames){
  show(c(length(mynames), length(unique(mynames)),
         sort(mynames)[c(1, length(mynames))]))
}


table2html <- function(resultstable){

  colheader <- paste(collapse = "</th><th>", c("idx", colnames(resultstable)))
  theader <- paste0("<thead>\n<tr bgcolor='#AAAAAA';><th>", colheader,
                    "</th></tr>\n</thead>\n<tbody>\n")
  head_str <- paste0("<script src='sorttable.js'></script>\n<TABLE class=",
                    "'sortable' border =1 >\n", theader)

  numtable <- cbind(1:dim(resultstable)[1], resultstable)
  rowinnerstrs <- apply(numtable, 1, paste, collapse = "</td><td>")
  rowstrs <- paste0("<tr><td>", rowinnerstrs, "</td></tr>")
  rows_str <- paste(collapse = "\n", rowstrs)

  return(paste0(collapse = "\n", head_str, rows_str, "\n</tbody></table><br>"))
}


cumSumTable <- function(mydataframe){
  dftable <- table(as.data.frame(as.matrix(mydataframe)))
  idxnums <- sort(as.numeric(rownames(dftable)), decreasing = T)
  dftable <- dftable[as.character(idxnums), ]
  if (is.null(dim(dftable))){ # if only one column or row
    nCount <- dftable
    cumCount <- cumsum(dftable)
    #show(tmptab)
    return(cbind(nCount, cumCount))
  } else {
    cumul_table <- apply(dftable, 2, cumsum)
    nCount <- rowSums(dftable)
    cumCount <- rowSums(cumul_table)
    dftable_final <- cbind(nCount,
                          cumCount,
                          signif(cumul_table / rowSums(cumul_table) * 100, 3)
    )
    #show(dftable_final)
    return(dftable_final)
  }
}

summarize_rungraphs <- function(rungraphs = NULL, iso_table = NULL,
                                weighted_chip_evidence = NULL,
                                graphstr = "data",
                                htmlinfo = list(htmldir = "html/",
                                                indexpath = "index.html",
                                                txtstr = "txts/"),
                                runinfo = list(nregs = "1",
                                               ntargets = "5",
                                               nboots = "10")){

  show(paste0("Processing ", graphstr, " graph method ..."))
  write(paste0("<br><br><b>", graphstr, "</b><br><br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = T)

  # summarize edge results
  edgesinfo <- edgeinfo_from_graphs(rungraphs = rungraphs,
                                    iso_table = iso_table,
                                    weighted_chip_evidence = weighted_chip_evidence)
  edgesinfo <<- edgesinfo

  # write edge results summary
  write_tables_all(edgesinfo, tabletype = "edges",
                   html_cols = c("edgekey", "weight", "reg-origid",
                                 "target-origid", "reg-geneid",
                                 "target-geneid", "chip-evidence",
                                 "num-chip-peaks"),
                   html_idxs = 1:min(1000, dim(edgesinfo)[1]),
                   filestr = graphstr,
                   htmlinfo = htmlinfo)

  # regulators summary
  regsinfo <- genetable_summary(filter = "reg-origid",
                                filterNeigh = "target-origid",
                                tabletype = "regulators",
                                nboots = as.integer(runinfo$nboots),
                                minsupport = MINSUPPORT,
                                edgesinfo = edgesinfo,
                                graphstr = graphstr,
                                htmlinfo = htmlinfo)

  # targets summary
  targetsinfo <- genetable_summary(filter = "target-origid",
                                   filterNeigh = "reg-origid",
                                   tabletype = "targets",
                                   nboots = as.integer(runinfo$nboots),
                                   minsupport = MINSUPPORT,
                                   edgesinfo = edgesinfo,
                                   graphstr = graphstr,
                                   htmlinfo = htmlinfo)

  # summarize edge results
  show(paste0("Summarizing ", graphstr, " graph method ..."))
  ngraphmods <- length(rungraphs)
  myweights <- as.numeric(edgesinfo[, "weight"])
  myevidence <- as.numeric(edgesinfo[, "chip-evidence"])

  write(paste0("<br>", ngraphmods, " graphs with an average of ",
               signif(sum(myweights) / ngraphmods, 3), " edges each.<br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = T)
  write(paste0(sum(myweights), " (",
               signif(sum(myweights) / ( runinfo$nregs * runinfo$ntargets *
                                         runinfo$nboots) * 100, 3),
               "%) of ", runinfo$nregs, "*", runinfo$ntargets, "*",
               runinfo$nboots, " possible edges across bootstraps.<br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = T)
  write(paste0(length(myweights), " (",
               signif(length(myweights) / (runinfo$nregs *
                                           runinfo$ntargets) * 100, 3),
               "%) of ", runinfo$nregs, "*", runinfo$ntargets,
               " unique edges found.<br><br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = T)

  evid_weight_raw <- cbind(c(myweights), c(myevidence))
  if (length(unique(myweights)) == 1 || length(unique(myweights)) == 1){
    evid_weight_raw <- rbind(evid_weight_raw, c(0, -1))
  }
  evidence_by_weight_final <- cumSumTable(evid_weight_raw)
  show(evidence_by_weight_final)

  summary_final <- cbind(rownames(evidence_by_weight_final),
                        evidence_by_weight_final)
  if (dim(summary_final)[2] == 6){
    colnames(summary_final) <- c("Support", "nEdges", "cumEdges",
                              "%NA", "%NoPeak", "%Peaks")
  } else {
    colnames(summary_final) <- c("Support", "nEdges", "cumEdges")
  }

  htmlstr <- table2html(summary_final)
  write(htmlstr, file = paste0(htmlinfo$htmldir, htmlinfo$indexpath),
        append = T)
  return(summary_final)
}

edgeinfo_from_graphs <- function(rungraphs, iso_table, weighted_chip_evidence){

  show("Extracting edges from all graphs...")
  edgelist_list <- lapply(rungraphs, extract_edge_strings_from_graph)
  edgestable <- table(unlist(edgelist_list))

  reg_origids <- unlist(lapply(strsplit(names(edgestable), "\\|\\|"), "[[", 2))
  tar_origids <- unlist(lapply(strsplit(names(edgestable), "\\|\\|"), "[[", 1))
  show("reg_origids: ")
  showfirstlast(reg_origids)
  show("tar_origids: ")
  showfirstlast(tar_origids)

  # get gene ids
  reg_gids <- as.character(iso_table[reg_origids, "iso_ensgs"])
  tar_gids <- as.character(iso_table[tar_origids, "iso_ensgs"])
  show("reg_gids: ")
  showfirstlast(reg_gids)
  show("tar_gids: ")
  showfirstlast(tar_gids)

  # get chip evidence
  raw_evidence <- rep(-1, length(edgestable))
  if (!is.null(weighted_chip_evidence)){
    show(c("...Extracting chip evidence"))
    raw_evidence <- weighted_chip_evidence[cbind(gsub("-", ".", reg_gids),
                                                 gsub("-", ".", tar_gids))]
  }
  binary_evidence <- raw_evidence
  binary_evidence[binary_evidence > 1] <- 1
  show(table(binary_evidence))
  sortval <- edgestable * 1000000 + raw_evidence

  # combine edge summary table
  edgesinfo <- cbind(names(edgestable), edgestable, reg_origids, tar_origids,
                    reg_gids, tar_gids, binary_evidence, raw_evidence,
                    sortval)
  colnames(edgesinfo) <- c("edgekey", "weight", "reg-origid", "target-origid",
                          "reg-geneid", "target-geneid", "chip-evidence",
                          "num-chip-peaks", "sort-val")

  return(edgesinfo[sort(as.numeric(edgesinfo[, "sort-val"]),
                        decreasing = TRUE, index.return = T)$ix, ])
}

extract_edge_strings_from_graph <- function(mygraph){
  library("igraph")
  elist <- get.edgelist(mygraph)
  # small graphs lose node names when extracted
  if (is.numeric(elist)){
    return()
  } else {
    return(apply(elist, 1, paste, collapse = "||"))
  }
}

write_tables_all <- function(mytab, tabletype="table",
                             html_idxs=1:dim(mytab)[1],
                             html_cols=colnames(mytab),
                             filestr="data",
                             htmlinfo=list(htmldir = "html/",
                                           indexpath = "index.html",
                                           txtstr = "txts/")){
  htmlpath <- paste0(filestr, "_", tabletype, ".html")
  resultspath <- paste0(htmlinfo$txtstr, filestr, "_", tabletype, ".txt")
  show(paste0("Writing table: ", resultspath))
  write.table(mytab, paste0(htmlinfo$htmldir, resultspath), sep = "\t",
              row.names = F, col.names = T, quote = F)
  write(paste0('<a href = "', htmlpath, '" target="_blank">',
               tabletype, "</a><br>"),
        file = paste0(htmlinfo$htmldir, htmlinfo$indexpath), append = T)
  write_html_table_page(resultstable = mytab[html_idxs, html_cols],
                        htmlpagefile = paste0(htmlinfo$htmldir, htmlpath),
                        resultsrelpath = resultspath,
                        indexpath = htmlinfo$indexpath)
}


write_html_table_page <- function(resultstable, htmlpagefile, resultsrelpath,
                                  indexpath = "index.html",
                                  glossarypath = "glossary.html"){

  write(paste0("<table border = 1 width = '100%'><tr bgcolor = ",
               "'#AAAAAA'><th><a href = '", indexpath, "' target='_blank'>",
               "Index</a></th><th><a href = '", glossarypath, "' target=",
               "'_blank'>Glossary</a></th><th><a href = '", resultsrelpath,
               "' target='_blank'>Download</a></th></tr></table><br><br>"),
        file = htmlpagefile)
  htmlstr <- table2html(resultstable)
  write(htmlstr, file = htmlpagefile, append = T)
}



genetable_summary <- function(filter, filterNeigh, tabletype,
                              nboots = 10, minsupport = MINSUPPORT, edgesinfo,
                              graphstr, htmlinfo){
  show(paste0(tabletype, " processing..."))
  mygids <- unique(edgesinfo[, filter])
  # only extract top 1000 targets
  if (tabletype == "targets"){
    tmptable <- table(as.data.frame(edgesinfo[, c("target-origid",
                                                  "weight")] ))
    topgenes <- rownames(tmptable)
    names(topgenes) <- rownames(tmptable)
    if (nboots > 1){
      scaledtable <- t(t(as.matrix(tmptable)) * as.numeric(colnames(tmptable)))
      keepsupport <- nboots:minsupport
      colstrs <- intersect(colnames(scaledtable), as.character(keepsupport))
      topgenes <- rowSums(scaledtable[, as.character(colstrs)])
    }
    sortidxs <- sort(topgenes, decreasing = T, index.return = T)$ix
    mygids <- rownames(tmptable)[sortidxs[1:min(1000, length(topgenes))]]
  }
  ###
  genetable <- t(sapply(mygids, extract_gene_row, filter, filterNeigh,
                       edgesinfo, nboots))
  colnames(genetable) <- c("gid", "nEdges", "nNeigh", "nMax-Conf-Neigh",
                          "topNeigh", "suppTable")
  sortval <- (as.numeric(genetable[, "nMax-Conf-Neigh"]) * 1000000000 +
               as.numeric(genetable[, "nEdges"]))
  sortidxs <- sort(sortval, decreasing = T,
                   index.return = T)$ix[1:min(1000, length(sortval))]
  write_tables_all(genetable, tabletype = tabletype,
                   html_idxs = sortidxs, filestr = graphstr,
                   htmlinfo = htmlinfo)
  return(genetable)
}


extract_gene_row <- function(mygid, myfield, myfield2, edgesinfo,
                             nboots=10, maxneigh=3){
  #    mygid = "ENSG00000196092"
  #    myfield = "reg-geneid"
  #    myfield2 = "target-origid"
  #    mygid = "ENSG00000181804"
  #    myfield2 = "reg-origid"
  #    myfield = "target-geneid"

  gididxs <- which(edgesinfo[, myfield] == mygid)
  mysubtab <- edgesinfo[gididxs, c(myfield2, "weight", "num-chip-peaks",
                                   "chip-evidence")]
  if (length(gididxs) == 1){
    mysubtab <- t(as.matrix(mysubtab))
  }

  evid_weight_tab <- ""
  if (length(gididxs) > 1){
    evid_weight_tab <- cumSumTable(mysubtab[, c("weight", "chip-evidence")])
    evid_weight_tab <- evid_weight_tab[c(1, dim(evid_weight_tab)[1]), ]
  }

  return(c(mygid,
           sum(as.numeric(mysubtab[, "weight"])),
           length(gididxs),
           sum(as.numeric(mysubtab[, "weight"]) == nboots),
           tab2cell(mysubtab[1:min(maxneigh, length(gididxs)), ],
                    keeprnames = F),
           tab2cell(evid_weight_tab) ))
}


tab2cell <- function(mytab, keepheader=T, keeprnames=T){
  if (keepheader){ mytab <- rbind(colnames(mytab), mytab) }
  if (keeprnames){ mytab <- cbind(rownames(mytab), mytab) }
  return(paste(collapse = "<br>", apply(mytab, 1, paste, collapse = " | ")))
}
