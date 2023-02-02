

# make a function which takes the original full data and gives back the counts and sampke info in a nmaed list based on selection statement
specify_comparison<-function(samples, counts, selection_statement){
  samples_return<-samples %>% filter(rlang::eval_tidy(rlang::parse_expr(selection_statement)))
  counts_return <- counts %>% dplyr::select(rownames(samples_return))
  return(list(samples=samples_return, counts=counts_return))
}


run_diffexp<-function(comparison, design_term, gene_lengths){
  formula_parsed<-paste("~", design_term)
  dds <- DESeqDataSetFromMatrix(countData = comparison[["counts"]],
                                colData = comparison[["samples"]],
                                design = as.formula(formula_parsed))
  mcols(dds)$basepairs<-gene_lengths
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= 8 ) >= 5
  dds <- dds[idx,]
  dds <- DESeq(dds)
  res <- results(dds)
  return(list(dds=dds, results=res))
}

get_dds_object<-function(comparison, design_term, gene_lengths, cpm, min_samples, vst=FALSE){
  formula_parsed<-paste("~", design_term)
  dds <- DESeqDataSetFromMatrix(countData = comparison[["counts"]],
                                colData = comparison[["samples"]],
                                design = as.formula(formula_parsed))
  mcols(dds)$basepairs<-gene_lengths
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= cpm ) >= min_samples
  dds <- dds[idx,]
  if (vst == TRUE) {
    dds<-varianceStabilizingTransformation(dds)
    print("running vst")
  }
  return(dds)
}


get_enriched_terms<-function(gene_list, mappings){
  # use the gene 2 GOterms mapping provided for D. incarnata
  geneID2GO<-mappings
  # the input genes form the input, use these to annotate all genes, 1 is present in input list, 0 is absent
  geneSel<-gene_list
  geneSel<-factor(as.integer(names(geneID2GO) %in% geneSel))
  names(geneSel)<-names(geneID2GO)
  
  # set up the topGO object
  sampleGOdata <- new("topGOdata",
                      ontology = "BP",
                      allGenes = geneSel, 
                      nodeSize = 10,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO)
  
  # run three tests, fisher, Kol-Smirn, and Kol-Smirn with elimination
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
  
  # generate summary tane and return it
  allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                     classicKS = resultKS, elimKS = resultKS.elim,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 30,
                     numChar=1000 )
  #allRes<-GenTable(sampleGOdata, Fis = resultFisher, topNodes = 20)
  return(allRes)
}

return_enrichment_table<-function(de_object){
  mp<-readMappings("/Users/katieemelianova/Desktop/Dactylorhiza/data/all_annotations_justGO.txt")
  top10_enriched<-get_significant_genes(de_object, mappings_format=TRUE) %>% 
    get_enriched_terms(mp) %>% 
    data.frame() %>% 
    filter(classicFisher < 0.05) %>% 
    dplyr::select(Term, Annotated, Significant, Expected)
  return(top10_enriched)
}

# a way to make a data frame of values to annotate a heatmap with on the fly
construct_heatmap_annotation_df<-function(dds_object, value_list){
  sapply(value_list, function(x) dds_object[["dds"]][[x]]) %>% data.frame()
}

# annotation_values here are heatmpa annotations, adding a coloured bar 
# along the top indicating different groups
draw_heatmap<-function(dds_object, annotation_values = NULL){
  dds_fpkm<-fpkm(dds_object$dds)
  new_column_order<-dds_fpkm %>% colnames %>% sort()
  dds_fpkm <- dds_fpkm %>% data.frame() %>% dplyr::select(new_column_order)
  dds_significant<-dds_object$results %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% rownames()
  dds_toplot<-dds_fpkm[rownames(dds_fpkm) %in% dds_significant,]
  if(!(is.null(annotation_values))){
    annotation_df <- construct_heatmap_annotation_df(dds_object, annotation_values)
    rownames(annotation_df) <- colnames(dds_toplot)
  } else {
    annotation_df <- NULL
  }
  pheatmap(dds_toplot, 
           scale = "row", 
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           treeheight_row = 0, treeheight_col = 0,
           annotation_col = annotation_df)
}


get_significant_genes<-function(results_object, fold_change=2, pvalue=0.05, mappings_format=FALSE){
  de_genes<-results_object$results %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.05) %>% rownames()
  # to make the gene names match with the names in the GO term mappings, we need to cut off the last bit (optional)
  if (mappings_format == TRUE){
    de_genes<-sapply(de_genes, function(x) str_split(x, ":")[[1]][1]) %>% unname()
  }
  return(de_genes)
}


