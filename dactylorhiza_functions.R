

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
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)
  #allRes<-GenTable(sampleGOdata, Fis = resultFisher, topNodes = 20)
  return(allRes)
}

draw_heatmap<-function(dds_object){
  dds_fpkm<-fpkm(dds_object$dds)
  new_column_order<-dds_fpkm %>% colnames %>% sort()
  dds_fpkm <- dds_fpkm %>% data.frame() %>% dplyr::select(new_column_order)
  dds_significant<-dds_object$results %>% data.frame() %>% filter(log2FoldChange > 2 & padj < 0.05) %>% rownames()
  pheatmap(dds_fpkm[rownames(dds_fpkm) %in% dds_significant,], 
           scale = "row", 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_rownames = FALSE)
}

get_significant_genes<-function(results_object, fold_change=2, pvalue=0.05){
  de_genes<-results_object$results %>% data.frame() %>% filter(log2FoldChange > fold_change & padj < pvalue) %>% rownames()
  return(de_genes)
}



