

# make a function which takes the original full data and gives back the counts and sampke info in a nmaed list based on selection statement
specify_comparison<-function(samples, counts, selection_statement){
  samples_return<-samples %>% filter(rlang::eval_tidy(rlang::parse_expr(selection_statement)))
  counts_return <- counts %>% dplyr::select(rownames(samples_return))
  return(list(samples=samples_return, counts=counts_return))
}


run_diffexp<-function(comparison, design_term, gene_lengths, cpm_threshold=5, min_count_per_sample=8){
  formula_parsed<-paste("~", design_term)
  dds <- DESeqDataSetFromMatrix(countData = comparison[["counts"]],
                                colData = comparison[["samples"]],
                                design = as.formula(formula_parsed))
  mcols(dds)$basepairs<-gene_lengths
  dds <- estimateSizeFactors(dds)
  idx <- rowSums( counts(dds, normalized=TRUE) >= cpm_threshold ) >= min_count_per_sample
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
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100,
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
  #sapply(value_list, function(x) dds_object[[x]]) %>% data.frame()
}

# annotation_values here are heatmpa annotations, adding a coloured bar 
# along the top indicating different groups
draw_heatmap<-function(dds_object, annotation_values = NULL, custom=FALSE){
  # use custom if youre just passing in an FPKM custom data frame
  if (custom == FALSE){
    dds_fpkm<-fpkm(dds_object$dds)
    new_column_order<-dds_fpkm %>% colnames %>% sort()
    dds_fpkm <- dds_fpkm %>% data.frame() %>% dplyr::select(new_column_order)
    dds_significant<-dds_object$results %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.0005) %>% rownames()
    dds_toplot<-dds_fpkm[rownames(dds_fpkm) %in% dds_significant,]
  } else if (custom == TRUE) {
    new_column_order<-dds_object %>% colnames %>% sort()
    dds_object <- dds_object %>% data.frame() %>% dplyr::select(new_column_order)
    dds_toplot<-dds_object
  }
  if(!(is.null(annotation_values))){
    annotation_df <- construct_heatmap_annotation_df(dds_object, annotation_values)
    print(annotation_df)
    rownames(annotation_df) <- colnames(dds_fpkm)
  } else {
    annotation_df <- NULL
  }
  pheatmap::pheatmap(dds_toplot, 
           scale = "row", 
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           treeheight_row = 0, treeheight_col = 0,
           annotation_col = annotation_df)
}




draw_heatmap2 <-function(dds){
  # get DE results
  de_genes<-dds %>% results() %>% data.frame() %>% filter(abs(log2FoldChange) > 2 & padj < 0.0005) %>% rownames() 
  gene_counts<-counts(dds, normalized=TRUE)
  de_gene_counts<-gene_counts[de_genes,]
  
  # make the columns alphabetical to sort conditions
  new_column_order<-de_gene_counts %>% colnames %>% sort()
  de_gene_counts <- de_gene_counts %>% data.frame() %>% dplyr::select(new_column_order)
  
  # z score normalise data
  #heatmap_data <- t(apply(de_gene_counts, 1, function(x) x/sd(x)))
  heatmap_data <- t(apply(de_gene_counts, 1, function(x) x/mean(x)))
  
  # log normalise, remove NAs
  heatmap_data <- log2(heatmap_data)
  # get data range for breaks
  max_data <- max(heatmap_data, na.rm = TRUE)
  min_data <- -min(heatmap_data, na.rm = TRUE)
  range <- min(max_data, min_data)
  heatmap_data[is.infinite(heatmap_data)] <- NA
  heatmap_data[is.nan(heatmap_data)] <- NA
  heatmap_data %<>% data.frame() %>% drop_na()
  
  # plot heatmap
  pheatmap(heatmap_data,
           breaks = seq(-range, range, length.out = 100),
           cluster_rows = TRUE, cluster_cols = FALSE,
           treeheight_row = 0, treeheight_col = 0,
           show_rownames = F, show_colnames = T, scale="none")
}



get_significant_genes<-function(results_object, fold_change=2, pvalue=0.0005){
  de_genes<-results_object$results %>% data.frame() %>% filter(abs(log2FoldChange) > !!fold_change & padj < !!pvalue) %>% rownames()
  return(de_genes)
}


