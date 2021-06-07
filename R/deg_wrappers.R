#' Wrapper function to create dge from txi
#'
#' Wrapper which returns a dge that has been filtered
#' and is ready for downstream processing. I've
#' added a bunch of errors here which should prevent
#' mistakes in creating/using this object
#'
#' @param txi 'lengthScaledTPM' txi object
#' @param p data.frame of sample phenodata
#' @param gene_annots gene annotations
#'
#' @return DGE object
#' @export
#'
#' @examples
dge_wrapper <- function(txi, p, gene_annots){
  #TODO Error messages don't output if condition
  # need to fix that
  if(all(colnames(txi$counts) != row.names(p))){
    stop(call. = TRUE)
  }
  if(! all(row.names(txi$counts) %in% row.names(gene_annots))){
    stop(call. = TRUE)
  }

  dge <- edgeR::DGEList(counts =  txi$counts, samples = p, genes = gene_annots[row.names(txi$counts),])
  dge <- dge[edgeR::filterByExpr(dge, group = dge$samples$group),]
  dge <- dge[(rowSums(dge$counts > 1) >= .75*dim(dge)[2]),]
  dge <- dge[!is.na(dge$genes$entrezgene),] %>%
    edgeR::calcNormFactors(.)
  return(dge)

}

#' Runs voom to sva pipeline
#'
#' This runs the voom to sva pipeline then voom
#' again pipeline. Returns a voom transformed
#' object fitted with the surrogate variables
#'
#' @param dge DGE object
#' @param mod model.matrix with variables of interest
#' @param mod0 model.matrix without variables of interest
#'
#' @return voom object
#' @export
#'
#' @examples
voom_sva_wrapper <- function(dge, mod, mod0){
  v <- voom(dge, design = mod)
  # v_qa <- voomWithQualityWeights(dge,design = mod)
  svobj <- sva(v$E, mod=mod, mod0=mod0)
  modSV <- cbind(mod, svobj$sv)
  v_sv <- voom(dge, design = modSV)
  return(v_sv)

}

# DEG related functions

#' Over Representation analysis (ORA)
#'
#' ORA wrapper which splits the long res list into
#' combined, up-regulated, and down-regulated list of
#' genes. It then performs ORA using clusterProfiler
#'
#' @param res_object res object in long(ish) format
#' @param lfc_thresh lfc threshold  absolute value
#'
#' @return data.frame of enriched pathways
#' @export
#'
#' @examples
enrich_wrapper <- function(res_object, lfc_thresh=0){
  res_object2 <- res_object[abs(res_object$logFC) >= lfc_thresh,] %>%
    dplyr::filter(., adj.P.Val < .05)
  temp_reslist <- list(
    combined=res_object2,
    upreg=res_object2[res_object2$logFC > 0,],
    downreg=res_object2[res_object2$logFC < 0,]
  )
  for(i in temp_reslist){
    print( paste(lfc_thresh,paste0(dim(i), collapse = ' ')))
  }
  ol <- lapply(temp_reslist, function(temp_res){
    temp_list <-  list(
      try(DOSE::enrichDGN(temp_res$entrezgene,  readable = TRUE, pvalueCutoff = .1, qvalueCutoff = .1)),
      try(DOSE::enrichDO(temp_res$entrezgene, readable = TRUE, pvalueCutoff = .1, qvalueCutoff = .1)),
      try(ReactomePA::enrichPathway(temp_res$entrezgene,readable = TRUE, pvalueCutoff = .1, qvalueCutoff = .1 )),
      try(clusterProfiler::enrichKEGG(temp_res$entrezgene, pvalueCutoff = .1, qvalueCutoff = .1))) %>% setNames(., c('DGN', 'DO','Reactome', 'KEGG' )) %>%
      lapply(., data.frame)
  }) %>% lapply(., function(x){
    rbind_named_df_list(x, 'database')
  })# %>%
  #rbind_named_df_list(., 'regulation')
  return(ol)
}
