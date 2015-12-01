#' Summarize results of beta-Poisson generalized linear models of a dataset
#'
#' @param object An object of BPglm class resulted from \code{\link{BPglm}} function
#' @param ... Inherited parameters of summary object
#' @return An summary.BPglm object contains a table of coefficients sorted by estimated false discovery rate (FDR)
#' @export
#' @examples
#' #See the example in the \code{\link{BPglm}} function
summary.BPglm <- function(object, ...){
  FDR=p.adjust(object$PVAL,method="BH")
  topTable=NULL;
  if (!object$keepFit){
    topTable=cbind(object$TVAL,object$PVAL,FDR)
    rownames(topTable)=names(object$PVAL)
    colnames(topTable)=c("t value","Pr(>|t|)","fdr")
  }else{
    fitRes=object$fitRes
    coef=object$coef
    fit.ids=which(names(fitRes)=="fit")
    topTable=NULL;
    for (i in 1:length(fit.ids)){
      coefficients.table=summary(fitRes[fit.ids[i]]$fit)$coefficients[2,]
      topTable=rbind(topTable,coefficients.table)
    }
    rownames(topTable)=names(object$PVAL)
    topTable=cbind(topTable,FDR)
    colnames(topTable)=c("Estimate", "Std. Error","t value", "Pr(>|t|)","fdr")
  }
  
  topTable=topTable[order(FDR),]
  cat("Top five fdr: \n")
  print( topTable[1:5,])
  
  ans=list()
  ans$topTable=topTable  
  class(ans) <- "summary.BPglm"
  return(ans)
}
