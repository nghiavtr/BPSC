#' Get Monte-Carlo p-values from the BP models and the Monte-Carlo null distributions
#'
#' @param bp.model.list List of beta-Poisson models that are results from \code{\link{estimateBPMatrix}} function
#' @param MCdis.list List of Monte-Carlo null distributions that are results of \code{\link{getBPMCnullmatrix}} function
#' @return A list of Monte-Carlo p-values
#' @export
#' @examples
#' #See \code{\link{getBPMCnullmatrix}} function
getMCpval<-function(bp.model.list=NULL,MCdis.list=NULL){
  observed.X2.ids=which(names(bp.model.list)=="X2")        
  MCdis.X2.ids=which(names(MCdis.list)=="X2")

  MCpval.list=NULL;
  for (i in 1:length(observed.X2.ids)){
    pval=0.0
    if (!is.null(MCdis.list[[MCdis.X2.ids[i]]][1]) && !is.na(MCdis.list[[MCdis.X2.ids[i]]][1])){
      pval=sum(MCdis.list[[MCdis.X2.ids[i]]] >= bp.model.list[observed.X2.ids[i]])/length(MCdis.list[[MCdis.X2.ids[i]]])
    }
    MCpval.list=c(MCpval.list,pval)
  }
  return(MCpval.list)
}