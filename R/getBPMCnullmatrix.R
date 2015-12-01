#' Generate Monte-Carlo null distributions for a list of beta-Poisson models
#'
#' @param bp.model.list List of beta-Poisson models that are results from \code{\link{estimateBPMatrix}} function
#' @param fout A *.RData file name to export results
#' @param sim.num A number of simulation of each model
#' @param useParallel An option for using parallel (=TRUE)
#' @param cpu.num The number of cpus if using parallel
#' @param ran.num The number of data points generated from the beta-Poisson model to approximate the theoretical model
#' @param E.esp An small value added to expected value when computing X2, E.esp=0.0 by default
#' @param tbreak.num Number of breaks for binning
#' @param useDebug A parameter setting of \code{\link{getBPMCnull}} function that is just used for debug and checking, so useDebug=FALSE by default
#' @return A list of Monte-Carlo null distributions from the input models (MCdis.list) and setting values of parameters sim.num, ran.num and E.esp
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @examples
#' set.seed(2015)
#' #create random data matrix from a beta-poisson model
#' N=10
#' alp=sample(100,N,replace=TRUE)*0.1;
#' bet=sample(100,N,replace=TRUE)*0.1;
#' lam1=sample(100,N,replace=TRUE)*10;
#' lam2=sample(100,N,replace=TRUE)*0.01;
#' n=100
#' bp.mat=NULL
#' for (i in 1:N)
#'   bp.mat=rbind(bp.mat,rBP(n,alp=alp[i],bet=bet[i],lam1=lam1[i],lam2=lam2[i]))
#' #Estimate parameters from the data set
#' mat.res=estimateBPMatrix(bp.mat,para.num=4,fout=NULL,estIntPar=FALSE,useParallel=FALSE)
#' MCnullmatrix.res=getBPMCnullmatrix(bp.model.list=mat.res$bp.model.list,fout=NULL,
#'                                    sim.num=100,useParallel=FALSE)
#' #Get Monte-Carlo p-values
#' MC.pval=getMCpval(bp.model.list=mat.res$bp.model.list,
#'                   MCdis.list=MCnullmatrix.res$MCdis.list)
#' MC.pval
getBPMCnullmatrix<-function(bp.model.list,fout=NULL,sim.num=1000,useParallel=FALSE,cpu.num=16,ran.num=100000,E.esp=0.0,tbreak.num=10,useDebug=FALSE){
  if (useParallel) registerDoParallel(cores=cpu.num)
  pval.model.ids=which(names(bp.model.list)=="PVAL")
  param.model.ids=which(names(bp.model.list)=="par")
  n.model.ids=which(names(bp.model.list)=="n")
  tbreak.model.ids=which(names(bp.model.list)=="tbreak")
  MCdis.list=NULL  
  if (!useParallel){
    #serial version
    for (i in 1:length(param.model.ids)){
      param=bp.model.list[param.model.ids[i]]$par
      PVAL=bp.model.list[pval.model.ids[i]]$PVAL
      if (!is.na(PVAL)){
        n=bp.model.list[n.model.ids[i]]$n
        tbreak=bp.model.list[tbreak.model.ids[i]]$tbreak
        MC.dis=getBPMCnull(param,n=n,tbreak=tbreak,tbreak.num=tbreak.num,sim.num=sim.num,ran.num=ran.num,E.esp=E.esp,useDebug=useDebug)
      } else {    
        MC.dis=list(PVAL=NA,X2=NA,df=NA,param=param)
      }
      MCdis.list[[i]]=MC.dis
    }
    MCdis.list=unlist(MCdis.list,recursive=FALSE)
  }else{
    ### paralell version    
    MCdis.list=foreach(i=1:length(param.model.ids),.combine=c) %dopar% {
      param=bp.model.list[param.model.ids[i]]$par
      MC.dis=list(PVAL=NA,X2=NA,df=NA,param=param)
      tryCatch({      
      PVAL=bp.model.list[pval.model.ids[i]]$PVAL
      if (!is.na(PVAL)){
        n=bp.model.list[n.model.ids[i]]$n
        tbreak=bp.model.list[tbreak.model.ids[i]]$tbreak
        MC.dis=getBPMCnull(param,n=n,tbreak=tbreak,tbreak.num=tbreak.num,sim.num=sim.num,ran.num=ran.num,E.esp=E.esp,useDebug=useDebug)
      } else {    
        MC.dis=list(PVAL=NA,X2=NA,df=NA,param=param)
      }
      }, error=function(e){
        cat("ERROR at ",i, " : ",conditionMessage(e), "\n")
        })
      MC.dis
    }
  }
  if (!is.null(fout)) save(MCdis.list,sim.num,ran.num,E.esp, file=fout)
  return (list(MCdis.list=MCdis.list,sim.num=sim.num,ran.num=ran.num,E.esp=E.esp))
}


