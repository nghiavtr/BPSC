#' Estimate parameters of beta-Poisson models for a data matrix
#'
#' @param dataMat Data matrix that needs to be modelled
#' @param para.num Mode of beta-Poisson model: 3, 4 (default) or 5 parameters
#' @param tbreak.num Number of breaks for binning
#' @param fout A *.RData file name to export results
#' @param break.thres A parameter setting of \code{\link{getTbreak}} function
#' @param estIntPar An option to allow estimating initial parameters for the model from only expressed values
#' @param extreme.quant A quantile probability to remove extrem values (outliers) higher than the quantile. If extreme.quant=NULL, no elimination of outliers is done
#' @param useExt A parameter setting of \code{\link{getTbreak}} function that allows to extend the last bin to infinity or not
#' @param min.exp A threshold for minimum expressed values. If a expression is less than min.exp, it is set to be zero
#' @param useParallel An option to allow using parallel computing
#' @return A list of optimal models corresponding to the rows of the matrix. Each model consists of optimal parameters, X2 test results (X2 and PVAL), etc..
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
#' #In this function, user can also set estIntPar=TRUE to have better estimated beta-Poisson 
#' #models for the generalized linear model. However, a longer computational time is required.
estimateBPMatrix<-function(dataMat,para.num=4,tbreak.num=10,fout=NULL,break.thres=10,estIntPar=TRUE,extreme.quant=NULL,useExt=FALSE,min.exp=1e-4,useParallel=FALSE){
  ind.set=NULL  
  bp.model.list=NULL

  if (useParallel){
    bp.model.list=foreach(i=1:nrow(dataMat),.combine=c) %dopar% {
      obpOk=NULL;
      ind=NULL;
      bp.sub.count=dataMat[i,]
      bpstat=bp.sub.count
      if (para.num==5) bpstat=c(bp.sub.count,0) # add one zero to data for the fifth parameter
      bpstat[bpstat<min.exp]=0 # make sure expression not too small
      if (!is.null(extreme.quant)){
        bpstat=bpstat[bpstat <= quantile(bpstat,prob=extreme.quant)]
      }
      if (sum(bpstat>0) > 0.05*length(bpstat)){                
        obp=estimateBP(bpstat,para.num=para.num,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt)

        if (estIntPar){
          bpstat2=bpstat[bpstat>0]

          if (para.num == 5){
            mypar2=getInitParam(bpstat2,para.num=4)
            oo2=estimateBP(bpstat2,para.num=4,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,param0=mypar2)
            mypar=c(oo2$par,0)
            if (para.num==5) mypar[5]=sum(bpstat==0)/length(bpstat)
          } else {            
            mypar2=getInitParam(bpstat2,para.num=para.num)
            oo2=estimateBP(bpstat2,para.num=para.num,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,param0=mypar2)
            mypar=oo2$par            
          }

          oo=estimateBP(bpstat,para.num=para.num,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,param0=mypar)
          # select the better model
          if (is.na(obp$PVAL)) obp=oo
          else {
              if (!is.na(oo$PVAL) && oo$PVAL > obp$PVAL && oo$X2>=0) obp=oo              
          }
        }
          
        # store parameters and model in lists if the results are valid
        if (obp$X2>=0 && !is.na(obp$PVAL)){
#          ind=i
          obpOk=obp
          obpOk[["ind"]]=i
        }
      }
      obpOk
    }
    ind.set=unlist(bp.model.list[which(names(bp.model.list)=="ind")])
    names(ind.set)=NULL
    
    } else{ 
    for (i in 1:nrow(dataMat)){
      bp.sub.count=dataMat[i,]    
      bpstat=bp.sub.count
      if (para.num==5) bpstat=c(bp.sub.count,0) # add one zero to data for the fifth parameter
      bpstat[bpstat<min.exp]=0 # make sure expression not too small
      if (!is.null(extreme.quant)){
        bpstat=bpstat[bpstat <= quantile(bpstat,prob=extreme.quant)]
      }
      if (sum(bpstat>0) > 0.05*length(bpstat)){
        obp=estimateBP(bpstat,para.num=para.num,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt)
        
        if (estIntPar){
          bpstat2=bpstat[bpstat>0]

          if (para.num == 5){
            mypar2=getInitParam(bpstat2,para.num=4)
            oo2=estimateBP(bpstat2,para.num=4,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,param0=mypar2)
            mypar=c(oo2$par,0)
            if (para.num==5) mypar[5]=sum(bpstat==0)/length(bpstat)
          } else {            
            mypar2=getInitParam(bpstat2,para.num=para.num)
            oo2=estimateBP(bpstat2,para.num=para.num,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,param0=mypar2)
            mypar=oo2$par            
          }
          

          oo=estimateBP(bpstat,para.num=para.num,tbreak.num=tbreak.num,break.thres=break.thres,useExt=useExt,param0=mypar)
          # select the better model
          if (is.na(obp$PVAL)) obp=oo
          else {
              if (!is.na(oo$PVAL) && oo$PVAL > obp$PVAL && oo$X2>=0) obp=oo              
          }
        }
        
        # store parameters and model in lists if the results are valid
        if (obp$X2>=0 && !is.na(obp$PVAL)){
          ind.set=c(ind.set,i)
          obp[["ind"]]=i
          bp.model.list=c(bp.model.list,obp)
        }
      }
    }
  }
  if (!is.null(fout)) save(ind.set,tbreak.num,bp.model.list,break.thres,para.num,file=fout)
  return(list(ind.set=ind.set,tbreak.num=tbreak.num,bp.model.list=bp.model.list,break.thres=break.thres,para.num=para.num,estIntPar=estIntPar,extreme.quant=extreme.quant,useExt=useExt))
}
