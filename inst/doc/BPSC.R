## ----include = FALSE-----------------------------------------------------
library(BPSC)

## ----include = TRUE------------------------------------------------------
hist(log2(rBP(10000,alp=0.2160247,bet=0.5989889,lam1=171.9596728,lam2=0.9094114)+1),prob=TRUE,col='grey',breaks=10, xlab="log2(expression+1)", main="")

## ----include = TRUE------------------------------------------------------
library(BPSC)
set.seed(2015)
###Generate a random data matrix from a beta-poisson model
#Set the number of genes
N=100
#Generate randomly the parameters of BP models
alp=sample(100,N,replace=TRUE)*0.1;
bet=sample(100,N,replace=TRUE)*0.1;
lam1=sample(100,N,replace=TRUE)*10;
lam2=sample(100,N,replace=TRUE)*0.01

#Generate a control group
n1=100
control.mat=NULL
for (i in 1:N) control.mat=rbind(control.mat,rBP(n1,alp=alp[i],
	bet=bet[i],lam1=lam1[i],lam2=lam2[i]))

#To create biological variation, we randomly set 10% as differentially expressed genes 
#by simply replacing the parameter lam1 in treated group by a fold-change fc
DE.ids= sample(N,N*0.1)
fc=2.0
lam1[DE.ids]=lam1[DE.ids] * fc

#Generate a treated group
n2=100
treated.mat=NULL
for (i in 1:N)treated.mat=rbind(treated.mat,rBP(n2,alp=alp[i],
	bet=bet[i],lam1=lam1[i],lam2=lam2[i]))

#Create a data set by merging the control group and the treated group
bp.mat=cbind(control.mat,treated.mat)
rownames(bp.mat)=c(1:nrow(bp.mat));
colnames(bp.mat)=c(1:ncol(bp.mat))
group=c(rep(1,ncol(control.mat)),rep(2,ncol(treated.mat)))

## ----include = TRUE------------------------------------------------------
#First, choose IDs of all cells of the control group for estimating parameters of BP models
controlIds=which(group==1)

#Create a design matrix including the group labels. All batch effects can be also added here if they are available
design=model.matrix(~group) 
#Select the column in the design matrix corresponding to the coefficient (the group label) for the GLM model testing
coef=2 

#Run BPglm for differential expression analysis
res=BPglm(data=bp.mat, controlIds=controlIds, design=design, coef=coef, estIntPar=FALSE, useParallel=FALSE) 
#In this function, user can also set estIntPar=TRUE to have better estimated beta-Poisson models for the generalized linear model. However, a longer computational time is required.

#Plot the p-value distribution
hist(res$PVAL, breaks=20)
#Summarize the resutls
ss=summary(res)
#Compare the discovered DE genes and the true DE genes predefined beforeward
fdr=p.adjust(res$PVAL, method="BH")
bpglm.DE.ids=which(fdr<=0.05)
#Print the indices of the true DE genes:
cat(sort(DE.ids))
#Print the indices of the DE genes discovered by BPglm:
cat(sort(bpglm.DE.ids))

## ----include = TRUE------------------------------------------------------
library(doParallel)
registerDoParallel(cores=16)

## ----include = TRUE------------------------------------------------------
library(BPSC)
set.seed(2015)
#Create a simulated gene expression by randomly generating 100 data points 
#from a beta-poisson model
alp=0.6;bet=1.5;lam1=20;lam2=0.05
par0=c(alp,bet,lam1,lam2)
bp.vec=rBP(100,par0)

#Estimate parameters of the four-parameter beta-Poisson model from the data points
res=estimateBP(bp.vec,para.num=4)
#Print the goodness-of-fit of the model and the optimal parameters
res$X2
res$par
#Generate Monte-Carlo null distribution of the model.
#Due to time limit, the number of simulations (sim.num) here is set by 100.
MCnull.res=getBPMCnull(res$par,n=100,tbreak=res$tbreak,sim.num=100)
#Compute Monte-Carlo p-value
MCpval=sum(MCnull.res$X2 >= res$X2)/length(MCnull.res$X2)
MCpval
#Plot the Monte-Carlo null distribution and the goodness-of-fit from the model (blue line).
hist(MCnull.res$X2,xlab="goodness-of-fit statistic", 
main="Monte-Carlo null distribution",breaks=20, prob=TRUE)
lines(c(res$X2,res$X2),c(0,par("usr")[4]),
      col="blue",lwd=2.0)

## ----include = TRUE------------------------------------------------------
set.seed(2015)
#create random data matrix from a beta-poisson model
N=10
alp=sample(100,N,replace=TRUE)*0.1;
bet=sample(100,N,replace=TRUE)*0.1;
lam1=sample(100,N,replace=TRUE)*10;
lam2=sample(100,N,replace=TRUE)*0.01;
n=100
bp.mat=NULL
for (i in 1:N)
  bp.mat=rbind(bp.mat,rBP(n,alp=alp[i],bet=bet[i],lam1=lam1[i],lam2=lam2[i]))

## ----include = TRUE------------------------------------------------------
#Estimate parameters from the data set
mat.res=estimateBPMatrix(bp.mat,para.num=4,fout=NULL,estIntPar=FALSE,useParallel=FALSE)

## ----include = TRUE------------------------------------------------------
MCnullmatrix.res=getBPMCnullmatrix(bp.model.list=mat.res$bp.model.list,fout=NULL,
                                   sim.num=100,useParallel=FALSE)
#Get Monte-Carlo p-values
MCpval=getMCpval(bp.model.list=mat.res$bp.model.list,
                  MCdis.list=MCnullmatrix.res$MCdis.list)
MCpval

