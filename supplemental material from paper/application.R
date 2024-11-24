
########### The tests ##########################################

library(glmnet)

## MCL test
#### Input:
# X: n x p matrix of covariates.
# Y: n dimensional vector.
#  nfolds: number of folds used for cross-validation of Lasso/logistic Lasso, default = 5.
#  lambda: tuning parameter for the logistic lasso, default = NA. If specified, nfolds will be overwritten.
#  tune.1, tune.2 are tuning parameters in the nodewise lasso. Should be values between 0 and 1.
#  fdr: false discovery rate for multiple testing.
#  please refer to Ma, Cai, and Li (2020) for details.
logistic.test <- function(X,Y,nfolds=5, lambda = 0, 
                          tune.1 = 1.5, tune.2 = 1, 
                          intercept = F, fdr = 0.05){
  f = function(x){
    exp(x)/(1+exp(x))
  }
  p = dim(X)[2]
  n = dim(X)[1]
  
  if(lambda == 0){
    logistic.cv = cv.glmnet(x = X, y = Y, family = "binomial", alpha = 1,  
                            intercept=intercept, nfolds = nfolds, type.measure = "class")
    lambda = logistic.cv$lambda.min
  }
  
  my.logistic.fit = glmnet(x = X, y = Y, family = "binomial", alpha = 1,  
                           intercept=intercept, lambda = lambda, standardize = F)
  
  b.hat = coef(my.logistic.fit)
  
  W.n1 = c(exp(X%*% as.matrix(b.hat)[-1,])/(1+exp(X%*% as.matrix(b.hat)[-1,]))^2)^(-1)
  zeta.try = matrix(nrow =5,ncol = p)
  tau.try = matrix(nrow =5,ncol = p)
  
  V = matrix(ncol=n, nrow = p)
  tau = c()
  for(i in 1:p){
    nodewise.try = glmnet(x= X[,-i], y = X[,i], family = "gaussian", alpha = 1, intercept = F, nlambda = 5, standardize = F)
    for(lambda.i in 1:5){
      V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.i]
      zeta.try[lambda.i,i] = max(abs(as.vector(V[i,]%*%X[,-i])/sqrt(sum((V[i,])^2*W.n1))))
      tau.try[lambda.i,i] = sqrt(sum((V[i,])^2*W.n1))/(V[i,]%*% X[,i])
    }
    zeta0 = sqrt(2*log(p))
    if(min(zeta.try[,i])>sqrt(2*log(p))) zeta0 = tune.1*min(zeta.try[,i])
    lambda.chosen = order(nodewise.try$lambda[which(zeta.try[,i] <= zeta0)], decreasing=T)[1]
    tau[i] = tau.try[lambda.chosen,i]
    lambda.chosen = order(nodewise.try$lambda[tau.try[,i]<=tune.2*tau[i]],decreasing = F)[1]
    tau[i] = tau.try[lambda.chosen,i]
    V[i,] = X[,i]-X[,-i]%*%nodewise.try$beta[,lambda.chosen]
  }
  
  V2 = t((t(V)*W.n1))
  b.check = c()
  for(j in 1:p){
    b.check[j] = b.hat[j+1]+(V2[j,]%*%(Y-f(X %*% as.matrix(b.hat)[-1])))/(V[j,] %*% X[,j])
  }
  test.value = max((b.check/tau)^2)
  
  return(test.value)
}

logi<-function(x){
  y = NULL
  y = exp(x)/(1+exp(x))
  return(y) 
}

logi.deri<-function(x){
  y = NULL
  y = exp(x)/((1+exp(x))^2)
  return(y)
}

 
## GC test
#### Input:
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  g.assume, g.assume.deri:  model assumption for mean and its first derivative.
##  V.assume:  model assumption for variance.
##  please refer to Guo and Chen (2016) for details.
CSX.test<-function(X,Y, g.assume, g.assume.deri, V.assume){
  
  n = nrow(X)
  p = ncol(X)
  mu.0 = rep(g.assume(0), n)  
  var.ratio=g.assume.deri(0)/ V.assume(0)
  Phi.0 = rep(var.ratio, n)            
  Un.tempt = 0
  tr.Var.tempt = 0
  
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      tempt.value = NULL
      tempt.value = (Y[i]- mu.0[i]) * (Y[j]- mu.0[j]) * Phi.0[j] * Phi.0[i] * crossprod(X[i,],X[j,])
      Un.tempt = Un.tempt + tempt.value
      tr.Var.tempt = tr.Var.tempt + tempt.value^2
    }
  }
  Un = NULL
  Un = (1/n) * 2 * Un.tempt
  tr.Var = NULL
  tr.Var = (1/(n*(n-1))) * 2 * tr.Var.tempt
  CSX.test = Un/sqrt(2*tr.Var) 
  
  return (CSX.test)
}


## RP test
#### Input:
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  rho: projection ratio.
##  D: random projection times.
##  S: parameter of sparse projection matrix
RP.test<-function(X,Y,rho,D=10,S=400){
  
  n = nrow(X)
  p = ncol(X)
  k = ceiling(rho*n)
  I.P1 = diag(n) - (1/n)*matrix(1,n,n)
  P.array=array(NA,dim = c(k,p,D))
  
  for(i in 1:D){
    data.p=NULL
    data.p=rnorm((k*p),0,1)
    P.tempt=NULL
    P.tempt=matrix(data.p,nrow=k,ncol=p)
    P.array[,,i]=P.tempt
  }
  
  P=NULL
  P=as.matrix(apply(P.array,c(1,2),mean))
  P.new=NULL
  P.new=as.matrix(P.array[,,1])
  
  prob=1/sqrt(S)
  data.p=NULL
  data.p=sample(c(-1,0,1),(k*p),replace = T,prob =c((1/2)*prob,(1-prob),(1/2)*prob) )
  
  P_sparse=NULL
  P_sparse=matrix(data.p,nrow=k,ncol=p)

  X.P.new=X%*%t(P.new)
  Z.new=I.P1%*%X.P.new
  PX.new=Z.new%*%solve(t(Z.new)%*%Z.new)%*%t(Z.new)
  T.new=((t(Y)%*%(PX.new)%*%Y)/k)/((t(Y)%*%(I.P1-PX.new)%*%Y)/(n-k-1))
  
  X.P.multi=X%*%t(P)
  Z.multi=I.P1%*%X.P.multi
  PX.multi=Z.multi%*%solve(t(Z.multi)%*%Z.multi)%*%t(Z.multi)
  T.multi=((t(Y)%*%PX.multi%*%Y)/k)/((t(Y)%*%(I.P1-PX.multi)%*%Y)/(n-k-1))
  
  X.P.sparse=X%*%t(P_sparse)
  Z.sparse=I.P1%*%X.P.sparse
  PX.sparse=Z.sparse%*%solve(t(Z.sparse)%*%Z.sparse)%*%t(Z.sparse)
  T.sparse=((t(Y)%*%PX.sparse%*%Y)/k)/((t(Y)%*%(I.P1-PX.sparse)%*%Y)/(n-k-1))
  
  deno<-sqrt(2/(n*rho*(1-rho)))
  
  multi<-(T.multi-1)/deno
  new<-(T.new-1)/deno
  sparse<-(T.sparse-1)/deno
  
  return(list(multi=multi,new=new,sparse=sparse))
  
}

## RP test for testing partial regression coefficients
##  n: sample size
##  p1, p2: dimensions of covariates x1 and x2 respectively. Here x^T=(x1^T,x2^T).
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  rho2: projection ratio.
##  D: random projection times.
##  S: parameter of sparse projection matrix

RP.test.sub<-function(n,p1,p2,X,Y,rho2,D,S){
  p=p1+p2
  
  X1=X[,1:p1]
  X2=X[,(p1+1):p]
  
  k2=ceiling(rho2*n)
  rho1=p1/n
  
  P1=(1/n)*matrix(1,n,n)
  I.P1=diag(n)-P1
  
  P.array=array(NA,dim = c(k2,p2,D))
  
  for(i in 1:D){
    data.p=NULL
    data.p=rnorm((k2*p2),0,1)
    P.tempt=NULL
    P.tempt=matrix(data.p,nrow=k2,ncol=p2)
    P.array[,,i]=P.tempt
  }
  
  P=NULL
  P=as.matrix(apply(P.array,c(1,2),mean))
  P.new=NULL
  P.new=as.matrix(P.array[,,1])
  
  prob=1/sqrt(S)
  data.p=NULL
  data.p=sample(c(-1,0,1),(k2*p2),replace = T,prob =c((1/2)*prob,(1-prob),(1/2)*prob) )
  
  P_sparse=NULL
  P_sparse=matrix(data.p,nrow=k2,ncol=p2)
  
  UX1=I.P1%*%X1
  PX1=UX1%*%solve(t(UX1)%*%UX1)%*%t(UX1)
  
  X2.P.new=X2%*%t(P.new)
  Z.new=I.P1%*%cbind(X1,X2.P.new)
  PX.new=Z.new%*%solve(t(Z.new)%*%Z.new)%*%t(Z.new)
  T.new=((t(Y)%*%(PX.new-PX1)%*%Y)/k2)/((t(Y)%*%(I.P1-PX.new)%*%Y)/(n-k2-1-p1))
  
  X2.P.multi=X2%*%t(P)
  Z.multi=I.P1%*%cbind(X1,X2.P.multi)
  PX.multi=Z.multi%*%solve(t(Z.multi)%*%Z.multi)%*%t(Z.multi)
  T.multi=((t(Y)%*%(PX.multi-PX1)%*%Y)/k2)/((t(Y)%*%(I.P1-PX.multi)%*%Y)/(n-k2-1-p1))
  
  X2.P.sparse=X2%*%t(P_sparse)
  Z.sparse=I.P1%*%cbind(X1,X2.P.sparse)
  PX.sparse=Z.sparse%*%solve(t(Z.sparse)%*%Z.sparse)%*%t(Z.sparse)
  
  T.sparse=((t(Y)%*%(PX.sparse-PX1)%*%Y)/k2)/((t(Y)%*%(I.P1-PX.sparse)%*%Y)/(n-k2-1-p1))
  
  deno<-sqrt((2*(1-rho1))/(n*rho2*(1-rho1-rho2)))
  multi<-(T.multi-1)/deno
  new<-(T.new-1)/deno
  sparse<-(T.sparse-1)/deno
  return(list(multi=multi,new=new,sparse=sparse))
  
}

###########Processing data##########################################
library(Biobase)
library(GEOquery)
library(stringr)
 
##---------------Download data--------------------
ori_gset <- getGEO("GSE50948", GSEMatrix =TRUE, getGPL=FALSE)
oi_gp<-getGEO("GPL570")	

save(ori_gset, file = 'ori_gset.Rdata' )
save(oi_gp, file = 'oi_gp.Rdata' )

##---------------load data& prepocessed--------------------
load("ori_gset.Rdata")
load("oi_gp.Rdata")

#ori_gset=read.csv("GSE50948_series_matrix.txt")
if (length(ori_gset) > 1) idx <- grep("GPL570", attr(ori_gset, "names")) else idx <- 1
gset <- ori_gset[[idx]]

data=as.data.frame(exprs(gset))
data=data.frame(rownames(data),data)
colnames(data)[1]="ID_REF"
data[,c((ncol(data)+1):(ncol(data)+2))]=Table(oi_gp)[match(data$ID_REF,Table(oi_gp)$ID),c("ENTREZ_GENE_ID","Gene Symbol")]
colnames(data)[(ncol(data)-1):ncol(data)]=c("GeneID","GeneSymbol")
data = data[!is.na(data[,"GeneID"]),]
data = data[data[,"GeneID"]!="",]
data=data[!grepl(pattern="///",data[,ncol(data)-1]),]

ids=data[,c(1,(ncol(data)-1),(ncol(data)))]
ids$OR=apply(data[,c(2:(ncol(data)-2))],1,IQR)

## each gene corresponding to one signature picked.
ids=ids[order(ids$GeneSymbol,ids$OR,decreasing = T),] 
ids=ids[!duplicated(ids$GeneSymbol),]
processed_data=data[as.vector(ids$ID_REF),]
## delet ID refer, gene ID, gene symbol
final_data=processed_data[,-c(1,(ncol(processed_data)-1),(ncol(processed_data)))]
final_data=as.matrix(sapply(final_data, as.numeric)) 
rownames(final_data)=as.vector(processed_data$GeneSymbol)

##========= clinicopathologically variables=====================================
## delete some patients for incomplete information
final_patient=pData(ori_gset[[1]])[-c(18,57,155,156),]
names_patient=row.names(final_patient)
genes_patient=final_data[,names_patient]

##(1) PCR
PCR=as.character(final_patient$"characteristics_ch1.28")
PCR[which(PCR=="pcr: pCR")]=1
PCR[which(PCR=="pcr: RD")]=0
PCR=as.numeric(PCR)

##(2) age at diagnosis 
age=NULL
age=as.character(final_patient$"characteristics_ch1.24")
library(stringr)
age<-as.numeric((str_replace(age, "age:","")))

##(3) Histological grade
histological_grade=NULL
histological_grade=as.numeric(final_patient$"invasive_tumor_grade:ch1")

##(4) status of ER: ER+:1; ER-:0
ER=NULL
ER=as.character(final_patient$"characteristics_ch1.19")
ER[which(ER=="er: ER+")]=1
ER[which(ER=="er: ER-")]=0
ER=as.numeric(ER)

##(5) status of PR: PR+:1; PR-:0
PR=NULL
PR=as.character(final_patient$"characteristics_ch1.20")
PR[which(PR=="pr: PR+")]=1
PR[which(PR=="pr: PR-")]=0
PR=as.numeric(PR)

##(6) status of HER2
HER2=NULL
HER2=as.character(final_patient$"characteristics_ch1.21")
patient_HER2_plus=which(HER2=="her2: HER2+")
patient_HER2_minus=which(HER2=="her2: HER2-")
HER2[patient_HER2_plus]=1
HER2[patient_HER2_minus]=0
HER2=as.numeric(HER2)

##(7) inflammatory breast cancer: yes:1 no:0
IBC=NULL
IBC=as.character(final_patient$"inflammatory.brca:ch1")
IBC[which(IBC=="yes")]=1
IBC[which(IBC=="no")]=0
IBC=as.numeric(IBC)

##(8) "invasive_tumor_area_size1 [mm]:ch1"
Size1=NULL
Size1=as.numeric(final_patient$"invasive_tumor_area_size1 [mm]:ch1")
 
##(9) "invasive_tumor_area_size2 [mm]:ch1"
Size2=NULL
Size2=as.numeric(final_patient$"invasive_tumor_area_size2 [mm]:ch1")

##(10) treatment
treatment=NULL
treatment=as.character(final_patient$"treatment:ch1" )
treatment[which(treatment=="neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF)")]=0 
treatment[which(treatment=="neoadjuvant doxorubicin/paclitaxel (AT) followed by cyclophosphamide/methotrexate/fluorouracil (CMF) + Trastuzumab for 1 year")]=1
treatment=as.numeric(treatment)

## all the variavles considered
## genes,age,histological grade,ER,PR,IBC,HER2,Size1,Size2,treatment
covatiate=NULL
covatiate=t(rbind(genes_patient,age,histological_grade,ER,PR,IBC,HER2,Size1,Size2,treatment))
names_covariate=as.character(colnames(covatiate))
length(names_covariate)

##======= PAM50 ===================================
PAM_tring="ACTR3B, ANLN, BAG1, BCL2, BIRC5, BLVRA, CCNB1, CCNE1, CDC20, CDC6, CDH3, CENPF, CEP55, CXXC5, EGFR, ERBB2, ESR1, EXO1, FGFR4, FOXA1, FOXC1, GPR160, GRB7, KIF2C, KRT14, KRT17, KRT5, MAPT, MDM2, MELK, MIA, MKI67, MLPH, MMP11, MYBL2, MYC, NAT1, NDC80, NUF2, ORC6L, PGR, PHGDH, PTTG1, RRM2, SFRP1, SLC39A6, TMEM45B, TYMS, UBE2C, UBE2T"
PAM50_gene=str_replace(unlist(strsplit(PAM_tring, ","))," ","")
#test_PAM50_gene=PAM50_gene[!grepl("CEP55|EXO1|GRB7|MIA|ORC6L|UBE2C|UBE2T",PAM50_gene)]

diff_gene=setdiff(PAM50_gene,colnames(covatiate))
test_PAM50_gene=PAM50_gene[-which(PAM50_gene==diff_gene[1])]
test_PAM50_gene=test_PAM50_gene[-which(test_PAM50_gene==diff_gene[2])]
##======= Start test================================

criterion=c("ER","PR","HER2")
sign=0.95
MCL=rep(NA,3)

Gambel<-function(x){
  y=NULL
  y=exp(-(1/sqrt(pi)*exp(-x/2)))
  return(y)
}

for (i in 1:3){
print(paste0("consider ",criterion[i]))
X=NULL
X=covatiate[,!grepl("criterion[i]|treatment",colnames(covatiate))] 
Y=NULL
Y=covatiate[,criterion[i]]

n=nrow(X); p=ncol(X)
criti_Gambel= 2*log(p)-log(log(p))-log(pi)-2*log(log(1/(sign)))

## testing for all the parameters:
##--------------------------------
## RP:  
RP.value = RP.test(X,Y,rho=0.2,D=10,S=400)
RP=as.numeric(RP.value$new) 
multi_RP=as.numeric(RP.value$multi) 
S_RP=as.numeric(RP.value$sparse) 
##(1)RP:
print(paste0("RP:", RP," multi-RP: ",multi_RP," S-RP: ",S_RP))
print(paste0("RP test: ",RP>qnorm(sign,0,1)," p_value: ", 1-pnorm(RP,0,1)))
##(1)multi RP:
print(paste0("multi_RP test: ",multi_RP>qnorm(sign,0,1)," p_value: ",1-pnorm(multi_RP,0,1)))
##(1)S RP:
print(paste0("S_RP test: ",S_RP>qnorm(sign,0,1), " p_value: ",1-pnorm(S_RP,0,1)))

## GC:
CSX.value = CSX.test(X,Y,g.assume=logi,g.assume.deri=logi.deri, V.assume=logi.deri)
print(paste0("CSX test: ",CSX.value>qnorm(sign,0,1), " CSX-p_value: ", 1-pnorm(CSX.value)))

## MCL:
MCL.value = logistic.test(X,Y, nfolds=5, lambda = 0, 
                         tune.1 = 1.5, tune.2 = 1.01, 
                       intercept = F, fdr = 0.05)
MCL[i]=MCL.value
print(paste0("MCL test: ",MCL.value > criti_Gambel, " TC-p-value: ", 1-Gambel(MCL.value-2*log(p)+log(log(p)))))

## testing for partial parameters:
##--------------------------------
delete_name=c("age","histological_grade","ER","PR","HER2","IBC","Size1","Size2", "treatment", test_PAM50_gene)
X1=covatiate[,delete_name]
X2=covatiate[,-which(colnames(covatiate) %in% delete_name)]
p1=ncol(X1);  p2=ncol(X2)

X.RE=cbind(X1,X2)
sub_RP_test=RP.test.sub(n,p1,p2,X.RE,Y,rho2=0.4,D=10,S=400)

sub_RP=as.numeric(sub_RP_test$new) 
sub_multi_RP=as.numeric(sub_RP_test$multi) 
sub_S_RP=as.numeric(sub_RP_test$sparse) 


print(paste0("sub_RP:", sub_RP," sub_multi-RP: ",sub_multi_RP," sub_S-RP: ",sub_S_RP))
print(paste0("sub_RP test: ",sub_RP>qnorm(sign,0,1)," sub_p_value: ", 1-pnorm(sub_RP,0,1)))
##(1)multi RP:
print(paste0("sub_multi_RP test: ",sub_multi_RP>qnorm(sign,0,1)," sub_p_value: ",1-pnorm(sub_multi_RP,0,1)))
##(1)S RP:
print(paste0("sub_S_RP test: ",sub_S_RP>qnorm(sign,0,1), " sub_p_value: ",1-pnorm(sub_S_RP,0,1)))
print(paste0("finish for ", criterion[i]))
}
 