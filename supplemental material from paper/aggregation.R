## This is the code for Appendix D.3 of the Supplementary material 

## Parallel computing
#######################################################################

args <- commandArgs(trailingOnly=TRUE)

core <- as.numeric(args[1])

#######################################################################

Q_agg_P<-function(test_seqs_pre, rho, gamma){
    deno = sqrt(2/(n*rho*(1-rho)))
    p_seqs = 1-pnorm(((test_seqs_pre-1)/deno))
    pg_seqs =sort(p_seqs)/gamma
    p_cdf = 1:length(pg_seqs)/length(pg_seqs)
    Q=NULL
    Q=min(1,pg_seqs[which(p_cdf >= gamma)[1]])
    return(Q)
}

Qadap_agg_P<-function(test_seqs_pre, rho, gammamin=0.05){
    gamma_seq=seq(gammamin,1,length=50)[-50]
    Q_seq=rep(NA,length(gamma_seq))
    for(i in 1:length(gamma_seq)){
        Q_seq[i]=Q_agg_P(test_seqs_pre, rho, gamma_seq[i])
    }
    P=NULL
    P=min(1,((1-log(gammamin))*min(Q_seq)))
    return(P)
}

#### Input:
##  X:  n x p matrix of covariates.
##  Y:  n dimensional vector.
##  rho: projection ratio.
##  D: random projection times.
##  S: parameter of sparse projection matrix
## mixRP: whether generate mix of sparse and normal random projection or not

RP.test<-function(X,Y,rho,D=100,S=400,mixRP=FALSE){
  
  n = nrow(X)
  p = ncol(X)
  k = ceiling(rho*n)
  
  I.P1 = diag(n)-(1/n)*matrix(1,n,n)
  pro_pre=1/sqrt(S)
  P.array = array(NA,dim = c(k,p,D))
 
  if(mixRP==FALSE){
      P.array = (1/sqrt(p))*array(rnorm((k*p*D),0,1),dim=c(k,p,D))
  }else{
      for(j in 1:(D/2)){
          data.p=NULL
          data.p=sample(c(-1,0,1),(k*p),replace = T,prob =c((1/2)*pro_pre,(1-pro_pre),(1/2)*pro_pre))
          P.array[,,(2*j-1)]=(1/sqrt(p))*matrix(rnorm(k*p), k,p)
          P.array[,, 2*j]= (1/sqrt(p))*matrix(data.p, k,p)
  }}
  
  H_seq = array(NA,dim = c(n,n,D))
  T_single = rep(NA,D)
  
  for (i in 1:D){
    U = NULL
    U = P.array[,,i]%*%t(X)%*% I.P1
    H = NULL
    H = t(U)%*%solve((U%*%t(U)))%*%U
    H_seq[,,i] = H
    T_single[i] = ((t(Y)%*%H%*%Y)/k)/((t(Y)%*%(I.P1-H)%*%Y)/(n-k-1))
  }
  
  T.ave_hat_seq=rep(NA,4)
  T.ave_all_seq=rep(NA,4)
  T.quantile_seq=rep(NA,4)
  T.quantile_seq_fix=rep(NA,4)
  
  number=c(10,20,50,100)
  
  for (i in 1:4){
    H_ave=NULL
    H_ave=as.matrix(apply(H_seq[,,1:number[i]],c(1,2),mean))
    T.ave_hat_seq[i]= ((t(Y)%*%H_ave%*%Y)/k)/((t(Y)%*%(I.P1-H_ave)%*%Y)/(n-k-1))
    T.ave_all_seq[i]=mean(T_single[1:number[i]])
    T.quantile_seq[i]=Qadap_agg_P(T_single[1:number[i]], rho, gammamin=0.05)
    T.quantile_seq_fix[i]=Q_agg_P(T_single[1:number[i]], rho, gamma=0.5)
  }
  
  deno = sqrt(2/(n*rho*(1-rho)))
  
  T.ave_hat_seq_out = (T.ave_hat_seq-1)/deno
  T.ave_all_seq_out = (T.ave_all_seq-1)/deno
  
  return(list(T.ave_hat_seq_out=T.ave_hat_seq_out, 
              T.ave_all_seq_out=T.ave_all_seq_out,
              T.quantile_seq=T.quantile_seq,
              T.quantile_seq_fix=T.quantile_seq_fix))
  
}

## =============================================================================

#### Input:
##  n, p, m : data dimensions.
##  beta 
##  Gamma: required by generating X.
##  data.gen:  for data generation.
##  rho: projection ratio.

simulation<-function(n,p,m,beta,Gamma,data.gen,rho,sign=0.95,index){
    k=ceiling(rho*n)
    count.RP=matrix(0,4,4)
    L=length(index)
    record.RP = matrix(NA,4,L) 
     
    for(i in 1:L){
        set.seed(index[i])
        data=NULL
        data=data.gen(n,m,p,Gamma,beta)
        X=as.matrix(data$X)
        Y=as.vector(data$Y)
        
        RP.value_list = RP.test(X,Y,rho,D=100,S=400,mixRP=FALSE)
        RP_ave_hat_seq_out = rbind(as.vector(RP.value_list$T.ave_hat_seq_out), 
                                       as.vector(RP.value_list$T.ave_all_seq_out),
                                       as.vector(RP.value_list$T.quantile_seq),
                                       as.vector(RP.value_list$T.quantile_seq_fix))
        for (k in 1:2){
                for (j in 1:4){
                    if ( RP_ave_hat_seq_out[k,j]> qnorm(sign,0,1) ){
                        count.RP[k,j] = count.RP[k,j] + 1
                    }
                }
            }
            for (k in 3:4){
                for (j in 1:4){
                    if ( RP_ave_hat_seq_out[k,j]<= (1-sign) ){
                        count.RP[k,j] = count.RP[k,j] + 1
                    }}}
            
            if (sum(beta^2)==0){
                record.RP[,i]=c(RP_ave_hat_seq_out[1,1],
                                RP_ave_hat_seq_out[2,1],
                                RP_ave_hat_seq_out[3,1],
                                RP_ave_hat_seq_out[4,1])
            }
        }
       
    count=list(count.RP=count.RP)
    return(list(count=count,record.RP=record.RP))
}

## Generate data following the logistic model --------

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

data.bino<-function(n,m,p,Gamma,beta){
    
    Z1 = NULL;   Z1 = matrix(sqrt(3/5)*rt((n*(m/2)),5,ncp=0), n, (m/2))  
    Z2 = NULL;   Z2 = matrix(((rbinom((n*(m/2)), 1, prob=0.5)-0.5)*2), n, (m/2))
    Z = NULL;    Z=cbind(Z1,Z2)
    X = NULL;   X = Z %*% Gamma 
    pro = NULL; pro = logi(X %*% beta)
    Y = NULL;   Y = rbinom(n, 1, prob=pro)
    
    return(list(X=X,Y=Y))
}

setting<-function(n,p,sp_small,omega2,seed){
    
    Gamma=NULL;m=NULL;beta_sp=NULL;beta_ran=NULL
    set.seed(seed)
    s=ceiling(n^1)
    M=ceiling(n^0.8)
    diag.t=rep(0,(p-M))
    v=4
    for (i in 1:(p-M)){
        diag.t[i]=i^(-v)
    }
    sum.v=sum(diag.t)
    dia.tempt=((s-M)/sum.v)*diag.t
    dia=sqrt(c(rep(1,M),dia.tempt))
    dia=sqrt(dia)
    
    save=rnorm(p^2,0,1)
    orignal=matrix(save,p,p)
    symtric=orignal+t(orignal)        
    O=as.matrix(eigen(symtric)$vectors)  
    
    Gamma=NULL
    Gamma=O%*%diag(dia)%*%t(O) 
    Sigma= Gamma%*%t(Gamma)
    m=p
    
    beta_sp_tempt = rep(0,p)
    beta_sp_tempt[sample(1:p,sp_small)] =  rep(1, sp_small)
     
    beta_ran_tempt=O[,1:100]%*%rnorm(100)
    
    t1=as.numeric(t(beta_sp_tempt)%*%Sigma%*%beta_sp_tempt)
    t2=as.numeric(t(beta_ran_tempt)%*%Sigma%*%beta_ran_tempt)
    
    beta_sp=sqrt(omega2/t1)* beta_sp_tempt
    beta_ran = sqrt(omega2/t2)*beta_ran_tempt
    
    return(list(Gamma=Gamma,m=m,beta_sp=beta_sp,beta_ran=beta_ran))
}


n=400;p=1000 
setting_tempt=setting(n,p,sp_small=10, omega2=0.4,seed=20)  
## omega2=b^2 takes value 0.4 or 0.8

Gamma=as.matrix(setting_tempt$Gamma)
m=as.numeric(setting_tempt$m)
beta_sp=as.vector(setting_tempt$beta_sp)
beta_ran=as.vector(setting_tempt$beta_ran)

index=(20*core-19):(20*core)  
true_m=array(NA,dim=c(4,4,3))

simu_null=simulation(n,p,m,beta=rep(0,p),Gamma,data.gen=data.bino,rho=0.4,sign=0.95,index)
simu_sp=simulation(n,p,m,beta=beta_sp,Gamma,data.gen=data.bino,rho=0.4,sign=0.95, index)
simu_ran=simulation(n,p,m,beta=beta_ran,Gamma,data.gen=data.bino,rho=0.4,sign=0.95,index)

true_m[,,1]=as.matrix(simu_null$count$count.RP)
true_m[,,2]=as.matrix(simu_sp$count$count.RP)
true_m[,,3]=as.matrix(simu_ran$count$count.RP)

filename <- paste0("true_m",core,".txt")

saveRDS(true_m, 
        file = filename)


 
