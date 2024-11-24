
library(parallel)

name="model5_B5"      ## give name for the file:
                      ##   the first term specifies the model, and the second term specifies the type of Sigma by the number of block 
                      ##   eg: model5_B5, model5_B100, model4_B5 ....
###############################################################
workerFunc <- function(core){
  system(paste("Rscript --vanilla model_power.R",core))  
}

combinX = list()
no_cores <- 50  ## number of cores

for(i in 1:no_cores){
  combinX[[i]]<- i}

c1 <- makeCluster(no_cores)
parLapply(c1,combinX,workerFunc)
stopCluster(c1)

########## output  combination ##########
T=1000   ## Total simulation times  
length=T/no_cores   

omega_seq=seq(from=0,to=2,by=0.05)
true_array=array(0,dim=c(2,2,length(omega_seq)))

for(i in 1:no_cores){  
    true_array.tempt=NULL
    filename <- paste0("true_array",i,".txt")
    true_array.tempt <- readRDS(filename) 
    true_array=true_array+true_array.tempt
  }

true_pro_array= true_array/T
save.image(paste0(name, ".RData"))
  