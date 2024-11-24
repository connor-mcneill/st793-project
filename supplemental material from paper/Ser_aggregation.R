
library(parallel)

name="aggregation"      ## give name for the generated file 

###############################################################

workerFunc <- function(core){
  system(paste("Rscript --vanilla aggregation.R",core))  
}
combinX = list()
no_cores <- 50  ## number of cores

for(i in 1:no_cores){
  combinX[[i]]<- i}

c1 <- makeCluster(no_cores)
parLapply(c1,combinX,workerFunc)
stopCluster(c1)

########## output  combination ##########
T=1000

true_m=array(0,dim=c(4,4,3))

for(i in 1:no_cores){  
    true_m_tempt=NULL
    filename = paste0("true_m",i,".txt")
    true_m_tempt = readRDS(filename ) 
    true_m = true_m + true_m_tempt
}

true_pro=true_m/T

saveRDS(true_pro,file = paste0(name, ".txt"))
save.image(paste0(name, ".RData"))
  