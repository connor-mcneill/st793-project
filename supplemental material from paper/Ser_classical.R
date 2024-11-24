library(parallel)

name="classical"      ## give name for the file 

###############################################################

workerFunc <- function(core){
  system(paste("Rscript --vanilla classical.R",core)) 
}
combinX = list()
no_cores <- 50  ## number of cores

for(i in 1:no_cores){
  combinX[[i]]<- i}

c1 <- makeCluster(no_cores)
parLapply(c1,combinX,workerFunc)
stopCluster(c1)

########## output  combination ##########
T=40000   ## Total simulation times  
length=T/no_cores   

p_seq=c(1,5,10,50,100,200)
true_null=array(NA,dim=c(length(p_seq),3,T))

for(i in 1:no_cores){  
    true_null_tempt=NULL
    filename <- paste0("true_null",i,".txt")
    true_null_tempt<- readRDS(filename) 
    true_null[,,(i*length-length+1):(i*length)]=true_null_tempt

}

save.image(paste0(name, ".RData"))
  