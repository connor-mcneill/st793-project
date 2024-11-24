library(parallel)

name="Sub_p2"      ## give name to the generated file 
###############################################################
workerFunc <- function(core){
  system(paste("Rscript --vanilla sub_p2.R",core))  
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
cor_array=array(0,dim = c(2,3,3))
 
for(i in 1:no_cores){  
    cor_array_tempt=NULL
    filename<- paste0("matrix_cor",i,".txt")
    cor_array_tempt <- readRDS(filename)
    cor_array= cor_array+cor_array_tempt
}

cor_pro=cor_array/T

saveRDS(cor_pro, file=paste0(name,".txt"))
save.image(paste0(name, ".RData"))
  