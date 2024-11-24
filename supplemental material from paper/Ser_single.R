library(parallel)

name="Single"      ## give name for the generated file
###############################################################
workerFunc <- function(core){
  system(paste("Rscript --vanilla single.R",core))  
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

no_seq=c(1,100,5)
cor_array=array(0,dim = c(3,3,length(no_seq)))
for(i in 1:no_cores){  
    cor_array_tempt = NULL
    filename = paste0("matrix_cor",i,".txt")
    cor_array_tempt = readRDS(filename)
    cor_array  = cor_array+ cor_array_tempt
}
pro_array=cor_array/T

saveRDS(pro_array, file=paste0(name,".txt"))
save.image(paste0(name, ".RData"))
  