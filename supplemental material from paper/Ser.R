
library(parallel)

name="n400p1000logiSigma2"      ## give name for the file 
###############################################################
workerFunc <- function(core){
  system(paste("Rscript --vanilla Main_code.R",core))
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

true_matrix=matrix(0,3,3)
true_record=rep(0,T)
sign_seq=seq(from=0.01,to=0.99,by=0.04)
true_array=array(0,dim=c(3,3,length(sign_seq)))

for(i in 1:no_cores){  
  true_matrix.tempt=NULL
  true_record.tempt=NULL
  true_array.tempt=NULL
  filename1 <- paste0("true_m",i,".txt")
  true_matrix.tempt <- as.matrix(read.table(filename1))
  true_matrix=true_matrix+true_matrix.tempt
  filename2 <- paste0("true_null",i,".txt")
  true_record.tempt<- as.matrix(read.table(filename2))
  true_record[(i*length-length+1):(i*length)]=true_record.tempt
  ## uncomment the following if one wants to generate figure such as Figure 1 in Section 5.1
  # filename3 <- paste0("true_array",i,".txt")  
  # true_array.tempt <- readRDS(filename3)
  # true_array=true_array+true_array.tempt
}

true.pro=true_matrix/T
## uncomment the following if one wants to generate figure such as Figure 1 in Section 5.1
# true.pro.array=true_array/T    
save.image(paste0(name, ".RData"))
  