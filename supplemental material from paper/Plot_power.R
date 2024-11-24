library(ggplot2)
library(extrafont)
 
setwd(" ")  ## set working directory to the place where the generated data is saved

omega_seq=seq(from=0,to=2,by=0.05); ll=length(omega_seq)
class=c("B1","B100","B5")

##-- load the generated data--------------------------------

cos_array<-array(NA,dim=c(2,2,length(omega_seq),3))

for(i in 1:3){  
    filename <- NULL
    filename <- paste0("model4_",class[i],".txt")
    cos_array_tempt <- readRDS(filename) 
    cos_array[,,,i]=cos_array_tempt
}

square_array<-array(NA,dim=c(2,2,length(omega_seq),3))

for(i in 1:3){  
    filename <- NULL
    filename <- paste0("model5_",class[i],".txt")
    square_array_tempt <- readRDS(filename) 
    square_array[,,,i]=square_array_tempt
}

## plot-----------------------------------------------------------

cbbPalette <- c("#009E73", "#D55E00")  

df1<-data.frame(
    omega_seq=c(omega_seq,omega_seq,omega_seq,omega_seq,omega_seq,omega_seq),
    cos=c(cos_array[1, 1, ,1],cos_array[1, 1, ,2],cos_array[1, 1, ,3], cos_array[2, 1, ,1],cos_array[2, 1, ,2],cos_array[2, 1, ,3]),
    cov=rep(c("B1","B100","B5","B1","B100","B5"),c(ll,ll,ll,ll,ll,ll)),
    beta=rep(c("s","s","s","r","r","r"),c(ll,ll,ll,ll,ll,ll))
)

p1 <-ggplot(df1, aes(x=omega_seq, y=cos, color=beta,linetype=cov)) +
  geom_line(size=1.5)+ 
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 0.5)) +
  theme(legend.position ="right", legend.text = element_text(size=20),legend.title = element_text(size=20), legend.key.size = unit(1, 'cm'))+
  theme(
    plot.title = element_text(size=20),
    axis.text.x = element_text(size=15), 
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20))+
  xlab(expression(b)) + ylab("Empirical power") +
  scale_color_manual(values = cbbPalette[1:2],
                     labels = c(expression(delta^1),expression( delta^2)), name=expression(beta))+
  scale_linetype_manual(values=c(1,6,2),
                        labels = c(expression(Sigma[1]), 
                                  expression(Sigma[2]),expression(Sigma[3])), name="Covariance")
p1
   
 
df2<-data.frame(
  omega_seq=c(omega_seq,omega_seq,omega_seq,omega_seq,omega_seq,omega_seq),
  sq=c(square_array[1, 1, ,1],square_array[1, 1, ,2],square_array[1, 1, ,3], square_array[2, 1, ,1],square_array[2, 1, ,2],square_array[2, 1, ,3]),
  cov=rep(c("B1","B100","B5","B1","B100","B5"),c(ll,ll,ll,ll,ll,ll)),
  beta=rep(c("s","s","s","r","r","r"),c(ll,ll,ll,ll,ll,ll))
)

p2 <-ggplot(df2, aes(x=omega_seq, y=sq, color=beta,linetype=cov)) +
  geom_line(size=1.5)+ 
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 0.5)) +
  theme(legend.position ="right", legend.text = element_text(size=20),legend.title = element_text(size=20), legend.key.size = unit(1, 'cm'))+
  theme(
    plot.title = element_text(size=20),
    axis.text.x = element_text(size=15), 
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20))+
  xlab(expression(b)) + ylab("Empirical power") +
  
  scale_color_manual(values = cbbPalette[1:2],
                     labels = c(expression(delta^1),expression( delta^2)), name=expression(beta))+
  scale_linetype_manual(values=c(1,6,2),
                        labels = c(expression(Sigma[1]), 
                                   expression(Sigma[2]),expression(Sigma[3])), name="Covariance")
p2
 