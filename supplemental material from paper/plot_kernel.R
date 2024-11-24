library(ggplot2)
library(extrafont)
font_import()
loadfonts(device="win")       
fonts()
 
B1=as.vector(unlist(read.table(" ")))   ## load data w.r.t. Sigma1
B100=as.vector(unlist(read.table(" "))) ## load data w.r.t. Sigma2
B5=as.vector(unlist(read.table(" ")))   ## load data w.r.t. Sigma3

n=1000
df1<-data.frame(
  B1=B1, B100=B100, B5=B5
)

p1 <- ggplot(df1) + 
  stat_density(aes(x=B1, colour="B1", linetype="B1"),   size = 1.5,   geom="line", position="identity")+
  stat_density(aes(x=B100, colour="B100", linetype="B100"), size = 1.5,  geom="line", position="identity")+
  stat_density(aes(x=B5, colour="B5", linetype="B5"),   size = 1.5,  geom="line", position="identity")+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color="black", size=2.0)+
  xlab("") + ylab("Density") + xlim(-3,3)+ylim(0,0.6)+
  theme(text=element_text(family="Times New Roman", face="bold", size=25), legend.position = c(0.86,0.78),legend.text = element_text(size=25),legend.title = element_text(size=25), legend.key.size = unit(1.25, 'cm')) +  
  scale_color_manual(values=c("B1"="#FF9900", "B100"="#00CC66","B5"= "#0066FF"),
                     labels = c(expression(Sigma[1]), 
                                expression(Sigma[2]),expression(Sigma[3])), name="Covariance")+
  scale_linetype_manual(values=c("B1"=1,"B100"=2,"B5"=6),
                        breaks=c("B1", "B100", "B5"),
                        labels = c(expression(Sigma[1]), 
                                   expression(Sigma[2]),expression(Sigma[3])), name="Covariance")

p1
