
## load the gerenerated data

name="ROCn4p1" ## give name to the files

sign_seq=seq(from=0.01,to=0.99,by=0.04)

pdf(file=paste0(name,"null","_fig.pdf"))

plot(sign_seq, sign_seq, type= "l", lty= 1, col= "black", lwd = 2,
   xlab="", ylab="", ann=F)
axis(1)
axis(2)
lines(sign_seq, as.vector(true.pro.array[1,1,]), type = "l", lty=1, xlim = c(0,1), ylim = c(0,1), lwd = 2, col="blue")
lines(sign_seq, as.vector(true.pro.array[1,2,]), type = "l", lty=6, lwd = 2, col="darkred")
lines(sign_seq, as.vector(true.pro.array[1,3,]), type = "l", lty=6, lwd = 2, col="darkgreen")
mtext(side = 1, line = 2.2, "Nominal level", cex=1.5)
mtext(side = 2, line = 2.0, "False positive rate", las = 0, cex=1.5)
mtext(side = 3, "")
legend("bottomright",c("RP","GC","MCL"), lty = c(1,6,6), lwd=2, col= c("blue","darkred","darkgreen"), bty = "o", cex=1.2)
grid()
dev.off()
dev.off()


pdf(file=paste0(name,"Sparse","_fig.pdf"))

plot( sign_seq,as.vector(true.pro.array[2,1,]), type= "l", lty=1, lwd=2, col="blue",
      xlab="",ylab="")
axis(1)
axis(2)
lines(sign_seq,as.vector(true.pro.array[2,2,]),lty=6, lwd=2, col="darkred")
lines(sign_seq,as.vector(true.pro.array[2,3,]),lty=6, lwd=2, col="darkgreen")
mtext(side = 1, line = 2.2, "False positive rate", cex = 1.5)
mtext(side = 2, line = 2.0, "True positive rate", las = 0, cex = 1.5)
mtext(side = 3, "")
legend("bottomright",c("RP","GC","MCL"), lty = c(1,6,6), lwd=2, col= c("blue","darkred","darkgreen"), bty = "o", cex=1.2)
grid() 
dev.off()
dev.off()
 

pdf(file=paste0(name, "Dense","_fig.pdf"))

plot( sign_seq,as.vector(true.pro.array[3,1,]), type= "l", lty=1, lwd=2,col="blue",
      xlab="",ylab="")
axis(1)
axis(2)
lines(sign_seq,as.vector(true.pro.array[3,2,]),lty=6,lwd=2, col="darkred")
lines(sign_seq,as.vector(true.pro.array[3,3,]),lty=6, lwd=2,col="darkgreen")
mtext(side = 1, line = 2.2, "False positive rate", cex = 1.5)
mtext(side = 2, line = 2.0, "True positive rate", las = 0, cex = 1.5)
mtext(side = 3, "")
legend("bottomright",c("RP","GC","MCL"), lty = c(1,6,6), lwd=2, col= c("blue","darkred","darkgreen"), bty = "o", cex=1.2)
grid() 

dev.off()
dev.off()

