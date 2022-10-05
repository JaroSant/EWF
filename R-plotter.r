rm(list = ls()) 
library(ggplot2) #Needs ggplot2 - run install.packages(ggplot2) if not already installed !!

cols <- hcl.colors(30, "Spectral")

Ne = 10000 #Quantities needed to recale time back into years before present 
g = 5

ImpT <- t(read.table("D:/MCMC4WF-ThetaZero/ImpT.txt", quote = "\"", comment.char = "")) #Input location of file ImpT.txt here! 
ImpT <-(((ImpT * 2 * Ne * g) - 20000)) #Converting from diffusion time to years before present
ImpHT <-(t(read.table("D:/MCMC4WF-ThetaZero/HorseTrajectories.txt", quote = "\"", comment.char = ""))) #Input location of file ImpHT.txt here!
OGT <-t(read.table("D:/MCMC4WF-ThetaZero/OGT.txt", quote = "\"", comment.char = "")) #Input location of file OGT.txt here! 
OGT <-(((OGT * 2 * Ne * g) - 20000)) #Converting from diffusion time to years before present
OGHT <-(t(read.table("D:/MCMC4WF-ThetaZero/OGHT.txt", quote = "\"", comment.char = ""))) #Input location of file OGHT.txt here

setEPS()
postscript(file = "horseTrajectories.eps") #Selection of figure format
matplot(ImpT, ImpHT, type = "l", col = cols, lwd = 1, lty = 1, xlab = "Time in years before present", ylab = "Frequency", ylim = c(0, 0.85)) #Plots all the generated paths
points(OGT, OGHT, pch = 4, bg = "black", col = "black", lwd = 4, cex = 1.5) #Superimposes the original observations dev.off()
