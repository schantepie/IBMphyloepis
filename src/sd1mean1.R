a=commandArgs()[4]

#setwd("/media/mnhn/Leca/Gmatrix_project/1_Data_Gmatrix/Barn_swallow_spain_moller")
data <- read.csv("Gmat_final.csv", sep=",")
PedR <- read.csv("PedR.csv")

library(MCMCglmm)


nb_ind=100
prob_mut=0.0025
loci=20
nb_repeat=1
length.burn=10000
cpus=1
tips=100
times=100000
intratrait=TRUE
#
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#
sd=1
EpiIntra=matrix(rnorm(loci*loci,mean = 1,sd = sd),loci,loci)
diag(EpiIntra)=0


# Epiinter=matrix(rnorm(loci*loci,mean = 0,sd = sd),loci,loci)
# diag(Epiinter)= 0
#
#  #note all replications use the same tree when nb_repeat>1
#
 results_simu=IBMepis::VCVevol_brownian(length.burn=length.burn,nb_ind=nb_ind,Vm=Vm,prob_mut=prob_mut,loci=loci,nb_repeat=nb_repeat,cpus=cpus,tips=tips,times=times,intratrait = TRUE,EpiIntra=EpiIntra)#,length.burn=length.burn
#
            
                  

save.image(paste(a,"epis.Rdata",sep=""))











