
plot_ibm<-function(result_simu,traits,rep=1){

  PhenoAlongTreeEpi=result_simu$Pheno_along_tree_Epi
  PhenoAlongTreeNoEpi=result_simu$Pheno_along_tree_NoEpi

min_plot_t1_Epi=  min(do.call(rbind,PhenoAlongTreeEpi[,rep])[,2])
min_plot_t2_Epi=  min(do.call(rbind,PhenoAlongTreeEpi[,rep])[,3])
max_plot_t1_Epi=  max(do.call(rbind,PhenoAlongTreeEpi[,rep])[,2])
max_plot_t2_Epi=  max(do.call(rbind,PhenoAlongTreeEpi[,rep])[,3])


x_plot=  max(do.call(rbind,PhenoAlongTreeEpi[,rep])[,1])

par(mfrow=c(2,1))

plot(x=c(min_plot_t1,0),type='n',xlim=c(0,x_plot),ylim=c(min_plot_t1_Epi,max_plot_t1_Epi),xlab = "time",ylab="trait 1")
for (i in 1:dim(summary(PhenoAlongTreeEpi[,rep]))[1]) {
  lines(x=PhenoAlongTreeEpi[,rep][[i]][,1],y=PhenoAlongTreeEpi[,rep][[i]][,2])
  lines(x=PhenoAlongTreeEpi[,rep][[i]][,1],y=PhenoAlongTreeNoEpi[,rep][[i]][,2],col="red")
}

plot(x=c(min_plot_t2,0),type='n',xlim=c(0,x_plot),ylim=c(min_plot_t2_Epi,max_plot_t2_Epi),xlab = "time",ylab="trait 2")
for (i in 1:dim(summary(PhenoAlongTreeEpi[,rep]))[1]) {
   lines(x=PhenoAlongTreeEpi[,rep][[i]][,1],y=PhenoAlongTreeEpi[,rep][[i]][,3])
  lines(x=PhenoAlongTreeEpi[,rep][[i]][,1],y=PhenoAlongTreeNoEpi[,rep][[i]][,3],col="red")

}
}

