######## FUNCTION AND EXEMPLE OF A NEUTRAL EVOLUTION ALONG A PHYLOGENETIC TREE



VCVevol_brownian<-function(nb_ind,Vm,prob_mut,loci,nb_repeat=1,cpus=1,length.burn,tips,times,intratrait,EpiIntra){

  require(ape)
  require(picante)
  require(phytools)
  # require(snowfall) # to comment for cluster
  # require(cluster) # to comment for cluster
  require(mvMORPH)


  ##########################################################################
  ######################" PART 1 : FUNCTIONS TO LOAD #######################
  ##########################################################################"

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #------------------------------ Core function to drift the VCV matrix ----------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  evol_vcv_tree<-function(rep,phy,nb_ind,Vm,prob_mut,loci,traits,intratrait,EpiIntra){

    ##burning phase
    burned = IBMepis::ibm(nbInd = nb_ind,length = phy$lengthburned,CovarMatrix = Vm,prob_mut = prob_mut, loci = loci,
                        pursue = FALSE, ParentGenetValue = matrix(0),intratrait = intratrait,matrixEpistasisIntra = EpiIntra,
                        intertraits = FALSE,matrixEpistasisInter = matrix(0))


    phy$node.m0[phy$edge[1,1],,] <- burned$Parents_geneticValues
    phy$node.meanPhenoEpi[phy$edge[1,1],,]<- burned$OffspringsValueMeanEpi[phy$lengthburned,]
    phy$node.meanPhenoNoEpi[phy$edge[1,1],,]<- burned$OffspringsValueMeanNoEpi[phy$lengthburned,]



    ## propagation along tree
    edges = list(NULL)
    edges_pheno_epi= list(NULL)
    edges_pheno_noepi= list(NULL)

    for (i in 1:phy$Nedge) {
      edges[[i]] = matrix(NA,phy$edge.length[i]+1,(1+(traits+1)))
      edges_pheno_epi[[i]] = matrix(NA,phy$edge.length[i]+1,(traits+1))
      edges_pheno_noepi[[i]] = matrix(NA,phy$edge.length[i]+1,(traits+1))

      if (phy$edge[i,1]==phy$Nterm+1) {

        anc.age = 0
        anc.m0= phy$node.m0[phy$edge[1,1],,]
        anc.meanPhenoEpi= phy$node.meanPhenoEpi[phy$edge[1,1],,]
        anc.meanPhenoNoEpi= phy$node.meanPhenoNoEpi[phy$edge[1,1],,]

        anc.vcv.sum= burned$Gmat_sum[phy$lengthburned,]
        anc.vcv.bd= burned$Gmat_bd[phy$lengthburned,]

        # phy$node.vcv.sum[phy$edge[1,1],] = anc.vcv.sum
        phy$node.vcv.bd[phy$edge[1,1],] = anc.vcv.bd

      } else {
        anc.edge = match(phy$edge[i,1],phy$edge[,2])
        anc.age = phy$ages[match(phy$edge[i,1],phy$edge[,2])]
        anc.vcv.bd = phy$node.vcv.bd[phy$edge[i,1],]
        anc.m0= phy$node.m0[phy$edge[i,1],,]
        anc.meanPhenoEpi= phy$node.meanPhenoEpi[phy$edge[i,1],,]
        anc.meanPhenoNoEpi= phy$node.meanPhenoNoEpi[phy$edge[i,1],,]
      }
      edges[[i]][,1]<- anc.age:phy$ages[i]
      edges[[i]][1,2:(1+(traits+1))] <- anc.vcv.bd
      edges_pheno_epi[[i]][,1]<- anc.age:phy$ages[i]
      edges_pheno_epi[[i]][1,2:(traits+1)] <- anc.meanPhenoEpi
      edges_pheno_noepi[[i]][,1]<- anc.age:phy$ages[i]
      edges_pheno_noepi[[i]][1,2:(traits+1)] <- anc.meanPhenoNoEpi

      #--call for "ind_based_pleio" function   ---
      model_evo=IBMepis::ibm(length=phy$edge.length[i],nbInd = nb_ind,CovarMatrix = Vm,prob_mut = prob_mut,loci = loci,
                             pursue=TRUE,ParentGenetValue=anc.m0,intratrait = intratrait,matrixEpistasisIntra = EpiIntra,
                             intertraits = FALSE,matrixEpistasisInter = matrix(0))


      edges[[i]][-1,2:((traits+1)+1)]<- model_evo$Gmat_bd
      edges_pheno_epi[[i]][-1,2:(traits+1)] <-  model_evo$OffspringsValueMeanEpi
      edges_pheno_noepi[[i]][-1,2:(traits+1)] <-  model_evo$OffspringsValueMeanNoEpi

      phy$node.vcv.bd[phy$edge[i,2],] <- edges[[i]][nrow(edges[[i]]),2:(1+(traits+1))]
      phy$node.m0[phy$edge[i,2],,] <- model_evo$Parents_geneticValues
      phy$node.meanPheno[phy$edge[i,2],,] <- model_evo$OffspringsValueMean[phy$edge.length[i],,drop=FALSE]
      phy$node.meanPhenoEpi[phy$edge[i,2],,] <- model_evo$OffspringsValueMeanEpi[phy$edge.length[i],,drop=FALSE]
      phy$node.meanPhenoNoEpi[phy$edge[i,2],,] <- model_evo$OffspringsValueMeanNoEpi[phy$edge.length[i],,drop=FALSE]


    }

    results=list(vcv_mat_bd=edges, PhenoAlongTreeEpi=edges_pheno_epi,PhenoAlongTreeNoEpi=edges_pheno_noepi,meanPhenoEpi=phy$node.meanPhenoEpi)
    return(results)
  }


  ##########################################################################
  ######################" PART 2 : parallelise the simulation ##############
  #########################################################################"
  traits=dim(Vm)[1]

  #------set random seed

  tseed=as.numeric(Sys.time())
  set.seed((tseed - floor(tseed)) * 1e8 )

  #------phy construction following Revell et al 2007

  b<-exp((log(tips)-log(2))/times)-1
  phy<-pbtree(b=b,n=tips,t=times,type="discrete")

  #----- reorder phy

  phy <- reorder(phy, 'cladewise')
  phy$Nterm = length(phy$tip.label)
  phy$Nedge = nrow(phy$edge)

  phy$node.vcv.sum = matrix(NA,nrow=(phy$Nterm+phy$Nnode),ncol=(traits*traits)) #only work for 2 traits
  phy$node.vcv.bd= matrix(NA,nrow=(phy$Nterm+phy$Nnode),ncol=(traits+1)) #only work for 2 traits
  phy$node.m0 =  array(NA, c((phy$Nterm+phy$Nnode),loci*2*nb_ind,traits)) #2 pour le nombre d'alleles
  phy$node.meanPhenoEpi =  array(NA, c((phy$Nterm+phy$Nnode),1,traits)) #2 pour le nombre d'alleles
  phy$node.meanPhenoNoEpi =  array(NA, c((phy$Nterm+phy$Nnode),1,traits)) #2 pour le nombre d'alleles

  # phy$node.m0[phy$edge[1,1],,] <-array(NA,c(traits,loci*2*nb_ind))


  phy$ages = node.age(phy)$ages
  phy$lengthburned = length.burn


  #----- Parallel evaluation of the evol_vcv_tree of the function

  # sfInit(parallel=TRUE, cpus=cpus)
  # sfExportAll()


  simu<- lapply(1:nb_repeat,function(rep,phy,nb_ind,Vm,prob_mut,loci,traits,intratrait,EpiIntra)  evol_vcv_tree(rep,phy,nb_ind,Vm,prob_mut,loci,traits,intratrait,EpiIntra),
                  phy=phy,
                  nb_ind=nb_ind,
                  Vm=Vm,
                  prob_mut=prob_mut,
                  loci=loci,
                  traits=traits,
                  intratrait=intratrait,
                  EpiIntra=EpiIntra)
  # sfStop()


  ##### version for cluster without parallelisation
  #   nb_repeat=1
  #   simu<-lapply(1:nb_repeat,function(rep,phy,nb_ind,Vm,prob_mut,loci,traits)   evol_vcv_tree(rep,phy,nb_ind,Vm,prob_mut,loci,traits),
  #                phy=phy,
  #                nb_ind=nb_ind,
  #                Vm=Vm,
  #                prob_mut=prob_mut,
  #                loci=loci,
  #                traits=traits)
  ################

  Simu_arraylist=simplify2array(simu, higher = TRUE)
  Simu_array=apply(Simu_arraylist,1,function(x) simplify2array(x, higher = TRUE))


  # ~~~ estimate a Gmatrix as mean of Gmatrices estimated along tree

  sim_vcv_mat=apply(Simu_array$vcv_mat,1,function(y) array(unlist(y), dim = c(dim(y[[1]]), length(y))))

   rearrange<-function(x){
    x=do.call(rbind,x)
    x=x[-duplicated(x),,drop=FALSE] #remove dupication introduced in the propagation process for two sister branches
    x=apply(x,2,mean)
  }
  Gmean_simul_bd_raw=as.matrix(do.call(rbind, lapply(apply(Simu_array$vcv_mat,2,function(x) simplify2array(x,higher=TRUE)),rearrange))[,-1])
  Gmean_simul_bd=matrix(NA,nb_repeat,traits*traits)


  if(nb_repeat==1) {

    Gmean_simul_bd[,1]=Gmean_simul_bd_raw[1,]##ONLY FOR TWO TRAITS
    Gmean_simul_bd[,4]=Gmean_simul_bd_raw[2,]##ONLY FOR TWO TRAITS
    Gmean_simul_bd[,2:3]=Gmean_simul_bd_raw[3,]##ONLY FOR TWO TRAITS

  }else{

    Gmean_simul_bd[,1]=Gmean_simul_bd_raw[,1]##ONLY FOR TWO TRAITS
    Gmean_simul_bd[,4]=Gmean_simul_bd_raw[,2]##ONLY FOR TWO TRAITS
    Gmean_simul_bd[,2:3]=Gmean_simul_bd_raw[,3]##ONLY FOR TWO TRAITS

  }



  #  ~~~ estimate a theoretic Gmatrix from Lynch et al. 1986
  Gtheo=matrix(rep(c(2*nb_ind*(2*loci*prob_mut*Vm)),nb_repeat),ncol=traits*traits,byrow=T)

  #  ~~~ estimate phenotypic values at the tips

  # terms=which(phy$ages==times) #get the edge number corresponding to the tips ; or terms=which(phy$edge[,2] <= Ntip(phy))
  Pheno_mean=Simu_array$meanPhenoEpi[1:tips,,,,drop=FALSE] # to keep only info from tips
  Pheno_mean=array(Pheno_mean,c(tips,traits,nb_repeat))
  row.names(Pheno_mean)=phy$tip.label

  ############## Uncomment to have an unbiased estimate of sigma (biais is due to ML estimate)
  #    multBM<-function(x,y) mvBM(data=x,tree=y, model="BM1",method="pic",diagnostic = FALSE,echo = FALSE)$sigma
  #     correct=length(phy$tip.label)/(length(phy$tip.label)-1)
  #     Rate=t(apply(Pheno_mean,3,multBM ,y=phy))
  #     Rate_correct=Rate*correct
  #     R=t(apply(Pheno_mean,3,multBM,phy)) #to get sigma estimate
  #     print(rbind(R,(Gmean_simul/nb_ind),(Gtheo/nb_ind))) # to check difference in sigma rate)
  ###########################"

  likmultBM<-function(x,y) mvBM(data=x,tree=phy ,method="pic", scale.height = TRUE, model="BM1",diagnostic = FALSE,echo = FALSE)$LogLik


  multLL<-function(rep,phy,sigma,Pheno_mean,traits){
    sigma_matrix=matrix(sigma[rep,],traits,traits)
    Pheno_mean_byrep=Pheno_mean[,,rep]
    logl=mvLL(data=Pheno_mean_byrep,tree=phy,method="pic",param=list(estim=FALSE, sigma= sigma_matrix))$logl
    return(logl)
  }

  #estimate of different LogLikelihood

  LogLik_R=apply(Pheno_mean,3,likmultBM,phy) #estimation of Rate and LogLik (only LogLik is keet)
  LogLik_R0_Gmean_simul=sapply(1:nb_repeat,function(rep,phy,sigma,Pheno_mean,traits) multLL(rep,phy,sigma,Pheno_mean,traits),phy=phy,sigma=(Gmean_simul_bd/nb_ind),Pheno_mean=Pheno_mean,traits=traits)
  LogLik_R0_Gtheo=sapply(1:nb_repeat,function(rep,phy,sigma,Pheno_mean,traits) multLL(rep,phy,sigma,Pheno_mean,traits),phy=phy,sigma=(Gtheo/nb_ind),Pheno_mean=Pheno_mean,traits=traits)

  ddl=((traits^2-traits)/2)+traits
  Signi_chi=cbind(LogLik_R,LogLik_R0_Gmean_simul,LogLik_R0_Gtheo,ddl,pchisq(2*(LogLik_R-LogLik_R0_Gmean_simul),ddl ,lower.tail =FALSE),pchisq(2*(LogLik_R-LogLik_R0_Gtheo),ddl ,lower.tail =FALSE))
  colnames(Signi_chi)=c("LogLik_estim_rate","LogLik_Gmean_simul_rate","LogLik_Gtheo_rate","ddl","CHI2_simul","CHI2_theo")

  #simu_results=list(Signi_chi=Signi_chi,pheno.termbranch=sim_node.m0,Gtheo=Gtheo,Gmean_simul=Gmean_simul,VCVmat=sim_vcv_mat,phy=phy)
  simu_results=list(Signi_chi=Signi_chi,phy=phy,Gtheo=Gtheo,Gmean_simul_bd=Gmean_simul_bd,Pheno_mean_tips_Epi=Pheno_mean,Pheno_along_tree_Epi=Simu_array$PhenoAlongTreeEpi,Pheno_along_tree_NoEpi=Simu_array$PhenoAlongTreeNoEpi)
  return(simu_results)
}


################### Example


# Vm=matrix(c(0.05,0,0,0.20),2,2)
#
#
# nb_ind=100
# prob_mut=0.0025
# loci=20
# nb_repeat=1
# length.burn=1000
# cpus=1
# tips=50
# times=100
# intratrait=TRUE
#
# # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #
# #
# sd=1
# EpiIntra=matrix(rnorm(loci*loci,mean = -1,sd = sd),loci,loci)
# diag(EpiIntra)=0


# Epiinter=matrix(rnorm(loci*loci,mean = 0,sd = sd),loci,loci)
# diag(Epiinter)= 0
#
#
#   #note all replications use the same tree when nb_repeat>1
#
# results_simu=VCVevol_brownian(length.burn=length.burn,nb_ind=nb_ind,Vm=Vm,prob_mut=prob_mut,loci=loci,nb_repeat=nb_repeat,cpus=cpus,tips=tips,times=times,intratrait = TRUE,EpiIntra=EpiIntra)#,length.burn=length.burn
#
# results_simu$Signi_chi
#
#
# mvEB(results_simu$phy,results_simu$Pheno_mean_tips[,,3])
