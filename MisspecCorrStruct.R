#########################################################################
# Functions to accompany the Shiny app for 
# "Estimating variance components when random effect correlation structure is misspecified"
# by J Kasza and A Forbes.
#
# 2018-02-27: vartreat_1v3 modified so that model 1 is the crossed model with no interaction
#########################################################################
#Functions to generate design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}
CRXOdesmat<- function(T) {
  if((T-1)%%2 == 0) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T-1)
    Xcrxo[1:(T-1)/2, seq(1,T,2)] <- 1
    Xcrxo[((T-1)/2 + 1):(T-1), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
  
  if((T-1)%%2 == 1) {
    
    Xcrxo <- matrix(data=0, ncol = T, nrow = T)
    Xcrxo[1:(T)/2, seq(1,T,2)] <- 1
    Xcrxo[((T)/2 + 1):(T), seq(2,T,2)] <- 1
    return(Xcrxo)    
  }
  
}
plleldesmat <- function(T) {
  if((T-1)%%2 == 0) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T -1)
    Xpllel[1:(T-1)/2,] <- 1
    return(Xpllel)
  }
  
  
  if((T-1)%%2 == 1) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T)
    Xpllel[1:(T)/2,] <- 1
    return(Xpllel)
  }
  
}

pllelBLdesmat <- function(T) {
  #For now assume 50% of periods are baseline measurements.
  
  if((T-1)%%2 == 0) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T -1)
    Xpllel[1:(T-1)/2,((T+1)/2):T] <- 1
    return(Xpllel)
  }
  
  
  if((T-1)%%2 == 1) {
    Xpllel <- matrix(data=0, ncol = T, nrow =T)
    Xpllel[1:(T)/2,(T/2 + 1):T] <- 1
    return(Xpllel)
  }
  
}

#########################################################################
#Functions to calculate variance of the treatment effect estimator
TreatEffVar <- function(T, m, sig2E, sig2C, sig2CP, design) {
  
  sig2 <- sig2E/m
  
  if(design == 1)  Xmat <- SWdesmat(T)
  else if (design == 2)  Xmat <- plleldesmat(T)
  else if (design == 4)  Xmat <- CRXOdesmat(T)
  else if (design == 5)  Xmat <- pllelBLdesmat(T)
  
  K <- nrow(Xmat) 
  Xvec <-  as.vector(t(Xmat))
  
  Vi <-diag(sig2 + sig2CP, T) + matrix(data=sig2C, nrow=T, ncol=T) 
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(Xmat)%*%solve(Vi)%*%(matrix(colSums(Xmat),nrow=T, ncol=1))/K )
  return(vartheta)
}
TreatEffVarEXP <- function(T, m, r0, sig2E, sig2C, sig2CP, design) {
  
  sig2 <- sig2E/m
  
  if(design == 1)  Xmat <- SWdesmat(T)
  else if (design == 2)  Xmat <- plleldesmat(T)
  else if (design == 4)  Xmat <- CRXOdesmat(T)
  else if (design == 5)  Xmat <- pllelBLdesmat(T)
  
  K <- nrow(Xmat) 
  Xvec <-  as.vector(t(Xmat))
  
  Vi <- diag(sig2,T) + (sig2C + sig2CP)*(r0^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  
  
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(Xmat)%*%solve(Vi)%*%(matrix(colSums(Xmat),nrow=T, ncol=1))/K )
  return(vartheta)
}
#########################################################################
#Functions to calculate the ratio of the variances when the within-cluster 
#correlation structure is incorrectly specified as model 2 or 1.
vartreat_2v3 <- function(T, m, rho0, r, totalvar=1, design){
  
  sig2alpha <- rho0*totalvar
  sig2E <- totalvar - sig2alpha
  
  #Expected values of variance components estimated using HG
  rterms <-  (r^(abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE)- matrix(1:T,nrow=T, ncol=T, byrow=TRUE))))
  exp_sigmagamma <- sig2alpha*(T/(T-1) - sum(rterms)/(T*(T-1)))
  exp_sigalpha <- sig2alpha*(sum(rterms)/(T*(T-1)) - 1/(T-1))
  #Variance of the treatment effect using these terms:
  varModel2 <- TreatEffVar(T, m, sig2E, exp_sigalpha, exp_sigmagamma, design)
  
  #Variance of the treatment effect using the exp decay model:
  varModel3 <- TreatEffVarEXP(T, m, r, sig2E, sig2alpha, 0, design)
  
  return(varModel2/varModel3)
  
}

vartreat_1v3 <- function(T, m, rho0, r, totalvar=1, design){
  
  sig2alpha <- rho0*totalvar
  sig2E <- totalvar - sig2alpha
  
  #Expected values of variance components estimated using HH
  #OLD VERSIONS
  #{
  #rterms <-  (r^(abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE)- matrix(1:T,nrow=T, ncol=T, byrow=TRUE))))
  #exp_sigeps <- sig2E + sig2alpha*(T*m/(T*m-1) - m*sum(rterms)/(T*(T*m-1)))
  #exp_sigalpha <- sig2alpha*(m*sum(rterms)/(T*(T*m-1)) - 1/(T*m-1))
  #}
  
  #NEW VERSIONS (two-way crossed model with no random interaction)
  #{
  if(design == 1)  Xmat <- SWdesmat(T)
  else if (design == 2)  Xmat <- plleldesmat(T)
  else if (design == 4)  Xmat <- CRXOdesmat(T)
  K <- nrow(Xmat) 
  
  rterms <-  (r^(abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE)- matrix(1:T,nrow=T, ncol=T, byrow=TRUE))))
  exp_sigeps <- sig2E + sig2alpha*( (T*m*(K-1))/(K*T*m - K -T + 1) - (K-1)*m*sum(rterms)/(T*(K*T*m-K-T+1)))
  exp_sigalpha <- sig2alpha*((K*m-1)*sum(rterms)/(T*(K*T*m-K-T+1)) - (K-1)/(K*T*m-K-T+1))
  #}
  
  #Variance of the treatment effect using these terms:
  varModel1 <- TreatEffVar(T, m, exp_sigeps, exp_sigalpha, 0, design)
  
  #Variance of the treatment effect using the exp decay model:
  varModel3 <- TreatEffVarEXP(T, m, r, sig2E, sig2alpha, 0, design)
  
  return(varModel1/varModel3)
  
}


#Plot the relative variances for a range of Ts and rs
vartreat_versus3_plot <- function(design, m, rho0, model){
  
  #Sequence of periods considered
  Tseq <- seq(4, 20 ,1) 
  #Sequence of decay values considered
  rseq <- seq(0, 1, 0.05)
  
  Trvars <- matrix(data=NA, nrow= length(Tseq), ncol=length(rseq))
  for(t in 1:length(Tseq)) {
    for(rind in 1:length(rseq)) {
      if(model==2) Trvars[t,rind] <- vartreat_2v3(Tseq[t], m, rho0, rseq[rind], totalvar=1, design)
      else if(model==1) Trvars[t,rind] <- vartreat_1v3(Tseq[t], m, rho0, rseq[rind], totalvar=1, design)
    }
  }
  
  #want to highlight cells with rel var between 0.99 and 1.01
  Trvarsind <- matrix(data=NA, nrow= length(Tseq), ncol=length(rseq))
  for(t in 1:length(Tseq)) {
    for(rind in 1:length(rseq)) {
      if(Trvars[t,rind]> 0.99 && Trvars[t,rind]< 1.01) Trvarsind[t,rind] <- 1
    }
  }
  
  #Plot the results using a contour plot 
  Trvars<-round(Trvars, 4)
  melted_Trvars <- melt(Trvars)
  names(melted_Trvars)[names(melted_Trvars)=="Var1"] <- "T"
  names(melted_Trvars)[names(melted_Trvars)=="Var2"] <- "r"
  
  Tvec <- as.vector(matrix(data= Tseq, nrow=length(Tseq), ncol=length(rseq), byrow=FALSE))
  rvec <- as.vector(matrix(data= rseq, nrow=length(Tseq), ncol=length(rseq), byrow=TRUE))
  melted_Trvars$Tseq <- Tvec
  melted_Trvars$rseq <- rvec
  
  melted_Trvarsind <- melt(Trvarsind)
  names(melted_Trvarsind)[names(melted_Trvarsind)=="Var1"] <- "T"
  names(melted_Trvarsind)[names(melted_Trvarsind)=="Var2"] <- "r"
  melted_Trvarsind$Tseq <- Tvec
  melted_Trvarsind$rseq <- rvec 
  
  
  #A workaround if the max or min value is 1, to ensure colours displayed appropriately
  mycolours <- c("blue","white","red")
  myvals <- rescale(c(min(melted_Trvars$value), 1, max(melted_Trvars$value)))
  
  if(max(melted_Trvars$value) == 1) {
    mycolours <- c("blue","white")
    myvals <- rescale(c(min(melted_Trvars$value), 1)) 
  } 
  
  if(min(melted_Trvars$value) == 1) {
    mycolours <- c("white", "red")
    myvals <- rescale(c(1, max(melted_Trvars$value) ) )  } 
  
  #scale_fill_gradientn(colours=c("blue","white","red"),
  #values  = rescale(c(min(melted_Trvars$value), 1, max(melted_Trvars$value)))) +
  
  myplot <- ggplot(melted_Trvars, aes(x=rseq, y=Tseq)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_gradientn(colours= mycolours,
                         values  = myvals) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme(aspect.ratio=3/5, legend.position ="none", legend.key.size = unit(1, "cm"), 
          legend.text=element_text(size=12), 
          legend.background = element_rect(fill="grey95")) +
    coord_fixed() + xlab("Decay parameter, r") +  ylab("Number of periods, T") +
    guides(fill=guide_legend(nrow=1, keywidth=2, unit="cm")) +
    geom_text(aes(rseq, Tseq, label = round(value,2)), color = "black", size = 3) 
  
  myplot
}


#########################################################
# Functions that take a user-input design matrix and 
# plot the variance ratio for a range of decays for 
# that particular design matrix

#########################################################################
#Functions to calculate the ratio of the variances when the within-cluster 
#correlation structure is incorrectly specified as model 2 or 1.
TreatEffVar_DESMAT <- function(T, m, sig2E, sig2C, sig2CP, desmat) {
  
  sig2 <- sig2E/m
  
  Xmat <- desmat
  
  K <- nrow(Xmat) 
  Xvec <-  as.vector(t(Xmat))
  
  Vi <-diag(sig2 + sig2CP, T) + matrix(data=sig2C, nrow=T, ncol=T) 
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(Xmat)%*%solve(Vi)%*%(matrix(colSums(Xmat),nrow=T, ncol=1))/K )
  return(vartheta)
}
TreatEffVarEXP_DESMAT <- function(T, m, r0, sig2E, sig2C, sig2CP, desmat) {
  
  sig2 <- sig2E/m
  
  Xmat <- desmat
  
  K <- nrow(Xmat) 
  Xvec <-  as.vector(t(Xmat))
  
  Vi <- diag(sig2,T) + (sig2C + sig2CP)*(r0^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  
  
  vartheta <-  1/(t(Xvec)%*%(diag(1,K)%x%solve(Vi))%*%Xvec -colSums(Xmat)%*%solve(Vi)%*%(matrix(colSums(Xmat),nrow=T, ncol=1))/K )
  return(vartheta)
}


vartreat_2v3_DESMAT <- function(T, m, rho0, r, totalvar=1, desmat){
  
  sig2alpha <- rho0*totalvar
  sig2E <- totalvar - sig2alpha
  
  #Expected values of variance components estimated using HG
  rterms <-  (r^(abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE)- matrix(1:T,nrow=T, ncol=T, byrow=TRUE))))
  exp_sigmagamma <- sig2alpha*(T/(T-1) - sum(rterms)/(T*(T-1)))
  exp_sigalpha <- sig2alpha*(sum(rterms)/(T*(T-1)) - 1/(T-1))
  #Variance of the treatment effect using these terms:
  varModel2 <- TreatEffVar_DESMAT(T, m, sig2E, exp_sigalpha, exp_sigmagamma, desmat)
  
  #Variance of the treatment effect using the exp decay model:
  varModel3 <- TreatEffVarEXP_DESMAT(T, m, r, sig2E, sig2alpha, 0, desmat)
  
  return(varModel2/varModel3)
  
}

vartreat_1v3_DESMAT <- function(T, m, rho0, r, totalvar=1, desmat){
  
  sig2alpha <- rho0*totalvar
  sig2E <- totalvar - sig2alpha
  
  K <- nrow(desmat) 
  
  rterms <-  (r^(abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE)- matrix(1:T,nrow=T, ncol=T, byrow=TRUE))))
  exp_sigeps <- sig2E + sig2alpha*( (T*m*(K-1))/(K*T*m - K -T + 1) - (K-1)*m*sum(rterms)/(T*(K*T*m-K-T+1)))
  exp_sigalpha <- sig2alpha*((K*m-1)*sum(rterms)/(T*(K*T*m-K-T+1)) - (K-1)/(K*T*m-K-T+1))
  
  
  #Variance of the treatment effect using these terms:
  varModel1 <- TreatEffVar_DESMAT(T, m, exp_sigeps, exp_sigalpha, 0, desmat)
  
  #Variance of the treatment effect using the exp decay model:
  varModel3 <- TreatEffVarEXP_DESMAT(T, m, r, sig2E, sig2alpha, 0, desmat)
  
  return(varModel1/varModel3)
  
}


#Plot the relative variances for a range of Ts and rs
vartreat_versus3_plotDESMAT <- function(desmat, m, rho0, model){
  
  #Sequence of decay values considered
  rseq <- seq(0, 1, 0.01)
  
  Tperiods<- ncol(desmat)
  
  Trvars <- matrix(data=NA, nrow=length(rseq), ncol=1)
    for(rind in 1:length(rseq)) {
      if(model==2) Trvars[rind,1] <- vartreat_2v3_DESMAT(Tperiods, m, rho0, rseq[rind], totalvar=1, desmat)
      else if(model==1) Trvars[rind,1] <- vartreat_1v3_DESMAT(Tperiods, m, rho0, rseq[rind], totalvar=1, desmat)
    }
  
  
  #Plot the results
  Trvars <- as.data.frame(Trvars)
  #names(Trvars)[names(melted_Trvars)=="V1"] <- "Ratio"
  Trvars$rseq <- rseq
  
  myplot<- plot_ly(Trvars, x = ~rseq, y = ~V1, type = 'scatter', mode = 'lines',
                   line = list(color = "black", width = 4)) %>%
    layout(xaxis = list(title = "Decay parameter, r"), yaxis = list (title = "Variance Ratio"),
           legend=list(orientation="h", xanchor="center", y=1.1, x=0.5)) 
  
  myplot
}

