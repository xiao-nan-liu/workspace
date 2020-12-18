options(useFancyQuotes = FALSE) 
library("Hmisc")
library("ADNIMERGE")
library("mice")
library(mvtnorm)
library(Matrix)
library(clusterGeneration)
library(ggplot2)
library(ROCR)
library(glmnet)
library("pROC")
library(lmtest)
help(package="ADNIMERGE")
#########################################################################
####################screening############################################
scr=function(X, y, M){
  p_value=rep(0,ncol(X))
  trld=0.1
  for(i in 1:nrow(X)){
    for(j in 1:ncol(X)){
      if(M[i,j]==0) X[i,j]=NA
    }
  }
  for(i in 1:ncol(X)){
    glm0=glm(y~X[,i], family="binomial", maxit=50)
    p_value[i]=lrtest(glm0)[2,5]
  }
  X1=X[,which(p_value<=trld)]
  p1=sum(which(p_value<=trld)<=12)
  p2=sum(which(p_value<=trld)>12&which(p_value<=trld)<=16)
  p3=sum(which(p_value<=trld)>=17)
  select=rep(0,17)
  select[p_value<=0.10]=1
  return(list(data.frame(X1, y),c(p1,p2,p3),which(p_value<=trld), select))
}
##############################################################################
###########################fold function######################################
fold=function(X1, X2, X3, X4, y1, y2, y3, y4, nfold, k){
  n1=nrow(X1);n2=nrow(X2);n3=nrow(X3);n4=nrow(X4);
  n=n1+n2+n3+n4;
  n1_vld=floor(n1/nfold);n2_vld=floor(n2/nfold);n3_vld=floor(n3/nfold);n4_vld=floor(n4/nfold);
  if (k<=nfold-1){
    val1=c((k-1)*n1_vld+1:n1_vld);val2=(k-1)*n2_vld+1:n2_vld;val3=(k-1)*n3_vld+1:n3_vld;val4=(k-1)*n4_vld+1:n4_vld;
    X_train=rbind(X1[-val1,], X2[-val2,],X3[-val3,], X4[-val4,])
    y_train=c(y1[-val1],y2[-val2],y3[-val3],y4[-val4])
    X_validate=rbind(X1[val1,], X2[val2,],X3[val3,], X4[val4,])
    y_validate=c(y1[val1],y2[val2],y3[val3],y4[val4])
  }
  if (k==nfold){
    val1=((k-1)*n1_vld+1):n1;val2=((k-1)*n2_vld+1):n2;val3=((k-1)*n3_vld+1):n3;val4=((k-1)*n4_vld+1):n4;
    X_train=rbind(X1[-val1,], X2[-val2,],X3[-val3,], X4[-val4,])
    y_train=c(y1[-val1],y2[-val2],y3[-val3],y4[-val4])
    X_validate=rbind(X1[val1,], X2[val2,],X3[val3,], X4[val4,])
    y_validate=c(y1[val1],y2[val2],y3[val3],y4[val4])
  }
  return(list(X_train=X_train, y_train=y_train, X_vld=X_validate, y_vld=y_validate, n_train=c(n1-length(val1), n2-length(val2),n3-length(val3),n4-length(val4)),
              n_vld=c(length(val1), length(val2),length(val3), length(val4))))
}
##############################################################################
#######################function estimation####################################
Estimate_fun=function(X0_train, y_train, X0_vld, y_vld, Miss_matrix, n1,n2,n3,n4,n1_vld,n2_vld,n3_vld,n4_vld){
  sect=rep(0,17)
  C1=scr(X0_train, y_train, Miss_matrix)
  sect=C1[[4]]
  X0_vld=X0_vld[,C1[[3]]];
  p_num=sum(C1[[2]])
  X0_train=as.matrix(C1[[1]][,1:p_num]); y_train=C1[[1]][,p_num+1]
  #X0_predictor=as.matrix(C1[[1]][,1:p_num]);
  p1=C1[[2]][1];p2=C1[[2]][2];p3=C1[[2]][3]
  C=glm(y_train~.,data=data.frame(y_train,X0_train[,1:p1]),family=binomial)$coef
  A_initial=matrix(0,nrow=p1,ncol=p2+p3)
  A_update=matrix(1,nrow=p1,ncol=p2+p3)
  b1_initial=matrix(0,nrow=1,ncol=p2+p3)
  b1_update=matrix(1,nrow=1,ncol=p2+p3)
  b2_initial=matrix(0,nrow=1,ncol=p2+p3)
  b2_update=matrix(1,nrow=1,ncol=p2+p3)
  Sigma_initial=as.matrix(diag(1,nrow=p2+p3, ncol=p2+p3))
  Sigma_update=as.matrix(diag(0.5,nrow=p2+p3, ncol=p2+p3))
  X0_est=as.matrix(X0_train)
  while(sum(abs(A_update-A_initial))>1e-7||sum(abs(b1_update-b1_initial))>1e-7||sum(abs(Sigma_update-Sigma_initial))>1e-7){
    A_initial=A_update;
    b1_initial=b1_update;
    b2_initial=b2_update;
    Sigma_initial=Sigma_update;
    ####E step####
    for(i in 1:n1){
      X0_est[i,(p1+1):p_num]=y_train[i]*(X0_est[i,1:p1]%*%A_update+b1_update)+
        (1-y_train[i])*(X0_est[i,1:p1]%*%A_update+b2_update)
    }
    for(i in n1+1:n2){
      if(p2>=1){
        X0_est[i,(p1+1):(p1+p2)]=y_train[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])+
        (1-y_train[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2])+
        (X0_est[i,(p1+p2+1):p_num]-
           y_train[i]*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b1_update[(p2+1):(p2+p3)])-
           (1-y_train[i])*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b2_update[(p2+1):(p2+p3)]))%*%
        solve(Sigma_update[(p2+1):(p2+p3),(p2+1):(p2+p3)])%*%Sigma_update[(p2+1):(p2+p3),1:p2]}
                                                                                              
    }
    for(i in (n1+n2)+1:n3){
      if(p2>=1){
        X0_est[i,(p1+p2+1):p_num]=y_train[i]*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b1_update[p2+1:p3])+
          (1-y_train[i])*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b2_update[p2+1:p3])+
          (X0_est[i,p1+1:p2]-y_train[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])-
             (1-y_train[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2]))%*%
          solve(Sigma_update[1:p2,1:p2])%*%Sigma_update[1:p2,(p2+1):(p2+p3)]
      }
      if(p2==0){
        X0_est[i,(p1+1):p_num]=y_train[i]*(X0_est[i,1:p1]%*%A_update+b1_update)+
          (1-y_train[i])*(X0_est[i,1:p1]%*%A_update+b2_update)
      }
    }
    ####M step####
    A_update=solve(t(X0_est[,1:p1])%*%X0_est[,1:p1])%*%t(X0_est[,1:p1])%*%(X0_est[,(p1+1):p_num]-as.matrix(y_train)%*%((b1_update))-as.matrix(1-y_train)%*%((b2_update)))
    X1=X0_est[y_train==1,];X2=X0_est[y_train==0,]
    b1_update=apply(X1[,(p1+1):p_num]-X1[,1:p1]%*%A_update,2,mean)
    b2_update=apply(X2[,(p1+1):p_num]-X2[,1:p1]%*%A_update,2,mean)
    Temp=matrix(0,p2+p3,p2+p3)
    Penal=matrix(0,p2+p3,p2+p3)
    for(i in 1:n){
      Temp=Temp+y_train[i]*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)+
        (1-y_train[i])*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)
    }
    Penal[1:p2,1:p2]=n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:p3]%*%solve(Sigma_initial[p2+1:p3,p2+1:p3])%*%Sigma_initial[p2+1:p3,1:p2])
    Penal[p2+1:p3,p2+1:p3]=n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,1:p2]%*%solve(Sigma_initial[1:p2,1:p2])%*%Sigma_initial[1:p2,p2+1:p3])
    Sigma_update=(Temp+n1*Sigma_initial+Penal)/n
  }
  ###############estimation####################
  if(p2>=1){
    beta1_est=C[2:(p1+1)]
    int1_est=C[1]
    beta2_1est=C[2:(p1+1)]-A_update[,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%
      (b1_update[p2+1:p3]-b2_update[p2+1:p3])
    beta2_3est=solve(Sigma_update[p2+1:p3,p2+1:p3])%*%(b1_update[p2+1:p3]-b2_update[p2+1:p3])
    int2_est=C[1]+(b2_update[p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b2_update[p2+1:p3]-
                     b1_update[p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b1_update[p2+1:p3])/2
    beta3_1est=C[2:(p1+1)]-A_update[,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%
      (b1_update[1:p2]-b2_update[1:p2])
    beta3_2est=solve(Sigma_update[1:p2,1:p2])%*%(b1_update[1:p2]-b2_update[1:p2])
    int3_est=C[1]+(b2_update[1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%b2_update[1:p2]-
                     b1_update[1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%b1_update[1:p2])/2
    beta4_1est=C[2:(p1+1)]-A_update%*%solve(Sigma_update)%*%(b1_update-b2_update)
    beta4_23est=solve(Sigma_update)%*%(b1_update-b2_update)
    int4_est=C[1]+(b2_update%*%solve(Sigma_update)%*%b2_update-
                     b1_update%*%solve(Sigma_update)%*%b1_update)/2
    eta1=X0_vld[1:n1_vld,1:p1]%*%beta1_est+int1_est
    eta2=X0_vld[n1_vld+1:n2_vld,1:p1]%*%beta2_1est+X0_vld[n1_vld+1:n2_vld,(p1+p2+1):p_num]%*%beta2_3est+rep(int2_est,n2_vld)
    eta3=X0_vld[(n1_vld+n2_vld)+1:n3_vld,1:p1]%*%beta3_1est+X0_vld[(n1_vld+n2_vld)+1:n3_vld,p1+1:p2]%*%beta3_2est+rep(int3_est,n3_vld)
    eta4=X0_vld[(n1_vld+n2_vld+n3_vld)+1:n4_vld,1:p1]%*%beta4_1est+X0_vld[(n1_vld+n2_vld+n3_vld)+1:n4_vld,p1+1:(p2+p3)]%*%beta4_23est+rep(int4_est,n4_vld)
    y1_predict=1/(1+exp(-eta1))
    y2_predict=1/(1+exp(-eta2))
    y3_predict=1/(1+exp(-eta3))
    y4_predict=1/(1+exp(-eta4))
    y_predict=c(y1_predict,y2_predict,y3_predict,y4_predict)
    
    glm_1=glm(y_train[1:n1]~X0_train[1:n1,1:p1],family="binomial",control = list(maxit = 50))
    glm_2=glm(y_train[n1+1:n2]~X0_train[n1+1:n2,c(1:p1,(p1+p2+1:p3))],family="binomial",control = list(maxit = 50))
    glm_3=glm(y_train[(n1+n2)+1:n3]~X0_train[(n1+n2)+1:n3,1:(p1+p2)],family="binomial",control = list(maxit = 50))
    glm_4=glm(y_train[(n1+n2+n3+1):n]~X0_train[(n1+n2+n3+1):n,1:p_num],family="binomial",control = list(maxit = 50))
    beta1=coef(glm_1)
    beta2=coef(glm_2)
    beta3=coef(glm_3)
    beta4=coef(glm_4)
    eta_1=cbind(1,X0_vld[1:n1_vld,1:p1])%*%beta1
    eta_2=cbind(1,X0_vld[n1_vld+1:n2_vld,c(1:p1,(p1+p2+1:p3))])%*%beta2
    eta_3=cbind(1,X0_vld[(n1_vld+n2_vld)+1:n3_vld,1:(p1+p2)])%*%beta3
    eta_4=cbind(1,X0_vld[(n1_vld+n2_vld+n3_vld+1):n_vld,1:p_num])%*%beta4
    y1_spredict=1/(1+exp(-eta_1))
    y2_spredict=1/(1+exp(-eta_2))
    y3_spredict=1/(1+exp(-eta_3))
    y4_spredict=1/(1+exp(-eta_4))
    y_spredict=c(y1_spredict,y2_spredict,y3_spredict,y4_spredict)
    
    beta_a1=C
    glm_a2=glm(y_train[c(n1+1:n2,(n-n4+1):n)]~X0_train[c(n1+1:n2,(n-n4+1):n),c(1:p1, (p1+p2)+1:p3)], family="binomial", control = list(maxit = 50))
    glm_a3=glm(y_train[n1+n2+1:(n3+n4)]~X0_train[n1+n2+1:(n3+n4),1:(p1+p2)], family="binomial", control=list(maxit=50))
    beta_a4=beta4
    beta_a2=coef(glm_a2)
    beta_a3=coef(glm_a3)
    eta_a1=cbind(1,X0_vld[1:n1_vld,1:p1])%*%beta_a1
    eta_a2=cbind(1,X0_vld[n1_vld+1:n2_vld,c(1:p1,(p1+p2+1:p3))])%*%beta_a2
    eta_a3=cbind(1,X0_vld[(n1_vld+n2_vld)+1:n3_vld,1:(p1+p2)])%*%beta_a3
    eta_a4=cbind(1,X0_vld[(n1_vld+n2_vld+n3_vld+1):n_vld,1:p_num])%*%beta_a4
    y1_apredict=1/(1+exp(-eta_a1))
    y2_apredict=1/(1+exp(-eta_a2))
    y3_apredict=1/(1+exp(-eta_a3))
    y4_apredict=1/(1+exp(-eta_a4))
    y_apredict=c(y1_apredict,y2_apredict,y3_apredict,y4_apredict)
    
    return(list(y_predict=y_predict, y_spredict=y_spredict, y_apredict=y_apredict, y_vld=y_vld, sect=sect))
  }
  if(p2==0){
    beta1_est=C[2:(p1+1)]
    int1_est=C[1]
    beta2_1est=C[2:(p1+1)]-A_update[,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%
      (b1_update[p2+1:p3]-b2_update[p2+1:p3])
    beta2_3est=solve(Sigma_update[p2+1:p3,p2+1:p3])%*%(b1_update[p2+1:p3]-b2_update[p2+1:p3])
    int2_est=C[1]+(b2_update[p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b2_update[p2+1:p3]-
                     b1_update[p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b1_update[p2+1:p3])/2
    beta3_1est=beta1_est
    int3_est=int1_est
    beta4_1est=beta2_1est
    beta4_3est=beta2_3est
    int4_est=int2_est
    
    eta1=X0_vld[1:n1_vld,1:p1]%*%beta1_est+int1_est
    eta2=X0_vld[n1_vld+1:n2_vld,1:p1]%*%beta2_1est+X0_vld[n1_vld+1:n2_vld,(p1+p2+1):p_num]%*%beta2_3est+rep(int2_est,n2_vld)
    eta3=X0_vld[(n1_vld+n2_vld)+1:n3_vld,1:p1]%*%beta3_1est+rep(int3_est,n3_vld)
    eta4=X0_vld[(n1_vld+n2_vld+n3_vld)+1:n4_vld,1:p1]%*%beta4_1est+X0_vld[(n1_vld+n2_vld+n3_vld)+1:n4_vld,p1+1:(p2+p3)]%*%beta4_23est+rep(int4_est,n4_vld)
    y1_predict=1/(1+exp(-eta1))
    y2_predict=1/(1+exp(-eta2))
    y3_predict=1/(1+exp(-eta3))
    y4_predict=1/(1+exp(-eta4))
    y_predict=c(y1_predict,y2_predict,y3_predict,y4_predict)
    
    glm_1=glm(y_train[1:n1]~X0_train[1:n1,1:p1],family="binomial",control = list(maxit = 50))
    glm_2=glm(y_train[n1+1:n2]~X0_train[n1+1:n2,c(1:p1,(p1+p2+1:p3))],family="binomial",control = list(maxit = 50))
    glm_3=glm(y_train[(n1+n2)+1:n3]~X0_train[(n1+n2)+1:n3,1:(p1+p2)],family="binomial",control = list(maxit = 50))
    glm_4=glm(y_train[(n1+n2+n3+1):n]~X0_train[(n1+n2+n3+1):n,1:p_num],family="binomial",control = list(maxit = 50))
    beta1=coef(glm_1)
    beta2=coef(glm_2)
    beta3=coef(glm_3)
    beta4=coef(glm_4)
    eta_1=cbind(1,X0_vld[1:n1_vld,1:p1])%*%beta1
    eta_2=cbind(1,X0_vld[n1_vld+1:n2_vld,c(1:p1,(p1+p2+1:p3))])%*%beta2
    eta_3=cbind(1,X0_vld[(n1_vld+n2_vld)+1:n3_vld,1:(p1+p2)])%*%beta3
    eta_4=cbind(1,X0_vld[(n1_vld+n2_vld+n3_vld+1):n_vld,1:p_num])%*%beta4
    y1_spredict=1/(1+exp(-eta_1))
    y2_spredict=1/(1+exp(-eta_2))
    y3_spredict=1/(1+exp(-eta_3))
    y4_spredict=1/(1+exp(-eta_4))
    y_spredict=c(y1_spredict,y2_spredict,y3_spredict,y4_spredict)
    
    beta_a1=C
    glm_a2=glm(y_train[c(n1+1:n2,(n-n4+1):n)]~X0_train[c(n1+1:n2,(n-n4+1):n),c(1:p1, (p1+p2)+1:p3)], family="binomial", control = list(maxit = 50))
    glm_a3=glm(y_train[n1+n2+1:(n3+n4)]~X0_train[n1+n2+1:(n3+n4),1:(p1+p2)], family="binomial", control=list(maxit=50))
    beta_a4=beta4
    beta_a2=coef(glm_a2)
    beta_a3=coef(glm_a3)
    eta_a1=cbind(1,X0_vld[1:n1_vld,1:p1])%*%beta_a1
    eta_a2=cbind(1,X0_vld[n1_vld+1:n2_vld,c(1:p1,(p1+p2+1:p3))])%*%beta_a2
    eta_a3=cbind(1,X0_vld[(n1_vld+n2_vld)+1:n3_vld,1:(p1+p2)])%*%beta_a3
    eta_a4=cbind(1,X0_vld[(n1_vld+n2_vld+n3_vld+1):n_vld,1:p_num])%*%beta_a4
    y1_apredict=1/(1+exp(-eta_a1))
    y2_apredict=1/(1+exp(-eta_a2))
    y3_apredict=1/(1+exp(-eta_a3))
    y4_apredict=1/(1+exp(-eta_a4))
    y_apredict=c(y1_apredict,y2_apredict,y3_apredict,y4_apredict)
    
    return(list(y_predict=y_predict, y_spredict=y_spredict, y_apredict=y_apredict, y_vld=y_vld,sect=sect))
  }
}
Accuracy=function(y1,ytrue){
  cutoff=seq(0.01,0.99,length=99)
  Acc=rep(0,99)
  Sen=rep(0,99)
  Spe=rep(0,99)
  cost=rep(0,99)
  y_temp=y1
  for(i in 1:99){
    y_temp[y1>=cutoff[i]]=1
    y_temp[y1<cutoff[i]]=0
    TP=length(which(y_temp==1&ytrue==1))
    FP=length(which(y_temp==0&ytrue==1))
    FN=length(which(y_temp==1&ytrue==0))
    TN=length(which(y_temp==0&ytrue==0))
    Acc[i]=(TP+TN)/length(y1)
    Sen[i]=TP/(TP+FN)
    Spe[i]=TN/(FP+TN)
  }
  cost=(1-Sen)^2+(1-Spe)^2
  i=which(cost==min(cost))[1]
  return(c(Acc[i],Sen[i],Spe[i],cutoff[i]))
}
##############################################################################
###############data preparation###############################################
#MCIdata=read.csv("C:/Users/Xiaonan/Dropbox/missing modality simulation/MCI.csv",header=T)
MCIdata=read.csv("C:/Users/xliu203/Dropbox/missing modality simulation/MCI.csv",header=TRUE)
X0=MCIdata[,c(6,11:14,22:24,25, 27, 29,30,32,33,43:48,64,67,69,70)]
X1=X0[,15:20]
lmba=eigen(cor(data.matrix(X1)))$values[1]
eigen1=data.matrix(X1)%*%eigen(cor(data.matrix(X1)))$vectors[,1]
ape=X0[,11]
APOEHT=rep(1,214)
APOEHT[ape=="NC"]=0
X00=data.matrix(X0[,c(3:5)])/X0[,2]
X0_predictor=data.frame(X00,X0[,c(6:8)],APOEHT, X0[,c(9,10,12:14,21:24)],eigen1)
for(i in 1:17){
  if(i!=4&i!=5&i!=7){
    X0_predictor[,i]=(X0_predictor[,i]-mean(X0_predictor[,i],na.rm=TRUE))/sd(X0_predictor[,i],na.rm=TRUE)
  }
}

tempData=mice(X0_predictor,m=5,maxit=50,meth='pmm',seed=500)
X0_predictor <- complete(tempData,1)
X0_predictor=as.matrix(X0_predictor)
y_response=rep(0,214)
y_response[which(X0$AV45.v1.Pos118=="+")]=1

AUC=matrix(0,nrow=50, ncol=3)
Acc=matrix(0,nrow=50, ncol=3)
Sen=matrix(0,nrow=50, ncol=3)
Spe=matrix(0,nrow=50, ncol=3)
colnames(AUC)=c("AUC_proposed", "AUC_sm", "AUC_avail")
colnames(Acc)=c("AUC_proposed", "AUC_sm", "AUC_avail")
colnames(Sen)=c("AUC_proposed", "AUC_sm", "AUC_avail")
colnames(Spe)=c("AUC_proposed", "AUC_sm", "AUC_avail")
select=rep(0,17)
for(l in 1:50){
a1=sample(1:214)
X0_predictor=X0_predictor[a1,]
y_response=y_response[a1]
n.grp=4; n.mod=3;
n.r=matrix(0,4,n.grp)
n.r[1,]=c(2.5,2.5,2.5,2.5)/10
n=214;p1=12;p2=4;p3=1;p_num=p1+p2+p3;q1=1
n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);n4=n-n1-n2-n3;

Miss_matrix=matrix(1,n,p1+p2+p3)
Miss_matrix[1:n1,p1+1:(p2+p3)]=0;
Miss_matrix[n1+1:n2,p1+1:p2]=0;
Miss_matrix[n1+n2+1:n3,p1+p2+1:p3]=0;
y1_response=y_response[1:n1];y2_response=y_response[n1+1:n2];y3_response=y_response[n1+n2+1:n3];y4_response=y_response[n1+n2+n3+1:n4];
X1_predictor=X0_predictor[1:n1,];X2_predictor=X0_predictor[n1+1:n2,];X3_predictor=X0_predictor[n1+n2+1:n3,];X4_predictor=X0_predictor[n1+n2+n3+1:n4,];
nfold=5

for(k in 1:nfold){
  fold1=fold(X1_predictor,X2_predictor,X3_predictor,X4_predictor,y1_response, y2_response, y3_response, y4_response, nfold, k)
  X0_train=fold1$X_train
  X0_vld=fold1$X_vld
  y_train=fold1$y_train
  y_vld=fold1$y_vld
  n1=fold1$n_train[1];n2=fold1$n_train[2];n3=fold1$n_train[3];n4=fold1$n_train[4];
  n1_vld=fold1$n_vld[1];n2_vld=fold1$n_vld[2];n3_vld=fold1$n_vld[3];n4_vld=fold1$n_vld[4];
  n=n1+n2+n3+n4;n_vld=n1_vld+n2_vld+n3_vld+n4_vld
  Miss_matrix=matrix(1,n,p1+p2+p3)
  Miss_matrix[1:n1,p1+1:(p2+p3)]=0;
  Miss_matrix[n1+1:n2,p1+1:p2]=0;
  Miss_matrix[n1+n2+1:n3,p1+p2+1:p3]=0;
  es1=Estimate_fun(X0_train, y_train, X0_vld, y_vld, Miss_matrix, n1,n2,n3,n4,n1_vld,n2_vld,n3_vld,n4_vld)
  select=select+es1$sect
  y_predict=NULL;
  y_spredict=NULL;
  y_apredict=NULL;
  y_test=NULL;
  y_predict=c(y_predict, es1$y_predict)
  y_spredict=c(y_spredict, es1$y_spredict)
  y_apredict=c(y_apredict, es1$y_apredict)
  y_test=c(y_test, es1$y_vld)
}

    
    pred1=prediction(y_predict, y_test)
    perf1=performance(pred1,"tpr","fpr")
    auc.tmp1=performance(pred1,"auc")
    auc1=as.numeric(auc.tmp1@y.values)
    ci.auc(y_vld,y_predict)*1
    
    
    pred2=prediction(y_spredict, y_test)
    perf2=performance(pred2,"tpr","fpr")
    auc.tmp2=performance(pred2,"auc")
    auc2=as.numeric(auc.tmp2@y.values)
    ci.auc(y_vld,y_spredict)*1
    
    pred3=prediction(y_apredict, y_test)
    perf3=performance(pred3,"tpr","fpr")
    auc.tmp3=performance(pred3,"auc")
    auc3=as.numeric(auc.tmp3@y.values)
    ci.auc(y_vld,y_apredict)*1
    AUC[l,1]=auc1
    AUC[l,2]=auc2
    AUC[l,3]=auc3
    Acc[l,1]=Accuracy(y_predict, y_test)[1]
    Acc[l,2]=Accuracy(y_spredict, y_test)[1]
    Acc[l,3]=Accuracy(y_apredict, y_test)[1]
    Sen[l,1]=Accuracy(y_predict, y_test)[2]
    Sen[l,2]=Accuracy(y_spredict, y_test)[2]
    Sen[l,3]=Accuracy(y_apredict, y_test)[2]
    Spe[l,1]=Accuracy(y_predict, y_test)[3]
    Spe[l,2]=Accuracy(y_spredict, y_test)[3]
    Spe[l,3]=Accuracy(y_apredict, y_test)[3]
}
apply(AUC,2,mean)
apply(AUC,2,sd)
apply(AUC,2,min)
apply(AUC,2,max)
t.test(AUC[,1],AUC[,2])
t.test(AUC[,1],AUC[,3]) 

apply(Sen,2,mean)
apply(Sen,2,sd)
apply(Sen,2,min)
apply(Sen,2,max)
t.test(Sen[,1],Sen[,2])
t.test(Sen[,1],Sen[,3])   

apply(Spe,2,mean)
apply(Spe,2,sd)
apply(Spe,2,min)
apply(Spe,2,max)
t.test(Spe[,1],Spe[,2])
t.test(Spe[,1],Spe[,3])   

apply(Acc,2,mean)
apply(Acc,2,sd)
apply(Acc,2,min)
apply(Acc,2,max)
t.test(Acc[,1],Acc[,2])
t.test(Acc[,1],Acc[,3])   
    



roc.test(y_vld[42:54],as.vector(y4_predict), as.vector(y4_spredict))
    


apply(AUCTL,2,mean)
apply(AUCSM,2,mean)
Model=c(rep("Transfer Learning",10),rep("Separate Model",10))
stdall=c(std,std)
Mn=c(apply(AUCTL,2,mean),apply(AUCSM,2,mean))
SD=c(apply(AUCTL,2,sd),apply(AUCSM,2,sd))
maxauc=Mn+SD
minauc=Mn-SD
Data=data.frame(inctr,stdall, Mn, maxauc, minauc )

ggplot(Data, aes(x=stdall, y=Mn,color=Model))+geom_point()+ geom_line()+
  scale_x_continuous(name="Standard deviation", breaks=seq(from=0,to=1,by=0.1))+
  scale_y_continuous(name="AUC", limits=c(0.05,1),breaks=seq(from=0,to=1,by=0.1))+
  geom_errorbar(aes(ymin=minauc, ymax=maxauc),width=0.05)+
  theme(plot.title = element_text(size = 20,face="plain"))+
  theme(axis.text.x=element_text(colour="black",size=14),axis.text.y=element_text(colour="black",size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_rect(fill = NA, colour = "black"))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16),legend.title=element_text(size=16))
pvalue=rep(0,10)
for(i in 1:10){
  pvalue[i]=t.test(AUCSM[,i], mu=0.792, alternative="less")$p.value
}
