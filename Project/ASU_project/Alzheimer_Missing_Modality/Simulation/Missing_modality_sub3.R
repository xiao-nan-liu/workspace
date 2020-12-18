library(mvtnorm)
library(Matrix)
library(clusterGeneration)
library(ggplot2)
library(ROCR)
library(glmnet)
library("pROC")
###################################################################################
#############################logistic regression###################################
AUC=NULL
for(l in 1:100){
  n=300;p1=10;p2=5;p3=5;p_num=p1+p2+p3
  Sigma_X1=matrix(0.6,p1,p1)
  diag(Sigma_X1)=1
  Sigma_X2=matrix(0.6,p2,p2)
  diag(Sigma_X2)=1
  Sigma_X3=matrix(0.6,p3,p3)
  diag(Sigma_X3)=1
  Sigma_X=matrix(0,p_num+1,p_num+1)
  Sigma_X[1:p1,1:p1]=Sigma_X1
  Sigma_X[(p1+1):(p1+p2),(p1+1):(p1+p2)]=Sigma_X2
  Sigma_X[(p1+p2+1):(p_num),(p1+p2+1):(p_num)]=Sigma_X3
  Sigma_X[p_num+1,p_num+1]=1
  Sigma_X[p_num+1,1:p_num]=0
  Sigma_X[1:p_num,p_num+1]=0
  beta=c(rep(0.2,p1),rep(0.2,p2),rep(0.2,p3),1)
  mu=c(rep(0,p1), rep(0,p2),rep(0,p3),0)
  X0_predictor=NULL
  for(i in 1:n){
    X0_predictor=rbind(X0_predictor,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  lim=X0_predictor%*%beta
  pai=1/(1+exp(lim))
  
  y_response=rbinom(n,size=1,prob=pai)
  
  X0=cbind(X0_predictor,y_response)
  Truesigma=Sigma_X[(p1+1):(p_num),(p1+1):(p_num)]-Sigma_X[(p1+1):(p_num),1:p1]%*%solve(Sigma_X[1:p1,1:p1])%*%Sigma_X[1:p1,(p1+1):(p_num)]
  ############################################################################
  ####################validation set##########################################
  n_vld=100
  X0_vld=NULL
  for(i in 1:n_vld){
    X0_vld=rbind(X0_vld,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  lim_vld=X0_vld%*%beta
  pai_vld=1/(1+exp(lim_vld))
  y_vld=rbinom(n_vld,size=1,prob=pai_vld)
  X_vld=cbind(X0_vld,y_vld)
  summary(glm(y_response~X0_predictor[,1:20], family=binomial, maxit=50))
  ###########################################################################
  ###########################################################################
  n.grp=4; n.mod=3;
  n.r=matrix(0,4,n.grp)
  n.r[1,]=c(1,1,1,0)/3
  n1=0;n2=0;n3=0;n4=0
  q1=1
  n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);
  ##n1_vld=round(n_vld*n.r[q1,1]);n2_vld=round(n_vld*n.r[q1,2]);n3_vld=round(n_vld*n.r[q1,3]);n4_vld=n_vld-n1_vld-n2_vld-n3_vld;
  Miss_matrix=matrix(1,n.grp,n.mod)
  Miss_matrix[1,2:3]=0;
  Miss_matrix[2,2]=0;
  Miss_matrix[3,3]=0;
  C=glm(y_response~X0_predictor[,1:p1],data=data.frame(y_response,X0_predictor),family=binomial)$coef
  A_initial=matrix(0,nrow=p1,ncol=p2+p3)
  A_update=matrix(1,nrow=p1,ncol=p2+p3)
  b1_initial=matrix(0,nrow=1,ncol=p2+p3)
  b1_update=matrix(1,nrow=1,ncol=p2+p3)
  b2_initial=matrix(0,nrow=1,ncol=p2+p3)
  b2_update=matrix(1,nrow=1,ncol=p2+p3)
  Sigma_initial=diag(rep(1,p2+p3))
  Sigma_update=diag(rep(0.5,p2+p3))
  X0_est=X0_predictor
  while(sum(abs(A_update-A_initial))>1e-7||sum(abs(b1_update-b1_initial))>1e-7||sum(abs(Sigma_update-Sigma_initial))>1e-7){
    A_initial=A_update;
    b1_initial=b1_update;
    b2_initial=b2_update;
    Sigma_initial=Sigma_update;
    ####E step####
    for(i in 1:n1){
      X0_est[i,(p1+1):p_num]=y_response[i]*(X0_est[i,1:p1]%*%A_update+b1_update)+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update+b2_update)
    }
    for(i in n1+1:n2){
      X0_est[i,(p1+1):(p1+p2)]=y_response[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2])+
        (X0_est[i,(p1+p2+1):p_num]-
           y_response[i]*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b1_update[(p2+1):(p2+p3)])-
           (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b2_update[(p2+1):(p2+p3)]))%*%
        solve(Sigma_update[(p2+1):(p2+p3),(p2+1):(p2+p3)])%*%Sigma_update[(p2+1):(p2+p3),1:p2]                                                                                        
    }
    for(i in (n1+n2)+1:n3){
      X0_est[i,(p1+p2+1):p_num]=y_response[i]*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b1_update[p2+1:p3])+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b2_update[p2+1:p3])+
        (X0_est[i,p1+1:p2]-y_response[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])-
           (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2]))%*%
        solve(Sigma_update[1:p2,1:p2])%*%Sigma_update[1:p2,(p2+1):(p2+p3)]
    }
    ####M step####
    A_update=solve(t(X0_est[,1:p1])%*%X0_est[,1:p1])%*%t(X0_est[,1:p1])%*%(X0_est[,(p1+1):p_num]-as.matrix(y_response)%*%((b1_update))-as.matrix(1-y_response)%*%((b2_update)))
    X1=X0_est[y_response==1,];X2=X0_est[y_response==0,]
    b1_update=apply(X1[,(p1+1):p_num]-X1[,1:p1]%*%A_update,2,mean)
    b2_update=apply(X2[,(p1+1):p_num]-X2[,1:p1]%*%A_update,2,mean)
    Temp=matrix(0,p2+p3,p2+p3)
    Penal=matrix(0,p2+p3,p2+p3)
    for(i in 1:n){
      Temp=Temp+y_response[i]*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)+
        (1-y_response[i])*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)
    }
    Penal[1:p2,1:p2]=n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:p3]%*%solve(Sigma_initial[p2+1:p3,p2+1:p3])%*%Sigma_initial[p2+1:p3,1:p2])
    Penal[p2+1:p3,p2+1:p3]=n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,1:p2]%*%solve(Sigma_initial[1:p2,1:p2])%*%Sigma_initial[1:p2,p2+1:p3])
    Sigma_update=(Temp+n1*Sigma_initial+Penal)/n
  }
  ###############estimation####################
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
  eta4=X0_vld[,1:p1]%*%beta4_1est+X0_vld[,p1+1:(p2+p3)]%*%beta4_23est+rep(int4_est,n_vld)
  y4_predict=1/(1+exp(-eta4))
  y_predict=c(y4_predict)
  
  
  pred1=prediction(y_predict, y_vld)
  perf1=performance(pred1,"tpr","fpr")
  auc.tmp1=performance(pred1,"auc")
  auc1=as.numeric(auc.tmp1@y.values)
  ci.auc(y_vld,y_predict)
  auc1
  AUC=c(AUC, auc1)
}
AUC
savePlot <- function(myPlot,number) {
  png(paste("C:/Users/xliu203/Desktop/compare performance plots/myPlot",number,".png"),width =746, height = 552)
  print(myPlot)
  dev.off()
}

x=c(0.2,0.4,0.6,0.8)
error=as.vector(t(AUC))
Method=c(rep("EM",4),rep("Sep",4))
mydata=data.frame(x,error,Method)
p=ggplot(mydata,aes(x=x,y=error,color=Method))+geom_line(size=0.8)
p_1=p+xlab("Missing Percentage")+ylab("Area Under the Curve")+ggtitle("Performance Comparison")
savePlot(p_1,9)
plot(p_1)

#######################################################################################################################
########################################Different correlation##########################################################
AUC2=NULL
for(l in 1:100){
  n=300;p1=10;p2=5;p3=5;p_num=p1+p2+p3
  Sigma_X1=matrix(0.6,p1,p1)
  diag(Sigma_X1)=1
  Sigma_X2=matrix(0.6,p2,p2)
  diag(Sigma_X2)=1
  Sigma_X3=matrix(0.6,p3,p3)
  diag(Sigma_X3)=1
  Sigma_X=matrix(0.6,p_num+1,p_num+1)
  Sigma_X[1:p1,1:p1]=Sigma_X1
  Sigma_X[(p1+1):(p1+p2),(p1+1):(p1+p2)]=Sigma_X2
  Sigma_X[(p1+p2+1):(p_num),(p1+p2+1):(p_num)]=Sigma_X3
  Sigma_X[p_num+1,p_num+1]=1
  Sigma_X[p_num+1,1:p_num]=0
  Sigma_X[1:p_num,p_num+1]=0
  beta=c(rep(0.2,p1),rep(0.2,p2),rep(0.2,p3),1)
  mu=c(rep(0,p1), rep(0,p2),rep(0,p3),0)
  X0_predictor=NULL
  for(i in 1:n){
    X0_predictor=rbind(X0_predictor,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  lim=X0_predictor%*%beta
  pai=1/(1+exp(lim))
  
  y_response=rbinom(n,size=1,prob=pai)
  
  X0=cbind(X0_predictor,y_response)
  Truesigma=Sigma_X[(p1+1):(p_num),(p1+1):(p_num)]-Sigma_X[(p1+1):(p_num),1:p1]%*%solve(Sigma_X[1:p1,1:p1])%*%Sigma_X[1:p1,(p1+1):(p_num)]
  ############################################################################
  ####################validation set##########################################
  n_vld=100
  X0_vld=NULL
  for(i in 1:n_vld){
    X0_vld=rbind(X0_vld,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  lim_vld=X0_vld%*%beta
  pai_vld=1/(1+exp(lim_vld))
  y_vld=rbinom(n_vld,size=1,prob=pai_vld)
  X_vld=cbind(X0_vld,y_vld)
  summary(glm(y_response~X0_predictor[,1:20], family=binomial, maxit=50))
  ###########################################################################
  ###########################################################################
  n.grp=4; n.mod=3;
  n.r=matrix(0,4,n.grp)
  n.r[1,]=c(1,1,1,0)/3
  n1=0;n2=0;n3=0;n4=0
  q1=1
  n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);
  ##n1_vld=round(n_vld*n.r[q1,1]);n2_vld=round(n_vld*n.r[q1,2]);n3_vld=round(n_vld*n.r[q1,3]);n4_vld=n_vld-n1_vld-n2_vld-n3_vld;
  Miss_matrix=matrix(1,n.grp,n.mod)
  Miss_matrix[1,2:3]=0;
  Miss_matrix[2,2]=0;
  Miss_matrix[3,3]=0;
  C=glm(y_response~X0_predictor[,1:p1],data=data.frame(y_response,X0_predictor),family=binomial)$coef
  A_initial=matrix(0,nrow=p1,ncol=p2+p3)
  A_update=matrix(1,nrow=p1,ncol=p2+p3)
  b1_initial=matrix(0,nrow=1,ncol=p2+p3)
  b1_update=matrix(1,nrow=1,ncol=p2+p3)
  b2_initial=matrix(0,nrow=1,ncol=p2+p3)
  b2_update=matrix(1,nrow=1,ncol=p2+p3)
  Sigma_initial=diag(rep(1,p2+p3))
  Sigma_update=diag(rep(0.5,p2+p3))
  X0_est=X0_predictor
  while(sum(abs(A_update-A_initial))>1e-7||sum(abs(b1_update-b1_initial))>1e-7||sum(abs(Sigma_update-Sigma_initial))>1e-7){
    A_initial=A_update;
    b1_initial=b1_update;
    b2_initial=b2_update;
    Sigma_initial=Sigma_update;
    ####E step####
    for(i in 1:n1){
      X0_est[i,(p1+1):p_num]=y_response[i]*(X0_est[i,1:p1]%*%A_update+b1_update)+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update+b2_update)
    }
    for(i in n1+1:n2){
      X0_est[i,(p1+1):(p1+p2)]=y_response[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2])+
        (X0_est[i,(p1+p2+1):p_num]-
           y_response[i]*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b1_update[(p2+1):(p2+p3)])-
           (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b2_update[(p2+1):(p2+p3)]))%*%
        solve(Sigma_update[(p2+1):(p2+p3),(p2+1):(p2+p3)])%*%Sigma_update[(p2+1):(p2+p3),1:p2]                                                                                        
    }
    for(i in (n1+n2)+1:n3){
      X0_est[i,(p1+p2+1):p_num]=y_response[i]*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b1_update[p2+1:p3])+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b2_update[p2+1:p3])+
        (X0_est[i,p1+1:p2]-y_response[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])-
           (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2]))%*%
        solve(Sigma_update[1:p2,1:p2])%*%Sigma_update[1:p2,(p2+1):(p2+p3)]
    }
    ####M step####
    A_update=solve(t(X0_est[,1:p1])%*%X0_est[,1:p1])%*%t(X0_est[,1:p1])%*%(X0_est[,(p1+1):p_num]-as.matrix(y_response)%*%((b1_update))-as.matrix(1-y_response)%*%((b2_update)))
    X1=X0_est[y_response==1,];X2=X0_est[y_response==0,]
    b1_update=apply(X1[,(p1+1):p_num]-X1[,1:p1]%*%A_update,2,mean)
    b2_update=apply(X2[,(p1+1):p_num]-X2[,1:p1]%*%A_update,2,mean)
    Temp=matrix(0,p2+p3,p2+p3)
    Penal=matrix(0,p2+p3,p2+p3)
    for(i in 1:n){
      Temp=Temp+y_response[i]*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)+
        (1-y_response[i])*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)
    }
    Penal[1:p2,1:p2]=n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:p3]%*%solve(Sigma_initial[p2+1:p3,p2+1:p3])%*%Sigma_initial[p2+1:p3,1:p2])
    Penal[p2+1:p3,p2+1:p3]=n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,1:p2]%*%solve(Sigma_initial[1:p2,1:p2])%*%Sigma_initial[1:p2,p2+1:p3])
    Sigma_update=(Temp+n1*Sigma_initial+Penal)/n
  }
  ###############estimation####################
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
  eta4=X0_vld[,1:p1]%*%beta4_1est+X0_vld[,p1+1:(p2+p3)]%*%beta4_23est+rep(int4_est,n_vld)
  y4_predict=1/(1+exp(-eta4))
  y_predict=c(y4_predict)
  
  
  pred1=prediction(y_predict, y_vld)
  perf1=performance(pred1,"tpr","fpr")
  auc.tmp1=performance(pred1,"auc")
  auc1=as.numeric(auc.tmp1@y.values)
  ci.auc(y_vld,y_predict)
  auc1
  AUC2=c(AUC2, auc1)
}
AUC2
#######################################################################################################################
#############################################Different sample size#####################################################
AUC3=NULL
for(l in 1:100){
  n=150;p1=10;p2=5;p3=5;p_num=p1+p2+p3
  Sigma_X1=matrix(0.6,p1,p1)
  diag(Sigma_X1)=1
  Sigma_X2=matrix(0.6,p2,p2)
  diag(Sigma_X2)=1
  Sigma_X3=matrix(0.6,p3,p3)
  diag(Sigma_X3)=1
  Sigma_X=matrix(0.6,p_num+1,p_num+1)
  Sigma_X[1:p1,1:p1]=Sigma_X1
  Sigma_X[(p1+1):(p1+p2),(p1+1):(p1+p2)]=Sigma_X2
  Sigma_X[(p1+p2+1):(p_num),(p1+p2+1):(p_num)]=Sigma_X3
  Sigma_X[p_num+1,p_num+1]=1
  Sigma_X[p_num+1,1:p_num]=0
  Sigma_X[1:p_num,p_num+1]=0
  beta=c(rep(0.2,p1),rep(0.2,p2),rep(0.2,p3),1)
  mu=c(rep(0,p1), rep(0,p2),rep(0,p3),0)
  X0_predictor=NULL
  for(i in 1:n){
    X0_predictor=rbind(X0_predictor,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  lim=X0_predictor%*%beta
  pai=1/(1+exp(lim))
  
  y_response=rbinom(n,size=1,prob=pai)
  
  X0=cbind(X0_predictor,y_response)
  Truesigma=Sigma_X[(p1+1):(p_num),(p1+1):(p_num)]-Sigma_X[(p1+1):(p_num),1:p1]%*%solve(Sigma_X[1:p1,1:p1])%*%Sigma_X[1:p1,(p1+1):(p_num)]
  ############################################################################
  ####################validation set##########################################
  n_vld=100
  X0_vld=NULL
  for(i in 1:n_vld){
    X0_vld=rbind(X0_vld,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  lim_vld=X0_vld%*%beta
  pai_vld=1/(1+exp(lim_vld))
  y_vld=rbinom(n_vld,size=1,prob=pai_vld)
  X_vld=cbind(X0_vld,y_vld)
  summary(glm(y_response~X0_predictor[,1:20], family=binomial, maxit=50))
  ###########################################################################
  ###########################################################################
  n.grp=4; n.mod=3;
  n.r=matrix(0,4,n.grp)
  n.r[1,]=c(1,1,1,0)/3
  n1=0;n2=0;n3=0;n4=0
  q1=1
  n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);
  ##n1_vld=round(n_vld*n.r[q1,1]);n2_vld=round(n_vld*n.r[q1,2]);n3_vld=round(n_vld*n.r[q1,3]);n4_vld=n_vld-n1_vld-n2_vld-n3_vld;
  Miss_matrix=matrix(1,n.grp,n.mod)
  Miss_matrix[1,2:3]=0;
  Miss_matrix[2,2]=0;
  Miss_matrix[3,3]=0;
  C=glm(y_response~X0_predictor[,1:p1],data=data.frame(y_response,X0_predictor),family=binomial)$coef
  A_initial=matrix(0,nrow=p1,ncol=p2+p3)
  A_update=matrix(1,nrow=p1,ncol=p2+p3)
  b1_initial=matrix(0,nrow=1,ncol=p2+p3)
  b1_update=matrix(1,nrow=1,ncol=p2+p3)
  b2_initial=matrix(0,nrow=1,ncol=p2+p3)
  b2_update=matrix(1,nrow=1,ncol=p2+p3)
  Sigma_initial=diag(rep(1,p2+p3))
  Sigma_update=diag(rep(0.5,p2+p3))
  X0_est=X0_predictor
  while(sum(abs(A_update-A_initial))>1e-7||sum(abs(b1_update-b1_initial))>1e-7||sum(abs(Sigma_update-Sigma_initial))>1e-7){
    A_initial=A_update;
    b1_initial=b1_update;
    b2_initial=b2_update;
    Sigma_initial=Sigma_update;
    ####E step####
    for(i in 1:n1){
      X0_est[i,(p1+1):p_num]=y_response[i]*(X0_est[i,1:p1]%*%A_update+b1_update)+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update+b2_update)
    }
    for(i in n1+1:n2){
      X0_est[i,(p1+1):(p1+p2)]=y_response[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2])+
        (X0_est[i,(p1+p2+1):p_num]-
           y_response[i]*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b1_update[(p2+1):(p2+p3)])-
           (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,(p2+1):(p2+p3)]+b2_update[(p2+1):(p2+p3)]))%*%
        solve(Sigma_update[(p2+1):(p2+p3),(p2+1):(p2+p3)])%*%Sigma_update[(p2+1):(p2+p3),1:p2]                                                                                        
    }
    for(i in (n1+n2)+1:n3){
      X0_est[i,(p1+p2+1):p_num]=y_response[i]*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b1_update[p2+1:p3])+
        (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,p2+1:p3]+b2_update[p2+1:p3])+
        (X0_est[i,p1+1:p2]-y_response[i]*(X0_est[i,1:p1]%*%A_update[,1:p2]+b1_update[1:p2])-
           (1-y_response[i])*(X0_est[i,1:p1]%*%A_update[,1:p2]+b2_update[1:p2]))%*%
        solve(Sigma_update[1:p2,1:p2])%*%Sigma_update[1:p2,(p2+1):(p2+p3)]
    }
    ####M step####
    A_update=solve(t(X0_est[,1:p1])%*%X0_est[,1:p1])%*%t(X0_est[,1:p1])%*%(X0_est[,(p1+1):p_num]-as.matrix(y_response)%*%((b1_update))-as.matrix(1-y_response)%*%((b2_update)))
    X1=X0_est[y_response==1,];X2=X0_est[y_response==0,]
    b1_update=apply(X1[,(p1+1):p_num]-X1[,1:p1]%*%A_update,2,mean)
    b2_update=apply(X2[,(p1+1):p_num]-X2[,1:p1]%*%A_update,2,mean)
    Temp=matrix(0,p2+p3,p2+p3)
    Penal=matrix(0,p2+p3,p2+p3)
    for(i in 1:n){
      Temp=Temp+y_response[i]*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b1_update)+
        (1-y_response[i])*t(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)%*%(X0_est[i,(p1+1):p_num]-X0_est[i,1:p1]%*%A_update-b2_update)
    }
    Penal[1:p2,1:p2]=n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:p3]%*%solve(Sigma_initial[p2+1:p3,p2+1:p3])%*%Sigma_initial[p2+1:p3,1:p2])
    Penal[p2+1:p3,p2+1:p3]=n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,1:p2]%*%solve(Sigma_initial[1:p2,1:p2])%*%Sigma_initial[1:p2,p2+1:p3])
    Sigma_update=(Temp+n1*Sigma_initial+Penal)/n
  }
  ###############estimation####################
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
  eta4=X0_vld[,1:p1]%*%beta4_1est+X0_vld[,p1+1:(p2+p3)]%*%beta4_23est+rep(int4_est,n_vld)
  y4_predict=1/(1+exp(-eta4))
  y_predict=c(y4_predict)
  
  
  pred1=prediction(y_predict, y_vld)
  perf1=performance(pred1,"tpr","fpr")
  auc.tmp1=performance(pred1,"auc")
  auc1=as.numeric(auc.tmp1@y.values)
  ci.auc(y_vld,y_predict)
  auc1
  AUC3=c(AUC3, auc1)
}
AUC3