library(mvtnorm)
library(Matrix)
library(clusterGeneration)
library(ggplot2)
library(ROCR)
library(glmnet)
library("pROC")
#######################################################################################
###################Adding covariates#####################################################
MSE=NULL
for (l in 1:100){
  n=300;p1=10;p2=5;p3=5;p_num=p1+p2+p3
  Sigma_X1=matrix(0.6,p1,p1)
  diag(Sigma_X1)=1
  Sigma_X2=matrix(0.6,p2,p2)
  diag(Sigma_X2)=1
  Sigma_X3=matrix(0.6,p3,p3)
  diag(Sigma_X3)=1
  Sigma_X=matrix(0.6,p_num,p_num)
  Sigma_X[1:p1,1:p1]=Sigma_X1
  Sigma_X[(p1+1):(p1+p2),(p1+1):(p1+p2)]=Sigma_X2
  Sigma_X[(p1+p2+1):(p_num),(p1+p2+1):(p_num)]=Sigma_X3
  beta=c(rep(0.2,p1),rep(0.2,p2),rep(0.2,p3))
  mu=c(rep(0,p1), rep(0,p2),rep(0,p3))
  X0_predictor=NULL
  for(i in 1:n){
    X0_predictor=rbind(X0_predictor,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  intcpt=2
  std=1
  y_response=X0_predictor%*%beta+intcpt+rnorm(n,mean=0,sd=std)
  X0=cbind(X0_predictor,y_response)
  ####################################################################
  ####################Validation data#################################
  n_vld=100
  X0_vld=rmvnorm(n_vld, mean=mu, sigma=Sigma_X, method="eigen")
  y_vld=X0_vld%*%beta+intcpt+rnorm(n_vld,mean=0,sd=std)
  X0_validation=cbind(X0_vld,y_vld)
  ####################################################################
  ####################################################################
  n.grp=4; n.mod=3;
  n.r=matrix(0,4,n.grp)
  n.r[1,]=c(1,1,1,0)/3
  n1=0;n2=0;n3=0;n4=0
  q1=1
  n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);
  
  #####################################################################
  ####################ECM estimation on regression#####################
  A_initial=matrix(0,p1+1,p2+p3+1)
  A2_initial=matrix(0,p1,p2)
  b2_initial=matrix(0,1,p2)
  A3_initial=matrix(0,p1,p3)
  b3_initial=matrix(0,1,p3)
  Ay_initial=matrix(0,p1,1)
  by_initial=matrix(0,1,1)
  Sigma_initial=diag(rep(0.5,p2+p3+1))
  A_update=matrix(1,p1+1,p2+p3+1)
  A2_update=matrix(1,p1,p2)
  b2_update=matrix(1,1,p2)
  A3_update=matrix(1,p1,p3)
  b3_update=matrix(1,1,p3)
  Ay_update=matrix(0,p1,1)
  by_update=matrix(0,1,1)
  Sigma_update=diag(rep(1,p2+p3+1))
  C=array(0,dim=(c(n,p2+p3+1,p2+p3+1)))
  X0_est=X0
  while(sum(abs(A_initial-A_update))>1e-7||sum(abs(Sigma_initial-Sigma_update))>1e-7){
    A_initial=A_update
    A2_initial=A2_update
    b2_initial=b2_update
    A3_initial=A3_update
    b3_initial=b3_update
    Ay_initial=Ay_update
    by_initial=by_update
    Sigma_initial=Sigma_update;
    ####E step#### 
    X0_est[1:n1,(p1+1):p_num]=X0_est[1:n1,1:p1]%*%cbind(A2_update, A3_update)+matrix(rep(1,n1),ncol=1)%*%c(b2_update, b3_update)+(X0_est[1:n1,p_num+1]-X0_est[1:n1,1:p1]%*%Ay_update-matrix(rep(1,n1),ncol=1)%*%by_update)%*%Sigma_update[1:(p2+p3),p2+p3+1]/Sigma_update[p2+p3+1,p2+p3+1]
    X0_est[n1+1:n2,(p1+1):(p1+p2)]=X0_est[n1+1:n2,1:p1]%*%A2_update+matrix(rep(1,n1),ncol=1)%*%b2_update+(X0_est[n1+1:n2,(p1+p2+1):(p_num+1)]-X0_est[n1+1:n2,1:p1]%*%cbind(A3_update,Ay_update)-matrix(rep(1,n2),ncol=1)%*%c(b3_update, by_update))%*%solve(Sigma_update[p2+1:(p3+1),p2+1:(p3+1)])%*%Sigma_update[p2+1:(p3+1),1:p2]
    X0_est[n1+n2+1:n3,(p1+p2+1):p_num]=X0_est[n1+n2+1:n3,1:p1]%*%A3_update+matrix(rep(1,n3),ncol=1)%*%b3_update+(X0_est[n1+n2+1:n3,c(p1+1:p2,p_num+1)]-X0_est[n1+n2+1:n3,1:p1]%*%cbind(A2_update,Ay_update)-matrix(rep(1,n3),ncol=1)%*%c(b2_update,by_update))%*%solve(Sigma_update[c(1:p2,p2+p3+1),c(1:p2,p2+p3+1)])%*%Sigma_update[c(1:p2,p2+p3+1),p2+1:p3]
    
    ####M step####
    A_update=solve(t(cbind(1,X0_est[,1:p1]))%*%cbind(1,X0_est[,1:p1]))%*%t(cbind(1,X0_est[,1:p1]))%*%X0_est[,(p1+1):(p_num+1)]
    A2_update=A_update[2:(p1+1),1:p2]
    b2_update=A_update[1,1:p2]
    A3_update=A_update[2:(p1+1),p2+1:p3]
    b3_update=A_update[1,p2+1:p3]
    Ay_update=A_update[2:(p1+1),p2+p3+1]
    by_update=A_update[1,p2+p3+1]
    Sigma_update=matrix(0,p2+p3+1,p2+p3+1)
    Penal=matrix(0,p2+p3+1,p2+p3+1)
    Penal[1:(p2+p3),1:(p2+p3)]=n1*(Sigma_initial[1:(p2+p3),1:(p2+p3)]-Sigma_initial[1:(p2+p3),p2+p3+1]%*%solve(Sigma_initial[p2+p3+1,p2+p3+1])%*%Sigma_initial[p2+p3+1, 1:(p2+p3)])
    Penal[1:p2,1:p2]=Penal[1:p2,1:p2]+n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:(p3+1)]%*%solve(Sigma_initial[p2+1:(p3+1),p2+1:(p3+1)])%*%Sigma_initial[p2+1:(p3+1),1:p2])
    Penal[p2+1:p3,p2+1:p3]=Penal[p2+1:p3,p2+1:p3]+n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,c(1:p2,p2+p3+1)]%*%solve(Sigma_initial[c(1:p2,p2+p3+1),c(1:p2,p2+p3+1)])%*%Sigma_initial[c(1:p2,p2+p3+1),p2+1:p3])
    Sigma_update=t(X0_est[,(p1+1):(p_num+1)]-cbind(1,X0_est[,1:p1])%*%A_update)%*%(X0_est[,(p1+1):(p_num+1)]-cbind(1,X0_est[,1:p1])%*%A_update)/n+Penal/n
    
  }
  
  beta1_est=Ay_update
  intcpt1_est=by_update
  
  beta1_1est=(Ay_update-Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%t(A3_update))
  beta1_3est=Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])
  intcpt1_1est=by_update-Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b3_update
  
  beta2_1est=(Ay_update-Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%t(A2_update))
  beta2_2est=Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])
  intcpt2_1est=by_update-Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%b2_update
  
  beta3_1est=(Ay_update-Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])%*%t(cbind(A2_update, A3_update)))
  beta3_23est=Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])
  intcpt3_1est=by_update-Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])%*%c(b2_update, b3_update)
  
  
  
  y_pred=X0_vld[,1:p1]%*%t(beta3_1est)+X0_vld[,p1+1:(p2+p3)]%*%t(beta3_23est)+rep(intcpt3_1est,n_vld)
  
  MSE=c(MSE,mean((y_vld-y_pred)^2))
}


MSE




beta1=lm(X0[(n1+1):n,p_num+1]~X0[(n1+1):n,1:p_num])$coef
beta2=lm(X0[1:n1,p_num+1]~X0[1:n1,1:p1])$coef
y1_spred=cbind(1,X0_vld[1:n1,1:p1])%*%beta2
y2_spred=cbind(1,X0_vld[(n1+1):n,])%*%beta1
MSE[2,q1]=mean((y_vld-c(y1_spred,y2_spred))^2, na.rm=TRUE)


#######################################################################################
###################Changing correlation#####################################################
MSE2=NULL
for (l in 1:100){
  n=300;p1=10;p2=5;p3=5;p_num=p1+p2+p3
  Sigma_X1=matrix(0.6,p1,p1)
  diag(Sigma_X1)=1
  Sigma_X2=matrix(0.6,p2,p2)
  diag(Sigma_X2)=1
  Sigma_X3=matrix(0.6,p3,p3)
  diag(Sigma_X3)=1
  Sigma_X=matrix(0,p_num,p_num)
  Sigma_X[1:p1,1:p1]=Sigma_X1
  Sigma_X[(p1+1):(p1+p2),(p1+1):(p1+p2)]=Sigma_X2
  Sigma_X[(p1+p2+1):(p_num),(p1+p2+1):(p_num)]=Sigma_X3
  beta=c(rep(0.2,p1),rep(0.2,p2),rep(0.2,p3))
  mu=c(rep(0,p1), rep(0,p2),rep(0,p3))
  X0_predictor=NULL
  for(i in 1:n){
    X0_predictor=rbind(X0_predictor,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  intcpt=2
  std=1
  y_response=X0_predictor%*%beta+intcpt+rnorm(n,mean=0,sd=std)
  X0=cbind(X0_predictor,y_response)
  ####################################################################
  ####################Validation data#################################
  n_vld=100
  X0_vld=rmvnorm(n_vld, mean=mu, sigma=Sigma_X, method="eigen")
  y_vld=X0_vld%*%beta+intcpt+rnorm(n_vld,mean=0,sd=std)
  X0_validation=cbind(X0_vld,y_vld)
  ####################################################################
  ####################################################################
  n.grp=4; n.mod=3;
  n.r=matrix(0,4,n.grp)
  n.r[1,]=c(1,1,1,0)/3
  n1=0;n2=0;n3=0;n4=0
  q1=1
  n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);
  
  #####################################################################
  ####################ECM estimation on regression#####################
  A_initial=matrix(0,p1+1,p2+p3+1)
  A2_initial=matrix(0,p1,p2)
  b2_initial=matrix(0,1,p2)
  A3_initial=matrix(0,p1,p3)
  b3_initial=matrix(0,1,p3)
  Ay_initial=matrix(0,p1,1)
  by_initial=matrix(0,1,1)
  Sigma_initial=diag(rep(0.5,p2+p3+1))
  A_update=matrix(1,p1+1,p2+p3+1)
  A2_update=matrix(1,p1,p2)
  b2_update=matrix(1,1,p2)
  A3_update=matrix(1,p1,p3)
  b3_update=matrix(1,1,p3)
  Ay_update=matrix(0,p1,1)
  by_update=matrix(0,1,1)
  Sigma_update=diag(rep(1,p2+p3+1))
  C=array(0,dim=(c(n,p2+p3+1,p2+p3+1)))
  X0_est=X0
  while(sum(abs(A_initial-A_update))>1e-7||sum(abs(Sigma_initial-Sigma_update))>1e-7){
    A_initial=A_update
    A2_initial=A2_update
    b2_initial=b2_update
    A3_initial=A3_update
    b3_initial=b3_update
    Ay_initial=Ay_update
    by_initial=by_update
    Sigma_initial=Sigma_update;
    ####E step#### 
    X0_est[1:n1,(p1+1):p_num]=X0_est[1:n1,1:p1]%*%cbind(A2_update, A3_update)+matrix(rep(1,n1),ncol=1)%*%c(b2_update, b3_update)+(X0_est[1:n1,p_num+1]-X0_est[1:n1,1:p1]%*%Ay_update-matrix(rep(1,n1),ncol=1)%*%by_update)%*%Sigma_update[1:(p2+p3),p2+p3+1]/Sigma_update[p2+p3+1,p2+p3+1]
    X0_est[n1+1:n2,(p1+1):(p1+p2)]=X0_est[n1+1:n2,1:p1]%*%A2_update+matrix(rep(1,n1),ncol=1)%*%b2_update+(X0_est[n1+1:n2,(p1+p2+1):(p_num+1)]-X0_est[n1+1:n2,1:p1]%*%cbind(A3_update,Ay_update)-matrix(rep(1,n2),ncol=1)%*%c(b3_update, by_update))%*%solve(Sigma_update[p2+1:(p3+1),p2+1:(p3+1)])%*%Sigma_update[p2+1:(p3+1),1:p2]
    X0_est[n1+n2+1:n3,(p1+p2+1):p_num]=X0_est[n1+n2+1:n3,1:p1]%*%A3_update+matrix(rep(1,n3),ncol=1)%*%b3_update+(X0_est[n1+n2+1:n3,c(p1+1:p2,p_num+1)]-X0_est[n1+n2+1:n3,1:p1]%*%cbind(A2_update,Ay_update)-matrix(rep(1,n3),ncol=1)%*%c(b2_update,by_update))%*%solve(Sigma_update[c(1:p2,p2+p3+1),c(1:p2,p2+p3+1)])%*%Sigma_update[c(1:p2,p2+p3+1),p2+1:p3]
    
    ####M step####
    A_update=solve(t(cbind(1,X0_est[,1:p1]))%*%cbind(1,X0_est[,1:p1]))%*%t(cbind(1,X0_est[,1:p1]))%*%X0_est[,(p1+1):(p_num+1)]
    A2_update=A_update[2:(p1+1),1:p2]
    b2_update=A_update[1,1:p2]
    A3_update=A_update[2:(p1+1),p2+1:p3]
    b3_update=A_update[1,p2+1:p3]
    Ay_update=A_update[2:(p1+1),p2+p3+1]
    by_update=A_update[1,p2+p3+1]
    Sigma_update=matrix(0,p2+p3+1,p2+p3+1)
    Penal=matrix(0,p2+p3+1,p2+p3+1)
    Penal[1:(p2+p3),1:(p2+p3)]=n1*(Sigma_initial[1:(p2+p3),1:(p2+p3)]-Sigma_initial[1:(p2+p3),p2+p3+1]%*%solve(Sigma_initial[p2+p3+1,p2+p3+1])%*%Sigma_initial[p2+p3+1, 1:(p2+p3)])
    Penal[1:p2,1:p2]=Penal[1:p2,1:p2]+n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:(p3+1)]%*%solve(Sigma_initial[p2+1:(p3+1),p2+1:(p3+1)])%*%Sigma_initial[p2+1:(p3+1),1:p2])
    Penal[p2+1:p3,p2+1:p3]=Penal[p2+1:p3,p2+1:p3]+n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,c(1:p2,p2+p3+1)]%*%solve(Sigma_initial[c(1:p2,p2+p3+1),c(1:p2,p2+p3+1)])%*%Sigma_initial[c(1:p2,p2+p3+1),p2+1:p3])
    Sigma_update=t(X0_est[,(p1+1):(p_num+1)]-cbind(1,X0_est[,1:p1])%*%A_update)%*%(X0_est[,(p1+1):(p_num+1)]-cbind(1,X0_est[,1:p1])%*%A_update)/n+Penal/n
    
  }
  
  beta1_est=Ay_update
  intcpt1_est=by_update
  
  beta1_1est=(Ay_update-Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%t(A3_update))
  beta1_3est=Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])
  intcpt1_1est=by_update-Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b3_update
  
  beta2_1est=(Ay_update-Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%t(A2_update))
  beta2_2est=Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])
  intcpt2_1est=by_update-Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%b2_update
  
  beta3_1est=(Ay_update-Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])%*%t(cbind(A2_update, A3_update)))
  beta3_23est=Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])
  intcpt3_1est=by_update-Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])%*%c(b2_update, b3_update)
  
  
  
  y_pred=X0_vld[,1:p1]%*%t(beta3_1est)+X0_vld[,p1+1:(p2+p3)]%*%t(beta3_23est)+rep(intcpt3_1est,n_vld)
  
  MSE2=c(MSE2,mean((y_vld-y_pred)^2))
}


MSE2






savePlot <- function(myPlot,number) {
  png(paste("C:/Users/xliu203/Desktop/compare performance plots/myPlot",number,".png"),width =746, height = 552)
  print(myPlot)
  dev.off()
}

x=c(0.2,0.4,0.6,0.8)
error=as.vector(t(MSE))
Method=c(rep("EM",4),rep("Sep",4))
mydata=data.frame(x,error,Method)
p=ggplot(mydata,aes(x=x,y=error,color=Method))+geom_line(size=0.8)
p1=p+xlab("Missing Percentage")+ylab("Mean Square Error")+ggtitle("Performance Comparison")
savePlot(p1,7)
plot(p1)


MSE3=NULL
for (l in 1:100){
  n=150;p1=10;p2=5;p3=5;p_num=p1+p2+p3
  Sigma_X1=matrix(0.6,p1,p1)
  diag(Sigma_X1)=1
  Sigma_X2=matrix(0.6,p2,p2)
  diag(Sigma_X2)=1
  Sigma_X3=matrix(0.6,p3,p3)
  diag(Sigma_X3)=1
  Sigma_X=matrix(0.6,p_num,p_num)
  Sigma_X[1:p1,1:p1]=Sigma_X1
  Sigma_X[(p1+1):(p1+p2),(p1+1):(p1+p2)]=Sigma_X2
  Sigma_X[(p1+p2+1):(p_num),(p1+p2+1):(p_num)]=Sigma_X3
  beta=c(rep(0.2,p1),rep(0.2,p2),rep(0.2,p3))
  mu=c(rep(0,p1), rep(0,p2),rep(0,p3))
  X0_predictor=NULL
  for(i in 1:n){
    X0_predictor=rbind(X0_predictor,rmvnorm(1, mean=mu, sigma=Sigma_X, method="eigen"))
  }
  intcpt=2
  std=1
  y_response=X0_predictor%*%beta+intcpt+rnorm(n,mean=0,sd=std)
  X0=cbind(X0_predictor,y_response)
  ####################################################################
  ####################Validation data#################################
  n_vld=100
  X0_vld=rmvnorm(n_vld, mean=mu, sigma=Sigma_X, method="eigen")
  y_vld=X0_vld%*%beta+intcpt+rnorm(n_vld,mean=0,sd=std)
  X0_validation=cbind(X0_vld,y_vld)
  ####################################################################
  ####################################################################
  n.grp=4; n.mod=3;
  n.r=matrix(0,4,n.grp)
  n.r[1,]=c(1,1,1,0)/3
  n1=0;n2=0;n3=0;n4=0
  q1=1
  n1=round(n*n.r[q1,1]);n2=round(n*n.r[q1,2]);n3=round(n*n.r[q1,3]);
  
  #####################################################################
  ####################ECM estimation on regression#####################
  A_initial=matrix(0,p1+1,p2+p3+1)
  A2_initial=matrix(0,p1,p2)
  b2_initial=matrix(0,1,p2)
  A3_initial=matrix(0,p1,p3)
  b3_initial=matrix(0,1,p3)
  Ay_initial=matrix(0,p1,1)
  by_initial=matrix(0,1,1)
  Sigma_initial=diag(rep(0.5,p2+p3+1))
  A_update=matrix(1,p1+1,p2+p3+1)
  A2_update=matrix(1,p1,p2)
  b2_update=matrix(1,1,p2)
  A3_update=matrix(1,p1,p3)
  b3_update=matrix(1,1,p3)
  Ay_update=matrix(0,p1,1)
  by_update=matrix(0,1,1)
  Sigma_update=diag(rep(1,p2+p3+1))
  C=array(0,dim=(c(n,p2+p3+1,p2+p3+1)))
  X0_est=X0
  while(sum(abs(A_initial-A_update))>1e-7||sum(abs(Sigma_initial-Sigma_update))>1e-7){
    A_initial=A_update
    A2_initial=A2_update
    b2_initial=b2_update
    A3_initial=A3_update
    b3_initial=b3_update
    Ay_initial=Ay_update
    by_initial=by_update
    Sigma_initial=Sigma_update;
    ####E step#### 
    X0_est[1:n1,(p1+1):p_num]=X0_est[1:n1,1:p1]%*%cbind(A2_update, A3_update)+matrix(rep(1,n1),ncol=1)%*%c(b2_update, b3_update)+(X0_est[1:n1,p_num+1]-X0_est[1:n1,1:p1]%*%Ay_update-matrix(rep(1,n1),ncol=1)%*%by_update)%*%Sigma_update[1:(p2+p3),p2+p3+1]/Sigma_update[p2+p3+1,p2+p3+1]
    X0_est[n1+1:n2,(p1+1):(p1+p2)]=X0_est[n1+1:n2,1:p1]%*%A2_update+matrix(rep(1,n1),ncol=1)%*%b2_update+(X0_est[n1+1:n2,(p1+p2+1):(p_num+1)]-X0_est[n1+1:n2,1:p1]%*%cbind(A3_update,Ay_update)-matrix(rep(1,n2),ncol=1)%*%c(b3_update, by_update))%*%solve(Sigma_update[p2+1:(p3+1),p2+1:(p3+1)])%*%Sigma_update[p2+1:(p3+1),1:p2]
    X0_est[n1+n2+1:n3,(p1+p2+1):p_num]=X0_est[n1+n2+1:n3,1:p1]%*%A3_update+matrix(rep(1,n3),ncol=1)%*%b3_update+(X0_est[n1+n2+1:n3,c(p1+1:p2,p_num+1)]-X0_est[n1+n2+1:n3,1:p1]%*%cbind(A2_update,Ay_update)-matrix(rep(1,n3),ncol=1)%*%c(b2_update,by_update))%*%solve(Sigma_update[c(1:p2,p2+p3+1),c(1:p2,p2+p3+1)])%*%Sigma_update[c(1:p2,p2+p3+1),p2+1:p3]
    
    ####M step####
    A_update=solve(t(cbind(1,X0_est[,1:p1]))%*%cbind(1,X0_est[,1:p1]))%*%t(cbind(1,X0_est[,1:p1]))%*%X0_est[,(p1+1):(p_num+1)]
    A2_update=A_update[2:(p1+1),1:p2]
    b2_update=A_update[1,1:p2]
    A3_update=A_update[2:(p1+1),p2+1:p3]
    b3_update=A_update[1,p2+1:p3]
    Ay_update=A_update[2:(p1+1),p2+p3+1]
    by_update=A_update[1,p2+p3+1]
    Sigma_update=matrix(0,p2+p3+1,p2+p3+1)
    Penal=matrix(0,p2+p3+1,p2+p3+1)
    Penal[1:(p2+p3),1:(p2+p3)]=n1*(Sigma_initial[1:(p2+p3),1:(p2+p3)]-Sigma_initial[1:(p2+p3),p2+p3+1]%*%solve(Sigma_initial[p2+p3+1,p2+p3+1])%*%Sigma_initial[p2+p3+1, 1:(p2+p3)])
    Penal[1:p2,1:p2]=Penal[1:p2,1:p2]+n2*(Sigma_initial[1:p2,1:p2]-Sigma_initial[1:p2,p2+1:(p3+1)]%*%solve(Sigma_initial[p2+1:(p3+1),p2+1:(p3+1)])%*%Sigma_initial[p2+1:(p3+1),1:p2])
    Penal[p2+1:p3,p2+1:p3]=Penal[p2+1:p3,p2+1:p3]+n3*(Sigma_initial[p2+1:p3,p2+1:p3]-Sigma_initial[p2+1:p3,c(1:p2,p2+p3+1)]%*%solve(Sigma_initial[c(1:p2,p2+p3+1),c(1:p2,p2+p3+1)])%*%Sigma_initial[c(1:p2,p2+p3+1),p2+1:p3])
    Sigma_update=t(X0_est[,(p1+1):(p_num+1)]-cbind(1,X0_est[,1:p1])%*%A_update)%*%(X0_est[,(p1+1):(p_num+1)]-cbind(1,X0_est[,1:p1])%*%A_update)/n+Penal/n
    
  }
  
  beta1_est=Ay_update
  intcpt1_est=by_update
  
  beta1_1est=(Ay_update-Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%t(A3_update))
  beta1_3est=Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])
  intcpt1_1est=by_update-Sigma_update[p2+p3+1,p2+1:p3]%*%solve(Sigma_update[p2+1:p3,p2+1:p3])%*%b3_update
  
  beta2_1est=(Ay_update-Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%t(A2_update))
  beta2_2est=Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])
  intcpt2_1est=by_update-Sigma_update[p2+p3+1,1:p2]%*%solve(Sigma_update[1:p2,1:p2])%*%b2_update
  
  beta3_1est=(Ay_update-Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])%*%t(cbind(A2_update, A3_update)))
  beta3_23est=Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])
  intcpt3_1est=by_update-Sigma_update[p2+p3+1,1:(p2+p3)]%*%solve(Sigma_update[1:(p2+p3),1:(p2+p3)])%*%c(b2_update, b3_update)
  
  
  
  y_pred=X0_vld[,1:p1]%*%t(beta3_1est)+X0_vld[,p1+1:(p2+p3)]%*%t(beta3_23est)+rep(intcpt3_1est,n_vld)
  
  MSE3=c(MSE3,mean((y_vld-y_pred)^2))
}


MSE3