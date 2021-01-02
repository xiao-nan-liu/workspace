#############################################################
################Generate DATA################################
setwd("C:/Users/xliu203/Desktop/new Vector")
rec=read.csv("C:/Users/xliu203/Desktop/rectum volume.csv")
fileName <- dir()
datalist <- vector("list",79)
recvol=rec$X*rec$Y*rec$Z
recvol=recvol[-c(8,13,34,40,43,64,74)]
for(i in 1:79){
  datalist[[i]]=read.table(fileName[i],header=T)
}

D0=matrix(0,79,101)
for(i in 1:79){
  for(j in 0:100){
    d=datalist[[i]][[1]]
    D0[i,j+1]=quantile(datalist[[i]][[1]],0.01*j)
  }
}
name=rep(0,101)
for(i in 1:101){
  name[i]=paste("D", i-1, sep="")
}
colnames(D0)=name
mean_dose=apply(D0, 2, mean)
cov_dose=cov(D0)

N=1000; ####number of patients
genDVH=function(N){
  x=cov_dose
  DVH=abs(matrix(rmvnorm(N, mean=mean_dose, sigma=x, method="svd"), nrow=N, ncol=101))
  return(DVH)
}
DVH=genDVH(N)

rho=0.2;
beta_ps=c(rnorm(2,0.6,0.1), rnorm(2,-0.6,0.1), rep(0,10))
sig=c(1, rho^seq(1:(length(beta_ps)-1)))
sig=toeplitz(sig)
X_ps=rmvnorm(N, mean=rep(0,14),sigma=sig)
alpha=0.3;
n0=0.10;
D=NULL
for(i in 1:100){
  D=cbind(D,(DVH[,i]+DVH[,i+1])/2)
}
Geud=rep(0,100)
for(i in 1:N){
  Geud[i]=sum(0.01*D[i,]^(1/n0))^n0
}
mean(pnorm(X_ps%*%beta_ps+alpha*Geud-18.277135))
y_ps=rep(0,N)
for(i in 1:N){
  y_ps[i]=rbinom(1,1,pnorm(X_ps[i,]%*%beta_ps+alpha*Geud[i]-18.277135))  ####ratio approximately 1:4####
}
sum(y_ps)
prob=glm(y_ps~., family=binomial(link="probit"), data=data.frame(X_ps[,1:4], y_ps, Geud))
summary(prob)
################################################################################
########################sample ratio 2:3########################################
int_1=which(y_ps==1)
int_0=which(y_ps==0)
int_sample=c(sample(int_1)[1:40],sample(int_0)[1:60])
y_sample=y_ps[int_sample]
X_sample=X_ps[int_sample,]
D_sample=D[int_sample,]
prob=glm(y_sample~., family=binomial(link="probit"), data=data.frame(X_sample[,1:4], y_sample, Geud[int_sample]))
summary(prob)
###############################################################################
########################Test data##############################################
N_test=100; ####number of patients
genDVH=function(N){
  x=cov_dose
  DVH=abs(matrix(rmvnorm(N, mean=mean_dose, sigma=x, method="svd"), nrow=N, ncol=101))
  return(DVH)
}
DVH_test=genDVH(N_test)

sig=c(1, rho^seq(1:(length(beta_ps)-1)))
sig=toeplitz(sig)
X_test=rmvnorm(N_test, mean=rep(0,14),sigma=sig)
D_test=NULL
for(i in 1:100){
  D_test=cbind(D_test,(DVH_test[,i]+DVH_test[,i+1])/2)
}
Geud_test=rep(0,100)
for(i in 1:N_test){
  Geud_test[i]=sum(0.01*D_test[i,]^(1/n0))^n0
}
mean(pnorm(X_test%*%beta_ps+alpha*Geud_test-18.277135))
y_test=rep(0,N_test)
for(i in 1:N_test){
  y_test[i]=rbinom(1,1,pnorm(X_test[i,]%*%beta_ps+alpha*Geud_test[i]-18.277135))  ####ratio approximately 1:4####
}
sum(y_test)

###############################################################################
############################LKB model using sample data########################
dntcpl_sp <- function(param){
  
  
  TD_50=param[1];
  m0=param[2];
  n=param[3];
  lntcp=0;
  if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1){
    
    for(i in 1:100){
      
      d=D[i,]
      eud=(sum(d^(1/n))/length(d))^(n)
      t=(eud-TD_50)/(m0*TD_50)
      lntcp=lntcp+y_ps[i]*log(pnorm(t))+(1-y_ps[i])*log(1-pnorm(t))
    }
  }else {lntcp=-100}
  return(lntcp)
}
dntcpl_sp(c(60,0.1,0.1))
model=maxLik(dntcpl_sp,start=c(60,0.1,0.1),method="nm")
summary(model)
coef=coef(model)
theta.start <- coef
md<-function(param) -dntcpl_sp(param)
out <- nlm(md, theta.start,hessian=TRUE)
theta.hat <- out$estimate
fish=out$hessian
theta.hat
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
inv.fish <- solve(fish)
theta.hat[1] + c(-1, 1) * crit * sqrt(inv.fish[1, 1])/sqrt(100)
theta.hat[2] + c(-1, 1) * crit * sqrt(inv.fish[2, 2])/sqrt(100)
ni=theta.hat[3] + c(-1, 1) * crit * sqrt(inv.fish[3, 3])/sqrt(100)

###############################################################################
##########################mLKB model###########################################
eps=0.001
lambda.max=100
lambda.min=eps*lambda.max
lambda=seq(from=log(lambda.max), to=log(lambda.min), length.out=100)
lambda=exp(lambda)
tau=0.2
mlikelihood=function(X, y, beta,tau){
  
  w1=tau/mean(y)
  w0=(1-tau)/(1-mean(y))
  X1=X
  llh=rep(0, nrow(X))
  for(i in 1:nrow(X)){
    if(pnorm(X1[i,]%*%beta)>1-1e-7){llh[i]=w0*(1-y[i])*(-16.12)}
    else if(pnorm(X1[i,]%*%beta)<1e-7){llh[i]=w1*y[i]*(-16.12)}
    else{llh[i]=w1*y[i]*log(pnorm(X1[i,]%*%beta))+w0*(1-y[i])*log(1-pnorm(X1[i,]%*%beta))}
  }
  
  return(sum(llh))
}

mglmnet_probit=function(X, y, penalty.factor,lambda, tau){
  c=(1-tau)*mean(y)/tau/(1-mean(y))
  w1=tau/mean(y)
  w0=(1-tau)/(1-mean(y))
  X=cbind(1,X)
  beta.initial=rep(0, ncol(X))
  beta0=beta.initial
  eps=0.001
  penalty.factor=c(0,penalty.factor)
  beta=matrix(0,ncol=length(lambda), nrow=ncol(X))
  deviance=rep(0,length(lambda))
  BIC=rep(0,length(lambda))
  AIC=rep(0,length(lambda))
  AICC=rep(0,length(lambda))
  for(i in 1:length(lambda)){
    s=lambda[i]
    beta.hat=beta.initial
    beta.final=c(1, rep(0.5, ncol(X)-1))
    while(sum(abs(beta.hat-beta.final))>eps){
      beta.final=beta.hat
      z=X%*%beta.hat+(w1*y-w0*pnorm(X%*%beta.hat)-(w1-w0)*y*pnorm(X%*%beta.hat))/dnorm(X%*%beta.hat)/(w0*pnorm(X%*%beta.hat)+w1*(1-pnorm(X%*%beta.hat)))
      w=dnorm(X%*%beta.hat)^2*(w0*pnorm(X%*%beta.hat)+w1*(1-pnorm(X%*%beta.hat)))/pnorm(X%*%beta.hat)/(1-pnorm(X%*%beta.hat))
      for(j in 1:nrow(X)){
        if(dnorm(X[j,]%*%beta.hat)<=1e-3){
          z[j]=X[j,]%*%beta.hat+1e-3
          w[j]=1e-3*w1
        }
      }
      
      beta.temp=rep(1,length(beta.final))
      while(sum(abs(beta.hat-beta.temp))>eps){
        beta.temp=beta.hat
        for(index in 1:ncol(X)){
          beta.lasso=sum(w*X[,index]*(z-X[,-index]%*%beta.hat[-index]))
          if(penalty.factor[index]==1){
            if(beta.lasso>0&&s<abs(beta.lasso)) beta.hat[index]=(beta.lasso-s)/sum(w*X[,index]^2)
            if(beta.lasso<0&&s<abs(beta.lasso)) beta.hat[index]=(beta.lasso+s)/sum(w*X[,index]^2)
            if(s>=abs(beta.lasso)) beta.hat[index]=0
            
          }
          if(penalty.factor[index]==0){
            beta.hat[index]=beta.lasso/sum(w*X[,index]^2)
          }
        }##end for
        
      }##end while
    }            
    beta.initial=beta.final
    beta[,i]=beta.final
    deviance[i]=-2*(mlikelihood(X,y,beta0,tau)-mlikelihood(X,y, beta.final,tau))
    AIC=-2*mlikelihood(X,y,beta.final,tau)+sum(beta.final!=0)*2
    BIC[i]=-2*mlikelihood(X,y,beta.final,tau)+sum(beta.final!=0)*log(nrow(X))
    AICC[i]=-2*mlikelihood(X,y,beta.final,tau)+sum(beta.final!=0)*2+2*sum(beta.final!=0)*
      (sum(beta.final!=0)+1)/(nrow(X)-1-sum(beta.final!=0))
  }
  result=list(lambda=lambda, beta=beta, deviance=deviance, BIC=BIC, AICC=AICC)
  return(result)
}
mglmnet_probit(X_sample, y_sample, penalty.factor=rep(1,14),lambda,tau=0.2)

predict_glmnet_probit=function(X, y,newx, penalty.factor, lambda){
  
  object1=glmnet_probit(X,y,penalty.factor,lambda)
  object2=cv_glmnet_probit(X,y,penalty.factor,lambda)
  lambda=object2[[2]]
  cv=object2[[1]]
  k=which(cv==min(cv))[1]
  s=lambda[k]
  beta=object1[[2]][,k]
  predict=pnorm(cbind(1,newx)%*%beta)
  return(list(predict,cv[k]))
}

predict_glmnet_probit(X,y,newx,rep(1,11),lambda)
###############################################################
##############Simulation mLKB##################################
nc=seq(0, 1, length.out=50)
geud_sp=rep(0,100)
BIC=rep(0,length(nc))
AIC=rep(0,length(nc))
AICC=rep(0,length(nc))
Beta=matrix(0,ncol=50, nrow=16)

for(k in 1:length(nc)){
  
  n=nc[k]
  for(i in 1:100){
    geud_sp[i]=(sum(D_sample[i,]^(1/n))*0.01)^(n)  ##calculate gEUD
  }
  
  Xa=cbind(X_sample,geud_sp);
  pre=mglmnet_probit(Xa, y_sample,c(rep(1,14),0), lambda,tau)
  j=which.min(pre$BIC)[1]
  BIC[k]=pre$BIC[j]
  Beta[,k]=pre$beta[,j]
}

k=which.min(BIC)
which.min(AIC)
which.min(AICC)

n=nc[k]
for(i in 1:100){
  geud_sp[i]=(sum(D_sample[i,]^(1/n))*0.01)^(n)  ##calculate gEUD
}
Xa=cbind(X_sample,geud_sp);
X1=Xa[,c(1:4,15)]
Data=data.frame(X1,y_sample)
predict=rep(0,100)
for(j in 1:100){
  Data1=Data[-j,]
  prob1=glm(y_sample~., family=binomial(link="probit"), data=Data1)
  beta=coef(prob1)
  X2=as.matrix(cbind(1,Data1[,1:5]))
  eta=X2%*%beta
  w=dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta))*(w0*pnorm(eta)+w1*(1-pnorm(eta)))
  Hes=matrix(0,nrow=6, ncol=6)
  for(i in 1:99){
    Hes=Hes+w[i]*X2[i,]%*%t(X2[i,])
  }
  Q=(X2)%*%solve(Hes)%*%t(X2)
  ksi=eta/2*diag(Q)
  b=solve(Hes)%*%t(X2)%*%diag(as.vector(w))%*%ksi ####bias
  beta.hat=beta-b
  predict[j]=pnorm(as.matrix(cbind(1,Data[j,1:5]))%*%beta.hat)
}


pred=prediction(predict, y_sample)
perf1_sp=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc1_sp=as.numeric(auc.tmp@y.values)
auc1_sp

prob=glm(y_sample~., family=binomial(link="probit"), data=Data)
beta=coef(prob)
X2=as.matrix(cbind(1,Data[,1:5]))
eta=X2%*%beta
w=dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta))*(w0*pnorm(eta)+w1*(1-pnorm(eta)))
Hes=matrix(0,nrow=6, ncol=6)
for(i in 1:78){
  Hes=Hes+w[i]*X2[i,]%*%t(X2[i,])
}
Q=(X2)%*%solve(Hes)%*%t(X2)
ksi=eta/2*diag(Q)
b=solve(Hes)%*%t(X2)%*%diag(as.vector(w))%*%ksi ####bias
beta.hat=beta-b

###########################################################################
###########################Statistical Model###############################
D_pt=matrix(0,100,17)
for(i in 1:100){
  for(j in 1:17){
    
    D_pt[i,j]=quantile(D[i,],0.05*j+0.05)
  }
}
colnames(D_pt)=c("D90","D85","D80","D75","D70","D65","D60","D55","D50",
                 "D45","D40","D35","D30","D25","D20","D15","D10")

Xc=cbind(X_sample, D_pt)
pred_stat=rep(0,100)
for(i in 1:100){
  lasso1=glmnet(Xc[-i,], y[-i], family="binomial")
  AICC=(1-lasso1$dev.ratio)*lasso1$nulldev+2*colSums(lasso1$beta!= 0)+
    2*colSums(lasso1$beta!= 0)*(colSums(lasso1$beta!= 0)+1)/(nrow(X)-1-colSums(lasso1$beta!= 0))
  j=which.min(AICC)
  pred_stat[i]=predict(lasso1, Xc, s=lasso1$lambda[j], type="response")[i]
}
pred=prediction(pred_stat, y_sample)
perf3_sp=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc3_sp=as.numeric(auc.tmp@y.values)

lasso0=glmnet(Xc, y_sample, family="binomial")
AICC=(1-lasso0$dev.ratio)*lasso0$nulldev+2*colSums(lasso0$beta!= 0)+
  2*colSums(lasso0$beta!= 0)*(colSums(lasso0$beta!= 0)+1)/(nrow(Xc)-1-colSums(lasso0$beta!= 0))
j=which.min(AICC)
lasso0$beta[,32]

###########################################################################################
###########################Simulation containing error bar#################################
beta_ps=c(rnorm(1,0.4,0.05), rnorm(1,-0.4,0.05),rnorm(1,0.4,0.05), rnorm(1,-0.4,0.05),rep(0,10))
alpha=0.25;
n0=0.20;
N_test=1000; ####number of patients
genDVH=function(N){
  x=cov_dose
  DVH=abs(matrix(rmvnorm(N, mean=mean_dose, sigma=x, method="svd"), nrow=N, ncol=101))
  return(DVH)
}
DVH_test=genDVH(N_test)

sig=c(1, rho^seq(1:(length(beta_ps)-1)))
sig=toeplitz(sig)
X_test=rmvnorm(N_test, mean=rep(0,14),sigma=sig)
D_test=NULL
for(i in 1:100){
  D_test=cbind(D_test,(DVH_test[,i]+DVH_test[,i+1])/2)
}
Geud_test=rep(0,N_test)
for(i in 1:N_test){
  Geud_test[i]=sum(0.01*D_test[i,]^(1/n0))^n0
}
mean(pnorm(X_test%*%beta_ps+alpha*Geud_test-12.60))
y_test=rep(0,N_test)
for(i in 1:N_test){
  y_test[i]=rbinom(1,1,pnorm(X_test[i,]%*%beta_ps+alpha*Geud_test[i]-12.60))  ####ratio approximately 1:4####
}
sum(y_test)
prob=glm(y_test~., family=binomial(link="probit"), data=data.frame(X_test[,1:4], y_test, Geud_test))
summary(prob)



lambda=exp(seq(log(0.1),log(50),length.out=20))
trainAUC=AUC=rep(0,50)
BETA=matrix(NA, ncol=50,nrow=14)
Alpha=rep(0,50)
n_sp=rep(0,50)
trainAUC_stat=AUC_stat=rep(0,50)
BETA_stat=matrix(NA, ncol=50, nrow=14)
Alpha_stat=matrix(NA, ncol=50,nrow=17)
lambda=exp(seq(log(0.1),log(50),length.out=20))
AUC=rep(0,50)
BETA=matrix(NA, ncol=50,nrow=14)
Alpha=rep(0,50)
n_sp=rep(0,50)
AUC_stat=rep(0,50)
BETA_stat=matrix(NA, ncol=50, nrow=14)
Alpha_stat=matrix(NA, ncol=50,nrow=17)
AUC_LKB=rep(0,50)
n_LKB=rep(0,50)
for(l in 1:50){
  N=1000; ####number of patients
  DVH=genDVH(N)
  
  rho=0.2;
  sig=c(1, rho^seq(1:(length(beta_ps)-1)))
  sig=toeplitz(sig)
  X_ps=rmvnorm(N, mean=rep(0,14),sigma=sig)
  D=NULL
  for(i in 1:100){
    D=cbind(D,(DVH[,i]+DVH[,i+1])/2)
  }
  Geud=rep(0,100)
  for(i in 1:N){
    Geud[i]=sum(0.01*D[i,]^(1/n0))^n0
  }
  mean(pnorm(X_ps%*%beta_ps+alpha*Geud-12.60))
  y_ps=rep(0,N)
  for(i in 1:N){
    y_ps[i]=rbinom(1,1,pnorm(X_ps[i,]%*%beta_ps+alpha*Geud[i]-12.60))  ####ratio approximately 1:4####
  }
  sum(y_ps)
  prob=glm(y_ps~., family=binomial(link="probit"), data=data.frame(X_ps[,1:5], y_ps, Geud))
  summary(prob)
  
  int_1=which(y_ps==1)
  int_0=which(y_ps==0)
  int_sample=c(sample(int_1)[1:80],sample(int_0)[1:120])
  y_sample=y_ps[int_sample]
  X_sample=X_ps[int_sample,]
  D_sample=D[int_sample,]
  prob=glm(y_sample~., family=binomial(link="probit"), data=data.frame(X_sample[,1:5], y_sample, Geud[int_sample]))
  summary(prob)
  
  nc=seq(0, 0.5, length.out=20)
  geud_sp=rep(0,200)
  geud_test=rep(0,N_test)
  geud_sample=rep(0,200)
  BIC=rep(0,length(nc))
  AIC=rep(0,length(nc))
  AICC=rep(0,length(nc))
  Beta=matrix(0,ncol=50, nrow=16)
  
  
  for(k in 1:length(nc)){
    
    n=nc[k]
    for(i in 1:200){
      geud_sp[i]=(sum(D_sample[i,]^(1/n))*0.01)^(n)  ####calculate gEUD
    }
    
    Xa=cbind(X_sample,geud_sp);
    pre=mglmnet_probit(Xa, y_sample,c(rep(1,14),0), lambda,tau)
    j=which.min(pre$BIC)[1]
    BIC[k]=pre$BIC[j]
    Beta[,k]=pre$beta[,j]
  }
  
  k=which.min(BIC)
  beta=Beta[,k]
  BETA[,l]=beta[-c(1,16)]
  n=nc[k]
  n_sp[l]=n
  for(i in 1:200){
    geud_sample[i]=(sum(D_sample[i,]^(1/n))*0.01)^(n)  ####calculate gEUD
  }
  Xa=cbind(X_sample,geud_sample);
  ind=which(beta!=0)[-1]
  X1=Xa[,ind-1]
  Data=data.frame(X1,y_sample)
  predict=rep(0,N_test)
  prob1=glm(y_sample~., family=binomial(link="probit"), data=Data)
  beta=coef(prob1)
  X2=as.matrix(cbind(1,Data[,-length(beta)]))
  eta=X2%*%beta
  c=(1-tau)*mean(y_sample)/tau/(1-mean(y_sample))
  w1=tau/mean(y_sample)
  w0=(1-tau)/(1-mean(y_sample))
  w=dnorm(eta)^2/pnorm(eta)/(1-pnorm(eta))*(w0*pnorm(eta)+w1*(1-pnorm(eta)))
  for(i in 1:200){
    if(pnorm(eta)[i]<=1e-10) w[i]=w0*pnorm(eta)[i]
    if(pnorm(eta)[i]>=1-1e-10) w[i]=w1*(1-pnorm(eta)[i])
  }
  Hes=matrix(0,nrow=length(beta), ncol=length(beta))
  for(i in 1:200){
    Hes=Hes+w[i]*X2[i,]%*%t(X2[i,])
  }
  Q=(X2)%*%solve(Hes)%*%t(X2)
  ksi=eta/2*diag(Q)
  b=solve(Hes)%*%t(X2)%*%diag(as.vector(w))%*%ksi ####bias
  beta.hat=beta-b
  predict_train=pnorm(as.matrix(cbind(1,X1)%*%beta.hat))
  pred=prediction(predict_train, y_sample)
  perf1_sp=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc1_sp=as.numeric(auc.tmp@y.values)
  trainAUC[l]=auc1_sp
  
  for(i in 1:N_test){
    geud_test[i]=(sum(D_test[i,]^(1/n))*0.01)^(n)  ####calculate gEUD
  }
  Xa=cbind(X_test,geud_test);
  X1=Xa[,ind-1]
  predict=pnorm(as.matrix(cbind(1,X1)%*%beta.hat))
  Alpha[l]=beta.hat[length(beta.hat)]
  
  pred=prediction(predict, y_test)
  perf1_sp=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc1_sp=as.numeric(auc.tmp@y.values)
  AUC[l]=auc1_sp
  
  #######################################################################################
  ###############Simulation with error bar, statistical model############################
  D_pt=matrix(0,200,17)
  for(i in 1:200){
    for(j in 1:17){
      
      D_pt[i,j]=quantile(D[i,],0.05*j+0.05)
    }
  }
  colnames(D_pt)=c("D90","D85","D80","D75","D70","D65","D60","D55","D50",
                   "D45","D40","D35","D30","D25","D20","D15","D10")
  D_tt=matrix(0,N_test,17)
  for(i in 1:N_test){
    for(j in 1:17){
      
      D_tt[i,j]=quantile(D_test[i,],0.05*j+0.05)
    }
  }
  colnames(D_tt)=c("D90","D85","D80","D75","D70","D65","D60","D55","D50",
                   "D45","D40","D35","D30","D25","D20","D15","D10")
  
  Xc=cbind(X_sample, D_pt)
  Xc_test=cbind(X_test, D_tt)
  pred_stat=rep(0,N_test)
  lasso1=glmnet(Xc, y_sample,family="binomial")
  AICC=(1-lasso1$dev.ratio)*lasso1$nulldev+2*colSums(lasso1$beta!= 0)+
    2*colSums(lasso1$beta!= 0)*(colSums(lasso1$beta!= 0)+1)/(nrow(X)-1-colSums(lasso1$beta!= 0))
  BIC=(1-lasso1$dev.ratio)*lasso1$nulldev+log(nrow(Xc))*colSums(lasso1$beta!= 0)
  j=which.min(BIC)
  beta=lasso1$beta[,j]
  BETA_stat[,l]=beta[1:14]
  Alpha_stat[,l]=beta[15:length(beta)]
  pred_train=predict(lasso1,Xc, s=lasso1$lambda[j],type="response")
  pred=prediction(pred_train, y_sample)
  perf3_sp=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc3_sp=as.numeric(auc.tmp@y.values)
  trainAUC_stat[l]=auc3_sp
  
  pred_stat=predict(lasso1, Xc_test, s=lasso1$lambda[j], type="response")
  pred=prediction(pred_stat, y_test)
  perf3_sp=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc3_sp=as.numeric(auc.tmp@y.values)
  AUC_stat[l]=auc3_sp
  ######################LKB model with error bar#######################
  pred_LKB=rep(0,N_test)
  dntcpl <- function(param){
    
    TD_50=param[1];
    m0=param[2];
    n=param[3];
    lntcp=0;
    if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1){
      
      for(i in 1:200){
        
        d=D_sample[i,]
        eud=(sum(d^(1/n))/length(d))^(n)
        t=(eud-TD_50)/(m0*TD_50)
        lntcp=lntcp+y_sample[i]*log(pnorm(t))+(1-y_sample[i])*log(1-pnorm(t))
      }
    }else {lntcp=-300}
    return(lntcp)
  }
  dntcpl(c(62,0.2,0.1))
  model=maxLik(dntcpl,start=c(60,0.1,0.1),method="nm")
  summary(model)
  coef=coef(model)
  TD_50=coef[1]
  m0=coef[2]
  n=coef[3]
  n_LKB[l]=n
  for(i in 1:N_test){
    d.pre=D_test[i,]
    eud=(sum(d.pre^(1/n))/length(d.pre))^(n)
    t=(eud-TD_50)/(m0*TD_50)
    pred_LKB[i]=pnorm(t)
  }
  
  pred=prediction(pred_LKB, y_test)
  perf2=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc2=as.numeric(auc.tmp@y.values)
  AUC_LKB[l]=auc2
}


TP=rep(0,50)  ####accuracy
FP=rep(0,50)  ####sensitivity
FN=rep(0,50) ####specificity
TN=rep(0,50) ####precision

TP_stat=rep(0,50)  ####accuracy
FP_stat=rep(0,50)  ####sensitivity
FN_stat=rep(0,50) ####specificity
TN_stat=rep(0,50) ####precision

MA1=matrix(0,ncol=50, nrow=14)
MA2=matrix(0, ncol=50, nrow=14)
beta_ind=rep(0,14)
MA1[BETA!=0]=1
MA2[BETA_stat!=0]=1
beta_ind[beta_ps!=0]=1
for(i in 1:50){
  tem=MA1[,i]
  TP[i]=length(which(tem==1&beta_ind==1))
  FP[i]=length(which(tem==1&beta_ind==0))
  FN[i]=length(which(tem==0&beta_ind==1))
  TN[i]=length(which(tem==0&beta_ind==0))
}

for(i in 1:50){
  
  TP_stat[i]=length(which(MA2[,i]==1&beta_ind==1))
  FP_stat[i]=length(which(MA2[,i]==1&beta_ind==0))
  FN_stat[i]=length(which(MA2[,i]==0&beta_ind==1))
  TN_stat[i]=length(which(MA2[,i]==0&beta_ind==0))
}

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)

TPR_stat=TP_stat/(TP_stat+FN_stat)
FPR_stat=FP_stat/(FP_stat+TN_stat)

quantile(n_sp, c(0.025,0.975))
quantile(n_LKB, c(0.025,0.975))
quantile(TPR, c(0.025,0.975))
quantile(FPR, c(0.025,0.975))
quantile(TPR_stat, c(0.025,0.975))
quantile(FPR_stat, c(0.025,0.975))


mean(TPR); sd(TPR)
mean(1-FPR); sd(1-FPR)
mean(TPR_stat); sd(TPR_stat)
mean(1-FPR_stat); sd(1-FPR_stat)
mean(AUC); sd(AUC)
mean(AUC_stat); sd(AUC_stat);
mean(AUC_LKB); sd(AUC_LKB)

wilcox.test(AUC_stat, AUC)
wilcox.test(TPR, TPR_stat)
wilcox.test(FPR, FPR_stat)





#AUCmatrix=matrix(0,nrow=50,ncol=4)
#AUC_statmatrix=matrix(0, nrow=50,ncol=4)
#AUC_LKBmatrix=matrix(0,nrow=50,ncol=4)
#TPRmatrix=matrix(0,nrow=50,ncol=4)
#FPRmatrix=matrix(0,nrow=50,ncol=4)
#TPR_statmatrix=matrix(0,nrow=50,ncol=4)
#FPR_statmatrix=matrix(0,nrow=50,ncol=4)

AUCmatrix[,4]=AUC
AUC_statmatrix[,4]=AUC_stat
AUC_LKBmatrix[,4]=AUC_LKB
TPRmatrix[,4]=TPR
TPR_statmatrix[,4]=TPR_stat
FPRmatrix[,4]=FPR
FPR_statmatrix[,4]=FPR_stat

