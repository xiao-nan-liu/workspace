install.packages("glmnet")
install.packages("pROC")
install.packages("maxLik")
install.packages("ROCR")
library(maxLik)
library(ROCR)
library(glmnet)
library("pROC")
###Initialization###
setwd("C:/Users/xliu203/Desktop/new Vector")
rec=read.csv("C:/Users/xliu203/Desktop/rectum volume.csv")
fileName <- dir()
datalist <- vector("list",79)
recvol=rec$X*rec$Y*rec$Z
recvol=recvol[-c(8,13,34,40,43,64,74)]
for(i in 1:79){
datalist[[i]]=read.table(fileName[i],header=T)
}

for(i in 1:79){
recvol[i]=recvol[i]*length(datalist[[i]][[1]])/1000
}

patient=read.csv("C:/Users/xliu203/Desktop/new patient.csv")
patient1=read.csv("C:/Users/xliu203/Desktop/300 patients.csv")
Acute=patient$Acute.GI.Toxicity
Late=patient$Max.Chronic.GI.Toxicity  ##should use final
Acute1=patient1$Acute.GI.Toxicity
Late1=patient1$Max.Chronic.GI.Toxicity
R=R2=rep(0,79)
R[Acute<=1]=0
R[Acute>1]=1
R[Late<=1]=0
R[Late>1]=1
####Effective Dose Reduction Scheme 2Gy fraction####
dntcpl <- function(param){
	
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	lntcp=0;
	if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1){

	for(i in 1:79){

		d=datalist[[i]][[1]]
		eud=(sum(d^(1/n))/length(d))^(n)
		t=(eud-TD_50)/(m0*TD_50)
		lntcp=lntcp+R[i]*log(pnorm(t))+(1-R[i])*log(1-pnorm(t))
	}
	}else {lntcp=-100}
	return(lntcp)
}

####Effective Volume Reduction Scheme 2Gy fraction####
vntcpl <- function(param){

	TD_50=param[1];
	m0=param[2];
	n=param[3];
	lntcp=0;
	for(i in 1:length){
		
		D=datalist[[i]][[1]]/100*0.96
		d=D/max(D)
		V=sum(d^(1/n)/length(d))
		t=(max(D)-TD_50/(V^n))/(m0*TD_50/(V^n))
		lntcp=lntcp+R[i]*log(pnorm(t))+(1-R[i])*log(1-pnorm(t))
	}
	return(lntcp)
}
####Estimation on Parameters#####
vntcpl(c(50,1,1))
dntcpl(c(65,0.2,0.3))
model=maxLik(dntcpl,start=c(65,0.2,0.3),method="nm")
summary(model)
coef=coef(model)
n=coef[3]
geud=rep(0,79)
for(i in 1:79){
d=datalist[[i]][[1]]/100*0.96
geud[i]=(sum(d^(1/n))/length(d))^(n)  ##calculate gEUD
}

####Error####
entcp <- function(param){
	
	lntcp0=rep(0,79)
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	for(i in 1:79){

		d=datalist[[i]][[1]]
		eud=(sum(d^(1/n))/length(d))^(n)
		t=(eud-TD_50)/(m0*TD_50)
		lntcp0[i]=pnorm(t)
	}

	pred=prediction(lntcp0,R)
	perf1=performance(pred,"tpr","fpr")
	auc.tmp=performance(pred,"auc")
	auc=as.numeric(auc.tmp@y.values)
	return(c(auc,lntcp0))
}
auc1=entcp(coef(model))       ###model with no cv
auc2=entcp(c(76.9,0.13,0.09)) ###Quantec 2Gy
auc3=entcp(c(81,0.14,0.068)) ###NTCP modeling of late rectal
auc4=entcp(c(68.5,0.15,0.13)) ###Parameters for LKB

####Interval estimation####
theta.start <- coef
md<-function(param) -dntcpl(param)
out <- nlm(md, theta.start,hessian=TRUE)
theta.hat <- out$estimate
fish=out$hessian
theta.hat
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
inv.fish <- solve(fish)
theta.hat[1] + c(-1, 1) * crit * sqrt(inv.fish[1, 1])/sqrt(79)
theta.hat[2] + c(-1, 1) * crit * sqrt(inv.fish[2, 2])/sqrt(79)
theta.hat[3] + c(-1, 1) * crit * sqrt(inv.fish[3, 3])/sqrt(79)

###sampling using LKB###
length1=NULL;length2=NULL
for (i in 1:79){
if(R[i]==1) length1=c(length1,i)
else length2=c(length2,i)
}
length1_sample=sample(length1,4*length(length1),replace=TRUE)
id=c(length1_sample,length2)
sntcpl <- function(param){
	
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	lntcp=0;
	if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1){
	for(i in 1:length(id)){

		d=datalist[[id[i]]][[1]]/100*0.96
		eud=(sum(d^(1/n))/length(d))^(n)
		t=(eud-TD_50)/(m0*TD_50)
		lntcp=lntcp+R[id[i]]*log(pnorm(t))+(1-R[id[i]])*log(1-pnorm(t))
	}
	}
	else{ 
		lntcp=-1000
      }
	return(lntcp)
}
sntcpl(c(56,0.3,0.2))
smodel=maxLik(sntcpl,start=c(65,0.3,0.4),method="nm")
summary(smodel)
scoef=coef(smodel)

####Interval estimation####
theta.start <- scoef
md<-function(param) -sntcpl(param)
out <- nlm(md, theta.start,hessian=TRUE)
theta.hat <- out$estimate
fish=out$hessian
theta.hat
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
inv.fish <- solve(fish)
theta.hat[1] + c(-1, 1) * crit * sqrt(inv.fish[1, 1])/sqrt(79)
theta.hat[2] + c(-1, 1) * crit * sqrt(inv.fish[2, 2])/sqrt(79)
theta.hat[3] + c(-1, 1) * crit * sqrt(inv.fish[3, 3])/sqrt(79)

###Modified LKB Model###
T=patient$age.at.IMRT
T=patient$Diabetes.y.n.
T=patient$psa.prior.to.IMRT
T=patient$Statins.2
antcpl <- function(param){
	
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	alpha1=param[4];
	lntcp=0;
	if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1){
		for(i in 1:79){
				d=datalist[[i]][[1]]
				eud=(sum(d^(1/n)/length(d)))^(n)
				TD50=TD_50+alpha1*T[i]
				t=(eud-TD50)/(m0*TD50)
				lntcp=lntcp+R[i]*log(pnorm(t))+
				(1-R[i])*log(1-pnorm(t))
			}
		} else{
			lntcp=-100
	}
	return(lntcp)
}
antcpl(c(60,0.2,0.3,0))
model=maxLik(antcpl,start=c(c(60,0.1,0.1,0)),method="nm")
summary(model)
acoef=coef(model)
param=acoef

###Modified LKB error####
entcp2 <- function(param){
	
	lntcp=NULL
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	alpha1=param[4];
	for(i in 1:79){
		
		d=datalist[[i]][[1]]
		eud=(sum(d^(1/n)/length(d)))^(n)
		TD50=TD_50+alpha1*T[i]
		t=(eud-TD50)/(m0*TD50)
		lntcp=c(lntcp,pnorm(t))
	}

	pred=prediction(lntcp,R)
	perf2=performance(pred,"tpr","fpr")
	auc.tmp=performance(pred,"auc")
	auc=as.numeric(auc.tmp@y.values)
	return(c(auc,lntcp))
}
auc6=entcp2(acoef)
###Modified LKB confidence interval###
theta.start <- acoef
md<-function(param) -antcpl(param)
out <- nlm(md, theta.start,hessian=TRUE)
theta.hat <- out$estimate
fish=out$hessian
theta.hat
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
inv.fish <- solve(fish)
theta.hat[1] + c(-1, 1) * crit * sqrt(inv.fish[1, 1])/sqrt(79)
theta.hat[2] + c(-1, 1) * crit * sqrt(inv.fish[2, 2])/sqrt(79)
theta.hat[3] + c(-1, 1) * crit * sqrt(inv.fish[3, 3])/sqrt(79)
theta.hat[4] + c(-1, 1) * crit * sqrt(inv.fish[4, 4])/sqrt(79)

####checking the grade 1 toxicity####

pntcpl <- function(param){
	
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	p=param[4]
	lntcp=0;
	if(TD_50>50&&TD_50<100&&m0>0&&m0<1&&n>0&&n<1&&p>-0.0001&&p<1){
		for(i in 1:length){

			d=datalist[[i]][[1]]/100*0.96
			eud=(sum(d^(1/n)/length(d)))^(n)
			t=(eud-TD_50)/(m0*TD_50)
                  if(Acute[i]==2){
                  lntcp=lntcp+log(pnorm(t))
				}
			if(Acute[i]==0){
                  lntcp=lntcp+log(1-pnorm(t))
				}
			if(Acute[i]==1){
                  lntcp=lntcp+p*log(pnorm(t))+(1-p)*log(1-pnorm(t))
				}

			}
		} else{
			lntcp=-100
	}
	lntcp
}
pntcpl(c(coef,0))

model=maxLik(pntcpl,start=c(coef,0),method="nm")
summary(model)
acoef=coef(model)
entcp2 <- function(param){
	
	lntcp=NULL
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	for(i in 1:length(Acute[Acute!=1])){
		
		d=datalist[[i]][[1]]/100*0.96
		eud=(sum(d^(1/n)/length(d)))^(n)
		t=(eud-TD_50)/(m0*TD_50)
		lntcp=c(lntcp,pnorm(t))
	}
      R=Acute[Acute!=1]
      R[R>1]=1
	pred=prediction(lntcp,Acute[Acute!=1])
      perf=performance(pred,"tpr","fpr")
	auc.tmp=performance(pred,"auc")
	auc=as.numeric(auc.tmp@y.values)
	return(auc)
}
auc6=entcp2(acoef)

###confidence interval for grade 1 ###
theta.start <- acoef
md<-function(param) -antcpl(param)
out <- nlm(md, theta.start,hessian=TRUE)
theta.hat <- out$estimate
fish=out$hessian
theta.hat
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
inv.fish <- solve(fish)
theta.hat[1] + c(-1, 1) * crit * sqrt(inv.fish[1, 1])/sqrt(79)
theta.hat[2] + c(-1, 1) * crit * sqrt(inv.fish[2, 2])/sqrt(79)
theta.hat[3] + c(-1, 1) * crit * sqrt(inv.fish[3, 3])/sqrt(79)
theta.hat[4] + c(-1, 1) * crit * sqrt(inv.fish[4, 4])/sqrt(79)

####binomial model####

bntcpl <- function(param){
	
	TD1_50=param[1]; TD2_50=param[3];
	m1=param[2]; m2=param[4]
	n=param[5]; mu=param[6]
	lntcp=0;		
	for(i in 1:length){

		d=datalist[[i]][[1]]/100*0.96
		eud=(sum(d^(1/n)/length(d)))^(n)
		t1=(eud-TD1_50)/(m1*TD1_50)
		t2=(eud-TD2_50)/(m2*TD2_50)
            if(Acute[i]!=2&&Late[i]!=2){
              lntcp=lntcp+log(1-pnorm(t1)-pnorm(t2)+mu)
			}
		if(Acute[i]==2&&Late[i]!=2){
              lntcp=lntcp+log(pnorm(t1)-mu)
			}
		if(Acute[i]!=2&&Late[i]==2){
               lntcp=lntcp+log(pnorm(t2)-mu)
			}
		if(Acute[i]==2&&Late[i]==2){
               lntcp=lntcp+log(mu)
			}
		}
	lntcp
}
bntcpl(c(60,0.06,70,0.11,0.08,0.02))
model=maxLik(bntcpl,start=c(61,0.06,72,0.11,0.08,0.02),method="nm")
summary(model)
bcoef=coef(model)

####prior correction####
dprior <- function(param){
	
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	lntcp=0;
	for(i in 1:length){

		d=datalist[[i]][[1]]/100*0.96
		eud=(sum(d^(1/n))/length(d))^(n)
		t=(eud-TD_50)/(m0*TD_50)
		p=1/(1+(1/pnorm(t)-1)*64/79/15/79*3/79/(76/79))
		lntcp=lntcp+R[i]*log(p)+(1-R[i])*log(1-p)
	}
	return(lntcp)
}
dprior(c(65,0.2,0.1))
model=maxLik(dprior,start=c(65,0.2,0.1),method="nm")
summary(model)
coef=coef(model)
entcp <- function(param){
	
	lntcp=rep(0,79)
	TD_50=param[1];
	m0=param[2];
	n=param[3];
	for(i in 1:length){

		d=datalist[[i]][[1]]/100*0.96
		eud=(sum(d^(1/n))/length(d))^(n)
		t=(eud-TD_50)/(m0*TD_50)
		lntcp[i]=p=1/(1+(1/pnorm(t)-1)*0.99/0.01*3/79/(76/79))
	}

	pred0=prediction(lntcp,R)
	auc.tmp=performance(pred0,"auc")
	auc=as.numeric(auc.tmp@y.values)
	return(auc)
}
auc1=entcp(coef(model))       ###model with no cv
auc2=entcp(c(76.9,0.13,0.09)) ###Quantec 2Gy
auc3=entcp(c(81,0.14,0.068)) ###NTCP modeling of late rectal
auc4=entcp(c(68.5,0.15,0.13)) ###Parameters for LKB
##############plot DVH#############
d=datalist[[1]][[1]]
h <- hist(d, plot=FALSE, breaks=100)
h$counts     <- 20170-cumsum(h$counts)
h$density    <- 1-cumsum(h$density)
plot(h, freq=F, xlab="Dose(Gy)",ylab="Fraction Volume", main="Culmulative DVH for patient 1", col="white", border="black", xlim=c(0,80))
plot(density(d), xlab="Dose(Gy)",ylab="Fraction Volume", main="Differential DVH for patient 1", type="h" )
hist(d, breaks=100, xlab="Dose(Gy)",ylab="Fraction Volume", main="Differential DVH for patient 1", freq=F, xlim=c(0,80))
plot(sort(d), 1-(1:length(d))/length(d), type="h", xlab="Dose(Gy)",ylab="Fraction Volume", main="Culmulative DVH for patient 1")

####new model using logistic regression###
x03=rep(0,79)
for(i in 1:79){
x03[i]=quantile(datalist[[i]][[1]],1-0.001)
}
x02=rep(0,79)
for(i in 1:79){
x02[i]=quantile(datalist[[i]][[1]],1-0.07)
}
x01=rep(0,79)
for(i in 1:79){
x01[i]=quantile(datalist[[i]][[1]],1-0.21)
}
x00=rep(0,79)
for(i in 1:79){
x00[i]=quantile(datalist[[i]][[1]],1-0.29)
}
x_D40=patient$rectal.D40..gy.
sum((x_D40-x00)^2)
x_D30=patient$rectal.d30..gy.
sum((x_D30-x01)^2)
x_D10=patient$rectal.d10..gy.
sum((x_D10-x02)^2)
x_D1.8=patient$rectal.d1.8.gy.
sum((x_D1.8-x03)^2)
x_Dmean=patient$Mean.Rectal.Dose
x04=rep(0,79)
for(i in 1:79){
x04[i]=mean(datalist[[i]][[1]])
}
x0=rep(0,79)
for(i in 1:79){
d=datalist[[i]][[1]]
x0[i]=length(d[d>=50])/length(d)
}
x0=rep(0,79)
for(i in 1:79){
d=datalist[[i]][[1]]
x0[i]=quantile(datalist[[i]][[1]],1-0.25)
}

z=matrix(0,79,13)
for(i in 1:79){
for(j in 1:13){
d=datalist[[i]][[1]]
z[i,j]=length(d[d>=5*j+5])/length(d)
}
}
colnames(z)=c("V10","V15","V20","V25","V30","V35","V40","V45","V50","V55","V60","V65","V70")

z0=matrix(0,79,17)
for(i in 1:79){
for(j in 1:17){
d=datalist[[i]][[1]]
z0[i,j]=quantile(datalist[[i]][[1]],0.05*j+0.05)
}
}
colnames(z0)=c("D90","D85","D80","D75","D70","D65","D60","D55","D50","D45","D40","D35","D30","D25","D20","D15","D10")

##############multiple indices#############
data_1=data.frame(R,z)
logit=glm(R~1,data=data_1,family="binomial")
logit1=glm(R~., data=data_1, family="binomial")
step(logit1, k=2, direction="backward")
step(logit, k=2, direction="forward", scope=formula(logit1))

leave_one_out_cv=function(data_1){
ntcp_cv=rep(0,79)
for(i in 1:79){
data_1f=data_1[-i,]
logit2=glm(R~., data=data_1f, family="binomial")
ntcp_cv[i]=predict(logit2,newdata=data_1[i,],type="response")
}
pred=prediction(ntcp_cv, R)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
}
AUC_max=NULL
lAUC=rep(0,13)
for(j in 1:13){
  data_1m=data_1[,c(1,j+1)]
  lAUC[j]=leave_one_out_cv(data_1m)
}
AUC_max=c(AUC_max, lAUC[which.max(lAUC)])
lAUC2=rep(0,13)
for(j in 1:13){
  if(j!=which.max(lAUC)){
    data_1m=data_1[,c(1,which.max(lAUC)+1,j+1)]
    lAUC2[j]=leave_one_out_cv(data_1m)
  }
}
AUC_max=c(AUC_max, lAUC2[which.max(lAUC2)])
lAUC3=rep(0,13)
for(j in 1:13){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)){
    data_1m=data_1[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,j+1)]
    lAUC3[j]=leave_one_out_cv(data_1m)
  }
}
which.max(lAUC3)
AUC_max=c(AUC_max, lAUC3[which.max(lAUC3)])
lAUC4=rep(0,13)
for(j in 1:13){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)&&j!=which.max(lAUC3)){
    data_1m=data_1[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,which.max(lAUC3)+1,j+1)]
    lAUC4[j]=leave_one_out_cv(data_1m)
  }
}
which.max(lAUC4)
AUC_max=c(AUC_max, lAUC4[which.max(lAUC4)])
lAUC5=rep(0,13)
for(j in 1:13){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)&&j!=which.max(lAUC3)&&j!=which.max(lAUC4)){
    data_1m=data_1[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,which.max(lAUC3)+1,which.max(lAUC4)+1,j+1)]
    lAUC5[j]=leave_one_out_cv(data_1m)
  }
}
which.max(lAUC5)
AUC_max=c(AUC_max, lAUC5[which.max(lAUC5)])

lAUC6=rep(0,13)
for(j in 1:13){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)&&j!=which.max(lAUC3)&&j!=which.max(lAUC4)&&j!=which.max(lAUC5)){
    data_1m=data_1[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,which.max(lAUC3)+1,which.max(lAUC4)+1,which.max(lAUC5)+1,j+1)]
    lAUC6[j]=leave_one_out_cv(data_1m)
  }
}
which.max(lAUC6)
AUC_max=c(AUC_max, lAUC6[which.max(lAUC6)])

colnames(z)[c(which.max(lAUC),which.max(lAUC2),which.max(lAUC3),which.max(lAUC4),which.max(lAUC5),which.max(lAUC6))]
AUC_max

data_2=data.frame(R,z0)
logit=glm(R~1,data=data_2,family="binomial")
logit1=glm(R~., data=data_2, family="binomial")
step(logit1, k=2, direction="backward")
step(logit, k=2, direction="forward", scope=formula(logit1))
logit3=glm(R~D20+D25+D60, data=data_2, family="binomial")
ntcp=predict(logit1,type="response")
pred=prediction(ntcp, R)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp)

leave_one_out_cv=function(data_2m){
  ntcp_cv=rep(0,79)
  for(i in 1:79){
    data_2f=data_2m[-i,]
    data_2g=data_2m[i,]
    logit2=glm(R~., data=data_2f, family="binomial")
    ntcp_cv[i]=predict(logit2,newdata=data_2g,type="response")
  }
  pred=prediction(ntcp_cv, R)
  perf3=performance(pred,"tpr","fpr")
  auc.tmp=performance(pred,"auc")
  auc=as.numeric(auc.tmp@y.values)
  auc
}
ntcp_cv2=rep(0,79)
for(i in 1:79){
  data_2f=data_2[-i,]
  logit3=glm(R~D70 + D55 + D50 + D40 + D20 + D15 + D10, data=data_2f, family="binomial")
  ntcp_cv2[i]=predict(logit3,newdata=data_2[i,],type="response")
}
pred=prediction(ntcp_cv2, R)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp_cv2)

AUC_max2=NULL
lAUC=rep(0,17)
for(j in 1:17){
  data_2m=data_2[,c(1,j+1)]
  lAUC[j]=leave_one_out_cv(data_2m)
}
AUC_max2=c(AUC_max2, lAUC[which.max(lAUC)])
lAUC2=rep(0,17)
for(j in 1:17){
  if(j!=which.max(lAUC)){
    data_2m=data_2[,c(1,which.max(lAUC)+1,j+1)]
    lAUC2[j]=leave_one_out_cv(data_2m)
  }
}
AUC_max2=c(AUC_max2, lAUC2[which.max(lAUC2)])
lAUC3=rep(0,17)
for(j in 1:17){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)){
    data_2m=data_2[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,j+1)]
    lAUC3[j]=leave_one_out_cv(data_2m)
  }
}
which.max(lAUC3)
AUC_max2=c(AUC_max2, lAUC3[which.max(lAUC3)])
lAUC4=rep(0,17)
for(j in 1:17){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)&&j!=which.max(lAUC3)){
    data_2m=data_2[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,which.max(lAUC3)+1,j+1)]
    lAUC4[j]=leave_one_out_cv(data_2m)
  }
}
which.max(lAUC4)
AUC_max2=c(AUC_max2, lAUC4[which.max(lAUC4)])
lAUC5=rep(0,17)
for(j in 1:17){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)&&j!=which.max(lAUC3)&&j!=which.max(lAUC4)){
    data_2m=data_2[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,which.max(lAUC3)+1,which.max(lAUC4)+1,j+1)]
    lAUC5[j]=leave_one_out_cv(data_2m)
  }
}
which.max(lAUC5)
AUC_max2=c(AUC_max2, lAUC5[which.max(lAUC5)])

lAUC6=rep(0,17)
for(j in 1:17){
  if(j!=which.max(lAUC)&&j!=which.max(lAUC2)&&j!=which.max(lAUC3)&&j!=which.max(lAUC4)&&j!=which.max(lAUC5)){
    data_2m=data_2[,c(1,which.max(lAUC)+1,which.max(lAUC2)+1,which.max(lAUC3)+1,which.max(lAUC4)+1,which.max(lAUC5)+1,j+1)]
    lAUC6[j]=leave_one_out_cv(data_2m)
  }
}
which.max(lAUC6)
AUC_max2=c(AUC_max2, lAUC6[which.max(lAUC6)])

colnames(z0)[c(which.max(lAUC),which.max(lAUC2),which.max(lAUC3))],which.max(lAUC4),which.max(lAUC5),which.max(lAUC6))]
AUC_max2
########################################################
y=rep(0,79)
y[Acute<=1]=0
y[Acute>=2]=1
x1=patient$age.at.IMRT
x2=patient$Hormones.1.yes.2.no
x3=patient$T.stage.greater.than.2a...1.
x4=patient$gleason.score
x5=patient$psa.prior.to.IMRT
x6=patient$Diabetes.y.n.
x7=patient$Prostate.vol.cc.
x8=patient$Statins.2
x9=patient$Neoadjuvant.1.yes.2.no
x10=patient$Concurrent.1.yes.2.no
x11=patient$Adjuvant.1.yes.2.no
X_79=data.frame(x00,x01,x02,x03,x04,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y)
colnames(X_79)=c("D40","D30","D10","D1.8", "Dmean","age","Hormones","T.stage","gleason.score",
"psa","Diabetes","Prostate","Neoadjuvant","Concurrent","Adjuvant",
"Statins","toxicity")

logit=glm(R~x_D1.8, data=data.frame(R, x_D1.8),family="binomial")
ntcp=predict(logit,type="response")
pred=prediction(ntcp, R)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp)
y_lin1=logit$coef[2]*x03_new+logit$coef[1]
y_pre1=1/(1+exp(-y_lin1))
logit_new=glm(y~x03_new, family=binomial)
summary(logit_new)
ntcp=predict(logit_new,type="response")
pred=prediction(ntcp, y)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(y,ntcp)

logit=glm(R~x03, data=data.frame(R, x03),family="binomial")
ntcp=predict(logit,type="response")
pred=prediction(ntcp, R)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp)
y_lin2=logit$coef[2]*x02_new+logit$coef[1]
y_pre2=1/(1+exp(-y_lin2))
logit_new=glm(y~x02_new, family=binomial)
summary(logit_new)
ntcp=predict(logit_new,type="response")
pred=prediction(ntcp, y)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(y,ntcp)

logit=glm(R~x01, data=data.frame(R, x01),family="binomial")
ntcp=predict(logit,type="response")
pred=prediction(ntcp, R)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp)
y_lin3=logit$coef[2]*x01_new+logit$coef[1]
y_pre3=1/(1+exp(-y_lin3))
logit_new=glm(y~x01_new, family=binomial)
summary(logit_new)
ntcp=predict(logit_new,type="response")
pred=prediction(ntcp, y)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(y,ntcp)





X=data.frame(x0,x1,x5,x6,x8,R)
colnames(X)=c("dose", "age", "psa", "diabetes", "Statins","R")
logit=glm(R~dose+psa+Statins,data=X,family="binomial")
summary(logit)
ntcp=predict(logit,type="response")
pred=prediction(y_pre,y)
perf3=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp)
exp(confint.default(logit))
exp(coef(logit)[2])


logit=glm(R~z$D90,family="binomial")
summary(logit)
ntcp=predict(logit,type="response")
pred=prediction(ntcp,R)
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
ci.auc(R,ntcp)


ci.auc(R,ntcp)
#+age+psa+diabetes+Statins
logit=glm(acute~., data=X_79,family="binomial")
summary(logit)
logit0=glm(R~1, data=X,family="binomial")
back=step(logit,k=2)
fore=step(logit0,scope=list(lower=formula(logit0),upper=formula(logit))
          ,direction="forward")
both=step(logit0,scope=list(lower=formula(logit0),upper=formula(logit))
          ,direction="both",trace=0)
summary(glm(toxicity ~ age + D10 + D30, family = "binomial",data = X_79))

###
plot(perf1,col="1")
par(new=TRUE)
plot(perf2,col="2")
par(new=TRUE)
plot(perf3,col="3")
par(new=TRUE)
plot(perf4,col="4")
par(new=TRUE)
plot(perf5,col="5")
par(new=TRUE)
plot(perf6,col="6")
legend(0.8,0.75,c("Acute","Late"),
       lty=c(1,1),lwd=c(1.5,1.5), col=c(1,2))
legend(0.55,0.7,c("LKB","mLKB (psa)","Unrestricted multivariate\nlogistic regression"),
       lty=c(1,1,1),lwd=c(1.5,1.5,1.5), col=c(1,2,3))



####verify####
patient_new=read.csv("C:/Users/xliu203/Desktop/new data.csv")
y=rep(0,295)
x00_new=patient_new$rectal.D40..gy.*0.96
x01_new=patient_new$rectal.d30..gy.*0.96
x02_new=patient_new$rectal.d10..gy.*0.96
x03_new=patient_new$rectal.d1.8.gy.*0.96
x04_new=patient_new$Mean.Rectal.Dose*0.96
x1=patient_new$age.at.IMRT
x2=patient_new$Hormones.1.yes.2.no
x3=patient_new$T.stage.greater.than.2a...1.
x4=patient_new$gleason.score
x5=patient_new$psa.prior.to.IMRT
x6=patient_new$Diabetes.y.n.
x7=patient_new$Prostate.vol.cc.
x8=patient_new$Neoadjuvant.1.yes.2.no
x9=patient_new$Concurrent.1.yes.2.no
x10=patient_new$Adjuvant.1.yes.2.no
x11=patient_new$Statins.2
Acute1=patient_new$Acute.GI.Toxicity
Late1=patient_new$Max.Chronic.GI.Toxicity
y[Acute1<=1]=0
y[Acute1>=2]=1
y[Late1<=1]=0
y[Late1>=2]=1
X_295=data.frame(x00,x01,x02,x03,x04,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y)
colnames(X_295)=c("D40","D30","D10","D1.8", "Dmean","age","Hormones","T.stage","gleason.score",
"psa","Diabetes","Prostate","Neoadjuvant","Concurrent","Adjuvant",
"Statins","toxicity")
X_295=na.omit(X_295)

logit=glm(toxicity~., data=X_295,family="binomial")
summary(logit)
logit0=glm(toxicity~1, data=X_295,family="binomial")
back=step(logit,k=2)
fore=step(logit0,scope=list(lower=formula(logit0),upper=formula(logit))
          ,direction="forward")
both=step(logit0,scope=list(lower=formula(logit0),upper=formula(logit))
          ,direction="both",trace=0)
log=glm(toxicity ~ age + D10 + D30, family = "binomial",data = X_295)
log=glm(toxicity ~ Concurrent, family = "binomial",data = X_295)
ntcp=predict(log,type="response")
pred=prediction(ntcp,X_295$toxicity)
perf5=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc

predict=NULL
for(i in 1:79){
predict=c(predict,predict(log,newdata=X_79[i,],type="response"))
}
pred=prediction(predict,X_79$toxicity)
perf5=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc

####cross validation####
X=cbind(z[,9], x5, x8,R)
X=as.data.frame(X)
summary(glm(R~., data=X, family="binomial"))
Q=rep(0,79)
for(i in 1:79){
gmodel=glm(R~., data=X[-i,], family="binomial")
Q[i]=predict(gmodel, newdata=X[i,],type="response")
}
pred=prediction(Q, R)
perf5=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc


X=cbind(z0[,14], x5, x8,R)
X=as.data.frame(X)
summary(glm(R~., data=X, family="binomial"))
Q=rep(0,79)
for(i in 1:79){
gmodel=glm(R~., data=X[-i,], family="binomial")
Q[i]=predict(gmodel, newdata=X[i,],type="response")
}
pred=prediction(Q, R)
perf5=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc

#############################################################################
#######################power analysis########################################
t(cor(X[,-c(2,4,6)], R))%*%solve(cor(X[,-c(2,4,6)]))%*%cor(X[,-c(2,4,6)], R)
d=function(theta){
(1+(1+theta^2)*exp(1.25*theta^2))/(1+exp(-theta^2/4))
}
delta=d(1.7)
N=function(theta){
(qnorm(0.1)+exp(-theta^2)/4*qnorm(0.3))^2/(0.2*theta^2)*(1+0.4*delta)
}
N(1.7)

#############################################################################
########################Modified LKB#########################################
nc=seq(ni[1],ni[2], length.out=20)
X=data.frame(x1,x2,x3,x4,x5,as.numeric(x6),x7,x8,x9,x10,x11)
colnames(X)=c("age","Hormones","T.stage","gleason.score",
"psa","Diabetes","Prostate","statins", "Neoadjuvant","Concurrent","Adjuvant")
X=as.matrix(X)

est=function(nc, X, y){
cvlist=list(100); temp=NULL; Xlist=list()
for(l in 1:length(nc)){
n=nc[l]
geud=rep(0,79)
for(i in 1:79){
d=datalist[[i]][[1]]
geud[i]=(sum(d^(1/n))/length(d))^(n)  ##calculate gEUD
}
Xa=cbind(X,geud);Xlist[[l]]=Xa
cv.lasso=cv.glmnet(Xa,y,alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0)
,family=binomial(link="probit"), nfolds=10)
if(coef(cv.lasso,s="lambda.min")[13]>=0){
lambda=cv.lasso$lambda.min
dev=cv.lasso$glmnet.fit$dev[which(cv.lasso$lambda==lambda)]
temp=c(temp, dev)
cvlist[[l]]=cv.lasso
}else{
cv.lasso=cv.glmnet(Xa,y,alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0),
exclude=12,family=binomial, nfolds=10)
lambda=cv.lasso$lambda.min
dev=cv.lasso$glmnet.fit$dev[which(cv.lasso$lambda==lambda)]
temp=c(temp, dev)
cvlist[[l]]=cv.lasso

}

}
k=which(temp==min(temp))
Q=rep(0,79)
for(i in 1:79){
cv.lasso=cv.glmnet(Xlist[[k]][-i,],y[-i],alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0)
,family="binomial", nfolds=10)
Q[i]=predict(cv.lasso, Xlist[[k]], s="lambda.min")[i]
}
pred=prediction(Q, y)
perf5=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
return(c(nc[k], auc))
}
MLKB=NULL
for(i in 1:50){
MLKB=rbind(MLKB,est(nc, X, R))
}
colnames(MLKB)=c("n","auc")
####statistical model####
auc=NULL
for(i in 1:50){
Xc=cbind(X, z0[,14])
colnames(Xc)=c("age","Hormones","T.stage","gleason.score",
"psa","Diabetes","Prostate","Neoadjuvant","Concurrent","Adjuvant",
"Statins", "D25")
Q=rep(0,79)
for(i in 1:79){
cv.lasso=cv.glmnet(Xc[-i,],R[-i],alpha=1,penalty.factor=c(1,1,1,1,1,1,1,1,1,1,1,0)
,family="binomial", nfolds=10)
Q[i]=predict(cv.lasso, Xc, s="lambda.min")[i]
}
pred=prediction(Q, R)
perf5=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=c(auc,as.numeric(auc.tmp@y.values))
}


#####################################################################
ntcp_cv=rep(0,79)
for(i in 1:79){
  data_2f=data_2m[-i,]
  logit2=glm(R~., data=data_2f, family="binomial")
  ntcp_cv[i]=predict(logit2,newdata=data_2m[i,],type="response")
}
pred=prediction(ntcp_cv, R)
perf=performance(pred,"tpr","fpr")
auc.tmp=performance(pred,"auc")
auc=as.numeric(auc.tmp@y.values)
auc
