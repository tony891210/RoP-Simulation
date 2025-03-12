library(snpStats)
library(doMC)
library(foreach)
library(lmtest)
registerDoMC(100)

arg<-commandArgs(T)
parameters<-expand.grid(c("additive","dominant","recessive"),c("cis","trans"),c(0,0.8))
m<-parameters[as.numeric(arg[1]),1]
e<-parameters[as.numeric(arg[1]),2]
Dp<-parameters[as.numeric(arg[1]),3]

N=1000
Niter=1000
PA<-0.2
b0=-2
b1=0.5
FreqB<-round(seq(0.05,0.95,length.out=20),2)

Phase.s<-function(PA,PB,Dp,N){
  D<-Dp*min(PA*(1-PB),PB*(1-PA))
  PAB<-PA*PB+D
  Pab<-(1-PA)*(1-PB)+D
  PAb<-PA*(1-PB)-D
  PaB<-(1-PA)*PB-D
  
  Haps<-sample(c("11","10","01","00"),size=2*N,replace=T,prob=c(PAB,PAb,PaB,Pab))
  Hap1<-data.frame(snp1=as.numeric(substring(Haps[c(1:N)],1,1)),
                   snp2=as.numeric(substring(Haps[c(1:N)],2,2)))
  Hap2<-data.frame(snp1=as.numeric(substring(Haps[c((N+1):(2*N))],1,1)),
                   snp2=as.numeric(substring(Haps[c((N+1):(2*N))],2,2)))
  
  G1<-Hap1[,1]+Hap2[,1]
  G2<-Hap1[,2]+Hap2[,2]
  Cis11<-Hap1[,1]*Hap1[,2]+Hap2[,1]*Hap2[,2]
  Cis10<-Hap1[,1]*(1-Hap1[,2])+Hap2[,1]*(1-Hap2[,2])
  Cis01<-(1-Hap1[,1])*Hap1[,2]+(1-Hap2[,1])*Hap2[,2]
  Cis00<-(1-Hap1[,1])*(1-Hap1[,2])+(1-Hap2[,1])*(1-Hap2[,2])
  Trans11<-Hap1[,1]*Hap2[,2]+Hap2[,1]*Hap1[,2]
  Trans10<-Hap1[,1]*(1-Hap2[,2])+Hap2[,1]*(1-Hap1[,2])
  Trans01<-(1-Hap1[,1])*Hap2[,2]+(1-Hap2[,1])*Hap1[,2]
  Trans00<-(1-Hap1[,1])*(1-Hap2[,2])+(1-Hap2[,1])*(1-Hap1[,2])
  
  phase<-data.frame(G1,G2,Cis11,Cis10,Cis01,Cis00,Trans11,Trans10,Trans01,Trans00)
  return(phase)
}
Test.wu<-function(Cis11,Cis10,Cis01,Cis00,Index){
  cis11a<-Cis11[which(Index==1)];cis10a<-Cis10[which(Index==1)];cis01a<-Cis01[which(Index==1)];cis00a<-Cis00[which(Index==1)]
  cis11n<-Cis11[which(Index==0)];cis10n<-Cis10[which(Index==0)];cis01n<-Cis01[which(Index==0)];cis00n<-Cis00[which(Index==0)]
  na<-sum(Index)
  nn<-sum(1-Index)
  
  p11a<-(sum(cis11a)+0.5)/(2*na);p10a<-(sum(cis10a)+0.5)/(2*na);p01a<-(sum(cis01a)+0.5)/(2*na);p00a<-(sum(cis00a)+0.5)/(2*na)
  p11n<-(sum(cis11n)+0.5)/(2*nn);p10n<-(sum(cis10n)+0.5)/(2*nn);p01n<-(sum(cis01n)+0.5)/(2*nn);p00n<-(sum(cis00n)+0.5)/(2*nn)
  
  lambda1<-log((p11a*p00a)/(p10a*p01a));lambda2<-log((p11n*p00n)/(p10n*p01n))
  va<-(1/p11a+1/p10a+1/p01a+1/p00a)/(2*na);vn<-(1/p11n+1/p10n+1/p01n+1/p00n)/(2*nn)
  
  chi<-(lambda1-lambda2)^2/(va+vn)
  p<-1-pchisq(chi,df=1)
  
  return(p)
  
}


power_r_trans<-as.numeric()
power_r_cis<-as.numeric()
power_s<-as.numeric()
power_i<-as.numeric()
#power_je<-as.numeric()
power_r_j<-as.numeric()
power_wu<-as.numeric()
power_h<-as.numeric()
Power<-data.frame()

  
for(PB in FreqB){
  pvalues<-as.numeric()
  pvalues<-foreach(i=c(1:Niter),.combine = "rbind")%dopar%{
    
    set.seed(i)
    Nsample<-round(N/(exp(-2)/(1+exp(-2)))*2)
    
    P<-Phase.s(PA,PB,Dp,Nsample)
    P$G1.D<-(P$G1==1)+0
    P$G2.D<-(P$G2==1)+0
    P$Cis11.D<-(P$Cis11>0)+0
    P$Trans11.D<-(P$Trans11>0)+0
    P$v<-(P$Cis11==1&P$G1==1&P$G2==1)+0
    
    if(e=="cis" & m=="additive"){lg<- b0+b1*P$Cis11}
    if(e=="cis" & m=="recessive"){lg<- b0+2*((P$Cis11>1)+0)}
    if(e=="cis" & m=="dominant"){lg<- b0+b1*((P$Cis11>0)+0)}
    if(e=="trans" & m=="additive"){lg<- b0+b1*P$Trans11}
    #if(e=="trans" & m=="recessive"){lg<- b0+b1*((P$Trans11>1)+0)}
    if(e=="trans" & m=="dominant"){lg<- b0+b1*((P$Trans11>0)+0)}
    
    p<-exp(lg)/(1+exp(lg))
    Y<-rbinom(n=Nsample,size=1,p)
    Index<-c(sample(which(Y==1),N),sample(which(Y==0),N))
    
    Y<-Y[Index]
    P<-P[Index,]
    
    #CASSI JE test
    if(F){
    X<-matrix(c(P$G1+1,P$G2+1),ncol=2)
    colnames(X)<-c("SNP1","SNP2")
    rownames(X)<-paste0("Sample",c(1:dim(X)[1]))
    snpMatrix<-new("SnpMatrix",X)
    pheno<- Y+1
    file_out<-paste0("temp_",e,m,"_",PB,"_",Dp,"_",i)
    write.plink(file.base=file_out,snps=snpMatrix,phenotype=pheno)
    
    system(paste0("./cassi -so -je -je-cc-only -je-th 1 -je-cellmin 0 -i ",file_out,".bed"," -o cassi_",file_out,".out"))
    cassi.out<-read.table(paste0("cassi_",file_out,".out"),header=T)
    p_je<-cassi.out$JE_CC_P
    system(paste0("rm ",file_out,".bim"))
    system(paste0("rm ",file_out,".bed"))
    system(paste0("rm ",file_out,".fam"))
    system(paste0("rm cassi_",file_out,".log"))
    system(paste0("rm cassi_",file_out,".out"))
    }
    #Phase regression
    fit0<-glm(Y~G1+G1.D+G2+G2.D,data=P,family=binomial())
    fit<-glm(Y~G1+G1.D+G2+G2.D+Cis11+Trans11,data=P,family = binomial())
    fit.cis<-glm(Y~G1+G1.D+G2+G2.D+Trans11,data=P,family = binomial())
    p_r_cis<-lrtest(fit.cis,fit)$`Pr(>Chisq)`[2]
    fit.trans<-glm(Y~G1+G1.D+G2+G2.D+Cis11,data=P,family = binomial())
    p_r_trans<-lrtest(fit.trans,fit)$`Pr(>Chisq)`[2]
    p_r_joint<-lrtest(fit,fit0)$`Pr(>Chisq)`[2]
    #Saturate
    fits<-glm(Y~G1+G1.D+G2+G2.D,data=P,family = binomial())
    fits1<-glm(Y~G1+G1.D+G2+G2.D+G1*G2+G1*G2.D+G2*G1.D+G1.D*G2.D+v,data=P,family = binomial())
    p_s<-lrtest(fits1,fits)$`Pr(>Chisq)`[2]
    
    #Interaction
    fiti<-glm(Y~G1+G1.D+G2+G2.D,data=P,family = binomial())
    fiti1<-glm(Y~G1+G1.D+G2+G2.D+G1*G2+G1*G2.D+G2*G1.D+G1.D*G2.D,data=P,family = binomial())
    p_i<-lrtest(fiti1,fiti)$`Pr(>Chisq)`[2]
    #Wu test
    p_wu<-Test.wu(P$Cis11,P$Cis10,P$Cis01,P$Cis00,Index=Y)
    #Hap test
    fith<-glm(Y~1,data=P,family = binomial())
    fith1<-glm(Y~Cis11+Cis10+Cis01,data=P,family = binomial())
    p_h<-lrtest(fith1,fith)$`Pr(>Chisq)`[2]
    return(c(p_r_cis,p_r_trans,p_s,p_r_joint,p_wu,p_h,p_i))
  }
  
  power_r_cis<-append(power_r_cis,sum(pvalues[,1]<0.05)/dim(pvalues)[1])
  power_r_trans<-append(power_r_trans,sum(pvalues[,2]<0.05)/dim(pvalues)[1])
  power_s<-append(power_s,sum(pvalues[,3]<0.05)/dim(pvalues)[1])
  power_r_j<-append(power_r_j,sum(pvalues[,4]<0.05)/dim(pvalues)[1])
  power_wu<-append(power_wu,sum(pvalues[,5]<0.05)/dim(pvalues)[1])
  power_h<-append(power_h,sum(pvalues[,6]<0.05)/dim(pvalues)[1])
  power_i<-append(power_i,sum(pvalues[,7]<0.05)/dim(pvalues)[1])
}
P<-data.frame(FreqB,
              power_r_cis,power_r_trans,power_s,power_r_j,power_wu,power_h,power_i,
              Dp=rep(Dp,length(FreqB)),model=rep(m,length(FreqB)),effect=rep(e,length(FreqB)))

write.table(P,paste0("Power_","PA_",PA,"_",e,"_",m,"_",Dp,".txt"),quote=F,col.names = T,row.names = F)
