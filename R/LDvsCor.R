
LD<-function(vector){
  rlt<-list()
  table<-matrix(table(vector),2,2)
  pAB=table[1,1]/sum(table)
  pA=(2*table[1,1]+table[2,1]+table[1,2])/(2*sum(table))
  pB=(2*table[2,2]+table[2,1]+table[1,2])/(2*sum(table))
  pa=1-pA
  pb=1-pB
  D=pAB-pA*pB
  if(D>0){
    Dmax=min(pA*pb,pa*pB)
  } else{
    Dmax=max(-pA*pB,-pa*pb)
  } 
  Dp=D/Dmax
  r=Dp/sqrt(pA*pa*pB*pb)
  test<-chisq.test(table,correct = T)
  chisq<-test$statistic
  phi=as.numeric(sqrt(test$statistic/length(vector)))
  A1<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))
  A2<-as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))  
  fit<-cor.test(A1,A2)
  rlt$corr=as.numeric(fit$estimate)
  rlt$corr.p=as.numeric(fit$p.value)
  rlt$p=test$p.value
  rlt$Dp=Dp
  rlt$nobs<-sum(table)
  rlt$phi<-phi
  rlt$chisq<-as.numeric(chisq)
  return(rlt)
}
X<-c()
Y<-c()
for(j in 1:3000){
mlc<-c()
Vector<-c()
a<-round(runif(1,1,100))
b<-round(runif(1,1,100))
for(i in 1:100){
  vector<-sample(c(rep("CC",a),rep("CT",b),rep("TC",b),rep("TT",a)),100,replace=T)
  Vector<-c(Vector,vector)
  A1<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,1,1)))))-2))/length(vector)
  A2<-abs(sum(as.numeric(as.factor(unlist(lapply(vector,function(x) substr(x,2,2)))))-2))/length(vector)
  tmp<-data.frame(A1,A2)
  mlc<-rbind(mlc,tmp)
}
x<-LD(Vector)
y<-cor.test(mlc[,1],mlc[,2])
X<-c(X,x$phi)
Y<-c(Y,abs(as.numeric(y$estimate)))
print(j)
}
lm(Y~X)
pdf("LDvsCOR.pdf")
plot(x=X,y=Y,xlab="LD (phi)",ylab="Absolute pearson correlation coefficient (r)",pch=16,col="blue")
abline(a=0,b=1,lwd=3,col="red")
dev.off()

