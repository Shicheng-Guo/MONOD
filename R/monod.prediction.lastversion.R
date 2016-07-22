load("input.RData")
load("monod.image-1.RData")
# feature selection: same directory principle
x1<-paste(gsirlt1[,1],gsirlt1[,2],sep=".")
x2<-paste(gsirlt2[,1],gsirlt2[,2],sep=".")
x<-c(x1,x2)
y<-names(table(x)[which(table(x)==2)])
region=unique(unlist(lapply(strsplit(y,"[.]"),function(x) x[[1]])))
gSirlt2<-gsirlt2[match(region,gsirlt2[,1]),]  # test
GsiRlt1<-TopGSIByCategory(gsi=gsirlt1,top=300,thresHigh=0.3,thresLow=0.1,allmin=0,plotfigure=F,figprefix="test-tissuespecific")
GsiRlt2<-TopGSIByCategory(gsi=gSirlt2,top=300,thresHigh=0.3,thresLow=0.1,allmin=0,plotfigure=F,figprefix="test-tissuespecific")
GsiRlt1<-GsiRlt1[order(GsiRlt1[,3],decreasing=F),]
A=which(GsiRlt1$group=="Lung")[1:280]
B=which(GsiRlt1$group=="Colon")[1:280]
C=which(GsiRlt1$group=="WBC")[1:280]
GsiRlt1<-GsiRlt1[-c(A,B,C),]
rlt<-rbind(GsiRlt1[1:3],GsiRlt2[1:3])
rlt<-unique(rlt)
write.table(rlt,file="HumanTissueGSI-2.txt",sep="\t",quote=F,col.names=T,row.names=F)
bio<-rlt

# for train data
load("input.RData")
colnames(input)[grep("NCP",colnames(input))]="WBC-test"
colnames(input)[grep("LCP",colnames(input))]="Lung-test"
colnames(input)[grep("CCP",colnames(input))]="Colon-test"
train=input[,-grep("-test",colnames(input))]
test=input[,grep("-test",colnames(input))]
colnames(train)<-unlist(lapply(strsplit(colnames(train),"[.]"),function(x) x[[1]]))
colnames(test)<-unlist(lapply(strsplit(colnames(test),"[-]"),function(x) x[[1]]))

bio<-read.table("HumanTissueGSI-2.txt",head=T,sep="\t",as.is=T)
colnames(bio)<-c("ID","group","GSI")
bio<-topGSIByCategory(bio,top=30)
colnames(bio)<-c("ID","group","GSI")
label=bio[,2]

newtrain<-train[match(bio[,1],rownames(train)),]
newtrain=newtrain[,order(colnames(newtrain))]
label=bio[,2]
ACCU<-c()
for(i in seq(0,0.5,0.01)){
DisMatrix<-apply(newtrain,2,function(x) tapply(x,label,function(x) sum(x>=i)))
DisMatrix
acc<-apply(DisMatrix,2,function(x) which.max(x))
LEN<-table(colnames(newtrain))
LOC<-c(0,cumsum(table(colnames(newtrain))))
ACC<-c()
for(j in 1:(length(LEN))){
ACC1<-sum(acc[(LOC[j]+1):LOC[j+1]]==j)/LEN[j]
ACC<-c(ACC,ACC1)
}
ACCU<-cbind(ACCU,ACC)
}
colnames(ACCU)<-seq(0,0.5,0.01)
ACCU

# for test data (5-fold cross validation)
fold1 <- cut(seq(1,30),breaks=5,labels=FALSE)
fold2 <- cut(seq(1,29),breaks=5,labels=FALSE)
fold3 <- cut(seq(1,75),breaks=5,labels=FALSE)
newtest<-test[match(bio[,1],rownames(test)),]
fold<-c(fold1,fold2,fold3)
traintest<-newtest[,which(fold !=i)]
rlt<-c()
for(i in 1:5){
  ACC<-c()
  for(j in seq(20,300,10)){
    ACC1<-c()
    for(k in seq(0.1,0.4,0.01)){
      bio<-read.table("HumanTissueGSI-2.txt",head=T,sep="\t",as.is=T)
      colnames(bio)<-c("ID","group","GSI")
      bio<-topGSIByCategory(bio,top=j)
      label=bio[,2]
      newtest<-test[match(bio[,1],rownames(test)),]
      fold<-c(fold1,fold2,fold3)
      traintest<-newtest[,which(fold !=i)]
      DisMatrix<-apply(traintest,2,function(x) tapply(x,label,function(x) sum(x>=k)))
      DisMatrix
      acc<-apply(DisMatrix,2,function(x) which.max(x))
      acc1<-sum(acc[grep("Colon",names(acc))]==2)/length(grep("Colon",names(acc)))
      acc2<-sum(acc[grep("Lung",names(acc))]==6)/length(grep("Lung",names(acc)))
      acc3<-sum(acc[grep("WBC",names(acc))]==10)/length(grep("WBC",names(acc)))
      ACCC<-c(j,k,acc1,acc2,acc3)
      ACC1<-rbind(ACC1,ACCC)
    }
      ACC<-rbind(ACC,ACC1)
  }
      parameter<-ACC[which.max(rowSums(ACC[,3:5])),]
      parameter.feature<-parameter[1]
      parameter.thres<-parameter[2]
      bio<-read.table("HumanTissueGSI-2.txt",head=T,sep="\t",as.is=T)
      colnames(bio)<-c("ID","group","GSI")
      bio<-topGSIByCategory(bio,top=parameter.feature)
      label=bio[,2]
      newtest<-test[match(bio[,1],rownames(test)),]
      testtest<-newtest[,which(fold ==i)]
      DisMatrix<-apply(testtest,2,function(x) tapply(x,label,function(x) sum(x>=parameter.thres)))
      DisMatrix
      acc<-apply(DisMatrix,2,function(x) which.max(x))
      acc1<-sum(acc[grep("Colon",names(acc))]==2)/length(grep("Colon",names(acc)))
      acc2<-sum(acc[grep("Lung",names(acc))]==6)/length(grep("Lung",names(acc)))
      acc3<-sum(acc[grep("WBC",names(acc))]==10)/length(grep("WBC",names(acc)))
      ACC.test<-c(parameter[1:5],acc1,acc2,acc3)
      print(ACC.test)
      rlt<-rbind(rlt,ACC.test)
}
write.table(rlt,file="CrossValidation-5Fold.Accuray.result.txt",sep="\t",quote=F,col.names=T,row.names=F)

# for one-time test
  ACC<-c()
  for(j in seq(20,300,10)){
    ACC1<-c()
    for(k in seq(0.1,0.4,0.01)){
      bio<-read.table("HumanTissueGSI-2.txt",head=F,sep="\t",as.is=T)
      colnames(bio)<-c("ID","group","GSI")
      bio<-topGSIByCategory(bio,top=j)
      label=bio[,2]
      newtest<-test[match(bio[,1],rownames(test)),]
      fold<-c(fold1,fold2,fold3)
      traintest<-newtest
      DisMatrix<-apply(traintest,2,function(x) tapply(x,label,function(x) sum(x>=k)))
      DisMatrix
      acc<-apply(DisMatrix,2,function(x) which.max(x))
      acc1<-sum(acc[grep("Colon",names(acc))]==2)/length(grep("Colon",names(acc)))
      acc2<-sum(acc[grep("Lung",names(acc))]==6)/length(grep("Lung",names(acc)))
      acc3<-sum(acc[grep("WBC",names(acc))]==10)/length(grep("WBC",names(acc)))
      ACCC<-c(j,k,acc1,acc2,acc3)
      ACC1<-rbind(ACC1,ACCC)
    }
      ACC<-rbind(ACC,ACC1)
  }
write.table(ACC,file="NoneCrossValidation.Accuray.result.txt",sep="\t",quote=F,col.names=T,row.names=F)

