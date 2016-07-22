
TopGSIByCategory<-function(gsi,top=1,thresHigh=0.3,thresLow=0.1,plotfigure=T,figprefix="tissuespecific"){
  GSIRlt<-c()
  group<-as.character(unique(gsi$group))
  rank<-c(rep(top,length(group)))
  otf1<-paste(figprefix,"boxplot.pdf",sep="")
  otf2<-paste(figprefix,"heatmap.pdf",sep="")
  parnum<-ceiling(sqrt(length(group)))
  pdf(otf1)
  par(mfrow=c(parnum,parnum),oma = c(2,2,2,2) + 0.1,mar = c(2,2,2,2) + 0.1,cex.axis=0.75, cex.lab=0.75)
  for (i in 1:length(group)){
    # select tissue-specific (remove target group<0.2 or non-target group>0.1)
    subset=gsi[which(gsi$group==group[i]),]
    rexclhigh<-which(apply(subset,1,function(x) x[grep(group[i],colnames(gsi))]<0.2))
    xx<-subset[,-grep(group[i],colnames(gsi))]
    rexcllow<-which(apply(xx,1,function(x) any(as.numeric(x[4:length(x)])>0.1)))
    rexcl<-c(rexclhigh,rexcllow)
    subset=subset[-rexcl,]
    subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
    GSIRlt<-rbind(GSIRlt,subset)
    if(plotfigure==T){
      zz=subset[which(subset$group==group[i]),]
      if(nrow(zz)>=1){
      boxplot(na.omit(zz[,4:ncol(zz)]),horizontal=T,las=2,col="red")
      }else{
      print(paste(group[i],"do not have biomarker,please check raw data",sep=" "))
      }
    }
  }
  dev.off()
  
  if(plotfigure==T){
    HeatMap(data=data.matrix(na.omit(GSIRlt[,4:ncol(GSIRlt)])),phen=gsub("AVE.","",colnames(GSIRlt)[4:ncol(GSIRlt)]),figure=otf2)
  }
  GSIRlt<-na.omit(GSIRlt)
  return(GSIRlt)
}
