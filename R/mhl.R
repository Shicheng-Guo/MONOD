hap_comb_tmp<-function(hapinfo.single){
  # hap_com.single is a character of haplotype, such as "TCTCTT"
  hap_com.single<-c()
  for(i in 1:nchar(hapinfo.single)){
    for(j in i:(nchar(hapinfo.single))){
      hap_com.single<-c(hap_com.single,(substr(hapinfo.single,i,j)))
    } 
  }
  return(hap_com.single)
}

hap_comb<-function(hapinfo){
  # haplotype is array of observed methylation haplotype and return
  hap_comb<-unlist(lapply(hapinfo,hap_comb_tmp))
  return(hap_comb)
}

mhl<-function(hap_comb_array){
  mhl_value<-c()
  hap_comb_array<-unlist(hap_comb_array)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array))){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep("C",hap_comb_array[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array[input_tmp]))
    mhl_value<-c(mhl_value,i*nmeth/(nunmeth+nmeth))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}

mf<-function(hapinfo){
  mf_value<-c()
  hap_comb_array<-unlist(hapinfo)
  meth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="C"
  unmeth_hap=lapply(lapply(hap_comb_array,function(x) unique(unlist(strsplit(x,"")))),function(x) paste(x,collapse=""))=="T"
  hap_comb_array<-c(hap_comb_array[meth_hap],hap_comb_array[unmeth_hap])
  for(i in unique(nchar(hap_comb_array))){
    input_tmp<-which(nchar(hap_comb_array)==i)  
    nmeth<-length(grep("C",hap_comb_array[input_tmp]))
    nunmeth<-length(grep("T",hap_comb_array[input_tmp]))
    mf_value<-c(mhl_value,nmeth/(nunmeth+nmeth))
  }
  mhl_value<-sum(mhl_value)/sum(unique(nchar(hap_comb_array)))
  return(mhl_value)
}


mhlRlt<-c()
for(K in seq(0,100,1)){
hapinfo<-c(rep("CCCC",K),rep("TTTT",100-K))
mhlrlt<-mhl(hap_comb(unique(hapinfo)))
mhlRlt<-c(mhlRlt,mhlrlt)
}
mhlRlt



