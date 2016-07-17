require(ForImp)
require(ggplot2)
setwd("C:/SPB_Data/lymphoma")

nrmse<-function(sim,obs){
  numerator<-sum((sim-obs)^2) / length(sim)
  mean<-sum(obs)/length(obs)
  var<-sum((obs-mean)^2)/(length(obs)-1)
  return (sqrt(numerator/var))
}
# set up params
files = c("m_1_1.txt","m_5_1.txt","m_10_1.txt","m_15_1.txt","m_20_1.txt")
ans=read.table("ans.txt",sep="\t",na.strings = "NA")


do<-function(method){
  data = matrix(NA,nrow=1,ncol=5);
  methods <- c('KNN','SKNN','IKNN','ISKNN','LS',"LLS",'SLLS',
               'ILLS','ISLLS',"WLLS",'WSLLS','WILLS','HLLS',
               'shrLLS','shrSLLS','shrILLS','shrWLLS')
  for(i in 1:length(files)){
    m=as.matrix(read.table(files[i],sep="\t",na.strings = "NA"))
    switch(match.arg(method,methods),
      KNN={imputed <-KNN(m,k=15,sim='EuDist')},
      SKNN={imputed <-SKNN(m,k=15,sim='EuDist')},
      IKNN={imputed <-IKNN(m, k=15, sim="EuDist", iter=2)},
      ISKNN={imputed <-ISKNN(m,k=15,sim='EuDist',iter=2)},
      LS={imputed <- LS(m,k=15,sim='EuDist')},
      LLS={imputed <-LLS(m,k=15,sim='EuDist')},
      SLLS={imputed <-SLLS(m,k=15,sim='EuDist')},
      ILLS={imputed <-ILLS(m,k=15,sim='EuDist',iter=2)},
      ISLLS={imputed <-ISLLS(m,k=15,sim='EuDist',iter=2)},
      WLLS={imputed<-WLLS(m,k=15,sim='EuDist',order=1)},
      WSLLS={imputed<-WSLLS(m,k=15,sim='EuDist',order=1)},
      WILLS={imputed<-WILLS(m,k=15,sim='EuDist',iter=2,order=1)},
      HLLS={imputed<-HLLS(m,k=15,sim='EuDist',iter=2,order=1)},
      shrLLS={imputed<-shrLLS(m,k=15,sim='EuDist')},
      shrSLLS={imputed<-shrSLLS(m,k=15,sim='EuDist')},
      shrILLS={imputed<-shrILLS(m,k=15,sim='EuDist',iter=2)},
      shrWLLS={imputed<-shrWLLS(m,k=15,sim='EuDist',order=2)}
    )
    obs<-ans[is.na(m)]
    sim<-imputed[is.na(m)]
    data[1,i] = nrmse(sim,obs)
  }
  return(data)
}
