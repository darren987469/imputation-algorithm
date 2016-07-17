require(MASS)
similarityCal<-function(vec, mat, method="EuDist"){
  epsilon <- 1e-8
  methods <- c("EuDist","cor","Angle")
  switch(match.arg(method,methods),
         EuDist=1/sqrt(rowSums((mat-matrix(vec,nc=length(vec),nr=nrow(mat),byrow=TRUE))^2)+epsilon),
         cor=t(abs(cor(vec,t(mat),use="everything",method="pearson"))),
         Angle=apply(mat,1,function(i) abs(sum(i * vec)/sqrt(sum(i^2)*sum(vec^2))))
  )
}

# Zero
Zero<-function(M) {
  M[is.na(M)] = 0
  return(M)
}

# Row average
RowAverage<-function(M){
  x = M
  for(i in 1:nrow(x)) {
    x[i,which(is.na(x[i,]))] = mean(as.numeric(x[i,which(!is.na(x[i,]))]))
  }
  return(x)
}

# KNN
KNN<-function(M, k=15, sim="EuDist"){
  x = M
  x.missRowIdx <- which(!complete.cases(x))
  x.miss <- x[x.missRowIdx,]
  x.comp <- x[-x.missRowIdx,]
  x.imp <- t(apply(x.miss,1,function(target){
    missidx <- which(is.na(target))
    neighbor.pool <- x.comp[,-missidx,drop=F]
    similarities <- similarityCal(target[-missidx],neighbor.pool,method=sim)
    neighborhood <- order(similarities,decreasing=TRUE)[1:k]
    imp.values <- sapply(missidx,function(j) {
      weight.avg<-similarities[neighborhood] %*% x.comp[neighborhood,j,drop=F] / sum(similarities[neighborhood])
      weight.avg
    })
    target[missidx] <- imp.values
    return(target)
  }))
  x[x.missRowIdx,] <- x.imp
  return(x)
}

# SKNN
SKNN<-function(M, k=15, sim="EuDist"){
  x = M
  miss.RowIdx<-which(!complete.cases(x))
  x.completeRows<-x[-miss.RowIdx,]
  x.missingRows<-x[miss.RowIdx,]
  miss.list<-order(rowSums(is.na(x.missingRows)))
  for (i in seq(nrow(x.missingRows))) {
    target<-x.missingRows[miss.list[i],]
    missColIdx<-which(is.na(target))
    dist.list<-similarityCal(target[-missColIdx], x.completeRows[,-missColIdx], sim)
    neighborsIdx<-order(abs(dist.list),decreasing=T)[1:k]
    estimation<-sapply(missColIdx, function(j){
      weight<-dist.list[neighborsIdx]
      weightedAvg<-weight %*% x.completeRows[neighborsIdx,j] / sum(weight)
      return(weightedAvg)
    })
    target[missColIdx]<-estimation
    x.missingRows[miss.list[i],]<-target
    x.completeRows<-rbind(x.completeRows, target)
  }
  x[miss.RowIdx,]<-x.missingRows
  return(x)
}

# IKNN
IKNN<-function(M, k=15, sim="EuDist", iter=2) {
  x = M
  miss.RowIdx<-which(!complete.cases(x))
  x.ravged<-RowAverage(x)
  x.miss<-(cbind(1:nrow(x), x))[miss.RowIdx,]
  x.complete<-x.ravged
  for (r in 1:iter) {
    x.imputed<-t(apply(x.miss, 1, function(row) {
      rowIdx<-row[1]
      row.origin<-row[-1]
      neighbor.pool<-x.complete[-rowIdx,]
      target<-x.complete[rowIdx,]
      missColIdx<-which(is.na(row.origin))
      dist.list<-similarityCal(target[-missColIdx], neighbor.pool[,-missColIdx], method=sim)
      neighborsIdx<-order(dist.list,decreasing=T)[1:k]
      estimation<-sapply(missColIdx, function(h){
        weight<-dist.list[neighborsIdx]
        weightedAvg<-weight %*% neighbor.pool[neighborsIdx,h]/sum(weight)
        return(weightedAvg)
      })
      row.origin[missColIdx]<-estimation
      return(row.origin)
    }))
    x.complete[miss.RowIdx,]<-x.imputed
  }
  x<-x.complete
  return(x)
}

# ISKNN
ISKNN<-function(M, k=15, sim="EuDist", iter=2) {
  x = M
  x.missRowIdx<-which(!complete.cases(x))
  x.ravged<-RowAverage(x)
  x.miss<-(cbind(1:nrow(x), x))[x.missRowIdx,]
  x.complete<-SKNN(M=x, k=k,sim=sim)
  for (r in 1:iter) {
    x.imputed<-t(apply(x.miss, 1, function(row) {
      rowIdx<-row[1]
      row.origin<-row[-1]
      neighbor.pool<-x.complete[-rowIdx,]
      target<-x.complete[rowIdx,]
      missColIdx<-which(is.na(row.origin))
      dist.list<-similarityCal(target[-missColIdx], neighbor.pool[,-missColIdx], method=sim)
      neighborsIdx<-order(dist.list,decreasing=T)[1:k]
      estimation<-sapply(missColIdx, function(h){
        weight<-dist.list[neighborsIdx]
        weightedAvg<-weight %*% neighbor.pool[neighborsIdx,h]/sum(weight)
        return(weightedAvg)
      })
      row.origin[missColIdx]<-estimation
      return(row.origin)
    }))
    x.complete[x.missRowIdx,]<-x.imputed
  }
  x<-x.complete
  return(x)
}

# LS
LS<-function(M, k=15, sim="EuDist"){
  x = M
  x.missRowIdx<-which(!complete.cases(x))
  x.completeRows<-x[-x.missRowIdx,]
  x.missingRows<-x[x.missRowIdx,]
  x.imputed<-t(apply(x.missingRows,1,function(i){
    missColIdx<-which(is.na(i))
    dist.list<-similarityCal(i[-missColIdx], x.completeRows[,-missColIdx,drop=F], method=sim)
    neighborsIdx<-order(dist.list,decreasing=T)[1:k]
    neighbor.sim<-dist.list[neighborsIdx]
    e=1e-6
    weight<-((neighbor.sim^2)/(1-(neighbor.sim)^2+e))^2
    fit.list<-lapply(neighborsIdx, function(m){
      fit<-lm(i[-missColIdx] ~ x.completeRows[m,-missColIdx])
      if (is.na(fit$coefficients[[2]])) fit$coefficients[[2]] <- 0
      reg<-c(fit$coefficients[[1]],fit$coefficients[[2]])
      return(reg)
    })
    estimation<-sapply(missColIdx,function(j){
      regression<-rep(0,k)
      for (n in seq(k)) {
        a0<-fit.list[n][[1]][1]
        b0<-fit.list[n][[1]][2]
        regression[n]<-a0+b0*x.completeRows[neighborsIdx[n],j]
      }
      weightedAvg<-(weight %*% regression)/sum(weight)
      return(weightedAvg)
    })
    i[missColIdx]<-estimation
    return(i)
  }))
  x[x.missRowIdx,]<-x.imputed
  return(x)
}

# LLS
LLS<-function(M, k=15, sim="EuDist"){
  x = M
  missRowIdx<-which(!complete.cases(x))
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed<-t(apply(x.missingPart,1,function(i){
    missColIdx<-which(is.na(i))
    dist.list<-similarityCal(i[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    Cp<-ginv(t(A))
    w<-i[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    ans<-t(b) %*% X
    i[missColIdx]<-ans
    return(i)
  }))
  x[missRowIdx,]<-imputed
  return(x)
}

# SLLS
SLLS <- function(M, k=15, sim){
  x = M
  missRowIdx<-which(!complete.cases(x))
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed = x.missingPart
  miss.list<-order(rowSums(is.na(x.missingPart)))
  # thershold = average row missing gene num
  threshold = sum(is.na(x.missingPart)) / nrow(x.missingPart)
  
    for(j in 1: nrow(x.missingPart)){
    target = x.missingPart[miss.list[j],]
    missColIdx<-which(is.na(target))
    dist.list<-similarityCal(target[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    Cp<-ginv(t(A))
    w<-target[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    ans<-t(b) %*% X
    # impute value
    target[missColIdx]<-ans
    # add imputed row to x.completepart if miss gene num < threshold
    if(length(missColIdx) < threshold){
      x.completePart <- rbind(x.completePart, target)
    }
    imputed[miss.list[j],] = target 
  }
  x[missRowIdx,]<-imputed
  return(x)
}

# ILLS
ILLS<-function(M, k=15, sim="EuDist", iter=2){
  x = M
  missRowIdx<-which(!complete.cases(x))
  # x.comp imputed with rowaverage
  x.comp<-RowAverage(M=x)
  x.miss<-(cbind(1:nrow(x),x))[missRowIdx,]
  for(iteration in 1:iter){
    imputed<-t(apply(x.miss,1,function(i){
      rowIdx<-i[1]
      i.origin<-i[-1]
      candidates<-x.comp[-rowIdx,]
      target<-x.comp[rowIdx,]
      missColIdx<-which(is.na(i.origin))
      dist.list<-similarityCal(target[-missColIdx],candidates[,-missColIdx],sim)
      neighborIdx<-order(dist.list,decreasing=T)[1:k]
      A<-candidates[neighborIdx,-missColIdx,drop=FALSE]
      b<-candidates[neighborIdx,missColIdx,drop=FALSE]
      Cp<-ginv(t(A))
      w<-i.origin[-missColIdx,drop=FALSE]
      X<-Cp %*% w
      ans<-t(b) %*% X
      i.origin[missColIdx]<-ans
      return(i.origin)
    }))
    x.comp[missRowIdx,]<-imputed
  }
  x <- x.comp
  return(x)
}

# ISLLS
ISLLS<-function(M, k=15, sim="EuDist", iter=2){
  x = M
  missRowIdx<-which(!complete.cases(x))
  # x.comp imputed with SLLS
  x.comp<-SLLS(x,k,sim)
  x.miss<-(cbind(1:nrow(x),x))[missRowIdx,]
  for(iteration in 1:iter){
    imputed<-t(apply(x.miss,1,function(i){
      rowIdx<-i[1]
      i.origin<-i[-1]
      candidates<-x.comp[-rowIdx,]
      target<-x.comp[rowIdx,]
      missColIdx<-which(is.na(i.origin))
      dist.list<-similarityCal(target[-missColIdx],candidates[,-missColIdx],sim)
      neighborIdx<-order(dist.list,decreasing=T)[1:k]
      A<-candidates[neighborIdx,-missColIdx,drop=FALSE]
      b<-candidates[neighborIdx,missColIdx,drop=FALSE]
      Cp<-ginv(t(A))
      w<-i.origin[-missColIdx,drop=FALSE]
      X<-Cp %*% w
      ans<-t(b) %*% X
      i.origin[missColIdx]<-ans
      return(i.origin)
    }))
    x.comp[missRowIdx,]<-imputed
  }
  x <- x.comp
  return(x)
}

# WLLS
WLLS <- function(M, k=15, sim='EuDist', order=1){
  x = M
  missidx<-is.na(x)
  missRowIdx<-which(rowSums(missidx) != 0)
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed = x.missingPart
  for(j in 1: nrow(x.missingPart)){
    i = x.missingPart[j,]
    missColIdx<-which(is.na(i))
    dist.list<-similarityCal(i[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    similarity <-dist.list[neighborIdx]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    #weight matrix S
    S<-matrix(0,ncol=k,nrow=k)
    for(z in 1:k){
      S[z,z] = (similarity[z])^order 
    }
    A=S%*%A
    b=S%*%b
    Cp<-ginv(t(A))
    w<-i[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    ans<-t(b) %*% X
    # impute value
    i[missColIdx]<-ans
    imputed[j,] = i 
  }
  x[missRowIdx,]<-imputed
  return(x)
}

# WSLLS
WSLLS <- function(M, k=15, sim, order=1){
  x = M
  missidx<-is.na(x)
  missRowIdx<-which(rowSums(missidx) != 0)
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed = x.missingPart
  miss.list<-order(rowSums(is.na(x.missingPart)))
  # thershold = average row missing gene num
  threshold = sum(is.na(x.missingPart)) / nrow(x.missingPart)
  
  for(j in 1: nrow(x.missingPart)){
    target = x.missingPart[miss.list[j],]
    missColIdx<-which(is.na(target))
    dist.list<-similarityCal(target[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    similarity <-dist.list[neighborIdx]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    #weight matrix S
    S<-matrix(0,ncol=k,nrow=k)
    for(z in 1:k){
      S[z,z] = (similarity[z])^order 
    }
    A=S%*%A
    b=S%*%b
    Cp<-ginv(t(A))
    w<-target[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    ans<-t(b) %*% X
    # impute value
    target[missColIdx]<-ans
    # add imputed row to x.completepart if miss gene num < threshold
    if(length(missColIdx) < threshold){
      x.completePart <- rbind(x.completePart, target)
    }
    imputed[miss.list[j],] = target 
  }
  x[missRowIdx,]<-imputed
  return(x)
}

# WILLS
WILLS<-function(M, k=15, sim="EuDist", iter=2, order=1){
  x = M
  missRowIdx<-which(!complete.cases(x))
  # x.comp imputed with rowaverage
  x.comp<-RowAverage(M=x)
  x.miss<-(cbind(1:nrow(x),x))[missRowIdx,]
  for(iteration in 1:iter){
    imputed<-t(apply(x.miss,1,function(i){
      rowIdx<-i[1]
      i.origin<-i[-1]
      candidates<-x.comp[-rowIdx,]
      target<-x.comp[rowIdx,]
      missColIdx<-which(is.na(i.origin))
      dist.list<-similarityCal(target[-missColIdx],candidates[,-missColIdx],sim)
      neighborIdx<-order(dist.list,decreasing=T)[1:k]
      A<-candidates[neighborIdx,-missColIdx,drop=FALSE]
      b<-candidates[neighborIdx,missColIdx,drop=FALSE]
      #weight matrix S
      similarity <-dist.list[neighborIdx]
      S<-matrix(0,ncol=k,nrow=k)
      for(z in 1:k){
        S[z,z] = (similarity[z])^order 
      }
      A=S%*%A
      b=S%*%b
      Cp<-ginv(t(A))
      w<-i.origin[-missColIdx,drop=FALSE]
      X<-Cp %*% w
      ans<-t(b) %*% X
      i.origin[missColIdx]<-ans
      return(i.origin)
    }))
    x.comp[missRowIdx,]<-imputed
  }
  x <- x.comp
  return(x)
}

# HLLS
HLLS<-function(M, k=15, sim="EuDist", iter=2, order=1){
  x = M
  missRowIdx<-which(!complete.cases(x))
  # x.comp imputed with WSLLS
  x.comp<-WSLLS(x,k,sim,order)
  x.miss<-(cbind(1:nrow(x),x))[missRowIdx,]
  # iteration: WLLS
  for(iteration in 1:iter){
    imputed<-t(apply(x.miss,1,function(i){
      rowIdx<-i[1]
      i.origin<-i[-1]
      candidates<-x.comp[-rowIdx,]
      target<-x.comp[rowIdx,]
      missColIdx<-which(is.na(i.origin))
      dist.list<-similarityCal(target[-missColIdx],candidates[,-missColIdx],sim)
      neighborIdx<-order(dist.list,decreasing=T)[1:k]
      A<-candidates[neighborIdx,-missColIdx,drop=FALSE]
      b<-candidates[neighborIdx,missColIdx,drop=FALSE]
      #weight matrix S
      similarity <-dist.list[neighborIdx]
      S<-matrix(0,ncol=k,nrow=k)
      for(z in 1:k){
        S[z,z] = (similarity[z])^order 
      }
      A=S%*%A
      b=S%*%b
      Cp<-ginv(t(A))
      w<-i.origin[-missColIdx,drop=FALSE]
      X<-Cp %*% w
      ans<-t(b) %*% X
      i.origin[missColIdx]<-ans
      return(i.origin)
    }))
    x.comp[missRowIdx,]<-imputed
  }
  x <- x.comp
  return(x)
}

# shrLLS
ShrLLS <- function(M, k=15, sim='EuDist'){
  x = M
  missRowIdx<-which(!complete.cases(x))
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed = x.missingPart
  
  for(j in 1: nrow(x.missingPart)){
    i = x.missingPart[j,]
    missColIdx<-which(is.na(i))
    dist.list<-similarityCal(i[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    Cp<-ginv(t(A))
    w<-i[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    #shrinkage
    shrink = 1-((nrow(A)-2)*var(X))/(ncol(A)*sum(X^2))
    x.shrink = as.numeric(shrink)*X
    ans<-t(b) %*% x.shrink
    # impute value
    i[missColIdx]<-ans
    imputed[j,] = i 
  }
  x[missRowIdx,]<-imputed
  return(x)
}

# shrSLLS
ShrSLLS <- function(M, k=15, sim){
  x = M
  missRowIdx<-which(!complete.cases(x))
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed = x.missingPart
  miss.list<-order(rowSums(is.na(x.missingPart)))
  # thershold = average row missing gene num
  threshold = sum(is.na(x)) / nrow(x.missingPart)
  
  for(j in 1: nrow(x.missingPart)){
    target = x.missingPart[miss.list[j],]
    missColIdx<-which(is.na(target))
    dist.list<-similarityCal(target[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    Cp<-ginv(t(A))
    w<-target[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    #shrinkage
    shrink = 1-((nrow(A)-2)*var(X))/(ncol(A)*sum(X^2))
    x.shrink = as.numeric(shrink)*X
    ans<-t(b) %*% x.shrink
    # impute value
    target[missColIdx]<-ans
    # add imputed row to x.completepart if miss gene num < threshold
    if(length(missColIdx) < threshold){
      x.completePart <- rbind(x.completePart, target)
    }
    imputed[miss.list[j],] = target 
  }
  x[missRowIdx,]<-imputed
  return(x)
}

# shrILLS
ShrILLS<-function(M, k=15, sim="EuDist", iter=2){
  x = M
  missRowIdx<-which(!complete.cases(x))
  # x.comp imputed with rowaverage
  x.comp<-RowAverage(M=x)
  x.miss<-(cbind(1:nrow(x),x))[missRowIdx,]
  for(iteration in 1:iter){
    imputed<-t(apply(x.miss,1,function(i){
      rowIdx<-i[1]
      i.origin<-i[-1]
      candidates<-x.comp[-rowIdx,]
      target<-x.comp[rowIdx,]
      missColIdx<-which(is.na(i.origin))
      dist.list<-similarityCal(target[-missColIdx],candidates[,-missColIdx],sim)
      neighborIdx<-order(dist.list,decreasing=T)[1:k]
      A<-candidates[neighborIdx,-missColIdx,drop=FALSE]
      b<-candidates[neighborIdx,missColIdx,drop=FALSE]
      Cp<-ginv(t(A))
      w<-i.origin[-missColIdx,drop=FALSE]
      X<-Cp %*% w
      #shrinkage
      shrink = 1-((nrow(A)-2)*var(X))/(ncol(A)*sum(X^2))
      x.shrink = as.numeric(shrink)*X
      ans<-t(b) %*% x.shrink
      i.origin[missColIdx]<-ans
      return(i.origin)
    }))
    x.comp[missRowIdx,]<-imputed
  }
  x <- x.comp
  return(x)
}

# shrWLLS
ShrWLLS <- function(M, k=15, sim='EuDist', order=1){
  x = M
  missidx<-is.na(x)
  missRowIdx<-which(rowSums(missidx) != 0)
  x.completePart<-x[-missRowIdx,]
  x.missingPart<-x[missRowIdx,]
  imputed = x.missingPart
  for(j in 1: nrow(x.missingPart)){
    i = x.missingPart[j,]
    missColIdx<-which(is.na(i))
    dist.list<-similarityCal(i[-missColIdx],x.completePart[,-missColIdx],sim)
    neighborIdx<-order(dist.list,decreasing=T)[1:k]
    similarity <-dist.list[neighborIdx]
    A<-x.completePart[neighborIdx,-missColIdx,drop=FALSE]
    b<-x.completePart[neighborIdx,missColIdx,drop=FALSE]
    #weight matrix S
    S<-matrix(0,ncol=k,nrow=k)
    for(z in 1:k){
      S[z,z] = (similarity[z])^order 
    }
    A=S%*%A
    b=S%*%b
    Cp<-ginv(t(A))
    w<-i[-missColIdx,drop=FALSE]
    X<-Cp %*% w
    #shrinkage
    shrink = 1-((nrow(A)-2)*var(X))/(ncol(A)*sum(X^2))
    x.shrink = as.numeric(shrink)*X
    ans<-t(b) %*% x.shrink
    # impute value
    i[missColIdx]<-ans
    imputed[j,] = i 
  }
  x[missRowIdx,]<-imputed
  return(x)
}