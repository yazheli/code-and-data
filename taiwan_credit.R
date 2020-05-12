setwd("~/Desktop/EM PAPER/jcgs-template/code and data")
load("credit_data.rdata")

library(mefa)
library(mnlogit)
library(mlogit)
library(hmeasure)
library(pROC)
library(foreach)

set.seed(70)
index = sample.int(10, size = nrow(credit_data), replace = T)

loopfunction = function(k){
  auc_2 = matrix(0,nrow = 2, ncol = 1)
  
  train = credit_data[index != k,]
  test = credit_data[index == k,]
  
  default = train[which(train$default==1),]   
  nondefault = train[-which(train$default==1),]   
  train = rbind(default,nondefault)
  
  weight_c2 = runif(nrow(default), 0.45, 0.55)
  weight_c3 = 1 - weight_c2
  
  mean2 = sapply(default, weighted.mean,  w = weight_c2)
  mean3 = sapply(default, weighted.mean,  w = weight_c3)
  
  ####initial M step
  new_data = rbind(nondefault,mean2,mean3)
  new_data = as.data.frame(new_data)
  new_data$class = c(rep("c1",nrow(nondefault)), "c2", "c3")
  
  new_mldata = mlogit::mlogit.data(new_data, choice = "class",shape = "wide")
  f <- as.formula(paste("class ~ 1|", paste(names(default)[-31], collapse = " + ")));f
  new_mlogit.model<-mnlogit(f , data = new_mldata, reflevel="c1")
  newQ = logLik(new_mlogit.model)
  oldQ = 0.5*newQ
  coef = as.numeric(coef(new_mlogit.model))
  phi_2 = sum(weight_c2)/nrow(default)
  phi_3 = sum(weight_c3)/nrow(default)
  
  
  #############################################
  #E step updata weights and cluster center
  #test = mlogit::mlogit.data(data, choice = "class",shape = "wide")[c(45001:45120),]
  n_iteration = 50
  weight = matrix(0, nrow = n_iteration, ncol = nrow(default))
  phi = matrix(0,nrow = n_iteration, ncol = 2)
  
  #i=1
  #phi = rep(0,50)

  #for(j in 1 : n_iteration){
  j = 1
  while((abs((newQ-oldQ)/oldQ)>0.0002 | j < 20) & (j < 51)){  
    print((newQ-oldQ)/oldQ)
    print(j)
    oldQ=newQ
    coef = matrix(coef , nrow = 2, byrow = TRUE)
    
    coef1 = coef[1,]
    coef2 = coef[2,]
    
    design_matrix = cbind(rep(1,nrow(train)), train[,-31])
    
    u2 = t(t(design_matrix[1:nrow(default),]) * coef1)
    u3 = t(t(design_matrix[1:nrow(default),]) * coef2)
    u2 = rowSums(u2)
    u3 = rowSums(u3)
    
    p2 = rep(0,nrow(default))
    for(i in 1:nrow(default)){
      p2[i] = exp(u2[i])/(1+exp(u2[i])+exp(u3[i]))
    }
    
    p3 = rep(0,nrow(default))
    for(i in 1:nrow(default)){
      p3[i] = exp(u3[i])/(1+exp(u2[i])+exp(u3[i]))
    }
    
    new_weight_c2 = phi_2*p2/(phi_2*p2 + phi_3*p3)
    new_weight_c3 = phi_3*p3/(phi_2*p2 + phi_3*p3)
    
    new_mu_2 = sapply(default, weighted.mean,  w = new_weight_c2)
    new_mu_3 = sapply(default, weighted.mean,  w = new_weight_c3)
    
    ####M step solve likelihood function
    new_data = rbind(nondefault,new_mu_2,new_mu_3)
    new_data = as.data.frame(new_data)
    new_data$class = c(rep("c1",nrow(nondefault)), "c2", "c3")
    
    
    new_mldata = mlogit::mlogit.data(new_data, choice = "class",shape = "wide")
    f <- as.formula(paste("class ~ 1|", paste(names(default)[-31], collapse = " + ")));f
    new_mlogit.model<-mnlogit(f , data = new_mldata, reflevel="c1")  
    newQ = logLik(new_mlogit.model)
    coef = as.numeric(coef(new_mlogit.model))
    
    #print(logLik(new_mlogit.model))
    
    
    phi_2 = sum(new_weight_c2)/nrow(default)
    phi_3 = sum(new_weight_c3)/nrow(default)
    
    phi[j,] = c(phi_2,phi_3)
    weight[j,] = new_weight_c2
    j = j+1
  }
  
  index_weight = rbind(new_weight_c2,new_weight_c3)
  cluster_index = apply(index_weight, 2, which.max) 
  cluster2 = default[which(cluster_index==1),]
  cluster3 = default[which(cluster_index==2),]
  
  cluster = rbind(nondefault, cluster2 ,cluster3)
  if(length(unique(c(which(colMeans(cluster2) == 0),which(colMeans(cluster3) == 0)))) != 0){
    cluster = cluster[,-unique(c(which(colMeans(cluster2) == 0),which(colMeans(cluster3) == 0)))]
  }
  cluster$class = c(rep("c1",nrow(nondefault)), rep("c2",nrow(cluster2)), rep("c3",nrow(cluster3)))
  
  new_mldata = mlogit::mlogit.data(cluster, choice = "class",shape = "wide")
  f <- as.formula(paste("class ~ 1|", paste(names(cluster)[1:(ncol(cluster)-2)], collapse = " + ")));f
  mlogit.model<-mnlogit(f , data = new_mldata, reflevel="c1")
  
  test$class = c(rep("c1",nrow(test)-2), "c2", "c3")
  test_mldata = mlogit::mlogit.data(test, choice = "class",shape = "wide")
  
  results = predict(mlogit.model, newdata=test_mldata, probability=TRUE)
  result = rowSums(results[,-1])
  
  auc_2[1,1] = as.numeric(auc(roc(response = test$default, predictor = result)))
  
  glm.fit = glm(default ~. , data=train, family=binomial)
  result = predict(glm.fit, newdata=test, probability=TRUE)
  auc_2[2,1] = as.numeric(auc(roc(response = test$default, predictor = result)))
  print(k)
  return(auc_2)
}

library(foreach)
library(doMC)

ptm <- proc.time()
registerDoMC(8)  #change the 8 to your number of CPU cores  
auc_2 = foreach(k=1:10) %dopar% {
  
  loopfunction(k)
  
}
proc.time() - ptm

auc_2 = matrix(unlist(auc_2),nrow = 2)

rowMeans(auc_2)


##################################################################

set.seed(70)
index = sample.int(10, size = nrow(credit_data), replace = T)

loopfunction = function(k){
  auc_3 = matrix(0,nrow = 2, ncol = 1)
  
  train = credit_data[index != k,]
  test = credit_data[index == k,]
  #train_cluster_index = sample(1:20,14)
  #train = credit_data[(index %in% train_cluster_index),]
  #test = credit_data[!(index %in% train_cluster_index),]
  
  default = train[which(train$default==1),]   
  nondefault = train[-which(train$default==1),]   
  data = rbind(default,nondefault)
  
  weight_c2 = runif(nrow(default), 0.3, 0.4)
  weight_c3 = runif(nrow(default), 0.3, 0.4)
  weight_c4 = 1 - weight_c2 - weight_c3
  
  mean2 = sapply(default, weighted.mean,  w = weight_c2)
  mean3 = sapply(default, weighted.mean,  w = weight_c3)
  mean4 = sapply(default, weighted.mean,  w = weight_c4)
  
  
  #r2 = matrix(rep(mean2,times=100), ncol = length(mean2))
  #r3 = matrix(rep(mean2,times=100), ncol = length(mean3))
  ####initial M step
  new_data = rbind(nondefault,mean2,mean3,mean4)
  new_data = as.data.frame(new_data)
  new_data$class = c(rep("c1",nrow(nondefault)), "c2", "c3", "c4")
  
  new_mldata = mlogit::mlogit.data(new_data, choice = "class",shape = "wide")
  f <- as.formula(paste("class ~ 1|", paste(names(default)[-31], collapse = " + ")));f
  new_mlogit.model<-mnlogit(f , data = new_mldata, reflevel="c1")
  newQ = logLik(new_mlogit.model)
  oldQ = 0.5*newQ
  
  coef = as.numeric(coef(new_mlogit.model))
  phi_2 = sum(weight_c2)/nrow(default)
  phi_3 = sum(weight_c3)/nrow(default)
  phi_4 = sum(weight_c4)/nrow(default)
  
  
  #############################################
  ####E step updata weights and cluster center
  #test = mlogit::mlogit.data(data, choice = "class",shape = "wide")[c(45001:45120),]
  n_iteration = 50
  
  weight_c2 = matrix(0, nrow = n_iteration, ncol = nrow(default))
  weight_c3 = matrix(0, nrow = n_iteration, ncol = nrow(default))
  weight_c4 = matrix(0, nrow = n_iteration, ncol = nrow(default))
  
  phi = matrix(0,nrow = n_iteration, ncol = 3)
  
  #i=1
  #phi = rep(0,50)
  j = 1
  while((abs((newQ-oldQ)/oldQ)>0.0002 | j < 20) & (j < 51)){ 
    print((newQ-oldQ)/oldQ)
    print(j)
    oldQ=newQ
    
    coef = matrix(coef , nrow = 3, byrow = TRUE)
    
    coef1 = coef[1,]
    coef2 = coef[2,]
    coef3 = coef[3,]
    
    
    design_matrix = cbind(rep(1,nrow(train)), train[,-31])
    
    u2 = t(t(design_matrix[1:nrow(default),]) * coef1)
    u3 = t(t(design_matrix[1:nrow(default),]) * coef2)
    u4 = t(t(design_matrix[1:nrow(default),]) * coef3)
    
    u2 = rowSums(u2)
    u3 = rowSums(u3)
    u4 = rowSums(u4)
    
    p2 = rep(0,nrow(default))
    for(i in 1:nrow(default)){
      p2[i] = exp(u2[i])/(1+exp(u2[i])+exp(u3[i])+exp(u4[i]))
    }
    
    p3 = rep(0,nrow(default))
    for(i in 1:nrow(default)){
      p3[i] = exp(u3[i])/(1+exp(u2[i])+exp(u3[i])+exp(u4[i]))
    }
    
    p4 = rep(0,nrow(default))
    for(i in 1:nrow(default)){
      p4[i] = exp(u4[i])/(1+exp(u2[i])+exp(u3[i])+exp(u4[i]))
    }
    
    new_weight_c2 = phi_2*p2/(phi_2*p2 + phi_3*p3 + phi_4*p4)
    new_weight_c3 = phi_3*p3/(phi_2*p2 + phi_3*p3 + phi_4*p4)
    new_weight_c4 = phi_4*p4/(phi_2*p2 + phi_3*p3 + phi_4*p4)
    
    new_mu_2 = sapply(default, weighted.mean,  w = new_weight_c2)
    new_mu_3 = sapply(default, weighted.mean,  w = new_weight_c3)
    new_mu_4 = sapply(default, weighted.mean,  w = new_weight_c4)
    
    ####M step solve likelihood function
    new_data = rbind(nondefault,new_mu_2,new_mu_3,new_mu_4)
    new_data = as.data.frame(new_data)
    new_data$class = c(rep("c1",nrow(nondefault)), "c2", "c3","c4")
    
    
    new_mldata = mlogit::mlogit.data(new_data, choice = "class",shape = "wide")
    f <- as.formula(paste("class ~ 1|", paste(names(default)[-31], collapse = " + ")));f
    new_mlogit.model<-mnlogit(f , data = new_mldata, reflevel="c1")  
    newQ = logLik(new_mlogit.model)
    coef = as.numeric(coef(new_mlogit.model))
    
    #print(logLik(new_mlogit.model))
    
    
    phi_2 = sum(new_weight_c2)/nrow(default)
    phi_3 = sum(new_weight_c3)/nrow(default)
    phi_4 = sum(new_weight_c4)/nrow(default)
    
    phi[j,] = c(phi_2,phi_3,phi_4)
    weight_c2[j,] = new_weight_c2
    weight_c3[j,] = new_weight_c3
    weight_c4[j,] = new_weight_c4
    
    j = j+1
  }
  
  index_weight = rbind(new_weight_c2,new_weight_c3,new_weight_c4)
  cluster_index = apply(index_weight, 2, which.max) 
  cluster2 = default[which(cluster_index==1),]
  cluster3 = default[which(cluster_index==2),]
  cluster4 = default[which(cluster_index==3),]
  
  if(nrow(cluster2)*nrow(cluster3)*nrow(cluster4) == 0){
    #auc_3[1,1] = as.numeric(auc(roc(response = test$default, predictor = result)))
    glm.fit = glm(default ~. , data=train, family=binomial)
    result = predict(glm.fit, newdata=test, probability=TRUE)
    auc_3[2,1] = as.numeric(auc(roc(response = test$default, predictor = result)))
    auc_3[1,1] = auc_3[2,1]
  } else {
    cluster = rbind(nondefault, cluster2 ,cluster3, cluster4)
    if(length(unique(c(which(colMeans(cluster2) == 0),which(colMeans(cluster3) == 0),which(colMeans(cluster4) == 0)))) != 0){
      cluster = cluster[,-unique(c(which(colMeans(cluster2) == 0),which(colMeans(cluster3) == 0),which(colMeans(cluster4) == 0)))]   
    }
    
    cluster$class = c(rep("c1",nrow(nondefault)), rep("c2",nrow(cluster2)), rep("c3",nrow(cluster3)),rep("c4",nrow(cluster4)))
    
    new_mldata = mlogit::mlogit.data(cluster, choice = "class",shape = "wide")
    f <- as.formula(paste("class ~ 1|", paste(names(cluster)[1:(ncol(cluster)-2)], collapse = " + ")));f
    mlogit.model<-mnlogit(f , data = new_mldata, reflevel="c1")
    
    test$class = c(rep("c1",nrow(test)-3), "c2", "c3", "c4")
    test_mldata = mlogit::mlogit.data(test, choice = "class",shape = "wide")
    
    results = predict(mlogit.model, newdata=test_mldata, probability=TRUE)
    result = rowSums(results[,-1])
    
    auc_3[1,1] = as.numeric(auc(roc(response = test$default, predictor = result)))
    
    glm.fit = glm(default ~. , data=train, family=binomial)
    result = predict(glm.fit, newdata=test, probability=TRUE)
    auc_3[2,1] = as.numeric(auc(roc(response = test$default, predictor = result)))
  }
  return(auc_3)
}

library(foreach)
library(doMC)

ptm <- proc.time()
registerDoMC(8)  #change to your number of CPU cores  
auc_3 = foreach(k=1:10) %dopar% {
  
  loopfunction(k)
  
}
proc.time() - ptm

auc_3 = matrix(unlist(auc_3),nrow = 2)
rowMeans(auc_3)
