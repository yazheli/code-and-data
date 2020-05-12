load("credit_data.rdata")
data = credit_data
library(mefa)
library(mnlogit)
library(mlogit)
library(hmeasure)
library(pROC)
library(foreach)
library(doMC)
set.seed(70)
index1 = sample.int(10, size = nrow(data), replace = T)



cv_2 = list()
ptm <- proc.time()
for(j in 1:10){
  credit_data = data[index1 != j,]
  index = sample.int(10, size = nrow(credit_data), replace = T)
  registerDoMC(7)
  auc_2 = foreach(k=1:10) %dopar% {loopfunction(k)}
  #the loopfunction() is the same as the first loopfunction() in taiwan_credit.R
  auc_2 = matrix(unlist(auc_2),nrow = 2)
  cv_2[[j]] = auc_2
  print(j)
}
proc.time() - ptm

#save(cv_2, file = "cv2.rdata")

cv_3 = list()
ptm <- proc.time()
for(j in 1:1){
  credit_data = data[index1 != j,]
  index = sample.int(10, size = nrow(credit_data), replace = T)
  registerDoMC(7)
  auc_3 = foreach(k=1:10) %dopar% {loopfunction(k)}
  #the loopfunction() is the same as the second loopfunction() in taiwan_credit.R
  auc_3 = matrix(unlist(auc_3),nrow = 2)
  cv_3[[j]] = auc_3
  print(j)
}
proc.time() - ptm
#save(cv_3, file = "cv3.rdata")



