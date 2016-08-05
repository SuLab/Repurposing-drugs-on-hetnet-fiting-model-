#############################################################
#
# Rephetio: Repurposing drugs on a hetnet (https://thinklab.com/p/rephetio)
# data from : https://github.com/dhimmel/learn/tree/e5e7459532944cafcf542ef6ddfe9184548224d9/all-features/data/matrix/rephetio-v2.0
# this R code file wrote by Siyue Wu
# supervisor: Andrew Su, Benjamin Good
# special thanks to Nuria Queralt Rosinach
# Aug 4, 2016
#
#############################################################
#
#############################################################
# import library or packages
#############################################################
library(VIM)
library(data.table)
library(dplyr)
library(randomForest)
library(xgboost)
library(corrplot)
library(caret)
library(ROCR)
library(caTools)
library(doMC)
library(glmnet)
#############################################################
#
# function that used in caculate roc and pre 
# and plot the auroc and auprc
#
#############################################################
aurprc=function(testY,yhat) {
  ##### caculate auroc
  # tpr: true positive rate
  # fpr: false positive rate
  rocr_pred=prediction(predictions=yhat, labels=testY)
  auroc=performance(rocr_pred, 'auc')@y.values[[1]]
  plotRoc=performance(rocr_pred,"tpr","fpr")

  #####
  # caculate auprc
  # rec: True positive rate
  # prec:Positive predictive value
  recallprecision <- data.frame(
    'recall'=performance(rocr_pred, measure='rec')@y.values[[1]],
    'precision'=performance(rocr_pred, measure='prec')@y.values[[1]]
  )
  trapz_df <- na.omit(recallprecision[, c('recall', 'precision')])
  auprc=trapz(trapz_df$recall, trapz_df$precision)
  plotPrc=performance(rocr_pred,"rec","prec")
  matrics <- list(
    'auroc'=auroc,
    'auprc'=auprc,
    'plotRoc'=plotRoc,
    'plotPrc'=plotPrc
  )
  return(matrics)
}

## function
## glmnet_trainq function is changed a litter bit from 
## https://github.com/dhimmel/hetior/blob/master/R/glmnet.R
glmnet_train1 <- function(X, y, testX, testY, alpha=0,w, s='lambda.1se', cores=7, seed=0, ...) {
  # Fit a regularized logistic regression model using the glmnet package.
  # alpha is the regularization parameter (0 for ridge, 1 for lasso).
  
  fit <- list(X=X, y=y, w=w, s=s, alpha=alpha, seed=seed)
  
  # train model
  doMC::registerDoMC(cores=cores)
  set.seed(seed)
  fit$cv_model <- cv.glmnet(x = data.matrix(X), y = y, family='binomial',
                                    alpha=alpha, parallel=TRUE)
  fit$lambda <- fit$cv_model[[s]]
  # model information and performance
  fit$y_pred <- glmnet_predict(fit$cv_model, data.matrix(testX), s=s)
  fit$vtm <- aurprc(testY, fit$y_pred)
  
  return(fit)
}

glmnet_predict <- function(cv_glmnet, X, s = 'lambda.1se') {
  y_pred <- as.numeric(
    predict(cv_glmnet, s=cv_glmnet[[s]], newx=X, type='response'))
  return(y_pred)
}
 
#############################################################
#
# load matrix and do visulation for the data
#
#############################################################
# load features matrix which is the all feature matrix but haven't transformed
features <- fread("~/Downloads/features.tsv", stringsAsFactors = TRUE)
features=data.frame(features)
# load transfered-features matrix which is the all feature matrix and also have transformed
transformed.features <- fread("~/Downloads/transformed-features.tsv", stringsAsFactors=TRUE)
tran_df=data.frame(select(transformed.features,-matches('[CD][pt][CD]',ignore.case=F)))
# select the prior+degree+dwpc matrix from transform-features
X_dwpc = as.matrix(select(tran_df,prior_logit, starts_with('degree'), starts_with('dwpc')))
dim(X_dwpc)
# select the prior+degree+pdwpc matrix from transform-features
X_pdwpc = as.matrix(select(tran_df, prior_logit, starts_with('degree'), starts_with('pdwpc')) )
dim(X_pdwpc)
# select the prior+degree+rdwpc matrix from transform-features
X_rdwpc = as.matrix(select(tran_df, prior_logit, starts_with('degree'), starts_with('rdwpc')) )
dim(X_rdwpc)
# select the prior+degree+pdwpc and rdwpc matrix from transfered-features
X_split = as.matrix(select(tran_df, prior_logit, starts_with('degree'), starts_with('pdwpc'), starts_with('rdwpc')) )
dim(X_split)
# select 
X_all = data.frame(select(tran_df, prior_logit, starts_with('degree'), contains('dwpc')))
dim(X_all)

#########start do the visulation for this data matrix
## to see the pattern for 0 value in the matrix
#########
# plot the 0 pattern for features matrix(not transform)
simplefeatures=features[,-c(1:7)]
miceplot0 <- aggr(ifelse(simplefeatures==0, NA, 1), col=c("skyblue","pink"),
                  numbers=F,combined=T, varheight=F, border=NA,
                  sortVars=F, sortCombs=T, bars=F, ylabs=c("observation (3775 pairs of compound and disease)"),
                  labels=names(simplefeatures), cex.axis=.5)
legend("topleft", legend = c("value not equal 0", "value equal 0"), fill = c("skyblue","pink"),cex = 1,bty ="n")
# plot the 0 pattern for tranform features including prior+degree+X_dwpc
miceplot0 <- aggr(ifelse(X_dwpc==0, NA, 1), col=c("skyblue","pink"),
                  numbers=F,combined=T, varheight=F, border=NA,
                  sortVars=F, sortCombs=T, bars=F, ylabs=c("observation (3775 pairs of compound and disease)"),
                  labels=names(X_dwpc), cex.axis=.5)
legend("topleft", legend = c("value not equal 0", "value equal 0"), fill = c("skyblue","pink"),cex = 1,bty ="n")
# plot the 0 pattern for tranformed features matrix including prior+degree+X_pdwpc
miceplot0 <- aggr(ifelse(X_pdwpc==0, NA, 1), col=c("skyblue","pink"),
                  numbers=F,combined=T, varheight=F, border=NA,
                  sortVars=F, sortCombs=T, bars=F, ylabs=c("observation (3775 pairs of compound and disease)"),
                  labels=names(X_pdwpc), cex.axis=.5)
legend("topleft", legend = c("value not equal 0", "value equal 0"), fill = c("skyblue","pink"),cex = 1,bty ="n")
# plot the 0 pattern for tranformed features matrix including prior+degree+X_rdwpc
miceplot0 <- aggr(ifelse(X_rdwpc==0, NA, 1), col=c("skyblue","pink"),
                  numbers=F,combined=T, varheight=F, border=NA,
                  sortVars=F, sortCombs=T, bars=F, ylabs=c("observation (3775 pairs of compound and disease)"),
                  labels=names(X_rdwpc), cex.axis=.5)
legend("topleft", legend = c("value not equal 0", "value equal 0"), fill = c("skyblue","pink"),cex = 1,bty ="n")

############ do the correlation between column
# first for dwpc
########
par(mfrow=c(1,1))
#####
# draw the correlation between degree columns
p0=data.frame(select(data.frame(X_dwpc), starts_with('degree')))
c0=cor(p0, use = "complete")
# the cut off of high correlation is 0.80
rem0=findCorrelation(c0, cutoff = .80)
p0.rem=cor(p0[,rem0] , use="complete")
p0.plot=corrplot(p0.rem,  method = 'ellipse' , type="lower" , 
                 order="FPC",tl.col = "black")
#####
# draw the correlation between dwpc columns
p1=data.frame(select(data.frame(X_dwpc), starts_with('dwpc')))
c1=cor(p1, use = "complete")
# the cut off =0.90 means only get out of the correlation >=0.90
rem1=findCorrelation(c1, cutoff = .95)
# since there are lots of features have correlation >=0.90 (length(rem1)=586)
# I use for loop to plot all out, each plot has only 10 columns, so we total have
# ceiling(length(rem1)/10) output graphs for high correlation columns(correlation>0.9)
j=30
for(i in 1:ceiling(length(rem1)/j)){
  number=pmin(j*i,length(rem1))
  p1.rem=cor(p1[,rem1[((i-1)*j):number]] , use="complete")
  p1.plot=corrplot(p1.rem,  method = 'ellipse' , type="lower" , 
                   order="alphabet",tl.cex=0.6)
}
#############################################################
#
# try logistic regression for X_dwpc
#
#############################################################
# permute the data and divide 70% of data to be training data
# and 30% of data to be testing data
index=sample(1:nrow(X_dwpc),0.7*nrow(X_dwpc),replace=F)
trainingX=X_dwpc[index,]
trainingY=tran_df$status[index]
testX=X_dwpc[-index,]
testY=tran_df$status[-index]
# call the function above
penalty = ifelse(colnames(trainingX) == 'prior_logit', 0, 1)

fit = glmnet_train1(trainingX,trainingY,testX,testY, tran_df$status, alpha = 0, penalty.factor=penalty, cores=10)
yhat=fit$y_pred
yhat=pmax(pmin(yhat,1-10^-15),10^-15)
glmnetlogloss=-mean(testY*log(yhat,base=exp(1))+(1-testY)*log(1-yhat,base=exp(1)))
glmnetAuroc=fit$vtm['auroc'][[1]]
glmnetAuprc=fit$vtm['auprc'][[1]]
cat("AUROC for glmnet model= ",glmnetAuroc)
cat("AUPRC for glmnet model= ",glmnetAuprc)
glmnetPlotPrc=fit$vtm['plotPrc'][[1]]
glmnetPlotRoc=fit$vtm['plotRoc'][[1]]
#############################################################
#
# try random forest for X_dwpc
#
#############################################################
# permute the data and divide 70% of data to be training data
# and 30% of data to be testing data
index=sample(1:nrow(X_dwpc),0.7*nrow(X_dwpc),replace=F)
trainingX=X_dwpc[index,]
trainingY=tran_df$status[index]
testX=X_dwpc[-index,]
testY=tran_df$status[-index]
# start the model 
rf=randomForest(as.factor(trainingY) ~., data=trainingX, ntree=150)
# the runing time is around 3 min
Sys.time()-start
yhat=predict(rf, testX, type="prob")[,2]
yhat=pmax(pmin(yhat,1-10^-15),10^-15)
# use log loss function to see the test error
randmforestlogloss=-mean(testY*log(yhat,base=exp(1))+(1-testY)*log(1-yhat,base=exp(1)))
# the log loss=0.15
####### to see the AUROC and AUPRC
randomForestAuroc=aurprc(testY,yhat)$auroc
randomForestAuprc=aurprc(testY,yhat)$auprc
cat("AUROC for random forest model= ",aurprc(testY,yhat)$auroc)
cat("AUPRC for random forest model= ",aurprc(testY,yhat)$auprc)
randomForestPlotPrc=aurprc(testY,yhat)$plotPrc
randomForestPlotRoc=aurprc(testY,yhat)$plotRoc
#############################################################
#
# try xgboosting for X_dwpc
#
#############################################################
# permute the data and divide 70% of data to be training data
# and 30% of data to be testing data
index=sample(1:nrow(X_dwpc),0.7*nrow(X_dwpc),replace=F)
trainingX=X_dwpc[index,]
trainingY=tran_df$status[index]
testX=X_dwpc[-index,]
testY=tran_df$status[-index]
# create cgb.DMatrix for using at xgboost package
xgtrain = xgb.DMatrix(as.matrix(trainingX), label = trainingY, missing=NA)
xgtest = xgb.DMatrix(as.matrix(testX), missing=NA)
# cv function for gradience boosting
# https://github.com/dmlc/xgboost/blob/master/doc/parameter.md
# eval_metrix="logloss"
param0 = list("objective"  = "binary:logistic", "eval_metric" = "logloss",
              "eta" = 0.05, "subsample" = 0.9, "colsample_bytree" = 0.9, "max_depth" = 10)
iter=500
model_cv = xgb.cv(params = param0, nrounds = iter, nfold = 2, data = xgtrain, #early.stop.round = 10, 
                  maximize = FALSE)# nthread = 8)
gc()
bestIter <- which(model_cv$test.logloss.mean==min(model_cv$test.logloss.mean))


ensemble=rep(0, nrow(testX))
Iter=round((bestIter-1)[1] * 1.5)
for (i in 1:3) {
  print(i)
  set.seed(i + 2017)
  model = xgb.train(nrounds = Iter, params = param0, data = xgtrain, 
                    watchlist = list('train' = xgtrain), print.every.n = 10, nthread =8)
  p <- predict(model, xgtest)
  rm(model)
  gc()
  ensemble= ensemble + p
}
finalprob= ensemble/i
yhat=pmax(pmin(finalprob,1-10^-15),10^-15)
length(testY)
length(yhat)
##loss log
xgboostlogloss=-mean(testY*log(yhat,base=exp(1))+(1-testY)*log(1-yhat,base=exp(1)))
#######to see the AUROC and AUPRC
xgboostAuroc=aurprc(testY,yhat)$auroc
xgboostAuprc=aurprc(testY,yhat)$auprc
cat("AUROC for gradient boosting model= ",aurprc(testY,yhat)$auroc)
cat("AUPRC for gradient boosting model= ",aurprc(testY,yhat)$auprc)
xgboostPlotPrc=aurprc(testY,yhat)$plotPrc
xgboostPlotRoc=aurprc(testY,yhat)$plotRoc

#############################################################
#
# compare the aur and prc curve with logistic regression
#
#############################################################
#log loss
cat("logloss for glmnet model= ", glmnetlogloss)
cat("logloss for random forest model= ", randmforestlogloss)
cat("logloss for gradient boosting model= ", xgboostlogloss)
# auroc
cat("AUROC for glmnet model= ",glmnetAuroc)
cat("AUROC for random forest model= ",randomForestAuroc)
cat("AUROC for gradient boosting model= ",xgboostAuroc)
# auprc
cat("AUPRC for glmnet model= ",glmnetAuprc)
cat("AUPRC for random forest model= ",randomForestAuprc)
cat("AUPRC for gradient boosting model= ",xgboostAuprc)

# plot roc curve 
plot(glmnetPlotRoc,col="red",lty=1, lwd=2,main="ROC curves")
plot(randomForestPlotRoc,col="blue",lty=1, lwd=2,add=T)
plot(xgboostPlotRoc,col="green",lty=1, lwd=2,add=T)
legend(0.5,0.2, legend = c("logisticRegression", "randomforest","gradientBoosting"), lty=1,
       lwd=2,col=c("red","blue" ,"green"),cex = 1,bty ="n")
legend(0.3,0.45,c(paste(c("AUROC for logistic = ","AUROC for randomforest = ",
                         "AUROC for gradientboosting = "),
                       c(round(glmnetAuroc, digits = 2),round(randomForestAuroc,digits=2),
                         round(xgboostAuroc,digits=2)),sep=""),"\n"),
                        border="white",cex=0.8,box.col = "white")


# plot prc curve 
plot(glmnetPlotPrc,col="red",lty=1, lwd=2,main="PRC curves")
plot(randomForestPlotPrc,col="blue",lty=1, lwd=2,add=T)
plot(xgboostPlotPrc,col="green",lty=1, lwd=2,add=T)
legend(0.2,0.2, legend = c("logisticRegression", "randomforest","gradientBoosting"), lty=1,
       lwd=2,col=c("red","blue" ,"green"),cex = 1,bty ="n")
legend(0.2,0.45,c(paste(c("AUPRC for logistic = ","AURPC for randomforest = ",
                          "AUPRC for gradientboosting = "),
                        c(round(glmnetAuprc, digits = 2),round(randomForestAuprc,digits=2),
                          round(xgboostAuprc,digits=2)),sep=""),"\n"),
       border="white",cex=0.8,box.col = "white")