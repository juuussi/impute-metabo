
simulated_data <- simulate_data(data=reference_data, nrow=100, ncol=600)
miss_data <- simulate_missingness(data=simulated_data, mcar=0.3)
colnames(miss_data) <- sapply(X=1:600, FUN=function(x) paste0('V', x))  

## Conigrave package it doesnt work
library(Conigrave)
correlatrix(miss_data,colnames(miss_data))


# for time data
library(cutoffR)
impdata <- cutoff(miss_data)
missdf <- as.data.frame(miss_data)

library(ForImp)
library(homals)
# SUMMARY OF MISSIGNESS PER VARIABLE
miss_variables <- missingness(miss_data)$missing_values_per_variable
 
library(GenForImp)# it doesnt work
xForImpMahala <- ForImp.Mahala(missdf)

mi <- ForImp.PCA(miss_data)

# computing the Relative Mean Square Error
error <- sum(apply((x0-xForImpPCA)^2/diag(var(x0)),2,sum)) / n
error
## hot.deck imputation
library(hot.deck) # it doesnt work
affinity(miss_data)
out <- hot.deck(miss_data)
out2 <- hd2amelia(out)
###  imp4p package

library(imp4p)
logmiss <- log2(miss_data)
estim.bound(logmiss,c)


#Simulating data
res.sim=sim.data(nb.pept=2000,nb.miss=600,pi.mcar=0.2,para=10,nb.cond=2,nb.repbio=3,
                 nb.sample=5,m.c=25,sd.c=2,sd.rb=0.5,sd.r=0.2);
data=res.sim$dat.obs;
cond=res.sim$conditions;
#Estimation of lower and upper bounds for each missing value
res=estim.bound(data,conditions=cond)

### 
library(imputeLCMD)
binVec  <- sample(rep(1,each =100))

output <- impute.MAR(miss_data, binVec, method = "KNN")
###


exprsDataObj = generate.ExpressionData(nSamples1 = 6, nSamples2 = 6,
                                       meanSamples = 0, sdSamples = 0.2,
                                       nFeatures = 1000, nFeaturesUp = 50, nFeaturesDown = 50,
                                       meanDynRange = 20, sdDynRange = 1,
                                       meanDiffAbund = 1, sdDiffAbund = 0.2)
exprsData = exprsDataObj[[1]]
# insert 15% missing data with 100% missing not at random
m.THR = quantile(exprsData, probs = 0.15)
sd.THR = 0.1
MNAR.rate = 100
exprsData.MD.obj = insertMVs(exprsData,m.THR,sd.THR,MNAR.rate)
exprsData.MD = exprsData.MD.obj[[2]]
# run model.Selector
m.s = model.Selector(exprsData.MD)
# perform MAR/MCAR imputation
exprsData.MAR.imputed = impute.MAR (exprsData.MD, m.s, method = "MLE")
## The function is currently defined as
function (dataSet.mvs, model.selector, method = "MLE")
{
  if (length(which(model.selector[[1]] == 1)) == 0) {
    dataSet.imputed = dataSet.mvs
  }
  else {
    dataSet.MCAR = dataSet.mvs[which(model.selector[[1]] ==
                                       1), ]
    switch(method, MLE = {
      dataSet.MCAR.imputed = impute.wrapper.MLE(dataSet.MCAR)
    }, SVD = {
      dataSet.MCAR.imputed = impute.wrapper.SVD(dataSet.MCAR,
                                                K = 2)
    }, KNN = {
      dataSet.MCAR.imputed = impute.wrapper.KNN(dataSet.MCAR,
                                                K = 15)
    })
    dataSet.imputed = dataSet.mvs
    dataSet.imputed[which(model.selector[[1]] == 1), ] = dataSet.MCAR.imputed
  }
  return(dataSet.imputed)
}

###

library(hydroGOF)

###### package mi
library(mi)
# STEP 0: GET DATA
data(nlsyV, package = "mi")
# STEP 0.5 CREATE A missing_variable (you never need to actually do this)
income <- missing_variable(nlsyV$income, type = "continuous")
show(income)
# STEP 1: CONVERT IT TO A missing_data.frame
mdf <- missing_data.frame(nlsyV) # this calls missing_variable() internally
show(mdf)

#####################
library(mi)

# STEP 1: Convert to a missing_data.frame
mdf <- missing_data.frame(as.data.frame(miss_data2[,1:100])) # this calls missing_variable() internally

V1 <- missing_variable(miss_data2[,1], type = "continuous")

show(mdf)
# STEP 2: change things
mdf <- change(mdf, y = "V81", what = "transformation", to = "identity")
# STEP 3: look deeper
summary(mdf)

# STEP 4: impute
## Not run:
imputations <- mi(mdf)## End(Not run)
mi::image(miss_data)

###################
library(imputeLCMD)
miss_data2 <- simulate_missingness(data=simulated_data, mnar=0.3)
exprsData.imputed2 = impute.wrapper.MLE(miss_data2)
exprsData.imputed2 =impute.wrapper.SVD(miss_data2, 15)
hist(exprsData.imputed2[,32])


###########
library(corrplot)
M <- cor(simulated_data)
corrplot(M, method="color",order = "hclust")

### 

# Label <- rbinom(150, 1, 0.5)
# 
# for (l in Label) {
#   if (l == 1) {
#     pic_col <- unique(sample(x=1:100, size=2))
#     print(pic_col)
#     te <- simulated_data[,pic_col]
#     mvalues <-  mean(te)
#     te[te>mvalues] <- NA
#     simulated_data[,pic_col] <- te
#   
# }
# }
#######################################################################
# MAR missingness 

simulated_data <- simulate_data(data=reference_data, nrow=15, ncol=30)
miss_data <- simulate_missingness(simulated_data,mcar =0.05)
initial_nas <- sum(colSums(is.na(simulated_data)))
current_nas <- initial_nas
print(current_nas)
simulated_data_n <- ncol(simulated_data) * nrow(simulated_data)

percentage <- 0.3


mar <- percentage
    
    current_nan_percentage <- 0
    columns <- NULL
    support <- 1:(ncol(simulated_data)-1)
    if  (((current_nas - initial_nas) / simulated_data_n) < mar){
      while(length(support) > 0 && current_nan_percentage < mar ) {
        # try to take unique column each time
        column <- sample(x=support, size=1)
        support <- setdiff(1:(ncol(simulated_data)-1), columns)
        
        columns <- c(columns, column)
        #print(columns)
        # condition: 
        # calculate the mean of the choosen column
        mn <- mean(simulated_data[,column])
        # dependency : remove the values on the next column above the mean 
        simulated_data[simulated_data[,column+1] > mn, column+1] <- NA 
        # check the current missing porpotion
        current_nan_percentage <- length(which(is.na(simulated_data)))/length(simulated_data)
        print(current_nan_percentage)
      }
    }


print(columns)


#############
library(Amelia)
library(mlbench)
# load dataset
data(Soybean)
# create a missing map
missmap(Soybean, col=c("black", "grey"), legend=FALSE)
#### 


library(MissMech)

order1 <- OrderMissing(miss_data, del.lesscases = 0)

out <- TestMCARNormality(miss_data)
########################################################
library(Amelia)
simulated_data <- simulate_data(reference_data,100,10)
miss_data <- simulate_missingness(simulated_data,mnar = 0.1)
a.out <- amelia(miss_data, m = 100, idvars=2)
View(a.out$imputations)
############################################
library(PEMM)
PEM.result3 = PEMM_fun(miss_data, phi=1)
nrmse_error <- missForest:: nrmse(PEM.result3$Xhat, miss_data, simulated_data)
############################################################################
## Load a sample metabolite dataset with 5\% missig values (metaboliteData)e
data(metaboliteData)
## Perform Bayesian PCA with 2 components
pc <- pca(miss_data, method="bpca", nPcs=10)
## Get the estimated principal axes (loadings)
loadings <- loadings(pc)
## Get the estimated scores
scores <- scores(pc)
## Get the estimated complete observations
cObs <- completeObs(pc)
nrmse_error <- missForest:: nrmse(cObs, miss_data, simulated_data)
