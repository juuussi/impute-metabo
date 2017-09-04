### R code from vignette source 'pcaMethods.Rnw'

###################################################
### code chunk number 1: pcaMethods.Rnw:102-114
###################################################
library(pcaMethods)
x <- c(-4,7); y <- c(-3,4)
distX <- rnorm(100, sd=0.3)*3
distY <- rnorm(100, sd=0.3) + distX * 0.3
mat <- cbind(distX, distY)
res <- pca(mat, nPcs=2, method="svd", center=F)
loading <- loadings(res)[1,]
grad <- loading[2] / loading[1]
if (grad < 0)
  grad <- grad * -1
lx <- c(-4,7)
ly <- c(grad * -4, grad * 7)


###################################################
### code chunk number 2: pcaMethods.Rnw:118-125
###################################################
par(mar=c(2, 3, 2, 2))
plot(x,y, type="n", xlab="", ylab="")
abline(v=0, col="dark gray", lwd = 2); abline(h=0, col = "dark gray", lwd = 2)
points(distX, distY, type = 'p', col = "blue")
lines(lx,ly, lwd = 2)
points(-1, -1 * grad + 0.5, pch = 19, col = "red", lwd=4)
points(6, 6 * grad + 0.5, pch = 19, col = "red", lwd=4)


###################################################
### code chunk number 3: pcaMethods.Rnw:253-255
###################################################
library(lattice)
library(pcaMethods)


###################################################
### code chunk number 4: pcaMethods.Rnw:258-261
###################################################
library(pcaMethods)
data(metaboliteData)
data(metaboliteDataComplete)


###################################################
### code chunk number 5: pcaMethods.Rnw:264-266
###################################################
md  <- prep(metaboliteData, scale="none", center=TRUE)
mdC  <- prep(metaboliteDataComplete, scale="none", center=TRUE)


###################################################
### code chunk number 6: pcaMethods.Rnw:271-277
###################################################
resPCA  <- pca(mdC, method="svd", center=FALSE, nPcs=5)
resPPCA  <- pca(md, method="ppca", center=FALSE, nPcs=5)
resBPCA  <- pca(md, method="bpca", center=FALSE, nPcs=5)
resSVDI  <- pca(md, method="svdImpute", center=FALSE, nPcs=5)
resNipals  <- pca(md, method="nipals", center=FALSE, nPcs=5)
resNLPCA <- pca(md, method="nlpca", center=FALSE, nPcs=5, maxSteps=300)


###################################################
### code chunk number 7: pcaMethods.Rnw:293-296
###################################################
sDevs <- cbind(sDev(resPCA), sDev(resPPCA), sDev(resBPCA), sDev(resSVDI), sDev(resNipals), sDev(resNLPCA))
matplot(sDevs, type = 'l', xlab="Eigenvalues", ylab="Standard deviation of PC", lwd=3)
legend(x="topright", legend=c("PCA", "PPCA", "BPCA", "SVDimpute","Nipals PCA","NLPCA"), lty=1:6, col=1:6, lwd=3)

###################################################
### code chunk number 8: pcaMethods.Rnw:308-311
###################################################
par(mfrow=c(1,2))
plot(loadings(resBPCA)[,1], loadings(resPCA)[,1], xlab="BPCA", ylab="classic PCA", main = "Loading 1")
plot(loadings(resBPCA)[,2], loadings(resPCA)[,2], xlab="BPCA", ylab="classic PCA", main = "Loading 2")
#############################
## LLSimpute

lls <- llsImpute(miss_data, k = 30)
lls <- llsImpute(miss_data, k = 10, center = FALSE, completeObs = TRUE, correlation = "pearson", allVariables = TRUE, maxSteps = 100)
####

data1 <- simulate_data(reference_data,10,20)
miss_data <- simulate_missingness(data1,mcar = 0.5)
result <- pca(miss_data, method="svdImpute", nPcs=3, center = TRUE)
## Get the estimated complete observations
cObs <- completeObs(result)
##
esti <- kEstimateFast(miss_data, method = "bpca", evalPcs = 2:4, nruncv=1, em="nrmsep")
## KNNN
knn_esti <- impute.knn(as.matrix(miss_data))
knn_esti$data
