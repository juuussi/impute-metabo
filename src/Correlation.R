library(Hmisc)

## Checking correlation between your variables

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
miss_data <- data.frame(miss_data)

res2 <-rcorr(as.matrix(miss_data))
flattenCorrMatrix(res2$r, res2$P)
symnum(res2$r, abbr.colnames = FALSE)