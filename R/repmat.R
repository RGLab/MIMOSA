##' R equivalent of repmat
##'
##'@details replicates mx rows x m and nx columns by n
##'@param X input matrix
##'@param m number times to replicate the rows
##'@param n number of times to replicate the columns
repmat = function(X,m,n){
if (is.matrix(X)) {
mx = dim(X)[1]
nx = dim(X)[2]
} else {
   mx = 1
   nx = length(X)
}
matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}
