#' likelihoodSurfacePlot
#' @param bmr The BetaMixResult object
#' @param grid.x The upper limit of the x grid to calculate the likelihood surface
#' @param grid.y The upper limit of the y grid to calculate the likelihood surface
#' @param whch The indxex of the observation for which to plot the likelihood surface
#' @param nlevels The number of contour levels to plot
#' @param fdr The false discovery rate threshold used to color significant points
#' @returnType NULL 
#' @return Silently returns NULL
#' @author finak
#' @export

setGeneric("likelihoodSurfacePlot",
		function(bmr, grid.x, grid.y, whch, nlevels, fdr, ...)
		{
			standardGeneric("likelihoodSurfacePlot")
		})

setMethod("likelihoodSurfacePlot",
		signature(bmr = "BetaMixResult", grid.x = "integer", grid.y = "integer", whch = "integer", nlevels = "integer", fdr = "numeric"),
		function(bmr, grid.x, grid.y, whch, nlevels = 30, fdr = 0.01)
		{
#set up the grid
			if(bmr@alternative.model=="not equal"){
				alternative<-"two.sided"
			}else{
				alternative<-"greater"
			}
			
			gr<-expand.grid(1:grid.x,1:grid.y)
			N<-colSums(matrix(unlist(bmr@data[whch,]),2,byrow=T))
			gr<-cbind(gr,t(N-t(gr)))
			colnames(gr)<-c("ns","nu","Ns","Nu")
			
			z<-sqrt(((matrix(bmr@w[1]*exp(.Call("MarginalNULL",ns=gr[,"ns"],Ns=gr[,"Ns"],nu=gr[,"nu"],Nu=gr[,"Nu"],alpha0=rep(bmr@alpha0,nrow(gr)),beta0=rep(bmr@beta0,nrow(gr)),alphaS=rep(bmr@alphaS,nrow(gr)),betaS=rep(bmr@betaS,nrow(gr)),bmr@w,log=TRUE,package="MIMOSA"))+bmr@w[2]*exp(.Call("MarginalGT",ns=gr[,"ns"],Ns=gr[,"Ns"],nu=gr[,"nu"],Nu=gr[,"Nu"],alpha0=rep(bmr@alpha0,nrow(gr)),beta0=rep(bmr@beta0,nrow(gr)),alphaS=rep(bmr@alphaS,nrow(gr)),betaS=rep(bmr@betaS,nrow(gr)),bmr@w,log=TRUE,500,package="MIMOSA")),ncol=grid.x))))
			
#plot the likelihood surface
			image(z,x=1:grid.x,y=1:grid.y,col=rev(gray.colors(30)),xlab="ns",ylab="nu")
			contour(1:grid.x,1:grid.y,z,xlab="ns",ylab="nu",add=TRUE,nlevels=nlevels)
#plot the data points
			points(bmr@data[whch,1:2,FALSE],pch=3,col=c("black","red")[(bmr@fdr<fdr)+1][whch])    
			points(bmr@data[whch,1:2,FALSE],pch=2,col=c("black","green")[(apply(bmr@data,1,function(x)p.adjust(fisher.test(matrix(x,2,byrow=T),alternative=alternative)$p.value,"fdr")<fdr)+1)[whch]])    
			
#plot the identity line
			lines(t(t(matrix(seq(0,1,l=2))%*%N)),lwd=2,lty=2)   
			invisible(NULL)
		})
