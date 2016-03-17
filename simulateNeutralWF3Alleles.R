#
#	Numerical solution of stationary distribution of the neutral Wright-Fisher model 
#		for 3-alleles for plotting Figure 2 of "An approximate stationary solution for 
#		multi-allele neutral diffusion when q_{ab} << 1".  
#
#	(takes about 6 minutes for simulation to converge if 
#		the elements of QMatrix are O(1/10^3) and N = 30.  
#		on a MacBook air (1.7 GHz Intel Core i5)
#	 takes about 25 secs for simulation to converge if 
#		the elements of QMatrix are O(1/10^3) and N = 15.  
#		Analytic solution and plots happen instantly)
#
#	Version 1 11/03/16 Conrad Burden and Yurong Tang
#
	library(combinat)
#
#		The next line (and others like it further down) are to read 
# 			the current time so one can check 
# 			how long it takes to run when M is large
#  
	t1 <- Sys.time()    
#
#	Population size, alphabet size
#
	N <- 15
	K <- 3	# Only works for K=3; it will crash if you change this line
#
#	Parameterise Q-matrix as GTR + net flux
#
	alpha12 <- 1/1000
	alpha13 <- 2/1000
	alpha23 <- 5/1000
	Pi <- c(5, 3, 2)
	Pi <- Pi/sum(Pi)
	phi <- 0.2/1000
#
	QMatrix <- array(0, dim=c(3, 3))
	QMatrix[1, 2] <- alpha12*Pi[2] + phi/Pi[1]/2
	QMatrix[1, 3] <- alpha13*Pi[3] - phi/Pi[1]/2
	QMatrix[2, 1] <- alpha12*Pi[1] - phi/Pi[2]/2
	QMatrix[2, 3] <- alpha23*Pi[3] + phi/Pi[2]/2
	QMatrix[3, 1] <- alpha13*Pi[1] + phi/Pi[3]/2
	QMatrix[3, 2] <- alpha23*Pi[2] - phi/Pi[3]/2
	diag(QMatrix) <- -rowSums(QMatrix)
#
#	Set up indices of P matrix defining the neutral Wright-Fisher model
#
	iVectors <- t(xsimplex(K, N)) # each row of this array is a vector of the form (i1, i2, i3)
	Pdim <- dim(iVectors)[1]
#
#	Construct Pmatrix
#
	uMatrix <- diag(1,K,K) + QMatrix/N      						# Eq.(3)
	Pmatrix <- as.matrix(array(dim=c(Pdim, Pdim)))
	for(i in 1:Pdim){
		iThisRow <- iVectors[i,]
		psi_ia <- rowSums(t(iThisRow*uMatrix)/N)  					# Eq.(2)
		combFactors <- exp(logfact(N) - rowSums(logfact(iVectors))) # the logs are to avoid roundoff errors
		Pmatrix[i,] <- combFactors * exp(iVectors%*%log(psi_ia))  	# Eq.(1)
		}
#
	cat("\n Pmatrix constructed \n")
	t1.5 <- Sys.time()
#
#	Find stationary distribution:
#		Start with a uniform distribution and multiply by P until it converges
#
	vec <- rep(1,Pdim)/Pdim
	diff <- 1
	iter <- 0
	while (abs(diff) > 1.e-11) {
		iter <- iter + 1
		oldVec <- vec
		vec <- vec%*%Pmatrix
		diff <- max(abs(vec - oldVec))
		cat("\n iter =", iter, "diff =", diff)
		}
	SS <- c(vec)
#
#	Print out time taken to convergence
#
	t2 <- Sys.time()
	deltaT <- t2 - t1
	cat("\n \n time to Convergence:")
	print(deltaT)
	cat("\n")
#
#############
#
#	Plots
#
	if(.Platform$OS.type=="unix") {
		quartz(height=7.5, width=7.5)
		}
	if(.Platform$OS.type=="windows") {
		windows(height=7.5, width=7.5)
		}
	oldpar <- par(mfrow=c(2,2))
#
#############
#
#	2d plot
#
	xCorner <- c(-sqrt(3)/2, sqrt(3)/2, 0)
	yCorner <- c(0.5, 0.5,-1)
	xCorner <- c(0, sqrt(3)/2, -sqrt(3)/2)
	yCorner <- c(1, -0.5, -0.5)
#
	xx <- iVectors%*%xCorner
	yy <- iVectors%*%yCorner
	plot(xx, yy, pch=16, cex=4*SS^(1/3), 
			xlim=1.2*N*c(-sqrt(3)/2, sqrt(3)/2), ylim=1.2*N*c(-0.5, 1), 
			bty="n", col.axis="white", ann=F, fg="white") 
	text(1.2*N*xCorner, 1.2*N*yCorner, c(expression(A[1]), expression(A[2]), expression(A[3])))
#
#############
#
#	normalisations for continuum solutions along edges	Eq.(27)
#
	C12 <- 2*alpha12*Pi[1]*Pi[2]
	C23 <- 2*alpha23*Pi[2]*Pi[3] 
	C31 <- 2*alpha13*Pi[1]*Pi[3]
#
#	Line densities Eq.(32)
#
	f_ab <- function(x, C, phii){
		f_ab <- C/x/(1 - x) - phii*(1/x - 1/(1 - x))
		return(f_ab)
		}   
#
#	Corner points: integral form 0 to 1/N of 2-allele projections of f(x)	Eq.(28) & (29)
#
	ProbAllA1 <- Pi[1]*N^(-2*(alpha12*Pi[2] + alpha13*Pi[3]))
	ProbAllA2 <- Pi[2]*N^(-2*(alpha12*Pi[1] + alpha23*Pi[3]))
	ProbAllA3 <- Pi[3]*N^(-2*(alpha13*Pi[1] + alpha23*Pi[2]))
#
#	range of x for plots
#
	x <- (1:499)/500
#
##############
#
#
#	Plot along 1-2 edge
#
	indices12 <- which(iVectors[,3]==0)
	xDiscrete <- iVectors[indices12,1]/N
	fDiscrete <- SS[indices12]*N
#
	plot(xDiscrete, log10(fDiscrete), pch=4, col="blue", 
		main="1-2 edge", xlab = "proportion of type-1 alleles", 
		ylab = expression(log[10]*f[12](x)))
	points(x, log10(f_ab(x, C12, phi)), type = "l", col = "red")
	points(c(0, 1), log10(N*c(ProbAllA2, ProbAllA1)), pch=3, col = "red")
#
#############
#
#	Plot along 2-3 edge
#
	indices23 <- which(iVectors[,1]==0)
	xDiscrete <- iVectors[indices23,2]/N
	fDiscrete <- SS[indices23]*N
#
	plot(xDiscrete, log10(fDiscrete), pch=4, col="blue", 
		main="2-3 edge", xlab = "proportion of type-2 alleles", ylab = expression(log[10]*f[23](x)))
	points(x, log10(f_ab(x, C23, phi)), type = "l", col = "red")
	points(c(0, 1), log10(N*c(ProbAllA3, ProbAllA2)), pch=3, col = "red")
#
#############
#
#	Plot along 3-1 edge
#
	indices31 <- which(iVectors[,2]==0)
	xDiscrete <- iVectors[indices31,3]/N
	fDiscrete <- SS[indices31]*N
#
	plot(xDiscrete, log10(fDiscrete), pch=4, col="blue", 
		main="3-1 edge", xlab = "proportion of type-3 alleles", ylab = expression(log[10]*f[31](x)))
	points(x, log10(f_ab(x, C31, phi)), type = "l", col = "red")
	points(c(0, 1), log10(N*c(ProbAllA1, ProbAllA3)), pch=3, col = "red")
#
##############
#
	par(oldpar)
#
##############
#
#	Check how close corner probabilitites are to theory
#
	SAtCorner1 <- SS[which(iVectors[,2] + iVectors[,3]==0)]
	SAtCorner2 <- SS[which(iVectors[,3] + iVectors[,1]==0)]
	SAtCorner3 <- SS[which(iVectors[,1] + iVectors[,2]==0)]
#
	cat("\n     SAtCorner1 =", SAtCorner1, "    	ProbAllA1 =", ProbAllA1, 
		"\n     SAtCorner2 =", SAtCorner2, "	    ProbAllA2 =", ProbAllA2,
		"\n     SAtCorner3 =", SAtCorner3, "	    ProbAllA3 =", ProbAllA3, "\n\n")
#
##############
	
	