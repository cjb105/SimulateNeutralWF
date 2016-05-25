#
#	Numerical solution of stationary distribution of Wright-Fisher model 
#		for 3-alleles (takes about 6 minutes for simulation to converge if 
#		the elements of QMatrix are O(1/10^3) and M = 30.  
#		Analytic solution and plots happen instantly)
#
#	Version 1 06/11/15
#	Version 2 10/11/15 plotted distribution along edges
#	Version 3 11/11/15 added plots of marginal distributions
#	Version 4 13/11/15 added plot of finite population equivalent 2-allele model
#	Version 5 02/12/15 reparameterise Q matrix as reversible part + flux
#	Version 6 07/12/15 1st attempt at superimposing continuum solution along 2-3 edge (not correct)
#	Version 7 20/01/16 continuum solutions matched up with (normalised-by-hand) discrete solution
#	Version 8 21/01/16 failed attempt at normalisation
#	Version 9 21/01/16 normalisation sorted for reversible QMatrix
#	Version 10 22/10/16 first atempt at including 'flux' part - not working yet
#	Version 11 22/01/16 use function Ibeta for incomplete beta function - still not working
#	Version 12 25/01/16 normalisation fixed, trick is to write soln in terms of 1/x & 1/(1 - x)
#	Version 13 25/01/16 tidied up plots
#	Version 14 28/01/16 further cosmetic changes to final plot
#	Version 15 02/02/16 corner (i.e. non-segregating) probabilities included
#	Version 16 22/02/16 replaced normalisations with 1st order approximations
#	Version 18 25/05/16 same y-axis scale in all 3 plots
#
	library(combinat)
#
# 		Change the next line to the directory this file is in
#
	setwd("~/Documents/Population genetics/Paper 2") 
#
#		The next line (and others like it further down) are only to read 
# 			the current time so I can check 
# 			how long it takes to run when M is large
#  
	t1 <- Sys.time()    
#
#	Population size, alphabet size
#
	M <- 30
	d <- 3	# Currently only works for d=3; it will crash if you change this line
#
#	Parameterise Q-matrix as GTR + net flux
#
	alpha <- 1/1000
	beta <- 2/1000
	gamma <- 5/1000
	Pi <- c(5, 3, 2)
	Pi <- Pi/sum(Pi)
	phi <- 0.2/1000
#
	QMatrix <- array(0, dim=c(3, 3))
	QMatrix[1, 2] <- alpha*Pi[2] + phi/Pi[1]/2
	QMatrix[1, 3] <- beta*Pi[3] - phi/Pi[1]/2
	QMatrix[2, 1] <- alpha*Pi[1] - phi/Pi[2]/2
	QMatrix[2, 3] <- gamma*Pi[3] + phi/Pi[2]/2
	QMatrix[3, 1] <- beta*Pi[1] + phi/Pi[3]/2
	QMatrix[3, 2] <- gamma*Pi[2] - phi/Pi[3]/2
	diag(QMatrix) <- -rowSums(QMatrix)
#
#	Set up indices of P matrix 
#
	iVectors <- t(xsimplex(d, M)) # each row of this array is a vector of the form (i1, i2, i3)
	Pdim <- dim(iVectors)[1]
#
#	Construct Pmatrix
#
	uMatrix <- diag(1,d,d) + QMatrix/M      # Eq.(32) of 9/6/15 handwritten notes
	Pmatrix <- as.matrix(array(dim=c(Pdim, Pdim)))
	for(i in 1:Pdim){
		iThisRow <- iVectors[i,]
		psi_ia <- rowSums(t(iThisRow*uMatrix)/M)  # Eq.(12) of 9/6/15 handwritten notes
		combFactors <- exp(logfact(M) - rowSums(logfact(iVectors))) # the logs are to avoid roundoff errors
		Pmatrix[i,] <- combFactors * exp(iVectors%*%log(psi_ia))  # Eq.(13) of 9/6/15 handwritten notes
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
#	copy stationary distribution into a convenient array 
#
	Stationary1 <- array(dim=(M + 1)^(d - 1))
	for(iRow in 1:Pdim){
		posnInArrayMinus1 <- iVectors[iRow, 1:(d - 1)]
		Stationary1[sum(posnInArrayMinus1*(M + 1)^(0:(d - 2))) + 1] <- SS[iRow]
		}
	Stationary <- array(Stationary1, dim=rep(M + 1, d - 1))
#
#	The array element Stationary[i1 + 1, i2 + 1] is the stationary distribution 
#		probability at the point (i1, i2, i3) where i3 = M - i1 - i2Stationary[M + 1, 1]
#
	t2 <- Sys.time()
	deltaT <- t2 - t1
#
#############
#
#	Plots
#
	quartz(height=7.5, width=7.5)
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
			xlim=1.2*M*c(-sqrt(3)/2, sqrt(3)/2), ylim=1.2*M*c(-0.5, 1), 
			bty="n", col.axis="white", ann=F, fg="white") 
	text(1.2*M*xCorner, 1.2*M*yCorner, c(expression(A[1]), expression(A[2]), expression(A[3])))
#
#############
#
#	normalisations for continuum solutions along edges
#
	A12 <- 2*alpha*Pi[1]*Pi[2]
	A23 <- 2*gamma*Pi[2]*Pi[3] 
	A31 <- 2*beta*Pi[1]*Pi[3]  
#
#	Corner points: integral form 0 to 1/M of 2-allele projections of f(x)
#
	ProbAllA1 <- Pi[1]*M^(-2*(alpha*Pi[2] + beta*Pi[3]))
	ProbAllA2 <- Pi[2]*M^(-2*(alpha*Pi[1] + gamma*Pi[3]))
	ProbAllA3 <- Pi[3]*M^(-2*(beta*Pi[1] + gamma*Pi[2]))
#
##############
#
#
#	Plot along 1-2 edge
#
#
#	Infinite population solution
#
	x <- c((1:9)/10000, (1:9)/1000, (1:99)/100, (991:999)/1000, (9991:9999)/10000)
	x <- (1:499)/500
	q12 <- QMatrix[1, 2]
	q21 <- QMatrix[2, 1]
	f3OfX <- A12/x/(1 - x) - phi*(1/x - 1/(1 - x))
#
	indices12 <- which(iVectors[,3]==0)
	xDiscrete <- iVectors[indices12,1]/M
	fDiscrete <- SS[indices12]*M
#
	plot(xDiscrete, log10(fDiscrete), pch=4, col="blue", 
		main="1-2 edge", xlab = "proportion of type-1 alleles", 
		ylab = expression(log[10]*f[12](x)), ylim=c(-3, 1.3))
	points(x, log10(f3OfX), type = "l", col = "red")
	points(c(0, 1), log10(M*c(ProbAllA2, ProbAllA1)), pch=3, col = "red")
#
#############
#
#	Plot along 2-3 edge
#
#
#	Infinite population solution
#
	q23 <- QMatrix[2, 3]
	q32 <- QMatrix[3, 2]
	f1OfX <- A23/x/(1 - x) - phi*(1/x - 1/(1 - x))
#
	indices23 <- which(iVectors[,1]==0)
	xDiscrete <- iVectors[indices23,2]/M
	fDiscrete <- SS[indices23]*M
#
	plot(xDiscrete, log10(fDiscrete), pch=4, col="blue", 
		main="2-3 edge", xlab = "proportion of type-2 alleles", 
		ylab = expression(log[10]*f[23](x)), ylim=c(-3, 1.3))
	points(x, log10(f1OfX), type = "l", col = "red")
	points(c(0, 1), log10(M*c(ProbAllA3, ProbAllA2)), pch=3, col = "red")
#
#############
#
#	Plot along 3-1 edge
#
#
#	Infinite population solution
#
	q31 <- QMatrix[3, 1]
	q13 <- QMatrix[1, 3]
	f2OfX <- A31/x/(1 - x) - phi*(1/x - 1/(1 - x))
#
	indices31 <- which(iVectors[,2]==0)
	xDiscrete <- iVectors[indices31,3]/M
	fDiscrete <- SS[indices31]*M
#
	plot(xDiscrete, log10(fDiscrete), pch=4, col="blue", 
		main="3-1 edge", xlab = "proportion of type-3 alleles", 
		ylab = expression(log[10]*f[31](x)), ylim=c(-3, 1.3))
	points(x, log10(f2OfX), type = "l", col = "red")
	points(c(0, 1), log10(M*c(ProbAllA1, ProbAllA3)), pch=3, col = "red")
#
##############
#
	par(oldpar)
#
##############
#
#	Check how close corner probabilitites are to theory
#
	S1 <- Stationary[M + 1, 1]
	S2 <- Stationary[1, M + 1]
	S3 <- Stationary[1, 1]
#
	cat("\n     S1 =", S1, "    ProbAllA1 =", ProbAllA1, 
		"\n     S2 =", S2, "	    ProbAllA2 =", ProbAllA2,
		"\n     S3 =", S3, "	    ProbAllA3 =", ProbAllA3, "\n\n")
#
#	Check how well Continuum limit is normalised
#
#	Corners
#
	f1 <- function(x){return(A23 * x^(2*q32 - 1) * (1 - x)^(2*q23 - 1))}
	f2 <- function(x){return(A31 * x^(2*q13 - 1) * (1 - x)^(2*q31 - 1))}
	f3 <- function(x){return(A12 * x^(2*q21 - 1) * (1 - x)^(2*q12 - 1))}
#
	I1 <- integrate(f1, lower=1/M, upper=1 - 1/M)$value
	I2 <- integrate(f2, lower=1/M, upper=1 - 1/M)$value
	I3 <- integrate(f3, lower=1/M, upper=1 - 1/M)$value
#
	cat("\n ProbAllA1 =", ProbAllA1, "		ProbAllA2 =", ProbAllA2, "	ProbAllA3 =", ProbAllA3)
	cat("\n I1 =", I1, "		I2 =", I2, "	I3 =", I3)
	cat("\n    Total Corners =", Stot <- ProbAllA1 + ProbAllA2 + ProbAllA3)
	cat("\n  Total Integrals =", Itot <- I1 + I2 + I3)
	cat("\n Total Everything =", Stot + Itot, "\n")
#
##############
	
	