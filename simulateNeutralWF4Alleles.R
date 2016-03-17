#
#	Numerical solution of stationary distribution of the neutral Wright-Fisher model 
#		for 4-alleles for plotting Figures 4 and 5 of "An approximate stationary solution for 
#		multi-allele neutral diffusion when q_{ab} << 1".  
#
#	With alpha12 = 1/1000 etc and M = 30, takes 7.12 hours to run on a MacBook Air (1.7 GHz Intel Core i5)
#	With alpha12 = 1/100 etc and M = 30, takes 51 minutes to run
#
#	Version 1 11/03/16 Conrad Burden and Yurong Tang
#
	library(combinat)
#
#		The next line (and others like it further down) are to read 
# 			the current time so one can check 
# 			how long it takes to run when N is large
#  
	t1 <- Sys.time()    

#	Population size, alphabet size
#
	N <- 15
	K <- 4	# Currently only works for K=4; it will crash if you change this line
#
#	Parameterise Q-matrix as GTR + net flux
#
	alpha12 <- 1/100
	alpha13 <- 2/100
	alpha14 <- 3/100
	alpha23 <- 4/100
	alpha24 <- 5/100
	alpha34 <- 6/100
	Pi <- c(1, 2, 3, 4)
	Pi <- Pi/sum(Pi)
	phi <- c(0.4/100, 0.1/100, -0.03/100)
#
	QMatrix <- array(0, dim=c(4, 4))
	QMatrix[1, 2] <- alpha12*Pi[2] + phi[3]/Pi[1]/2
	QMatrix[1, 3] <- alpha13*Pi[3] - phi[2]/Pi[1]/2
	QMatrix[1, 4] <- alpha14*Pi[4] + phi[2]/Pi[1]/2 - phi[3]/Pi[1]/2
	QMatrix[2, 1] <- alpha12*Pi[1] - phi[3]/Pi[2]/2
	QMatrix[2, 3] <- alpha23*Pi[3] + phi[1]/Pi[2]/2
	QMatrix[2, 4] <- alpha24*Pi[4] - phi[1]/Pi[2]/2 + phi[3]/Pi[2]/2
	QMatrix[3, 1] <- alpha13*Pi[1] + phi[2]/Pi[3]/2
	QMatrix[3, 2] <- alpha23*Pi[2] - phi[1]/Pi[3]/2
	QMatrix[3, 4] <- alpha34*Pi[4] + phi[1]/Pi[3]/2 - phi[2]/Pi[3]/2
	QMatrix[4, 1] <- alpha14*Pi[1] - phi[2]/Pi[4]/2 + phi[3]/Pi[4]/2
	QMatrix[4, 2] <- alpha24*Pi[2] + phi[1]/Pi[4]/2 - phi[3]/Pi[4]/2
	QMatrix[4, 3] <- alpha34*Pi[3] - phi[1]/Pi[4]/2 + phi[2]/Pi[4]/2
	diag(QMatrix) <- -rowSums(QMatrix)
#
#	Set up indices of P matrix 
#
	iVectors <- t(xsimplex(K, N)) # each row of this array is a vector of the form (i1, i2, i3, i4)
	Pdim <- dim(iVectors)[1]
#
#	Construct Pmatrix
#
	uMatrix <- diag(1,K,K) + QMatrix/N      							# Eq.(3)
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
#	Calculation of theoretical distribution along edges of tetrahedron starts here
#
#	Eq.(43):
#
	CMatrix <- array(0, dim=c(4, 4))
	CMatrix[1, 2] <- CMatrix[2, 1] <- 2*alpha12*Pi[1]*Pi[2]
	CMatrix[1, 3] <- CMatrix[3, 1] <- 2*alpha13*Pi[1]*Pi[3]
	CMatrix[1, 4] <- CMatrix[4, 1] <- 2*alpha14*Pi[1]*Pi[4]
	CMatrix[2, 3] <- CMatrix[3, 2] <- 2*alpha23*Pi[2]*Pi[3]
	CMatrix[2, 4] <- CMatrix[4, 2] <- 2*alpha24*Pi[2]*Pi[4]
	CMatrix[3, 4] <- CMatrix[4, 3] <- 2*alpha34*Pi[3]*Pi[4]
#
#	Eq.(39)
#
	PhiMatrix <- array(dim=c(3, 4))
	PhiMatrix[1, 2] <- phi[3]
	PhiMatrix[1, 3] <- - phi[2]
	PhiMatrix[1, 4] <- phi[2] - phi[3]
	PhiMatrix[2, 3] <- phi[1]
	PhiMatrix[2, 4] <- phi[3] - phi[1]
	PhiMatrix[3, 4] <- phi[1] - phi[2]
#
#	Corner points: integral from 0 to 1/N of 2-allele projections of f(x), approx by Eq.(44)
#
	ProbCorner <- Pi - 2*rowSums(CMatrix)*log(N)
#
#	Labels for plots of 6 edges
#
	label <- array(dim=c(3, 4)); yLabel <- array(dim=c(3, 4))
	label[1, 2] <- "A1-A2 edge"; yLabel[1, 2] <- "12"
	label[1, 3] <- "A1-A3 edge"; yLabel[1, 3] <- "13"
	label[1, 4] <- "A1-A4 edge"; yLabel[1, 4] <- "14"
	label[2, 3] <- "A2-A3 edge"; yLabel[2, 3] <- "23"
	label[2, 4] <- "A2-A4 edge"; yLabel[2, 4] <- "24"
	label[3, 4] <- "A3-A4 edge"; yLabel[3, 4] <- "34"
#
#############
#
#	Plot distribution on all 6 edges
#
	x <- (1:499)/500
#
	quartz()
	oldpar <- par(mfrow=c(2,3))
	for(firstBase in 1:3){
		for(secondBase in (firstBase + 1):4){
			Stat <- array(dim=(N + 1))
			for (i in 0:N){
				index <- which(iVectors[,firstBase]==(N - i) & iVectors[,secondBase]==i)
				Stat[i + 1] <- SS[index]
				}
			plot((N:0)/N,log10(N*Stat), 
				pch=4, 
				col="blue", 	# want LH end of plot to be x=0, i.e. for f_12, all A2 at x=0
				xlab = expression(x), 
				ylab = substitute(log[10]*f[subscript](x), list(subscript=yLabel[firstBase, secondBase])), 
				main=label[firstBase, secondBase])
			fab0fx <- CMatrix[firstBase, secondBase]/x/(1 - x) - 
							PhiMatrix[firstBase, secondBase]*(1/x - 1/(1 - x))
			points(x, log10(fab0fx), type = "l", col = "red")
			points(c(1, 0), log10(N*c(ProbCorner[firstBase], ProbCorner[secondBase])), pch=1, col = "red")
			}
		}
	par(oldpar)
#