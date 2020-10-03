################## Function for O'Sullivan Splines ##################
#	Code to generate O'Sullivan Spline Penalty Matrix
#		using data from MATLAB and flat file communication
#
#	Created:	10/03/2018
#	Modified:	10/15/2018
#
#	By:			MJ Meyer
#####################################################################

library(R.matlab)
library(splines)

setwd('/Users/mjm556/Documents/Code')
source('formOmega.R')

setwd('/Users/mjm556/Documents/MATLAB')
Dlist	<- readMat('D.mat')
Klist	<- readMat('K.mat')

D 		<- as.numeric(Dlist$D)
K		<- as.numeric(Klist$K)

x				<- 1:D
numIntKnots		<- K

intKnots		<- quantile(unique(x), seq(0,1,length = numIntKnots+2))[-c(1,(numIntKnots+2))]
names(intKnots)	<- NULL

a				<- 0
b				<- D+1

Theta			<- bs(x, knots = intKnots, degree = 3, Boundary.knots = c(a, b), intercept = TRUE)
Omega			<- formOmega(a, b, intKnots)

writeMat('Theta.mat', Theta = Theta)
writeMat('Omega.mat', Omega = Omega)
