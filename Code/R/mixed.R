################## Function for Profile Analysis ##################
#	Code to CLMM Profile Analysis
#		using data from MATLAB and flat file communication
#
#	Created:	01/14/2020
#	Modified:	01/14/2020
#
#	By:			MJ Meyer
#####################################################################

library(ordinal)
library(R.matlab)
library(tidyr)

setwd('/Users/mjm556/Documents/MATLAB')
Ylist	    <- readMat('Y.mat')
Xlist	    <- readMat('X.mat')

Y         <- Ylist$Y
X         <- Xlist$x1

df        <- data.frame(Y = Y, X = X, ID = 1:nrow(Y))
dfl       <- gather(df, key = day, value = count, Y.1:Y.365)
dfls      <- dfl[order(dfl$ID),]
dfls$t    <- rep(1:ncol(Y), nrow(Y))

dfls$resp   <- ordered(as.factor(dfls$count))
dfls$stime  <- dfls$t/max(dfls$t)

stime     <- proc.time()
modelit   <- clmm(resp ~ X + stime + (1 + stime | ID), link = 'probit', data = dfls)
etime     <- proc.time() - stime

Xw        <- spread(dfls[,c('X', 'ID', 't')], key = t, value = X)
tw        <- spread(dfls[,c('stime', 'ID', 't')], key = t, value = stime)

eta       <- as.matrix(Xw[,-1]*modelit$beta['X'] + tw[,-1]*modelit$beta['stime'])
rtime     <- as.numeric(etime)[3]

writeMat('eta.mat', eta = eta)
writeMat('rtime.mat', rtime = rtime)
