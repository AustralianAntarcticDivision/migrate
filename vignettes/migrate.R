## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(migrate)

## ------------------------------------------------------------------------
yrs <- 2007:2015
relT <- rep( 1000, length(yrs))
# Annual Tag loss rate (for single-tag model) 
tagloss	<- rep(0.006, length(yrs))	
# Natural Mortality
M <- 0.155
# Release mortality
relM <- 0.1
# Proportion of fish that migrate
true_move	<- 0.01		
# Catch
C1 <- rep(10000, length(yrs))	
C2 <- rep(10000, length(yrs))
# Abundance in the two areas
B1 <- rep(1000000, length(yrs))	
B2 <- rep(1000000, length(yrs))

## ------------------------------------------------------------------------
# T: Number of tags available in Areas 1 and 2 (after half of M)
T1 <- T2	<- matrix(0,nrow=length(yrs),ncol=length(yrs),dimnames=list(yrs,yrs))
# Number of tags being recaptured
R1 <- R2 <- T1
# Number of tags moving
T_move <- T1

## ------------------------------------------------------------------------
for (x in 1:length(yrs)) {
  # Release tags
	T1[x,x]	<- relT[x]		
	# Apply release mortality
	T1[x,x]	<- (1-relM) * T1[x,x]		
	for (y in x:length(yrs)) {			
			# Available Tags after applying 0.5 M
			T1[x,y] <- T1[x,y] * exp(-0.5*M)
			T2[x,y] <- T2[x,y] * exp(-0.5*M)			 
			# Tagloss rate
			T1[x,y] <- T1[x,y] * exp(-tagloss[y])				 
			T2[x,y] <- T2[x,y] * exp(-tagloss[y])
			# Recaptures (by Petersen)
			R1[x,y] <- rbinom(n=1,size=round(T1[x,y],0),prob=C1[y]/B1[y])	
			R2[x,y]	<- rbinom(n=1,size=round(T2[x,y],0),prob=C2[y]/B2[y])	
			if((y+1) <= length(yrs)) {				# Prevent error
				T_move[x,y+1]	<- (T1[x,y] - R1[x,y])*true_move			# Movement after removals [= T*move*(1-move)**(y-x-1)]
				T1[x,y+1]	<- T1[x,y] - R1[x,y] - T_move[x,y+1]	# = T*(1-move)**(y-x)
				T2[x,y+1]	<- T2[x,y] - R2[x,y] + T_move[x,y+1]	# = 1-T*(1-move)**(y-x)
			}
			T1[x,y]		<- T1[x,y] * exp(-0.5*M)	# Apply 0.5 M - do not apply so T1&T2 remains mid-season numbers
			T2[x,y]		<- T2[x,y] * exp(-0.5*M)	# 
	}
}

## ------------------------------------------------------------------------
T1
T2

## ------------------------------------------------------------------------
R1
R2

## ------------------------------------------------------------------------
obs <- list(yrs=yrs,relT=relT,R1=R1,R2=R2,prec1=C1/B1,
            prec2=C2/B2,M=M,relM=relM,tagloss=tagloss)
para <- logitT(c(move=0.011694))	# move

## ------------------------------------------------------------------------
fit_ls <- move_sim(obs,para,rfit="R2",lltype="prop")

## ------------------------------------------------------------------------
cbind(Est=ilogitT(fit_ls$coef),
	    Lwr=ilogitT(fit_ls$coef-2*fit_ls$coef.se),
	    Upr=ilogitT(fit_ls$coef+2*fit_ls$coef.se))

## ------------------------------------------------------------------------
fit_ml <- move_sim(obs,para,rfit="R2",lltype="multinom")

## ------------------------------------------------------------------------
cbind(Est=ilogitT(fit_ml$coef),
	    Lwr=ilogitT(fit_ml$coef-2*fit_ml$coef.se),
	    Upr=ilogitT(fit_ml$coef+2*fit_ml$coef.se))

