#' Estimate migration between two areas from tag recaptures
#'
#' Estimate migration between two areas from tag recaptures, catches and
#' population sizes in each area.
#' @param obs A list that contains the following "yrs" years, "relT" vector
#' of annual tag releases, "R1" matrix of recaptures in area 1, "R2" matrix
#' of recaptures in area 2, "prec1" (catch numbers/vulnerable numbers) in area 1
#' , "prec2" (catch numbers/vulnerable numbers) in area 2, "M" instantanous
#' natural mortality rate, "relM" tag release mortality and "tagloss" annual tag
#' loss.
#' @param para starting value for the migration rate "move".
#' @param rfit The recaptures to fit to, "R1" fits to recaptures in area 1 only,
#'  "R2" fits the recaptures in area 2 only and "R12" fits to the recaptures in
#'  both areas.
#' @param lltype the fitting method: either "prop" (least squares) or "multinom"
#' (maximum likelihood).
#' @export
move_sim <- function(obs, para, rfit, lltype) {
  # extract the parameters and data
  yrs <- obs[["yrs"]]
  relT <- obs[["relT"]]
  R1 <- obs[["R1"]]
  R2 <- obs[["R2"]]
  prec1 <- obs[["prec1"]]
  prec2 <- obs[["prec2"]]
  M <- obs[["M"]]
  relM <- obs[["relM"]]
  tagloss <- obs[["tagloss"]]
  rfit <- rfit
  lltype <- lltype
  #move 	<- para[1]
  # Estimation
  nlogl	<- function(ps) {
    t_move <- matrix(0,nrow=length(yrs),ncol=length(yrs),dimnames=list(yrs,yrs))
    t1 <- t2 <- tt1 <- tt2 <- r1 <- r2 <- t_move

    ps 	<- ilogitT(ps)
    move <- ps[1]
    for (x in 1:length(yrs)) {
      t1[x,x]	<- relT[x]					# Release tags
      t1[x,x]	<- (1-relM) * t1[x,x]		# Apply release mortality
      for (y in x:length(yrs)) {
        t1[x,y]		<- t1[x,y] * exp(-0.5*M)	# Apply first 0.5 M
        t2[x,y]		<- t2[x,y] * exp(-0.5*M)
        t1[x,y]		<- t1[x,y] * exp(-tagloss[y-x+1])	# Tagloss: Probability of at least one tag
        t2[x,y]		<- t2[x,y] * exp(-tagloss[y-x+1])
        tt1[x,y]	<- t1[x,y]							# Store available tags in tt
        tt2[x,y]	<- t2[x,y]
        r1[x,y]		<- t1[x,y]*prec1[y]					# Recaptures (by Petersen), recapture probability = C/B
        r2[x,y]		<- t2[x,y]*prec2[y]
        if((y+1) <= length(yrs)) {						# Prevent error
          t_move[x,y+1]	<- (t1[x,y] - r1[x,y])*move			# Movement after removals [= T*move*(1-move)**(y-x-1)]
          t1[x,y+1]	<- t1[x,y] - r1[x,y] - t_move[x,y+1]	# = T*(1-move)**(y-x)
          t2[x,y+1]	<- t2[x,y] - r2[x,y] + t_move[x,y+1]	# = 1-T*(1-move)**(y-x)
        }
        t1[x,y]		<- t1[x,y] * exp(-0.5*M)	# Apply second 0.5 M
        t2[x,y]		<- t2[x,y] * exp(-0.5*M)
      }
    }
    for(x in 1:length(yrs)) r1[x,x] <- r2[x,x] <- 0	# Remove within-season recaptures
    nll <- 0
    if(lltype == "prop") {				# Simple proportional fit by least-square
      if(rfit=="R1") 	nll <- sum((R1-r1)**2)
      if(rfit=="R2") 	nll <- sum((R2-r2)**2)
      if(rfit=="R12") nll <- sum((R1-r1)**2) + sum((R2-r2)**2)
    }
    if(lltype == "multinom") {				# Multinomial likelihood: all recaptured & total un-recaptured
      for(x in 1:length(yrs)) {
        yy <- seq.int(x,length(yrs))
        # print(paste("year: ",x,sep=""))
        # print(c("seq.year: ",yy))
        if(rfit=="R1") 	obsN <- c(R1[x,yy],relT[x]-sum(R1[x,yy]))	# All obs recaptures & total un-recaptured
        if(rfit=="R2") 	obsN <- c(R2[x,yy],relT[x]-sum(R2[x,yy]))
        if(rfit=="R12") obsN <- c(R1[x,yy],R2[x,yy],relT[x]-sum(R1[x,yy])-sum(R2[x,yy]))
        if(rfit=="R1") 	expN <- c(r1[x,yy],relT[x]-sum(r1[x,yy]))	# All obs recaptures & total un-recaptured
        if(rfit=="R2") 	expN <- c(r2[x,yy],relT[x]-sum(r2[x,yy]))
        if(rfit=="R12") expN <- c(r1[x,yy],r2[x,yy],relT[x]-sum(r1[x,yy])-sum(r2[x,yy]))
        mm	<- dmultinom(x=obsN,size=relT[x],prob=expN)		# All obs recaptures & total un-recaptured
        if (mm != 0) 	nll <- nll - log(mm)
      }
    }
    #print(list(tt1=tt1,tt2=tt2,r1=r1,r2=r2,t_move=t_move,nll=nll))
    #print(paste("nll: ",nll,sep=""))
    #print("********************************************")
    nll
  }
  mn 	<- optim(para,nlogl,method="BFGS",hessian=TRUE)
  V 	<- try(solve(mn$hessian),silent=TRUE)
  if(!is.numeric(V)) V <- 0
  return(list(nll=mn$value,coef=mn$par,coef.se=sqrt(diag(V)),V=V,mn=mn))
  # print(paste("R1: ",sum((R1-r1)**2),"R2: ", sum((R2-r2)**2),sep=" "))
}
