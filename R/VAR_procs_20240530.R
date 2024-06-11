#****************************************************************
Init.var <- function(grid=1001,order.max.init=10,inv.roots.def=NULL,a.niter=5000,
		a.eps.max.for.UIA=1E-10,a.eps.for.roots=1E-5,a.eps.for.spectra=1E-4)
#********************************
{
a.var <- list()
a.var$grid <- grid
a.var$pmax.init <- order.max.init
a.var$niter<-a.niter
a.var$eps.max.for.UIA<-a.eps.max.for.UIA
a.var$eps.for.roots<-a.eps.for.roots
a.var$eps.for.spectra<-a.eps.for.spectra


if ((a.eps.for.roots<0) | (a.eps.for.roots>1) | (a.eps.max.for.UIA<0) | (a.eps.max.for.UIA>1) | (a.eps.for.spectra<0) | (a.eps.for.spectra>1))
		{
		stop("Arguments a.eps.for.roots, a.eps.max.for.UIA and a.eps.for.spectra should be in (0,1)")
		}
if (a.niter<=100) 
		{ 
		stop("Argument a.niter should be in >100")
		}
if (order.max.init<=0) 
		{ 
		stop("Argument order.max.init should be in >0")
		}


if (is.null(inv.roots.def))
	{
		a.matrix<- matrix(data=c(NA,NA,rep(1,13)),nrow=1,ncol=15)
		dimnames(a.matrix)<-list()
		dimnames(a.matrix)[[2]]<-c("radius","angle","det","cross","chi.1","chi.2","chi.1.prod.2","ma.1","ma.2","eta.1","eta.2","ksi.1","ksi.2","ksi.c","zeta")
		a.var$inv.roots<- as.data.frame(a.matrix)		
	} else if (is.data.frame(inv.roots.def))
		{
		a.var$inv.roots<- as.data.frame(inv.roots.def)
		} else if (file.exists(inv.roots.def))
			{
			a.var$inv.roots <- read.table(inv.roots.def,header=TRUE,na.strings="#N/A",fill=TRUE,comment.char="")
			} else
				{
				stop("Argument inv.roots.def should be an existing text file or a data.frame or NULL")
				}


a.var$inv.roots<-check.inv.roots(a.var$inv.roots,order.max.init)

class(a.var)<- "var"
a.var
}
#********************************************************************************************************************************
calculate.VAR <- function(a.var,calc.method="from.det.cross",M.fact=1.1,plot.spectra=TRUE,suppr.spec.check.warn=FALSE)
#********************************************************************************************************************************
{
if (is.null(a.var))
	{
	stop("First initialize a var calling init.var")
	} else if (is.null(a.var$inv.roots))
		{
		stop("First initialize a a.var$inv.roots passing to init.var the inv.roots data frame text file or NULL")
		} else
			{
			a.var$inv.roots<-check.inv.roots(a.var$inv.roots,a.var$pmax.init)
			method.ok<- FALSE
			if (calc.method=="from.det.cross")
				{
					a.var <- calculate.VAR.from.det.cross(a.var)
					method.ok<-TRUE
				} else if (calc.method=="from.eta.ksi.zeta")
				{
					if (M.fact<=1) 
						{
						stop("Argument M.fact should be >1")
						}
					a.var <- calculate.VAR.from.eta.ksi.zeta(a.var,M.fact)
					method.ok<-TRUE
				} else 
					{
					stop("Method of defining the VAR should be 'from.det.cross' or 'from.eta.ksi.zeta' ")
					}

			if (method.ok)
				{
				a.var$spec.1 <- a.var$det$inv.values$spec-a.var$chi.1$inv.values$spec-log(2*pi)
				a.var$spec.2 <- a.var$det$inv.values$spec-a.var$chi.2$inv.values$spec-log(2*pi)
				a.var$Coher<- a.var$chi.1.prod.2 $inv.values$spec-a.var$cross$inv.values$spec
				a.var <- calculate.Phase.etc(a.var)

				a.var <- calc.covs(a.var,a.var$order)
				a.var <- find.coefs(a.var)
				a.var <- calc.VAR.spec.from.coefs(a.var)

				if (plot.spectra) plot_VAR.spectra(a.var)
				a.var<-check.VAR.spectra.calculation(a.var,suppr.spec.check.warn)
				}
			}
class(a.var)<- "var"
a.var
}
#****************************************************************
plot_VAR.spectra <- function(a.var,both=TRUE)
#********************************
{
par.def<-par(no.readonly = TRUE)
on.exit(par(par.def))
par(mfrow=c(2,2))
par(ps=10)
par(mar=c(2,4, 1, 1) + 0.1)

freq.labels <- c("0","0.25 pi","0.5 pi","0.75 pi","pi")
freq.labels.where <-seq(from=0,to=1,length=5)*a.var$grid

if (both==TRUE)
{
plot.ts(as.ts(cbind(a.var$spec.1,
		log(Re(a.var$spec[(1:a.var$grid),1,1])))),plot.type="single",xaxt="n",ylab='log spectra',xlab='frequency',
		type="l", main="spec.1",col=c(4,2),lwd=c(2,1) )
axis(1,at=freq.labels.where,labels=freq.labels)
legend("topright",c("|chi.1|^2/|det|^2","from VAR coefs"),col=c(4,2),lwd=c(2,1) )



plot.ts(as.ts(cbind(a.var$spec.2,
		log(Re(a.var$spec[(1:a.var$grid),2,2])))),plot.type="single",xaxt="n",ylab='log spectra',xlab='frequency',
		type="l", main="spec.2",col=c(4,2), lwd=c(2,1) )
axis(1,at=freq.labels.where,labels=freq.labels)
legend("topright",c("|chi.2|^2/|det|^2","from VAR coefs"),col=c(4,2),lwd=c(2,1) )


help <- Mod(a.var$spec[(1:a.var$grid),1,2])^2 /(    Mod(a.var$spec[(1:a.var$grid),1,1])*Mod(a.var$spec[(1:a.var$grid),2,2])   )
plot.ts(    as.ts(cbind(exp(a.var$Coher),  help      )), plot.type="single",xaxt="n",ylab='in [0,1]',xlab='frequency',
		type="l", main="Coherency",col=c(4,2),lwd=c(2,1)  )
axis(1,at=freq.labels.where,labels=freq.labels)
legend("topright",c("|cross|^2/(|cross|^2+|det|^2)","from VAR coefs"),col=c(4,2),lwd=c(2,1) )


help <- -Arg(a.var$spec[(1:a.var$grid),1,2]) 
plot.ts(    as.ts(cbind((a.var$Phase),  help      )), plot.type="single",xaxt="n",ylab='in [-pi,pi]',xlab='frequency',
		type="l", main="Phase",col=c(4,2),lwd=c(2,1)  )
axis(1,at=freq.labels.where,labels=freq.labels)
legend("topright",c("Arg(cross)","from VAR coefs"),col=c(4,2),lwd=c(2,1) )
}

if (both==F)
{
plot.ts(as.ts(log(Re(a.var$spec[(1:a.var$grid),1,1]))),plot.type="single",xaxt="n",ylab='log spectra',xlab='frequency',
		type="l", main="spec.1",col=c(4,2),lwd=c(2,1) )
axis(1,at=freq.labels.where,labels=freq.labels)


plot.ts(as.ts(		log(Re(a.var$spec[(1:a.var$grid),2,2]))),plot.type="single",xaxt="n",ylab='log spectra',xlab='frequency',
		type="l", main="spec.2",col=c(4,2), lwd=c(2,1) )
axis(1,at=freq.labels.where,labels=freq.labels)


help <- Mod(a.var$spec[(1:a.var$grid),1,2])^2 /(    Mod(a.var$spec[(1:a.var$grid),1,1])*Mod(a.var$spec[(1:a.var$grid),2,2])   )
plot.ts(    as.ts(  help      ), plot.type="single",xaxt="n",ylab='in [0,1]',xlab='frequency',
		type="l", main="Coherency",col=c(4,2),lwd=c(2,1), ylim=c(0,1) )
axis(1,at=freq.labels.where,labels=freq.labels)


help <- -Arg(a.var$spec[(1:a.var$grid),1,2]) 
plot.ts(    as.ts(  help      ), plot.type="single",xaxt="n",ylab='in [-pi,pi]',xlab='frequency',
		type="l", main="Phase",col=c(4,2),lwd=c(2,1),ylim=c(-pi,pi)  )
axis(1,at=freq.labels.where,labels=freq.labels)
}


}

#****************************************************************
plot_VAR.Phase.details<- function(a.var)
#********************************
{
par.def<-par(no.readonly = TRUE)
on.exit(par(par.def))

par(mfrow=c(2,2))
freq.labels <- c("0","0.25 pi","0.5 pi","0.75 pi","pi")
freq.labels.where <-seq(from=0,to=1,length=5)*a.var$grid
plot.ts(as.ts(a.var$Phase),
		type="p", main="Phase lag",xaxt="n",xlab='frequency',ylab='in [-pi,pi]')
axis(1,at=freq.labels.where,labels=freq.labels)

plot.ts(as.ts(a.var$Phase.div.freq),
		type="p", main="Time lag = Phase.div.freq" ,xaxt="n",xlab='frequency',ylab='lag-lead time units')
axis(1,at=freq.labels.where,labels=freq.labels)

plot.ts(as.ts(cbind(a.var$group.delay,a.var$group.delay)),
		type="l", main="group.delay", plot.type="single",xaxt="n",xlab='frequency',ylab='')
axis(1,at=freq.labels.where,labels=freq.labels)

help <- Mod(a.var$spec[(1:a.var$grid),1,2])^2 /(    Mod(a.var$spec[(1:a.var$grid),1,1])*Mod(a.var$spec[(1:a.var$grid),2,2])   )
plot.ts(    as.ts(cbind(exp(a.var$Coher),  help      )), plot.type="single",xaxt="n",ylab='in [0,1]',xlab='frequency',
		type="l", main="Coherency",col=c(4,2),lwd=c(2,1)  )
axis(1,at=freq.labels.where,labels=freq.labels)

}



#****************************************************************
simulate.VAR<- function(a.var,sample.size, burn.in = 1000)
#********************************
{
p <- a.var$order
series.number <- dim(a.var$ar.list$ar)[3]
sigma <- a.var$ar.list$var.pred

innovations.iid <- array(data = NA, dim = c((burn.in+sample.size), series.number))
for (i in 1:series.number)
	innovations.iid[,i] <-rnorm((burn.in+sample.size))

e <- eigen(sigma)
V <- e$vectors
root <- V %*% diag((e$values)^(1/2)) %*% t(V)

innovations <- t( root %*% t(innovations.iid))
init.series <- array(data = 0, dim = c((burn.in+sample.size), series.number))
if (p==0) {init.series<-innovations} else
	{
	init.series[(1:p),]<- innovations [(1:p),]

	for ( s in (p+1):(burn.in+sample.size) )
		{
		for ( i in 1:p )
			{
			phi <- as.matrix(a.var$ar.list$ar[i,,])
			init.series[s,] <- init.series[s,] +  t(   phi %*% as.vector(t(init.series[s-i,]))  )
			}
		init.series[s,] <- init.series[s,]+ innovations [s,]
		}
	}

series <- array(data = NA, dim = c((sample.size), series.number))
series <- init.series[((burn.in+1):(burn.in+sample.size)),]
series 
}



#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************


#****************************************************************
check.inv.roots<- function(a.inv.roots,a.pmax.init)
#*****************************************************************
{

rows<-nrow(a.inv.roots)
cols<-ncol(a.inv.roots)
target.rows<- 6*a.pmax.init+1
if (rows < target.rows)
	for (i in (1: (target.rows-rows) ))
		a.inv.roots<-rbind(a.inv.roots, c(NA,NA,rep(0,(cols-2))) )


for (what in c("radius","angle","det","cross","chi.1","chi.2","chi.1.prod.2","ma.1","ma.2","eta.1","eta.2","ksi.1","ksi.2","ksi.c","zeta"))
	{
	which.column.from <-eval(parse(text=paste("a.inv.roots$",what,sep="")))
	if (is.null(which.column.from)) 
		{
		rows<-nrow(a.inv.roots)
		cols<-ncol(a.inv.roots)
		if ((what=="radius") | (what=="angle") )
			{a.inv.roots<-cbind(a.inv.roots,c(rep(NA,(rows))))} else
		a.inv.roots<-cbind(a.inv.roots,c(1,rep(0,(rows-1))))
		names(a.inv.roots)[cols+1]<-what
		}
	}

rows<-nrow(a.inv.roots)
cols<-ncol(a.inv.roots)
for (k in (2:rows))
	{
	if (!is.na(a.inv.roots$angle[k]))
	{
	if (a.inv.roots$angle[k]<0)
		{ 
		stop("Please change root",k-1,"! Angle should be in [0,pi]. Complex conjugates will be added automatically.")
		}
	if (a.inv.roots$angle[k]>pi) 
		{ 
		stop("Please change root",k-1,"! Angle should be in [0,pi].")
		}
	}

	if (!is.na(a.inv.roots$radius[k]))
	{
	if (a.inv.roots$radius[k]<=0) 
		{
		stop("Please change root",k-1,"! Radius should be in (0,1) (except for cross), since inverse roots are entered.")
		}

	if (a.inv.roots$radius[k]==1) 
		{
		stop("Please change root",k-1,"! Radius should be in (0,1) (except for cross), since inverse roots are entered.")

		}
	
	if ((a.inv.roots$radius[k]>1) & ((a.inv.roots$cross[k]==0)|(a.inv.roots$cross[k]!= sum(a.inv.roots[k,3:cols])))) 
		{
		stop("Please change root",k-1,". Radius should be in (0,1) (except for cross), since inverse roots are entered.")
		}

	}
	}

a.inv.roots$chi.1.prod.2<-a.inv.roots$chi.1+a.inv.roots$chi.2
a.inv.roots$chi.1.prod.2[1]<-1


min.ksi.1.2<-pmin(a.inv.roots$ksi.1,a.inv.roots$ksi.2)
a.inv.roots$ksi.1<-a.inv.roots$ksi.1-min.ksi.1.2
a.inv.roots$ksi.2<-a.inv.roots$ksi.2-min.ksi.1.2
a.inv.roots$ksi.c<-a.inv.roots$ksi.c+min.ksi.1.2
a.inv.roots$ksi.1[1]<-1
a.inv.roots$ksi.2[1]<-1
a.inv.roots$ksi.c[1]<-1

a.inv.roots
}

#****************************************************************
check.VAR.spectra.calculation<- function(a.var,suppr.spec.check.warn)
#*****************************************************************
{

v.1 <- -a.var$det$inv.values$spec
v.2 <- -a.var$cross$inv.values$spec
msg<-vector(mode = "list", length = 0)

msg[[1]] <- a.var$UIA.msg


eps2<-max(abs(a.var$chi.1.prod.2$inv.values$spec+log(exp(v.1)+exp(v.2))))
msg[[2]]<-paste("In log scale, |chi.1|^2*|chi.2|^2 differs from |det|^2+|cross|^2","by eps=",eps2," at most.")
if (!suppr.spec.check.warn)if (abs(eps2)>a.var$eps.for.spectra)warning(msg[[2]])

eps3<-max(abs(a.var$spec.1-log(a.var$spec[(1:a.var$grid),1,1])))
msg[[3]] <- paste("In log scale, spec.1 calculated from chi.1 & det differs from the one",
			       "calculated from the VAR coefficients by eps=",eps3," at most.")
if (!suppr.spec.check.warn) if (abs(eps3)>a.var$eps.for.spectra)warning(msg[[3]])


eps4<-max(abs(a.var$spec.2-log(a.var$spec[(1:a.var$grid),2,2])))
msg[[4]] <- paste("In log scale, spec.2 calculated from chi.2 & det differs from the one",
					"calculated from the VAR coefficients by eps=",eps4," at most.")
if (!suppr.spec.check.warn) if (abs(eps4)>a.var$eps.for.spectra) warning(msg[[4]])


help <- Mod(a.var$spec[(1:a.var$grid),1,2])^2 /(    Mod(a.var$spec[(1:a.var$grid),1,1])*Mod(a.var$spec[(1:a.var$grid),2,2])   )
eps5<-max(abs(exp(a.var$Coher)/help)-1)
msg[[5]] <- paste("In relative scale, the coherency calculated from cross & det differs from the one",
			"calculated from the VAR coefficients by eps=",eps5," at most.")
if (!suppr.spec.check.warn) if (abs(eps5)>a.var$eps.for.spectra) warning(msg[[5]])

help <- -Arg(a.var$spec[(1:a.var$grid),1,2]) 
eps6<-max((sin(a.var$Phase)-sin(help))^2)
msg[[6]] <- paste("In sin^2 scale, the Phase calculated from cross differs from the one","calculated from the VAR coefficients by eps=",eps6," at most.")
if (!suppr.spec.check.warn) if (abs(eps6)>a.var$eps.for.spectra) warning(msg[[6]])

a.var$validity.msg<-msg
a.var
}


#*******************************************************************
calculate.Phase.etc <- function(a.var)
#*******************************************************************
{
lamdas <- a.var$det$inv.values$freq*pi*2
b.lamdas <- complex(length.out=a.var$grid,modulus=1,argument=-lamdas)
values <- complex(length.out=a.var$grid,modulus=0,argument=lamdas)
p <- a.var$cross$order
coefs <- a.var$cross$coefs
for (i in 1:(p+1))
	{
	c.lamdas <- b.lamdas^(i-1)
	values <- values + coefs[i]*c.lamdas
	}
d.lamdas <- Conj(b.lamdas^(a.var$order))

a.var$Phase	<- Arg(d.lamdas*values) 
a.var$Phase.div.freq	<- a.var$Phase/lamdas

a.var$group.delay<- rep(a.var$order,a.var$grid)
lamdas <- a.var$cross$inv.values$freq*pi*2
b.lamdas <- complex(length.out=a.var$grid,modulus=1,argument=lamdas)

j<-0
if (a.var$cross$inv.roots.number>0)
for (i in 1:a.var$cross$inv.roots.number)
	{
	j <-j+1
	root <- complex(modulus=a.var$cross$inv.roots[i,1],argument=a.var$cross$inv.roots[i,2])
	c.lamdas <- rep(root,a.var$grid)
	d.lamdas <- b.lamdas*c.lamdas
	a.var$group.delay <- a.var$group.delay + Re(d.lamdas/(1-d.lamdas))
	if (Im(root)!=0)
		{
		j <-j+1
		root <- Conj(root)
		c.lamdas <- rep(root,a.var$grid)
		d.lamdas <- b.lamdas*c.lamdas
		a.var$group.delay <- a.var$group.delay +Re(d.lamdas/(1-d.lamdas))
		}
	}
a.var
}


#****************************************************************

#****************************************************************
find.coefs <- function(a.var)
#********************************
{
	p <- a.var$order
	gammas <- array(data = NA, dim = c(p + 1, 2, 2))
	for (u in (0:p)) 
		{
		gammas[(u + 1), 1, 1] <- a.var$cov.1[(u + 1)]
		gammas[(u + 1), 2, 2] <- a.var$cov.2[(u + 1)]
		gammas[(u + 1), 2, 1] <- a.var$cov.cross[(u + 1)]
		gammas[(u + 1), 1, 2] <- a.var$cov.cross.neg[(u + 1)]
		#print(gammas[(u + 1), ,] )
		}

if (p>0)
{
	big.gamma.matrix <- matrix(data=NA, nrow=2*p,ncol=2*p)
	for (i in (1:(p))) 
	for (j in (1:(p))) 
		{
		big.gamma.matrix[(2*(i-1)+1),(2*(j-1)+1)] <- gammas[(abs(i-j) + 1), 1, 1] 
		big.gamma.matrix[(2*(i-1)+2),(2*(j-1)+2)] <- gammas[(abs(i-j) + 1), 2, 2] 
		if (i>=j)
			{
			big.gamma.matrix[(2*(i-1)+1),(2*(j-1)+2)] <- gammas[(abs(i-j) + 1), 2, 1] 
			big.gamma.matrix[(2*(i-1)+2),(2*(j-1)+1)] <- gammas[(abs(i-j) + 1), 1, 2] 
			}
		else
			{
			big.gamma.matrix[(2*(i-1)+1),(2*(j-1)+2)] <- gammas[(abs(i-j) + 1), 1, 2] 
			big.gamma.matrix[(2*(i-1)+2),(2*(j-1)+1)] <- gammas[(abs(i-j) + 1), 2, 1] 
			}
		}
	#print("big.gamma.matrix")
	#print(big.gamma.matrix)


	big.gamma.vector <- matrix(data=NA, nrow=2*p,ncol=2)
	for (i in (1:p)) 
		{
		big.gamma.vector [(2*(i-1)+1),1] <- gammas[(i + 1), 1, 1] 
		big.gamma.vector [(2*(i-1)+1),2] <- gammas[(i + 1), 1, 2] 
		big.gamma.vector [(2*(i-1)+2),1] <- gammas[(i + 1), 2, 1] 
		big.gamma.vector [(2*(i-1)+2),2] <- gammas[(i + 1), 2, 2] 
		}

	big.gamma.vector.neg <- matrix(data=NA, nrow=2*p,ncol=2)
	for (i in (1:p)) 
		{
		big.gamma.vector.neg  [(2*(i-1)+1),1] <- gammas[(i + 1), 1, 1] 
		big.gamma.vector.neg  [(2*(i-1)+1),2] <- gammas[(i + 1), 2, 1] 
		big.gamma.vector.neg  [(2*(i-1)+2),1] <- gammas[(i + 1), 1, 2] 
		big.gamma.vector.neg  [(2*(i-1)+2),2] <- gammas[(i + 1), 2, 2] 
		}
	#print("big.gamma.vector")
	#print(big.gamma.vector)
	phi.vector.tr <-  solve(t(big.gamma.matrix)) %*% big.gamma.vector.neg
	#print("phi.vector.tr")
	#print(phi.vector.tr)
}
	
	a.var$ar.list <-list()
	a.var$ar.list$order <- a.var$order
if (p>0)
{
	a.var$ar.list$ar <- array(data = NA, dim = c(p, 2, 2))
	for (i in (1:p)) 
		{
		a.var$ar.list$ar [i , 1, 1] <-  phi.vector.tr[(2*(i-1)+1),1] 
		a.var$ar.list$ar [i , 2, 1] <-  phi.vector.tr[(2*(i-1)+1),2] 
		a.var$ar.list$ar [i , 1, 2] <-  phi.vector.tr[(2*(i-1)+2),1] 
		a.var$ar.list$ar [i , 2, 2] <-  phi.vector.tr[(2*(i-1)+2),2] 
		#print(a.var$ar.list$ar [i , , ] )
		}

	sigma <- gammas[1, , ] -  t(t(big.gamma.vector.neg)  %*% phi.vector.tr)
}
else
{
	a.var$ar.list$ar <- array(data = 0, dim = c(1, 2, 2))
	sigma <- gammas[1, , ] 
}
	a.var$ar.list$var.pred <- sigma
	#print("sigma")

	#print(sigma)
a.var
}

#****************************************************************
calc.covs<- function(a.var,umax)
#********************************
{
a.var$cov.1 <- rep(NA,(umax+1))
a.var$cov.2 <- rep(NA,(umax+1))
a.var$cov.cross <- rep(NA,(umax+1))
a.var$cov.cross.neg <- rep(NA, (umax + 1))

N <-a.var$grid


help <-  seq(from=0,to=0.5,length=a.var$grid)
lamdas <- help*pi*2
for (u in 0:(umax))
	{
	values <- cos(u*lamdas)*exp(a.var$spec.1) 
	a.var$cov.1[(u+1)] <- 2*pi* (sum(values)-(values[1]+values[N])/2)/(N-1)

	values <- cos(u*lamdas)*exp(a.var$spec.2)
	a.var$cov.2[(u+1)] <- 2*pi* (sum(values)-(values[1]+values[N])/2)/(N-1)

	help <- (   exp(a.var$Coher) * (exp(a.var$spec.1) )* (exp(a.var$spec.2) )   )^(1/2) 
	b.lamdas <- complex(length.out=N,modulus=1,argument=u*lamdas)
	help.values <- complex(length.out=N,modulus= help,argument=a.var$Phase)
	values <- Re(help.values * b.lamdas )
	a.var$cov.cross[(u+1)] <- 2*pi* (sum(values)-(values[1]+values[N])/2)/(N-1)

	b.lamdas <- complex(length.out = N, modulus = 1, argument =  - u * lamdas)
	help.values <- complex(length.out = N, modulus=help,argument = a.var$Phase)
	values <- Re(help.values * b.lamdas)
	a.var$cov.cross.neg[(u + 1)] <- 2*pi* ( (sum(values) - (values[1] + values[N])/2))/(N - 1)
	}
a.var

}

#****************************************************************
calc.VAR.spec.from.coefs<- function(a.var)
#********************************
{
	series.number<- dim(a.var$ar.list$ar)[3]
	a.var$freq <- seq(from=0,to=0.5,length=a.var$grid)*pi*2
	a.var$spec <- array(complex(length.out=a.var$grid*series.number*series.number,modulus=0),
											dim=c(a.var$grid,series.number,series.number))

	lamdas <- a.var$freq 
	b.lamdas <- complex(length.out=a.var$grid,modulus=1,argument=-lamdas)
	phi <- array(complex(length.out=a.var$grid*series.number*series.number,modulus=0),
											dim=c(a.var$grid,series.number,series.number))
	p <- a.var$order
	for (i in 1:series.number)
		phi[,i,i] <- rep(1,a.var$grid)

	 if (p>0) for (i in 1:p)
		{
		coefs <- a.var$ar.list$ar[i,,]
		c.lamdas <- b.lamdas^i
		my.mult <- function(x) { x*c.lamdas }
		phi <- phi - apply(coefs,c(1,2),my.mult )
		}

	inv_sigma<- solve(a.var$ar.list$var.pred)
	my.mult.transp.phi.inv_sigma.phi <- function(x) { as.matrix(solve(t(Conj(as.matrix(x)))%*%inv_sigma%*%as.matrix(x)))/(2*pi) }
	for (k in 1:a.var$grid)
		a.var$spec[k,,] <- my.mult.transp.phi.inv_sigma.phi(phi[k,,])

	a.var
}

#********************************************************************************************************************************
calculate.VAR.from.det.cross <- function(a.var)
#********************************************************************************************************************************
{
#read and build det and cross from roots stored in data frame "a.var$inv.roots"
a.var$det <- Specify.Var.polynom(a.var$det,"det",a.var)
a.var$cross <- Specify.Var.polynom(a.var$cross,"cross",a.var)
check.orders.det.cross(a.var)

a.var <- Calculate.forced.Product.polynom(a.var)

#calculates the LCD of (det,cross) (they should also be roots of chi1*chi2),
#and store them in column "chi.1.prod.2" of a.var$inv.roots
#****************

a.var$chi.1.prod.2.rest <- Calculate.rest.product.polynom.fourier.coefs(
				a.var$det$fourier.coefs+a.var$cross$fourier.coefs,
				max(a.var$det$order,a.var$cross$order),
				a.var$chi.1.prod.2.forced$fourier.coefs,
				a.var$chi.1.prod.2.forced$order)

#takes the fourier coefs (covs) of det^2+cross^2 and chi.1.prod.2.forced^2, 
#writes each one of them as inverse fourier transform of its fourier coefs
#and calculates the quotient between them (nom/denom). It returns its fourier coefficients
#(its roots should also be roots of chi1 or of chi2. This be enforced below)
#****************

a.var <- Calculate.Product.polynom(a.var$chi.1.prod.2.rest,a.var)

#finds the coefficients of chi.1.prod.2.rest (MA polynom) given the fourier coefs of chi.1.prod.2.rest^2 (covs)
#remider: chi.1.prod.2.rest:=det^2+cross^2/chi.1.prod.2.forced^2, where chi.1.prod.2.forced:=LCD(det,cross).
#Then finds the roots of chi.1.prod.2.rest. Keeps only roots with positive Arg and forces them into the UC. 
#These are then added to chi.1.prod.2 column in a.var$inv.roots, which already contains the roots of chi.1.prod.2.forced
#****************

a.var$chi.1.prod.2 <- Specify.Var.polynom(a.var$chi.1.prod.2,"chi.1.prod.2",a.var)
#rebuild chi1*chi2 from stored roots, then adjust constant appropriatelly
#****************

a.var <- distribute.roots.to.chi1.2(a.var)

#define chi.1 and chi.2 once their product is defined
#adjusts in a.var$inv.roots the order of the roots of chi.1 and chi.2, 
#so that their sum matches the orders of chi.1.prod.2 
#priority is given to chi.1, unless its order is higher than of chi.1.prod.2, 
#in which case it is set to it, while chi.2 =0
#otherwise, missing roots are added to chi.2, the initial values of which are ignored

#****************

a.var$chi.1 <- Specify.Var.polynom(a.var$chi.1,"chi.1",a.var)
a.var$chi.2 <- Specify.Var.polynom(a.var$chi.2,"chi.2",a.var)
a.var$order <- max(a.var$chi.1$order,a.var$chi.2$order,
			ceiling(a.var$det$order/2),ceiling(a.var$cross$order/2))
a.var$det <- calc.values.new(a.var$det,a.var$grid)
a.var$cross <- calc.values.new(a.var$cross,a.var$grid)
a.var$chi.1.prod.2 <- calc.values.new(a.var$chi.1.prod.2,a.var$grid)
a.var$chi.1 <- calc.values.new(a.var$chi.1,a.var$grid)
a.var$chi.2 <- calc.values.new(a.var$chi.2,a.var$grid)

a.var
}
#********************************************************************************************************************************
calculate.VAR.from.eta.ksi.zeta <- function(a.var,M.fact)
#********************************************************************************************************************************
{
#read roots of eta, ksi, zeta and calculate their values
a.var$eta.1 <- Specify.Var.polynom(a.var$eta.1,"eta.1",a.var)
a.var$eta.2 <- Specify.Var.polynom(a.var$eta.2,"eta.2",a.var)
a.var$ksi.1 <- Specify.Var.polynom(a.var$ksi.1,"ksi.1",a.var)
a.var$ksi.2 <- Specify.Var.polynom(a.var$ksi.2,"ksi.2",a.var)
a.var$ksi.c <- Specify.Var.polynom(a.var$ksi.c,"ksi.c",a.var)
a.var$zeta <- Specify.Var.polynom(a.var$zeta,"zeta",a.var)

check.orders.eta.ksi.zeta(a.var)


a.var$eta.1 <- calc.values.new(a.var$eta.1,a.var$grid)
a.var$eta.2 <- calc.values.new(a.var$eta.2,a.var$grid)
a.var$ksi.1 <- calc.values.new(a.var$ksi.1,a.var$grid)
a.var$ksi.2 <- calc.values.new(a.var$ksi.2,a.var$grid)
a.var$ksi.c <- calc.values.new(a.var$ksi.c,a.var$grid)
a.var$zeta <- calc.values.new(a.var$zeta,a.var$grid)

#find constant M for zeta, chi.1 and chi.2, ensuring that Mod(cross)^2:=Mod(chi.1*chi2)^2-Mod(det)^2 is positive
M<-M.fact*sqrt(max(a.var$ksi.c$values$spec/(a.var$eta.1$values$spec*a.var$eta.2$values$spec*a.var$zeta$values$spec)))


#set constant M appropriatelly, define and recalculate chi.1,chi.2, det
#a.var$inv.roots[1,"zeta"]<-M
#a.var$zeta <- Specify.Var.polynom(a.var$zeta,"zeta",a.var)
#a.var$zeta <- calc.values.new(a.var$zeta,a.var$grid)

a.var$inv.roots[,"chi.1"]<-a.var$inv.roots[,"eta.1"]+a.var$inv.roots[,"ksi.2"]+a.var$inv.roots[,"zeta"]
a.var$inv.roots[1,"chi.1"]<-M
a.var$chi.1 <- Specify.Var.polynom(a.var$chi.1,"chi.1",a.var)
a.var$chi.1 <- calc.values.new(a.var$chi.1,a.var$grid)

a.var$inv.roots[,"chi.2"]<-a.var$inv.roots[,"eta.2"]+a.var$inv.roots[,"ksi.1"]+a.var$inv.roots[,"zeta"]
a.var$inv.roots[1,"chi.2"]<-M
a.var$chi.2 <- Specify.Var.polynom(a.var$chi.2,"chi.2",a.var)
a.var$chi.2 <- calc.values.new(a.var$chi.2,a.var$grid)

a.var$inv.roots[,"det"]<-a.var$inv.roots[,"ksi.1"]+a.var$inv.roots[,"ksi.2"]+a.var$inv.roots[,"ksi.c"]+a.var$inv.roots[,"zeta"]
a.var$inv.roots[1,"det"]<-M
a.var$det <- Specify.Var.polynom(a.var$det,"det",a.var)
a.var$det <- calc.values.new(a.var$det,a.var$grid)


#define and calculate chi.1*chi.2
a.var$inv.roots[,"chi.1.prod.2"]<-a.var$inv.roots[,"chi.1"]+a.var$inv.roots[,"chi.2"]
a.var$inv.roots[1,"chi.1.prod.2"]<-M^2
a.var$chi.1.prod.2 <- Specify.Var.polynom(a.var$chi.1.prod.2,"chi.1.prod.2",a.var)
a.var$chi.1.prod.2 <- calc.values.new(a.var$chi.1.prod.2,a.var$grid)

#find the Fourier coefficients of Mod(cross)^2:=Mod(chi.1*chi.2)^2-Mod(det)^2
a.var$tent.cross<-list()
a.var$tent.cross$fourier.coefs<-a.var$chi.1.prod.2$fourier.coefs-a.var$det$fourier.coefs
a.var$tent.cross$order<-max(a.var$chi.1.prod.2$order,a.var$det$order)

a.var$inv.roots[,"chi.1.prod.2"]<-0
a.var$inv.roots[1,"chi.1.prod.2"]<-1
a.var$chi.1.prod.2 <- Specify.Var.polynom(a.var$chi.1.prod.2,"chi.1.prod.2",a.var)

#find the roots of cross from the Fourier coefficients of Mod(cross)^2,
#using the Innovations Algorithm. These are temporarily stored in column "chi.1.prod.2" of VAR.inv.roots
a.var <- Calculate.Product.polynom(a.var$tent.cross,a.var,invert=TRUE)
#take over the roots of cross from column "chi.1.prod.2" of VAR.inv.roots and calculate it
a.var$inv.roots[,"cross"]<-a.var$inv.roots[,"chi.1.prod.2"]
a.var$cross <- Specify.Var.polynom(a.var$cross,"cross",a.var)
a.var$cross <- calc.values.new(a.var$cross,a.var$grid)

#set the order of the VAR
a.var$order <- max(a.var$chi.1$order,a.var$chi.2$order,
			ceiling(a.var$det$order/2),ceiling(a.var$cross$order/2))

#for compleetness recalculate chi.1*chi*2
a.var$inv.roots[,"chi.1.prod.2"]<-a.var$inv.roots[,"chi.1"]+a.var$inv.roots[,"chi.2"]
a.var$inv.roots[1,"chi.1.prod.2"]<-M^2
a.var$chi.1.prod.2 <- Specify.Var.polynom(a.var$chi.1.prod.2,"chi.1.prod.2",a.var)
a.var$chi.1.prod.2 <- calc.values.new(a.var$chi.1.prod.2,a.var$grid)

a.var
}

#****************************************************************
Calculate.forced.Product.polynom <- function(a.var)
#********************************
{
#forced.polynom contains the roots which are common to det and cross, i.e. is their least common denominator
#they should also be roots of chi1*chi2. 
#The roots are stored in column "chi.1.prod.2" of dataframe a.var$inv.roots

forced.polynom<-list()
a.var$inv.roots$chi.1.prod.2[1] <-1
row.of.last.root.in.frame <- 2
for (k in 2:(6*a.var$pmax.init+1))
	{
	if (is.na(a.var$inv.roots$radius[k])==FALSE)
		{
		row.of.last.root.in.frame <-k
		a.var$inv.roots$chi.1.prod.2[k] <-min(a.var$inv.roots$det[k],a.var$inv.roots$cross[k])
		}
	else
		a.var$inv.roots$chi.1.prod.2[k] <-0

	}
a.var$chi.1.prod.2.forced<- Specify.Var.polynom(forced.polynom,"chi.1.prod.2",a.var)
a.var
}


#****************************************************************
Calculate.rest.product.polynom.fourier.coefs <- function(nom.coefs,nom.order,denom.coefs,denom.order)
#********************************
{
#takes the fourier coefs (covs) of two positive polynomials on the UC (nom,denom), 
#writes each of them as inverse fourier transform of the fourier coefs
#and calculates the quotient between them (nom/denom). Then it returns its fourier coefficients

quotient <- list()
my.nom<-list()
my.nom$order<-2*nom.order
my.nom$coefs<-rep(0,(my.nom$order+1))
for (i in 1:(nom.order+1))
	{
	my.nom$coefs[nom.order+i]<-nom.coefs[i]
	my.nom$coefs[nom.order-i+2]<-nom.coefs[i]
	}
my.denom<-list()
my.denom$order<-2*denom.order
my.denom$coefs<-rep(0,(my.denom$order+1))
for (i in 1:(denom.order+1))
	{
	my.denom$coefs[denom.order+i]<-denom.coefs[i]
	my.denom$coefs[denom.order-i+2]<-denom.coefs[i]
	}
my.quot<-divide.polynom(my.nom,my.denom)

quotient$order<-(my.quot$order)/2
quotient$fourier.coefs<-rep(0,(quotient$order+1))
quotient$check.fourier.coefs<-rep(0,(quotient$order+1))
for (i in 1:(quotient$order+1))
	{
	quotient$fourier.coefs[i]<-my.quot$coefs[quotient$order+i]
	quotient$check.fourier.coefs[i]<-my.quot$coefs[quotient$order-i+2]
	}

quotient
}
#****************************************************************
reset.roots.stored.in.dataframe <- function(inv.roots,inv.roots.number,what,a.var)
#********************************
{
#called with what=chi.1.prod.2
#for each root in inv.roots checks whether in appears among the roots listed in a.var$inv.roots
#(including those added by this function call
# if it does, it augments its order (degree) by 1
# if it odes not, it adds the root to the list and sets the order to 1
# the new roots are added directly to a.var$inv.roots
# the new orders are returned in the value of the function and a.var$inv.roots should then be set to it.
# calling a.var$inv.roots$chi.1.prod.2 <- reset.roots.stored.in.dataframe(my.polynom$inv.roots,my.polynom$inv.roots.number,"chi.1.prod.2")
  

#cat("entered reset.roots.stored.in.dataframe ","\n")

which.column.to <-eval(parse(text=paste("a.var$inv.roots$",what,sep="")))
row.of.last.root.in.frame <- 2
for (k in 2:(6*a.var$pmax.init+1))
	{
	if (is.na(a.var$inv.roots$radius[k])==FALSE)
		row.of.last.root.in.frame <-k
	}
if (inv.roots.number>0)
for (j in 1:inv.roots.number)
	{
	found <- 0
	k <- 1
	while ((found==0) & (k<row.of.last.root.in.frame))
		{
		k <- k+1
		if ( (abs(inv.roots[j,1]-a.var$inv.roots$radius[k])<a.var$eps.for.roots) & (abs(inv.roots[j,2]-a.var$inv.roots$angle[k])<a.var$eps.for.roots) )
			{
			found <- 1
			which.column.to[k] <- which.column.to[k]+1
			}
		}
		if (found==0)
			{
			row.of.last.root.in.frame <- row.of.last.root.in.frame+1
			k <- row.of.last.root.in.frame
			a.var$inv.roots$radius[k] <- inv.roots[j,1]
			a.var$inv.roots$angle[k] <- inv.roots[j,2]
			which.column.to[k] <- 1
			}

	}
a.var$inv.roots$chi.1.prod.2<-which.column.to
a.var
}
#****************************************************************
Calculate.Product.polynom <- function(my.polynom,a.var,invert=TRUE)
#********************************
{
#finds the coefficients of my.polynom (MA polynom) given the fourier coefs of my.polynom^2 (covs), using the Univ. Innovations algor.
#here applied to chi.1.prod.2.rest:=det^2+cross^2/chi.1.prod.2.forced^2, where initially chi.1.prod.2.forced:=LCD(det,cross).
#Then finds the roots of my.polynom. Keeps only roots with positive Arg and forces them into the UC. 
#These are then added to chi.1.prod.2 column in a.var$inv.roots, which already contains the roots of chi.1.prod.2.forced
#


#cat("started calculating ch1*ch2-rest from its fourier coefficients","\n")
my.polynom <- Univariate.Innovations.Algorithm(my.polynom,my.polynom$fourier.coefs,my.polynom$order,a.var$niter,a.var$eps.max.for.UIA)
my.polynom$const <- my.polynom$coefs[1]

my.polynom$help.roots <- polyroot(my.polynom$coefs)
order<-length(my.polynom$help.roots)
my.polynom$inv.roots <- matrix(data=NA,nrow=order+1, ncol=2)
k<-0
if (order>0)
for (j in 1:order)
	{
	root <-my.polynom$help.roots[j]
#cat(j,order,root,1/Mod(root),Arg(root),"\n")
if (invert)
	{
	if (Arg(root)>=-a.var$eps.for.roots)
		{
		k <- k+1
		if (abs(Arg(root))<a.var$eps.for.roots)
			my.polynom$inv.roots[k,2] <- 0
		else
			my.polynom$inv.roots[k,2] <- Arg(root)
		if (Mod(root)>1)
			my.polynom$inv.roots[k,1] <- 1/Mod(root)
		else
			my.polynom$inv.roots[k,1] <- Mod(root)
		}
	}
	}

my.polynom$inv.roots.number <- k
a.var <- reset.roots.stored.in.dataframe(my.polynom$inv.roots,my.polynom$inv.roots.number,"chi.1.prod.2",a.var)
a.var$inv.roots$chi.1.prod.2[1] <- my.polynom$const

#cat("finished calculating rest polynomial","\n")
a.var$chi.1.prod.2.rest<-my.polynom 
a.var$UIA.msg<-my.polynom$msg
a.var
}



#****************************************************************
distribute.roots.to.chi1.2 <- function(a.var)
#********************************
#adjusts in a.var$inv.roots the order of the roots of chi.1 and chi.2, 
#so that their sum matches the oreders of chi.1.prod.2 
#priority is given to chi.1, unless its order i higher than of chi.1.prod.2, 
#in which case it is set to it, while chi.2 =0
#otherwise, missing roots are added to chi.2, the initial values of which are ignored
{
if (a.var$inv.roots$chi.1[1]==0)
	{
	if (a.var$inv.roots$chi.1.prod.2[1]!=0)
		a.var$inv.roots$chi.1[1] <-1	
	}
else 
	a.var$inv.roots$chi.2[1] <- a.var$inv.roots$chi.1.prod.2[1]/a.var$inv.roots$chi.1[1]

row.of.last.root.in.frame <- 2
for (k in 2:(6*a.var$pmax.init+1))
	{
	if (is.na(a.var$inv.roots$radius[k])==FALSE)
		row.of.last.root.in.frame <-k
	}
count <-0 
for (k in 2:row.of.last.root.in.frame)
	{
	a.var$inv.roots$chi.2[k]<-0
	if (a.var$inv.roots$chi.1[k]>a.var$inv.roots$chi.1.prod.2[k])
		a.var$inv.roots$chi.1[k] <- a.var$inv.roots$chi.1.prod.2[k]

	diff <-a.var$inv.roots$chi.1.prod.2[k]-a.var$inv.roots$chi.1[k]
	while ((count<a.var$pmax.init) & (diff>0))
		{
		a.var$inv.roots$chi.2[k] <- a.var$inv.roots$chi.2[k]+1
		count <-count+1
		diff<-diff-1
		}
	a.var$inv.roots$chi.1[k] <- a.var$inv.roots$chi.1[k]+diff
	}


a.var
}
#****************************************************************
Specify.Var.polynom <- function(my.polynom,what,a.var)
#********************************
{
#reads roots and constant (row=1) from column "what" from dataframe "a.var$inv.roots"
#and calculates the coefficients and the fourier coefficients of poynomial indicated by "what"
#order (degree) of polynomial should be <= 2*a.var$pmax.init

my.polynom <- list()
which.column.from <-eval(parse(text=paste("a.var$inv.roots$",what,sep="")))
my.polynom$const <- which.column.from[1]


my.polynom$inv.roots <- matrix(data=NA,nrow=sum(which.column.from[2:length(which.column.from)]), ncol=2)



my.polynom$inv.roots.number <- 0
for (k in 2:(6*a.var$pmax.init+1))
	{
	if ((which.column.from[k]>0) & (is.na(a.var$inv.roots$radius[k])==FALSE) & (is.na(a.var$inv.roots$angle[k])==FALSE))

		{
		for (j in 1:(which.column.from[k]))
			{
			my.polynom$inv.roots.number <- my.polynom$inv.roots.number+1
			r <- my.polynom$inv.roots.number
			my.polynom$inv.roots[r,1] <- a.var$inv.roots$radius[k]
			my.polynom$inv.roots[r,2] <- a.var$inv.roots$angle[k]
			}
		}
	}
my.polynom  <- calculate.coefficients(my.polynom)
my.polynom  <- calculate.fourier.coefficients(my.polynom,a.var)
if (my.polynom$order>2*a.var$pmax.init)
	{
	stop("Order of",what,"is too high. order=",my.polynom $order,", while order.max=",2*a.var$pmax.init)
	}
my.polynom 
}


#****************************************************************
divide.polynom <- function(nominator,denominator)
#********************************
{
quotient<-list()
pol.nom <- nominator$coefs
p.nom <- nominator$order

pol.denom <- denominator$coefs
p.denom <- denominator$order

p <-p.nom-p.denom
quotient$order <- p
pol.quot <- rep(0,(p+1))


#cat("division:nominator",p.nom,pol.nom,"\n")
while (p>=0)
	{
	
	pol.quot[p+1] <-pol.nom[p.nom+1]/pol.denom[p.denom+1]
	for (i in (p.nom-p.denom+1):(p.nom+1))
		{
		pol.nom[i] <-pol.nom[i] -pol.quot[(p+1)]*pol.denom[(p.denom+i-p.nom)]
		}
	#cat("division:quotient",p,pol.quot,"\n")
	#cat("division:nominator",p.nom,pol.nom,"\n")
	p.nom<-p.nom-1
	p <- p-1
	}

quotient$coefs <- pol.quot
quotient
}	
#****************************************************************
calculate.coefficients <- function(a.polynom)
#********************************
{
order <- 0
coefs <- c(1, rep(0,2*a.polynom$inv.roots.number+1) )
k <- 0
while (k<a.polynom$inv.roots.number)
	{
	k <- k+1
	
	if ((a.polynom$inv.roots[k,2]==0) || (a.polynom$inv.roots[k,2]==pi))
		{
		order <- order + 1
		j <- order+1
		
		if (a.polynom$inv.roots[k,1]==0) 
			{
			while (j>1)
				{
				coefs[j] <- coefs[j-1]

				j <- j-1
				}
				coefs[1] <- 0
			} else
		
			{
			while (j>1)
				{
				coefs[j] <- coefs[j] -
						coefs[j-1]*a.polynom$inv.roots[k,1]*
							cos(a.polynom$inv.roots[k,2])
				j <- j-1
				}
			}
		} else
	
		{
		order <- order + 2
		j <- order+1
		while (j>2)
			{
			coefs[j] <- coefs[j] -
					 2*coefs[j-1]*a.polynom$inv.roots[k,1]*
								cos(a.polynom$inv.roots[k,2])+
					 coefs[j-2]*((a.polynom$inv.roots[k,1])^2)
			j <- j-1
			}
		coefs[2] <- coefs[2] -
					 2*a.polynom$inv.roots[k,1]*cos(a.polynom$inv.roots[k,2])
		} 

	} 
a.polynom$order <- order
a.polynom$coefs <- a.polynom$const*coefs
a.polynom
}
#****************************************************************
Univariate.Innovations.Algorithm <- function(my.polynom,my.coefs,q,a.niter,a.eps.max.for.UIA)
#********************************
{
#cat("started Univ. Innov. algorithm","\n")
eps <- 1
is.a.zero.polynom <- 0
if ( my.coefs[1]==0)
	{
	eps <- 0
	is.a.zero.polynom<-1
	msg<-"Univ. Innov. Algor.: four.coef[1]=0"
	for (j in (2:(q+1)))
		{
		if ( my.coefs[j]!=0)
			{
			msg<-paste(msg,", but (invalid input): four.coef[")
			msg<-paste(msg,j)
			msg<-paste(msg,"]=")
			msg<-paste(msg,my.coefs[j])
			stop(msg)
			}
		}
	}

v <-rep(0,a.niter+1)
new.coefs <- rep(0,a.niter+1)
for (i in 1:(q+1))
	new.coefs[i] <- my.coefs[i]
if (q<1) q <- 1
theta <- matrix(data=0,nrow=a.niter,ncol=a.niter)
v[1] <- new.coefs[1]



n <- 0
while ((n<a.niter)&(eps>a.eps.max.for.UIA))
	{
	n <- n+1
	for (k in 0:(n-1))
		{
		if (k>(n-q-1))
			{
			theta[n,(n-k)] <- new.coefs[(n-k+1)]/v[(k+1)]
			start.j <- n-q
			if (start.j<0) 
				start.j <-0 
			if (k>start.j )
			for (j in start.j:(k-1))
				{
				theta[n,(n-k)] <- theta[n,(n-k)]-theta[k,(k-j)]*theta[n,(n-j)]*v[j+1]/v[k+1]
				}
			}
		}
	v[n+1] <- new.coefs[1]
	start.j <- n-q
	if (start.j<0) 
		start.j <-0 

	for (j in start.j :(n-1))
		{
		v[n+1] <- v[n+1]-theta[n,(n-j)]*theta[n,(n-j)]*v[j+1]
		}
	if (n>2) 
		eps <- abs(v[n+1]/v[n-1]-1)
	}
if (n==a.niter)
	{
	msg <- paste("WARNING: Univ. Innov. Algor. stopped without achieving convergence after n=",n,"with eps=",eps)
	#warning(msg)
	}
if (n<a.niter)
	{
	msg <- paste("Univ. Innov. Algor. achieved convergence after n=",n,"with eps=",eps)
	}

my.polynom$coefs <- rep(0,q+1)
if (is.a.zero.polynom==0)
	{
	std.deviation<-sqrt(v[n])
	my.polynom$coefs[1] <- std.deviation
	for (j in 2:(q+1))
		my.polynom$coefs[j] <- std.deviation*theta[n,j-1]
	my.polynom$v <- v
	#my.polynom$theta <- theta
	}
my.polynom$msg<-msg
#cat("finished Univ. Innov. algorithm","\n")
my.polynom 
}

#****************************************************************
calculate.fourier.coefficients <- function(a.polynom,a.var)
#********************************
{
fourier.coefs <- rep(0,2*a.var$pmax.init+1)
for (k in 1:(a.polynom$order+1))    
	{
	fourier.coefs[k] <- 0
	for (j in 1:(a.polynom$order+2-k))    
		{
		fourier.coefs[k] <- fourier.coefs[k] + a.polynom$coefs[j]*a.polynom$coefs[j+k-1]
		}
	}
a.polynom$fourier.coefs <- fourier.coefs
a.polynom
}
#****************************************************************
calc.values.new <- function(polynom,a.grid)
#********************************
{
polynom$inv.values <- list()
polynom$inv.values$freq <- seq(from=0,to=0.5,length=a.grid)
lamdas <- polynom$inv.values$freq*pi*2
b.lamdas <- complex(length.out=a.grid,modulus=1,argument=-lamdas)
values <- complex(length.out=a.grid,modulus=0,argument=lamdas)
p <- polynom$order
coefs <- polynom$coefs
for (i in 1:(p+1))
	{
	c.lamdas <- b.lamdas^(i-1)
	values <- values + coefs[i]*c.lamdas
	}
polynom$inv.values$spec <- -log(Mod(values)^2)
polynom$values$spec <- (Mod(values)^2)
polynom$values$freq<-polynom$inv.values$freq
polynom
}

#****************************************************************
check.orders.det.cross <- function(a.var)
#********************************
{
if (a.var$det$order>2*a.var$pmax.init)
	{
	stop("Order of det is too high. order(det)=",a.var$det$order,", while order.max=",2*a.var$pmax.init)
	}
if (a.var$cross$order>2*a.var$pmax.init)
	{
	stop("Order of cross is too high. order(cross)=",a.var$cross$order,", while order.max=",2*a.var$pmax.init)
	}
}


#****************************************************************
check.orders.eta.ksi.zeta <- function(a.var)
#********************************
{
if (a.var$eta.1$order>2*a.var$pmax.init)
	{
	stop("Order of eta.1 is too high. order(eta.1)=",a.var$eta.1$order,", while order.max=",2*a.var$pmax.init)
	}
if (a.var$eta.2$order>2*a.var$pmax.init)
	{
	stop("Order of eta.2 is too high. order(eta.2)=",a.var$eta.2$order,", while order.max=",2*a.var$pmax.init)
	}
if (a.var$ksi.1$order>2*a.var$pmax.init)
	{
	stop("Order of ksi.1 is too high. order(ksi.1)=",a.var$ksi.1$order,", while order.max=",2*a.var$pmax.init)
	}
if (a.var$ksi.2$order>2*a.var$pmax.init)
	{
	stop("Order of ksi.2 is too high. order(ksi.2)=",a.var$ksi.2$order,", while order.max=",2*a.var$pmax.init)
	}
if (a.var$ksi.c$order>2*a.var$pmax.init)
	{
	stop("Order of ksi.c is too high. order(ksi.c)=",a.var$ksi.c$order,", while order.max=",2*a.var$pmax.init)
	}
if (a.var$zeta$order>2*a.var$pmax.init)
	{
	stop("Order of zeta is too high. order(zeta)=",a.var$zeta$order,", while order.max=",2*a.var$pmax.init)
	}
}

#****************************************************************
check.VAR.basic.relation<- function(a.var)
#********************************
{
par.def<-par(no.readonly = TRUE)
on.exit(par(par.def))

par(mfrow=c(1,2))
freq.labels <- c("0","0.25 pi","0.5 pi","0.75 pi","pi")
freq.labels.where <-seq(from=0,to=1,length=5)*a.var$grid

v.1 <- -a.var$det$inv.values$spec
v.2 <- -a.var$cross$inv.values$spec
plot.ts(as.ts(cbind(a.var$chi.1.prod.2$inv.values$spec,-log(exp(v.1)+exp(v.2)))),
		type="l", main="1/(chi1^2*chi2^2),1/(det^2+cross^2)" , plot.type="single",xaxt="n",ylab='log spectra',xlab='frequency',
 		col=c(4,2),lwd=c(2,1) )
axis(1,at=freq.labels.where,labels=freq.labels)

plot.ts(as.ts(cbind(a.var$chi.1.prod.2$inv.values$spec+log(exp(v.1)+exp(v.2)))), xaxt="n",ylab='log spectra',xlab='frequency',
		plot.type="single",type="l", main="1/(chi1^2*chi2^2-[det^2+cross^2])", col=c(4,2),lwd=c(2,1) )
axis(1,at=freq.labels.where,labels=freq.labels)
}

