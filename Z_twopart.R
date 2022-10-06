### set workspace
setwd("your workspace")

library(INLA)

low.resol <- TRUE

### load the data:
load("data_train.RData")

newd <- data_train_DF

na2 <- is.na(newd[c('CNT', 'BA')])

iszero2 <- newd[c('CNT', 'BA')]==0

### where BA is zero and CNT is NA, set to zero
newd$cnt <- newd$CNT
newd$cnt[na2[,1]&iszero2[,2]] <- 0
### where CNT is zero and BA is NA, set to zero
newd$ba <- newd$BA
newd$ba[na2[,2]&iszero2[,1]] <- 0

na2n <- is.na(newd[c('cnt', 'ba')])
is0n <- newd[c('cnt', 'ba')]==0

newd$z <- NA # define latent binary variable
bernoulli.var <- function(cnt,ba){
  if (is.na(cnt) & is.na(ba)) {
    Z = NA
  } else {
    temp = min(cnt, ba, na.rm=TRUE)
    Z = ifelse(temp>0, 1, 0)
  }
  return(Z)
}
newd$z <- mapply(bernoulli.var, newd$cnt, newd$ba)

### number of observations and covariate names
(nd <- nrow(newd))
(nx <- length(xnams <- colnames(newd)[c(6:37)]))
### mean and SD vectors
x.mean <- colMeans(newd[xnams])
x.sd <- sapply(newd[xnams], sd)
### standardize covariates
xx.std <- sweep(sweep(newd[xnams], 2, x.mean), 2, x.sd, '/')


### unique spatial locations
loc.area <- unique(newd[c("lon", "lat", 'area')])
locs <- as.matrix(loc.area[,1:2])
(n.locs <- nrow(locs))

### create an area index
newd$id.area <- pmatch(
  paste(newd$lon, newd$lat),
  paste(locs[,1], locs[,2]),
  duplicates.ok=TRUE)

### compute sd (over time) by area
locs.sd <- aggregate(newd[c('z')],
                     newd['id.area'], sd, na.rm=TRUE)


### build a mesh for spatial model
mesh1 <- inla.mesh.2d(locs, max.edge=3, offset=5, cutoff=0.3)
mesh1$n

### build a mesh for spatio-temporal model
mesh2 <- inla.mesh.2d(locs, max.edge=5, offset=5, cutoff=0.7)

if(low.resol) {
  mesh2 <- inla.mesh.2d(locs, max.edge=5, offset=5, cutoff=1)
}
mesh2$n

### IDEA
### non-stationary spatial model for each outcome
### the variance depends on the empirical variance
### which is computed over time for each location

### build a SPDE model with variance depending on the SD of CNT and BA
str(mesh1$idx) ### locs as part of the mesh nodes
str(mesh1id.unique <- unique(mesh1$idx$loc))
str(locs.id.unique <- which(!duplicated(mesh1$idx$loc)))

### base parameters
nu <- 1
alpha <- nu + 2/2
lkappa0 <- log(8*nu)/2
ltau0 <- (lgamma(nu)-lgamma(alpha) -log(4*pi))/2 -lkappa0

mesh1.sd.z <- rep(0, mesh1$n)
mesh1.sd.z[mesh1$idx$loc] <- log(1+locs.sd$z)


spde.z <- inla.spde2.matern(
  mesh=mesh1, alpha=alpha,
  B.tau=cbind(ltau0, nu, cbind(-1, -mesh1.sd.z)), 
  B.kappa=cbind(lkappa0, -1, cbind(0.0, 0.0*mesh1.sd.z)))


if(FALSE) { ### checking the model
  
  qq1 <- inla.spde2.precision( ### build precision 
    spde.z, theta=c(1, 0.5, 0.3)) ### for a given theta
  vv1 <- inla.qinv(qq1)
  summary(diag(vv1))
  
  par(mfrow=c(2,1), mar=c(0,0,1,0))
  plot(locs, asp=1, axes=FALSE, pch=19, cex=log(1.2+locs.sd$z)/4)
  plot(mesh1$loc[,1:2], asp=1, axes=FALSE, pch=19, cex=0.1+diag(vv1)/2)
  
}

### build projection matrices
A.s <- inla.spde.make.A(mesh1, cbind(newd$lon, newd$lat))
stopifnot(all.equal(sum(A.s),nrow(newd)))

### year as replications and month as group
(nrepl <- max(newd$id.year <- 1+newd$year-min(newd$year)))
(ngroup <- max(newd$id.month <- 1+newd$month-min(newd$month)))
A.st <- inla.spde.make.A(
  mesh=mesh2,
  loc=cbind(newd$lon, newd$lat),
  group=newd$id.month, repl=newd$id.year,
  n.group=ngroup, n.repl=nrepl)
stopifnot(all.equal(sum(A.st),nrow(newd)))

### spatio-temporal index set and replicate (each outcome + copy)
st1.idx <- inla.spde.make.index(
  name='st1', n.spde=mesh2$n,
  n.group=ngroup, n.repl=nrepl)


### z stack
(nd <- nrow(newd))
table(ipred.z <- is.na(newd$z))
### stack the observed data on Z
z.obs.stack <- inla.stack(
  tag='z.obs',
  data=list(Y=newd$z[!ipred.z], 
            ilink=rep(1, sum(!ipred.z))),
  effects=list(data.frame(z.b0=1, z=xx.std[!ipred.z,]), 
               list(z.s=1:spde.z$n.spde),        
               st1.idx),
  A=list(1, A.s[!ipred.z,], A.st[!ipred.z, ]))

### prediction stack for Z
z.pred.stack <- inla.stack(
  tag='z.pred',
  data=list(Y=rep(NA, sum(ipred.z)),
            ilink=rep(1, sum(ipred.z))),
  effects=list(data.frame(z.b0=1, z=xx.std[ipred.z,]), 
               list(z.s=1:spde.z$n.spde),        
               st1.idx),
  A=list(1, A.s[ipred.z,], A.st[ipred.z, ]))


### model formula
zterms <- paste(paste0('z.', c('b0', xnams)), collapse='+')
f0 <- update(Y ~ 0, paste0('.~.+', zterms))
f0

spde.st0 <- inla.spde2.pcmatern(
  mesh=mesh2, alpha=2,
  prior.range=c(1, 0.1),
  prior.sigma=c(0.5, 0.1))

pccor1 <- list(rho=list(prior='pc.cor1', param=c(0, 0.95)))
cgroup <- list(model='ar1', hyper=pccor1)

ff <- update(
  f0, .~.+ f(z.s, model=spde.z) +
    f(st1, model=spde.st0, replicate=st1.repl,
      group=st1.group, control.group=cgroup))
ff

inla.setOption(
  pardiso.license='~/.pardiso.lic',
  smtp='pardiso',
  inla.mode='experimental',
  num.threads='6:-8') ### to run on a 128 threads machine

gc(reset=TRUE)

result.z <- inla(
  ff, family=c('binomial'),
  data=inla.stack.data(inla.stack(z.obs.stack,z.pred.stack)),
  control.predictor=list(
    A=inla.stack.A(inla.stack(z.obs.stack,z.pred.stack)),
    compute=TRUE, link=ilink),
  verbose=TRUE)

save('result.z',
     file='spdetwos_z_4_new.RData',
     compress='xz')
