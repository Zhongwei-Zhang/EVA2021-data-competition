### set workspace
setwd("your workspace")

library(INLA)
#library(maps)

low.resol <- TRUE

### load the data:
load("data_train.RData")

head(data_train_DF,1)
newd <- data_train_DF

na2 <- is.na(newd[c('CNT', 'BA')])
table(cnt=na2[,1], ba=na2[,2])

iszero2 <- newd[c('CNT', 'BA')]==0
table(cnt=iszero2[,1], ba=iszero2[,2])

table(cnt=na2[,1], ba=na2[,2],
      any0=iszero2[,1] | iszero2[,2])
table(cnt=na2[,1], ba=na2[,2],
      two0=iszero2[,1]&iszero2[,2])

### where BA is zero and CNT is NA, set to zero
newd$cnt <- newd$CNT
newd$cnt[na2[,1]&iszero2[,2]] <- 0
### where CNT is zero and BA is NA, set to zero
newd$ba <- newd$BA
newd$ba[na2[,2]&iszero2[,1]] <- 0

na2n <- is.na(newd[c('cnt', 'ba')])
is0n <- newd[c('cnt', 'ba')]==0

table(cnt=na2n[,1], ba=na2n[,2])
table(cnt=is0n[,1], ba=is0n[,2])
table(is0n[,1]&is0n[,2])
prop.table(table(is0n[,1]&is0n[,2]))*100

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
locs.sd <- aggregate(newd[c('cnt', 'ba')],
                     newd['id.area'], sd, na.rm=TRUE)

#par(mfrow=c(2,1), mar=c(0,0,0,0))
#plot(locs, cex=log(1.2+locs.sd$cnt)/4, pch=19, axes=F)
#plot(locs, cex=log(1.2+locs.sd$ba)/15, pch=19, axes=F)

### build a mesh for spatial model
mesh1 <- inla.mesh.2d(locs, max.edge=3, offset=5, cutoff=0.3)
mesh1$n

### build a mesh for spatio-temporal model
mesh2 <- inla.mesh.2d(locs, max.edge=5, offset=5, cutoff=0.7)
mesh2$n

if(low.resol) {
  mesh2 <- inla.mesh.2d(locs, max.edge=5, offset=5, cutoff=1)
}

print(c(nmesh1=mesh1$n, nmesh2=mesh2$n))

### build a SPDE model with variance depending on the SD of CNT and BA
str(mesh1$idx) ### locs as part of the mesh nodes
str(mesh1id.unique <- unique(mesh1$idx$loc))
str(locs.id.unique <- which(!duplicated(mesh1$idx$loc)))

### base parameters
nu <- 1
alpha <- nu + 2/2
lkappa0 <- log(8*nu)/2
ltau0 <- (lgamma(nu)-lgamma(alpha) -log(4*pi))/2 -lkappa0

mesh1.sd.cnt <- rep(0, mesh1$n)
mesh1.sd.cnt[mesh1$idx$loc] <- log(1+locs.sd$cnt)
mesh1.sd.ba <- rep(0, mesh1$n)
mesh1.sd.ba[mesh1$idx$loc] <- log(1+locs.sd$ba)

spde.cnt <- inla.spde2.matern(
  mesh=mesh1, alpha=alpha,
  B.tau=cbind(ltau0, nu, cbind(-1, -mesh1.sd.cnt)), 
  B.kappa=cbind(lkappa0, -1, cbind(0.0, 0.0*mesh1.sd.cnt)))

spde.ba <- inla.spde2.matern(
  mesh=mesh1, alpha=alpha,
  B.tau=cbind(ltau0, nu, cbind(-1, -mesh1.sd.ba)),
  B.kappa=cbind(lkappa0, -1, cbind(0.0, 0.0*mesh1.sd.ba)))

if(FALSE) { ### checking the model
  
  qq1 <- inla.spde2.precision( ### build precision 
    spde.cnt, theta=c(1, 0.5, 0.3)) ### for a given theta
  vv1 <- inla.qinv(qq1)
  summary(diag(vv1))
  
  par(mfrow=c(2,1), mar=c(0,0,1,0))
  plot(locs, asp=1, axes=FALSE, pch=19, cex=log(1.2+locs.sd$cnt)/4)
  plot(mesh1$loc[,1:2], asp=1, axes=FALSE, pch=19, cex=0.1+diag(vv1)/2)
  
}

### build projection matrices
A.s <- inla.spde.make.A(mesh1, cbind(newd$lon, newd$lat))
stopifnot(all.equal(sum(A.s),nrow(newd)))

### year as replication and month as group
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
st2.idx <- inla.spde.make.index(
  name='st2', n.spde=mesh2$n,
  n.group=ngroup, n.repl=nrepl)
stc.idx <- inla.spde.make.index(
  name='stc', n.spde=mesh2$n,
  n.group=ngroup, n.repl=nrepl)

sapply(st1.idx, range)
sapply(st2.idx, range)
sapply(stc.idx, range)

### cnt stack
(nd <- nrow(newd))
ipos.cnt <- newd$CNT>0
iipos.cnt <- which(ipos.cnt)
ipred.cnt <- is.na(newd$cnt)
### stack the observed data on CNT
cnt.obs.stack <- inla.stack(
  tag='cnt.obs',
  data=list(Y=cbind(newd$cnt[iipos.cnt]-1, NA),  #CNT-1
            ilink=rep(1, length(iipos.cnt))),
  effects=list(data.frame(cnt.b0=1, cnt=xx.std[iipos.cnt,]), 
               list(cnt.s=1:spde.cnt$n.spde),        
               st1.idx, stc.idx),
  A=list(1, A.s[iipos.cnt,], A.st[iipos.cnt, ], A.st[iipos.cnt, ]))
### prediction stack for CNT
cnt.pred.stack <- inla.stack(
  tag='cnt.pred',
  data=list(Y=cbind(rep(NA, sum(ipred.cnt)), NA), 
            ilink=rep(1, sum(ipred.cnt))),
  effects=list(data.frame(cnt.b0=1, cnt=xx.std[ipred.cnt,]), 
               list(cnt.s=1:spde.cnt$n.spde),        
               st1.idx, stc.idx),
  A=list(1, A.s[ipred.cnt,], A.st[ipred.cnt, ], A.st[ipred.cnt, ]))


### consider min positive BA to treat zero
ipos.ba <- newd$BA>0
iipos.ba <- which(ipos.ba)
ipred.ba <- is.na(newd$ba)
### stack the observed data on BA
ba.obs.stack <- inla.stack(
  tag='ba.obs',
  data=list(Y=cbind(NA, log(newd$ba[iipos.ba])),  #log(BA)
            ilink=rep(2, length(iipos.ba))),
  effects=list(data.frame(ba.b0=1, ba=xx.std[iipos.ba,]), 
               ba.s=1:spde.ba$n.spde,        
               st2.idx),
  A=list(1, A.s[iipos.ba,], A.st[iipos.ba, ]))
### prediction stack for BA
ba.pred.stack <- inla.stack(
  tag='ba.pred',
  data=list(Y=cbind(rep(NA, sum(ipred.ba)), NA), 
            ilink=rep(2, sum(ipred.ba))),
  effects=list(data.frame(ba.b0=1, ba=xx.std[ipred.ba,]), 
               ba.s=1:spde.ba$n.spde,        
               st2.idx),
  A=list(1, A.s[ipred.ba,], A.st[ipred.ba, ]))


### model formula
cntterms <- paste(paste0('cnt.', c('b0', xnams)), collapse='+')
baterms <- paste(paste0('ba.', c('b0', xnams)), collapse='+')
f0 <- update(Y ~ 0, paste0('.~.+', cntterms,'+', baterms))
f0

spde.st0 <- inla.spde2.pcmatern(
  mesh=mesh2, alpha=2,
  prior.range=c(1, 0.1),
  prior.sigma=c(0.5, 0.1))

pccor1 <- list(rho=list(prior='pc.cor1', param=c(0, 0.95)))
cgroup <- list(model='ar1', hyper=pccor1)

ff <- update(
  f0, .~.+ f(cnt.s, model=spde.cnt) +
    f(ba.s, model=spde.ba) +
    f(st1, model=spde.st0, replicate=st1.repl,
      group=st1.group, control.group=cgroup) +
    f(st2, model=spde.st0, replicate=st2.repl,
      group=st2.group, control.group=cgroup) +
    f(stc, copy='st2', group=stc.group, replicate=stc.repl, fixed=FALSE))
ff

inla.setOption(
  pardiso.license='~/.pardiso.lic',
  smtp='pardiso',
  inla.mode='experimental',
  num.threads='8:-16') ### to run on a 48 threads machine

gc(reset=TRUE)

result.cntba <- inla(
  ff, family=c('poisson', 'gaussian'),
  data=inla.stack.data(inla.stack(cnt.obs.stack,cnt.pred.stack, ba.obs.stack,ba.pred.stack)),
  control.predictor=list(
    A=inla.stack.A(inla.stack(cnt.obs.stack,cnt.pred.stack, ba.obs.stack,ba.pred.stack)),
    compute=TRUE, link=ilink),
  verbose=TRUE)

save('result.cntba',
     file='spdetwos_cntba.RData',
     compress='xz')
