

delta.estimate= function(xone,xzero, deltaone, deltazero, t, std= FALSE, conf.int = FALSE, weight.perturb = NULL) {
	delta.f = delta.estimate.RMST(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, weight = weight.perturb, delta.only=F)
	delta = delta.f$delta
	rmst.1 = delta.f$s1
	rmst.0 = delta.f$s0
	if(std | conf.int)	{
		n1 = length(xone)
		n0 = length(xzero)
		if(is.null(weight.perturb)) {weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)} 
		delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, delta.only = T)
		delta.p.vec = unlist(delta.p.vec)
		sd.delta = sd(delta.p.vec)
		mad.delta = mad(delta.p.vec)
	}
	if(conf.int){
		conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	}
	if(!std & !conf.int) {return(list("delta" = delta, "rmst.1" = rmst.1, "rmst.0" = rmst.0))}
	if(std & !conf.int) {return(list("delta" = delta, "rmst.1" = rmst.1, "rmst.0" = rmst.0, "delta.sd" = sd.delta, "delta.mad" = mad.delta))}
	if(conf.int) {return(list("delta" = delta, "rmst.1" = rmst.1, "rmst.0" = rmst.0, "delta.sd" = sd.delta, "delta.mad" = mad.delta, "conf.int" = conf.quantile.delta))}
}

delta.estimate.RMST= function(xone,xzero, deltaone, deltazero, t, weight = NULL, delta.only=F) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	censor1.t = censor.weight(xone, deltaone, t, weight = weight[1:length(xone)])
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight[(1+length(xone)):(length(xone)+length(xzero))])
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight[1:length(xone)])
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight[(1+length(xone)):(length(xone)+length(xzero))])
		M.1 = rep(0,length(xone))
		M.0 = rep(0,length(xzero))
		M.1[xone>t] = t/censor1.t
		M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 	
		M.0[xzero>t] = t/censor0.t
		M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
		delta = 1/(sum(weight[1:length(xone)])) * sum(M.1*weight[1:length(xone)]) - 1/(sum(weight[(1+length(xone)):(length(xone)+length(xzero))])) * sum(M.0*weight[(1+length(xone)):(length(xone)+length(xzero))])
		s1 = 1/(sum(weight[1:length(xone)])) * sum(M.1*weight[1:length(xone)])
		s0 = 1/(sum(weight[(1+length(xone)):(length(xone)+length(xzero))])) * sum(M.0*weight[(1+length(xone)):(length(xone)+length(xzero))])
	if(delta.only) { return(list("delta" = delta)) }
	if(!delta.only) { return(list("delta" = delta, "s1" = s1, "s0" = s0)) } 

}

R.q.event = function(xone,xzero, deltaone, deltazero, sone, szero, t, landmark, number = 40, transform = FALSE, extrapolate = TRUE, std = FALSE, conf.int = FALSE, weight.perturb = NULL, type = "np") {
	if(!(type %in% c("np","semi"))) {warning("Warning: Invalid type, default `np' for nonparametric estimator being used", call. = FALSE); type = "np"}
	warn.te = FALSE
	warn.support = FALSE
	n1 = length(xone)
	n0 = length(xzero)
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)
	}
	delta = delta.estimate.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t)
	delta.estimate = delta$delta
	delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, delta.only = T)
	delta.p.vec = unlist(delta.p.vec)
	sd.delta  = sd(delta.p.vec)
	mad.delta = mad(delta.p.vec)
	conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	if(0>conf.quantile.delta[1] & 0< conf.quantile.delta[2]) {warning("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting", call. = FALSE)
		warn.te = TRUE}
	if(delta.estimate < 0) {warning("Warning: it looks like you need to switch the treatment groups", call. = FALSE)}
	range.1 = range(pmin(sone,landmark))
	range.0 = range(pmin(szero,landmark))
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if(range.ind & !extrapolate & !transform) {
		warning("Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation", call. = FALSE)
		warn.support = TRUE
	}
	
	if(type == "np"){
		delta.q.estimate = delta.q.event.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	if(type == "semi"){
		delta.q.estimate = delta.q.event.semi.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	if(std | conf.int){
		if(type == "np"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number, deltaslist=FALSE, warn.extrapolate = FALSE))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
		if(type == "semi"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.semi.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, deltaslist=FALSE, number = number))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
	}
	
	if(conf.int) {	

		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.quantile.delta.q = quantile(delta.q.p.vec, 0.025)
		conf.u.quantile.delta.q = quantile(delta.q.p.vec, 0.975)
		
		conf.l.quantile.R.q = quantile(R.p, 0.025)
		conf.u.quantile.R.q = quantile(R.p, 0.975)
}
	if(!std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q))}
	if(std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p)))}
	if(conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p), "conf.int.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.delta.q" = as.vector(c(conf.l.quantile.delta.q, conf.u.quantile.delta.q)),
"conf.int.R.q" = as.vector(c(conf.l.quantile.R.q, conf.u.quantile.R.q))))
}	
}




delta.q.event.RMST = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight = NULL, landmark=landmark, deltaslist = TRUE, transform = FALSE, extrapolate = TRUE, number, warn.extrapolate = TRUE) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
	
	#first term
	t.vec = c(landmark, landmark+c(1:(number-1))*(t-landmark)/number, t)
	#number surrogate events in control group before landmark
	number.s.o = length(szero[xzero>landmark & szero < landmark])
	phi.a.tst0 = matrix(nrow = number.s.o, ncol = number+1)
	warn.temp = FALSE
	for(i in 1:(number+1)) {
	temp.phi = pred.smooth.surv(xone.f=xone[xone>landmark & sone < landmark], deltaone.f = deltaone[xone>landmark & sone <  landmark], sone.f=log(sone[xone>landmark & sone <  landmark]), szero.one = log(szero[xzero>landmark & szero < landmark]), myt=t.vec[i], weight = weight.group1[xone>landmark & sone <  landmark], transform = transform, extrapolate = extrapolate)
	phi.a.tst0[,i] =  temp.phi$Phat.ss
	if(temp.phi$warn.flag == 1) {warn.temp = TRUE}
	}
	if(warn.temp == TRUE & warn.extrapolate == TRUE) {if(warn.extrapolate) {warning("Note: Values were extrapolated.", call. = FALSE)}}
	phi.int = vector(length = number.s.o)
	for(j in 1:number.s.o) {
		phi.int[j] = landmark + (t-landmark)/number*(phi.a.tst0[j,1]/2 + sum(phi.a.tst0[j,2:number]) + phi.a.tst0[j,number+1]/2)
	}
	first.term = sum(weight.group0[xzero>landmark & szero <landmark]*phi.int)/(sum(weight.group0)*censor0.landmark) 
	
	
	#second term
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight.group1)
	M.1 = rep(0,length(xone))
	M.1[xone>t] = t/censor1.t
	M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 
	denom = sum(1*(sone > landmark & xone > landmark)*weight.group1)	
	psi.a.tt0 = sum(censor1.landmark*M.1*weight.group1*(sone > landmark & xone > landmark)/(denom))

	second.term = sum(weight.group0*(szero > landmark & xzero > landmark)*psi.a.tt0)/(sum(weight.group0)*censor0.landmark) 
	
	#third.term
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight.group0)
	M.0 = rep(0,length(xzero))
	M.0[xzero>t] = t/censor0.t
	M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
	denom = sum(1*( xzero > landmark)*weight.group0)	
	nu.term = sum(censor0.landmark*M.0*weight.group0*(xzero > landmark)/(denom))
	third.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term

	
	delta.q = first.term+second.term - third.term
	
	if(deltaslist) {return(list("delta.q" = delta.q))}
	if(!deltaslist) {return(c(first.term, second.term, third.term, delta.q))}


}

delta.q.event.semi.RMST = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight = NULL, landmark=landmark, deltaslist = TRUE, number) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
	
	#first term
	s.predictor = sone[xone>landmark & sone <  landmark]
	x.adjust = xone[xone>landmark & sone < landmark] - landmark
	sum.model = coxph(Surv(x.adjust, deltaone[xone>landmark & sone <  landmark])~s.predictor, weights = weight.group1[xone>landmark & sone <  landmark])
	s.new = as.data.frame(szero[xzero>landmark & szero < landmark])
	names(s.new) = "s.predictor"
	predict.score = predict(sum.model, type = "risk", newdata = s.new)
	baseline.hazard = basehaz(sum.model, centered = TRUE); 
	#yes, we really do want centered = TRUE
	
	t.vec = c(landmark, landmark+c(1:(number-1))*(t-landmark)/number, t)
	#number surrogate events in control group before landmark
	number.s.o = length(szero[xzero>landmark & szero < landmark])
	phi.a.tst0 = matrix(nrow = number.s.o, ncol = number+1)
	for(i in 1:(number+1)) {
		gap.time = t.vec[i]-landmark
		baseline.t <- approx(baseline.hazard$time,baseline.hazard$hazard,gap.time, rule = 2)$y
		phi.a.tst0[,i] =  exp(-baseline.t*predict.score)
	}
	phi.int = vector(length = number.s.o)
	for(j in 1:number.s.o) {
		phi.int[j] = landmark + (t-landmark)/number*(phi.a.tst0[j,1]/2 + sum(phi.a.tst0[j,2:number]) + phi.a.tst0[j,number+1]/2)
	}
	first.term = sum(weight.group0[xzero>landmark & szero <landmark]*phi.int)/(sum(weight.group0)*censor0.landmark) 
	
	
	#second term
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight.group1)
	M.1 = rep(0,length(xone))
	M.1[xone>t] = t/censor1.t
	M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 
	denom = sum(1*(sone > landmark & xone > landmark)*weight.group1)	
	psi.a.tt0 = sum(censor1.landmark*M.1*weight.group1*(sone > landmark & xone > landmark)/(denom))

	second.term = sum(weight.group0*(szero > landmark & xzero > landmark)*psi.a.tt0)/(sum(weight.group0)*censor0.landmark) 
	
	#third.term
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight.group0)
	M.0 = rep(0,length(xzero))
	M.0[xzero>t] = t/censor0.t
	M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
	denom = sum(1*( xzero > landmark)*weight.group0)	
	nu.term = sum(censor0.landmark*M.0*weight.group0*(xzero > landmark)/(denom))
	third.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term

	
	delta.q = first.term+second.term - third.term
	
	if(deltaslist) {return(list("delta.q" = delta.q))}
	if(!deltaslist) {return(c(first.term, second.term, third.term, delta.q))}


}

censor.weight = function(data.x, data.delta, t, weight=NULL) {
	if(is.null(weight)) {weight = rep(1,length(data.x))}
	S.KM = survfit(Surv(data.x,1-data.delta)~1, weights = weight)
	S.t.KM = approx(S.KM$time,S.KM$surv,t, rule=2)$y
	return(S.t.KM)
}

pred.smooth.surv <- function(xone.f, deltaone.f, sone.f, szero.one, myt, bw = NULL, weight, transform, extrapolate = T)
  { 
  	warn.flag = 0
    if(transform){	 	
    	mean.o= mean(c(sone.f, szero.one))
  		sd.o = sd(c(szero.one, sone.f))
    	sone.f.new = pnorm((sone.f - mean.o)/sd.o)
    	szero.one.new = pnorm((szero.one - mean.o)/sd.o)
		sone.f = sone.f.new
		szero.one = szero.one.new
	} 

	if(is.null(bw))
      {
        bwini = bw.nrd(sone.f)
        n.s = length(sone.f)
        bw <- bwini/(n.s^0.11)
      }      
    kerni.ss = Kern.FUN(zz=sone.f,zi=szero.one,bw)           
    tmpind = (xone.f<=myt)&(deltaone.f==1); tj = xone.f[tmpind]; 
    kerni.1 = t(weight*t(kerni.ss))
    pihamyt0.tj.ss = helper.si(tj, "<=", xone.f, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
	#dLamhat.tj.ss[is.na(dLamhat.tj.ss)] = 0
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	warn.flag = 1
    	c.mat = cbind(szero.one, Phat.ss)
    	for(o in 1:length(Phat.ss)) {
    		if(is.na(Phat.ss[o])){
    			distance = abs(szero.one - szero.one[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			Phat.ss[o] = new.est[1]   #in case there are multiple matches
    	}
    }
	}
    return(list("Phat.ss"=Phat.ss, "warn.flag" = warn.flag))
    }

VTM<-function(vc, dm){
	#takes vc and makes it the repeated row of a matrix, repeats it dm times
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    

Kern.FUN <- function(zz,zi,bw=NULL,kern0="gauss") ## returns an (n x nz) matrix ##
  { 
    if(is.null(bw))
      {
        bwini = bw.nrd(zz)
        n.s = length(zz)
        bw <- bwini/(n.s^0.11)
      } 
     out = (VTM(zz,length(zi))- zi)/bw
    switch(kern0,
            "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
            "gauss"= dnorm(out)/bw
           )
  }
  
  
  cumsum2 <- function(mydat)     #cumsum by row, col remains the same
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }

helper.si <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  }
 

R.t.estimate = function(xone,xzero, deltaone, deltazero, t, landmark, std = FALSE, conf.int = FALSE, weight.perturb = NULL) {
	warn.te = FALSE
	n1 = length(xone)
	n0 = length(xzero)
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)
	}
	delta = delta.estimate.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t)
	delta.estimate = delta$delta
	delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, delta.only = T)
	delta.p.vec = unlist(delta.p.vec)
	sd.delta  = sd(delta.p.vec)
	mad.delta = mad(delta.p.vec)
	conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	if(0>conf.quantile.delta[1] & 0< conf.quantile.delta[2]) {warning("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting", call. = FALSE)
		warn.te = TRUE}
	if(delta.estimate < 0) {warning("Warning: it looks like you need to switch the treatment groups", call. = FALSE)}
	delta.t.estimate = delta.t.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark)$delta.t
	R.t = 1-delta.t.estimate/delta.estimate	
	if(std | conf.int){
		delta.t.p.vec.temp = t(apply(weight.perturb, 2, delta.t.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark))
		delta.t.p.vec = unlist(delta.t.p.vec.temp)
		R.t.p = 1- delta.t.p.vec/delta.p.vec
	}
	if(conf.int) {	

		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.quantile.delta.t = quantile(delta.t.p.vec, 0.025)
		conf.u.quantile.delta.t = quantile(delta.t.p.vec, 0.975)
		
		conf.l.quantile.R.t = quantile(R.t.p, 0.025)
		conf.u.quantile.R.t = quantile(R.t.p, 0.975)
}
	if(!std & !conf.int) {return(list("delta" = delta.estimate, "delta.t" =delta.t.estimate, "R.t" = R.t))}
	if(std & !conf.int) {return(list("delta" = delta.estimate, "delta.t" =delta.t.estimate, "R.t" = R.t, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p)))}
	if(conf.int) {return(list("delta" = delta.estimate, "delta.t" =delta.t.estimate, "R.t" = R.t, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p), "conf.int.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.delta.t" = as.vector(c(conf.l.quantile.delta.t, conf.u.quantile.delta.t)),
"conf.int.R.t" = as.vector(c(conf.l.quantile.R.t, conf.u.quantile.R.t))))
}	
}

 
delta.t.RMST = function(xone,xzero, deltaone, deltazero, t, weight = NULL, landmark=landmark) {
	if(is.null(weight)) {weight = rep(1,length(xone)+length(xzero))}
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)	
	
	#a term
	censor1.times = censor.weight(xone, deltaone, xone, weight = weight.group1)
	M.1 = rep(0,length(xone))
	M.1[xone>t] = t/censor1.t
	M.1[xone<=t] = xone[xone<=t]*deltaone[xone<=t]/censor1.times[xone<=t] 
	denom = sum(1*( xone > landmark)*weight.group1)	
	nu.term.a = sum(censor1.landmark*M.1*weight.group1*(xone > landmark)/(denom))
	a.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term.a
		
	#b.term
	censor0.times = censor.weight(xzero, deltazero, xzero, weight = weight.group0)
	M.0 = rep(0,length(xzero))
	M.0[xzero>t] = t/censor0.t
	M.0[xzero<=t] = xzero[xzero<=t]*deltazero[xzero<=t]/censor0.times[xzero<=t] 
	denom = sum(1*( xzero > landmark)*weight.group0)	
	nu.term.b = sum(censor0.landmark*M.0*weight.group0*(xzero > landmark)/(denom))
	b.term = (sum(weight.group0*(xzero > landmark))/(sum(weight.group0)*censor0.landmark))*nu.term.b

	
	delta.t = a.term - b.term
	
	return(list("delta.t" = delta.t))


}


IV.event = function(xone,xzero, deltaone, deltazero, sone, szero, t, landmark, number = 40, transform = FALSE, extrapolate = TRUE, std = FALSE, conf.int = FALSE, weight.perturb = NULL, type = "np") {
	if(!(type %in% c("np","semi"))) {warning("Warning: Invalid type, default `np' for nonparametric estimator being used", call. = FALSE); type = "np"}
	warn.te = FALSE
	warn.support = FALSE
	n1 = length(xone)
	n0 = length(xzero)
	if(is.null(weight.perturb)){
		weight.perturb = matrix(rexp(500*(n1+n0), rate=1), ncol = 500)
	}
	delta = delta.estimate.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t)
	delta.estimate = delta$delta
	delta.p.vec = apply(weight.perturb, 2, delta.estimate.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, delta.only = T)
	delta.p.vec = unlist(delta.p.vec)
	sd.delta  = sd(delta.p.vec)
	mad.delta = mad(delta.p.vec)
	conf.quantile.delta = c(new.q(delta.p.vec, 0.025), new.q(delta.p.vec, 0.975))
	if(0>conf.quantile.delta[1] & 0< conf.quantile.delta[2]) {warning("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting", call. = FALSE)
		warn.te = TRUE}
	if(delta.estimate < 0) {warning("Warning: it looks like you need to switch the treatment groups", call. = FALSE)}
	range.1 = range(pmin(sone,landmark))
	range.0 = range(pmin(szero,landmark))
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if(range.ind & !extrapolate & !transform) {
		warning("Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation", call. = FALSE)
		warn.support = TRUE
	}
	
	if(type == "np"){
		delta.q.estimate = delta.q.event.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	if(type == "semi"){
		delta.q.estimate = delta.q.event.semi.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, number = number)$delta.q
		R.q = 1-delta.q.estimate/delta.estimate	
	}
	delta.t.estimate = delta.t.RMST(xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark)$delta.t
	R.t = 1-delta.t.estimate/delta.estimate	
	IV = R.q-R.t
	if(std | conf.int){
		delta.t.p.vec.temp = t(apply(weight.perturb, 2, delta.t.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, t = t, landmark=landmark))
		delta.t.p.vec = unlist(delta.t.p.vec.temp)
		R.t.p = 1- delta.t.p.vec/delta.p.vec
		if(type == "np"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, extrapolate = extrapolate, transform = transform, number = number, deltaslist=FALSE, warn.extrapolate = FALSE))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
		if(type == "semi"){
			delta.q.p.vec.temp = t(apply(weight.perturb, 2, delta.q.event.semi.RMST, xone = xone, xzero = xzero, deltaone = deltaone, deltazero = deltazero, sone = sone, szero = szero, t = t, landmark=landmark, deltaslist=FALSE, number = number))
			delta.q.p.vec = delta.q.p.vec.temp[,4]
			R.p = 1-delta.q.p.vec/delta.p.vec
		}
		IV.p = R.p-R.t.p
	}
	
	if(conf.int) {	

		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.quantile.delta.q = quantile(delta.q.p.vec, 0.025)
		conf.u.quantile.delta.q = quantile(delta.q.p.vec, 0.975)
		
		conf.l.quantile.R.q = quantile(R.p, 0.025)
		conf.u.quantile.R.q = quantile(R.p, 0.975)
		
		conf.l.quantile.delta.t = quantile(delta.t.p.vec, 0.025)
		conf.u.quantile.delta.t = quantile(delta.t.p.vec, 0.975)
		
		conf.l.quantile.R.t = quantile(R.t.p, 0.025)
		conf.u.quantile.R.t = quantile(R.t.p, 0.975)
		
		conf.l.quantile.IV = quantile(IV.p, 0.025)
		conf.u.quantile.IV = quantile(IV.p, 0.975)
}
	if(!std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q, "delta.t" =delta.t.estimate, "R.t" = R.t, "IV" = IV))}
	if(std & !conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q,  "delta.t" =delta.t.estimate, "R.t" = R.t, "IV" = IV, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p), "IV.sd" = sd(IV.p), "IV.mad" = mad(IV.p)))}
	if(conf.int) {return(list("delta" = delta.estimate, "delta.q" =delta.q.estimate, "R.q" = R.q,  "delta.t" =delta.t.estimate, "R.t" = R.t, "IV" = IV, "delta.sd" = sd(delta.p.vec),"delta.mad" = mad(delta.p.vec), "delta.q.sd" = sd(delta.q.p.vec),  "delta.q.mad" = mad(delta.q.p.vec), "R.q.sd" = sd(R.p), "R.q.mad" = mad(R.p), "delta.t.sd" = sd(delta.t.p.vec),  "delta.t.mad" = mad(delta.t.p.vec), "R.t.sd" = sd(R.t.p), "R.t.mad" = mad(R.t.p), "IV.sd" = sd(IV.p), "IV.mad" = mad(IV.p), "conf.int.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.delta.q" = as.vector(c(conf.l.quantile.delta.q, conf.u.quantile.delta.q)),
"conf.int.R.q" = as.vector(c(conf.l.quantile.R.q, conf.u.quantile.R.q)), "conf.int.delta.t" = as.vector(c(conf.l.quantile.delta.t, conf.u.quantile.delta.t)),
"conf.int.R.t" = as.vector(c(conf.l.quantile.R.t, conf.u.quantile.R.t)), "conf.int.IV" = as.vector(c(conf.l.quantile.IV, conf.u.quantile.IV))))
}	
}


new.q = function(x, p) {
	if(sum(is.na(x))>0) {return(NA)}
	if(sum(is.na(x))==0) {return(quantile(x,p))}
}


resam<- function(vv,t,t.0,tt,data,data1,data2,indexindex){
  
  nn=200 #number grid points for numeric integration

  n1=nrow(data1)
  n2=nrow(data2)

  causal=rep(NA,length(tt)); causal2=rep(NA,length(tt))
  causals=rep(NA,length(tt)); causals2=rep(NA,length(tt))
  causalind=rep(NA,length(tt)); causalind2=rep(NA,length(tt))
  causalsind=rep(NA,length(tt)); causalsind2=rep(NA,length(tt))
  for (j in 1:length(tt)){
  t=tt[j] 
  ################ pte2 given data1  
  xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1;v=vv[indexindex]
  
  from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
  s=seq(from, to, by = step)
  
  #### optimal function 
  {
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v, t.0=t.0, t=t)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
  kern = Kern.FUN(zz=s,zi=sob,bw)
  
  f0.s.t0.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  f0.s.t.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t)
  f1.s.t0.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  f1.s.t.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t)
  temp=(sob>t.0); temp[is.na(temp)]=1
  p0.t0.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  p0.t.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t)
  p1.t0.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  p1.t.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
  
  mu0t=sum(as.numeric(v)*(aob==0)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==0)*wei.t)
  mu1t=sum(as.numeric(v)*(aob==1)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
  
  lambda=nu/de
  
  gs.hat=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
  g3.hat=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
  }
  
  #### pte
  xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2;v=vv[-indexindex]
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v, t.0=t.0, t=t)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  # causal=mu1t-mu0t
  causal.t0=sum(as.numeric(v)*(xob>t.0)*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0)-sum(as.numeric(v)*(xob>t.0)*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0)
  causal[j]=sum(as.numeric(v)*(xob>t)*aob*wei.t)/sum(as.numeric(v)*aob*wei.t)-sum(as.numeric(v)*(xob>t)*(1-aob)*wei.t)/sum(as.numeric(v)*(1-aob)*wei.t)
  tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
  temp=(sob>t.0); temp[is.na(temp)]=1
  causals[j]=sum(as.numeric(v)*aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum(as.numeric(v)*aob*wei.t0)+
    sum(as.numeric(v)*(xob>t.0)*temp*g3.hat*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0) -
    (sum(as.numeric(v)*(1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum(as.numeric(v)*(1-aob)*wei.t0)+
       sum(as.numeric(v)*(xob>t.0)*temp*g3.hat*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0))
  pte.es=causals[j]/causal[j]
  

  ################ pte1 given data2  
  xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2;v=vv[-indexindex]
  
  from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
  s=seq(from, to, by = step)
  
  #### optimal function 
  {
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v, t.0=t.0, t=t)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
  kern = Kern.FUN(zz=s,zi=sob,bw)
  
  f0.s.t0.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  f0.s.t.t0.hat=apply(as.numeric(v)*(aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==0)*wei.t)
  f1.s.t0.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  f1.s.t.t0.hat=apply(as.numeric(v)*(aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum(as.numeric(v)*(aob==1)*wei.t)
  temp=(sob>t.0); temp[is.na(temp)]=1
  p0.t0.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t0)
  p0.t.t0.hat=sum(as.numeric(v)*(aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==0)*wei.t)
  p1.t0.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t0)
  p1.t.t0.hat=sum(as.numeric(v)*(aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
  
  mu0t=sum(as.numeric(v)*(aob==0)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==0)*wei.t)
  mu1t=sum(as.numeric(v)*(aob==1)*(xob>t)*wei.t)/sum(as.numeric(v)*(aob==1)*wei.t)
  
  integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
  temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
  nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
  
  lambda=nu/de
  
  gs.hat2=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
  g3.hat2=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
  }
  
  #### pte
  xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1;v=vv[indexindex]
  temp=WEIGHT.p(xob,deltaob,aob,n=n,v, t.0=t.0, t=t)
  wei.t0=temp[,1];wei.t=temp[,2]
  
  # causal=mu1t-mu0t
  causal.t0=(causal.t0+sum(as.numeric(v)*(xob>t.0)*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0)-
               sum(as.numeric(v)*(xob>t.0)*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0))/2
  causal2[j]=sum(as.numeric(v)*(xob>t)*aob*wei.t)/sum(as.numeric(v)*aob*wei.t)-sum(as.numeric(v)*(xob>t)*(1-aob)*wei.t)/sum(as.numeric(v)*(1-aob)*wei.t)
  tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
  temp=(sob>t.0); temp[is.na(temp)]=1
  causals2[j]=sum(as.numeric(v)*aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum(as.numeric(v)*aob*wei.t0)+
    sum(as.numeric(v)*(xob>t.0)*temp*g3.hat2*aob*wei.t0)/sum(as.numeric(v)*aob*wei.t0) -
    (sum(as.numeric(v)*(1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum(as.numeric(v)*(1-aob)*wei.t0)+
       sum(as.numeric(v)*(xob>t.0)*temp*g3.hat2*(1-aob)*wei.t0)/sum(as.numeric(v)*(1-aob)*wei.t0))
  pte.es=(pte.es+causals2[j]/causal2[j])/2
  

  g3.hat=(g3.hat+g3.hat2)/2
  gs.hat=(gs.hat+gs.hat2)/2
  }
  
#### output
  out=c(pte.es,g3.hat,gs.hat)#,pteind.es,ptermst,ptermstind
  
}


WEIGHT<-function(xob,deltaob,aob,n, t.0,t){
  x=xob[aob==0]; delta=1-deltaob[aob==0]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk,1,sum)
  Lam=cumsum(1/risk.n)
  s0=exp(-Lam)
  sur0=data.frame("time"=xsort,"surv"=s0)
  x=xob[aob==1]; delta=1-deltaob[aob==1]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk,1,sum)
  Lam=cumsum(1/risk.n)
  s1=exp(-Lam)
  sur1=data.frame("time"=xsort,"surv"=s1)
  ind0=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur0$time))}))
  G0=sur0$surv[ind0]
  ind1=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur1$time))}))
  G1=sur1$surv[ind1]
  G=(1-aob)*G0+aob*G1
  G.t0=(1-aob)*G0[which.min(abs(t.0-xob))]+aob*G1[which.min(abs(t.0-xob))]
  G.t=(1-aob)*G0[which.min(abs(t-xob))]+aob*G1[which.min(abs(t-xob))]
  wei.t0=(xob<=t.0)*deltaob/G+(xob>t.0)/G.t0;   wei.t=(xob<=t)*deltaob/G+(xob>t)/G.t;   
  out=cbind(wei.t0,wei.t,G.t,G.t0,G1,G0,G)#,pte2
}


WEIGHT.p<-function(xob,deltaob,aob,n,v, t.0,t){
  x=xob[aob==0];delta=1-deltaob[aob==0]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk*VTM(v[aob==0],length(xsort)),1,sum)
  N=(VTM(x,length(xsort))<=xsort)*VTM(delta,length(xsort))
  dN=rbind( N[1,],N[-1,]-N[-length(xsort),] )
  nu=apply(VTM(v[aob==0],length(xsort))*dN, 1, sum)
  Lam=cumsum(nu/risk.n)
  s0=exp(-Lam)
  sur0=data.frame("time"=xsort,"surv"=s0)
  x=xob[aob==1];delta=1-deltaob[aob==1]
  xsort=sort(x[delta==1])
  risk=VTM(x,length(xsort))>=xsort
  risk.n=apply(risk*VTM(v[aob==1],length(xsort)),1,sum)
  N=(VTM(x,length(xsort))<=xsort)*VTM(delta,length(xsort))
  dN=rbind( N[1,],N[-1,]-N[-length(xsort),] )
  nu=apply(dN*VTM(v[aob==1],length(xsort)), 1, sum)
  Lam=cumsum(nu/risk.n)
  s1=exp(-Lam)
  sur1=data.frame("time"=xsort,"surv"=s1)
  ind0=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur0$time))}))
  G0=sur0$surv[ind0]
  ind1=c(sapply(1:n, function(kk){which.min(abs(xob[kk]-sur1$time))}))
  G1=sur1$surv[ind1]
  G=(1-aob)*G0+aob*G1
  G.t0=(1-aob)*G0[which.min(abs(t.0-xob))]+aob*G1[which.min(abs(t.0-xob))]
  G.t=(1-aob)*G0[which.min(abs(t-xob))]+aob*G1[which.min(abs(t-xob))]
  wei.t0=(xob<=t.0)*deltaob/G+(xob>t.0)/G.t0;   wei.t=(xob<=t)*deltaob/G+(xob>t)/G.t; 
  out=cbind(wei.t0,wei.t,G.t,G.t0,G1,G0,G)#,pte2
}


R.opt.event<- function(xone,xzero, deltaone, deltazero, sone, szero, t, landmark, std = FALSE, conf.int = FALSE, gopt = FALSE, ind=FALSE){
  t.0=landmark
  re=500 #resample number
  nn=200 #number grid points for numeric integration
  stept=0.1
  tt = seq(t.0,t,stept) #set of time points
  allx = c(xone, xzero)
  alldelta = c(deltaone,deltazero)
  alls = c(sone, szero)
  alla = c(rep(1, length(xone)), rep(0, length(xzero)))
  data = cbind(allx, alldelta, alla, alls)
  
  
  n=nrow(data)
  indexindex=sample(n, n/2, replace = FALSE)
  data1=data[indexindex,]
  data2=data[-indexindex,]
  n1=nrow(data1)
  n2=nrow(data2)
  
  causal=rep(NA,length(tt)); causal2=rep(NA,length(tt))
  causals=rep(NA,length(tt)); causals2=rep(NA,length(tt))
  causalind=rep(NA,length(tt)); causalind2=rep(NA,length(tt))
  causalsind=rep(NA,length(tt)); causalsind2=rep(NA,length(tt))
  for (j in 1:length(tt)){
    t=tt[j]
    ################ pte2 given data1  
    xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1
    
    ####
    from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
    s=seq(from, to, by = step)
       
    #### optimal function 
    {
      temp=WEIGHT(xob,deltaob,aob,n=n, t.0=t.0, t=t)
      wei.t0=temp[,1];wei.t=temp[,2]
      
      bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
      kern = Kern.FUN(zz=s,zi=sob,bw)
      
      f0.s.t0.t0.hat=apply((aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==0)*wei.t0)
      f0.s.t.t0.hat=apply((aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==0)*wei.t)
      f1.s.t0.t0.hat=apply((aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==1)*wei.t0)
      f1.s.t.t0.hat=apply((aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==1)*wei.t)
      temp=(sob>t.0); temp[is.na(temp)]=1
      p0.t0.t0.hat=sum((aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==0)*wei.t0)
      p0.t.t0.hat=sum((aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==0)*wei.t)
      p1.t0.t0.hat=sum((aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==1)*wei.t0)
      p1.t.t0.hat=sum((aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==1)*wei.t)
      
      integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
      
      mu0t=sum((aob==0)*(xob>t)*wei.t)/sum((aob==0)*wei.t)
      mu1t=sum((aob==1)*(xob>t)*wei.t)/sum((aob==1)*wei.t)
      mu0t0=sum((aob==0)*(xob>t.0)*wei.t0)/sum((aob==0)*wei.t0)
      mu1t0=sum((aob==1)*(xob>t.0)*wei.t0)/sum((aob==1)*wei.t0)
      
      integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
      
      lambda=nu/de
      lambda2=(mu0t-mu0t0*mu1t/mu1t0)/(mu0t0^2/mu1t0)
      
      gs.hat=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
      g3.hat=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
      g.hat=(lambda2*mu0t0+mu1t)/mu1t0
    }
    
    #### pte
    xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2
    temp=WEIGHT(xob,deltaob,aob,n=n, t.0=t.0, t=t)
    wei.t0=temp[,1];wei.t=temp[,2]
    
    # causal=mu1t-mu0t
    causal.t0=sum((xob>t.0)*aob*wei.t0)/sum(aob*wei.t0)-sum((xob>t.0)*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
    causal[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
    temp=(sob>t.0); temp[is.na(temp)]=1
    causals[j]=sum(aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum(aob*wei.t0)+
      sum((xob>t.0)*temp*g3.hat*aob*wei.t0)/sum(aob*wei.t0) -
      (sum((1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat[tempind],na.rm = T)/sum((1-aob)*wei.t0)+
         sum((xob>t.0)*temp*g3.hat*(1-aob)*wei.t0)/sum((1-aob)*wei.t0))
    pte.es=causals[j]/causal[j]
    
    causalind[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    temp=(xob>t.0); temp[is.na(temp)]=1
    causalsind[j]=sum((xob>t.0)*temp*g.hat*aob*wei.t0)/sum(aob*wei.t0) -
      sum((xob>t.0)*temp*g.hat*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
    pteind.es=causalsind[j]/causalind[j]
    
    ################ pte1 given data2  
    xob=data2[,1];deltaob=data2[,2];aob=data2[,3];sob=data2[,4];n=n2
    
    ####
    from = min(sob[sob!=0],na.rm = T); to = quantile(sob[sob!=0],.95,na.rm = T); step=((to - from)/nn)
    s=seq(from, to, by = step)
    
    #### optimal function 
    {
      temp=WEIGHT(xob,deltaob,aob,n=n, t.0=t.0, t=t)
      wei.t0=temp[,1];wei.t=temp[,2]
      
      bw = 1.06*sd(sob,na.rm=T)*n^(-1/5)/(n^0.06)
      kern = Kern.FUN(zz=s,zi=sob,bw)
      
      f0.s.t0.t0.hat=apply((aob==0)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==0)*wei.t0)
      f0.s.t.t0.hat=apply((aob==0)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==0)*wei.t)
      f1.s.t0.t0.hat=apply((aob==1)*(xob>t.0)*(sob<=t.0)*kern*wei.t0,2,sum,na.rm=T)/sum((aob==1)*wei.t0)
      f1.s.t.t0.hat=apply((aob==1)*(xob>t)*(sob<=t.0)*kern*wei.t,2,sum,na.rm=T)/sum((aob==1)*wei.t)
      temp=(sob>t.0); temp[is.na(temp)]=1
      p0.t0.t0.hat=sum((aob==0)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==0)*wei.t0)
      p0.t.t0.hat=sum((aob==0)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==0)*wei.t)
      p1.t0.t0.hat=sum((aob==1)*(xob>t.0)*temp*wei.t0, na.rm = T)/sum((aob==1)*wei.t0)
      p1.t.t0.hat=sum((aob==1)*(xob>t)*temp*wei.t, na.rm = T)/sum((aob==1)*wei.t)
      
      integrand<-f0.s.t0.t0.hat^2/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      de=temp+p0.t0.t0.hat^2/p1.t0.t0.hat
      
      mu0t=sum((aob==0)*(xob>t)*wei.t)/sum((aob==0)*wei.t)
      mu1t=sum((aob==1)*(xob>t)*wei.t)/sum((aob==1)*wei.t)
      mu0t0=sum((aob==0)*(xob>t.0)*wei.t0)/sum((aob==0)*wei.t0)
      mu1t0=sum((aob==1)*(xob>t.0)*wei.t0)/sum((aob==1)*wei.t0)
      
      integrand<-f0.s.t0.t0.hat*f1.s.t.t0.hat/f1.s.t0.t0.hat
      temp=(integrand[1] + integrand[nn+1] + 2*sum(integrand[seq(2,nn,by=2)]) + 4 *sum(integrand[seq(3,nn-1, by=2)]) )*step/3
      nu=mu0t-temp-p0.t0.t0.hat*p1.t.t0.hat/p1.t0.t0.hat
      
      lambda=nu/de
      lambda2=(mu0t-mu0t0*mu1t/mu1t0)/(mu0t0^2/mu1t0)
      
      gs.hat2=(lambda*f0.s.t0.t0.hat+f1.s.t.t0.hat)/f1.s.t0.t0.hat
      g3.hat2=(lambda*p0.t0.t0.hat+p1.t.t0.hat)/p1.t0.t0.hat
      g.hat2=(lambda2*mu0t0+mu1t)/mu1t0
      g3.es=(g3.hat+g3.hat2)/2
      gs.es=(gs.hat+gs.hat2)/2
      g.es=(g.hat+g.hat2)/2
    }
    
    #### pte
    xob=data1[,1];deltaob=data1[,2];aob=data1[,3];sob=data1[,4];n=n1
    temp=WEIGHT(xob,deltaob,aob,n=n, t.0=t.0, t=t)
    wei.t0=temp[,1];wei.t=temp[,2]
    
    # causal=mu1t-mu0t
    causal.t0=(causal.t0+sum((xob>t.0)*aob*wei.t0)/sum(aob*wei.t0)-sum((xob>t.0)*(1-aob)*wei.t0)/sum((1-aob)*wei.t0))/2
    causal2[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    tempind=c(sapply(1:n, function(kk){which.min(abs(sob[kk]-s))})); tempind=as.numeric(tempind)
    temp=(sob>t.0); temp[is.na(temp)]=1
    causals2[j]=sum(aob*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum(aob*wei.t0)+
      sum((xob>t.0)*temp*g3.hat2*aob*wei.t0)/sum(aob*wei.t0) -
      (sum((1-aob)*wei.t0*(xob>t.0)*(sob<=t.0)*gs.hat2[tempind],na.rm = T)/sum((1-aob)*wei.t0)+
         sum((xob>t.0)*temp*g3.hat2*(1-aob)*wei.t0)/sum((1-aob)*wei.t0))
    pte.es=(pte.es+causals2[j]/causal2[j])/2
    
    causalind2[j]=sum((xob>t)*aob*wei.t)/sum(aob*wei.t)-sum((xob>t)*(1-aob)*wei.t)/sum((1-aob)*wei.t)
    temp=(xob>t.0); temp[is.na(temp)]=1
    causalsind2[j]=sum((xob>t.0)*temp*g.hat2*aob*wei.t0)/sum(aob*wei.t0) -
      sum((xob>t.0)*temp*g.hat2*(1-aob)*wei.t0)/sum((1-aob)*wei.t0)
    pteind.es=(pteind.es+causalsind2[j]/causalind2[j])/2
  }
  
  mm=length(causal)-1
  integrand=causal
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causals
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermst=(rmsts+causal.t0)/(rmst+causal.t0)
  
  integrand=causal2
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causals2
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermst=(ptermst+(rmsts+causal.t0)/(rmst+causal.t0))/2
  
  integrand=causalind
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causalsind
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermstind=(rmsts+causal.t0)/(rmst+causal.t0)
  
  integrand=causalind2
  rmst=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  integrand=causalsind2
  rmsts=(integrand[1] + integrand[mm+1] + 2*sum(integrand[seq(2,mm,by=2)]) + 4 *sum(integrand[seq(3,mm-1, by=2)]) )*stept/3
  ptermstind=(ptermstind+(rmsts+causal.t0)/(rmst+causal.t0))/2
  
  causal.est=(causal[length(causal)]+causal2[length(causal2)])/2
  causals.est=(causals[length(causals)]+causals2[length(causals2)])/2
    
  if(std | conf.int) {
  	
  ######## variance of proposed estimators
  vv=matrix(rexp((n1+n2)*re),nrow=n1+n2)
  temp=apply(vv,2,resam,t,t.0,tt=5,data,data1,data2,indexindex)
  aa=which(temp[1,]>1 | temp[1,]<0)
  if (length(aa)>= re-1){
    pte.se=NA
      g3.se=NA
    gs.se=rep(NA,nn+1)
    conf.int.R = c(NA, NA)
  } else {
    pte.se=sd(temp[1,temp[1,]<=1 & temp[1,]>=0])
      g3.se=sd(temp[2,temp[1,]<=1 & temp[1,]>=0])
    gs.se=apply( temp[-(1:2),temp[1,]<=1 & temp[1,]>=0],1,sd)
    R.p = temp[1,temp[1,]<=1 & temp[1,]>=0]
    conf.l.quantile.R = quantile(R.p, 0.025)
	conf.u.quantile.R = quantile(R.p, 0.975)
    conf.int.R = as.vector(c(conf.l.quantile.R, conf.u.quantile.R))
  }
  }
  
 out=list("R.opt"=pte.es)
 if(ind) {out = c(out, list("R.opt.ind" =pteind.es ))}
 if(gopt) {out = c(out, list("g1.opt"=gs.es, "g2.opt"=as.numeric(g3.es) ))}
 if(std) {out = c(out, list("R.opt.sd"=pte.se ))}
  if(conf.int) {out = c(out, list("conf.int.R"=conf.int.R))}
 if(std & gopt) {out = c(out, list("g1.opt.sd"=as.vector(gs.se), "g2.opt.sd" = g3.se ))}
 return(out)
 
}
