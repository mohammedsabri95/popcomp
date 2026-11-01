# Independent analytic and contract checks; no test framework dependency.
# Fabricated posterior arrays below are deterministic test fixtures, not fitted estimates.
library(popcomp)
results <- list()
check <- function(name, code) {
  err <- tryCatch({ force(code); NULL }, error = function(e) conditionMessage(e))
  results[[length(results)+1L]] <<- data.frame(test=name, pass=is.null(err),
                                            detail=if (is.null(err)) "OK" else err)
}
near <- function(x,y,tol=1e-12) stopifnot(isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tol)))
rejects <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error=identity),"error"))

rejects_message <- function(expr, text) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(err, "error"), grepl(text, conditionMessage(err), fixed = TRUE))
}

a <- c(1,2,4,8)
q <- rbind(c(.1,.2,.3,.4), c(.4,.3,.2,.1))
Y <- rbind(c(0,0,0,1),c(1,0,0,0),c(0,1,1,0))
unit <- c(1L,2L,1L); time <- c(1L,1L,2L); X <- cbind(x=c(-1,0,1))
make <- function(y=Y, ref=q, tt=time, xx=X) pc_data(y,unit,tt,ref,xx,a)
obj <- make()
check("softmax normalisation and finite extreme logits", {
  z <- pc_softmax(rbind(c(1000,1001,999,0),c(-1000,-1001,-999,-1002)))
  stopifnot(all(is.finite(z)),all(z>=0)); near(rowSums(z),c(1,1))
})
check("softmax additive-shift invariance", near(pc_softmax(rbind(a)),pc_softmax(rbind(a+10000))))
check("softmax one-category composition", near(pc_softmax(matrix(c(2,3),2,1)),c(1,1)))
check("W1 and mean for point masses", {
  near(pc_W1(rbind(c(1,0,0,0)),c(0,0,0,1),a),7)
  near(pc_M(rbind(c(1,0,0,0)),c(0,0,0,1),a),-7)
})
check("W1 detects shape change with zero mean", {
  p <- rbind(c(.5,0,.5)); ref <- c(0,1,0)
  near(pc_M(p,ref,c(0,1,2)),0); near(pc_W1(p,ref,c(0,1,2)),1)
})
check("row-specific references and identity", {
  near(pc_W1(q,q,a),c(0,0)); near(pc_M(q,q,a),c(0,0))
  near(pc_W1(q,q[2:1,],a),rep(sum(abs(cumsum(q[1,]-q[2,]))[-4]*diff(a)),2))
})
check("orthonormal and centred basis", {
  b <- pc_basis(a); B <- cbind(b$B,b$Bh)
  near(crossprod(B),diag(3)); near(colSums(B),rep(0,3))
})
check("invalid ordered support rejected by basis", rejects(pc_basis(c(1,2,2,4))))
check("too few categories rejected", rejects(pc_basis(c(1,2,3))))
check("zero-total rows explicitly rejected", {z<-Y;z[1,]<-0;rejects_message(make(y=z), "rows with n = 0 carry no information")})
check("one-case rows and zero cells accepted with correct offsets", {
  stopifnot(identical(obj$n,c(1,1,2)),all(is.finite(obj$data$off)))
  near(obj$data$y,as.vector(t(Y)))
  near(obj$data$off,as.vector(t(log(q[unit,])))+rep(log(rowSums(Y)),each=4))
})
check("unnormalised reference rejected", rejects(make(ref=q*2)))
check("unnamed covariates rejected", rejects(make(xx=unname(X))))
check("BYM2 requires graph", rejects(pc_formula(obj)))
check("missing predictor rejected", rejects(pc_regionyear(list(predictor=NULL),obj)))
check("posterior summaries and strict signed tolerance", {
  v<-c(-2,-1,0,1,2); s<-pc_tolerance(cbind(E=v),delta=2)
  near(s$mean,0);near(s$sd,sqrt(2.5));near(s$median,0)
  near(c(s$q025,s$q975),c(-1.9,1.9));near(s$p_gt0,.4)
  near(s$ub975_abs,2);near(s$p_within_2,.6)
})
check("distance tolerance uses strict inequality", {
  s<-pc_tolerance(cbind(W1=c(0,1,2,3)),delta=2,distance=TRUE)
  near(s$p_within_2,.5)
})
S<-3L
dr<-list(S=S,beta1=matrix(c(.2,.4,.6),1),beta2=matrix(0,1,S),
         gamma1=matrix(0,2,S),gamma2=matrix(0,2,S),
         u1=matrix(0,2,S),u2=matrix(0,2,S),z=NULL,sigma_eps=NULL)
pts<-list(hi=c(x=1),lo=c(x=0))
cn<-list(dM=list(type="M",weights=c(hi=1,lo=-1)),W1=list(type="W1",pair=c("hi","lo")))
ref<-colMeans(q)
expected<-t(vapply(c(.2,.4,.6),function(b) {
  w<-ref*exp(obj$B[,1]*b);p<-w/sum(w)
  c(sum((p-ref)*a),sum(abs(cumsum(p-ref))[-4]*diff(a)))
},numeric(2)))
check("typical standardised contrasts against direct formula", {
  got<-pc_contrasts(dr,obj,pts,cn,reference=ref,version="typical")
  stopifnot(identical(dim(got),c(S,2L)));near(got,expected)
})
check("marginal agrees when random and period effects vanish", {
  near(pc_contrasts(dr,obj,pts,cn,reference=ref,version="marginal",R=20,seed=8),expected)
})
check("default reference weights observed rows, not case totals", {
  r<-colSums(q*tabulate(unit,2))/3
  near(pc_contrasts(dr,obj,pts,cn,version="typical"),
       pc_contrasts(dr,obj,pts,cn,reference=r,version="typical"))
})
check("fixed-seed marginal Monte Carlo reproducibility", {
  dr2<-dr;dr2$u1[,]<-c(-.2,.3);dr2$sigma_eps<-rep(.2,S)
  r1<-pc_contrasts(dr2,obj,pts,cn,reference=ref,R=20,seed=44)
  r2<-pc_contrasts(dr2,obj,pts,cn,reference=ref,R=20,seed=44)
  stopifnot(identical(r1,r2))
})
check("region-year summaries against direct compositions", {
  pp<-rbind(c(.4,.3,.2,.1),c(.1,.2,.3,.4),rep(.25,4))
  dp<-list(S=2,predictor=matrix(rep(as.vector(t(log(pp))),2),ncol=2))
  rr<-pc_regionyear(dp,obj)
  near(rr$W1[,1],pc_W1(pp,q[unit,],a));near(rr$M[,1],pc_M(pp,q[unit,],a))
  stopifnot(all(rr$W1>=abs(rr$M)-1e-12))
})
# These public API contract tests must remain failures until approved fixes exist.
check("one contrast returns S by 1 matrix", {
  got<-pc_contrasts(dr,obj,pts,cn["dM"],reference=ref,version="typical")
  stopifnot(identical(dim(got),c(S,1L)));near(got[,1],expected[,1])
})
check("negative counts rejected", {z<-Y;z[1,]<-c(-1,0,0,2);rejects_message(make(y=z), "non-negative integer counts")})
check("fractional counts rejected", {z<-Y;z[1,]<-c(.5,0,0,.5);rejects_message(make(y=z), "non-negative integer counts")})
check("compatible structural-zero reference remains supported", {
  z <- q; z[1, ] <- c(0, .2, .3, .5)
  got <- make(ref=z)
  is_zero <- as.vector(t(z[unit, , drop=FALSE] == 0))
  stopifnot(all(got$data$off[is_zero] == -Inf),
            all(is.finite(got$data$off[!is_zero])))
  prob <- pc_softmax(matrix(got$data$off, ncol=length(a), byrow=TRUE))
  near(rowSums(prob), rep(1, nrow(Y)))
  stopifnot(all(prob[unit == 1, 1] == 0))
})
check("invalid time index rejected", rejects_message(make(tt=c(0L,1L,2L)), "positive integer indices"))
check("nonfinite exposure rejected", {z<-X;z[1,1]<-Inf;rejects_message(make(xx=z), "finite numeric covariates")})
check("unordered support rejected by W1", rejects_message(pc_W1(rbind(c(1,0,0,0)),c(0,0,0,1),rev(a)), "strictly increasing"))
check("unknown contrast type rejected", {
  bad<-cn;bad$W1$type<-"TYPO";rejects_message(pc_contrasts(dr,obj,pts,bad,reference=ref,version="typical"), "exactly 'M' or 'W1'")
})
check("positive counts at a zero reference probability rejected as incompatible", {
  ref0 <- q; ref0[1, ] <- c(0, .2, .3, .5)
  y0 <- Y; y0[1, ] <- c(1, 0, 0, 0)
  rejects_message(make(y=y0, ref=ref0), "incompatible with the model")
})
check("fractional time indices rejected", {
  rejects_message(make(tt=c(.5, 1, 2)), "positive integer indices")
})
check("valid zero case cells remain accepted", {
  got <- make()
  near(got$data$y, as.vector(t(Y)))
  stopifnot(all(is.finite(got$data$off)))
})
check("W1 symmetry and zero-mass boundary", {
  p <- rbind(c(.5, 0, .5)); ref <- c(0, 1, 0); ages <- c(0, 1, 2)
  near(pc_W1(p, ref, ages), 1)
  near(pc_W1(p, ref, ages), pc_W1(rbind(ref), as.vector(p), ages))
  near(pc_W1(p, p, ages), 0)
})
check("W1 rejects malformed support and probabilities", {
  rejects_message(pc_W1(rbind(c(.5, .5)), c(.5, .5), c(1, 1)), "strictly increasing")
  rejects_message(pc_W1(rbind(c(.5, .5)), c(1, 1), c(1, 2)), "summing to one")
  rejects_message(pc_W1(rbind(c(.5, .5)), 1, c(1, 2)), "one probability")
})
r<-do.call(rbind,results)
print(r,row.names=FALSE)
cat(sprintf("Release tests: %d passed; %d failed; %d total\n",sum(r$pass),sum(!r$pass),nrow(r)))
# An optional external output location keeps package checks free of source writes.
out<-Sys.getenv("POPCOMP_TEST_RESULTS","")
if (nzchar(out)) write.csv(r,out,row.names=FALSE)
if (any(!r$pass)) stop("Release contract failures: see printed results; scientific implementation unchanged")
