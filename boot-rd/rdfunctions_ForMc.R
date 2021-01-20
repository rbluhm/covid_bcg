## Copyright (c) 2016, Ottavio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md at https://github.com/grayclhn/boot-rd

## Additions by Richard Bluhm, December, 2020
## i) cluster bootstrap, 
## ii) partialling out covariates, 
## iii) different wild weights
## iv) parallel computation, and
## v) different BWs on each side using new rdbwselect
## Removed support for residual bootstrap
## Also available under same MIT "Expat" License

kweight <- rdrobust:::rdrobust_kweight

# Mammen's 2-point distribution
mammen_values <- c(1 - sqrt(5), 1 + sqrt(5)) / 2
mammen_weights <- c(sqrt(5) + 1, sqrt(5) - 1) / (2 * sqrt(5))

# Webb's 6-point distribution
webb_values <- c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2))
webb_weights <- rep(1/6, 6)

# Rademacher's 2-point distribution
rademacher_values <- c(-1,1)
rademacher_weights <- rep(1/2, 2)

boot_ci_basic <- function(estimate, boots, boot_parameter, a = 0.05) {
 ci <- as.vector(estimate) - quantile(boots, c(1 - a/2, a/2)) + as.vector(boot_parameter)
 names(ci) <- rev(names(ci))
 ci
}

## this function generates bias-corrected estimator
boot_estimator_wild <- function(ypl, ypr, yql, yqr, 
                                 coef.ql, coef.qr, coef.pl, coef.pr,
                                 Xql, Xqr, Xpl, Xpr,
                                 WXql, WXqr, e.ql.adj, e.qr.adj,
                                 wqr, wql, ihr, ihl, Nbc, gen.wild) {
  
  b.ql       <- crossprod(WXql, yql)
  fitted.l   <- Xql %*% b.ql
  residual.l <- (yql - fitted.l) * e.ql.adj
  
  b.qr       <- crossprod(WXqr, yqr)
  fitted.r   <- Xqr %*% b.qr
  residual.r <- (yqr - fitted.r) * e.qr.adj
  
  boot_parameter <- coef.qr %*% yqr - coef.ql %*% yql
  estimate       <- coef.pr %*% ypr - coef.pl %*% ypl
  
  foo <- function(gen.wild) {
     wild.e <- gen.wild()
     return(coef.pr %*% (fitted.r[ihr] + residual.r[ihr]* wild.e$wild.q.e.R[ihr]) -
                         coef.pl %*% (fitted.l[ihl] + residual.l[ihl]* wild.e$wild.q.e.L[ihl]))
  }
    
  boots <- replicate(Nbc,foo(gen.wild))

  return(as.numeric(estimate - mean(boots) + boot_parameter))
}

## this function generates bootstrap distribution of bias-corrected estimator
boot_dist_wild <- function(ypl, ypr, yql, yqr, 
                            coef.ql, coef.qr, coef.pl, coef.pr,
                            Xql, Xqr, Xpl, Xpr,
                            WXql, WXqr, e.ql.adj, e.qr.adj,
                            wqr, wql, ihr, ihl, Nbc, Nci, gen.wild,
                            parallel = FALSE, ncores = NULL){
  
  b.ql       <- crossprod(WXql, yql)
  fitted.l   <- Xql %*% b.ql
  residual.l <- (yql - fitted.l) * e.ql.adj
  
  b.qr       <- crossprod(WXqr, yqr)
  fitted.r   <- Xqr %*% b.qr
  residual.r <- (yqr - fitted.r) * e.qr.adj
  
  # use parallel for speed gains   
   if (parallel==TRUE) {
   
      if (is.null(ncores)) ncores=parallel::detectCores()
                           
       boots <- mclapply(1:Nci,function(t) { 
       
          wild.e <- gen.wild()
          newyqr <- fitted.r + residual.r* wild.e$wild.q.e.R
          newyql <- fitted.l + residual.l* wild.e$wild.q.e.L
        
          boot_estimator_wild(newyql[ihl], newyqr[ihr], newyql, newyqr,
                               coef.ql, coef.qr, coef.pl, coef.pr,
                               Xql, Xqr, Xpl, Xpr,
                               WXql, WXqr, e.ql.adj, e.qr.adj,
                               wqr, wql, ihr, ihl, Nbc, gen.wild)
                               
       }, mc.cores=ncores)
       
       boots <- unlist(boots)    
   
   } else {
   
        boots <- replicate(Nci, {

        wild.e <- gen.wild()
        newyqr <- fitted.r + residual.r* wild.e$wild.q.e.R
        newyql <- fitted.l + residual.l* wild.e$wild.q.e.L
    
        boot_estimator_wild(newyql[ihl], newyqr[ihr], newyql, newyqr,
                             coef.ql, coef.qr, coef.pl, coef.pr,
                             Xql, Xqr, Xpl, Xpr,
                             WXql, WXqr, e.ql.adj, e.qr.adj,
                             wqr, wql, ihr, ihl, Nbc, gen.wild)})
   
   }                   
                           
  return(boots)
}

## a wrapper for both point and interval estimator
rdboot_wild <- function(y, x, a = 0.05, Nbc = 500, Nci = 999, p = 1, q = 2, 
                         kernel = c("triangular", "uniform", "epanechnikov"),
                         residual = c("hc0","hc1", "hc2", "hc3"), bwselect = "mserd", h = NULL, b = NULL,
                         cluster = NULL, covs = NULL, 
                         parallel = FALSE, ncores = NULL, wilddist = NULL){
  
  # to improve speed, (1) create all necessary objects only once and pass
  # them into functions for double bootstrap. (2) avoid multiplication of
  # matrix with large size.
  
  kernel <- match.arg(kernel)
  residual <- match.arg(residual)

  if ( wilddist == "rademacher" ) {
        wild_values <- rademacher_values
        wild_weights <- rademacher_weights
    } else if ( wilddist == "webb" ) {
        wild_values <- webb_values
        wild_weights <- webb_weights
    } else {
        wild_values <- mammen_values
        wild_weights <- mammen_weights
    }
  
  if (!is.null(cluster)) {
      vce = "hc0"
      } else {
      vce = "nn"
    }
      
  if (parallel == TRUE) {
        require(parallel)
        RNGkind("L'Ecuyer-CMRG")
    }

  if (!is.null(h) & is.null(b)) b=h

  # allow fixed BWs too
  if (!is.null(h)) {
    hl <- hr <- h; bl <- br <- b
  } else {
  bw <- rdbwselect(y, x, p=p, q=q, kernel=kernel, bwselect=bwselect, cluster=cluster, vce=vce, covs=covs)$bws
  hl <- bw[1]; hr <- bw[2]; bl <- bw[3]; br <- bw[4]  
  }
  
  ## partial out covariates 
  if (!is.null(covs)) {
  
      d <- x > 0
      temp.df <- cbind(y=y, x=x, d=as.numeric(d), covs)
      f4m <- as.formula(paste("y~ d + x + x:d +", paste( names(covs), collapse="+")))
      
      ## add weights
      wl <- kweight(abs(x), 0, abs(hl), kernel = kernel)
      wr <- kweight(abs(x), 0, abs(hr), kernel = kernel)
      wl[ x >= 0 ] <- wr[ x >= 0 ]
      w <- wl; rm(wr,wl)
      
      m_covs <- lm(f4m, data = temp.df, subset = (x > -hl & x < hr), weights=w)

      pr.df <- data.frame(x = x, d = d)
      pr.df[] <- 0 
      pr.df <- cbind(pr.df, covs)
      y <- y - predict(m_covs, pr.df)    
      rm(d, temp.df, f4m, w, m_covs, pr.df)

  }
  
  yql <- y[x > -max(hl,bl) & x < 0]
  xql <- x[x > -max(hl, bl) & x < 0]
  yqr <- y[x >= 0 & x < max(hr, br)]
  xqr <- x[x >= 0 & x < max(hr, br)]
  
  ihl <- xql > -hl
  ihr <- xqr < hr
  
  ypl <- y[x > -hl & x < 0]
  xpl <- x[x > -hl & x < 0]
  ypr <- y[x >= 0 & x < hr]
  xpr <- x[x >= 0 & x < hr]
  
  #### clustering
  
    if (!is.null(cluster)) {

      i.hb <- (x > -max(hl, bl) & x < max(hr, br))
      i.hb.L <- x[i.hb] < 0; i.hb.R <- !i.hb.L  
      cluster <- cluster[i.hb]

      gen.wild <- function() {
        e <- vector(length = length(cluster))
        for (i in unique(cluster)) {
          e[cluster == i] <- sample(wild_values, 1, prob = wild_weights)
        }
        
      return(list(wild.q.e.L = e[i.hb.L], wild.q.e.R = e[i.hb.R]))
      }
    } else {
      
      gen.wild <- function() {
        return(list(wild.q.e.L = sample(wild_values, length(yql), T, wild_weights),
                    wild.q.e.R = sample(wild_values, length(yqr), T, wild_weights)))
      }
    }
      
  ####
  
  # a vector of weight for residual bootstrap
  wql <- kweight(xql, 0, bl, kernel)
  wqr <- kweight(xqr, 0, br, kernel)
  
  # orthogonal polynomials
  xql.poly <- poly(xql, q)
  xqr.poly <- poly(xqr, q)
  xpl.poly <- poly(xpl, p)
  xpr.poly <- poly(xpr, p)
  
  # design matrix
  Xql <- cbind(1, poly(xql, q))
  Xqr <- cbind(1, poly(xqr, q))
  Xpl <- cbind(1, poly(xpl, p))
  Xpr <- cbind(1, poly(xpr, p))
  
  KXql <- kweight(xql, 0, bl, kernel) * Xql
  KXqr <- kweight(xqr, 0, br, kernel) * Xqr
  KXpl <- kweight(xpl, 0, hl, kernel) * Xpl
  KXpr <- kweight(xpr, 0, hr, kernel) * Xpr
  
  # parameter maker
  WXql <- t(solve(crossprod(Xql, KXql), t(KXql)))
  WXqr <- t(solve(crossprod(Xqr, KXqr), t(KXqr)))
  WXpl <- t(solve(crossprod(Xpl, KXpl), t(KXpl)))
  WXpr <- t(solve(crossprod(Xpr, KXpr), t(KXpr)))
  
  # intercept maker
  coef.ql <- tcrossprod(c(1, predict(xql.poly, 0)), WXql)
  coef.qr <- tcrossprod(c(1, predict(xqr.poly, 0)), WXqr)
  coef.pl <- tcrossprod(c(1, predict(xpl.poly, 0)), WXpl)
  coef.pr <- tcrossprod(c(1, predict(xpr.poly, 0)), WXpr)
  
  # HC adjustment: diagnal vector of the projection matrix
  KXql.sqrt <- sqrt(kweight(xql, 0, bl, kernel)) * Xql
  KXqr.sqrt <- sqrt(kweight(xqr, 0, br, kernel)) * Xqr
  
  h_l <- diag(KXql.sqrt %*% solve(crossprod(Xql, KXql), t(KXql.sqrt)))
  h_r <- diag(KXqr.sqrt %*% solve(crossprod(Xqr, KXqr), t(KXqr.sqrt)))
  
  e.ql.adj <- switch (residual,
                      hc0 = 1,
                      hc1 = sqrt(sum(kweight(xql, 0, bl, kernel))/
                                   (sum(kweight(xql, 0, bl, kernel)) - q - 1)),
                      hc2 = 1/sqrt(1 - h_l),
                      hc3 = 1/(1 - h_l)
  )
  
  e.qr.adj <- switch (residual,
                      hc0 = 1,
                      hc1 = sqrt(sum(kweight(xqr, 0, br, kernel))/
                                   (sum(kweight(xqr, 0, br, kernel)) - q - 1)),
                      hc2 = 1/sqrt(1 - h_r),
                      hc3 = 1/(1 - h_r)
  )
  
  
  
  estimate <- boot_estimator_wild(ypl, ypr, yql, yqr, 
                                   coef.ql, coef.qr, coef.pl, coef.pr,
                                   Xql, Xqr, Xpl, Xpr,
                                   WXql, WXqr, e.ql.adj, e.qr.adj,
                                   wqr, wql, ihr, ihl, Nbc, gen.wild)
  
  boots <- boot_dist_wild(ypl, ypr, yql, yqr, 
                           coef.ql, coef.qr, coef.pl, coef.pr,
                           Xql, Xqr, Xpl, Xpr,
                           WXql, WXqr, e.ql.adj, e.qr.adj,
                           wqr, wql, ihr, ihl, Nbc, Nci, gen.wild, parallel, ncores)
  
  boot_parameter <- coef.qr %*% yqr - coef.ql %*% yql
  
  ci_basic <- boot_ci_basic(estimate, boots, boot_parameter, a)
  ci_percentile <- quantile(boots, c(a/2, 1-a/2))
  se_basic <- sqrt(var(boots))
  
  result <- rbind(c(estimate, se_basic, ci_basic, hl,hr,bl,br), c(estimate, se_basic, ci_percentile, hl,hr,bl,br))
  row.names(result) <- c("Basic CI", "Percentile CI")
  dimnames(result)[[2]] <-  c("Coef.", "Std. Err.",  "CI LB", "CI UB", "h(l)", "h(r)", "b(l)", "b(r)") 
  return(result)
}
