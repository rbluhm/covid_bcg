# required
lop <- c("rdrobust", "tidyverse", "haven", "parallel", "foreach", "miscFuncs")
newp <- lop[!(lop %in% installed.packages()[,"Package"])]
if(length(newp)) install.packages(newp)
lapply(lop, require, character.only = TRUE)

# code for rd boot
source("./boot-rd/rdfunctions_ForMc.R")

# grab data
bcg.df <- read_dta("./intermediate/fulldata_for_R_otherjumps.dta") %>% select(-starts_with("mob_"))
cov.df <- bcg.df %>%  select(starts_with("Fseg_fix_")) %>% select(-Fseg_fix_sp_5)
cov.df

# simulation options
RNGkind("L'Ecuyer-CMRG")
set.seed(10101)
reps.tau <- 1e5-1
reps.ci <- 1e3

# helper function
stars.pval <- function(p.value) {
  unclass(symnum(p.value, corr = FALSE, na = FALSE, 
                 cutpoints = c(0,0.01, 0.05, 0.1, 1), symbols = c("***", "**",  "*", "")))
}

# function to fit CCT, run bootstrap and generate table cells
build.tbl <- function(depvar=NULL, deplabel = NULL, distvar=NULL, covs=NULL, cluster=NULL, h=NULL, bwselect=NULL) {

  if(is.na(h)) h <- NULL
  if(is.na(bwselect)) bwselect <- "mserd" # have to provide some value, even if fixed h is used
  
  collab <- h
  if(is.null(h)) collab <- bwselect
  
  cct <- rdrobust(y=depvar, x=distvar, covs = covs, kernel = "uniform", all=T,  
                  cluster =cluster, vce="hc0", bwselect = bwselect, h=h)
  summary(cct)
  
  rdboot.ci <- rdboot_wild(depvar, distvar,  covs = covs, a = 0.05, Nbc = reps.tau, Nci = reps.ci,  
                           residual = "hc0", kernel = "uniform", bwselect = bwselect, h = h,
                           cluster = cluster, parallel= TRUE, wilddist = "webb")
  
  print(rdboot.ci)
  
  coef <- round(cct$coef[3], 2)
  se <- round(cct$se[3], 2)
  strs <- stars.pval(cct$pv[3])
  N_l <- cct$N_h[1]
  N_r <- cct$N_h[2]
  h <- round(rdboot.ci[1,5], 1)
  b <- round(rdboot.ci[1,7], 1)
  
  est <- round(rdboot.ci[1,1], 2)
  ci_lb <- round(rdboot.ci[1,3], 2)
  ci_ub <- round(rdboot.ci[1,4], 2)
  
  e1 <- paste0(coef, strs, " (", se, ")")
  e2 <- paste0("[",ci_lb, ", ", ci_ub, "]")
  e3 <- paste0(h, ", ", b)
  e4 <- paste0(N_l, " $|$ ", N_r)
  ret <- rbind("",e1,est,e2,e3,e4)
  colnames(ret) <- collab
  rownames(ret) <- c(deplabel, "\\sc{East}", "Wild Est.", "Wild 95\\% CI", "BWs $(h, b)$", "N (West $|$ East)")
  return(ret)  
}

# rows
dvs <- c("ldisp",  "lpopdens", "larea", "percent_age_over60", "percent_age_under35", "elapsed_2020_04_26")

stopifnot(dvs %in% names(bcg.df))

dv_labs <- c("Panel A. Disposable income per capita", 
             "Panel B. Population density", 
             "Panel C. Area", 
             "Panel D. Percent older than 60", 
             "Panel E. Percent younger than 35",
             "Panel F. Days since first case")

add_ctrls <- rep(NA, length(dvs))

# columns
hs <- c(NA, NA, 30, 60)
bwsels <- c("cerrd", "mserd", NA, NA)

## start timer for all operations
start.time <- Sys.time()

# build table
full.tbl <- foreach(r = 1:length(dvs), .combine = rbind) %:%
              foreach(c = 1:length(hs), .combine = cbind) %do% {
                
                if(!is.na(add_ctrls[r])) {
                  cov.c.df <- cbind(cov.df, bcg.df[add_ctrls[r]])
                } else {
                  cov.c.df <- cov.df
                }
              
                build.tbl(depvar=bcg.df[[dvs[r]]], deplabel = dv_labs[r], 
                          distvar=bcg.df[["d"]], covs=cov.c.df, cluster=bcg.df[["statenum"]], 
                          h=hs[c], bwselect=bwsels[c])
              
              }

full.tbl

## time code
time.taken <- Sys.time() - start.time
print("the complete code took")
print(time.taken)

# create tex output
sink("./output/Supplement_Table_S4_OtherRDs.tex")
latextable(full.tbl)
sink()

#read and modify

fin.tex <- readLines("./output/Supplement_Table_S4_OtherRDs.tex")
fin.tex <- fin.tex[which(fin.tex == "\\begin{table}[htbp]"):which(fin.tex == "\\end{table}")]

for(i in 1:length(dv_labs)) {
  lb.ind <- grepl(dv_labs[i], fin.tex, fixed = TRUE)
  fin.tex[lb.ind] <- paste0("\\multicolumn{5}{l}{\\it ", dv_labs[i], "} \\\\ \\noalign{\\smallskip}")  
}

fin.tex

lb.ind <- grepl("Wild Est.", fin.tex, fixed = TRUE)
fin.tex[lb.ind] <- "\\noalign{\\smallskip}"
lb.ind <- grepl("N (West $|$ East)", fin.tex, fixed = TRUE)
fin.tex[lb.ind] <- "\\noalign{\\smallskip}"

fileConn<-file("./output/Supplement_Table_S4_OtherRDs.tex")
writeLines(fin.tex, fileConn)
close(fileConn)
