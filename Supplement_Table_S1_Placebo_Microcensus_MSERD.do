** options first
set seed 10101
global opt "mserd"
global covs_rd "covs(F*)"
global covs_rdd "i.seg_fix_split i.seg_fix_splitXold"
global vce_rd "vce(cluster statenum)" //"vce(hc1)"
global vce_rdd "vce(cluster statenum)"
global boot_reps 99999

** table options next
cap global table "Supplement_Table_S1_Microcensus_MSERD"
cap global title "Balance tests: Microcensus 2017"
cap global tlabel "tab:placebo_fdz_mserd"

** table labels
global Lab1 "Panel A. Log labour force participation"
global Lab2 "Panel B. Log unemployed"
global Lab3 "Panel C. Log full time employment"
global Lab4 "Panel D. Log public employment"
global Lab5 "Panel E. Log commuters"


use ./data/export_kreise_nonoverlapping, clear

ren state state_fdz

gen id = string(county)
drop if id=="11000"

destring id, replace force
tab id

drop county

save ./intermediate/fdz_microcensus, replace


use ./data/fulldata_kreise, clear

keep id county type state lat lon dist2border ddr seg_fix_*
egen statenum=group(state)
xi i.seg_fix_split, pref(F) noomit


clonevar d=dist2border
replace d=-d if ddr==0

*Polynomials in distance to DDR border
forvalues i=1/1 {
forvalues j=1/`i' {
gen d_`j'_`i'=d^`j'
gen dXddr_`j'_`i'=ddr*d^`j'
}
}

so id
merge 1:m id using ./intermediate/fdz_microcensus


drop _me
sort id cohort

local allvars in_lf_ilo unemployed_status fulltime public_employment commuter 

local f_cohort "2 1 5 4"
local l_cohort "3 4 6 7"



forv i = 1(1)4 {
global col=`i'
global pan=1
preserve

local f_c : word `i' of `f_cohort'
local l_c : word `i' of `l_cohort'

di "first cohort: `f_c'"
di "last cohort: `l_c'"
di "diff cohort: `=`l_c'-`f_c''"

gen old=(cohort==`f_c') if cohort>=`f_c' & cohort<=`l_c'
if (`=`l_c'-`f_c''>1) {
replace old = 1 if cohort==`=`f_c'+1'
}

tab old cohort

collapse (mean) `allvars' (first) statenum F* ddr d d*_1 seg_fix_split [aw=W], by(old id)


foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
} 


local shares in_lf_ilo unemployed_status fulltime public_employment commuter  
foreach v of local shares {
replace `v' = log(1+1e6*`v')

}



foreach cvar of local allvars {


rdrobust `cvar' d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust `cvar' d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b1l = e(h_l)
local b1r = e(h_r)

cap drop samp
gen samp = 0
replace samp = 1 if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | ///
		    (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.1f `=-`b0l''
local b0r: di %3.1f `=`b0r''
local b1l: di %3.1f `=-`b1l''
local b1r: di %3.1f `=`b1r''

noi reg `cvar' ddr ddrXold d*_1 d*_1Xold ${covs_rdd} if samp==1,  ${vce_rdd}

local b_ddrXold : di %3.2f `=round(_b[ddrXold],.01)'
local se_ddrXold : di %3.2f `=round(_se[ddrXold],.01)'

local pv=2*ttail(e(df_r),abs(_b[ddrXold]/_se[ddrXold]))
local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

global b_ddrXold_${col}_${pan}=`b_ddrXold'
global se_ddrXold_${col}_${pan}=`se_ddrXold'
global bpse_ddrXold_${col}_${pan}="`b_ddrXold'`strs' (`se_ddrXold')"
global b0_ddrXold_${col}_${pan} ="[`b0l', `b0r']"
global b1_ddrXold_${col}_${pan} ="[`b1l', `b1r']"
global bwboth_ddrXold_${col}_${pan} ="`b0r', `b1r'"
	
cap boottest ddrXold,  boot(wild) nograph weight(webb) reps(${boot_reps})
mat CI = r(CI)
local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'


global ci_ddrXold_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"
global pan=1+${pan}

}

restore

}

	
cap global col1=3
global pan=${pan}-1

/* End of Computation */


*Here we construct the table itself


#delimit ;
cap global stringc=""; forvalues i=1/$col {; cap global stringc="$stringc"+"c"; };
#delimit cr

#delimit;
cap global stringnum=""; forvalues i=1/$col {; cap global stringnum="$stringnum"+"& (`i') "; };
#delimit cr

#delimit;
cap global stringsp=""; forvalues i=1/$col {; cap global stringsp="$stringsp"+"& "; };
#delimit cr
	
	
quietly {
cap log close
log using ./output/l${table}.txt, replace

noi di "\documentclass[9pt]{report}"
noi di "\pagestyle{plain}"
noi di "\usepackage{amsmath}"
noi di "\usepackage{geometry}"
*noi di "\usepackage{tabularx}"
noi di "\geometry{right=.5in,left=0.5in,top=.5in,bottom=.5in}"
noi di "\begin{document}"
noi di "\pagenumbering{gobble}"

noi di "\begin{table}[h!]"
noi di "\centering"
noi di "\caption{${title}}\label{${tlabel}}"

noi di "\begin{tabular*}{\textwidth}{@{\extracolsep\fill} l $stringc}"
noi di "\hline\hline \noalign{\smallskip}"
noi di "& \multicolumn{4}{c}{\it Policy change} \\"
noi di "& \multicolumn{2}{c}{1974 West ends universal} & \multicolumn{2}{c}{1990 East ends mandatory}  \\ \noalign{\smallskip} \cline{2-5}\noalign{\smallskip}"
noi di "& \multicolumn{4}{c}{\it Age interval around policy change } \\" 
noi di "&  5 years & 10 years & 5 years &  10 years  \\" 

noi di "$stringnum \\ \noalign{\smallskip} \hline\noalign{\smallskip} "


forvalues p=1(1)$pan {
global p=`p'
global V "ddrXold"
foreach v in $V {
cap global v `v'

foreach S in b bpse ci bwboth b0 b1 {
cap global S `S'
cap global string${S}_${v}_${p}=""
forvalues c=1(1)$col {
cap global c=`c'
cap global string${S}_${v}_${p}="${string${S}_${v}_${p}}"+"& ${${S}_${v}_${c}_${p}} "
}
}
}


noi di "\multicolumn{`=$col+1'}{l}{\it ${Lab${p}}} \\ \noalign{\smallskip}"

foreach v in $V {
global v `v'

noi di "$\text{\sc East}\times \text{\sc Treated}$   ${stringbpse_${v}_${p}} \\ \noalign{\smallskip}"
noi di "WCR 95\% CI	${stringci_${v}_${p}} \\"
noi di "BWs $(h_{NT}, h_{T})$            ${stringbwboth_${v}_${p}} \\ \noalign{\smallskip}"
nois di ""


}

}

	foreach S in  N {
	cap global S `S'
	cap global string${S}=""
	forvalues c=1/$col {
	cap global c=`c'
	cap global string${S}="${string${S}}"+"& ${${S}_${c}} "
	}
	}

noi di "\hline\hline"
noi di "\end{tabular*}"
noi di "\end{table}"

noi di "\end{document}"
log close
}

** strip tex out of log
qui shell perl -s tablemod.pl -tablenew="./output/${table}" -table="./output/l${table}"

** mp's tex
*shellout using ./intermediate/"${table}.tex"
qui !texify -p -c -b --run-viewer ${table}.tex

** rb's tex 
cd ./output/
!pdflatex -interaction=nonstopmode ${table}.tex
!rm ${table}.aux
!rm ${table}.log
!rm l${table}.txt

/*End of Table Construction*/
