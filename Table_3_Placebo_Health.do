** options first
set seed 10101
global opt "cerrd"
global covs_rd "covs(S*)"
global covs_rdd "i.seg_fix_split i.seg_fix_splitXold"
global vce_rd "vce(cluster statenum)" //"vce(hc1)"
global vce_rdd "vce(cluster statenum)"
global boot_reps 99999

** table options next
cap global table "Table_2_Placebo_Health"
cap global title "Balance tests: Health Status 2016"
cap global tlabel "tab:placebo"

** table labels
global Lab1 "Panel A. Log all-cause mortality (2016)"
global Lab2 "Panel B. Log mortality from infectious diseases (2016)"
global Lab3 "Panel C. Log mortality from respiratory diseases (2016)"
global Lab4 "Panel D. Log hospitalizations (2016)"
global Lab5 "Panel E. Log hospitalizations for infectious diseases (2016)"
global Lab6 "Panel F. Log hospitalizations for respiratory diseases (2016) "

use ./data/fulldata_kreise, clear

keep id county type state lat lon dist2border ddr seg_fix_*
egen state_num=group(state)

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
save ./intermediate/distance_stub, replace
global pan=1
use ./data/all_cause_mortality_by_age_withpop_2016, clear
drop if id=="11000"
destring id, replace force
so id
merge id using ./intermediate/distance_stub
ta _me
drop _me
forvalues i=1/4 {
global col=`i'
if `i'==1 | `i'==2 {
global a_`i'=40
global d_`i'=`i'
}
if `i'==3 | `i'==4 {
global a_`i'=25
global d_`i'=`i'-2
}
preserve
drop if f_age==${a_`i'}
keep if (f_age<=${a_`i'}+${d_`i'}*5 & f_age>=${a_`i'}-${d_`i'}*5)
gen old=(f_age<=${a_`i'}+${d_`i'}*5 & f_age>${a_`i'})
collapse (sum) deaths pop, by(id ddr lat lon d*_1 seg_fix_split old d state)
gen mort=ln(1+1e6*deaths/pop)
foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
} 

xi i.seg_fix_split, pref(S) noomit
cap drop statenum
egen statenum=group(state)

macro drop b0l b0r b1l b1r

rdrobust mort d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust mort d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b1l = e(h_l)
local b1r = e(h_r)

keep if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.1f `=-`b0l''
local b0r: di %3.1f `=`b0r''
local b1l: di %3.1f `=-`b1l''
local b1r: di %3.1f `=`b1r''

noi reg mort ddr ddrXold d*_1 d*_1Xold ${covs_rdd}, ${vce_rdd}

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

restore
}


global pan=1+${pan}
use ./data/infectious_disease_mortality_by_age_withpop_2016, clear
drop if id=="11000"
destring id, replace force
so id
merge id using ./intermediate/distance_stub
ta _me
drop _me
forvalues i=1/4 {
global col=`i'
if `i'==1 | `i'==2 {
global a_`i'=40
global d_`i'=`i'
}
if `i'==3 | `i'==4 {
global a_`i'=25
global d_`i'=`i'-2
}
preserve
drop if f_age==${a_`i'}
keep if (f_age<=${a_`i'}+${d_`i'}*5 & f_age>=${a_`i'}-${d_`i'}*5)
gen old=(f_age<=${a_`i'}+${d_`i'}*5 & f_age>${a_`i'})
collapse (sum) deaths pop, by(id ddr lat lon d*_1 seg_fix_split old d state)
gen mort=ln(1+1e6*deaths/pop)

foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
} 

xi i.seg_fix_split, pref(S) noomit
cap drop statenum
egen statenum=group(state)

macro drop b0l b0r b1l b1r

rdrobust mort d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust mort d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd}
local b1l = e(h_l)
local b1r = e(h_r)

keep if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.1f `=-`b0l''
local b0r: di %3.1f `=`b0r''
local b1l: di %3.1f `=-`b1l''
local b1r: di %3.1f `=`b1r''

noi reg mort ddr ddrXold d*_1 d*_1Xold ${covs_rdd}, ${vce_rdd}

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

restore
}


global pan=1+${pan}
use ./data/respiratory_disease_mortality_by_age_withpop_2016, clear
drop if id=="11000"
destring id, replace force
so id
merge id using ./intermediate/distance_stub
ta _me
drop _me
forvalues i=1/4 {
global col=`i'
if `i'==1 | `i'==2 {
global a_`i'=40
global d_`i'=`i'
}
if `i'==3 | `i'==4 {
global a_`i'=25
global d_`i'=`i'-2
}
preserve
drop if f_age==${a_`i'}
keep if (f_age<=${a_`i'}+${d_`i'}*5 & f_age>=${a_`i'}-${d_`i'}*5)
gen old=(f_age<=${a_`i'}+${d_`i'}*5 & f_age>${a_`i'})
collapse (sum) deaths pop, by(id ddr lat lon d*_1 seg_fix_split old d state)
gen mort=ln(1+1e6*deaths/pop)
foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
} 

xi i.seg_fix_split, pref(S) noomit
cap drop statenum
egen statenum=group(state)

macro drop b0l b0r b1l b1r

rdrobust mort d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust mort d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd}
local b1l = e(h_l)
local b1r = e(h_r)

keep if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.1f `=-`b0l''
local b0r: di %3.1f `=`b0r''
local b1l: di %3.1f `=-`b1l''
local b1r: di %3.1f `=`b1r''

noi reg mort ddr ddrXold d*_1 d*_1Xold ${covs_rdd}, ${vce_rdd}

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

restore
}

use ./data/hospitalizations_by_type_and_age_withpop_2016, clear
drop if disease==""
drop ICD
reshape wide hosp, i(id f_age l_age pop) j(disease) string
drop if id=="11000"
destring id, replace force
so id
merge id using ./intermediate/distance_stub
ta _me
drop _me
foreach m in all infectious respiratory {
global pan=1+${pan}
forvalues i=1/4 {
global col=`i'
if `i'==1 | `i'==2 {
global a_`i'=40
global d_`i'=`i'
}
if `i'==3 | `i'==4 {
global a_`i'=25
global d_`i'=`i'-2
}
preserve
drop if f_age==${a_`i'}
keep if (f_age<=${a_`i'}+${d_`i'}*5 & f_age>=${a_`i'}-${d_`i'}*5)
gen old=(f_age<=${a_`i'}+${d_`i'}*5 & f_age>${a_`i'})
collapse (sum) hosp`m' pop, by(id ddr lat lon d*_1 seg_fix_split old d state)
gen depvar=ln(1+1e6*hosp`m'/pop)

foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
} 

xi i.seg_fix_split, pref(S) noomit
cap drop statenum
egen statenum=group(state)

macro drop b0l b0r b1l b1r

rdrobust depvar d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust depvar d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd}
local b1l = e(h_l)
local b1r = e(h_r)

keep if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.1f `=-`b0l''
local b0r: di %3.1f `=`b0r''
local b1l: di %3.1f `=-`b1l''
local b1r: di %3.1f `=`b1r''

noi reg depvar ddr ddrXold d*_1 d*_1Xold ${covs_rdd}, ${vce_rdd}

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
	
restore
}
}
	
	
cap global col1=3


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
