** options first
set seed 10101
global opt "cerrd"
global covs_rd "covs(S*)"
global covs_rdd "i.seg_fix_split i.seg_fix_splitXold"
global vce_rd "vce(cluster statenum)" //"vce(hc1)"
global vce_rdd "vce(cluster statenum)"
global boot_reps 99999

** table options next
cap global table "Table_4_RD_DD_1974"
cap global title "RD-DD results: 1974 experiment"
cap global tlabel "tab:main_1974"

global Lab1 "Panel A. 1974 baseline (cases/million)"
global Lab2 "Panel B. 1974 differences-in-differences (cases/million)"
global Lab3 "Panel C. 1974 severity measure (symptomatic cases/million)"
global Lab4 "Panel D. 1974 differences-in-differences (symptomatic cases/million)"

use ./data/fulldata_kreise, clear
drop mob* pop_age_* 
so id
save ./intermediate/basedata, replace

foreach cat in all c {

use ./data/cases_pop_cat_`cat'_singleages_week17, clear
cap ren age_*_cases_cat_`cat' cases_`cat'_age_*
drop *_80p *_unk*
destring id, replace
so id
save ./intermediate/cases_`cat'_pop_singleages, replace

use ./intermediate/basedata, clear
so id
merge id using ./intermediate/cases_`cat'_pop_singleages
ta _me
drop if _me==2
drop _me
so id
save ./intermediate/basedata, replace
}

use ./intermediate/basedata, clear
foreach f in 21 29 45 {
foreach v in cases_all_age cases_c_age pop_age {
local fm10=`f'-10
local fm1=`f'-1
local fp1=`f'+1
local fp10=`f'+10
forvalues i=`fm10'/`fm1' {
cap egen `v'_`i'_`fm1'=rowtotal(`v'_`i'-`v'_`fm1')
}
forvalues i=`fp1'/`fp10' {
cap egen `v'_`fp1'_`i'=rowtotal(`v'_`fp1'-`v'_`i')
}
cap egen `v'_total=rowtotal(`v'_0-`v'_74)
}


foreach cat in all c {
forvalues i=`fm10'/`fm1' {
local s `i'_`fm1'
cap gen lcp_`cat'_`s'=ln(1+1e6*cases_`cat'_age_`s'/pop_age_`s')
}

forvalues i=`fp1'/`fp10' {
local s `fp1'_`i'
cap gen lcp_`cat'_`s'=ln(1+1e6*cases_`cat'_age_`s'/pop_age_`s')
}
cap gen lcp_`cat'_total=ln(1+1e6*cases_`cat'_age_total/pop_age_total)
}
}


gen cp= 1e6*cum_cases/pop_latest
gen dp= 1e6*cum_deaths/pop_latest

gen ldp=ln(1+dp)
gen lcp=ln(1+cp)




gen lpop_latest=ln(pop_latest)

clonevar d=dist2border
replace d=-d if ddr==0

gen lpopdens=ln(pop_latest/area)


*Polynomials in distance to DDR border
forvalues i=1/1 {
forvalues j=1/`i' {
gen d_`j'_`i'=d^`j'
gen dXddr_`j'_`i'=ddr*d^`j'
}
}


*Baseline 

global col=1
forv i = 1(2)9 {

local f=45
local fm11=`f'-11
local fm1=`f'-1
local fp1=`f'+1
local fp11=`f'+11
local fm10=`f'-10
local fp10=`f'+10
local a `=`fm11'+`i''_`fm1'
local b `fp1'_`=`fp11'-`i''

preserve
global pan=1
keep id state seg_fix_split lat lon ddr d*_1 d lcp_all_`b' lcp_all_`a'
reshape long lcp_all_ , i(id state seg_fix_split lat lon ddr d*_1 d) j(age_group) str
gen old=(age_group=="`b'")
foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
}
lab var ddrXold "\\sc{East} $\times$ \\{Treated}"

xi i.seg_fix_split, pref(S) noomit
cap drop statenum
egen statenum=group(state)

rdrobust lcp_all_ d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all 
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust lcp_all_ d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd}
local b1l = e(h_l)
local b1r = e(h_r)

gen samp = 0
replace samp = 1 if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.2f `=-`b0l''
local b0r: di %3.2f `=`b0r''
local b1l: di %3.2f `=-`b1l''
local b1r: di %3.2f `=`b1r''

reg lcp_all_ ddr ddrXold d*_1 d*_1Xold ${covs_rdd} if samp==1, ${vce_rdd}

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

* DD version 1974 baseline		

global pan=${pan}+1

reg lcp_all_ ddr ddrXold old, ${vce_rdd}

sum d
local b0l = r(min)
local b0r = r(max)
local b1l = r(min)
local b1r = r(max)

local b0l: di %3.2f `=-`b0l''
local b0r: di %3.2f `=`b0r''
local b1l: di %3.2f `=-`b1l''
local b1r: di %3.2f `=`b1r''

local b_ddrXold : di %3.2f `=round(_b[ddrXold],.01)'
local se_ddrXold : di %3.2f `=round(_se[ddrXold],.01)'

local pv=2*ttail(e(df_r),abs(_b[ddrXold]/_se[ddrXold]))
local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

global b_ddrXold_${col}_${pan}=`b_ddrXold'
global se_ddrXold_${col}_${pan}=`se_ddrXold'
global bpse_ddrXold_${col}_${pan}="`b_ddrXold'`strs' (`se_ddrXold')"
global b0_ddrXold_${col}_${pan} ="[`b0l', `b0r']"
global b1_ddrXold_${col}_${pan} ="[`b1l', `b1r']"
global bwboth_ddrXold_${col}_${pan} ="`b0l', `b1l'"

boottest ddrXold,  boot(wild) nograph weight(webb) reps(${boot_reps})
mat CI = r(CI)
local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'

global ci_ddrXold_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

sum d

local b0l: di %3.2f `=-abs(r(min))'
local b0r: di %3.2f `=r(max)'

global b0_ddrXold_${col}_${pan} ="[`b0l', `b0r']"
global b1_ddrXold_${col}_${pan} ="[`b0l', `b0r']"

global pan=${pan}-1

	
global col=1+$col
restore
}
global pan=1+${pan}

*Category C

global pan=1+${pan}
global col=1
forv i = 1(2)9 {

local f=45
local fm11=`f'-11
local fm1=`f'-1
local fp1=`f'+1
local fp11=`f'+11
local fm10=`f'-10
local fp10=`f'+10
local a `=`fm11'+`i''_`fm1'
local b `fp1'_`=`fp11'-`i''

preserve
keep id state seg_fix_split lat lon ddr d*_1 d lcp_c_`b' lcp_c_`a'
reshape long lcp_c_ , i(id state seg_fix_split lat lon ddr d*_1 d) j(age_group) str
gen old=(age_group=="`b'")
foreach v of varlist ddr d*_1 seg_fix_split {
gen `v'Xold=`v'*old
}

xi i.seg_fix_split, pref(S) noomit
cap drop statenum
egen statenum=group(state)

rdrobust lcp_c_ d if old==0, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd} all 
local b0l = e(h_l)
local b0r = e(h_r)

rdrobust lcp_c_ d if old==1, kernel(uni) bwselect(${opt}) ${covs_rd} ${vce_rd}
local b1l = e(h_l)
local b1r = e(h_r)

gen samp = 0
replace samp = 1 if (old==0 & ddr==0 & abs(d)<`b0l') | (old==0 & ddr==1 & abs(d)<`b0r') | (old==1 & ddr==0 & abs(d)<`b1l') | (old==1 & ddr==1 & abs(d)<`b1r')

local b0l: di %3.2f `=-`b0l''
local b0r: di %3.2f `=`b0r''
local b1l: di %3.2f `=-`b1l''
local b1r: di %3.2f `=`b1r''

reg lcp_c_ ddr ddrXold d*_1 d*_1Xold ${covs_rdd} if samp==1, ${vce_rdd}

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

* DD version 1974 severity	

global pan=${pan}+1

reg lcp_c_ ddr ddrXold old, ${vce_rdd}

sum d
local b0l = r(min)
local b0r = r(max)
local b1l = r(min)
local b1r = r(max)

local b0l: di %3.2f `=-`b0l''
local b0r: di %3.2f `=`b0r''
local b1l: di %3.2f `=-`b1l''
local b1r: di %3.2f `=`b1r''

local b_ddrXold : di %3.2f `=round(_b[ddrXold],.01)'
local se_ddrXold : di %3.2f `=round(_se[ddrXold],.01)'

local pv=2*ttail(e(df_r),abs(_b[ddrXold]/_se[ddrXold]))
local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

global b_ddrXold_${col}_${pan}=`b_ddrXold'
global se_ddrXold_${col}_${pan}=`se_ddrXold'
global bpse_ddrXold_${col}_${pan}="`b_ddrXold'`strs' (`se_ddrXold')"
global b0_ddrXold_${col}_${pan} ="[`b0l', `b0r']"
global b1_ddrXold_${col}_${pan} ="[`b1l', `b1r']"
global bwboth_ddrXold_${col}_${pan} ="`b0l', `b1l'"	

boottest ddrXold,  boot(wild) nograph weight(webb) reps(${boot_reps})
mat CI = r(CI)
local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'

global ci_ddrXold_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

sum d

local b0l: di %3.2f `=-abs(r(min))'
local b0r: di %3.2f `=r(max)'

global b0_ddrXold_${col}_${pan} ="[`b0l', `b0r']"
global b1_ddrXold_${col}_${pan} ="[`b0l', `b0r']"

global pan=${pan}-1

	
global col=1+$col
restore
}
global pan=1+${pan}

	
cap global col1=${col}
cap global col=${col1}-1




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
noi di "\usepackage{tabularx}"
noi di "\geometry{right=.5in,left=0.5in,top=.5in,bottom=.5in}"
noi di "\begin{document}"
noi di "\pagenumbering{gobble}"

noi di "\begin{table}[h!]"
noi di "\centering"
noi di "\caption{${title}}\label{${tlabel}}"

noi di "\begin{tabular*}{\textwidth}{@{\extracolsep\fill} l $stringc}"
noi di "\hline\hline \noalign{\smallskip}"
noi di "& \multicolumn{5}{c}{\it Age interval around policy change } \\"
noi di "&  10 years & 8 years & 6 years &  4 years &  2 years \\"

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



*noi di "\hline"
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
