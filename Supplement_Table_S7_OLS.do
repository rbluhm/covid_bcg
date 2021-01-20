** boot options
set seed 10101
global boot_reps 99999

** table options next
cap global table "Supplement_Table_S7_OLS"
cap global title "Bivariate OLS regressions of log(1+cases/million) on control variables"
cap global tlabel "tab:app_ols"


cap confirm file ./intermediate/total_stub.dta
if _rc!=0 {
do death_rate_adjustment.do
do hosp_rate_adjustment.do
}

use ./data/fulldata_kreise, clear


* age adjusted death rates
merge 1:1  id using ./intermediate/total_stub
drop _merge

gen ln_adj_total_rate = ln(1+adj_total_rate)

merge 1:1  id using ./intermediate/infectious_stub
drop _merge

gen ln_adj_infectious_rate = ln(1+adj_infectious_rate)


merge 1:1  id using ./intermediate/respiratory_stub
drop _merge

gen ln_adj_respiratory_rate = ln(1+adj_respiratory_rate)

* age adjusted hosp rates
merge 1:1  id using ./intermediate/hosp_stub
drop _merge

gen ln_adj_hosp_rate = ln(1+adj_hosp_rate)

merge 1:1  id using ./intermediate/hosp_infectious_stub
drop _merge

gen ln_adj_infectious_hosp_rate = ln(1+adj_infectious_hosp_rate)

merge 1:1  id using ./intermediate/hosp_respiratory_stub
drop _merge

gen ln_adj_respiratory_hosp_rate = ln(1+adj_respiratory_hosp_rate)

drop if id == 11000

ren pop_2018  pop
gen cp= 1e6*cum_cases/pop
gen dp= 1e6*cum_deaths/pop

gen ldp=ln(1+dp)
gen lcp=ln(1+cp)

gen lpopdens=ln(pop_latest/area)
gen larea = ln(area)
gen ldisp=ln(disp_income_2017)

gen percent_age_over60 = 100* (pop_age_60_79 + pop_age_80p) / pop
gen percent_age_under35 = 100* (pop_age_00_04 + pop_age_05_14 + pop_age_15_35) / pop

clonevar d=dist2border
replace d=-d if ddr==0


xi i.seg_fix_state, pref(S) noomit
xi i.seg_fix_split, pref(F) noomit

egen statenum=group(state)

ren ln_adj_infectious_rate ln_adj_inf_rate
ren ln_adj_respiratory_rate ln_adj_res_rate
ren ln_adj_infectious_hosp_rate ln_adj_inf_hrate
ren ln_adj_respiratory_hosp_rate ln_adj_res_hrate


global col=1
foreach l in 0 1 All {
global pan=1


global Lab${pan} "Panel A. Disposable income per capita"

cap global varlist ldisp

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}	

global Lab${pan} "Panel B. Population density"

cap global varlist lpopdens

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}	


global Lab${pan} "Panel C. Percent older than 60" 

cap global varlist percent_age_over60

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N) 
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}	

global Lab${pan} "Panel D. Percent younger than 35"

cap global varlist percent_age_under35

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
	

global Lab${pan} "Panel E. Age-adjusted overall death rate per million"

cap global varlist ln_adj_total_rate

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N) 
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}


global Lab${pan} "Panel G. Age-adjusted infectious diseases death rate per million"

cap global varlist ln_adj_inf_rate

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}

global Lab${pan} "Panel H. Age-adjusted respiratory diseases death rate per million"

cap global varlist ln_adj_res_rate

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}


global Lab${pan} "Panel I. Age-adjusted hospitalization rate per million"

cap global varlist ln_adj_hosp_rate

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}


global Lab${pan} "Panel J. Age-adjusted infectious diseases hospitalization rate per million"

cap global varlist ln_adj_inf_hrate

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N) 
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N) 
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}


global Lab${pan} "Panel K. Age-adjusted respiratory diseases hospitalization rate per million"

cap global varlist ln_adj_res_hrate

if "`l'"!="All" {
	reg lcp ${varlist} if ddr==`l', cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}
if "`l'"=="All" {
	reg lcp ${varlist}, cluster(state)
	local b : di %3.2f `=round(_b[${varlist}],.01)'
	local se : di %3.2f `=round(_se[${varlist}],.01)'
	local pv=2*ttail(e(df_r),abs(_b[${varlist}]/_se[${varlist}]))
	local strs "`=cond(`pv'>0.1,"",cond(`pv'>0.05,"*",cond(`pv'>0.01,"**","***")))'"

	global b_${varlist}_${col}_${pan}=`b'
	global se_${varlist}_${col}_${pan}=`se'
	global bpse_${varlist}_${col}_${pan}="`b'`strs' (`se')"

	boottest ${varlist},  boot(wild) nograph weight(webb) reps(${boot_reps})
	mat CI = r(CI)
	local wild_ci_l : di %3.2f `=round(CI[1,1],.01)'
	local wild_ci_h : di %3.2f `=round(CI[1,2],.01)'
	global ci_${varlist}_${col}_${pan}="[`wild_ci_l', `wild_ci_h']"

	cap global N_${col}=e(N)
	cap global VV${pan} "${varlist}"
	global pan=1+${pan}
}

global col=1+${col}
}	
	
	
cap global col1=${col}
cap global col=${col1}-1
cap global pan=${pan}-1




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
noi di "& West & East & All  \\" 

noi di "$stringnum \\ \noalign{\smallskip} \hline\noalign{\smallskip} "


forvalues p=1/$pan {
global p=`p'
global V "${VV`p'}"
foreach v of varlist $V {
cap global v `v'

foreach S in b bpse ci {
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

noi di "Coefficient (Std. Err.)  ${stringbpse_${v}_${p}} \\ \noalign{\smallskip}"
noi di "WCR 95\% CI	${stringci_${v}_${p}} \\ \noalign{\smallskip}"


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

noi di "Observations ${stringN} \\"
noi di "\hline"
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


