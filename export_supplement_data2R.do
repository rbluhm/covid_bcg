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

*Polynomials in distance to DDR border
forvalues i=1/1 {
forvalues j=1/`i' {
gen d_`j'_`i'=d^`j'
gen dXddr_`j'_`i'=ddr*d^`j'
}
}

xi i.seg_fix_split, pref(F) noomit

egen statenum=group(state)

save ./intermediate/fulldata_for_R_otherjumps, replace

