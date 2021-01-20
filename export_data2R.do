cap global D=50

cap confirm file ./intermediate/cases_simulation_donut_${D}.dta
if _rc!=0 {
do simulation_donut.do
}

cap confirm file ./intermediate/facebook_stub_noself.dta
if _rc!=0 {
do facebook.do
}


use ./data/fulldata_kreise, clear

cap drop mob*
so id
merge 1:1 id using ./intermediate/cases_simulation_donut_${D}
ta _me
keep if _me==3
drop _me
gen lcp_sim=ln(1+1e6*cases_simulation/pop_latest)


preserve
use ./data/cases_pop_cat_c_singleages_week17, clear
egen severe_cases_total=rowtotal(age_*_cases_cat_c)
destring id, replace
keep id severe_cases_total
so id
save ./intermediate/severe_cases_stub, replace
restore
so id
merge 1:1 id using ./intermediate/severe_cases_stub
ta _me
keep if _me==3
drop _me
gen lcp_c_total=ln(1+1e6*severe_cases_total/pop_latest)

preserve
foreach s in all c {
foreach wk in 17 50 {
use ./data/cases_pop_cat_`s'_singleages_week`wk', clear
egen cases_`s'_week`wk'_total=rowtotal(age_*_cases_cat_`s')
keep id cases_`s'_week`wk'_total
destring id, replace
so id
save ./intermediate/stub_`s'_`wk', replace
}
}
use ./intermediate/stub_all_17, clear
keep id
foreach s in all c {
foreach wk in 17 50 {
so id
merge 1:1 id using ./intermediate/stub_`s'_`wk'
ta _me
drop if _me==2
drop _me
save ./intermediate/cases_time_stub, replace
}
}
foreach s in all c {
gen new_cases_`s'=cases_`s'_week50-cases_`s'_week17
}
so id
save ./intermediate/cases_time_stub, replace
restore

so id
merge 1:1 id using ./intermediate/cases_time_stub
ta _me
drop if _me==2
drop _me
foreach s in all c {
gen lncp_`s'_total=ln(1+1e6*new_cases_`s'/pop_latest)
}

so nuts_id
merge 1:1 nuts_id using ./intermediate/facebook_stub_noself
ta _me
keep if _me==3
drop _me
so id

preserve
use ./data/kreise_counts_panel_daily, clear
keep if date>date("26apr2020","DMY") & date<=date("13dec2020","DMY") 
collapse (sum) new_deaths, by(id)
destring id, replace
so id
save ./intermediate/new_deaths_stub, replace
restore

so id
merge 1:1 id using ./intermediate/new_deaths_stub
ta _me
keep if _me==3
drop _me
so id

so id
merge 1:1 id using ./intermediate/commute_donut_${D}_stub
ta _me
keep if _me==3
drop _me
gen ratio_incoming=incoming0/(incoming0+incoming1)
gen ratio_outgoing=outgoing0/(outgoing0+outgoing1)


gen cp= 1e6*cum_cases/pop_latest
gen lcp_total=ln(1+cp)
gen dp= 1e6*cum_deaths/pop_latest
gen ldp_total=ln(1+dp)
gen ndp=1e6*new_deaths/pop_latest
gen lndp_total=ln(1+ndp)

clonevar d=dist2border
replace d=-d if ddr==0

gen lpopdens=ln(pop_latest/area)
gen larea = ln(area)
gen ldisp=ln(disp_income_2017)

gen percent_age_over60 = 100* (pop_age_60_79 + pop_age_80p) / pop_latest
gen percent_age_under35 = 100* (pop_age_00_04 + pop_age_05_14 + pop_age_15_35) / pop_latest

*Polynomials in distance to DDR border
forvalues i=1/1 {
forvalues j=1/`i' {
gen d_`j'_`i'=d^`j'
gen dXddr_`j'_`i'=ddr*d^`j'
}
}

xi i.seg_fix_split, pref(F) noomit

egen statenum=group(state)

save ./intermediate/fulldata_for_R, replace
 
