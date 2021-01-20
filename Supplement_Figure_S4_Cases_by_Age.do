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


keep id state seg_fix_state lat lon ddr cases_all_*


reshape long cases_all_age_ , i(id state seg_fix_state lat lon ddr) j(age) 
ren cases_all_age_ cases


collapse (sum) cases, by(age ddr)

reshape wide cases , i(age) j(ddr) 

tw bar cases0 age, ///
	graphregion(color(white)) bgcolor(white) ///
	ytitle("Cases per age group") xtitle("Age")
graph export ./output/Supplement_Figure_S4_a.pdf, replace
tw bar cases1 age, color(maroon)  ///
	graphregion(color(white)) bgcolor(white) ///
	ytitle("Cases per age group") xtitle("Age")
graph export ./output/Supplement_Figure_S4_b.pdf, replace


use ./intermediate/basedata, clear


keep id state seg_fix_state lat lon ddr cases_c_*


reshape long cases_c_age_ , i(id state seg_fix_state lat lon ddr) j(age) 
ren cases_c_age_ cases


collapse (sum) cases, by(age ddr)

reshape wide cases , i(age) j(ddr) 

tw bar cases0 age, ///
	graphregion(color(white)) bgcolor(white) ///
	ytitle("Severe cases per age group") xtitle("Age")
graph export ./output/Supplement_Figure_S4_c.pdf, replace
tw bar cases1 age, color(maroon)  ///
	graphregion(color(white)) bgcolor(white) ///
	ytitle("Severe cases per age group") xtitle("Age")
graph export ./output/Supplement_Figure_S4_d.pdf, replace

