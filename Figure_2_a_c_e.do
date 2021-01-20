use ./data/fulldata_kreise, clear
so id

keep id-date dist2border ddr seg_fix_split
drop date

** merge in from daily panel, so that data is really week 17 ending in apr 26
merge 1:m nuts_id  using ./data/kreise_counts_panel_daily, force update

drop if _me == 2
keep if date==td(26apr2020)
drop _me

merge 1:1 nuts_id using ./data/cases_pop_cat_all_singleages_week17, force
ta _me
drop if _me==2
drop _me

*drop cases_age_*
ren age_*_cases_cat_all cases_age_*

local k=10
local j0=45-`k'
local j1=45+`k'

foreach v in cases_age pop_age {
egen `v'_15_28=rowtotal(`v'_15-`v'_28)
egen `v'_30_44=rowtotal(`v'_30-`v'_44)
egen `v'_`j0'_44=rowtotal(`v'_`j0'-`v'_44)
egen `v'_46_`j1'=rowtotal(`v'_46-`v'_`j1')
egen `v'_46_59=rowtotal(`v'_46-`v'_59)
egen `v'_total=rowtotal(`v'_0-`v'_74)
}

foreach s in total `j0'_44 46_`j1' 15_28 30_44 46_59 {
gen lcp_`s'=ln(1+1e6*cases_age_`s'/pop_age_`s')
}

local list "lcp_total lcp_`j0'_44 lcp_46_`j1'"


global Tlcp_total "Log(1+cases/million)"
global Tlcp_`j0'_44 "Log(1+cases/million) `j0'-44 Year Olds"
global Tlcp_46_`j1' "Log(1+cases/million)  46-`j1' Year Olds"

clonevar d=dist2border
replace d=-d if ddr==0

*RD Graphs
foreach v of varlist `list' {
preserve
rdplot `v' d, binselect(qsmv) genvars hide
keep rdplot_mean_y rdplot_mean_bin
save ./intermediate/binstub, replace
restore
append using ./intermediate/binstub

foreach b in 30 {
twoway lpolyci `v' d if (d>0 ), kernel(rectangle) bwidth(`b') degree(1) ///
    || lpolyci `v' d if (d<0 ), kernel(rectangle) bwidth(`b') degree(1) ///
    || scatter rdplot_mean_y rdplot_mean_bin, sort msize(med) xline(0) mcolor(black)  ///
	xtitle("Distance to former border, kilometers") ///
	ytitle("${T`v'}") ylabel(, angle(horiz)) ///
	legend(off) graphregion(color(white)) bgcolor(white)
	
graph export ./output/Figure_2_all_`v'.pdf, replace
drop rdplot_*
}
}

