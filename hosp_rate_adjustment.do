use ./data/hospitalizations_by_type_and_age_withpop_2016, clear

keep if disease=="all"

by f_age, so: egen pop_total=total(pop)
gen hosp_rate=hosp/pop
levelsof id, local(levels)
gen adj_hosp_rate_kreis=.
foreach l of local levels {
su hosp_rate [aw=pop_total] if id=="`l'"
replace adj_hosp_rate_kreis=1e6*r(mean) if id=="`l'"
}
ren adj adj_hosp_rate
keep id adj_hosp_rate
duplicates drop
destring id, replace
so id
save ./intermediate/hosp_stub, replace

use ./data/hospitalizations_by_type_and_age_withpop_2016, clear

keep if disease=="infectious"

by f_age, so: egen pop_total=total(pop)
gen hosp_rate=hosp/pop
levelsof id, local(levels)
gen adj_hosp_rate_kreis=.
foreach l of local levels {
su hosp_rate [aw=pop_total] if id=="`l'"
replace adj_hosp_rate_kreis=1e6*r(mean) if id=="`l'"
}
ren adj adj_infectious_hosp_rate
keep id adj_infectious_hosp_rate
duplicates drop
destring id, replace
so id
save ./intermediate/hosp_infectious_stub, replace

use ./data/hospitalizations_by_type_and_age_withpop_2016, clear

keep if disease=="respiratory"

by f_age, so: egen pop_total=total(pop)
gen hosp_rate=hosp/pop
levelsof id, local(levels)
gen adj_hosp_rate_kreis=.
foreach l of local levels {
su hosp_rate [aw=pop_total] if id=="`l'"
replace adj_hosp_rate_kreis=1e6*r(mean) if id=="`l'"
}
ren adj adj_respiratory_hosp_rate
keep id adj_respiratory_hosp_rate
duplicates drop
destring id, replace
so id
save ./intermediate/hosp_respiratory_stub, replace
