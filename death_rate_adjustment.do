use ./data/all_cause_mortality_by_age_withpop_2016, clear
by f_age, so: egen pop_total=total(pop)
gen death_rate=deaths/pop
levelsof id, local(levels)
gen adj_death_rate_kreis=.
foreach l of local levels {
su death_rate [aw=pop_total] if id=="`l'"
replace adj_death_rate_kreis=1e6*r(mean) if id=="`l'"
}
ren adj adj_total_rate
keep id adj_total_rate
duplicates drop
destring id, replace
so id
save ./intermediate/total_stub, replace

use ./data/infectious_disease_mortality_by_age_withpop_2016, clear
by f_age, so: egen pop_total=total(pop)
gen death_rate=deaths/pop
levelsof id, local(levels)
gen adj_death_rate_kreis=.
foreach l of local levels {
su death_rate [aw=pop_total] if id=="`l'"
replace adj_death_rate_kreis=1e6*r(mean) if id=="`l'"
}
ren adj adj_infectious_rate
keep id adj_infectious_rate
duplicates drop
destring id, replace
so id
save ./intermediate/infectious_stub, replace

use ./data/respiratory_disease_mortality_by_age_withpop_2016, clear
by f_age, so: egen pop_total=total(pop)
gen death_rate=deaths/pop
levelsof id, local(levels)
gen adj_death_rate_kreis=.
foreach l of local levels {
su death_rate [aw=pop_total] if id=="`l'"
replace adj_death_rate_kreis=1e6*r(mean) if id=="`l'"
}
ren adj adj_respiratory_rate
keep id adj_respiratory_rate
duplicates drop
destring id, replace
so id
save ./intermediate/respiratory_stub, replace


/*
use ./data/fulldata_kreise, clear
so id
merge id using total_stub
ta _me
keep if _me==3
drop _me
so id
merge id using infectious_stub
ta _me
keep if _me==3
drop _me
save ./intermediate/fulldata_kreise_age, replace
*/
