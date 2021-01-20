use ./data/kreise_distances_panel, clear
destring id_kreis_a id_kreis_b, replace
save ./data/kreise_distances_panel_int, replace

use ./data/commuters, clear
so id_kreis_a id_kreis_b
merge 1:1 id_kreis_a id_kreis_b using ./data/kreise_distances_panel_int
ta _me
drop _me
replace incoming=0 if dist<${D}
replace outgoing=0 if dist<${D}

preserve
use ./data/fulldata_kreise, clear
keep id ddr
duplicates drop
so id
save ./intermediate/ddr_stub, replace
restore

preserve
clonevar id=id_kreis_b
so id
merge id using ./intermediate/ddr_stub
ta _me
drop _me
replace ddr=1 if ddr==. /*This is Berlin*/
drop id
clonevar id=id_kreis_a
collapse (sum) incoming outgoing, by(id ddr)
reshape wide incoming outgoing, i(id) j(ddr)
so id
save ./intermediate/commute_donut_${D}_stub, replace
restore

*Housekeeping
global T=60
global R0=2.5
global gamma=1/7
global beta=${R0}*${gamma}
destring id_*, replace
ren incoming m_a
ren outgoing m_b

preserve
use ./data/fulldata_kreise, clear
keep id pop_latest
so id
save ./intermediate/pop_stub, replace
restore


foreach s in a b {
cap drop id
clonevar id=id_kreis_`s'
so id
merge id using ./intermediate/pop_stub
ta _me
drop _me
ren pop_latest pop_`s'
replace pop_`s'=3613495 if id==11000

so id
merge id using ./data/initial_cases
ta _me
drop _me
ren cases_de cases_init_`s'
replace cases_init_`s'=0 if cases_init_`s'==.
drop id
}

by id_kreis_a, so: egen m_a_sum=total(m_a)
by id_kreis_b, so: egen m_b_sum=total(m_b)

*Initialization
gen S_0_a=pop_a
gen I_0_a=cases_init_a
gen R_0_a=0

gen S_0_b=pop_b
gen I_0_b=cases_init_b
gen R_0_b=0

*Iterative step
local t=0
while `t'<$T {
noi di "`t'"
if `t'>22 {
replace m_a=0
replace m_b=0
}
cap drop x_a x_b mx_a mx_b mx_a_sum mx_b_sum Itilde_a Itilde_b
gen x_a=I_`t'_a/pop_a
gen x_b=I_`t'_b/pop_b
gen mx_a=m_a*x_a
gen mx_b=m_b*x_b
by id_kreis_a, so: egen mx_a_sum=total(mx_b)
by id_kreis_b, so: egen mx_b_sum=total(mx_a)
gen Itilde_a=I_`t'_a+pop_a*mx_a_sum/(pop_a+m_a_sum)
gen Itilde_b=I_`t'_b+pop_b*mx_b_sum/(pop_b+m_b_sum)
gen S_`=`t'+1'_a=S_`t'_a-${beta}*S_`t'_a*Itilde_a/pop_a
gen I_`=`t'+1'_a=I_`t'_a+${beta}*S_`t'_a*Itilde_a/pop_a-${gamma}*I_`t'_a
gen R_`=`t'+1'_a=R_`t'_a+${gamma}*I_`t'_a
gen S_`=`t'+1'_b=S_`t'_b-${beta}*S_`t'_b*Itilde_b/pop_b
gen I_`=`t'+1'_b=I_`t'_b+${beta}*S_`t'_b*Itilde_b/pop_b-${gamma}*I_`t'_b
gen R_`=`t'+1'_b=R_`t'_b+${gamma}*I_`t'_b
local t=`t'+1
drop x_a x_b mx_a mx_b mx_a_sum mx_b_sum Itilde_a Itilde_b 
}

gen cases_simulation=pop_a-S_${T}_a
ren id_*_a id
keep id cases_simulation
duplicates drop
so id
save ./intermediate/cases_simulation_donut_${D}, replace
