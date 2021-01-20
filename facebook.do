use ./data/facebook, clear
drop if id_kreis_a==id_kreis_b
gen nuts3_a=substr(id_kreis_a,3,1)
gen nuts3_b=substr(id_kreis_b,3,1)
gen ddr_a=inlist(nuts3_a,"4","8","D","E","G")
gen ddr_b=inlist(nuts3_b,"4","8","D","E","G")
by id_kreis_a, so: egen total=total(scaled_sci)
by id_kreis_a, so: egen total_ddr=total(scaled_sci) if ddr_b==1
keep id_kreis_a ddr_a total total_ddr
drop if total_ddr==.
duplicates drop
gen frac_ddr=total_ddr/total
ren id_kreis_a nuts_id
keep nuts_id frac_ddr
so nuts_id
save ./intermediate/facebook_stub_noself, replace
