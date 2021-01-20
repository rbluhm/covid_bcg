use ./data/kreise_counts_panel_daily


drop if date>date("13dec2020","DMY")
collapse (sum) new_cases, by(date)
tsset date

tssmooth ma new_cases_ma=new_cases, window(7,1,0)

label var new_cases "New cases"
label var new_cases_ma "7-day avg."

tsline new_cases new_cases_ma, tline("26apr2020", lc(gs6) lp(dash)) ///
	tlabel(, format(%tdm)) legend(pos(11) ring(0) row(2)) ///
	ttext(15000 26apr2020 "26 Apr 2020", placement(east) orientation(vertical))  ///
	tline("12sep2020", lc(gs6) lp(dash)) ///
	ttext(15000 12sep2020 "12 Sep 2020", placement(east) orientation(vertical)) ///
	tline("13dec2020", lc(gs6) lp(dash)) ///
	ttext(15000 13dec2020 "13 Dec 2020", placement(east) orientation(vertical)) ///
	ytitle("New daily cases") ttitle("") graphregion(color(white)) bgcolor(white)

graph export ./output/Supplement_Figure_S3_Time_Series.pdf, replace
