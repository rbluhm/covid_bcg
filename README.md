# README for "The Spread of COVID-19 and the BCG Vaccine" by Richard Bluhm and Maxim Pinkovskiy, 2021, Econometrics Journal

This replication is written for Stata version 15.2 and R version 4.0.3.

It uses the user-written Stata  commands "rdrobust", available via "net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace", and boottest, available via "ssc install  boottest".

Bootstrapping of the RD confidence intervals and the visualizations of spatial data were done in R using the user written packages 

    - "sf", version 0.9.6, 
    - "haven", version 2.3.1,  
    - "tidyverse", version 1.3.0, 
    - "lubridate", version 1.7.9.2, 
    - "RcolorBrewer", version 1.1-2,
    - "rdrobust", version 0.99.9, 
    - "parallel", version 4.0.3, 
    - "foreach", version 1.5.1, 
    - "miscFuncs", version 1.3.

The bootstrap code is based based on boot-rd by Bartalotti et al. (2017, from https://github.com/grayclhn/boot-rd) but has been modified and extended. 

There are three folders.

- /data/ contains the primitive files for the analysis, containing the raw data
- /intermediate/ contains files created during the analysis
- /output/ contains the tables and the figures created during the analysis.

To replicate the entire analysis at once, please run 00_RUN_EVERYTHING.sh (a Unix/OSX shell script which includes the successive calls to R and Stata, on Windows, just run each line of the file calling these programs in the command prompt or install Linux shell support) 

- The Stata files export_data2R.do and export_supplement_data2R.do have to be run prior to running the R-scripts generating the tables for the main text or supplement.
- Each table is replicated by a Stata do-file or R-script beginning with the table number, e.g. Table_2_Placebo_Microcensus.do replicates Table 2 in the main text and Supplement_Table_S4_OtherRDs.R replicates Table S.4. in the Online Appendix. 
- The same goes for the figures in the main text and supplement.
- Each of these files is stand-alone and can be run independently of the others once all the necessary data files have been created. 
- Each file creates the corresponding output (.tex, .pdf, and/or .png)  in the ./output/ folder.

Several do-files exist as utilities or to produce intermediate inputs but do not need to be run independently, instead they are called with the appropriate parameters by the other files:

- death_rate_adjustment.do – performs age-adjustment for death rates.
- hosp_rate_adjustment.do – performs age-adjustment for hospitalization rates.
- simulation_donut.do – computes the dataset of simulated cases.
- facebook.do – computes the social connectedness measure used in the main text.

The .dta files (for every .dta file there is a corresponding .csv ASCII file providing the data without need of purchasing Stata) contain the variables used in the analysis. Below is a description of what each .dta file in the /data/ folder contains (intermediate .dta files produced by the analysis are not described).
 
- all_cause_mortality_by_age_withpop_2016.dta – county-age level data on all cause mortality rates as of 2016.
- cases_pop_cat_Y_singleages_weekX.dta – county by single age bin level data on COVID-19 cases for week X of 2020 (week 17 corresponds to April 26 and week 50 to December 13) for case category Y (Y can be "all" for all cases and "c" for symptomatic cases).
- commuters.dta. commuter_panel_filled_dec2019.dta, commuters.dta and commute_stub.dta – data on county-to-county commuting flows in different formats.
- ddr_border.gpkg – geospatial line features containing the GDR (DDR) border
- export_kreise_nonoverlapping.dta – Microcensus 2017 results for country-by-age group, incl. full labels to the underlying microcensus variables.
- facebook.dta – county-to-county connections from the social connectedness index (SCI) made available by the Facebook Data for Good Program
- fulldata_kreise.dta – county-level characteristics, including cumulative case and death counts as of April 26th, population, covariates such as distance to the border, dates of first case, etc.
- hospitalizations_by_type_and_age_withpop_2016.dta – county-age level data on hospitalizations from all causes, infectious diseases, and respiratory diseases as of 2016.
- infectious_disease_mortality_by_age_withpop_2016.dta – county-age level data on mortality rates from infectious diseases as of 2016.
- initial_cases.dta – county distribution of COVID-19 cases on 29 February 2020.
- kreise_counts_panel_daily.dta – data on new case counts by county and day from Jan 1 2020 until Dec 13 2020.
- kreise_distances_panel.dta – bilateral county-to-county distances.
- nuts_id_stub.dta – helper file with NUTS-3 IDs and numeric IDs.
- respiratory_disease_mortality_by_age_withpop_2016.dta – county-age level data on mortality rates from respiratory diseases as of 2016.
- RKI_Corona_Landkreise.gpkg – geospatial polygon features containing the shape of each German county in 2020.
