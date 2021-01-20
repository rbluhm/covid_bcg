#!/bin/bash
  
# This bash scripts run the complete replication kit for Bluhm and Pinkovskiy 2021 "The Spread of COVID-19 and the BCG Vaccine", requires Linux/OSX or a Linux subsystem in Windows

mkdir intermediate
echo "Creating intermediate data files with Stata ..."
stata-se -b export_data2R.do
stata-se -b export_supplement_data2R.do

mkdir output
echo "Running main figures with R and Stata ..."
R CMD BATCH Figure_1_map_cases.R
stata-se -b Figure_2_a_c_e.do
stata-se -b Figure_2_b_d_f.do

echo "Running main RD-DD tables with Stata ..."
stata-se -b Table_2_Placebo_Microcensus.do
stata-se -b Table_3_Placebo_Health.do
stata-se -b Table_4_RD_DD_1974.do
stata-se -b Table_5_RD_DD_1990.do


read -p "Do you want to run the main RD tables with R (*bootstrapping takes hours*)? [Y/y]es or [N/n]o " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Running main RD tables with R (*bootstrapping takes hours*) ..."
    R CMD BATCH Table_6_Cases_RDs.R
    R CMD BATCH Table_7_Deaths_RDs.R
fi

echo "Running supplement figures with R and Stata ..."
R CMD BATCH Supplement_Figure_S1_Flows.R
R CMD BATCH Supplement_Figure_S2_Map.R
stata-se -b Supplement_Figure_S3_Time_Series.do
stata-se -b Supplement_Figure_S4_Cases_by_Age.do

echo "Running supplement RD-DD and OLS tables with Stata ..."
stata-se -b Supplement_Table_S1_Placebo_Microcensus_MSERD.do
stata-se -b Supplement_Table_S2_Placebo_Health_MSERD.do
stata-se -b Supplement_Table_S3_RD_DD_MSERD.do
stata-se -b Supplement_Table_S7_OLS.do

read -p "Do you want to run the supplement RD tables with R (*bootstrapping takes hours*)? [Y/y]es or [N/n]o " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Running supplement RD tables with R ..."
    R CMD BATCH Supplement_Table_S4_OtherRDs.R
    R CMD BATCH Supplement_Table_S5_OtherRDsII.R
    R CMD BATCH Supplement_Table_S6_RDs_controls.R
fi


