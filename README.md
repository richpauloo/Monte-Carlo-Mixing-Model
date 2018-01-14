# Deterministic Monte Carlo Mixing Model of Solute Transport in the Tulare Basin, California

created July 7, 2017 by [Rich Pauloo](richpauloo@gmail.com) at the University of California Davis

## This repository contains 5 files:
 - README.json - a description of the files and how to run the Model
 - MM_working_MC.R - a .R file to be used by the computing language R <https://www.r-project.org/>
 - GW.csv - The 40 year groundwater budget*
 - LB.csv - The 40 year Land Zone Budget*
 - RZ.csv - The 40 year Root Zone Budget*

#### Notes: 
 - The 40 year period for all data is the same (10/31/1961-9/30/2001), and is based on the water year, which begins in October.
 - Results are from C2VSim Version 3.02-CG (R374) <http://baydeltaoffice.water.ca.gov/modeling/hydrology/C2VSim/index_C2VSIM.cfm>

## Running the model:
 - 1. Open MM_working_MC.R in R or R Studio. R Studio is recommended for ease of use and visualization of graphics.
 - 2. Change the working directory in line 11 to this repository, which contains the necessary water budget files: GW.csv, LB.csv, RZ.csv
 - 3. Run the code to visualize output

[To view output for one set of 1000 Monte Carlo samples click here.](http://rpubs.com/richpauloo/mcmm)
