# Deterministic Monte Carlo Mixing Model of Solute Transport in the Tulare Basin, California  

created July 7, 2017 by [Rich Pauloo](richpauloo@gmail.com) at the University of California Davis   


~[Groundwater salinization with depth is an output of the mixing model.](salinization.gif)  


***  

## This repository contains 5 files   
 - README.json - a description of the files and how to run the Model  
 - MCMM_no_RWI.Rmd - a .Rmd file to be run with the [computing language R](https://www.r-project.org/)  
 - GW.csv - The 40 year groundwater budget  
 - LB.csv - The 40 year Land Zone Budget  
 - RZ.csv - The 40 year Root Zone Budget  
 
***  

## Notes  
 - The 40 year period for all data is 10/31/1961--9/30/2001, and begins on October, the start of the water year.  
 - Water budgets are dervied from [C2VSim Version 3.02-CG (R374)](http://baydeltaoffice.water.ca.gov/modeling/hydrology/C2VSim/index_C2VSIM.cfm)  

***   

## Running the model  
 - 1. Open MCMM_no_RWI.Rmd (without rock water interactions) in R or R Studio.  
 - 2. Change the working directory in line 11 to your local repository, which contains the necessary water budget files: GW.csv, LB.csv, RZ.csv  
 - 3. Run the code to visualize output  
 
 
