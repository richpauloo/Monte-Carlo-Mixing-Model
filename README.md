# Mixing Cell Model of Solute Transport in the Tulare Basin, California  

Created July 7, 2017 by Rich Pauloo at the University of California Davis   


![Groundwater TDS-depth profile across a grid of timesteps is an output of the mixing model.](salinization.gif)  


***  

## Getting Started 

1. Clone this repository  
2. Ensure you have [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/) installed on your machine.  
3. Open MCMM_no_RWI.Rmd in Studio. An R Markdown (`.Rmd`) file is a notebook version of an R file ([more details here](https://rmarkdown.rstudio.com/)), and must be run in RStudio.   
4. Run the code to evaluate the model.  


## Contents

The model depends on 4 input files   
 - `boundary_dat.rds` - initial TDS-depth profile at $t_0$  
 - `GW.csv` - C2VSim 40 year groundwater budget  
 - `LB.csv` - C2VSim 40 year Land Zone budget  
 - `RZ.csv` - C2VSim 40 year Root Zone budget  
 
***  

## Notes  
 - The 40 year period for all data is from 1961-10-31 : 2001-09-30, and begins on October, the start of the water year.  
 - Water budgets are derived from [C2VSim Version 3.02-CG (R374)](http://baydeltaoffice.water.ca.gov/modeling/hydrology/C2VSim/index_C2VSIM.cfm)  

***  

## Contact

Please contact me at **rpauloo at ucdavis dot edu** or **richpauloo at gmail dot com** with any questions.   