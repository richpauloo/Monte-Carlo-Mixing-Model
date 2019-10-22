####################################################################################
# Deterministic Monte Carlo Mixing Model of ABCSAL
# Rich Pauloo
# 7/4/2017
####################################################################################

# if these packages are not already installed, use the function
# `install.packages()` to install them
library(dplyr)
library(extrafont)
library(scales)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)
library(kableExtra)
library(knitr)
library(xtable)


####################################################################################
# Read data
####################################################################################

# Groundwater budget (10/31/1961-9/30/2001)
GW = read.csv(file = "~/GitHub/Monte-Carlo-Mixing-Model/data/GW.csv", stringsAsFactors = FALSE, header = TRUE) 

# Land Budget budget (10/31/1961-9/30/2001)
LB = read.csv(file = "~/GitHub/Monte-Carlo-Mixing-Model/data/LB.csv", stringsAsFactors = FALSE, header = TRUE)

# Root Zone budget (10/31/1961-9/30/2001)
RZ = read.csv(file = "~/GitHub/Monte-Carlo-Mixing-Model/data/RZ.csv", stringsAsFactors = FALSE, header = TRUE) 

# Bring in boundary condition data and RWI from dissertation/code/02_reanalyze_gw_tds.R 

boundary_dat <- readRDS("~/GitHub/Monte-Carlo-Mixing-Model/data/boundary_dat.rds")
bd <- boundary_dat %>% dplyr::select(x,y) %>% 
  mutate(x = abs(x) * 3.28084) %>% # convert m to ft for model
  arrange(x) %>% 
  as.matrix()


####################################################################################
# Define Common Unit Conversion Functions
####################################################################################

# converts acre-feet to liters
af_to_L = function(a){
  a * 43560 * 28.3168
}

# converts miligrams to million tons
mg_to_Mton = function(a){
  a / (1000*907185*1000000)
}

# converts Mtons to metric metric tons (NOT metric Mtons)
Mtons_to_metric_tons <- function(x){return(x * 907184.74)}

# acre feet to km3
af_to_km3 <- function(a){
  a * 43560 * 28.3168 * 1e-12
}

# liters to km3
l_to_km3 <- function(l){
  l * 1e-12
}


####################################################################################
## Define Boundary Conditions 
####################################################################################

# Rock Water Interaction Coefficients dervied from groundwater 
# salinity literature of the Tulare Basin [@Kang2016; @Williamson1989]. 
# Groundwater velocity was computed with Z Budget and the C2VSim model, 
# Version 3.02-CG (R374) [@Brush2013].

# Rock water interactions
# (Kang, 2016) 3000 mg/L at 1719 ft. in Fresno co., 2205 ft. in Kern co.
Kang       = 3000/mean(c(1719,2205)) 
Williamson = 2000/2500


# velocity field from C2VSim
dv12 = 0.76106  # Darcy velocity between centroids of layers 1-2 (ft/yr)
dv23 = -0.01772 # Darcy velocity between centroids of layers 2-3 (ft/yr)

# C2VSim Groundwater budget (units = AF)
DP     = sum(GW$DP) # deep percolation
Si     = sum(GW$Si) # beginning storage
Sf     = sum(GW$Sf) # ending storage
NDP    = sum(GW$DP) # net deep percolation
Stream = sum(GW$Stream) # stream gain
R      = sum(GW$R) # recharge
Lake   = sum(GW$Lake) # lake gain
BI     = sum(GW$BI) # sum(GW$Si) boundary inflow
Sub    = sum(GW$Sub) # subsidence flow
P      = sum(GW$P) # GW pumping
SI     = sum(GW$Net.Subsurface.Inflow....) # subsurface inflow

# sanity check: average annual Delta S of groundwater system:
af_to_km3(R + Stream + Lake + BI + Sub + SI + NDP - P)/40

# spatial/temporal dimensions
Dim = data.frame(subregion = seq(14,21,1), 
                 # acre feet
                 SA = c(670228.87, 904472.46, 302449.08, 372889.22, 
                        897090.81, 801419.62, 423713.29, 652847.11), 
                 # feet: from C2VSim model
                 thickness = c(1915.155, 1990.01, 1509.194, 1868.821, 
                               1904.244, 1715.117, 1969.524, 3244.262))

# total surface area (acres)
Tot_SA         = sum(Dim$SA)           

# Volume (AF)
V              = Dim$SA * Dim$thickness 

# total volume (L)
Tot_V          = af_to_L(sum(V))        

# surface area weighted thickness
SA_w_thickness = sum(Dim$SA/Tot_SA*Dim$thickness) 

# percent aquifer
#pa = 0.989

# time step length (years)
t = 50

# length of C2VSim groundwater budget (years)
l = 40 


# salt concentration of natural SW and thus GWW from (Cismowski, 2006)
# within report search for: Annual Salt Load (thousand tons/year)
sw <- data.frame(term = c("sac","yolo","sjr"),
                 ttons = c(1945, 405,922),
                 taf_yr = c(16953,2980,3082))

# concentration of delivered sw
sw$c <- sw$ttons / sw$taf_yr

# convert ktons to mg and thousand acre feet to L
sw$c <- sw$c * 9.072e+11 / 1.233e+9 

# volume weighted concetration
sw <- sw %>% 
  mutate(pv  = taf_yr / (sum(sw$taf_yr)), # percent volume
         vwc = pv * c) # volume weighted concentration
# sum(sw$vwc) # sanity check: final concentration of surface water input

# volume weighted concetration
sw <- sw %>% 
  mutate(pv  = taf_yr / (sum(sw$taf_yr)), # percent volume
         vwc = pv * c) # volume weighted concentration

# final concentration of natural surface water and thus GW
gwc <- sum(sw$vwc)

####################################################################################

# C2VSim average annual fluxes from GW budget
annual_fluxes = data.frame(term = c("NDP", "Stream", "R", "Lake", "BI", "Sub", "P", "SI", "P_alt", "M"),
                           "L/yr" = c(af_to_L(NDP/l), af_to_L(Stream/l), af_to_L(R/l), 
                                      af_to_L(Lake/l), af_to_L(BI/l), af_to_L(Sub/l), 
                                      af_to_L(P/l), af_to_L(SI/l), NA, NA))


# eliminate overdraft (delta S = 0)
# first calculate the volume to eliminate: (1) pumping, and 
# (2) subsidence flow, which should not exist in a steady state system)
od  <- sum(annual_fluxes$L.yr[c(1:6,8)]) - sum(annual_fluxes$L.yr[7])
odc <- abs(od) + annual_fluxes$L.yr[6] # overdraft + C
l_to_km3(odc)

# managed aquifer recharge following (Hanak, 2019) has a theoretical
# max of 550 TAF = 0.68 km3. Thus we add M proportional to that amount.
# first we, solve for the proportion of (|dS| + C) that this amounts to
pp <- af_to_km3(550000) / l_to_km3(odc) 
M  <- pp * odc 

# sanity check that M = max recharge from Hanak (2019):
l_to_km3(M)

# new Pumping (P_alt) = R + B + M + I + N 
P_alt <- sum(annual_fluxes[c(1:5,8), 2]) + abs(M)

# percentage reduction in Pumping to stop overdraft,
# used to adjust vertical groundwater velocity
P_alt_prop <- P_alt / annual_fluxes[7,2] 

# percentage increase in M relative to all other recharge terms
# (R,B,I,N) to stop overdraft.
# used to adjust vertical groundwater velocity
#M_prop <- 1 + (M / sum(annual_fluxes$L.yr[c(1:5,8)]))

# sanity check for steady state conditions (no overdraft)
sum(annual_fluxes$L.yr[c(1:5,8)]) + M - P_alt

# add P_alt and M to `layer_fluxes`
annual_fluxes[9,  2] = P_alt
annual_fluxes[10, 2] = M


# calculate flux volume into top layer to enforce steady state conditions
# 50 yr Q from N, R, B, M (R = recharge from streams, lakes, and watersheds)
top_q  <- sum(annual_fluxes[c(1:5, 10), 2]) * 50 

# calculate flux volume into layers 2:m 
# only subsurface inflow because we eliminate overdraft. minimal term
side_q <- annual_fluxes[8,2] * 50


# calculate percentage reduction in vertical groundwater velocity
# which is the ratio of alternate downwards flux to original downwards flux
# alternate downwards flux is given by P_alt
# original  downwards flux is given by P (minus subsidence flow)
# P_alt / (P - C)
v_alt_prop <- P_alt/(annual_fluxes$L.yr[7] - annual_fluxes$L.yr[6])


####################################################################################

# assign TDS to static terms (i.e. - constant mass flux throught simulation)
# note that the TDS of M is the last element in annual_fluxes$TDS
annual_fluxes$TDS = c(0,gwc,gwc,gwc,gwc,gwc,0,gwc,0,gwc) # TDS
annual_fluxes$mass_flux = annual_fluxes$L.yr * annual_fluxes$TDS # mass flux (mg/yr)

# Internal mass flux: 
# only from subsurface inflow because subsidence flow eliminated
internal_mf = sum(annual_fluxes[8, 4]) * t # mg/50 yr 

# Top down Mass Flux: same terms as `top_q`, excluding NDP, 
# the mass of which is dynamically calculated in
# the loop to account for evapoconcentration
top_mf = sum(annual_fluxes[c(2:5,10), 4]) * t # mg/50 yr 

# Root Zone Budget (includes Agriculture, Urban, and Native Vegetation)
I  = sum(RZ$AG_INF) + sum(RZ$URB_INF) + sum(RZ$NAT_INF) # Net infiltration
ET = sum(RZ$AG_ET) + sum(RZ$URB_ET) + sum(RZ$NAT_ET) # Net ET

# Percentages
perc_I   = NDP/I # percentage of infiltration that becomes NDP
#pGWP = P/I # Percentage of groundwater pumping in Total Applied Water
pGWP = (sum(LB$Ag..Pumping) + sum(LB$Urban.Pumping)) / 
  (sum(LB$Ag..Pumping)   + sum(LB$Ag..Diversion) + 
   sum(LB$Urban.Pumping) + sum(LB$Urban.Diversion))
pSW  = 1 - pGWP # percentage of surface water diversions and reuse in Total Applied Water

# Average Annual Surface Water TDS, volume, and mass flux
# TDS_SW = .38*((247.22540+254.26799)/2) + .62*(87.24187) # TDS of surface water

####################################################################################
# salt concentration of SW imports from (Cismowski, 2006)
# within report search for: Annual Salt Load (thousand tons/year)
####################################################################################

# period of record 1985 - 1994
sw <- data.frame(term = c("swp","cvp"),
                 ttons = c(1004, 900),
                 taf_yr = c(2471,2328))

# concentration of delivered sw
sw$c <- sw$ttons / sw$taf_yr
sw$c <- sw$c * 9.072e+11 / 1.233e+9 # convert ktons to mg and thousand acre feet to L

# volume weighted concetration
sw <- sw %>% 
  mutate(pv  = taf_yr / (sum(sw$taf_yr)), # percent volume
         vwc = pv * c) # volume weighted concentration
TDS_SW <- sum(sw$vwc) # final concentration of surface water input

####################################################################################

# sanity check: matches 5.54 km3/yr mass balance in Table S1
V_SW = (sum(LB$Ag..Diversion) + sum(LB$Urban.Diversion))/40 
SW_annual_mass = mg_to_Mton(af_to_L(V_SW) * TDS_SW) # Mton
Mtons_to_metric_tons(SW_annual_mass) / 1e6

## Number of Monte Carlo Samples
###### This number determines how many random samples are drawn from the input parameters modeled as distributions of a random variable. 1000 samples is standard.
N = 1000

## Compute all output  

run_model <- function(irg_eff, RWI_on, N){
  
  # set up model arrays
  TDS                  = array(0, dim=c(8,8,N)) # TDS array
  layer_depths_array   = array(0, dim=c(8,1,N)) # layer depths array
  layer_velocity_array = array(0, dim=c(8,1,N)) # layer velocity array
  mass                 = array(0, dim=c(8,1,N)) # array of sum of masses per time step
  TDS_GWP              = array(0, dim=c(8,1,N)) # array of TDS of pumped groundwater per time step
  TDS_AW               = array(0, dim=c(8,1,N)) # array of TDS of applied water (SW+GW) per time step
  EC                   = array(0, dim=c(8,1,N)) # array of evaporative concentration of NDP per time step (TDS)
  NDP_mf               = array(0, dim=c(8,1,N)) # array of NDP mass flux per time step (mg/yr)
  RWI_cont             = array(0, dim=c(8,1,N)) # array of Rock water interactions
  lay_vol              = array(0, dim=c(8,1,N)) # array of layer volumes
  cbc                  = vector("list", length=N)
  
  # loop over N realizations
  for(k in 1:N){
    
    set.seed(k)
    
    # average regional porosity
    #n = 0.30
    
    # Random Variables
    n    = 0.3 #runif(1, min = porosity[1], max = porosity[2]) # porosity
    pa   = 0.990 # percent aquifer
    aVol = Tot_V * n * pa # saturated aquifer volume (L)
    pI   = runif(1, min = 1-irg_eff[2], max = 1-irg_eff[1]) # irrigation efficiency
    RWI  = runif(1, min = Williamson, max = Kang)  # Rock water interactions (TDS/ft)
    
    
    # velocity field from C2VSim
    v12    = 0.76106/(n*pa) # linear pore velocity between centroids of layers 1-2 (ft/yr)
    v23    = -0.01772/(n*pa) # linear pore velocity between centroids of layers 2-3 (ft/yr)
    vel    = c(v12,v23) # vector of velocities
    #0.5*(mean(dd$c1) + mean(dd$c2)); 0.5*(mean(dd$c2) + mean(dd$c3))
    #c(447,924)#
    depths = c(-631.5,-2077) # C2VSim midpoints between layer centroids over which velcoity is calculated
    vel    = as.data.frame(cbind(vel,depths)) # bind vel and depth vectors into one data frame for analysis
    l      = lm(depths~vel, data = vel) # create a linear model between velocity and depth
    
    
    # vertical velocity profile
    # initalize vectors to hold:
    #   - computed velocities at various depths
    #   - layer thicknesses, computed from groundwater velocity
    v    = matrix(0,8,1) 
    b    = matrix(0,8,1) 
    
    # compute first velocity and first layer thickness
    v[1] = ( (0 + -coef(l)[[1]]) / coef(l)[[2]] ) * v_alt_prop
    b[1] = -(v[1] * t) 
    
    # compute all remaining layer velocities and thicknesses
    for(i in 1:7){
      v[i+1] = ( (sum(b) -coef(l)[[1]]) / coef(l)[[2]] ) * v_alt_prop
      b[i+1] = -(v[i+1] * t) 
    }
    
    
    # layer thickness (ft)
    b[8] = -(SA_w_thickness + sum(b[1:7]))
    pb   = abs(b)/SA_w_thickness # percentage of total thickness per layer
    
    
    # Calculate the layer depths
    layer_depths    = matrix(0,8)
    layer_depths[1] = abs(b[1])
    for(i in 2:8){
      layer_depths[i,1] = abs(layer_depths[i-1,1]) + abs(b[i]) # vector of layer depths (ft)
    }
    
    # layer volumes (L)
    lay_vol[,1,k] = aVol * pb # saturated volume of all layers (L)
    
    
    # Calculate proportion Pumping in each layer (using P_alt)
    # most pumping occurs in the deep aquifer
    lay_pump = matrix(0,8,1)
    for(i in 1:8){
      lay_pump[i,1] = b[i]/sum(b) * annual_fluxes[9,2] 
    }
    
    # repeat for layer subsurface inflow (proportional `side_q`) 
    # in cells 2:m
    lay_si = matrix(0,8,1)
    lay_si[1] = NA
    for(i in 2:8){
      lay_si[i,1] = b[i]/sum(b[2:length(b)]) * annual_fluxes[8,2] 
    }
    
    
    # Vertical Flux per layer (L/yr)
    # based on velocity alone, these fluxes would occur with NO OTHER
    # inputs and outputs to the system
    Q = v * n * pa * Tot_SA * 43560 * 28.3168 # L/yr
    layer_fluxes = cbind(v, Q) # units of "ft/yr" and "L/yr" respectively
    
    # to enforce steady state conditions, the layer to layer fluxes 
    # are the residual of all other fluxes between the layers
    # for layer 1, Q_{in} from the top is `top_q` and must equal Q_{out},
    # which = P + Q_{1,2}. Q_{1,2} is the closure term.
    # Thus Q_{1,2} is Q_{in} - P_{1} (pumping in layer 1):
    q12 <- top_q - (lay_pump[1] * t)      # 50 year in minus out
    z <- vector()                         # vector of layer to layer fluxes
    z[1] <- q12                           # initalize with Q_{1,2}
    for(i in 2:8){                        # solve for Q_{m-1,m}
      z[i] <- (z[i-1] + lay_si[i] * t) -  # in minus...
        (lay_pump[i] * t)                 # out, over 50 yr time step
    }
    
    # store layer to layer fluxes
    layer_fluxes <- cbind(layer_fluxes, z)
    
    # arrange table of layer to layer fluxes
    cbc[[k]] <- data.frame(p = -lay_pump, 
                           qc = z/t, 
                           r = c(sum(annual_fluxes[2:4, 2]), rep(NA, 7)),
                           n = c(annual_fluxes[1,2], rep(NA, 7)),
                           b = c(annual_fluxes[5,2], rep(NA, 7)),
                           m = c(annual_fluxes[10,2], rep(NA, 7)),
                           i = lay_si) %>% 
      mutate_all(l_to_km3)
    
    # Initalize TDS array with baseline TDS
    #t0 = layer_depths * RWI # inital TDS-depth profile from RWI
    
    bc_match <- matrix(0, nrow = 8, ncol = 1)
    for(i in 1:nrow(layer_depths)){
      temp  <- layer_depths[i,1] - bd[,1]
      index <- which.min(abs(temp))
      bc_match[i] = bd[index,2]
    }
    
    t0 = bc_match
    TDS[,1,k] = t0
    
    # RWI == FALSE : no rock-water interactions
    if(RWI_on == FALSE){RWI_cont[,1,k] = 0}
    
    # # RWI == TRUE
    if(RWI_on == TRUE){
      RWI_cont[,1,k] = abs(b) * RWI # add RWI cont proportional to volume, depth
    }
    
    # Record layer depth/velocity for all simulations
    layer_depths_array[,,k]   = layer_depths
    layer_velocity_array[,,k] = v
    
    # layer basline mass
    lbm = bc_match * lay_vol[,,k]
    
    # Compute TDS for all layers and times
    for(j in 1:7){
      
      # calculate mass flux of NDP
      mass[j,1,k]    = sum(TDS[1:7,j,k] * lay_vol[1:7,1,k]) # sum of mass in pumping layers 1-7 (mg)
      TDS_GWP[j,1,k] = mass[j,1,k]/sum(lay_vol[1:7,1,k]) # TDS of pumped GW
      TDS_AW[j,1,k]  = TDS_GWP[j,1,k] * pGWP + TDS_SW * pSW # TDS of applied water for irrigation
      EC[j,1,k]      = TDS_AW[j,1,k] * (1/pI) # TDS of NDP Evaporative concentration
      NDP_mf[j,1,k]  = EC[j,1,k] * annual_fluxes[1,2] # mass flux of NDP
      
      # TDS in layer 1 is...
      # initial mass, PLUS
      # top down mass flux, PLUS
      # ...
      # mass from NDP, MINUS
      # mass flux to layer below (q_{1,2}), MINUS
      # Pumping (P_alt)... ALL DIVIDED BY
      # the Volume of layer 1   
      TDS[1,j+1,k] = (
        (
          TDS[1,j,k]    * lay_vol[1,1,k] +       
          top_mf +                         
          #internal_mf   * pb[1] +              
          NDP_mf[j,1,k] * t -                 
          TDS[1,j,k]    * layer_fluxes[1,3]  -  
          TDS[1,j,k]    * lay_pump[1,1] * t
          ) /
          lay_vol[1,1,k]
      )                  
      
      # TDS in layers 2-7 is...
      # initial mass, PLUS
      # internal mass flux, PLUS  
      # mass flux from layer above, MINUS  
      # Pumping, MINUS
      # mass flux into layer below, ALL DIVDED BY  
      # the layer volume
      for(i in 2:7){
        TDS[i,j+1,k] = (
          (
            TDS[i,j,k]   * lay_vol[i,1,k] +         
            internal_mf  * pb[i] +                 
            TDS[i-1,j,k] * layer_fluxes[i-1,3] -   
            TDS[i,j,k]   * lay_pump[i,1] * t -        
            TDS[i,j,k]   * layer_fluxes[i,3]
            ) / 
            lay_vol[i,1,k]
        )                     
      }
      
      # baseline correction from rock water interactions 
      mean_bl <- mean((TDS[,j,k] - TDS[,j+1,k])[(j+1):7])
      
      TDS[,j+1,k][(j+1):7] <- TDS[,j+1,k][(j+1):7] + (TDS[,j,k] - TDS[,j+1,k])[(j+1):7]
      TDS[1:j,j+1,k] <- TDS[1:j,j+1,k] + mean_bl
      
      # when RWI == TRUE, add the RWI zero order source term
      if(RWI_on == TRUE){
        TDS[,j+1,k] <- TDS[,j+1,k] + RWI_cont[,,k] 
      }
      
    }
    
  }
  
  # compile model arrays to output
  res <- list(TDS                  = TDS,
              layer_depths_array   = layer_depths_array,
              layer_velocity_array = layer_velocity_array,
              mass                 = mass,
              TDS_GWP              = TDS_GWP,
              TDS_AW               = TDS_AW,
              EC                   = EC,
              NDP_mf               = NDP_mf,
              RWI_cont             = RWI_cont,
              lay_vol              = lay_vol,
              cbc                  = cbc) 
  
  return(res)
  
}



##########################################################################
# evaluate the model with and without RWI
##########################################################################
# with RWI
z1 <- run_model(irg_eff = c(0.778, 1-perc_I),
                RWI_on = TRUE, 
                N = 1000)

# without RWI
z2 <- run_model(irg_eff = c(0.778, 1-perc_I),
                RWI_on = FALSE, 
                N = 1000)


##########################################################################
# cell by cell budget. sanity check - steady state
# since eta and f are fixed, so also is b, and thus all cbc budgets are =
# thus we take the first result
# the code is set up this way in case we want to vary the parameters that
# control aquifer thickness (eta, f)

cbc <- z1$cbc[[1]]

# layer 1
cbc1 <- tibble(term = colnames(cbc),
           q    = unlist(slice(cbc, 1), use.names = FALSE)) %>% 
  slice(-7) %>% 
  mutate(q = ifelse(term == "qc", q * -1, q))

# sanity check: steady state
pull(cbc1, q) %>% sum()
cbc1$term <- c("P{alt, 1}", "q{1,2}", "R", "N", "B", "M")
cbc1 <- cbc1[c(4,3,5,6,1,2), ]
cbc1 <- bind_rows(cbc1, data.frame(term = "dS", q = sum(cbc1$q)))

# layers 2:m
cbc_in <- select(cbc, p, qc, i)
cbc2 <- vector("list", length=6)
for(i in 1:length(cbc2)){
  counter <- i + 1
  cbc_i <- tibble(term = c(colnames(cbc_in), "qc_in"),
         q    = c(unlist(slice(cbc_in, counter), use.names = FALSE),
                  cbc_in$qc[i])) 
  cbc_i$q[2] <- cbc_i$q[2] * -1
  cbc_i <- bind_rows(cbc_i, data.frame(term = "dS", q = sum(cbc_i$q)))
  
  cbc_i$term <- c(paste0("P{alt, ", counter, "}"), 
                  paste0("q{", counter, ",", counter + 1, "}"),
                  paste0("I", counter),
                  paste0("q{", i, ",", counter, "}"),
                  "dS")
  cbc_i <- cbc_i[c(4,3,1,2,5), ]
  
  cbc2[[i]] <- cbc_i
}

print(
  xtable(
    bind_rows(cbc1, 
              bind_rows(cbc2)
              ),
    digits = 3, 
    display = c("d","s","e")
    ),
  include.rownames = FALSE
  )





##########################################################################
## Plot 1: time, layer depth, TDS, simulation df

# gather appraoch 


# make depth-time array
ld   <- lapply(1:1000, function(x){return(z1$layer_depths_array[1:7,,x])}) # all layer depths in a list
td   <- lapply(1:1000, function(x){return(z1$TDS[1:7,,x])}) # tds as a list of matrices
ldtd <- lapply(1:1000, function(x){return(data.frame(cbind(ld[[x]], td[[x]])))}) # bind

for(i in 1:1000){
  colnames(ldtd[[i]]) <- c("d", "0", "50","100","150","200","250","300","350")
  ldtd[[i]] <- gather(ldtd[[i]], time, tds, -d) # wide to long
  ldtd[[i]]$sim <- i # add simulation
}

df <- bind_rows(ldtd)
df$sim <- as.factor(df$sim)
df$time <- factor(paste0("t = ", df$time), 
                  levels = paste0("t = ", 
                                  c("0", "50","100","150","200","250","300","350")))
df$d <- df$d * 0.3048       # ft to meters
df$rwi <- "RWI"




# make depth-time array
ld   <- lapply(1:1000, function(x){return(z2$layer_depths_array[1:7,,x])}) # all layer depths in a list
td   <- lapply(1:1000, function(x){return(z2$TDS[1:7,,x])}) # tds as a list of matrices
ldtd <- lapply(1:1000, function(x){return(data.frame(cbind(ld[[x]], td[[x]])))}) # bind

for(i in 1:1000){
  colnames(ldtd[[i]]) <- c("d", "0", "50","100","150","200","250","300","350")
  ldtd[[i]] <- gather(ldtd[[i]], time, tds, -d) # wide to long
  ldtd[[i]]$sim <- i # add simulation
}

df2 <- bind_rows(ldtd)
df2$sim <- as.factor(df2$sim)
df2$time <- factor(paste0("t = ", df2$time), 
                   levels = paste0("t = ", 
                                   c("0", "50","100","150","200","250","300","350")))
df2$d <- df2$d * 0.3048       # ft to meters
df2$rwi <- "noRWI"




ll <- c(
  't = 0'  ="t = 0 yrs (1960)",
  't = 50' ="t = 50 yrs (2010)",
  't = 100'="t = 100 yrs (2060)",
  't = 150'="t = 150 yrs (2110)",
  't = 200'="t = 200 yrs (2160)",
  't = 250'="t = 250 yrs (2210)",
  't = 300'="t = 300 yrs (2260)",
  't = 350'="t = 350 yrs (2310)"
)



# p <- df %>% filter(time %in% c("t = 0", "t = 50", "t = 100", "t = 200", "t = 250", "t = 300") & sim %in% 251:500) %>%
#   ggplot(aes(-d, tds)) +
#   geom_line(aes(color = sim), alpha = 0.2) +
#   geom_smooth(se = FALSE, color = "red") +
#   coord_flip(ylim = c(0, 7500)) +
#   scale_color_grey() +
#   guides(color = FALSE) +
#   theme_minimal() +
#   facet_wrap(~time, labeller = as_labeller(ll)) +
#   scale_y_continuous(breaks = c(0, 1000, 2500, 5000, 7500),
#                      labels = c('0', '1,000', '2,500', '5,000', '7,500')) +
#   labs(x = "Depth (m)", y = "TDS (mg/L)") +
#   theme(panel.grid.minor = element_blank())
# 
# p

#ggsave(p, filename = "results/p_sim.pdf", device = cairo_pdf, height = 5, width = 7)



# df == RWI
df$tc <- rep(1:8, each = 7) # add time class for easy filtering (1 = 0 yrs, 2 = 50 yrs...)
df$dc <- 1:7                # depth class for easy filtering (1 = layer 1, 2 = layer 2, ...)

filter(df, dc == 1, tc == 1) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))
filter(df, dc == 1, tc == 2) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))
filter(df, dc == 1, tc == 5) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))

# df2 == no RWI
df2$tc <- rep(1:8, each = 7)
df2$dc <- 1:7

filter(df2, dc == 1, tc == 1) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))
filter(df2, dc == 1, tc == 2) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))
filter(df2, dc == 1, tc == 5) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))
filter(df2, dc == 1, tc == 7) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))

filter(df2, dc == 4, tc == 5) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))

filter(df2, dc == 7, tc == 7) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75))


# tables
t1 <- group_by(df2, dc, time) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75)) %>% 
  ungroup() %>% 
  select(-dc) %>% 
  filter(time != "t = 350") %>% 
  mutate_at(.vars = c("mean_depth","p50","p25","p75"), .funs = round) %>% 
  mutate(IQR = paste0(p25, "-", p75)) %>% 
  select(-c("p25", "p75")) %>% 
  rename(median = p50) %>% 
  mutate(mean_depth = as.character(mean_depth),
         median = as.character(median))

t2 <- group_by(df, dc, time) %>% 
  summarise(mean_depth = mean(d), p50 = median(tds), 
            p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75)) %>% 
  ungroup() %>% 
  select(-dc) %>% 
  filter(time != "t = 350") %>% 
  select(p50:p75) %>%   
  mutate_at(.vars = c("p50","p25","p75"), .funs = round) %>% 
  mutate(IQR = paste0(p25, "-", p75)) %>% 
  select(-c("p25", "p75")) %>% 
  rename(median = p50) %>% 
  mutate(median = as.character(median))

t3 <- cbind.data.frame(t1, t2)
t3$time <- stringr::str_sub(t3$time, 5,10)

print(xtable(t3), tabluar.environment = "longtable", include.rownames = FALSE)




# Combine RWI and no RWI

#df1 <- df # RWI
#df2 <- filter(df2, time != "t = 0") # no RWI: second run
#df1$rwi <- "RWI"; df2$rwi <- "no RWI"

df3 <- bind_rows(df, df2) %>% mutate(sim2 = paste0(as.character(sim), "_", rwi))

p <- df3 %>% 
  filter(time %in% c("t = 50", "t = 100", "t = 200", "t = 250", "t = 300") & 
           sim %in% 1:250) %>% 
  group_by(dc, time, rwi) %>% 
  summarise(md = median(d), p25 = quantile(tds, 0.25), p75 = quantile(tds, 0.75)) %>% 
  ggplot(aes(-md)) +
  geom_ribbon(aes(ymin = p25, ymax=p75, group = rwi), fill = "grey70", alpha = 0.8) + 
  #ggplot(aes(-d, tds)) + 
  #geom_line(aes(-d, tds, color = sim2), alpha = 0.2, lwd = 0.4) +
  # RWI == TRUE: purple
  geom_line(data = filter(df,
                            time %in% c("t = 0", "t = 50", "t = 100", "t = 200", "t = 250", "t = 300")) %>% 
                group_by(dc, time) %>% 
                summarise(mtds = median(tds),
                          md   = median(d)), 
              aes(-md, mtds),
              color = "#440154ff", 
            lwd = 0.8) +
  # RWI == FALSE: blue
  geom_line(data = filter(df2,
                          time %in% c("t = 0", "t = 50", "t = 100", "t = 200", "t = 250", "t = 300")) %>% 
              group_by(dc, time) %>% 
              summarise(mtds = mean(tds),
                        md   = median(d)) , 
              aes(-md, mtds), 
              color = "#21908dff", 
            lwd = 0.8) +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  scale_color_grey() + 
  guides(color = FALSE) + 
  coord_flip(ylim = c(-0, 3500))+#21000)) + 
  theme_minimal() + 
  facet_wrap(~time, labeller = as_labeller(ll)) +
  # scale_y_continuous(breaks = c(0, 5000, 10000),#, 15000, 20000), 
  #                    labels = c('0', '5,000', '10,000'))+#, '15,000', '20,000')) + 
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000), 
                     labels = c('0', '1,000', '2,000', '3,000')) + 
  scale_x_continuous(breaks = c(-50, -100, -150, -200, -250), 
                     labels = c('-50', '-100', '-150', '-200', '-250')) + 
  labs(x = "Depth (m)", y = "TDS (mg/L)")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.75, "cm"))
 
p
ggsave(p, filename = "~/GitHub/Monte-Carlo-Mixing-Model/results/p_sim_both2.pdf", device = cairo_pdf, height = 5, width = 7)



## Plot 2: Compare concentrations of GW, TAW, NDP
# initalize matricies
GWP_mean  = as.data.frame(matrix(0,7,1))
AW_mean   = as.data.frame(matrix(0,7,1))
EC_mean   = as.data.frame(matrix(0,7,1))
EC_median = as.data.frame(matrix(0,7,1))
GWP_5     = as.data.frame(matrix(0,7,1))
AW_5      = as.data.frame(matrix(0,7,1))
EC_5      = as.data.frame(matrix(0,7,1))
GWP_95    = as.data.frame(matrix(0,7,1))
AW_95     = as.data.frame(matrix(0,7,1))
EC_95     = as.data.frame(matrix(0,7,1))

for(i in 1:7){
  GWP_mean[i,1]  = mean(z1$TDS_GWP[i,1,])
  AW_mean[i,1]   = mean(z1$TDS_AW[i,1,])
  EC_mean[i,1]   = mean(z1$EC[i,1,])
  EC_median[i,1] = median(z1$EC[i,1,])
  GWP_5[i,1]     = quantile(z1$TDS_GWP[i,1,], c(0.25))
  GWP_95[i,1]    = quantile(z1$TDS_GWP[i,1,], c(0.75))
  AW_5[i,1]      = quantile(z1$TDS_AW[i,1,], c(0.25))
  AW_95[i,1]     = quantile(z1$TDS_AW[i,1,], c(0.75))
  EC_5[i,1]      = quantile(z1$EC[i,1,], c(0.25))
  EC_95[i,1]     = quantile(z1$EC[i,1,], c(0.75))
}

## compile stats for TDS of pumped GW, Total Applied Water (GW+SW), and Evaporative concentration of NDP

# no RWI
EC_dat_rwi = data.frame("legend" = rep(c("A","B","C"), each=7),
                        "time"   = factor(seq(0,300,50)),
                        "mean"   = rbind(GWP_mean, AW_mean, EC_median),
                        "p5"     = rbind(GWP_5, AW_5, EC_5),
                        "p95"    = rbind(GWP_95, AW_95, EC_95))
colnames(EC_dat_rwi) = c("legend","time","mean","p5", "p95")
EC_dat_rwi$class <- "Rock-water interactions"




# initalize matricies
GWP_mean  = as.data.frame(matrix(0,7,1))
AW_mean   = as.data.frame(matrix(0,7,1))
EC_mean   = as.data.frame(matrix(0,7,1))
EC_median = as.data.frame(matrix(0,7,1))
GWP_5     = as.data.frame(matrix(0,7,1))
AW_5      = as.data.frame(matrix(0,7,1))
EC_5      = as.data.frame(matrix(0,7,1))
GWP_95    = as.data.frame(matrix(0,7,1))
AW_95     = as.data.frame(matrix(0,7,1))
EC_95     = as.data.frame(matrix(0,7,1))

for(i in 1:7){
  GWP_mean[i,1]  = mean(z2$TDS_GWP[i,1,])
  AW_mean[i,1]   = mean(z2$TDS_AW[i,1,])
  EC_mean[i,1]   = mean(z2$EC[i,1,])
  EC_median[i,1] = median(z2$EC[i,1,])
  GWP_5[i,1]     = quantile(z2$TDS_GWP[i,1,], c(0.25))
  GWP_95[i,1]    = quantile(z2$TDS_GWP[i,1,], c(0.75))
  AW_5[i,1]      = quantile(z2$TDS_AW[i,1,], c(0.25))
  AW_95[i,1]     = quantile(z2$TDS_AW[i,1,], c(0.75))
  EC_5[i,1]      = quantile(z2$EC[i,1,], c(0.25))
  EC_95[i,1]     = quantile(z2$EC[i,1,], c(0.75))
}

# RWI: rerun with RWI
EC_dat_nrwi = data.frame("legend"=rep(c("A","B","C"), each=7),
                         "time" = factor(seq(0,300,50)),
                         "mean" = rbind(GWP_mean, AW_mean, EC_median),
                         "p5" = rbind(GWP_5, AW_5, EC_5),
                         "p95" = rbind(GWP_95, AW_95, EC_95))
colnames(EC_dat_nrwi) = c("legend","time","mean","p5", "p95")
EC_dat_nrwi$class <- "No rock-water interactions"
#EC_dat_nrwi$class <- ""

# combine
EC_dat2 <- bind_rows(EC_dat_nrwi, EC_dat_rwi)

# plot
p2 <- ggplot(filter(EC_dat2, time %in% c(0,50,100,150)), aes(x=time, y=mean, fill=legend)) +
  geom_bar(position=position_dodge(), stat="identity", 
           colour = "black", size=.3,
           alpha = 0.7) +
  geom_errorbar(aes(ymin=p5, ymax=p95),
                size=.3,    # thinner lines
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  xlab("Time (yrs)") +
  ylab("TDS (mg/L)") +
  scale_y_continuous(labels = comma) +
  scale_fill_viridis_d(name="Source",
                       breaks=c("A","B","C"),
                       labels=c("Pumped Groundwater", 
                                "Total Applied Water", 
                                "Net Deep Percolation")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = c(0.15, 0.77), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = "white", color = "transparent"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        legend.key.size = unit(.5, "cm")) +
  coord_cartesian(ylim = c(0,7100)) +
  scale_y_continuous(breaks = c(0,1000,3000,5000,7000), 
                     labels = formatC(c(0,1000,3000,5000,7000), big.mark=",")) +
  facet_wrap(~class, ncol = 2)

p2
# save
ggsave("~/GitHub/Monte-Carlo-Mixing-Model/results/p_ec2.pdf", p2, height= 4, width = 7, device = cairo_pdf)




# Plot 3: Salt Budget for Water Budget terms and NDP

## non:constant mass: extract all GW pumping TDS and convert to mass
gwp <- gwp2 <- vector("list", length = 7)

for(i in 1:7){
  gwp[[i]]  <- data.frame(tds = z1$TDS_GWP[i, 1, ], t = i, type = "GWP_RWI")
  gwp2[[i]] <- data.frame(tds = z2$TDS_GWP[i, 1, ], t = i, type = "GWP_nRWI")
}

gwp  <- dplyr::bind_rows(gwp)
gwp2 <- dplyr::bind_rows(gwp2)

gwp <- gwp %>% 
  mutate(m = Mtons_to_metric_tons(mg_to_Mton(tds * annual_fluxes[7,2])))

gwp2 <- gwp2 %>% 
  mutate(m = Mtons_to_metric_tons(mg_to_Mton(tds * annual_fluxes[7,2])))

# constant mass: surface water and allother budget terms

# surface water
swm <- Mtons_to_metric_tons(SW_annual_mass)

# all other budget terms
all_else <- Mtons_to_metric_tons(mg_to_Mton(sum(annual_fluxes[c(2,3,4,5,6,8),4])))

# combine into one df
df1 <- data.frame(m = rep(c(swm, all_else), each = 7), 
                  t = rep(1:7, times = 2), 
                  type = rep(c("surface_water_import", "all_else"), each = 7)) 
df2 <- select(gwp, m, t, type) %>% 
  group_by(t) %>% 
  summarise(median = median(m), p5 = quantile(m, 0.25), p95 = quantile(m, 0.75))

df2b <- select(gwp2, m, t, type) %>% 
  group_by(t) %>% 
  summarise(median = median(m), p5 = quantile(m, 0.25), p95 = quantile(m, 0.75))

df1$p5  <- NA
df1$p95 <- NA
df2  <- df2 %>% rename(m = median) %>% mutate(type = "gwp_rwi") %>% select(m, t, type, p5, p95)
df2b <- df2b %>% rename(m = median) %>% mutate(type = "gwp_no_rwi") %>% select(m, t, type, p5, p95)

# rwi: need to go back and re-run code with RWI
mg_to_metric_tons <- function(x){return(x * 1e-9)}

#RWI_cont[7,1,] <- RWI_cont[1,1,]
RWI_int <- mg_to_metric_tons(colSums((z1$RWI_cont[1:7,1,]*z1$lay_vol[1:7,1,])/50)) %>% 
  quantile(c(0.25, 0.5, 0.75))
df3 <- data.frame(m = RWI_int[2], t = 1:7, type = "rwi", p5 = RWI_int[1], p95 = RWI_int[3])

df4 <- rbind.data.frame(df1, df2, df3)
df4$t <- rep(seq(0,300, 50), times = 4) # format times
df4$type <- factor(rep(c("Surface Water Diversions", "I, C, R, B", "Pumped Groundwater", "Rock-water Interactions"), each = 7), # format labels
                   levels = c("I, C, R, B", "Surface Water Diversions", "Pumped Groundwater","Rock-water Interactions"))

# RWI
df4_rwi <- df4
df4_rwi$class <- "Rock-water interactions"

# no RWI: re-run without RWI
df3b <- df3
df3b$m <- 0; df3b$p5 <- 0; df3b$p95 <- 0 
df4 <- rbind.data.frame(df1, df2b, df3b)
df4$t <- rep(seq(0,300, 50), times = 4) # format times
df4$type <- factor(rep(c("Surface Water Diversions", "I, C, R, B", "Pumped Groundwater", "Rock-water Interactions"), each = 7), # format labels
                   levels = c("I, C, R, B", "Surface Water Diversions", "Pumped Groundwater","Rock-water Interactions"))
df4_nrwi <- df4 #%>% filter(m != 0)
df4_nrwi$class <- "No rock-water interactions"

# bind
df5 <- bind_rows(df4_rwi, df4_nrwi)

# plot
p3 <- ggplot(filter(df5, t %in% c(0,50,100,150)), aes(factor(t), m/1000000, fill = type)) + # convert m from tons to Mtons
  #p3 <- ggplot(filter(df4_rwi, t %in% c(0,50,100,150)), aes(factor(t), m/1000000, fill = type)) + # convert m from tons to Mtons
  geom_bar(position = position_dodge(), stat = "identity", 
           colour = "black", size=.3,
           alpha = 0.7) +
  geom_errorbar(aes(ymin = p5/1000000, ymax = p95/1000000),
                size=.3,    # thinner lines
                width=.2,   # Width of the error bars
                position=position_dodge(.9)) +
  scale_fill_viridis_d("Source", option = "E") +
  theme_minimal(base_size = 13) +
  theme(legend.position = c(0.16, 0.74), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = "white", color = "transparent"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        legend.key.size = unit(.5, "cm")) +
  #theme(legend.position = "bottom") +
  labs(x = "Time (yrs)", y = "Annual mass (Metric Mtons)") +
  facet_wrap(~class, ncol = 2) +
  scale_y_continuous(breaks = seq(0,10,2), labels = as.character(seq(0,10,2))) +
  coord_cartesian(ylim = c(0,10))

p3
ggsave("~/GitHub/Monte-Carlo-Mixing-Model/results/p_salt_budget2.pdf", p3, height= 4, width = 7, device = cairo_pdf)

##########################################################################
# combine plots 2 and 3
p3b <- p3 + labs(x=NULL) 
p2b <- p2 + theme(strip.background = element_blank(),
                  strip.text = element_blank())
p23 <- cowplot::plot_grid(p3b, p2b, align = "hv", ncol = 1, labels = c("(A)","(B)"))
ggsave("~/GitHub/Monte-Carlo-Mixing-Model/results/p23.pdf", p23, height= 8, width = 7, device = cairo_pdf)



##########################################################################
# GIF for github
# mean TDS
TDS_mean = as.data.frame(matrix(0,8,8))
for(i in 1:8){
  for(j in 1:8){
    TDS_mean[i,j] = mean(z2$TDS[i,j,])
  }
}

# mean depths
sd.top    = matrix(0,7,8)
sd.bottom = matrix(0,7,8)
for(d in 1:7){
  for(t in 1:8){
    sd.top[d,t]    = TDS_mean[d,t] - sd(z2$TDS[d,t,]) # sd
    sd.bottom[d,t] = TDS_mean[d,t] + sd(z2$TDS[d,t,]) # sd
  }
}

# 5th and 95th percentiles
TDS_5  = matrix(0,7,8)
TDS_95 = matrix(0,7,8)
for(d in 1:7){
  for(t in 1:8){
    TDS_5[d,t]  = quantile(z2$TDS[d,t,], c(0.25)) # 5% CI
    TDS_95[d,t] = quantile(z2$TDS[d,t,], c(.75))  # 95% CI
  }
}

## calculate layer depths/velocities of simulation means for main plot
layer_depths_mean   = as.data.frame(matrix(0,7,1))
layer_velocity_mean = as.data.frame(matrix(0,7,1))
for(j in 1:7){
  layer_depths_mean[j,1]   = mean(z2$layer_depths_array[j,1,]) * -1
  layer_velocity_mean[j,1] = mean(z2$layer_velocity_array[j,1,])
}

## reshape layer_depths_mean for plotting
layer_depths_mean = t.default(layer_depths_mean)

## trim 8th layer of TDS_mean matrix for plotting
TDS_mean = TDS_mean[1:7,]


# compile into one df
pl = data.frame(depth = layer_depths_mean, TDS1 = TDS_mean[,1], TDS2 = TDS_mean[,2],
                TDS3 = TDS_mean[,3], TDS4 = TDS_mean[,4], TDS5 = TDS_mean[,5],
                TDS6 = TDS_mean[,6], TDS7 = TDS_mean[,7], TDS8 = TDS_mean[,8],
                U1 = TDS_95[,1], L1 = TDS_5[,1],
                U2 = TDS_95[,2], L2 = TDS_5[,2],
                U3 = TDS_95[,3], L3 = TDS_5[,3],
                U4 = TDS_95[,4], L4 = TDS_5[,4],
                U5 = TDS_95[,5], L5 = TDS_5[,5],
                U6 = TDS_95[,6], L6 = TDS_5[,6],
                U7 = TDS_95[,7], L7 = TDS_5[,7],
                U8 = TDS_95[,8], L8 = TDS_5[,8])

# make gif
library(gganimate)
temp <- select(pl, depth, TDS1:TDS8)
temp2 <- tidyr::gather(temp, key, value, - depth)

# data reshaping
temp2 <- temp2 %>% 
  mutate(t = as.numeric(substr(key, 4, 5))) 
temp3 <- select(pl, depth, starts_with("U")) %>%
  tidyr::gather(key, value, - depth) %>% 
  mutate(t = as.numeric(substr(key, 2, 3))) %>% 
  rename(p95 = value) %>% 
  select(-key)
temp4 <- select(pl, depth, starts_with("L")) %>% 
  tidyr::gather(key, value, - depth) %>% 
  mutate(t = as.numeric(substr(key, 2, 3))) %>% 
  rename(p5 = value) %>% 
  select(-key)
temp5 <- left_join(temp3, temp4, by = c("depth", "t")) 
temp6 <- left_join(temp5, temp2, by = c("depth", "t"))

temp6$depth <- temp6$depth * 0.3048 # convert depth to m
temp6$t     <- rep(seq(0, 350, 50), each = 7) # time from index to years
temp6       <- filter(temp6, t <= 300) 

# make the animation
anim <- ggplot(temp6, aes(depth, value)) + 
  geom_ribbon(aes(ymin = p5, ymax = p95), alpha = 0.4) +
  geom_line(col = "red", size = 1) + 
  #geom_hline(yintercept = 1000, linetype = "dashed") +
  coord_flip() +
  #annotate(geom = "text", x = -250, y = 1600, label = "TDS upper MCL") +
  #geom_curve(data = data.frame(x1=-260,x2=-300,y1=1600,y2=1050),
  #aes(x = x1, y = y1, xend = x2, yend = y2),
  #arrow = arrow(length = unit(0.03, "npc")),
  #curvature = -0.4) +
  scale_y_continuous(labels = comma, limits = c(0,2500)) +
  theme_minimal(base_size = 16) +
  labs(title = 'Time: {round(frame_time, 2)} yrs ({1960 + round(frame_time, 0)})', 
       y = "TDS (mg/L)", x = "Depth (m)") + 
  transition_time(t) + 
  ease_aes("linear") +
  theme(panel.grid.minor = element_blank())

anim_save("~/GitHub/Monte-Carlo-Mixing-Model/results/salinization.gif", anim) # save to root




