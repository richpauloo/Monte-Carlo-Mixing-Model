# =================================================================================================== #
#                                                                                                     #
#                   Deterministic Monte Carlo Mixing Tank Model of Solute Transport                   #
#                                         Rich Pauloo, 2017                                           #
#                                                                                                     #
# =================================================================================================== #

# ======================================== #
# Load Data                                #
# ======================================== #
setwd("/Users/richpauloo/Desktop/School and Work/R/MM")
GW = read.csv(file = "GW.csv", stringsAsFactors = FALSE, header = TRUE) # Groundwater budget (10/31/1961-9/30/2001)
LB = read.csv(file = "LB.csv", stringsAsFactors = FALSE, header = TRUE) # Land Budget budget (10/31/1961-9/30/2001)
RZ = read.csv(file = "RZ.csv", stringsAsFactors = FALSE, header = TRUE) # Root Zone budget (10/31/1961-9/30/2001)

# ======================================== #
# Load Packages                            #
# ======================================== #
library(extrafont)
library(scales)
library(ggplot2)
library(grid)
library(gridExtra)

# ======================================== #
# Unit Conversion Functions                #
# ======================================== #
af_to_L = function(a){
  a * 43560 * 28.3168
}
mg_to_Mton = function(a){
  a / (1000*907185*1000000)
}

# ======================================== #
# Boundary Conditions                      #
# ======================================== #

# Rock water interactions
Kang = 3000/mean(c(1719,2205)) # (Kang, 2016) 300 mg/L at 1719 ft. in Fresno co., 2205 ft. in Kern co.
Williamson = 2000/2500

# velocity field from C2VSim
dv12 = 0.76106 # Darcy velocity between centroids of layers 1-2 (ft/yr)
dv23 = -0.01772 # Darcy velocity between centroids of layers 2-3 (ft/yr)

# C2VSim Groundwater budget (units = AF)
DP = sum(GW$DP) # deep percolation
Si = sum(GW$Si) # beginning storage
Sf = sum(GW$Sf) # ending storage
NDP = sum(GW$NDP) # net deep percolation
Stream = sum(GW$Stream) # stream gain
R = sum(GW$R) # recharge
Lake = sum(GW$Lake) # lake gain
BI = sum(GW$BI) # sum(GW$Si) boundary inflow
Sub = sum(GW$Sub) # subsidence flow
P = sum(GW$P) # GW pumping
SI = sum(GW$Net.Subsurface.Inflow....) # subsurface inflow

# spatial/temporal dimensions
Dim = data.frame(subregion = seq(14,21,1), 
                 SA = c(670228.87, 904472.46, 302449.08, 372889.22, 897090.81, 801419.62, 423713.29, 652847.11), # acre feet
                 thickness = c(1915.155, 1990.01, 1509.194, 1868.821, 1904.244, 1715.117, 1969.524, 3244.262)) # feet: from C2VSim model

Tot_SA = sum(Dim$SA) # total surface area (acres)
V = Dim$SA * Dim$thickness # Volume (AF)
Tot_V = af_to_L(sum(V)) # total volume (L)
SA_w_thickness = sum(Dim$SA/Tot_SA*Dim$thickness) # surface area weighted thickness

# time step length
t = 50 # years

# length of C2VSim simulation
l = 40 # years

# C2VSim average annual fluxes from GW budget # years
annual_fluxes = data.frame(term = c("NDP", "Stream", "R", "Lake", "BI", "Sub", "P", "SI"),
                           "L/yr" = c(af_to_L(NDP/l), af_to_L(Stream/l), af_to_L(R/l), af_to_L(Lake/l), af_to_L(BI/l), af_to_L(Sub/l), af_to_L(P/l),af_to_L(SI/l)))
annual_fluxes$TDS = c(0,200,200,200,200,200,0,400) # TDS
annual_fluxes$mass_flux = annual_fluxes$L.yr * annual_fluxes$TDS # mass flux (mg/yr)

# Internal mass flux
internal_mf = (annual_fluxes[5,4] + annual_fluxes[6,4] + annual_fluxes[8,4])*t # mg/50 yr 
# Top down Mass Flux
top_mf = (annual_fluxes[2,4] + annual_fluxes[3,4] + annual_fluxes[4,4])*t # mg/50 yr 

# Root Zone Budget (includes Agriculture, Urban, and Native Vegetation)
I = sum(RZ$AG_INF) + sum(RZ$URB_INF) + sum(RZ$NAT_INF) # Net infiltration
ET = sum(RZ$AG_ET) + sum(RZ$URB_ET) + sum(RZ$NAT_ET) # Net ET

pI = NDP/I # percentage of infiltration that becomes NDP

pGWP = P/I # Percentage of groundwater pumping in Total Applied Water
pSW = 1 - pGWP # percentage of surface water diversions and reuse in Total Applied Water

TDS_SW = 200 # TDS of surface water
V_SW = (I - P) / l # Average Annual SW diversion volume (AF/yr)
SW_annual_mass = mg_to_Mton(af_to_L(V_SW) * TDS_SW) # Mton

# ======================================== #
# Number of Monte Carlo Samples            #
# ======================================== #
N = 1000

# ======================================== #
# Initialize 3D output arrays              #
# ======================================== #
TDS = array(0, dim=c(8,8,N)) # TDS array
layer_depths_array = array(0, dim=c(8,1,N)) # layer depths array
layer_velocity_array = array(0, dim=c(8,1,N)) # layer velocity array

mass = array(0, dim=c(8,1,N)) # array of sum of masses per time step
TDS_GWP = array(0, dim=c(8,1,N)) # array of TDS of pumped groundwater per time step
TDS_AW = array(0, dim=c(8,1,N)) # array of TDS of applied water (SW+GW) per time step
EC = array(0, dim=c(8,1,N)) # array of evaporative concentration of NDP per time step (TDS)
NDP_mf = array(0, dim=c(8,1,N)) # array of NDP mass flux per time step (mg/yr)
RWI_cont = array(0, dim=c(8,1,N)) # array of Rock water interactions
lay_vol = array(0, dim=c(8,1,N)) # array of layer volumes

# ======================================== #
# Compute all output                       #
# ======================================== #
for(k in 1:N){
  
  # Random Variables
  ## Aquifer parameters
  n = runif(1, min = 0.25, max = 0.35) # porosity
  pa = runif(1, min = 0.35, max = 0.50) # percent aquifer
  aVol = Tot_V * n * pa # saturated aquifer volume (L)
  
  ## Rock water interactions (TDS/ft) 
  RWI = runif(1, min = Williamson, max = Kang) 
  
  # velocity field from C2VSim
  v12 = 0.76106/n # linear pore velocity between centroids of layers 1-2 (ft/yr)
  v23 = -0.01772/n # linear pore velocity between centroids of layers 2-3 (ft/yr)
  vel = c(v12,v23)
  depths = c(-631.5,-2077) # C2VSim midpoints between layer centroids over which velcoity is calculated
  vel = as.data.frame(cbind(vel,depths))
  lm = lm(vel$depths~vel$vel)
  
  # vertical velocity profile
  v = matrix(0,8,1)
  b = matrix(0,8,1)
  v[1] = (0 + -coef(lm)[[1]]) / coef(lm)[[2]]
  b[1] = -(v[1] * t)
  
  for(i in 1:7){
    v[i+1] = (sum(b) -coef(lm)[[1]]) / coef(lm)[[2]]
    b[i+1] = -(v[i+1] * t)
  }
  
  # Vertical Flux per layer (L/yr)
  Q = v * n * pa * Tot_SA * 43560 * 28.3168 # L/yr
  layer_fluxes = cbind(v, Q)
  #names(layer_fluxes) = c("ft/yr","L/yr")
  
  # layer thickness (ft)
  b[8] = -(SA_w_thickness + sum(b[1:7]))
  pb = abs(b)/SA_w_thickness # percentage of total thickness per layer
  
  # layer volumes (L)
  lay_vol[,1,k] = aVol * pb # saturated volume of all layers (L)
  
  # ======================================== #
  # Calculate Pumping in each layer          #
  # ======================================== #
  lay_pump = matrix(0,8,1)
  for(i in 1:7){
    lay_pump[i,1] = b[i]/sum(b) * annual_fluxes[7,2] # proportion pumping from each layer (L/yr)
  }
  
  # ======================================== #
  # Calculate the layer depths               #
  # ======================================== #
  layer_depths = matrix(0,8)
  layer_depths[1] = abs(b[1])
  for(i in 2:8){
    layer_depths[i,1] = abs(layer_depths[i-1,1]) + abs(b[i]) # vector of layer depths (ft)
  }
  
  # ======================================== #
  # Initalize TDS array with baseline TDS    #
  # ======================================== #
  t0 = layer_depths * RWI # inital TDS-depth profile from RWI
  TDS[,1,k] = t0
  # RWI
  RWI_cont[,1,k] = t0
  
  # ======================================== #
  # Record layer depth/velocity for all sim  #
  # ======================================== #
  layer_depths_array[,,k] = layer_depths
  layer_velocity_array[,,k] = v
  
  # ======================================== #
  # Compute TDS for all layers and times     #
  # ======================================== # 
  for(j in 1:7){
    # calculate mass flux of NDP
    mass[j,1,k] = sum(TDS[1:7,j,k] * lay_vol[1:7,1,k]) # sum of mass in pumping layers 1-7 (mg)
    TDS_GWP[j,1,k] = mass[j,1,k]/sum(lay_vol[1:7,1,k]) # TDS of pumped GW
    TDS_AW[j,1,k] = TDS_GWP[j,1,k] * pGWP + TDS_SW * pSW # TDS of applied water for irrigation
    EC[j,1,k] = TDS_AW[j,1,k] * (1/pI) # TDS of NDP Evaporative concentration
    NDP_mf[j,1,k] = EC[j,1,k] * annual_fluxes[1,2] 
    
    # TDS in layer 1
    TDS[1,j+1,k] = ((TDS[1,j,k]*lay_vol[1,1,k] + top_mf + NDP_mf[j,1,k]*t - TDS[1,j,k]*layer_fluxes[2,2]*t - 
                       TDS[1,j,k]*lay_pump[1,1]*t + internal_mf*pb[1]) / lay_vol[1,1,k]) + RWI_cont[1,1,k]
    # TDS in layers 2-7
    for(i in 2:7){
      TDS[i,j+1,k] = ((TDS[i,j,k]*lay_vol[i,1,k] + internal_mf*pb[i] + TDS[i-1,j,k]*layer_fluxes[i,2]*t -
                         TDS[i,j,k]*lay_pump[i,1]*t - TDS[i,j,k]*layer_fluxes[i+1,2]*t) / lay_vol[i,1,k]) + 
        RWI_cont[i,1,k]
    }
  }
}

# ======================================== #
# calculate statistics for output          #
# ======================================== #

# for Plot 1:
TDS_mean = as.data.frame(matrix(0,8,8))
for(i in 1:8){
  for(j in 1:8){
    TDS_mean[i,j] = mean(TDS[i,j,])
  }
}

sd.top = matrix(0,7,8)
sd.bottom = matrix(0,7,8)
for(d in 1:7){
  for(t in 1:8){
    sd.top[d,t] = TDS_mean[d,t] - sd(TDS[d,t,]) # sd
    sd.bottom[d,t] = TDS_mean[d,t] + sd(TDS[d,t,]) # sd
  }
}

TDS_5 = matrix(0,7,8)
TDS_95 = matrix(0,7,8)
for(d in 1:7){
  for(t in 1:8){
    TDS_5[d,t] = quantile(TDS[d,t,], c(0.05)) # 5% CI
    TDS_95[d,t] = quantile(TDS[d,t,], c(.95))  # 95% CI
  }
}

## calculate layer depths/velocities of simulation means for main plot
layer_depths_mean = as.data.frame(matrix(0,7,1))
layer_velocity_mean = as.data.frame(matrix(0,7,1))
for(j in 1:7){
  layer_depths_mean[j,1] = mean(layer_depths_array[j,1,]) * -1
  layer_velocity_mean[j,1] = mean(layer_velocity_array[j,1,])
}

## reshape layer_depths_mean for plotting
layer_depths_mean = t.default(layer_depths_mean)

## trim 8th layer of TDS_mean matrix for plotting
TDS_mean = TDS_mean[1:7,] 


# for Plot 2:
GWP_mean = as.data.frame(matrix(0,7,1))
AW_mean = as.data.frame(matrix(0,7,1))
EC_mean = as.data.frame(matrix(0,7,1))
GWP_5 = as.data.frame(matrix(0,7,1))
AW_5 = as.data.frame(matrix(0,7,1))
EC_5 = as.data.frame(matrix(0,7,1))
GWP_95 = as.data.frame(matrix(0,7,1))
AW_95 = as.data.frame(matrix(0,7,1))
EC_95 = as.data.frame(matrix(0,7,1))
for(i in 1:7){
  GWP_mean[i,1] = mean(TDS_GWP[i,1,])
  AW_mean[i,1] = mean(TDS_AW[i,1,])
  EC_mean[i,1] = mean(EC[i,1,])
  GWP_5[i,1] = quantile(TDS_GWP[i,1,], c(0.05))
  GWP_95[i,1] = quantile(TDS_GWP[i,1,], c(0.95))
  AW_5[i,1] = quantile(TDS_AW[i,1,], c(0.05))
  AW_95[i,1] = quantile(TDS_AW[i,1,], c(0.95))
  EC_5[i,1] = quantile(EC[i,1,], c(0.05))
  EC_95[i,1] = quantile(EC[i,1,], c(0.95))
}

## compile stats for TDS of pumped GW, Total Applied Water (GW+SW), and Evaporative concentration of NDP
EC_dat = data.frame("legend"=rep(c("A","B","C"), each=7),
                    "time" = factor(seq(50,350,50)),
                    "mean" = rbind(GWP_mean, AW_mean, EC_mean),
                    "5 CI" = rbind(GWP_5, AW_5, EC_5),
                    "95 CI" = rbind(GWP_95, AW_95, EC_95))
colnames(EC_dat) = c("legend","time","mean","5 CI", "95 CI")


# for Table 1:
RWI_mean = mean(colSums((RWI_cont[1:7,1,]*lay_vol[1:7,1,])/50)) # mg
RWI_5 = quantile(colSums(RWI_cont[1:7,1,]), c(0.05))
RWI_95 = quantile(colSums(RWI_cont[1:7,1,]), c(0.95))


# ======================================== #
# Plot Output                              #
# ======================================== #

# plot 1: salinization evolution: 0-350 yrs

## data frame of all variables to plot
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

ymax = 25000 # maximum TDS concentration (mg/L) to display on y axis

## 0 yrs
p1 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS1), color = 'red') +
  geom_ribbon(aes(ymin = pl$L1, ymax = pl$U1), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("0 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 50 yrs
p2 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS2), color = 'red') +
  geom_ribbon(aes(ymin = pl$L2, ymax = pl$U2), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("50 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 100 yrs
p3 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS3), color = 'red') +
  geom_ribbon(aes(ymin = pl$L3, ymax = pl$U3), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("100 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 150 yrs
p4 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS4), color = 'red') +
  geom_ribbon(aes(ymin = pl$L4, ymax = pl$U4), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("150 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 200 yrs
p5 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS5), color = 'red') +
  geom_ribbon(aes(ymin = pl$L5, ymax = pl$U5), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("200 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 250 yrs
p6 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS6), color = 'red') +
  geom_ribbon(aes(ymin = pl$L6, ymax = pl$U6), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("250 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 300 yrs
p7 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS7), color = 'red') +
  geom_ribbon(aes(ymin = pl$L7, ymax = pl$U7), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("300 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))

## 350 yrs
p8 = ggplot(pl, aes(depth)) +
  geom_line(aes(y = pl$TDS8), color = 'red') +
  geom_ribbon(aes(ymin = pl$L8, ymax = pl$U8), alpha=.4) + 
  aes(ymax = ymax) +
  xlab("TDS (mg/L)") +
  ylab("depth below land surface (ft.)")  + 
  ggtitle("350 yrs.") + 
  geom_hline(yintercept = 1000, linetype = "longdash") +
  theme_light() + 
  theme(text=element_text(family="Times New Roman", size = 12)) + 
  coord_flip() + 
  scale_y_continuous(label=comma, breaks = c(5000,15000,25000))


## arrange p1 - p8 into one plot
plot1 = grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=3, nrow=3, vp=viewport(width=0.95, height=0.95))

## visualize
plot1


# plot 2: evapoconcentration: compare concentrations of GW, TAW, NDP
plot2 = ggplot(EC_dat, aes(x=time, y=mean, fill=legend)) +
  geom_bar(position=position_dodge(), stat="identity", 
           colour = "black", size=.3) +
  geom_errorbar(aes(ymin=EC_dat$`5 CI`, ymax=EC_dat$`95 CI`),
                size=.3,    # thinner lines
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  xlab("Time (yrs)") +
  ylab("TDS (mg/L)") +
  scale_fill_hue(name="",
                 breaks=c("A","B","C"),
                 labels=c("Pumped Groundwater", 
                          "Total Applied Water", 
                          "Net Deep Percolation")) +
  ggtitle("Concentration of pumped GW, TAW, and NDP") +
  theme_bw() + scale_y_continuous(labels = comma) + 
  theme(text=element_text(family="Times New Roman", size=18)) +
  theme(legend.position = "bottom", legend.text = element_text(size=16)) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18))

## visualize
plot2


# plot 3: average vertical velocity profile
plot3data = data.frame(x = layer_velocity_mean[1:7,1], 
                  y = as.matrix(layer_depths_mean[1:7,1]))
plot3 = ggplot(plot3data, aes(x = x, y = y)) + 
  geom_smooth(method="lm", se= FALSE, color = c("gray40"), formula = y ~ x) +
  geom_point(shape=16, cex=3) +
  xlab("averge linear pore velocity (ft/yr)") +
  ylab("depth below land surface") +
  theme_bw() +
  ggtitle("Mean Linear Pore Velocity") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18)) + 
  theme(text=element_text(family="Times New Roman", size=18)) 

## visualize
plot3


# table 1: mass loading into entire GW system from various water budget terms and sources (Mtons)
tstep = 1 # time step of NDP
table1 = data.frame("term" = c("Stream","R","Lake","BI","Sub","SW_D",
                               "SI","EC_mean", "EC_5","EC_95", 
                               "RWI_mean","RWI_5", "RWI_95"),
                    "mass" = c(mg_to_Mton(annual_fluxes[2,4]), mg_to_Mton(annual_fluxes[3,4]),
                               mg_to_Mton(annual_fluxes[4,4]), mg_to_Mton(annual_fluxes[5,4]),
                               mg_to_Mton(annual_fluxes[6,4]), SW_annual_mass,
                               mg_to_Mton(annual_fluxes[8,4]),
                               mg_to_Mton(EC_mean[tstep,1] * annual_fluxes[1,2]),
                               mg_to_Mton(EC_5[tstep,1] * annual_fluxes[1,2]), 
                               mg_to_Mton(EC_95[tstep,1] * annual_fluxes[1,2]),
                               mg_to_Mton(RWI_mean), 
                               mg_to_Mton(RWI_5),
                               mg_to_Mton(RWI_95)))

## visualize
table1

sum(table1[1:8,2]) # new mass minus RWI (Steam, R, Lake, BI, Sub, SW_D, SI, EC)
table1[11,2] # mass from RWI 
  

## calculate annual mass loading in entire system from RWI and non-RWI sources per time step
table1 = matrix(0,4,7)
for(tstep in 1:7){
  table1[1,tstep] = sum(mg_to_Mton(annual_fluxes[2,4]), mg_to_Mton(annual_fluxes[3,4]),
                        mg_to_Mton(annual_fluxes[4,4]), mg_to_Mton(annual_fluxes[5,4]),
                        mg_to_Mton(annual_fluxes[6,4]), SW_annual_mass,
                        mg_to_Mton(annual_fluxes[8,4]),
                        mg_to_Mton(EC_mean[tstep,1] * annual_fluxes[1,2])) 
  table1[2,tstep] = sum(mg_to_Mton(annual_fluxes[2,4]), mg_to_Mton(annual_fluxes[3,4]),
                        mg_to_Mton(annual_fluxes[4,4]), mg_to_Mton(annual_fluxes[5,4]),
                        mg_to_Mton(annual_fluxes[6,4]), SW_annual_mass,
                        mg_to_Mton(annual_fluxes[8,4]),
                        mg_to_Mton(EC_mean[tstep,1] * annual_fluxes[1,2])) - 
                            ((mg_to_Mton(EC_mean[tstep,1] * annual_fluxes[1,2])) - (mg_to_Mton(EC_5[tstep,1] * annual_fluxes[1,2])))
  table1[3,tstep] = sum(mg_to_Mton(annual_fluxes[2,4]), mg_to_Mton(annual_fluxes[3,4]),
                        mg_to_Mton(annual_fluxes[4,4]), mg_to_Mton(annual_fluxes[5,4]),
                        mg_to_Mton(annual_fluxes[6,4]), SW_annual_mass,
                        mg_to_Mton(annual_fluxes[8,4]),
                        mg_to_Mton(EC_mean[tstep,1] * annual_fluxes[1,2])) + 
                            ((mg_to_Mton(EC_95[tstep,1] * annual_fluxes[1,2])) - (mg_to_Mton(EC_mean[tstep,1] * annual_fluxes[1,2])))
  table1[4,tstep] = mg_to_Mton(RWI_mean)
}

table1 = as.data.frame(table1)
rownames(table1) = c('NON-RWI mean','NON-RWI 5 CI', 'NON-RWI 95 CI',
                     'RWI mean')
colnames(table1) = paste(seq(50,350,50), "yrs", " ")
table1 = t.data.frame(table1)

mass_dat = data.frame("legend"=rep(c("A","B"), each=7),
                    "time" = factor(seq(50,350,50)),
                    "mean" = c(table1[,1], table1[,4]),
                    "5 CI" = c(table1[,2], rep(0,7)),
                    "95 CI" = c(table1[,3], rep(0,7)))
colnames(mass_dat) = c("legend","time","mean","5 CI", "95 CI")


# plot 4: Annual mass loading in entire system from RWI and non-RWI sources per time step
plot4 = ggplot(mass_dat, aes(x=time, y=mean, fill=legend)) +
  geom_bar(position=position_dodge(), stat="identity", 
           colour = "black", size=.3) +
  geom_errorbar(aes(ymin=mass_dat$`5 CI`, ymax=mass_dat$`95 CI`),
                size=.3,    # thinner lines
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  xlab("Time (yrs)") +
  ylab("mass (Million Tons)") +
  scale_fill_hue(name="",
                 breaks=c("A","B"),
                 labels=c("Water Budget + EC of NDP", 
                          "Rock Water Interactions")) +
  ggtitle("Annual Mass Loading in Entire System") +
  theme_bw() + scale_y_continuous(labels = comma) + 
  theme(text=element_text(family="Times New Roman", size=18)) +
  theme(legend.position = "bottom", legend.text = element_text(size=16)) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18))

## visualize
plot4



