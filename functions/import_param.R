### Functions for importing the parameters of the PBPK model -------------------------------------------------------------------

# Author: Nicola Melillo
# Date: 04/03/2021
# Version: v0.1

# example on how to use these functions is at the end of the file

### output of getPBPKParam
# The output is a list param.PBPK containing 14 elements:
# - organ_bf     [L/h]    organs blood flow
# - organ_w      [Kg]     organs weight
# - organ_d      [Kg/L]   organs density
# - organ_v      [L]      organs volumes
# - ACAT_v_lum   [L]      volumes of the ACAT compartments
# - ACAT_l       [cm]     length of the various ACAT sections
# - ACAT_d       [cm]     diameter of the various ACAT sections
# - ACAT_pH      []       pH in a given ACAT section
# - ACAT_v_ent   [L]      enterocytes volumes of a given ACAT section
# - ACAT_f_CO    []       fraction of the cardiac output (CO) directed to the given enterocytes section
# - general_p             - "weight"    [Kg]  weight of the subjects
#                         - "Htc"       []    haematocrit
#                         - "water_po"  [L]   volume of water drinked together with the oral dose (standard 250 mL)
# - param_drug            drug related parameters specified by the user before calling the function
# - comp_PBPK_names       names of PBPK compartments
# - comp_ACAT_names       names of the ACAT compartments
# - part_coeff   []       partition coefficients
# - physical_param

### PBPK organs for which we have the various parameters (in order)
# 1  lungs
# 2  brain
# 3  heart
# 4  kidneys
# 5  bone
# 6  muscle
# 7  stomach
# 8  spleen
# 9  liver
# 10 gut
# 11 pancreas
# 12 skin
# 13 fat
# 14 arterial_blood
# 15 venous_blood

### ACAT compartments
# 1 stomach
# 2 duodenum
# 3 jejunum1
# 4 jejunum2
# 5 ileum1
# 6 ileum2
# 7 ileum3
# 8 ileum4


### general notes
# - as now only a mean human adult female and male subject of 70kg is supported
# - as now only Poulin & Theil (PT) and Berezhkovsky (bere) partition coefficients are supported

### todos
# - RR partition coefficients

### function get partition coefficients ----------------------------------------------------------------------------------------
# https://github.com/metrumresearchgroup/PBPK_PC/tree/master
pc_sel <- function(logP, pKa, fup, BP, type, method){
  # PC calculation (import dataset)
  dat_uni <- read.csv("./data/unified_tissue_comp.csv")         # data for unified tissue composition
  
  if(method=="P&T"){
    uni <- calcKp_PT(logP, pKa, fup, BP, type, dat=dat_uni)        #prediction using unified tissue composition
  }
  else if(method=="Berezhkovskiy"){
    uni <- calcKp_Berez(logP, pKa, fup, BP, type, dat=dat_uni)     #prediction using unified tissue composition
  }
  else if(method=="R&R"){
    uni <- calcKp_RR(logP, pKa, fup, BP, type, dat=dat_uni)        #prediction using unified tissue composition
  }
  else if(method=="Schmitt"){
    uni <- calcKp_Schmitt(logP, pKa, fup, type, dat=dat_uni)       #prediction using unified tissue composition
  }
  else if(method=="PK-Sim"){
    uni <- calcKp_pksim(logP, fup, dat_uni)                    #prediction using unified tissue composition
  }
  drug_vals <- data.frame(kcp = unlist(uni))
  return(drug_vals)
}

getPartitionCoeff <- function(param_drug, type_part_coeff) {
  
  # derive drug parameters
  logPow <- log10(param_drug["Pow"])
  pKa    <- param_drug["pKa"]
  fup    <- param_drug["fup"]
  fut    <- param_drug["fut"] # *not used yet
  BP     <- param_drug["BP"]
  type   <- param_drug["type"]
  
  # Partition coefficient method as (n=5):
  # > Poulin & Theil (P&T), Berezhkovskiy, Rodgers & Rowland (R&R), Schmitt, PK-Sim standard
  # basic physiologic properties
  f_water_int <- 0.935 # fractional volume content of water in interstitial (human)
  f_water_pls <- 0.926 # fractional volume content of water in plasma (human)
  f_proteins_int_to_pls <- 0.37
  # test line:
  PC <- pc_sel(logP=logPow, pKa=pKa, fup=fup, BP=BP, type=type, method=type_part_coeff) # kcp: partition coefficient (intracellular/plasma)
  
  # PC Calculation in interfaces
  PC$kip <- (f_water_int + f_proteins_int_to_pls * (1/fup - f_water_pls)) * fup # kip: partition coefficient (interstitial/plasma)
  PC$kcw <- fup/PC$kcp # partition coefficient (intracellular/water)
  PC$kiw <- fup/PC$kip # partition coefficient (interstitial/water)
  
  return(PC)
}


### function get PBPK model parameters -----------------------------------------------------------------------------------------
getPBPKParam <- function(filename, param_drug, type_part_coeff, specie){
  
  ### notes:
  # - filename must be formatted as the file I origially used
  # - to read the excel files I used read_excel from readxl package. This return tibbles. I want as output named vectors, thus, I need to convert them.
  # - all the vectors relative to organs magnitude MUST have the same organs order
  
  library(readxl)
  
  ### read the files (sheet names hard coded - modify for rats)
  sheet_param_organs          <- "parameters_organs"
  sheet_param_general         <- "parameters_general"
  sheet_param_acat            <- "parameters_ACAT"
  sheet_param_part_coeff_RR   <- "parameters_partition_coeff_RR"
  sheet_param_part_coeff_PT   <- "parameters_partition_coeff_PT"
  sheet_param_dissolution     <- "physical_param_dissolution"
  
  param_organs        <- read_excel(filename, sheet=sheet_param_organs)
  param_general       <- read_excel(filename, sheet=sheet_param_general)
  param_ACAT          <- read_excel(filename, sheet=sheet_param_acat)
  param_part_coeff_RR <- read_excel(filename, sheet=sheet_param_part_coeff_RR)
  param_part_coeff_PT <- read_excel(filename, sheet=sheet_param_part_coeff_PT)
  param_diss          <- read_excel(filename, sheet=sheet_param_dissolution)
  
  ### get partition coefficients
  PC <- getPartitionCoeff(param_drug, type_part_coeff)
  # Define PCs
  PC_kcp <- structure(PC$kcp, names=rownames(PC)) # intracellular - plasma
  PC_kip <- structure(PC$kip, names=rownames(PC)) # interstitial - plasma
  PC_kcw <- structure(PC$kcw, names=rownames(PC)) # intracellular - water
  PC_kiw <- structure(PC$kiw, names=rownames(PC)) # interstitial - water
  
  ### readapt ORGANS parameters
  # 1) base organ param
  organ_bf <- structure(param_organs$blood_flow, names=param_organs$organ_name) # [L/h]
  organ_d  <- structure(param_organs$density, names=param_organs$organ_name)    # [Kg/L]
  organ_v  <- structure(param_organs$volume, names=param_organs$organ_name)     # whole volume [Kg]
  organ_w  <- organ_v*organ_d                                                   # [L]
  # 2) added organ param (detailed organ volume)
  organ_v_int <- organ_v*param_organs$f_int_vol                                 # interstitial volume [Kg]
  organ_v_vas <- organ_v*param_organs$f_vasc_vol                                # vascular volume [Kg]
  organ_v_cell <- organ_v - organ_v_int - organ_v_vas                           # cellular volume [Kg]
  
  ### surface area (SA)
  # 1) plasma <-> interstitial SA
  SA_factor <- 950 # [1/cm]
  SA_pls_int <- 1000*SA_factor*organ_v_vas # factor [1/cm] * interstitial volume [Kg to cm3, 1mg/mL density assumption]
  # > organ vascularization method by Niederalt et al. (doi: 10.1007/s10928-017-9559-4)
  
  # 2) interstitial <-> intracellular SA (PKsim based, fitted value)
  SA_int_cell <- sapply(names(organ_v), FUN = function(x){
    value <- organ_v[x] # total organ volume
    if(x=="lungs"){         return((value * 1000  / 2.2)^0.75 * 0.096 * 1E2 * 100)} # lungs
    if(x=="brain"){         return((value * 1000 / 1.671)^0.75 * 0.0006 * 1E2 * 100)} # brain
    if(x=="heart"){         return((value * 1000 / 1.2)^0.75 * 7.54 * 1E2 * 100)} # heart
    if(x=="kidneys"){       return((value * 1000 / 7)^0.75 * 1000 * 1E2 * 100)} # kidneys
    if(x=="bone"){          return((value * 1000 / 28.2)^0.75 * 10 * 1E2 * 100)} # bone
    if(x=="muscle"){        return((value * 1000 / 110.1)^0.75 * 7.54 * 1E2 * 100)} # muscle
    if(x=="stomach"){       return((value * 1000 / 1.1)^0.75 * 1000 * 1E2 * 100)} # stomach
    if(x=="spleen"){        return((value * 1000 / 1.3)^0.75 * 1000 * 1E2 * 100)} # spleen
    if(x=="liver"){         return((value * 1000 / 10)^0.75 * 82 * 1E2 * 100)} # liver
    if(x=="gut"){           return((value * 1000 / 11.1)^0.75 * 1000 * 1E2 * 100)} # gut = large + small
    if(x=="pancreas"){      return((value * 1000 / 1.3)^0.75 * 1000 * 1E2 * 100)} # pancreas
    if(x=="skin"){          return((value * 1000 / 43.4)^0.75 * 0.12 * 1E2 * 100)} # skin
    if(x=="fat"){           return((value * 1000 / 14.2)^0.75 * 5 * 1E2 * 100)} # fat
    if(x=="arterial_blood"){return(0)} # arterial_blood
    if(x=="venous_blood"){  return(0)} # venous_blood
  }) |> unlist(use.names = FALSE) |> structure(names=names(organ_v))
  
  ### readapt ACAT parameters
  ACAT_v_lum <- structure(param_ACAT$volume_lum, names=param_ACAT$compartment)  # [L]
  ACAT_l     <- structure(param_ACAT$length, names=param_ACAT$compartment)      # [cm]
  ACAT_d     <- structure(param_ACAT$diameter, names=param_ACAT$compartment)    # [cm]
  ACAT_pH    <- structure(param_ACAT$pH, names=param_ACAT$compartment)
  ACAT_v_ent <- structure(param_ACAT$volume_ent, names=param_ACAT$compartment)  # [L]
  ACAT_f_CO  <- structure(param_ACAT$fraction_CO, names=param_ACAT$compartment)
  
  ### readapt and define some dissolution parameters
  physical_param <- structure(param_diss$value, names=param_diss$parameter)
  
  ### readapt GENERAL parameters
  # weight [Kg]
  # Htc    
  # water_po
  general_p <- structure(param_general$value, names=param_general$parameter)
  
  ### get PBPK compartments names (add sink compartments)
  comp_PBPK_names <- names(organ_v)
  comp_PBPK_names <- c(comp_PBPK_names, "sink_CLh", "sink_CLr")
  
  ### define some parameter for the dissolution
  # for details of all these equations see supplementary materials of https://doi.org/10.1007/s10928-018-9615-8 
  kb    <- physical_param["kb"]            # [J/K]              Boltzmann constant
  an    <- physical_param["avogadro_num"]  # [mol^-1]           Avogadro number
  eta_w <- physical_param["eta_w"]         # [Pa*s]=[s*J/m^3]   dynamic viscosity of water at 37 ?C
  temp  <- physical_param["temp"]          # [K]                absolute temperature 
  
  # derive diffusivity parameter
  # diffusion coefficient calculated with the Stockes Einstein equation
  # hypothesis is that the drug molecule is spherical in shape
  Rh   <- ( (3*param_drug["mw"])/(4*pi*an*param_drug["rho"]) )^(1/3)*10^-1  # [m]                hydrodynamic radius of diffusing drug
  Diff <- (kb*temp/(6*pi*eta_w*Rh))*10^4*60*60;                             # [m^2/s]->[cm^2/h]  diffusion coefficient
  
  # calculate the effective thickness of the hydrodynamic diffusion layer
  # Hintz and Johnson model
  if(param_drug["r"]>30){
    ht <- 30               # [um]
  }else{
    ht <- param_drug["r"]  # [um]
  }
  
  # add Diff and h to the drug related parameters
  param_drug["ht"]   <- ht
  param_drug["Diff"] <- Diff
  
  ### define ACAT compartment names
  comp_ACAT_names_or <- names(ACAT_v_lum)
  comp_ACAT_names_s  <- comp_ACAT_names_or                                  # solid names
  comp_ACAT_names_d  <- comp_ACAT_names_or                                  # dissolved names
  comp_ACAT_names_e  <- comp_ACAT_names_or[2:(length(comp_ACAT_names_or))]  # enterocytes names
  comp_ACAT_names_ea <- comp_ACAT_names_e                                   # drug absorbed
  comp_ACAT_names_ec <- comp_ACAT_names_e                                   # drug cleared by the hepatocytes
  
  for(i in 1:length(comp_ACAT_names_or)){
    comp_ACAT_names_s[i] <- paste(comp_ACAT_names_s[i], "s", sep="_")
    comp_ACAT_names_d[i] <- paste(comp_ACAT_names_d[i], "d", sep="_")
  }
  for(i in 2:length(comp_ACAT_names_or)){
    comp_ACAT_names_e[i-1]  <- paste(comp_ACAT_names_e[i-1], "e", sep="_")
    comp_ACAT_names_ea[i-1] <- paste(comp_ACAT_names_ea[i-1], "ea", sep="_")
    comp_ACAT_names_ec[i-1] <- paste(comp_ACAT_names_ec[i-1], "ec", sep="_")
  }
  comp_ACAT_names <- c(comp_ACAT_names_s, comp_ACAT_names_d, comp_ACAT_names_e, comp_ACAT_names_ea, comp_ACAT_names_ec, "sink_s", "sink_d")
  
  ### define output list
  output_list <- list(organ_bf = organ_bf,
                      organ_w = organ_w,
                      organ_d = organ_d,
                      organ_v = organ_v, # whole volume
                      organ_v_int = organ_v_int, # interstitial volume
                      organ_v_vas = organ_v_vas, # vascular volume
                      organ_v_cell = organ_v_cell, # cellular volume
                      SA_pls_int = SA_pls_int, # plasma-interstitial surface area
                      SA_int_cell = SA_int_cell, # interstitial-intracellular surface area
                      ACAT_v_lum = ACAT_v_lum,
                      ACAT_l = ACAT_l,
                      ACAT_d = ACAT_d,
                      ACAT_pH = ACAT_pH,
                      ACAT_v_ent = ACAT_v_ent,
                      ACAT_f_CO = ACAT_f_CO,
                      general_p = general_p,
                      param_drug = param_drug,
                      PC_kcp = PC_kcp, # partition coefficient (cell-plasma)
                      PC_kip = PC_kip, # partition coefficient (interstitial-plasma)
                      PC_kcw = PC_kcw, # partition coefficient (cell-water)
                      PC_kiw = PC_kiw, # partition coefficient (interstitial-water)
                      comp_PBPK_names = comp_PBPK_names,
                      comp_ACAT_names = comp_ACAT_names,
                      physical_param = physical_param
                      )
  
  return(output_list)
  
}






### reorganize parameters for RxODE --------------------------------------------------------------------------------------------

reorganizeParam.rxode <- function(param.PBPK){

  # get names compartments and remove from param.PBPK list
  comp_PBPK_names <- param.PBPK$comp_PBPK_names
  comp_ACAT_names <- param.PBPK$comp_ACAT_names
  param.PBPK$comp_PBPK_names <- NULL
  param.PBPK$comp_ACAT_names <- NULL
  
  names.list <- names(param.PBPK)
  param.PBPK.v <- c()
  
  for(i in 1:length(names.list)){
    vect.i <- unname(param.PBPK[[names.list[i]]])
    names(vect.i) <- paste(names.list[i],names(param.PBPK[[names.list[i]]]),sep="_")
    param.PBPK.v <- c(param.PBPK.v, vect.i)
  }
  
  return(param.PBPK.v)
  
  }




