### functions describing the PBPK model

### Units
# x        - [mg]
# volumes  - [L]
# time     - [h]

### PBPK organs for which we have the various parameters (in order)

# PBPK distribution
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

# PBPK sink
# 16 sink liver
# 17 sink kidney

# ACAT solid
# 18 stomach_s
# 19 duodenum_s
# 20 jejunum1_s
# 21 jejunum2_s
# 22 ileum1_s
# 23 ileum2_s
# 24 ileum3_s

# ACAT dissolved
# 25 stomach_d
# 26 duodenum_d
# 27 jejunum1_d
# 28 jejunum2_d
# 29 ileum1_d
# 30 ileum2_d
# 31 ileum3_d

# ACAT enterocytes
# 32 duodenum_e
# 33 jejunum1_e
# 34 jejunum2_e
# 35 ileum1_e
# 36 ileum2_e
# 37 ileum3_e

# ACAT_absorbed
# 38 duodenum_ea
# 39 jejunum1_ea
# 40 jejunum2_ea
# 41 ileum1_ea
# 42 ileum2_ea
# 43 ileum3_ea

# ACAT cleared
# 44 duodenum_ec
# 45 jejunum1_ec
# 46 jejunum2_ec
# 47 ileum1_ec
# 48 ileum2_ec
# 49 ileum3_ec

# ACAT sink
# 50 sink_s
# 51 sink_d


PBPK.ACAT <- RxODE({
  
  ### ACAT model ===============================================================
  
  input_GI <- 0
  
  # term in front of the dissolution equation
  # > 앞서 만들어놓은 diffusion rate(by Noyes식)을 이용한 공식/ Thickness by Hintz model(Kd,st)
  Kd1 <- 10^5*(3*param_drug_Diff/(rho*param_drug_ht*r)) # [L/(h*mg)] 
  
  # define Henderson-Hasselback term for dissolution
  # particle movement: stomach(s) => duodenum(d) => jejunum(j, 1-2) => ileum(i, 1-3) ~ total 7 section
  if(type==0){           # neutral
    Ys  <- 1
    Yd  <- 1
    Yj1 <- 1
    Yj2 <- 1
    Yi1 <- 1
    Yi2 <- 1
    Yi3 <- 1
  }else if(type==1){     # acid
    Ys  <- 1 + 10^( ACAT_pH_stomach  - pKa )
    Yd  <- 1 + 10^( ACAT_pH_duodenum - pKa )
    Yj1 <- 1 + 10^( ACAT_pH_jejunum1 - pKa )
    Yj2 <- 1 + 10^( ACAT_pH_jejunum2 - pKa )
    Yi1 <- 1 + 10^( ACAT_pH_ileum1   - pKa )
    Yi2 <- 1 + 10^( ACAT_pH_ileum2   - pKa )
    Yi3 <- 1 + 10^( ACAT_pH_ileum3   - pKa )
  } else if (type==2){   # base
    Ys  <- 1 + 10^( -ACAT_pH_stomach  + pKa )
    Yd  <- 1 + 10^( -ACAT_pH_duodenum + pKa )
    Yj1 <- 1 + 10^( -ACAT_pH_jejunum1 + pKa )
    Yj2 <- 1 + 10^( -ACAT_pH_jejunum2 + pKa )
    Yi1 <- 1 + 10^( -ACAT_pH_ileum1   + pKa )
    Yi2 <- 1 + 10^( -ACAT_pH_ileum2   + pKa )
    Yi3 <- 1 + 10^( -ACAT_pH_ileum3   + pKa ) 
  }
  
  ### stomach equations
  
  # define parameters and useful indices
  Cs_st <- Csint * Ys # intrinsic  solubility * correction (by gut section)
  Kd_st <- Kd1 * (Cs_st - stomach_d/(ACAT_v_lum_stomach + general_p_water_po)) # drug dissolution rate in stomach
  kst   <- 1/general_p_GET # [1/h], Inversed Gastric Empyting Time
  
  # define differential equations
  output_s_s   <- kst * stomach_s
  output_d_s   <- kst * stomach_d
  
  # Drug movement in stomach (first term:solid movement, second term:dissolved movement)
  d/dt(stomach_s) <- -output_s_s - Kd_st * stomach_s
  d/dt(stomach_d) <- -output_d_s + Kd_st * stomach_s
  
  
  ###### equations for all the small intestine transit compartments
  
  length_small_int <- ACAT_l_duodenum + ACAT_l_jejunum1 + ACAT_l_jejunum2 + ACAT_l_ileum1 + ACAT_l_ileum2 + ACAT_l_ileum3
  CO <- organ_bf_venous_blood  # cardiac output [L/h]
  
  ### duodenum 
  Cs_it <- Csint * Yd
  Kd_it <- Kd1 * ( Cs_it - duodenum_d/ACAT_v_lum_duodenum)
  kit   <- 1/(general_p_SITT * ACAT_l_duodenum/length_small_int) # [1/h]
  ka    <- 2 * Peff / (ACAT_d_duodenum/2) # [1/h]
  Qe_i  <- ACAT_f_CO_duodenum * CO # [L/h]
  
  # define differential equations
  d/dt(duodenum_s)  <- +output_s_s - Kd_it * duodenum_s - kit * duodenum_s
  d/dt(duodenum_d)  <- +output_d_s + Kd_it * duodenum_s - kit * duodenum_d - ka * duodenum_d
  d/dt(duodenum_e)  <- ka * duodenum_d - Qe_i * duodenum_e/ACAT_v_ent_duodenum - CLent * duodenum_e/ACAT_v_ent_duodenum
  d/dt(duodenum_ea) <- Qe_i  * duodenum_e/ACAT_v_ent_duodenum
  d/dt(duodenum_ec) <- CLent * duodenum_e/ACAT_v_ent_duodenum
  
  output_s_d   <- kit * duodenum_s
  output_d_d   <- kit * duodenum_d
  input_GI     <- input_GI + Qe_i  * duodenum_e/ACAT_v_ent_duodenum
  
  
  ### jejunum 1 
  Cs_it_j1 <- Csint * Yj1
  Kd_it_j1 <- Kd1 * ( Cs_it_j1 - jejunum1_d/ACAT_v_lum_jejunum1)
  kit      <- 1/(general_p_SITT * ACAT_l_jejunum1/length_small_int) # [1/h]
  ka       <- 2 * Peff / (ACAT_d_jejunum1/2) # [1/h]
  Qe_i     <- ACAT_f_CO_jejunum1 * CO # [L/h]
  
  # define differential equations
  d/dt(jejunum1_s)  <- +output_s_d - Kd_it_j1 * jejunum1_s - kit * jejunum1_s
  d/dt(jejunum1_d)  <- +output_d_d + Kd_it_j1 * jejunum1_s - kit * jejunum1_d - ka * jejunum1_d
  d/dt(jejunum1_e)  <- ka * jejunum1_d - Qe_i * jejunum1_e/ACAT_v_ent_jejunum1 - CLent * jejunum1_e/ACAT_v_ent_jejunum1
  d/dt(jejunum1_ea) <- Qe_i  * jejunum1_e/ACAT_v_ent_jejunum1
  d/dt(jejunum1_ec) <- CLent * jejunum1_e/ACAT_v_ent_jejunum1
  
  output_s_j1   <- kit * jejunum1_s
  output_d_j1   <- kit * jejunum1_d
  input_GI      <- input_GI + Qe_i  * jejunum1_e/ACAT_v_ent_jejunum1
  
  ### jejunum 2 
  Cs_it_j2 <- Csint * Yj2
  Kd_it_j2 <- Kd1 * ( Cs_it_j2 - jejunum2_d/ACAT_v_lum_jejunum2)
  kit      <- 1/(general_p_SITT * ACAT_l_jejunum2/length_small_int) # [1/h]
  ka       <- 2 * Peff / (ACAT_d_jejunum2/2) # [1/h]
  Qe_i     <- ACAT_f_CO_jejunum2 * CO # [L/h]
  
  # define differential equations
  d/dt(jejunum2_s)  <- +output_s_j1 - Kd_it_j2 * jejunum2_s - kit * jejunum2_s
  d/dt(jejunum2_d)  <- +output_d_j1 + Kd_it_j2 * jejunum2_s - kit * jejunum2_d - ka * jejunum2_d
  d/dt(jejunum2_e)  <- ka * jejunum2_d - Qe_i * jejunum2_e/ACAT_v_ent_jejunum2 - CLent * jejunum2_e/ACAT_v_ent_jejunum2
  d/dt(jejunum2_ea) <- Qe_i  * jejunum2_e/ACAT_v_ent_jejunum2
  d/dt(jejunum2_ec) <- CLent * jejunum2_e/ACAT_v_ent_jejunum2
  
  output_s_j2   <- kit * jejunum2_s
  output_d_j2   <- kit * jejunum2_d
  input_GI      <- input_GI +  Qe_i * jejunum2_e/ACAT_v_ent_jejunum2
  
  ### ileum 1
  Cs_it_i1 <- Csint * Yi1
  Kd_it_i1 <- Kd1 * ( Cs_it_i1 - ileum1_d/ACAT_v_lum_ileum1)
  kit      <- 1/(general_p_SITT * ACAT_l_ileum1/length_small_int)      # [1/h]
  ka       <- 2 * Peff / (ACAT_d_ileum1/2)                             # [1/h]
  Qe_i     <- ACAT_f_CO_ileum1 * CO                                    # [L/h]
  
  # define differential equations
  d/dt(ileum1_s)  <- +output_s_j2 - Kd_it_i1 * ileum1_s - kit * ileum1_s
  d/dt(ileum1_d)  <- +output_d_j2 + Kd_it_i1 * ileum1_s - kit * ileum1_d - ka * ileum1_d
  d/dt(ileum1_e)  <- ka * ileum1_d - Qe_i * ileum1_e/ACAT_v_ent_ileum1 - CLent * ileum1_e/ACAT_v_ent_ileum1
  d/dt(ileum1_ea) <- Qe_i  * ileum1_e/ACAT_v_ent_ileum1
  d/dt(ileum1_ec) <- CLent * ileum1_e/ACAT_v_ent_ileum1
  
  output_s_i1   <- kit * ileum1_s
  output_d_i1   <- kit * ileum1_d
  input_GI      <- input_GI + Qe_i  * ileum1_e/ACAT_v_ent_ileum1
  
  
  ### ileum 2
  Cs_it_i2 <- Csint * Yi2
  Kd_it_i2 <- Kd1 * ( Cs_it_i2 - ileum2_d/ACAT_v_lum_ileum2)
  kit      <- 1/(general_p_SITT * ACAT_l_ileum2/length_small_int)      # [1/h]
  ka       <- 2 * Peff / (ACAT_d_ileum2/2)                             # [1/h]
  Qe_i     <- ACAT_f_CO_ileum2 * CO                                    # [L/h]
  
  # define differential equations
  d/dt(ileum2_s)  <- +output_s_i1 - Kd_it_i2 * ileum2_s - kit * ileum2_s
  d/dt(ileum2_d)  <- +output_d_i1 + Kd_it_i2 * ileum2_s - kit * ileum2_d - ka * ileum2_d
  d/dt(ileum2_e)  <- ka * ileum2_d - Qe_i * ileum2_e/ACAT_v_ent_ileum2 - CLent * ileum2_e/ACAT_v_ent_ileum2
  d/dt(ileum2_ea) <- Qe_i  * ileum2_e/ACAT_v_ent_ileum2
  d/dt(ileum2_ec) <- CLent * ileum2_e/ACAT_v_ent_ileum2
  
  output_s_i2   <- kit * ileum2_s
  output_d_i2   <- kit * ileum2_d
  input_GI      <- input_GI + Qe_i  * ileum2_e/ACAT_v_ent_ileum2
  
  ### ileum 3
  Cs_it_i3 <- Csint * Yi3
  Kd_it_i3 <- Kd1 * ( Cs_it_i3 - ileum3_d/ACAT_v_lum_ileum3)
  kit      <- 1/(general_p_SITT * ACAT_l_ileum3/length_small_int )     # [1/h]
  ka       <- 2 * Peff / (ACAT_d_ileum3/2)                             # [1/h]
  Qe_i     <- ACAT_f_CO_ileum3 * CO                                    # [L/h]
  
  # define differential equations
  d/dt(ileum3_s)  <- +output_s_i2 - Kd_it_i3 * ileum3_s - kit * ileum3_s
  d/dt(ileum3_d)  <- +output_d_i2 + Kd_it_i3 * ileum3_s - kit * ileum3_d - ka * ileum3_d
  d/dt(ileum3_e)  <- ka * ileum3_d - Qe_i * ileum3_e/ACAT_v_ent_ileum3 - CLent * ileum3_e/ACAT_v_ent_ileum3
  d/dt(ileum3_ea) <- Qe_i  * ileum3_e/ACAT_v_ent_ileum3
  d/dt(ileum3_ec) <- CLent * ileum3_e/ACAT_v_ent_ileum3
  
  output_s_i3   <- kit * ileum3_s
  output_d_i3   <- kit * ileum3_d
  input_GI      <- input_GI + Qe_i  * ileum3_e/ACAT_v_ent_ileum3
  
  ### equations for the sink compartments & model output
  
  # equations for the sink compartments
  d/dt(sink_s) <- output_s_i3
  d/dt(sink_d) <- output_d_i3
  
  
  
  ### PBPK model ===============================================================
  
  # preallocation of some support variable...
  bf_tot      <- 0
  outflow_tot <- 0
  outflow_spl <- 0
  
  ### parallel tissues equations...  c("brain","heart","kidneys","bone","muscle","skin","fat")
  
  # 1) brain -------------------------------------------------------------------
  # blood flows
  brain_inflow <- organ_bf_brain*(arterial_blood/organ_v_arterial_blood) # brain blood inflow
  brain_outflow <- organ_bf_brain*(brain_vas/organ_v_vas_brain) # brain blood outflow
  # concentrations
  brain_vas_conc <- brain_vas/organ_v_vas_brain # drug concentration in brain vascular
  brain_int_conc <- brain_int/organ_v_int_brain # drug concentration in brain interstitial
  brain_cell_conc <- brain_cell/organ_v_cell_brain # drug concentration in brain cells
  # diffusion gradient
  brain_intdiff <- fup*Pend*SA_pls_int_brain*(brain_vas_conc - brain_int_conc/(PC_kip_brain/BP)) # plasma-interstitial diffusion
  brain_celldiff <- Porg*SA_int_cell_brain*(brain_int_conc*PC_kiw_brain - brain_cell_conc*PC_kcw_brain) # interstitial-intracellular diffusion
  # diffential equations
  d/dt(brain_vas) <- brain_inflow - brain_outflow - brain_intdiff # vascular (blood flow + interstitial diffusion)
  d/dt(brain_int) <- brain_intdiff - brain_celldiff # interstitial (interstitial diffusion + intracellular diffusion)
  d/dt(brain_cell) <- brain_celldiff # cell (intracellular diffusion)
  
  bf_tot        <- bf_tot + organ_bf_brain # blood flow
  outflow_tot   <- outflow_tot + brain_outflow # total outflow
  
  # 2) heart -------------------------------------------------------------------
  heart_inflow <- organ_bf_heart * (arterial_blood/organ_v_arterial_blood) # blood inflow
  heart_outflow <- organ_bf_heart * (heart_vas/organ_v_vas_heart) # blood outflow
  # concentrations
  heart_vas_conc <- heart_vas/organ_v_vas_heart
  heart_int_conc <- heart_int/organ_v_int_heart
  heart_cell_conc <- heart_cell/organ_v_cell_heart
  # diffusion gradients
  heart_intdiff <- fup * Pend * SA_pls_int_heart * (heart_vas_conc - heart_int_conc/(PC_kip_heart/BP))
  heart_celldiff <- Porg * SA_int_cell_heart * (heart_int_conc * PC_kiw_heart - heart_cell_conc * PC_kcw_heart)
  # differential equations
  d/dt(heart_vas) <- heart_inflow - heart_outflow - heart_intdiff
  d/dt(heart_int) <- heart_intdiff - heart_celldiff
  d/dt(heart_cell) <- heart_celldiff
  
  bf_tot <- bf_tot + organ_bf_heart
  outflow_tot <- outflow_tot + heart_outflow
  
  # 3) kidneys -----------------------------------------------------------------
  clear_kid       <- CLr * fup * (kidneys_cell/organ_v_kidneys)/(PC_kcp_kidneys) + general_p_GFR * fup * arterial_blood/organ_v_arterial_blood/BP * GFR_flag
  kidneys_inflow <- organ_bf_kidneys * (arterial_blood/organ_v_arterial_blood)
  kidneys_outflow <- organ_bf_kidneys * (kidneys_vas/organ_v_vas_kidneys)
  # concentrations
  kidneys_vas_conc <- kidneys_vas/organ_v_vas_kidneys
  kidneys_int_conc <- kidneys_int/organ_v_int_kidneys
  kidneys_cell_conc <- kidneys_cell/organ_v_cell_kidneys
  # diffusion gradients
  kidneys_intdiff <- fup * Pend * SA_pls_int_kidneys * (kidneys_vas_conc - kidneys_int_conc/(PC_kip_kidneys/BP))
  kidneys_celldiff <- Porg * SA_int_cell_kidneys * (kidneys_int_conc * PC_kiw_kidneys - kidneys_cell_conc * PC_kcw_kidneys)
  # differential equations
  d/dt(kidneys_vas) <- kidneys_inflow - kidneys_outflow - kidneys_intdiff
  d/dt(kidneys_int) <- kidneys_intdiff - kidneys_celldiff
  d/dt(kidneys_cell) <- kidneys_celldiff - clear_kid
  
  bf_tot <- bf_tot + organ_bf_kidneys
  outflow_tot <- outflow_tot + kidneys_outflow
  
  # 4) bone --------------------------------------------------------------------
  bone_inflow <- organ_bf_bone * (arterial_blood/organ_v_arterial_blood)
  bone_outflow <- organ_bf_bone * (bone_vas/organ_v_vas_bone)
  # concentrations
  bone_vas_conc <- bone_vas/organ_v_vas_bone
  bone_int_conc <- bone_int/organ_v_int_bone
  bone_cell_conc <- bone_cell/organ_v_cell_bone
  # diffusion gradients
  bone_intdiff <- fup * Pend * SA_pls_int_bone * (bone_vas_conc - bone_int_conc/(PC_kip_bone/BP))
  bone_celldiff <- Porg * SA_int_cell_bone * (bone_int_conc * PC_kiw_bone - bone_cell_conc * PC_kcw_bone)
  # differential equations
  d/dt(bone_vas) <- bone_inflow - bone_outflow - bone_intdiff
  d/dt(bone_int) <- bone_intdiff - bone_celldiff
  d/dt(bone_cell) <- bone_celldiff
  
  bf_tot <- bf_tot + organ_bf_bone
  outflow_tot <- outflow_tot + bone_outflow
  
  # 5) muscle ------------------------------------------------------------------
  muscle_inflow <- organ_bf_muscle * (arterial_blood/organ_v_arterial_blood)
  muscle_outflow <- organ_bf_muscle * (muscle_vas/organ_v_vas_muscle)
  # concentrations
  muscle_vas_conc <- muscle_vas/organ_v_vas_muscle
  muscle_int_conc <- muscle_int/organ_v_int_muscle
  muscle_cell_conc <- muscle_cell/organ_v_cell_muscle
  # diffusion gradients
  muscle_intdiff <- fup * Pend * SA_pls_int_muscle * (muscle_vas_conc - muscle_int_conc/(PC_kip_muscle/BP))
  muscle_celldiff <- Porg * SA_int_cell_muscle * (muscle_int_conc * PC_kiw_muscle - muscle_cell_conc * PC_kcw_muscle)
  # differential equations
  d/dt(muscle_vas) <- muscle_inflow - muscle_outflow - muscle_intdiff
  d/dt(muscle_int) <- muscle_intdiff - muscle_celldiff
  d/dt(muscle_cell) <- muscle_celldiff
  
  bf_tot <- bf_tot + organ_bf_muscle
  outflow_tot <- outflow_tot + muscle_outflow
  
  
  # 6) skin --------------------------------------------------------------------
  skin_inflow <- organ_bf_skin * (arterial_blood/organ_v_arterial_blood)
  skin_outflow <- organ_bf_skin * (skin_vas/organ_v_vas_skin)
  # concentrations
  skin_vas_conc <- skin_vas/organ_v_vas_skin
  skin_int_conc <- skin_int/organ_v_int_skin
  skin_cell_conc <- skin_cell/organ_v_cell_skin
  # diffusion gradients
  skin_intdiff <- fup * Pend * SA_pls_int_skin * (skin_vas_conc - skin_int_conc/(PC_kip_skin/BP))
  skin_celldiff <- Porg * SA_int_cell_skin * (skin_int_conc * PC_kiw_skin - skin_cell_conc * PC_kcw_skin)
  # differential equations
  d/dt(skin_vas) <- skin_inflow - skin_outflow - skin_intdiff
  d/dt(skin_int) <- skin_intdiff - skin_celldiff
  d/dt(skin_cell) <- skin_celldiff
  
  bf_tot <- bf_tot + organ_bf_skin
  outflow_tot <- outflow_tot + skin_outflow
  
  
  # 7) fat ---------------------------------------------------------------------
  fat_inflow <- organ_bf_fat * (arterial_blood/organ_v_arterial_blood)
  fat_outflow <- organ_bf_fat * (fat_vas/organ_v_vas_fat)
  # concentrations
  fat_vas_conc <- fat_vas/organ_v_vas_fat
  fat_int_conc <- fat_int/organ_v_int_fat
  fat_cell_conc <- fat_cell/organ_v_cell_fat
  # diffusion gradients
  fat_intdiff <- fup * Pend * SA_pls_int_fat * (fat_vas_conc - fat_int_conc/(PC_kip_fat/BP))
  fat_celldiff <- Porg * SA_int_cell_fat * (fat_int_conc * PC_kiw_fat - fat_cell_conc * PC_kcw_fat)
  # differential equations
  d/dt(fat_vas) <- fat_inflow - fat_outflow - fat_intdiff
  d/dt(fat_int) <- fat_intdiff - fat_celldiff
  d/dt(fat_cell) <- fat_celldiff
  
  bf_tot <- bf_tot + organ_bf_fat
  outflow_tot <- outflow_tot + fat_outflow
  
  
  ### splanchnic organs equations c("stomach", "spleen", "pancreas", "gut")
  
  # 8) stomach -----------------------------------------------------------------
  stomach_inflow <- organ_bf_stomach * (arterial_blood/organ_v_arterial_blood)
  stomach_outflow <- organ_bf_stomach * (stomach_vas/organ_v_vas_stomach)
  # concentrations
  stomach_vas_conc <- stomach_vas/organ_v_vas_stomach
  stomach_int_conc <- stomach_int/organ_v_int_stomach
  stomach_cell_conc <- stomach_cell/organ_v_cell_stomach
  # diffusion gradients
  stomach_intdiff <- fup * Pend * SA_pls_int_stomach * (stomach_vas_conc - stomach_int_conc/(PC_kip_stomach/BP))
  stomach_celldiff <- Porg * SA_int_cell_stomach * (stomach_int_conc * PC_kiw_stomach - stomach_cell_conc * PC_kcw_stomach)
  # differential equations
  d/dt(stomach_vas) <- stomach_inflow - stomach_outflow - stomach_intdiff
  d/dt(stomach_int) <- stomach_intdiff - stomach_celldiff
  d/dt(stomach_cell) <- stomach_celldiff
  
  bf_tot <- bf_tot + organ_bf_stomach
  outflow_spl <- outflow_spl + stomach_outflow
  
  
  # 9) spleen ------------------------------------------------------------------
  spleen_inflow <- organ_bf_spleen * (arterial_blood/organ_v_arterial_blood)
  spleen_outflow <- organ_bf_spleen * (spleen_vas/organ_v_vas_spleen)
  # concentrations
  spleen_vas_conc <- spleen_vas/organ_v_vas_spleen
  spleen_int_conc <- spleen_int/organ_v_int_spleen
  spleen_cell_conc <- spleen_cell/organ_v_cell_spleen
  # diffusion gradients
  spleen_intdiff <- fup * Pend * SA_pls_int_spleen * (spleen_vas_conc - spleen_int_conc/(PC_kip_spleen/BP))
  spleen_celldiff <- Porg * SA_int_cell_spleen * (spleen_int_conc * PC_kiw_spleen - spleen_cell_conc * PC_kcw_spleen)
  # differential equations
  d/dt(spleen_vas) <- spleen_inflow - spleen_outflow - spleen_intdiff
  d/dt(spleen_int) <- spleen_intdiff - spleen_celldiff
  d/dt(spleen_cell) <- spleen_celldiff
  
  bf_tot <- bf_tot + organ_bf_spleen
  outflow_spl <- outflow_spl + spleen_outflow
  
  # 10) pancreas ---------------------------------------------------------------
  pancreas_inflow <- organ_bf_pancreas * (arterial_blood/organ_v_arterial_blood)
  pancreas_outflow <- organ_bf_pancreas * (pancreas_vas/organ_v_vas_pancreas)
  # concentrations
  pancreas_vas_conc <- pancreas_vas/organ_v_vas_pancreas
  pancreas_int_conc <- pancreas_int/organ_v_int_pancreas
  pancreas_cell_conc <- pancreas_cell/organ_v_cell_pancreas
  # diffusion gradients
  pancreas_intdiff <- fup * Pend * SA_pls_int_pancreas * (pancreas_vas_conc - pancreas_int_conc/(PC_kip_pancreas/BP))
  pancreas_celldiff <- Porg * SA_int_cell_pancreas * (pancreas_int_conc * PC_kiw_pancreas - pancreas_cell_conc * PC_kcw_pancreas)
  # differential equations
  d/dt(pancreas_vas) <- pancreas_inflow - pancreas_outflow - pancreas_intdiff
  d/dt(pancreas_int) <- pancreas_intdiff - pancreas_celldiff
  d/dt(pancreas_cell) <- pancreas_celldiff
  
  bf_tot <- bf_tot + organ_bf_pancreas
  outflow_spl <- outflow_spl + pancreas_outflow
  
  # 11) gut --------------------------------------------------------------------
  gut_inflow <- organ_bf_gut * (arterial_blood/organ_v_arterial_blood)
  gut_outflow <- organ_bf_gut * (gut_vas/organ_v_vas_gut)
  # concentrations
  gut_vas_conc <- gut_vas/organ_v_vas_gut
  gut_int_conc <- gut_int/organ_v_int_gut
  gut_cell_conc <- gut_cell/organ_v_cell_gut
  # diffusion gradients
  gut_intdiff <- fup * Pend * SA_pls_int_gut * (gut_vas_conc - gut_int_conc/(PC_kip_gut/BP))
  gut_celldiff <- Porg * SA_int_cell_gut * (gut_int_conc * PC_kiw_gut - gut_cell_conc * PC_kcw_gut)
  # differential equations
  d/dt(gut_vas) <- gut_inflow - gut_outflow - gut_intdiff
  d/dt(gut_int) <- gut_intdiff - gut_celldiff
  d/dt(gut_cell) <- gut_celldiff
  
  bf_tot <- bf_tot + organ_bf_gut
  outflow_spl <- outflow_spl + gut_outflow
  
  
  # 12) liver ------------------------------------------------------------------
  liver_inflow <- outflow_spl + organ_bf_liver * (arterial_blood/organ_v_arterial_blood)
  liver_outflow <- (organ_bf_liver + organ_bf_gut + organ_bf_pancreas + organ_bf_spleen + organ_bf_stomach) * (liver_vas/organ_v_vas_liver)
  # concentrations
  liver_vas_conc <- liver_vas/organ_v_vas_liver
  liver_int_conc <- liver_int/organ_v_int_liver
  liver_cell_conc <- liver_cell/organ_v_cell_liver
  # diffusion gradients
  liver_intdiff <- fup * Pend * SA_pls_int_liver * (liver_vas_conc - liver_int_conc/(PC_kip_liver/BP))
  liver_celldiff <- Porg * SA_int_cell_liver * (liver_int_conc * PC_kiw_liver - liver_cell_conc * PC_kcw_liver)
  # clearance
  clear_liv <- CLh * fup * liver_cell_conc/PC_kcp_liver
  # differential equations
  d/dt(liver_vas) <- liver_inflow - liver_outflow - liver_intdiff
  d/dt(liver_int) <- liver_intdiff - liver_celldiff
  d/dt(liver_cell) <- liver_celldiff - clear_liv + input_GI
  
  bf_tot <- bf_tot + organ_bf_liver
  outflow_tot <- outflow_tot + liver_outflow
  
  # 13) venous blood equation --------------------------------------------------
  outflow_ven <- bf_tot * venous_blood/organ_v_venous_blood
  d/dt(venous_blood) <- outflow_tot - outflow_ven
  
  # 14) lung equation ----------------------------------------------------------
  lungs_inflow <- outflow_ven # blood inflow to lungs
  lungs_outflow <- bf_tot * (lungs_vas/organ_v_vas_lungs) # blood outflow from lungs
  # concentrations
  lungs_vas_conc <- lungs_vas/organ_v_vas_lungs # vascular concentration
  lungs_int_conc <- lungs_int/organ_v_int_lungs # interstitial concentration
  lungs_cell_conc <- lungs_cell/organ_v_cell_lungs # intracellular concentration
  # diffusion gradients
  lungs_intdiff <- fup * Pend * SA_pls_int_lungs * (lungs_vas_conc - lungs_int_conc/(PC_kip_lungs/BP)) # vascular to interstitial diffusion
  lungs_celldiff <- Porg * SA_int_cell_lungs * (lungs_int_conc * PC_kiw_lungs - lungs_cell_conc * PC_kcw_lungs) # interstitial to intracellular diffusion
  # differential equations
  d/dt(lungs_vas) <- lungs_inflow - lungs_outflow - lungs_intdiff # vascular compartment
  d/dt(lungs_int) <- lungs_intdiff - lungs_celldiff # interstitial compartment
  d/dt(lungs_cell) <- lungs_celldiff # intracellular compartment
  
  # 15) arterial blood equation ------------------------------------------------
  d/dt(arterial_blood) <- lungs_outflow - bf_tot * arterial_blood/organ_v_arterial_blood
  
  # 16) sink compartments ------------------------------------------------------
  d/dt(sink_CLh) <- clear_liv # liver sink
  d/dt(sink_CLr) <- clear_kid # kidney sink
  
  plasma_conc <- venous_blood/BP/organ_v_venous_blood
})