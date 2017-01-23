#!/bin/csh
set echo
##################################
#
# Snow Albedo Feedback (SAF) script
# Release version 1.0
# September 2014
#
# <<Acknowledgement>>
# Created by Chris Fletcher (UWaterloo)
# with Toby Woerthle and Chad Thackeray.
#
# This script is provided freely to the community for the purpose of reproducibility and 
# model evaluation.  It can be modified and distributed freely, but we ask that proper
# citation of the article below be included with any future publication that uses the code.
#
# Reference: Fletcher, C.G., C.W. Thackeray, and T.M. Burgers (2014), Evaluating biases in simulated snow
#                albedo feedback in two generations of climate models, J. Geophys. Res. Atmos., Accepted.
# ** If you use or modify the code, please include a reference to this article **
#
# <<Background>>
# This script outputs NET SAF (dalpha/dT), and its components SNC and TEM, 
# computed for spatial grids of any input dataset containing albedo, snow cover fraction, temperature.
# The script is currently set up to work with monthly mean files in NetCDF format,
# but it can be modified to work with other formats, and data on different time frequencies.
#
# <<File structure>>
# The code assumes a file structure by default.
# INPUT_DATA = /path/to/dir/netcdf_files
# ROOT_DIR = /path/to/SAF_FTB2014_v1.0
# AUX_DIR = $ROOT_DIR/aux_data
# OUT_DIR = $ROOT_DIR/results
# WORK_DIR = $ROOT_DIR/work
# 
# <<Auxillary Files>>
# These files are required, and provided in the aux_data directory
# minimum_scf_extent_mask.nc -- EASE grid projection onto regular 2.5x2.5-deg lat-lon grid.
# dswrf.ntat.climo.mon.mean.nc -- NCEP-II-DOE climatological monthly mean incoming SW radiation at TOA
# filename.nc : land mask (minimum extent mask)
# filename.nc : TOA insolation for computing weights
#
# <<Input data requirements>>
# The code requires individual files containing grids of:
#    surface albedo (albedo == rsus/rsds), snow cover fraction (scf), and near-surface temperature (tas).
# All input files for a particular model, or observational product, are assumed to be on the same grid.
# If this is not the case, the user will need to interpolate all files to a common grid first.
#
# The input files are assumed be in the directory called INPUT_DATA. 
# File naming conventions are assumed to follow the CMIP5 standard, e.g.:
#     ${variable}_Amon_${model}_${expt}_r1i1p1_${time_period_str}.nc
# If your data are for a different realization than r1i1p1, you will need to modify the realization variable below.
#
# <<Units>>
# Input SCF and albedo data are assumed to be decimal fraction values (*not* percentages)
#
# <<Notes and modifications>>
# Sep 2014: This version: uses monthly maximum albedo to calculate alpha_snow
#
# << CONTACT >>
# For all code-related enquiries contact Chris Fletcher
# email: chris.fletcher@uwaterloo.ca
##################################


######################################################################################
 ################  ((((( START of user modification section ))))) ##################
###################################################################################### 
# Setup input/output directories: 
# Change these variables to match your data/preference
# The INPUT_DATA directory should contain monthly mean 
# albedo, snow cover fraction or swe, tas, and a land-sea mask
# for as many years as specified below.
 set INPUT_DATA = /path/to/CMIP5/data
 set convert_snw = 0  # change this option to 1 to convert SWE to SCF on the fly.
# Models - change to model of choice, or create model_list
 set model = "CanESM2"  # "CCSM4 CESM1-CAM5"
 set expt = historical
# Set years of study - should be changed to match years available in your data files
 set start_year = 1980
 set end_year = 2005
 set time_period_str = "${start_year}_${end_year}"
# Are we using an ensemble mean (default -- needs to be created by the user) or a single realization?
 set realization = "emi1p1"  ## em = ensemble mean; can be changed to "r1i1p1" "r2i1p1" etc.
 set albedo_file = $INPUT_DATA/albedo_Amon_${model}_${expt}_${realization}_${time_period_str}.nc
 set tas_file = $INPUT_DATA/tas_Amon_${model}_${expt}_${realization}_${time_period_str}.nc
# This is the scf file (to be created if convert_snw > 0)
 set scf_file = $INPUT_DATA/scf_LImon_${model}_${expt}_${realization}_${time_period_str}.nc
# This is the snw file
 set snw_file = $INPUT_DATA/snw_LImon_${model}_${expt}_${realization}_${time_period_str}.nc
# This is the land mask file
 set lmask_file = $INPUT_DATA/sftlf.historical.${model}.ea.fx.nc

# ROOT_DIR is the path to the aux_data, and also where the working dir and results will be stored by default.
 set ROOT_DIR = /path/to/SAF
 set outdir_root = $ROOT_DIR/results
# Set work directory
 set work = $ROOT_DIR/work
## Set locations of auxillary files
 set obs_mask = $ROOT_DIR/aux_data/EASE_grid_mask.nc
 set insol_file = $ROOT_DIR/aux_data/dswrf.ntat.climo.mon.mean.nc

# The seasonal transitions are computed as the difference between mon1 and mon2
# Here we can define the first and last transition to be analysed:
# startmon is the first mon1 (e.g., 3 (March), 4 (April) or 5 (May) )
 setenv startmon 3
# endmon is the last mon1 (e.g., 3 (March), 4 (April) or 5 (May) )
 setenv endmon 5
# Define analysis region (degrees; default: 30N-90N, 0-360E)
 set latmin = 30
 set latmax = 90
 set lonmin = 0
 set lonmax = 360
# Snow Threshold - the minimum amount of snow cover required to include a grid cell in the analysis
 set snow_thresh = "0.25"
######################################################################################
  ################  ((((( END of user modification section ))))) ##################
###################################################################################### 
#
#
#   WARNING: modifications not recommended below here.
#
#
#
# set up paths and time vars
if (! -d $outdir_root ) mkdir -p $outdir_root
mkdir -p $work
cd $work
set dates = "${start_year}-01-01,${end_year}-12-31"
set regstr = "${lonmin},${lonmax},${latmin},${latmax}"
# We can opt to set missing values to 0 in the calculation of K terms, but defauly is not to.
 set misstocstr = "" ## "-setmisstoc,0.0"
 set flagstr = "nomisstoc" ## "misstoc"

# Loop over models from model_list 
 set model_list = ( $model )
# outdir pattern is the identifier for all SAF terms computed using this script version 
 set outdir_pattern = ${snow_thresh}_${flagstr}
# Set the snow cover threshold - we will calculate SAF all cells with more than 0.25 SCF (previously used 0.1)
 set snow_present_thresh = $snow_thresh 

 echo ">>>>>>   PROCESSING: $model"

 # Locate land mask for this model
   if (! -f $lmask_file ) exit 1
   cdo sellonlatbox,${regstr} $lmask_file lmask.nc || exit 1
   set lmask_file = lmask.nc
 #
 # +++++++++++++++++++++++++++++++++++++++++++++
 # +++++++++++++++++++++++++++++++++++++++++++++
 # 	Multi-stage SAF calculation begins here
 # +++++++++++++++++++++++++++++++++++++++++++++
 # 1.  Read in input files for albedo tas and scf
 # Get native-resolution grids of alpha_s, Ts and Scf for all 12 months (climatological mean, NH30 domain + landmask).
 # +++++++++++++++++++++++++++++++++++++++++++++
 foreach var (albedo tas scf)
  #
  if ($var == albedo) set infile = $albedo_file
  if ($var == tas) set infile = $tas_file
  if ($var == scf) then
   set infile = $scf_file
   # Convert SNW to SCF on the fly if needed:   
   # Begin conversion if requested
   if ($convert_snw == 1) then
    echo ">>> Converting snw to scf"
    set gridstr = ""
    # Added check for using selgrid 
    set flag = `cdo infos $snw_file | grep "1 : generic"`
    if ($#flag != 0 ) set gridstr = "-selgrid,2"
    cdo -b f64 lec,0 $gridstr $snw_file msk0.nc   
    cdo -b f64 gtc,60 $gridstr $snw_file msk1.nc   
    cdo -b f64 divc,60 $gridstr $snw_file tmp1.nc
    # Create scf from snow mass and merge into same file     
    # Procedure using CDO: make 3 mask fields: msk1.nc = snw > 60, msk2.nc = 0 < snw < 60 and msk3.nc = snw = 0    
    # Then for msk1.nc apply scf = 1.0, and for msk0.nc apply scf = 0.0 
    cdo -b f64 setmisstoc,1.0 -ifnotthen msk1.nc tmp1.nc tmp2.nc || exit 1
    cdo -b f64 setname,scf -setmisstoc,0.0 -ifnotthen msk0.nc tmp2.nc $infile || exit 1
   endif
   if (! -f $infile) then
    echo "FATAL: scf file does not exist"
    echo "Either data not available, or convert from SWE by selecting convert_snw = 1."
    exit 1
   endif
  endif # are we SCF?
  # Processing the data into monthly climatologies for NH30 domain
  if ($var == albedo) then
    cdo ymonmean -selname,$var -setmissval,-999.0 -seldate,$dates -sellonlatbox,${regstr} $infile tmp.albedo.nc
  else if ($var == scf) then
    set gridstr = ""
  # Added check for using selgrid 
  set flag = `cdo infos $infile | grep "1 : generic"`
  if ($#flag != 0 ) set gridstr = "-selgrid,2"
    cdo ymonmean -selname,$var -setmissval,-999.0 -seldate,$dates -sellonlatbox,${regstr} $gridstr $infile tmp.scf.nc
  else if ($var == tas) then
    cdo ymonmean -setmissval,-999.0 -seldate,$dates -sellonlatbox,${regstr} -selname,tas $infile tmp.tas.nc
  endif
  # if albedo then multiply by 100% to get a percentage value
  if ($var == albedo) then 
   cdo mulc,100 tmp.albedo.nc tmp.nc && \mv tmp.nc tmp.albedo.nc
  endif
 end  # get var loop

 # ++++++++++++++++++++++++++++++++++++++++
 # 2. Compute max_albedo and alpha_land -- these are time-invariant quantities
 # max_albedo is the maximum albedo over all months (monmean instead of ymonmean)
 set allmons_albedo_file = $albedo_file
 cdo mulc,100 -timmax -selmon,2,3,4,5,6 -seldate,$dates -sellonlatbox,${regstr} -setmissval,-999.0 $allmons_albedo_file tmp.max_albedo.nc

 # ++++++++++++++++++++++++++++++++++++++++
 # Compute alpha_snow and alpha_land using Qu and Hall (2013) method:
 # i. Define threshold scf >= 0.9 (90%) / < 0.1 (10%) for "snow-covered" / "snow-free"
 # ii. Create mask from monthly climo (12-month) scf files using threshold
 cdo gec,0.9 tmp.scf.nc mask.scf90.nc
 cdo ltc,0.1 tmp.scf.nc mask.scf10.nc 
 # iii. Average albedo from monthly climo (same 12-months) for non-missing cells in the mask
 cdo setname,qh_asnow -timmean -ifthen mask.scf90.nc tmp.albedo.nc tmp.qh_asnow.nc
 cdo setname,qh_aland -timmean -ifthen mask.scf10.nc tmp.albedo.nc tmp.qh_aland.nc
 # iv. Save these new quantities: qh_aland_scf10 and qh_asnow_scf90
 # done later, after interpolating and applying obs mask. 

 # ++++++++++++++++++++++++++++++++++++++++
 # alpha_land is the surface albedo for the first month after snow melts
 @ mon = 1
 @ ind = 10
 while ($mon <= 11)  # final mon + 1 = 12
  @ mon2 = $mon + 1
  @ mon3 = $mon + 2
  cdo selmon,$mon tmp.scf.nc this_scf_before.nc
  cdo selmon,$mon2 tmp.scf.nc this_scf_after.nc
  cdo selmon,$mon2,$mon3 tmp.albedo.nc this_albedo.nc
  # 
  cdo gtc,0.1 this_scf_before.nc mask.nc	
  cdo ifthen mask.nc this_scf_after.nc tmp2.nc
  cdo ltc,0.1 tmp2.nc mask2.nc
 # Now mask albedo grids with mask2.nc
  cdo ifthen mask2.nc this_albedo.nc this_albedo_msk.nc
 # Now average to get this month's values of alpha_land
  cdo setname,alpha_land -timmean this_albedo_msk.nc tmp.$ind.nc 
  @ mon ++
  @ ind ++
 end 
 # Finally, take sum of the first two non-missing months = alpha_land
 cdo -O enssum tmp.??.nc tmp.alpha_land.nc

 # ++++++++++++++++++++++++++++++++++++++++
 # 3. For March, April, May and June only: compute alpha_snow
 # ++++++++++++++++++++++++++++++++++++++++ 
 @ mon = 3
 while ($mon <= 6)
  # Merge alpha_land, alpha_s and snow cover for this month
  cdo selmon,$mon tmp.albedo.nc tmp2.nc
  cdo selmon,$mon tmp.scf.nc tmp3.nc  
  \rm -rf tmp.allvars.nc
  cdo -O merge tmp2.nc tmp.alpha_land.nc tmp3.nc tmp.allvars.nc
  cdo setmon,$mon tmp.allvars.nc tmp.nc && \mv tmp.nc tmp.allvars.nc
  
  # Expression for alpha_snow (albedo in %, scf as fraction):
  cdo -b 64 setvrange,0,100 -setmissval,-999 -expr,"alpha_snow=(albedo-((1-scf)*alpha_land))/scf" tmp.allvars.nc tmp.alsnw.0$mon.nc

  # Here we add the "checks" imposed by F09 for small SCF 
  # (I will impose for all cells, not just those where SCF < 10%):
  if ($mon <= 3) then
   # alpha_snow(March) = min(alpha_snow(march), max_snow_albedo)
   cdo min tmp.alsnw.0$mon.nc tmp.max_albedo.nc tmp.nc
   \mv tmp.nc tmp.alsnw.0$mon.nc 
  else
   # Cases other than march: we compare this month with previous and take the minimum
   # alpha_snow(April) = min(alpha_snow(april), alpha_snow(march))
   # alpha_snow(May) = min(alpha_snow(may), alpha_snow(april))
   @ mm1 = $mon - 1
   cdo min tmp.alsnw.0$mon.nc tmp.alsnw.0$mm1.nc tmp.nc
   \mv tmp.nc tmp.alsnw.0$mon.nc
  endif
  # Final check: alpha_snow = max(alpha_snow, alpha_land), i.e. snow albedo cannot be < land albedo!
  cdo max tmp.alsnw.0$mon.nc tmp.alpha_land.nc tmp.nc
  \mv tmp.nc tmp.alsnw.0$mon.nc
  @ mon ++
 end 
 if (-f tmp.alpha_snow.nc ) \rm -rf tmp.alpha_snow.nc
 cdo mergetime tmp.alsnw.??.nc tmp.alpha_snow.nc 
 #

 # ++++++++++++++++++++++++++++++++++++++++
 # 4. Loop over monthly transitions and compute k4(NET), k1k2(SNC) and k3(TEM)
 # ++++++++++++++++++++++++++++++++++++++++
 @ mon1 = $startmon
 set monstr = "MA"
 while ($mon1 <= $endmon)
  if ($mon1 == 4) set monstr = "AM"
  if ($mon1 == 5) set monstr = "MJ"
  if ($mon1 > 5) exit 1
  @ mon2 = $mon1 + 1
 
# output location for k-coefficients interpolated to the obs grid
 set outdir = $outdir_root/${monstr}_${outdir_pattern}
  if (! -d $outdir ) mkdir -p $outdir


# Define the snow mask based on March, for all time periods.
 cdo gec,0.0 -selmon,3 tmp.scf.nc tmp.scf_gt_0_mask.nc

#
 # ++++++++++++++++++++++++++++++++++++++++
 # 4a. Insolation weight all albedo fields using observed TOA insolation weighting 
 # ++++++++++++++++++++++++++++++++++++++++
 cdo -b 64 -setmissval,-999.0 -interpolate,tmp.albedo.nc -setname,albedo -timmean -selmon,$mon1,$mon2 -sellonlatbox,${regstr} $insol_file tmp1.nc
 # land mask insolation (weighting only relative to land)
 cdo ifthen $lmask_file tmp1.nc tmp.nc && \mv tmp.nc tmp1.nc
 cdo fldmean tmp1.nc tmp2.nc
 cdo div tmp1.nc -enlarge,tmp1.nc tmp2.nc tmp.obs.insol_toa.$monstr.nc
 # multiply tmp.albedo, tmp.alpha_land and tmp.alpha_snow by insol weight field
 cdo mul tmp.albedo.nc tmp.obs.insol_toa.$monstr.nc tmp.albedo.wgt.nc
 cdo mul tmp.alpha_snow.nc tmp.obs.insol_toa.$monstr.nc tmp.alpha_snow.wgt.nc
 cdo mul tmp.alpha_land.nc tmp.obs.insol_toa.$monstr.nc tmp.alpha_land.wgt.nc
 # new: weight QH_asnow and QH_aland by insolation
 cdo mul tmp.qh_asnow.nc tmp.obs.insol_toa.$monstr.nc tmp.qh_asnow.wgt.nc
 cdo mul tmp.qh_aland.nc tmp.obs.insol_toa.$monstr.nc tmp.qh_aland.wgt.nc

 # ++++++++++++++++++++++++++++++++++++++++
 # 4b. Here we calculate the K terms
 # ++++++++++++++++++++++++++++++++++++++++
 # Calculate k4 (NET SAF)
 cdo sub -selmon,$mon2 tmp.tas.nc -selmon,$mon1 tmp.tas.nc tmp.delta_tas.nc    
 cdo ifthen $lmask_file tmp.delta_tas.nc tmp1.nc || exit 1
 cdo enlarge,tmp1.nc -fldmean -sellonlatbox,${regstr} tmp1.nc tmp.delta_ts_nh.nc
 # next compute Delta alpha_s and divide by tmp.delta_ts_nh, then mask, then area-avg
 cdo sub -selmon,$mon2 tmp.albedo.wgt.nc -selmon,$mon1 tmp.albedo.wgt.nc tmp.delta_albedo.nc    
 cdo -setname,k4 -div tmp.delta_albedo.nc tmp.delta_ts_nh.nc k4.$model.nc            
 #
 # k1k2 (SNC COMPONENT)
 # k1 is Delta Scf / <Delta_Ts>_NH
 cdo sub -selmon,$mon2 tmp.scf.nc -selmon,$mon1 tmp.scf.nc tmp.delta_scf.nc
 cdo setname,k1 -mulc,100 -div tmp.delta_scf.nc tmp.delta_ts_nh.nc k1.$model.nc   
 # k2 is mean albedo_snow - albedo_land
 cdo divc,2 -add -selmon,$mon2 tmp.alpha_snow.wgt.nc -selmon,$mon1 tmp.alpha_snow.wgt.nc tmp.mean_alpha_snow.wgt.nc
 cdo setvrange,0,100.0 -setname,k2 -sub tmp.mean_alpha_snow.wgt.nc tmp.alpha_land.wgt.nc k2.$model.nc
 cdo -setname,k1k2 -divc,100 -mul k2.$model.nc k1.$model.nc k1k2.$model.nc 

 # k3 (TEM component):
 # First, k3 as a residual (k4 - k1k2) [method of Fernandes et al. 2009, GRL] 
 cdo -setname,k3_res -sub k4.$model.nc k1k2.$model.nc k3_res.$model.nc
 
 # We use alpha_snow to compute k3 explicitly for this month pair
 # compute mean snow, Delta alpha_snow, multiply them and divide by dTs NH
 cdo divc,2 -add -selmon,$mon2 tmp.scf.nc -selmon,$mon1 tmp.scf.nc tmp.mean_snow.nc
 cdo sub -setname,delta_alpha_snow -selmon,$mon2 tmp.alpha_snow.wgt.nc -selmon,$mon1 tmp.alpha_snow.wgt.nc tmp.delta_alpha_snow.nc
 cdo mul tmp.mean_snow.nc tmp.delta_alpha_snow.nc tmp2.nc
 # Divide by NH avg Delta Ts to create TEM:
 cdo setname,k3 -div tmp2.nc tmp.delta_ts_nh.nc k3.$model.nc
 
 # ++++++++++++++++++++++++++++++++++++++++ 
 # 4c. Regrid all k fields to obs grid, apply obs land mask, then compute area-averages.
 # We need to save monthly mean snow albedo, Ts, S and alpha_s.
 # -- first save important diagnostics, like the different albedos
 set save_vars = ( max_albedo alpha_land.wgt delta_albedo delta_scf delta_tas mean_alpha_snow.wgt delta_alpha_snow albedo.wgt scf tas qh_aland.wgt qh_asnow.wgt mean_snow)
 foreach var ($save_vars)
  # Need to do a copy here, so that the constant vars (alpha_land) are preserved
  \cp tmp.$var.nc $var.$model.nc
 end
 foreach k (k1 k2 k1k2 k3_res k3 k4 alpha_land.wgt delta_albedo delta_scf delta_tas mean_alpha_snow.wgt delta_alpha_snow albedo.wgt scf tas qh_aland.wgt qh_asnow.wgt mean_snow)
  # Add explicit condition limiting domain to S > 0 cells.               
  # If misstocstr is set, then cells where S < threshold are set to zero.  
  # This means they will not be missing, and contribute zero to average SAF.                 
  # We need to reapply the land-sea mask so that ocean = missing,and not = 0 (this messes up interpolation).
  # -- APPLY ONLY TO k-fields, for all others we set cells with S < thresh to MISSING.  
  if ($k == "k1" || $k == "k2" || $k == "k1k2" || $k == "k3" || $k == "k3_res" || $k == "k4" || ) then
   cdo ifthen $lmask_file $misstocstr -ifthen tmp.scf_gt_0_mask.nc $k.$model.nc tmp1.nc || exit 1
   \mv tmp1.nc $k.$model.nc
  else
   cdo ifthen $lmask_file -ifthen tmp.scf_gt_0_mask.nc $k.$model.nc tmp1.nc
   \mv tmp1.nc $k.$model.nc
  endif
  cdo interpolate,$obs_mask $k.$model.nc tmp1.nc
   # The following should now be the only land-sea mask applied to the fields.
   # This ensures we have the correct missing data area for obs and all model grids.
  cdo ifthen $obs_mask tmp1.nc $outdir/$k.$model.EASE_grid.nc
  cdo fldmean $outdir/$k.$model.EASE_grid.nc $outdir/$k.$model.EASE_grid.aa.nc
 end
 #
 @ mon1 ++
 end # loop over months

 # clean up work dir
 set files = `ls tmp*.nc *mask*.nc tas_20.*.nc this_*.nc`
 \rm -rf $files


