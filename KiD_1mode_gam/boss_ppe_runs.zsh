#!/bin/zsh

# config of the run
mps=("boss_4m_3069")
config_name="condevpcoal"
caselist=(102) #(101 102 103 105 106 107)
case_num=${#caselist[@]}

KiD_path="/global/homes/a/arthurhu/KiD_AMP/KiD_1mode_gam/"

nikki=$(date +'%Y-%m-%d')
s_sample_dist="lhs"
# lhs_path="/Users/arthurhu/Library/Mobile Documents/com~apple~CloudDocs/storage/postdoc/KiD_AMP/KiD_1mode_gam/lhs_nc"
lhs_path="${KiD_path}lhs_nc"
echo $lhs_path
# param_val_fpath="/Users/arthurhu/Library/Mobile Documents/com~apple~CloudDocs/storage/postdoc/CloudBOSS/param_consolid_updated.csv"
param_val_fpath="/global/homes/a/arthurhu/Cloud_BOSS/param_consolid_updated.csv"
# custom_dens_path="/Users/arthurhu/github/BOSS_PPE/MCMC_posterior/condevp_withcoal_r1_param_density_RWM.csv"
# custom_bins_path="/Users/arthurhu/github/BOSS_PPE/MCMC_posterior/condevp_withcoal_r1_param_bins_RWM.csv"
custom_dens_path="/global/homes/a/arthurhu/BOSS_PPE/MCMC_posterior/fullmic_r1_param_density_RWM.csv"
custom_bins_path="/global/homes/a/arthurhu/BOSS_PPE/MCMC_posterior/fullmic_r1_param_bins_RWM.csv"

l_ppe_nevp=".true."
l_ppe_condevp=".true."
l_ppe_coal=".true."
l_ppe_sed=".false."

# initial condition for all cases
cim3=0.
cim0=0.

rim3=0.
rim0=0.

rs_dm=0. # mean-mass diameter (m), ignores the case once this is non-zero
rs_N=0. # number mixing ratio (#/kg)
ztop=6000. # top of the domain
zcb=1000. # cloud base height
zct=3000. # cloud bottom height
t1=3600.
t2=900.

# switches for nucleation/condensation, collision, sedimentation, and advection
l_nuc_cond_s=1
l_coll_s=1
l_sed_s=0
l_adv_s=1

# []==if, &&==then, ||=else
[ $l_nuc_cond_s -eq 1 ] && l_nuc_cond_f='.true.' || l_nuc_cond_f='.false.'
[ $l_coll_s -eq 1 ] && l_coll_f='.true.' || l_coll_f='.false.'
[ $l_sed_s -eq 1 ] && l_sed_f='.true.' || l_sed_f='.false.'

if [ "$l_adv_s" = 1 ]; then
   l_adv='.true.'
   l_noadv_qv='.false.'
   l_noadv_hyd='.false.'
else
   l_adv='.false.'
   l_noadv_qv='.true.'
   l_noadv_hyd='.true.'
fi

isp_c=4
isp_r=4

Na_min=$1
Na_max=$2
w_min=$3
w_max=$4

# # reset oscillation time based on updraft speed to prevent overshooting
# if [[ $((($ztop-$zct)/$w_max)) -lt $t2 && $l_adv_s -eq 1 ]]; then
#   t2=$((($ztop-$zct)/$w_max))
#   t1=$(($t2*4))
# fi

for ((imp=1; imp<=${#mps[@]}; imp=imp+1))
do
mp=${mps[$imp]}
echo $mp
   if [[ $mp = 'AMP4m_SBM' ]]; then
      isp_c=2  # shape parameter for cloud
      isp_r=2  # shape parameter for rain
      imc1=4 # II moment for cloud
      imc2=5 # III moment for cloud
      imr1=4 # II moment for rain
      imr2=5 # III moment for rain
   elif [[ $mp = *AMP4m_TAU* ]]; then
      imc1=${mp[-2]}
      imc2=${mp[-1]}
      imr1=${mp[-2]}
      imr2=${mp[-1]}
   elif [[ $mp = 'AMP2m_SBM' ]]; then
      isp_c=$isp_c  # shape parameter for cloud
      isp_r=$isp_r  # shape parameter for rain
      imc1=0 # II moment for cloud
      imc2=0 # III moment for cloud
      imr1=0 # II moment for rain
      imr2=0 # III moment for rain
   elif [[ $mp = 'AMP2m_TAU' ]]; then
      imc1=0 # II moment for cloud
      imc2=0 # III moment for cloud
      imr1=0 # II moment for rain
      imr2=0 # III moment for rain
   elif [[ $mp = *boss_4m* ]]; then
      imc1=${mp[-2]}
      imc2=${mp[-1]}
      imr1=${mp[-2]}
      imr2=${mp[-1]}
   elif [[ $mp = *boss_2m* ]]; then
      imc1=0
      imc2=0
      imr1=0
      imr2=0
   fi

   if [[ $mp = *ML_* ]]; then
     l_use_nn='.true.'
   else
     l_use_nn='.false.'
   fi

   if [[ $mp = *pres* ]]; then
     Pdep='pres_'
   else
     Pdep=''
   fi

   if [[ $mp = *boss* ]] then
     mp_id='boss'
     bintype='tau'
     isp_c=$isp_c
     isp_r=$isp_r
     ampORbin=''
   else
     mp_id='amp'
   fi

   if [[ $mp = *4m* ]]; then
     nhm='4,4'
     nhb='1,1'
     n_cat=1
     momnames="'M1','M2','M3','M4'"
   elif [[ $mp = *2m* ]]; then
     nhm='2,2'
     nhb='1,1'
     n_cat=2
     momnames="'M1','M2'"
   elif [[ $mp = 'BIN_SBM' ]]; then
     nhm='1,1'
     nhb='33,33'
   elif [[ $mp = 'BIN_TAU' ]]; then
     nhm='2,2'
     nhb='34,1'
     momnames="'M1','M2'"
     imc1=6
     imc2=9
   fi

   if [[ $mp = *AMP* ]]; then
     ampORbin='amp'
   elif [[ $mp = *BIN* ]]; then
     ampORbin='bin'
   fi

   if [[ $mp = *TAU* ]]; then
     bintype='tau'
     isp_c=$isp_c
     isp_r=$isp_r
   elif [[ $mp = *SBM* ]]; then
     bintype='sbm'
   fi


   # outdir=~/research/KiD_output/$nikki/$config_name/${mp}_ens${5}/
   outdir=/pscratch/sd/a/arthurhu/KiD_output/$nikki/$config_name/${mp}_ens${5}/
   for ((ic=1; ic<=case_num; ic++))
   do
      if [[ ${caselist[ic]} -gt 104 ]] && [[ ${caselist[ic]} -lt 200 ]]
      then
         zc=0
      else
         zc="$ztop,$zcb,$zct"
      fi
      if [ ! -d $outdir ]; then
         mkdir -p $outdir
      fi
      echo ${config_name}_${mp}
      nml_fn=${KiD_path}namelists/jobnml/${config_name}_${mp}_ens${5}.nml
      cat > $nml_fn << END
&mphys
! hydrometeor names
h_names='cloud','rain'

!Moment names
mom_names='M1','M2','M3','M4','M5','M6'

!Initial shape parameter
h_shape=${isp_c},${isp_r}

!Initial cloud mixing ratio (kg/kg) and number mixing ratio (#/kg)
cloud_init=${cim3},${cim0},${cimx},${cimy}  !.001,100.e6

!Initial rain mixing ratio (kg/kg) and number mixing ratio (#/kg)
rain_init=${rim3},${rim0}

!Constant rain source mean-mass diameter (m) and number mixing ratio (#/kg)
rain_source=${rs_dm},${rs_N}

! number of moments for each species
!To run AMP as the bin scheme, set num_h_moments = 1 and num_h_bins = 33
!To run AMP as AMP, set num_h_moments = 2 or 3 and num_h_bins = 1
num_h_moments=${nhm}
num_h_bins=${nhb}

!AMP control - which moments to predict
imomc1 = ${imc1}  !1st predicted cloud moment
imomc2 = ${imc2}  !2nd predicted cloud moment (if 3M)
imomr1 = ${imr1}  !1st predicted rain moment
imomr2 = ${imr2}  !2nd predicted rain moment (if 3M)

!Microphysics process control
donucleation = ${l_nuc_cond_f}
docondensation = ${l_nuc_cond_f}
docollisions = ${l_coll_f}
dosedimentation = ${l_sed_f}
log_predictNc = .true.

! Aerosol initialization
num_aero_moments=1
num_aero_bins=1
aero_N_init=0. !or CCN at 1% SS
aero_sig_init=1.4
aero_rd_init=0.05e-6

! Background values for each moment (assumed the same for all species)
!This is the original KiD method for initializing hydrometeors. Not
!recommended (Adele)
mom_init=0,0,0

! param_val_fpath="../../CloudBOSS/boss_slc_param_values_${Pdep}30${imc1}${imc2}.csv"
param_val_fpath="$param_val_fpath"
! param_infl_sigma_fpath="../../CloudBOSS/boss_slc_param_sigma_${Pdep}30${imc1}${imc2}.csv"
!param_val_fpath_2cat="../../CloudBOSS/boss_2cat_param_values.csv"
param_val_fpath_2cat="../../BOSS-drizzLES/params/boss_post_mcmcNUTS0p8_les_obsσ_rpn_d_covobsrun_lwprr.csv"

iautoq = 1,
ivTnc = 2,
ivTqc = 2,
ivTnr = 2,
ivTqr = 2,
dNc_min = 1e-6,
dNc_max = 100e-6,
dNr_min = 10e-6,
dNr_max = 2e-3,
vTncmax = 0.5,
vTqcmax = 0.5,
vTnrmax = 10.,
vTqrmax = 10.,

n_cat = ${n_cat}
idraw=1

nmom_diag = 14

/

&case
icase=${caselist[ic]}
/

&control
mphys_scheme='${mp_id}'
dt=0.5            !Timestep length (s)
dgstart=0.0       !When to start diagnostic output
dg_dt=5.0         !Timestep for diagnostic output
wctrl(1)=0.      !Updraft speed
tctrl(1)=${t1}    !Total length of simulation (s)
tctrl(2)=${t2}     !May not be used, depends on the case. Typically the period of w oscillation
tctrl(3)=1080.    !For cases 105-107
tctrl(4)=1200.    !For cases 105-107
zctrl=${zc} !zctrl(1) is the domain height, (2) and (3) specify the location to init. hydromets.
/

&ppe
l_ppe=.true.
s_sample_dist="$s_sample_dist"
custom_dens_path="$custom_dens_path"
custom_bins_path="$custom_bins_path"
n_ppe=$6
irealz=$5
deflation_factor=1.
Dm_min=0.
Dm_max=0.
Na_min=${Na_min}e6
Na_max=${Na_max}e6
w_min=$w_min
w_max=$w_max
lhs_path="$lhs_path"
l_ppe_nevp=$l_ppe_nevp
l_ppe_condevp=$l_ppe_condevp
l_ppe_coal=$l_ppe_coal
l_ppe_sed=$l_ppe_sed
/

&switch
l_advect=${l_adv}
l_noadv_qv=${l_noadv_qv}
l_noadv_hydrometeors=${l_noadv_hyd}
l_noadv_theta=.true.
l_diverge=.false.
l_nodiv_hydrometeors=.true.
l_fix_theta=.true.
l_diverge_advection=.false.
l_fix_aerosols=.true.
l_periodic_bound=.false.
l_truncated=.false.
l_init_test=.false.
l_use_nn=${l_use_nn} ! whether use NN based AMP or old AMP algo
l_boss_partition_liq=.false.
l_boss_save_dsd=.false.
l_getrates=.false.
/

&addcontrol
KiD_outdir='$outdir'
ampORbin='${ampORbin:l}'
bintype='${bintype:l}'
mp_proc_dg=.true.
initprof='i' ! 'i' for an increasing initial water profile wrt height, 'c' for constant
l_hist_run=.false.
extralayer=.false.
!l_diag_nu=.false.
moments_diag = -6, -3, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15
/
END
     ${KiD_path}bin/KiD_1D.exe $nml_fn
  done
done
