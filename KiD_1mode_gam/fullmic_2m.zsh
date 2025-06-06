#!/bin/zsh

# config of the run
conf_basename="fullmic_2m_fullg" # case/folder name. determined automatically if set empty
caselist=(101) #(101 102 103 105 106 107)
case_num=${#caselist[@]}
ampORbin=("AMP")
bintype=("TAU")
tests2run_num=$((${#ampORbin[@]}*${#bintype[@]}))

# initial condition for all cases
icimm=0. # initial cloud mass kg/kg
icinm=0. # initial cloud number 1/kg
irimm=0.
irinm=0.
rs_dm=0. # mean-mass diameter (m), ignores the case once this is non-zero
rs_N=0. # number mixing ratio (#/kg)
imc1=0 # II moment for cloud
imc2=0 # III moment for cloud
imr1=0 # II moment for rain
imr2=0 # III moment for rain
ztop=6000. # top of the domain
zcb=600. # cloud base height
zct=1200. # cloud bottom height
t1=3600.
t2=900.

# switches for nucleation/condensation, collision, sedimentation, and advection
l_nuc_cond_s=1
l_coll_s=0
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

ia=$1
iw=$2
var1str=Na$ia
var2str=w$iw

# reset oscillation time based on updraft speed to prevent overshooting
if [[ $((($ztop-$zct)/$iw)) -lt $t2 && $l_adv_s -eq 1 ]]; then
  t2=$((($ztop-$zct)/$iw))
  t1=$(($t2*4))
fi

config_fname=${conf_basename}
for ((iab=1; iab<=${#ampORbin[@]}; iab=iab+1))
do
  for ((ibt=1; ibt<=${#bintype[@]}; ibt=ibt+1))
  do
	echo "${ampORbin[$iab]}"-"${bintype[$ibt]}"
   if [[ ${bintype[$ibt]} = 'SBM' ]]; then
      isp_c=4  # shape parameter for cloud
      isp_r=4  # shape parameter for rain
   else
      isp_c=10
      isp_r=10
   fi
    if [[ ${ampORbin[$iab]} = 'AMP' ]]; then
      nhm='2,2'
      nhb='1,1'
      else
        if [[ ${bintype[$ibt]} = 'SBM' ]]; then
          nhm='1,1'
          nhb='33,33'
        else
          nhm='2,1'
          nhb='34,1'
        fi
      fi
      outdir=output/$(date +'%Y-%m-%d')/$config_fname/${ampORbin[$iab]}_${bintype[$ibt]}/$var1str/$var2str/
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
	    echo "${caselist[ic]}"
	    cat > namelists/jobnml/${config_fname}_${ampORbin[$iab]}_${bintype[$ibt]}_${var1str}_${var2str}.nml << END
&mphys
! hydrometeor names
h_names='cloud','rain'

!Moment names
mom_names='M1','M2'

!Initial shape parameter
h_shape=${isp_c},${isp_r}

!Initial cloud mixing ratio (kg/kg) and number mixing ratio (#/kg)
cloud_init=${icimm},${icinm}  !.001,100.e6

!Initial rain mixing ratio (kg/kg) and number mixing ratio (#/kg)
rain_init=${irimm},${irinm}

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
donucleation = .true.
docondensation = ${l_nuc_cond_f}
docollisions = ${l_coll_f}
dosedimentation = ${l_sed_f}

! Aerosol initialization
num_aero_moments=1
num_aero_bins=1
aero_N_init=${ia}.e6 !or CCN at 1% SS
aero_sig_init=1.4
aero_rd_init=0.05e-6

! Background values for each moment (assumed the same for all species)
!This is the original KiD method for initializing hydrometeors. Not
!recommended (Adele)
mom_init=0,0,0
/

&case
icase=${caselist[ic]}
/

&control
mphys_scheme='amp'
dt=0.5            !Timestep length (s)
dgstart=0.0       !When to start diagnostic output
dg_dt=5.0         !Timestep for diagnostic output
wctrl(1)=${iw}      !Updraft speed
tctrl(1)=${t1}    !Total length of simulation (s)
tctrl(2)=${t2}     !May not be used, depends on the case. Typically the period of w oscillation
tctrl(3)=1080.    !For cases 105-107
tctrl(4)=1200.    !For cases 105-107
zctrl=${zc} !zctrl(1) is the domain height, (2) and (3) specify the location to init. hydromets.
!rhctrl=${rh}
!pctrl_v=${pcpt}
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
l_periodic_bound=.False.
l_truncated=.true.
l_init_test=.false.
/

&addcontrol
KiD_outdir='$outdir'
ampORbin='${ampORbin[$iab]:l}'
bintype='${bintype[$ibt]:l}'
mp_proc_dg=.true.
initprof='i' ! 'i' for an increasing initial water profile wrt height, 'c' for constant
extralayer=.false.
l_hist_run=.false.
!l_diag_nu=.false.
/
END
     ./bin/KiD_1D.exe namelists/jobnml/${config_fname}_${ampORbin[$iab]}_${bintype[$ibt]}_${var1str}_${var2str}.nml
    done
  done
done
