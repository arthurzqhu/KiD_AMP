#!/bin/zsh

# config of the run
mconfig_temp='evaponly' # case/folder name. determined automatically if set empty
caselist=(101) #(101 102 103 105 106 107)
case_num=${#caselist[@]}
ampORbin=("BIN" "AMP")
bintype=("TAU" "SBM")
tests2run_num=$((${#ampORbin[@]}*${#bintype[@]}))

# initial condition for all cases
icimm=0. # initial cloud mass kg/kg
icinm=0. # initial cloud number 1/kg
irimm=0.
irinm=0.
rs_dm=0. # mean-mass diameter (m), ignores the case once this is non-zero
rs_N=0.  # number mixing ratio (#/kg)
isp_c=4  # shape parameter for cloud
isp_r=4  # shape parameter for rain
imc1=0 # II moment for cloud
imc2=6 # III moment for cloud
imr1=0 # II moment for rain
imr2=6 # III moment for rain
ztop=6000. # top of the domain
t1=100.
t2=5.
# switches
l_nuc_cond_s=1
l_coll_s=0
l_sed_s=0
l_adv_s=0


# set initial water if nucleation/condensation and/or adv is turned off 
#if [[ $l_nuc_cond_s -eq 0 || $l_adv_s -eq 0 ]]; then 
#   icimm=0.001     
#   icinm=100.e6  
#   irimm=0.5e-3
#   irinm=3.e3
#else 
#   icimm=0.
#   icinm=0.      
#fi

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


iw=2
ia=100

idm=$1
irh=$2

#icinm=100.e6
#icimm=$((($idm*1.e-6)**3*3.14159/6*1000.*$icinm))
icimm=0.001
icinm=$(($icimm/(($idm*1.e-6)**3*3.14159/6*1000.)))

#irinm=1.e4
#irimm=$((($idm*1.e-6)**3*3.14159/6*1000.*$irinm))

var1str=dm$1
var2str=rh$2

#idmr=400
#irinm=${inr}
#

#echo rainmass=$irimm
#echo rainnumber=$irinm
#
#idmc=20
#icinm=${inc}.e6
#icimm=$((($idmc*1.e-6)**3*3.14159/6*1000.*$icinm))
#
#echo cloudmass=$icimm
#echo cloudnumber=$icinm

mconfig=${mconfig_temp}

for ((iab=1; iab<=${#ampORbin[@]}; iab=iab+1))
do
  for ((ibt=1; ibt<=${#bintype[@]}; ibt=ibt+1))
  do
	echo "${ampORbin[$iab]}"-"${bintype[$ibt]}"
    if [[ ${ampORbin[$iab]} = 'AMP' ]]; then
      nhm='3,3'
      nhb='1,1'
      # changes nhm based on the input 
      if [ $imc1 = $imc2 ]; then
        nhm=${nhm//3,/$'2,'}
      fi
      if [ $imr1 = $imr2 ]; then
        nhm=${nhm//,3/$',2'}
      fi
      else
        if [[ ${bintype[$ibt]} = 'SBM' ]]; then
          nhm='1,1'
          nhb='33,33'
        else
          nhm='2,1'
          nhb='34,1'
        fi
      fi
      outdir=output/$(date +'%Y-%m-%d')/$mconfig/${ampORbin[$iab]}_${bintype[$ibt]}/${var1str}/${var2str}/
	  for ((ic=1; ic<=case_num; ic++))
	  do
	    if [[ ${caselist[ic]} -gt 104 ]] && [[ ${caselist[ic]} -lt 200 ]]
	    then
	      zc=0
        else
	      zc="$ztop,600.,1200."
        fi
	    if [ ! -d $outdir ]; then
	      mkdir -p $outdir
	    fi
	    echo "${caselist[ic]}"
	    cat > namelists/jobnml/${mconfig}_${ampORbin[$iab]}_${bintype[$ibt]}_${var1str}_${var2str}.nml << END
&mphys
! hydrometeor names
h_names='cloud','rain'

!Moment names
mom_names='M1','M2','M3'

!Initial shape parameter
h_shape=${isp_c},${isp_r}

!Initial cloud mixing ratio (kg/kg) and number mixing ratio (#/kg)
cloud_init=${icimm},${icinm}  !.001,100.e6

!Initial rain mixing ratio (kg/kg) and number mixing ratio (#/kg)
rain_init=${irimm},${irinm}

!Constant rain source mean-mass diameter (m) and number mixing ratio (#/kg)
rain_source=${rs_dm},${rs_N}

!Initial supersaturation ratio
!ss_init=0.1

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
donucleation = .false.
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
dt=0.05           !Timestep length (s)
dgstart=0.0       !When to start diagnostic output
dg_dt=0.1         !Timestep for diagnostic output
wctrl(1)=${iw}      !Updraft speed
tctrl(1)=${t1}    !Total length of simulation (s)
tctrl(2)=${t2}     !May not be used, depends on the case. Typically the period of w oscillation
tctrl(3)=1080.    !For cases 105-107
tctrl(4)=1200.    !For cases 105-107
zctrl=${zc} !zctrl(1) is the domain height, (2) and (3) specify the location to init. hydromets.
rhctrl=${irh}
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
/

&addcontrol
KiD_outdir='$outdir'
ampORbin='${ampORbin[$iab]:l}'
bintype='${bintype[$ibt]:l}'
mp_proc_dg=.true.
initprof='i' ! 'i' for a initial water profile increase wrt height, 'c' for constant
!l_diag_nu=.false.
/
END
     ./bin/KiD_1D.exe namelists/jobnml/${mconfig}_${ampORbin[$iab]}_${bintype[$ibt]}_${var1str}_${var2str}.nml
#            done
#          done
#        done
#      done
    done
  done
done
