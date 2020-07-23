#!/bin/bash

caselist=(101 102 103 105 106 107)
case_num=${#caselist[@]}


icimm=0.0
icinm=0.0
isp_c=4  # shape parameter for cloud
isp_r=4  # shape parameter for rain
imc1=0 # II moment for cloud
imc2=6 # III moment for cloud
imr1=0 # II moment for rain
imr2=6 # III moment for rain
#ia=50

#for icimm in 0.0 0.0003 0.001 0.003 #initial cloud mass mixing ratio
#do
	#for icinm in 0 30 100 300 #initial number mixing ratio
	#do
	
#for isp_c in 2 15
#do
#	for isp_r in 2 15
#	do

#mc1=(0 2 2 4 4 6)
#mc2=(2 0 4 2 6 4)
#mnum=${#mc1[@]}

#for ((imnum=0; imnum<mnum; imnum++))
#do
#	imc1=${mc1[imnum]}
#	imc2=${mc2[imnum]}
#for ((imc1=0; imc1<8; imc1=imc1+2))
#do
#	for ((imc2=imc1+2; imc2<=8; imc2=imc2+2))
#	do
for ((ia=50; ia<100; ia=ia+100))
do
echo $ia
		echo $imc1 $imc2
		outdir=output/SBM/$(date +'%Y-%m-%d')/a${ia}/
		for ((ic=0; ic<case_num; ic++))
		do
			if [ ${caselist[ic]} -gt 104 ] && [ ${caselist[ic]} -lt 200 ]
			then
				zc=0
			else
				zc="3000.,600.,900."
			fi
			if [ ! -d $outdir ]; then
			    mkdir -p $outdir
			fi
			echo “${caselist[ic]}”
			cat > namelists/BIA.nml << END

&mphys
! hydrometeor names
h_names='cloud','rain'

!Moment names
mom_names='M1','M2','M3'

!Initial shape parameter
h_shape=$isp_c,$isp_r

!Initial cloud mixing ratio (kg/kg) and number mixing ratio (#/kg)
cloud_init=$icimm,$icinm  !.001,100.e6

!Initial rain mixing ratio (kg/kg) and number mixing ratio (#/kg)
rain_init=0.0,0.0

! number of moments for each species
!To run AMP as the bin scheme, set num_h_moments = 1 and num_h_bins = 33
!To run AMP as AMP, set num_h_moments = 2 or 3 and num_h_bins = 1
num_h_moments=1,1
num_h_bins=33,33

!AMP control - which moments to predict
imomc1 = $imc1  !1st predicted cloud moment
imomc2 = $imc2  !2nd predicted cloud moment (if 3M)
imomr1 = $imr1  !1st predicted rain moment
imomr2 = $imr2  !2nd predicted rain moment (if 3M)

!Microphysics process control
donucleation = .true.
docondensation = .true.
docollisions = .true.
dosedimentation = .true.

! Aerosol initialization
num_aero_moments=0
num_aero_bins=1
aero_N_init=$ia.e6 !or CCN at 1% SS
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
dt=1.0            !Timestep length (s)
dgstart=0.0       !When to start diagnostic output
dg_dt=1.0         !Timestep for diagnostic output
wctrl(1)=2.0      !Updraft speed
tctrl(1)=3600.    !Total length of simulation (s)
tctrl(2)=600.     !May not be used, depends on the case. Typically the period of w oscillation
tctrl(3)=1080.	  !For cases 105-107
tctrl(4)=1200.    !For cases 105-107
zctrl=$zc !zctrl(1) is the domain height, (2) and (3) specify the location to init. hydromets.
/

&switch
l_advect=.true.
l_noadv_qv=.false.
l_noadv_hydrometeors=.false.
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
ampORbin='bin'
bintype='sbm'
/
END
		./bin/KiD_1D.exe namelists/BIA.nml
				#done
			#done
#		done
	done
done
