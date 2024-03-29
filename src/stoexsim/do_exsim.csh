#!/bin/tcsh

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

### the base implementation of EXSIM ###

rm -rf output
# Setup some variables
set pi=3.14159265
set dip=85	# The dip of the fault (faults are dipping 60 to 90 degrees)
set Zh=$argv[2]	# Depth to hypocentre (this will be 3 km)
set HypW=0 # Down-dip distance for Hypocentre (always start at top = 0)
set Ztop_min=1 # The top of the seismogenic zone (the reservoir top)
set SeisDepth=30 # The seimogenic depth (pers. comm. Dost: 10 - 13 km)
set mech=$argv[1]	# Mechanism (used to calculate fault geometry)
set RANDOM=`bash -c 'echo $RANDOM'`	# select a psuedo-random number
set sigma_L=0.15 # sigma of the fault length (ln-unit)
set sigma_W=`echo -$sigma_L` # sigma of the fault width (ln-unit) - note epsilon for sigma_L and sigma_W is assumed prefectly anti-correlated (to maintain fault area)
set niter=1 # number of simulations to perform (no. hypocentres)
set iter=1
# Stress-parameter model: Define hinge, above which SD is const, below hinge 
# SD decreases linearly in log-space to sd_low at sd_low_M, below which is constant.
set sd_M=5.0	# M hinge point beyond which stress-parameter is constant
set sd_scale=1.000	# A scaling factor for the stress-parameter
set sd_low=27	# Stress-parameter at M=sd_low_M
set sd_low_M=3 # M for sd_low
set sd_in=50	# Stress-parameter above M=sd_M

# min max and mo. response periods
set fmin=0.01
set fmax=100
set fn=99

# Calculate the required M-scaling between sd_M_low and sd_M
set sdi_in=`echo "" | awk '{print (((log('$sd_in')/log(10.))-(log('$sd_low')/log(10.)))/(('$sd_M'-'$sd_low_M')))}'`
while ( $iter <= $niter ) 
set DYN_RAN=`bash -c 'echo $RANDOM'` 
set epsilon=`./random/normdev 0 1 $RANDOM $DYN_RAN -5 5` #normal distibution
#set epsilon=0
# We will work and store the results in this directory:
mkdir output 
cd output 

foreach Mw ( $argv[3] )
set DYN_RAN2_Td575_interEV=`bash -c 'echo $RANDOM'` # inter-event epsilon (used for Td)
echo $epsilon >> LW_epsilon.dat
# Set the SD (stress-parameter) at the given M
set sd=`echo "" | awk '{print '$sd_scale'*'$sd_low'*10^('$sdi_in'*('$Mw'-'$sd_low_M'))}'`
# If above sd_M hinge
if ( `echo "" | awk '{if ('$Mw'>'$sd_M'){print 1}; if ('$Mw'<='$sd_M'){print 0}}'` == "1" ) then
	set sd=`echo "" | awk '{print '$sd_in'}'`
endif
# If below sd_M_low hinge
if ( `echo "" | awk '{if ('$Mw'<'$sd_low_M'){print 1}; if ('$Mw'>='$sd_low_M'){print 0}}'` == "1" ) then
	set sd=`echo "" | awk '{print '$sd_low'}'`
endif

# Now setup the fault dimensions 
# M < 5.25: square - area equiv to circular rupture model (Eshelby)
set HypL=0 # Along-strike distance for Hypocentre (for small M, point source = 0)
if ( `echo "" | awk '{if('$Mw'<5.25){print 1}; if('$Mw'>=5.25){print 0}; }'` == 1 ) then
	set W=`echo "" | awk '{print 2*((7/16)*10^(1.5*'$Mw'+9.1)/('$sd'*1e5))^(1/3)/1000}'`
	set W=`echo "" | awk '{print sqrt('$pi'*('$W'/2)^2)}'`
	# Adjustment to W&C94 stress reference - assume ref (25 bar: SS; 15 bar: N,R) for W&C94
	if ( "$mech" == "S" ) then
		set stress_factor=`echo ""|awk '{print ('$sd'/25)^(1.0/3.0)}'`
	else
		set stress_factor=`echo ""|awk '{print ('$sd'/15)^(1.0/3.0)}'`
	endif
	set W=`echo "" | awk '{print '$W'*'$stress_factor'}'`
	set L=`echo "" | awk '{print '$W'}'`
	set FL=`echo "" | awk '{print '$L'}'`
	set FW=`echo "" | awk '{print '$W'}'`
	set Wx=`echo "" | awk '{print '$W'*cos('$dip'*'$pi'/180)}'`
	set Wy=`echo "" | awk '{print '$W'*sin('$dip'*'$pi'/180)}'`
    echo "Fault Dims (mean): W = "$FW"; L = "$FL"; Wx = "$Wx"; Wy = "$Wy" (dip = "$dip")"
	# Add epsilon to fault dimensions
	set W=`echo ""| awk '{print '$W'*10^('$sigma_W'*'$epsilon')}'`
	set Wx=`echo ""| awk '{print '$Wx'*10^('$sigma_W'*'$epsilon')}'`
	set Wy=`echo ""| awk '{print '$Wy'*10^('$sigma_W'*'$epsilon')}'`
	set FW=`echo ""| awk '{print '$FW'*10^('$sigma_W'*'$epsilon')}'`
	set L=`echo ""| awk '{print '$L'*10^('$sigma_L'*'$epsilon')}'`
	set FL=`echo ""| awk '{print '$FL'*10^('$sigma_L'*'$epsilon')}'`
    echo "Fault Dims (+/- sigma): W = "$FW"; L = "$FL"; Wx = "$Wx"; Wy = "$Wy" (dip = "$dip")"
	set Ztop=`echo "" | awk '{print '$Zh'-'$Wy'/2}'` # Assume Zh at centre of fault
	set Ztop=`echo "" | awk '{if('$Ztop' > '$Ztop_min'){ print '$Ztop'};if('$Ztop'<= '$Ztop_min'){print '$Ztop_min' }}'`
endif

# for M >= 5.25: use Wells and Coppersmith '94)
# Distinguish S: strike-slip and N: normal
if ( `echo "" | awk '{if('$Mw'<5.25){print 1}; if('$Mw'>=5.25){print 0}; }'` == 0 ) then
	if ( "$mech" == "U" ) then
		set W=`echo "" | awk '{print 10^(-1.01+0.32*'$Mw')}'`
		set L=`echo "" | awk '{print 10^(-2.44+0.59*'$Mw')}'`
	else if ( "$mech" == "S" ) then
		set W=`echo "" | awk '{print 10^(-0.76+0.27*'$Mw')}'`
		set L=`echo "" | awk '{print 10^(-2.57+0.62*'$Mw')}'`
	else if ( "$mech" == "N" ) then
		set W=`echo "" | awk '{print 10^(-1.14+0.35*'$Mw')}'`
		set L=`echo "" | awk '{print 10^(-1.88+0.50*'$Mw')}'`
	else	
		echo "Unsupported mechanism"
		exit
	endif
	# Adjustment to W&C for different stress - assume ref 30 bar (central model)
	#set stress_factor=`echo ""|awk '{print (30/'$sd')^(1.0/3.0)}'`
	set stress_factor=1 # don't use this
	# Assuming since W&C is empirical - don't want to adjust: assume instead tha the L/W stress parameter relation needs adjustment (no physical constraint)
	set W=`echo ""| awk '{print '$W'*'$stress_factor'}'`
	set L=`echo ""| awk '{print '$L'*'$stress_factor'}'`

	# Caluculate Wx and Wy 	
	set Wx=`echo "" | awk '{print '$W'*cos('$dip'*'$pi'/180)}'`
	set Wy=`echo "" | awk '{print '$W'*sin('$dip'*'$pi'/180)}'`
	# This will force W&C'94
	#set FL=`echo "" | awk '{print 0}'`
	#set FW=`echo "" | awk '{print 0}'`
	set FL=$L
	set FW=$W
	# Add epsilon to fault dimensions
    echo "Fault Dims (mean): W = "$FW"; L = "$FL"; Wx = "$Wx"; Wy = "$Wy" (dip = "$dip")"
	set W=`echo ""| awk '{print '$W'*10^('$sigma_W'*'$epsilon')}'`
	set Wx=`echo ""| awk '{print '$Wx'*10^('$sigma_W'*'$epsilon')}'`
	set Wy=`echo ""| awk '{print '$Wy'*10^('$sigma_W'*'$epsilon')}'`
	set FW=`echo ""| awk '{print '$FW'*10^('$sigma_W'*'$epsilon')}'`
	set L=`echo ""| awk '{print '$L'*10^('$sigma_L'*'$epsilon')}'`
	set FL=`echo ""| awk '{print '$FL'*10^('$sigma_L'*'$epsilon')}'`
	set Ztop=`echo "" | awk '{print '$Zh'-'$Wy'/2}'`
	set Ztop=`echo "" | awk '{if('$Ztop' > '$Ztop_min'){ print '$Ztop'};if('$Ztop'<= '$Ztop_min'){print '$Ztop_min' }}'`
    echo "Fault Dims (+/- sigma): W = "$FW"; L = "$FL"; Wx = "$Wx"; Wy = "$Wy" (dip = "$dip")"
	if ( `echo "" | awk '{if(3+'$Wy'>'$SeisDepth'){print 1}else{print 0} }'` == 1 )  then
        # If so set to max width (SeisDepth-3km) and adjust FL to maintain fault area

        set Wy=`echo "" | awk '{print '$SeisDepth'-'$Ztop'}'`
        set FArea=`echo "" | awk '{print '$FL'*'$FW'}'`
        set FW=`echo "" | awk '{print '$Wy'/sin('$dip'*'$pi'/180)}'`
        set FL=`echo "" | awk '{print '$FArea'/'$FW'}'`
        echo "Fault exceeds seismogenic depth."
        echo "Revised fault Dims: W = "$FW"; L = "$FL"; Wx = "$Wx"; Wy = "$Wy" (dip = "$dip")"
    endif
	set HypL=`echo "" | awk '{srand('$DYN_RAN'); print '$FL'*rand()}'` # Along-strike distance for Hypocentre (random)
#	set HypL=`echo "" | awk '{srand('$DYN_RAN'); print '$FL'*0.5}'` # Along-strike distance for Hypocentre (middle)

endif

rm -f  Rjb.dat
foreach Rjb( \
$argv[4] \
)
	
	echo $Rjb "-90" >> Rjb.dat
	set Rhyp_e=`echo | awk '{print sqrt('$Rjb'^2+'$Zh'^2)}'` # Assume here Rjb=Repi
end


# Update the input file
cat ../exsim_dmb.params_base > exsim_dmb.params.tmp 

# Push the variables int othe template param file
sed s/"#MW#"/"$Mw"/g ./exsim_dmb.params.tmp | sed s/"#HypL#"/"$HypL"/g | sed s/"#HypW#"/"$HypW"/g | sed s/"#ZTOP#"/"$Ztop"/g | sed s/"#MECH#"/"$mech"/g |  sed s/"#DIP#"/"$dip"/g | sed s/"#FL#"/"$FL"/g |sed s/"#FW#"/"$FW"/g  | sed s/"#Fmin#"/"$fmin"/g | sed s/"#Fmax#"/"$fmax"/g | sed s/"#Fn#"/"$fn"/g | sed s/"#SD#"/"$sd"/g |sed s/"#SEED#"/"$DYN_RAN"/g  > exsim_dmb.params 

# get the template site and crustal amp files
cp ../site_amps.txt ../crustal_amps.txt .

cat Rjb.dat >> exsim_dmb.params
echo "stop" >> exsim_dmb.params

# Print out the inputs
echo "M = "$Mw
echo "Iteration = "$iter
echo "stress drop = "$sd
echo "Ztop = "$Ztop
echo "L,W Epsilon = "$epsilon

# Run the code (Dave Boores Version of EXSIM)
echo "" | ../exsim_dmb_171016/exsim_dmb_be 

end
#compress
foreach file ( *.out )
	gzip $file 
end
cd ../
@ iter+=1
end
