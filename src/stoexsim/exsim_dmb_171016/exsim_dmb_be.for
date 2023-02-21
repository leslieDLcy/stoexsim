      program EXSIM_DMB  
! Changes to this code by Ben Edwards indicated comments preceded by BE
! In addition all fault dimension arrays (200,200) have been extended to 1000,1000
      
! This is a modification of Motazedian and Atkinson's EXSIM.  See 
!    Boore, D. M. (2009). Comparing stochastic point-source and finite-source 
!    ground-motion simulations: SMSIM and EXSIM, Bull. Seismol. Soc. Am. 99, 3202-3216.
!    for details regarding the modifications.


!
!     EXSIM_DMB is designed to work for a list of observation sites for a
!     user-defined fault geometry and location. For each site the output consists of:
!       For the first of nsims simulations and the first of n_hypocenters hypocenters:
!         acceleration time series
!         <more files>
!       The average Fourier spectra and PSA over all of the simulations and hypocenters.

!^^^^^^^ FOR FUTURE REVISIONS ^^^^^^^^^^^^^^^^^^^^
!
! Allow an option not to write to screen (Andreas Skarlatoudis found that this
! will speed up the execution of the program by more than a factor of 10, although some
! trial runs that I did found this not to be true).
!
! Use more meangingful variable and subroutine names: e.g., subroutine sum4avg, 
! averagePSA4hypo rather than averagePSA
!
!^^^^^^^ FOR FUTURE REVISIONS ^^^^^^^^^^^^^^^^^^^^


! !!! NOTE: If the PSA and/or acceleration time series of each realization is
! desired, suitable modifications need to be made to the program

!Here is a sample parameter file for the program:

!!Parameter file for program exsim_dmb
!! Revision of program involving a change in the control file on this date:
!   11/22/11
!!
!!-------------------------------------------------------------------------
!! ******* Input parameters common to SMSIM and EXSIM (in the order in which 
!!         they appear in the SMSIM parameter file) *******
!!-------------------------------------------------------------------------
!!
!!Title
!  Runs for comparing EXSIM and SMSIM 
!!rho, beta, prtitn, radpat, fs:
!    2.8 3.7 0.707 0.55 2.0
!!gsprd: r_ref, nsegs, (rlow(i), a_s, b_s, m_s(i))  (Usually set
!!r_ref = 1.0 km)
!    1.0
!    3
!      1.0 -1.3 0.0 6.5
!     70.0 +0.2 0.0 6.5
!    140.0 -0.5 0.0 6.5
!!q: fr1, Qr1, s1, ft1, ft2, fr2, qr2, s2, c_q
!   1.0 1000.0 0.0 1.4242 1.4242 1.0 893.0 0.32 3.7   
!!path duration (ndur_hinges, 
!! (rdur(i), dur(i), i = 1, ndur_hinges), durslope)
!    4
!    0.0 0.0
!   10.0 0.0
!   70.0 9.6
!  130.0 7.8
!  0.04
!!site diminution parameters: fmax, kappa_0, dkappadmag, amagkref
!! (NOTE: fmax=0.0 or kappa_0=0.0 => fmax or kappa are not used.  I included this
!!  to prevent the inadvertent use of both fmax and kappa to control the diminution
!!  of high-frequency motion (it would be very unusual to use both parameters
!!  together.  Also note that if do not want to use kappa, dkappadmag must also
!!  be set to 0.0).
!    0.0 0.005 0.0 6.0
!!low-cut filter corner, nslope (0 ==> no filter)
! 0.05 8
!!window params: iwind(0=box,1=exp), taper(<1), eps_w, eta_w, f_tb2te, f_te_xtnd
!! (see SMSIM manual for the meaning of the parameters)
!! As of 11/25/11, I will not use the shape parameters, using the default
!! parameters in the call to wind2 instead.  The only parameters I use as of this
!! date are iwind and taper.
!! BUT: placeholders must be included for eps_w, eta_w, f_tb2te, f_te_xtnd, because some day
!! they may be used.
!    1 0.05 0.2 0.05 2.0 2.0
!!timing stuff: dur_fctr, dt, tshift, seed, nsims, iran_type (0=normal;1=uniform)
!! NOTE: these are the SMSIM parameters, but for now (11/25/11) I will read and use
!! the current EXSIM parameters, as given in the next uncommented line.
!! The reason not to change to the SMSIM parameters is that I do not have the time to
!! make sure that the program is revised correctly.   The tpadl and tpadt params do not
!! automatically account for magnitude, and if the values are not changed each time amag
!! is changed they may be unnecessarily long for small events and too short for large 
!! events.  It is up to the user to specify adequate values (these values are adjusted
!! to appropriate sizes automatically in SMSIM).
!!    1.3  0.002 60.0 123.0 100 0
!!tpadl, tpadt, dt, seed, nsims
! 50.0 20.0 0.002 123.0 10               
!!
!!-------------------------------------------------------------------------
!! ******* Input parameters specific to EXSIM *******
!!-------------------------------------------------------------------------
!!
!! SOURCE PARAMETERS:
!!
!!MW, Stress
!  7.0 140.0     
!!lat and lon of upper edge of fault
!  0.0 0.0  
!!strike,dip, depth of fault
!  0.0 50.0 2.0             
!!fault type (S=strikeslip; R=reverse; N=normal; U=undifferentiated) 
!! (Only used if Wells and Coppersmith is used to obtain FL and FW).
!  R                               
!!fault length and width, dl, dw, stress_ref
!!Note: Force program to use Wells and Coppersmith (WC) for FL and/or FW if
!! either entry = 0.0.
!! If Wells and Coppersmith are used to obtain FL and/or FW, the WC values are
!! modified to account for the scaling implied by differences in the stress
!! specified above and a stress that is assumed to be valid for the generic WC
!! relations; this stress is stress_ref. The value of 70 bars is an educated
!! guess for stress_ref, but it is not based on a quantitative analysis.
!! The WC values of FL and/or FW are multiplied by the factor
!! (stress_ref/stress)^(1/3).
!! Note that four entries on the following line are needed as placeholders,
!! even if not used)
!  0  0  1.5 1.5 70.0 !fault length and width, dl, dw, stress_ref
!!vrup/beta 
!  0.8              
!!hypo location in along fault and down dip distance from the fault 
!!reference point (an upper corner)(-1.0, -1.0 for a random location); 
!!number of iterations over hypocenter (need an entry, but only used if 
!!either of the first two values are -1.0, indicating a random location)
! -1.0 -1.0 2                           
!!Enter type of risetime (1=original, 2=1/f0)
! 2
!!DynamicFlag (0=no), PulsingPercent
!  1   50.0                   
!!iflagscalefactor (1=vel^2; 2=acc^2; 3=asymptotic acc^2 (dmb))
!  2                               
!!islipweight = -1  -> unity slip for all subfaults,
!!islipweight =  0  -> specify slips read from text file, 
!!islipweight =  1  -> random weights
!  -1                        
!! Text file containing matrix of slip weights.  The weights are a matrix of 
!! values, with rows corresponding to the nl subfaults along the length, 
!! and the columns corresponding to the nw subfaults down the width of 
!! the fault. The values entered are normalized by the sum of all of the 
!! values before computing the moment for each subfault, so it is only
!! the relative size of the slip weights that matter (need a placeholder
!! even if do not assign the slip weights
!  slip_weights.txt
!!!deterministic flag,gama,nu,t0, impulse peak
!  0   1.0  90.0  4.0  10.  
!!
!!-------------------------------------------------------------------------
!! PARAMETERS RELATED TO PATH AND SITE:
!!-------------------------------------------------------------------------
!!
!!Name of crustal amplification file:
!  crustal_amps_sample.txt
!!Name of site amplification file:
!  site_amps_sample.txt
!!
!!-------------------------------------------------------------------------
!! PARAMETERS RELATED TO COMPUTATIONS OF AVERAGES:
!!-------------------------------------------------------------------------
!!
!!iflagfas_avg (1=arithmetic, 2=geometric, 3=rms: USE 3!)
!  3                               
!!iflagpsa_avg_over_sims (1=arithmetic: USE 1!, 2=geometric, 3=rms)
!! NOTE on 22 November 2011.  I used to advise using the geometric mean, but in
!! the course of working on a paper with Eric Thompson on RV calculations in SMSIM
!! I found that my TD calculations have used geoemtric averages until 03 August 1994,
!! when I switched to arithmetic averages, apparently as a result of a recommendation
!! by Bill Joyner.
!  1                               
!!iflagpsa_avg_over_hypos (1=arithmetic, 2=geometric, 3=rms)
!! The program first computes the average ground-motion intensity measure over the number
!! of simulations for a given hypocenter, and then computes an average of these over the
!! hypocenters.  There might some justification to use a geometric mean for this.
!  2                               
!!
!!-------------------------------------------------------------------------
!! PARAMETERS RELATED TO THE OUTPUT:
!!-------------------------------------------------------------------------
!!
!!Write acc, psa, husid files for each site?
! Y
!!Output file names stem:
!  M7.0_dl_dw_1.5_140b_2hyp_5trials_pulse50_1_f0_rt
!! %damping of response spectra
! 5.0
!!# of f and Min and Max F for response spectra
!  100 0.1   99.                
!!no. of frequencies for summary output (10 max):
! 4 
!!frequency (-1.0, 99.0 for pgv, pga):
! -1.0 99.0 0.5 5.0
!!
!!-------------------------------------------------------------------------
!! PARAMETERS RELATED TO THE SITES AT WHICH MOTIONS ARE COMPUTED:
!! Put this last for convenience in editing the params file.  For example, the site list 
!! can be very long, but by inserting "stop" in the list it is easy to select a small subset 
!! of site at which motions will be computed.
!!-------------------------------------------------------------------------
!!
!!Site coord flag (1=lat,long; 2=R,Az; 3=N,E)
!  2                       
!!If "Y" below and strike = 0.0:
!!  if site coord flag = 2, move origin of the radial line to the midpoint of
!!                         the top edge of the fault
!!  if site coord flag = 3 and siteLocation(1) = 0, redefine 
!!                         siteLocation(1) = 0 to be the midpoint of the 
!!                         top edge of the fault (so that the sites will be
!!                         along a line normal to the midpoint)
!!  if site coord flag = 3 and siteLocation(2) = 0, redefine
!!                         siteLocation(1) = 0 to be the far end of the fault,
!!                         so that the sites are along a line along the
!!                         strike of the fault 
! Y
!!Coordinates of each site (siteLocation(1), siteLocation(2)):
! 1              45
! 100            45
! stop
! 1000           45
! 1              90
! 100            90
! 1000           90
   
 
! Dates: 08/17/08 - Modifications by D. M. Boore, based on EXSIM release version 1.0, 
!                   October 10, 2005 (for more details regarding EXSIM, see 
!                   Motazedian,Dariush and Gail M. Atkinson, (2005)
!                   "Stochastic Finite-Fault Modeling Based on a Dynamic Corner Frequency",
!                   Bulletin of the Seismological Society of America,95,995-1010.
!        08/17/08 - Corrected typo "ISlipWeigth" and "falg", and allow all
!                   unit slip weights if flag -1.
!        08/18/08 - Read in y (=vrup/beta), write R in PSA_FA output  
!                   Revise column headers of PSA_FA output.
!                   Lots of changes to output, including writing
!                   separate PSA_FA files for each site (to make it
!                   easier to import into a graphics program).
!                   Allow the site coordinates to be entered as 
!                        isitecoordflag   coords
!                                     1  lat,long 
!                                     2  R, Az
!                                     3  N, E
!        08/19/08 - Reverse order of fault input, and allow use of Wells
!                   and Coppersmith if FL=0 or FW=0, put magnitude input before
!                   FL, FW (because the later needs amag).  Also get fault
!                   type (needed for Wells and Coppermsith).
!        08/20/08 - Change averaging algorithm (log avg for PSA, E(FAS^2)^1/2 for FAS)
!                   Return E(FAS^2)^1/2, not log of the value, out of the 
!                   binning routine SAMPLE, called by compute FACCN.
!
!                   Add input variable iflagscalefactor to determine whether to use
!                   scale factor based on the following:

!                   iflagscalefactor   integrand of integral
!                           1            velocity spectrum squared
!                           2            acceleration spectrum squared
!                           3            D. Boore's factor based on acc. sq. high-frequency asymptotes
 
!                   Add input variable to determine if use arithmetic or 
!                   geometric average for PSA
!        08/21/08 - Added scalefactor (h) to computeStochasticWave output and print out 
!                   the scalefactor and corner frequency for each subfault for 
!                   the first simulation.  Also include an option for tapering the 
!                   scaling factor.
!        08/22/08 - Write out the time series for the first realization only.
!        08/25/08 - Read in qmin
!                   In checkFFTpoints, increased upper index so that 2**n=32768.  
!                   Also, doubled dimensions of wave, totalWave, and subwave.
!                   In many subroutines replaced dimensions of arrays passed 
!                   in through the parameter list with "(*)".
!        08/27/08 - Added dur_sub to common /sub/, in order to print it out at some point.
!        11/06/08 - Add option to compute RMS of PSA, and also use same options for computing the average of FAS 
!                   (for study purposes).
!        11/07/08 - Changed way that averages are computed.  Now the accumulated sum is computed
!                   in computeAverageResponse and computeAverageFA and this is used in the main program
!                   to compute an average.  But a modification was needed to allow
!                   for negative values of log10(y) for geometric means.  I set the average to -9.99
!                   if accum_fas = 0.0 or accum_psa = 0.0.
!        11/08/08 - I may be having problems with the rms average of FAS, so compute this differently now.
!        11/15/08 - Increase maximum number of nl, nw to 200.
!        11/25/08 - I am concerned about the calls to ran1.   Note that the distributed version of
!                   exsim uses urand, not ran1.  Apparently I substituted ran1 for urand, but I think
!                   I made a mistake.  If used to
!                   generate a series of "random" numbers, ran1 should be called the
!                   first time with a negative integer (idum) as a seed.  
!                   Ran1 performs some initializations, including assigning
!                   values to Saved variables iv and iy, as well as changing
!                   idum.  In exsim ran1 is called for a number of things,
!                   possibly including locateHypocentreRandomly, 
!                   createRandomWeight, computeStochasticWave, and
!                   computing a random delay for the time series from 
!                   each subsource.   I would think that the idum, iv, 
!                   and iy values would need to be in sync with one another,
!                   but this will not be the case if ran1 is called with
!                   a different value for idum that is not the same as in
!                   previous calls.   In the previous version of exsim, iseed
!                   is used as idum in the calls in all calls except for
!                   the call for computing the random delay.  But iseed was 
!                   not constrained to be less than 0.0 initially, so the ran1
!                   initialization was probably not done.  And also there was
!                   one call to ran1 that used islipweight for idum, which
!                   could have values less than, equal, to, or greater than 0.0
!                   (in my runs it was less than 0.0).  When less than 0.0, as
!                   it was in all of my runs, the iv and iy would have been
!                   initialzed and would have been used in the next call to
!                   ran1, but the idum would not have been the reset value
!                   of islipweight.  I do not know the consequences of all of
!                   this.  I have changed exsim to contrain iseed to be less 
!                   than 0 initially and have used iseed in all calls to ran1.
!        11/27/08 - Added duration calculation
!        11/28/08 - Use 1/f0 for risetime
!        11/28/08 - Change input file to have a comment line before the name of the scale
!                   factor file.
!        11/30/08 - Move computation of random slipweights 
!                   if Islipweight == 1.0 from main program to getInputParameters
!                   (a more logical place because the slipweights array is filled
!                   in the same place for all values of ISlipWeight).
!        11/30/08 - Move "findSubFaultMoment" out of loop over isite
!        12/01/08 - Use "if then" in loop over subsources in 
!                   finding NumberofActiveSubs
!        12/01/08 - Remove "NumberofActivesSubs" from the arguments passed
!                   to subroutines computeStochasticWave and sourceSpectra, as
!                   the variable is not used.
!        12/01/08 - Added writing of istart, istop, nptsTotalWave, maxPointsOfWave
!                   in WritePar (I'm trying to figure out if findindex is needed and if
!                   code setting istop is correct).
!        11/28/08 - Use 1/firstElementF0 for risetime
!        12/01/08 - Use original exsim risetime
!        12/01/08 - Add option to choose one of three ways of obtaining risetime
!        12/06/08 - Allow iteration of random hypocenters for each station
!        12/08/08 - Read initial seed.
!        12/09/08 - Add scaling_coefficient, following Dariush's suggestion.
!        12/09/08 - Increase length of output file names.
!        12/17/08 - Move computations that do not depend on the loop over
!                   nsims out of the loop.  This required defining f0, etc as
!                   arrays (one value per subsource).  This is computationally efficient.
!                   Another advantage is that
!                   I can calculate the maximum duration for the total wave before doing
!                   the subsource wave calculation, and therefore I can use dynamic
!                   allocation for the time series array dimension.  This required
!                   changing common sub not to include arrays of variables that 
!                   depend on the subsource, passing these things through argument lists.
!                   I renamed common /sub/ to common /params/.  I also removed dur_sub
!                   from this common block.
!                   I also compute dur_sub outside of computeStochasticWave and pass the
!                   value into that subroutine through the argument list.
!        12/17/08 - Major revision involving the time series (use tpads and front and back, do not truncate
!                   the motions from each subsource; use dynamic allocations). NOTE: 
!                   Using Mavroeidis and Papageorgiou has not been tested, and it will probably not work 
!                   because I have not taken into account the added npadl.  I comment out a section of code
!                   before the call to M&P that may have to do with determining the index for which the pulse
!                   starts.  It would be easy to fix this later using t_arrive(i,j).
!        12/22/08 - If nl = nw = 1, over-ride the computation of a random hypocenter and 
!                   set n_hypocenters, i0, j0 = 1, 1, 1
!        12/31/08 - Write sqrt(nl*nw)*fas to help in deciding if that factor can correct for
!                   the incoherent summation.
!        12/31/08 - replace "np" for the taper correction at low frequencies by "ftaper_f0", 
!                   the ratio of ftaper and f0 read from the input file.  Also, there is no
!                   need for "taper_scalefactor", because ftaper_f0 == 0.0 implies
!                   no tapering toward low frequencies.
!        01/04/09 - Apply a new taper that attempts to correct for the incoherent summation,
!                   the decay of the subsource spectra for frequencies less than the subsource
!                   corner frequency (following Frankel, 1995), and the improper long-period spectral
!                   level.
!        01/06/09 - Remove "scaling_coefficient" from input and program.  
!                 - Put computation of c4taper, fc4taper outside of the sourcespectra function.
!                 - Revise format of input, using a header line before each line containing input.
!                 - Include a low-cut filter.
!                 - Compute average pga and pgv.
!        01/10/09 - Increase max number of frequencies to 500, without changing anything else.
!        01/13/09 - "apply_taper" was read in from the parameter file, but it has not been used since
!                   I added the Frankel-like scale factor (on 01/04/09).  Thus I have formally removed it from this
!                   version of the program.
!                 - Add SD, PSV to FAS, PSA output file and change order of FAS and PSA.
!        02/06/09 - Compute and write out husid time series and 95%-5% duration
!                 - Change input to be closer to smsim input
!                 - Read only output file stem name
!        02/14/09 - Correct error in setting i0, j0 if n_hypocenters = 1 and nsites > 1
!        02/16/09 - Change input, use stress scaling to adjust Wells and Coppersmith FL, FW
!                 - Compute RJB, RCD, print all distances for each site.
!                 - Include an option to move the sites to positions relative to
!                   the midpoint of the fault or to the tip of the fault
!                   if move_site == .true. and isitecoordflag == 2 (r and epi)
!                   (reset isitecoordflag to 3 in this case) 
!                   or if move_site == .true.  and isitecoordflag == 3 and one of the 
!                   site locations = 0.0
!                 - Moved writePar before loops, and just after call to getInputParameters.
!        02/17/09 - Modify findDistanceAndAzimuth subroutine
!                 - Modify findSubfaultDistance subroutine
!                 - Change FaultStrikeDeg to FaultStrike
!                 - Delete convertDegreeToRadian subroutine
!                 - Delete findDistance subroutine
!                 - Delete shortestSubfaultDistance subroutine
!                 - Read in dl, dw, and compute nl, nw
!        02/18/09 - Replace calls to get_time with a call to the system routine date_time
!        02/23/09 - Hardwire in choice of 75% rather than 95% in determining the duration.
!                   75% is recommended by Ou and Herrmann (SRL) and is consistent with a few
!                   Husid plots I made.
!        05/28/09 - Specify a hypocenter location by distances along the fault and down the dip.
!                 - Add period column to psa_fs output.
!        12/05/09 - Allow for up to a six-segment path duration (Atkinson and Boore 
!                   use a four-segment function).
!                   The form of the input matches that in the smsim programs.
!        12/06/09 - Pass path duration parameters through calling arguments rather than common.
!                 - Obtain number of frequencies in the crustal and site amp files from the 
!                   first line of the files (not from the params file).
!                 - Move the slip weight params before the site coord params and obtain slip weights
!                   from a text file (if slip weights are to be specified).
!                 - Modify order of site coord params so that it is easy to edit if do not want
!                   to do sims at all sites in the list (just change nsites; the list of sites need
!                   not be edited if the first nsites are those for which simulations are desired).
!                 - changed "y" to "vrup_beta" and "v" to "vrup" (the former names makes it very
!                   difficult to search for those variables).
!        11/20/11 - Add some Fortran 95 enhancements (e.g., named loops).  
!                 - Save and write to ioPSA the minimum, maximum, and std of the ground-motion 
!                   intensity measures.
!        11/21/11 - Change name and content of subroutine for accumulating sums for use in
!                   computing averages (to avoid repeating if statements for different types of
!                   averages).
!                 - Reduce the amount of information written to the screen
!        11/22/11 - Write elapsed time for each site to "parameters" file.
!                 - Comment out most non-error-related "Print" and "write(*" and "write (*" statements
!                 - Started changes required by a complete revision of the parameter file.
!                   This included: 
!                      Reordering the parameters into a more logical order, with the first
!                        set of parameters duplicating, to the extent possible, the SMSIM params.
!                      Including partition, free surface, radpat values as params rather than
!                        being hardwired in computeStochasticWave, with these values:
!                              prtitn=1.0/sqrt(2.0)
!                              rtp=0.55
!                              fs=2.0
!        11/23/11 - Use "include 'common_blocks.fi'" for the common blocks.
!        11/25/11 - NOTE: The statement "iseed = iabs(iseed)" has been removed from the main body
!                   of the program.  iseed needs to be a negative number before the first call to
!                   the random number routines (ran1 or gasdev) (see also the 11/25/08 revision).
!                   As it was, this was not the case
!                   if random slip weights were specified.  iseed is now set to a negative
!                   integer in getInputParameters, using the equation iseed = -int(abs(seed)),
!                   where "seed" is read in as a real number.
!        04/12/13 - Changed "fi2*180./pi" in call to writePSA_FA to "fi2", because fi2 in degrees.
!        02/04/14 - Increased dimension of siteLocation from (300,2) to (3000,2)
!        10/14/14 - Increase length of InputFileName
!        04/28/15 - Various small edits to correct errors found when Bob Herrmann compiled it using gfortran
!                   (but I used ! for all comments, which is proper Fortran 95 style, whereas Bob
!                   replaced ! and * with C).  I also changed (I think) all tabs to spaces using the 
!                   Search Replace... in TextPad, checking "regular expression" and replacing "\t"
!                   with " ".
!        05/13/15 - Various cosmetic code changes, including reducing the number of continuation lines
!                   for a number of statements (19 are allowed, although the program compiled and ran
!                   with no problem even with more continuation lines; this was caught by Visual Analyzer).
!                 - Write slipWeights in the f_h_f0 output file
!                 - Write nsims and n_hypocenters at the top of the dist_psa output file.
!        05/15/15 - Deleted declarations of many variables that Visual Analyzer found not to be used.
!        05/17/15 - In an attempt to avoid calling routines with different entry points, include the q(f)
!                   function here and compute the transition variables qt1, st in getInputParameters, passing them through
!                   the getInputParameterss argument list and then to sourceSpectra via the 
!                   common block in common_blocks.fi.
!                 - Replace Q parameters s1, s2, st with s1q, s2q, stq
!                 - Replace fault parameters s1, w1, s2, w2 with s1f, w1f, s2f, w2f
!		 10/17/16 - Edwards: changed velocity used by Q filter from beta to c_q
!				  - Added some more output of individual waveform simulations			

      real wave(:), totalWave(:), vel_total(:), husid(:)
      allocatable :: wave, totalWave, vel_total, husid 
      dimension slipWeights(1000,1000),
     :           accum_fas(500), accum_psa(500),
     :           averagePSA(500),
     :           averageFA(500),nnFA(500)
      dimension weightedMoment(1000,1000)
      dimension freq(500),psa(500),fa(500)
      
      real subfaultMoment(1000,1000), f0(1000,1000), risetime(1000,1000)
      integer NumberOfActiveSubs(1000,1000) 
c BE 01/12/16: freq_out(10) --> freq_out(50)
      real subfaultDistance(1000,1000), actualSlip(1000,1000), 
     :     dur_sub(1000,1000), dur_path(1000,1000),
     :     t_arrive(1000,1000), delay(1000,1000), freq_out(50), 
     :     freq4intrp(500), psa4intrp(500)
     
      real dur_75_05, arias
      
      real minPGA, maxPGA
      
      real alogFASsum(500), alogPSAsum(500)
      integer nFASsum(500), nPSAsum(500)
c BE 01/12/16: avgavgPSA_out(10)  --> avgavgPSA_out(50)
      real avgavgFAS(500), avgavgPSA(500), avgavgPSA_out(50)
c BE 01/12/16: f4head(10)*5 --> f4head(50)*5 
      character f4head(50)*5
      
      real amag, r_ref, rlow(10), a_s(10), b_s(10), m_s(10)

      real stutter
      
      dimension siteLocation(3000,2)
      character InputFileName*200
      character f_stem*120
      character fpar*120, facc*120, fpsa*120, f_h_f0*120, fhusid*120,
     :          f_times*120, f_dur_75_05*120, f_dist_psa*120
     
      logical f_exist, rmv_trend, specify_length, specify_width,
     :        move_site, write_site_files

      character fault_type*10

      character datx*8, time_start*10, time_stop*10

!BE 25/11/16      real rpathdur(6), pathdur(6)      
      real rpathdur(200), pathdur(200)      

      include 'common_blocks.fi'

      integer nptsTotalWave
!     nptsTotalWave is the total number of points in the final timeseries.

      nu_Read_Params = 99
      
      nu_Write_Params = nu_Read_Params - 1
      ioPSA = nu_Read_Params - 2
      ioHUS = nu_Read_Params - 3
      ioACC = nu_Read_Params - 4
      nu_h = nu_Read_Params - 8
      nu_times = nu_Read_Params - 10
      nu_dur_75_05 = nu_Read_Params - 11
       
      nu_dist_psa = nu_Read_Params - 20
      


      pi=4.0*atan(1.0)
      twopi = 2.0 * pi

      f_exist = .false.
      do while (.not. f_exist)
        InputFileName = ' '
        write(*, '(a)') 
     :    ' Enter name of input parameter file '//
     :    '(Enter = exsim_dmb.params): '
        read(*, '(a)') InputFileName
        if (InputFileName(1:4) == '    ') 
     :                 InputFileName = 'exsim_dmb.params'
        call trim_c(InputFileName, nc_f_in)
        inquire(file=InputFileName(1:nc_f_in), exist=f_exist)
        if (.not. f_exist) then
          write(*,'(a)') ' ******* FILE '//
     :                     InputFileName(1:nc_f_in)//
     :                   ' DOES NOT EXIST ******* '
        end if
      end do
      open(unit=nu_Read_Params, file=InputFileName(1:nc_f_in), 
     :     status='unknown')

!     Read (and echo) the input parameters that are common to all grid points.
      CALL GETINPUTPARAMETERS(nu_Read_Params,rho,beta,prtitn,rtp,fs,
     :     r_ref, nsprd_segs, rlow, a_s, b_s, m_s,fr1, qr1, s1q, ft1, 
     :     ft2, fr2, qr2, s2q, c_q, 
     :     qt1, stq, ndur_hinges, rpathdur, pathdur, 
     :     durslope, fmax, akappa_0, dkappadmag, amagkref, akappa,
     :     flocut, nslope, iwind, taper, eps_w, eta_w, f_tb2te, 
     :     f_te_xtnd, tpadl, tpadt, dt, iseed, nsims, amag, stress,
     :     FaultLat,FaultLon,FaultStrike,FaultDip,h,fault_type,
     :     FaultLength,FaultWidth,dl,dw,stress_ref,nl,nw,nsubs, 
     :     specify_length, specify_width,vrup_beta,
     :     hyp_loc_fl, hyp_loc_fw, n_hypocenters,i0_in, j0_in, 
     :     i_rise_time, iDynamicFlag, pulsingPercent, 
     :     iflagscalefactor, islipweight, slipWeights,
     :     iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :     f_crustal_amps, f_site_amps,
     :     iflagfas_avg,iflagpsa_avg_over_sims,iflagpsa_avg_over_hypos,
     :     write_site_files, f_stem,
     :     damp, nfreq,freq1,freq2, nfout, freq_out,
     :     isitecoordflag, move_site, numberOfSites, siteLocation)

      close(nu_Read_Params)
      
!DEBUG
      print *,' amag, stress = ', amag, stress
!DEBUG
      
      vrup=vrup_beta * beta

      subFaultRadius=sqrt((dl*dw)/pi)

!%%%%%%%%%%
      riseTime_original=subFaultRadius/vrup
!%%%%%%%%%%

      CALL GETAMPS(nu_Read_Params, 
     :   f_crustal_amps,
     :   n_crustal_amps, freq_crustal_amps, amp_crustal_amps, 
     :   f_site_amps,
     :   n_site_amps, freq_site_amps, amp_site_amps) 
     
      call trim_c(f_stem, nc_f_stem)
      
      fpar = ' '
      fpar = f_stem(1:nc_f_stem)//'_parameters.out'
      call trim_c(fpar, nc_fpar)
      open (nu_Write_Params,file=fpar(1:nc_fpar),status='unknown')

      CALL WRITEPAR(nu_Write_Params,
     :  FaultStrike,FaultDip,h, FaultLat, FaultLon, 
     :  siteLocation,numberOfSites, move_site,
     :  fault_type, FaultLength, FaultWidth,nl,nw,dl,dw, 
     :  specify_length, specify_width, stress_ref,
     :  vrup_beta,
     :  hyp_loc_fl, hyp_loc_fw, i0_in,j0_in,n_hypocenters,
     :  nsims,
     :  amag,
     :  stress, 
     :  pulsingPercent,iDynamicFlag, i_rise_time,
     :  iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :  iflagscalefactor,   
     :  iflagfas_avg, iflagpsa_avg_over_hypos, iflagpsa_avg_over_sims,
     :  tpadl, tpadt,
     :  r_ref, nsprd_segs, rlow, a_s, b_s, m_s,
     :  rpathdur, pathdur, durslope, ndur_hinges 
     :                   )
      
! Initialize q_f
!      dummy = q_f_setup(fr1, qr1, s1q, ft1, ft2,
!     :                  fr2, qr2, s2q, c_q)

! Initialize gsprd_q_f:

      dummy = gsprd_f_setup(r_ref,nsprd_segs,rlow,
     :                  a_s,b_s,m_s,
     :                  amag)          ! remove numsource from argument list
     
      
      npadl = tpadl/dt
      npadt = tpadt/dt


      f_h_f0 = ' '
      f_h_f0 = f_stem(1:nc_f_stem)//'_h_f0.out'
      call trim_c(f_h_f0, nc_f_h_f0)
      open(nu_h, file=f_h_f0(1:nc_f_h_f0), status='unknown')
      write(nu_h,'(2x,a, 1x,a, 2x, a, 2x,a, 3x, a, 3x,a, 1x,a, 
     :  1x,a, 1x,a,
     :  1x,a, 1x,a,
     :  5x,a, 1x,a, 9x,a, 
     :  8x,a, 1x,a,
     :  3x,a, 4x,a, 7x,a, 
     :  1x,a, 4x,a,
     :  1x,a,
     :  2x,a)') 
     :      'isite', 'ihyp', 'i0', 'j0', 'i', 'j', 'nsubs', 
     :      'slipWeights', 'Diff_T_arrive',
     :      'NoOfEffectiveSubfaults', 'NumberOfActiveSubs',
     :      'f0main', '1stElmntF0', 'f0', 
     :      'hij', 'iflagscalefactor', 
     :      '1/f0main', '1/1stF0', '1/f0', 
     :      'riset_orig', 'dur_sub',
     :      'subFaultM0',       
     :      'Tot/subM0'
 
      f_dur_75_05 = ' '
      f_dur_75_05 = f_stem(1:nc_f_stem)//'_dur_75_05.out'
      call trim_c(f_dur_75_05, nc_f_dur_75_05)
      open(nu_dur_75_05, file=f_dur_75_05(1:nc_f_dur_75_05), 
     :             status='unknown')
      write(nu_dur_75_05,'(2x,a,
     :             1x,a,
     :             2x, a, 2x,a, 
     :             1x,a, 
     :             1x,a)') 
     :                 'isite',
     :                 'ihyp',
     :                 'i0', 'j0',
     :                 'dur_75_05',
     :                 'arias(cm/s)'
 
!     Generate logarithmically-spaced frequencies
      freq(1)=freq1
      if (nfreq>=2) then
         frinc=alog(freq2/freq1)/(nfreq-1)
         do k=2,nfreq
            freq(k)=freq1*exp((k-1)*frinc)
         enddo
      endif
!     End of Generate logarithmically-spaced frequencies

!      iseed = -iabs(iseed)  ! See 11/25/11 revision.

      if(i0_in==0.or.j0_in==0)
     :    CALL LOCATEHYPOCENTRERANDOMLY(i0,j0,nl,nw)

      NoOfEffectiveSubfaults=real(nl)*(pulsingPercent/100.0)
      NoOfEffectiveSubfaults=NoOfEffectiveSubfaults/2. !!double side propogation
      if(NoOfEffectiveSubfaults<=1) NoOfEffectiveSubfaults=1

      totalMoment=10.**(1.5*amag+16.05)
      avgSubFaultMoment=totalMoment/real(nw*nl)

      firstElementF0=  
     :                4.9e+6*beta*(stress/avgSubFaultMoment)**(1.0/3.0)

      F0main=4.9e+6*beta*(stress/totalMoment)**(1.0/3.0)

      CALL FINDSUBFAULTMOMENT(nl,nw,slipWeights
     :                 ,totalMoment,weightedMoment)

      f_dist_psa = ' '
      f_dist_psa = f_stem(1:nc_f_stem)//'_distances_psa.out'
      call trim_c(f_dist_psa, nc_f_dist_psa)
      open(nu_dist_psa, file=f_dist_psa(1:nc_f_dist_psa), 
     :     status='unknown')
    
      f4head = '00000'
      
      do i = 1, nfout
        if (freq_out(i) < 0.0 .or. freq_out(i) >= 10.0) then
           write(f4head(i)(1:5),'(f5.2)') freq_out(i)
        else  
           write(f4head(i)(2:5),'(f4.2)') freq_out(i)
        end if
      end do
        
c    BE 02/12/16 : 10(7x,a4, 3x,a8))' --> 50(7x,a4, 3x,a8))'
      write(nu_dist_psa,'(a, 1x,i3, a, 1x,i3, a )') 
     : ' Average motions, where for each site the motions are '//
     :  'averaged over', nsims, ' simulations and', n_hypocenters, 
     : ' hypocenters (more than 1 when random hypocenters are '//
     :  'specified in the input params file).'
      write(nu_dist_psa,'(1x,a, 
     :                1x,a, 1x,a, 1x,a, 
     :                1x,a, 1x,a, 
     :                5x,a, 3x,a,
     :                50(7x,a4, 3x,a8))')
     :  'isite',
     :  'sitecoord(1)', 'sitecoord(2)', 'isitecoordflag',
     :  'site_lat_degrees', 'site_lon_degrees', 
     :  'd_jb', 'd_cd2f',
     :  ('freq', 'psa'//f4head(i), i = 1, nfout)
     
      write(nu_Write_Params,*)
 
!___________________________________________________________________________________________________
!///////////////////////////////////////////|\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!---------------------------------------------------------------------------------------------------
      loop over number of sites: DO isite=1,numberOfSites   ! Loop on Sites around the Fault
      
! Standard Fortran 90 intrinsic Subroutine DATE_AND_TIME
        call DATE_AND_TIME( datx, time_start )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'
       
        SiteLat=siteLocation(isite,1)
        SiteLon= siteLocation(isite,2)
        print *, "Working on Site #",isite,SiteLat, SiteLon

! Compute some distances: 

        CALL SITE_COORDS_IN_DEGREES(
     :         FaultLat, FaultLon, 
     :         SiteLat, SiteLon, isitecoordflag, 
     :         site_lat_degrees, site_lon_degrees ) 

        h_min_c = 3.0 ! Campbell depth to seismogenic region
        w1f = 0.0
        w2f = faultwidth
        s1f = 0.0
        s2f = faultlength
     
        CALL DIST_3DF(
     :   site_lat_degrees, site_lon_degrees, 
     :   FaultLat, FaultLon, h, FaultStrike, FaultDip,
     :   w1f, w2f, s1f, s2f, 
     :   h_min_c, d_jb, az_jb, d_cd2f, az_cd2f, d_c, az_c,
     :   d_sta_n, d_sta_e, irgn_cd2f, irgn_c, irgn_jb)
 
        nnFA=0
        fi2=0
        CALL FINDDISTANCEANDAZIMUTH(FaultLat,FaultLon,SiteLat,
     :                           SiteLon,R,fi2,isitecoordflag)
!     Note that R, azm fi2 are the distance, azimuth (in degrees). w.r.t. origin (not epicenter),
!     Note also I must input faultstrike in degrees and its converted to rad.
 
        alogFASsum = 0.0
        alogPSAsum = 0.0
           nFASsum = 0
           nPSAsum = 0
        alogPGAsum = 0.0
        alogPGVsum = 0.0
        
        loop over hypocenters: DO ihypo = 1, n_hypocenters
        

          if (n_hypocenters > 1) then
            call locateHypocentreRandomly(i0,j0,nl,nw)
          else if (i0_in == 0 .or. j0_in == 0) then
            call locateHypocentreRandomly(i0,j0,nl,nw)
          else
            i0 = i0_in
            j0 = j0_in
          end if

          print *, ' Working on hypocenter # ',ihypo
          print *, ' located at i0, j0 = ',i0, j0
          
          nsubfaults = 0
          
          t_arrive_min=10000.
          t_end_max=0.
          risetime_max = 0.0
            
          loop on nl for non nsims quantities: DO i=1,nl
            loop on nw for non nsims quantities: DO j=1,nw
            
              nsubfaults = nsubfaults + 1
              
              subFaultMoment(i,j)=weightedMoment(i,j)
                
              if(iDynamicFlag==0) then
                NumberOfActiveSubs(i,j)=1
              else
                NumberOfActiveSubs(i,j)=NumberOfPulsingSubs(i0,j0,i,j
     :                 ,nl,nw,NoOfEffectiveSubfaults)
                if(NumberOfActiveSubs(i,j)==0) then
                  NumberOfActiveSubs(i,j)=1
                end if
              end if
                
              f0(i,j)=
     :           firstElementF0*(NumberOfActiveSubs(i,j)**(-1.0/3.0))
              if (i_rise_time == 1) then
                risetime(i,j) = riseTime_original                
              else if (i_rise_time == 2) then
                riseTime(i,j) = 1.0/f0(i,j)
              else
                print *,' ERROR: invalid value of i_rise_time, = ',
     :                    i_rise_time
                stop
              end if
              
              if(risetime(i,j) > risetime_max) then
                risetime_max = risetime(i,j)
              end if
 
              actualSlip(i,j)=1.e-22*subFaultMoment(i,j)/
     :                            (rho*beta**2.*dl*dw)
! The actual slip is not written to any output files in this version.

              subfaultDistance(i,j)=
     :                 findSubfaultDistance(R,h,FaultStrike,
     :                           fi2,FaultDip,dl,dw,i,j)
     
!     calculate duration of subsource time history at distance
!     subfaultDistance
              call dur_path_cmp(subfaultDistance(i,j),
     :                        rpathdur, pathdur, durslope, ndur_hinges,
     :                  dur_path(i,j))  ! computes path part of duration

              dur_sub(i,j) = dur_path(i,j) + risetime(i,j)
              
              if ((i-i0) /= 0 .or. (j-j0) /= 0) then
                delay(i,j)=sqrt((dl*(i-i0))**2.+(dw*(j-j0))**2.)/vrup
              else
                delay(i,j)=0.  ! delay from hypocenter to itself
              end if
              
              t_arrive(i,j) = delay(i,j) + subfaultDistance(i,j)/beta            
              if (t_arrive(i,j) < t_arrive_min) 
     :                    t_arrive_min = t_arrive(i,j)
               
              t_end= t_arrive(i,j) + dur_sub(i,j)            
              if (t_end > t_end_max) then
                imax = i
                jmax = j
                t_end_max = t_end
              end if
            
            END DO loop on nw for non nsims quantities 
          END DO loop on nl for non nsims quantities 

          nptsTotalWave=(t_end_max - t_arrive_min + tpadl + tpadt + 
     :                        risetime_max)/dt   ! need risetime_max because the stutter
                                                 ! delay can shift the subsource time series
                                                 ! corresponding to t_end_max by an amount up to
                                                 ! risetime of that subsource; to be safe I
                                                 ! use risetimemax rather than risetime(imax,jmax),
                                                 ! where imax, jmax are the indices of the subsource
                                                 ! corresponding to t_end_max.
                                                 
          signnpw2 = +1.0
          call get_npw2(nptsTotalWave,signnpw2,npw2TotalWave)
          
          if (write_site_files) then
          
            if (ihypo == 1) then
 
              f_times = ' '
              f_times = f_stem(1:nc_f_stem)//'_times_s000.out'
              write(f_times(nc_f_stem+9:nc_f_stem+11),'(i3.3)') isite
              call trim_c(f_times, nc_f_times)         
              open (nu_times,file=f_times(1:nc_f_times),
     :              status='unknown')
 
              write(nu_times,'(1x,a,i8, 1x,a,i8)') 
     :          'nptsTotalWave = ', nptsTotalWave, 
     :          'npw2TotalWave = ', npw2TotalWave
              write(nu_times,*)


              write(nu_times,'(1x,a, 2x,a, 2x,a, 3x,a, 3x,a,
     :                  1x,a, 1x,a, 
     :                  1x,a, 1x,a, 1x,a,
     :                  1x,a, 1x,a,
     :                  1x,a,
     :                  1x,a, 1x,a, 1x,a,
     :                  1x,a,
     :                  1x,a,
     :                  1x,a, 1x,a)')   
     :        'itr', 'i0', 'j0', 'i', 'j', 
     :        't_arrive(i0,j0)', 't_arrive_min', 
     :        't_arrive(i,j)', 'risetime(i,j)', 'risetime_max', 
     :        'dur_path(i,j)', 'dur_sub(i,j)', 
     :        't_arrive(i,j)+dur_sub(i,j)',
     :        'imx', 'jmx', 't_arrive(imax,jmax)', 
     :        'dur_sub(imax,jmax)', 
     :        't_arrive(imax,jmax)+dur_sub(imax,jmax)', 
     :        'risetime_max', 't_end_max' 
     
              do i = 1, nl
                do j = 1, nw
              
                  write(nu_times,'(5(1x,i3), 
     :                      8x,f8.1, 5x,f8.1,
     :                      6x,f8.1, 6x,f8.1, 5x,f8.1, 
     :                      5x,f9.2, 4x,f9.2,
     :                      19x,f8.1,
     :                      2(1x,i3), 12x,f8.1,
     :                      11x,f8.1,
     :                      31x,f8.1,
     :                      5x,f8.1, 2x,f8.1)') 
     :          ihypo, i0, j0, i, j, 
     :          t_arrive(i0, j0), t_arrive_min, 
     :          t_arrive(i,j), risetime(i,j), risetime_max, 
     :          dur_path(i,j), dur_sub(i,j), 
     :          t_arrive(i,j)+dur_sub(i,j),
     :          imax, jmax, t_arrive(imax,jmax), 
     :          dur_sub(imax,jmax), 
     :          t_arrive(imax,jmax)+dur_sub(imax,jmax), 
     :          risetime_max, t_end_max 
     
                end do
              end do
            
              close(nu_times)
       
              end if
              
          end if
          
          HypoDistance=findSubfaultDistance(R,h,FaultStrike,
     :                           fi2,FaultDip,dl,dw,i0,j0)
     
!     now do actual generation and summation of subfault time series ***

! Initialize arrays (index over frequency) containing cumulative sums for averages
          averagePSA=0.0
          averageFA=0.0
          accum_fas = 0.0
          accum_psa = 0.0
          nnFA = 0
          accum_pga = 0.0
          accum_pgv = 0.0
          
          accum_dur_75_05 = 0.0
          accum_arias = 0.0
          
          minPGA = 1.0e20
          maxPGA = 0.0
          
          loop over nsims for a hypocenter: DO isim=1,nsims
          
!            print *, ' Site isite, Hypocenter iteration ihypo, '//
!     :               'Simulation number isim = ', isite, ihypo, isim

            allocate (totalWave(npw2TotalWave))
            allocate (husid(npw2TotalWave))
          
            totalWave = 0.0

! For a given hypocenter and simulation, begin loop over subfaults:
            DO i=1,nl
              DO j=1,nw
            
                nsubsource = npadl + dur_sub(i,j)/dt + npadt
                signnpw2 = +1.0
                call get_npw2(nsubsource,signnpw2,iFFTsub) 

                allocate (wave(iFFTsub))
                
                wave=0
                
                call computeStochasticWave(
     :              wave,
     :              subFaultMoment(i,j), subfaultDistance(i,j),
     :              F0Main, f0(i,j), dur_sub(i,j),
     :              iflagscalefactor, hij)                            ! dmb
      
! NOTE: The length of wave is iFFTsub; this is passed through common /par/

                if (isim == 1) then
                  
                  write(nu_h,'(4x,i3,
     :            2x,i3, 1x,i3, 1x,i3, 
     :            1x,i3, 1x,i3, 1x,i5,1x,es11.4,
     :            4x,es10.3, 
     :            18x,i5,14x,i5,
     :            4(1x,es10.3), 
     :            16x,i1,
     :            7(1x,es10.3))') 
     :            isite, ihypo, i0, j0,
     :            i, j, nsubfaults,slipWeights(i,j),
     :            t_arrive(i,j)-t_arrive(i0,j0),
     :            NoOfEffectiveSubfaults,NumberOfActiveSubs(i,j), 
     :            f0main, firstElementF0, f0(i,j), 
     :            hij, iflagscalefactor,  
     :            1.0/f0main, 1.0/firstElementF0, 1.0/f0(i,j),
     :               risetime_original, dur_sub(i,j),
     :               subFaultMoment(i,j), 
     :               TotalMoment/subFaultMoment(i,j)
                end if
                
                
!               subsource accelerogram is randomly delayed to
!               simulate complexity in slip process
!               subtract minimum delay to eliminate static time shift

                 
                stutter = Ran1(iseed) * riseTime(i,j)
                ishift = ((t_arrive(i,j) - t_arrive_min) + stutter)/dt
                
                if (npadl+ishift+dur_sub(i,j)+npadt > 
     :                                             npw2TotalWave) then
!                   print *,'npadl+ishift+dur_sub(i,j)+'//
!     :                     'npadt>npw2TotalWave'
!                   print *,  npadl+ishift+dur_sub(i,j)+npadt, 
!     :                       npw2TotalWave
                end if   
                 
! No need to call shift---just figure out the shifted index and use it  
! to fill the totalWave vector.
                 
!               add current accelerogram to total wave field

                nhigh = min(iFFTsub, npw2TotalWave - ishift)

                DO k = 1, nhigh
                  totalWave(k+ishift)=totalWave(k+ishift)+wave(k)
                END DO
                
                deallocate (wave)
                
              END DO ! loop over nw
            END DO ! loop over nl
! For a given hypocenter and simulation, end loop over subfaults:
          
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     include Mavroeidis and Papageorgiou, 2003 approach
            if( int(iPapaFlag)==1)then
              call MavroPapa(totalWave,npw2TotalWave,
     *           PapaGama,PapaNu,PapaT0,amag,PapaPeak)
            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!     ************************** end of summation **********************

!     accelerogram and model parameters are saved only once

            pga=findPeakValue(npw2TotalWave,totalWave)
            if(pga < minPGA) then
              minPGA = pga
            end if
            if(pga > maxPGA) then
              maxPGA = pga
            end if
            call sum4AverageY(accum_pga, pga, iflagpsa_avg_over_sims)
!DEBUG
!      print *,' pga, minPGA, maxPGA = ', pga, minPGA, maxPGA
!DEBUG
      
            allocate (vel_total(npw2TotalWave))
            vel_total = 0.0

          
            call acc2v(totalWave, npw2TotalWave, dt, .false., vel_total)

            pgv=findPeakValue(npw2TotalWave,vel_total)
            call sum4AverageY(accum_pgv, pgv, iflagpsa_avg_over_sims)
            
            deallocate(vel_total)

            rmv_trend = .false.
            husid = 0.0 
            call accsqint(totalWave, npw2TotalWave, dt, rmv_trend, 
     :                    husid)
            husid_end = husid(npw2TotalWave)
            do i = 1, npw2TotalWave
              husid(i) = husid(i)/husid_end
            end do
      
            call locate(husid,npw2TotalWave,0.05,j05) 
            call locate(husid,npw2TotalWave,0.75,j75) 
        
            accum_dur_75_05 = accum_dur_75_05 + 
     :                        alog10(real(j75-j05)*dt)  ! for geometric mean
            
            accum_arias = accum_arias + 
     :                    alog10((pi/(2.0*981.0))*husid_end)  ! for geometric mean 
              

          if (write_site_files) then

c edwards edit 171016
c write everything
c            if (ihypo.eq.1 .and. isim .eq. 1) then
c end edit            
              facc = ' '
c edwards edit 171016
c              facc = f_stem(1:nc_f_stem)//'_acc_s000.out'
c              write(facc(nc_f_stem+7:nc_f_stem+9),'(i3.3)') isite
c add hypo and sim #
              facc = f_stem(1:nc_f_stem)//'_acc_hypo000_sim000_s000.out'
              write(facc(nc_f_stem+10:nc_f_stem+12),'(i3.3)') ihypo
              write(facc(nc_f_stem+17:nc_f_stem+19),'(i3.3)') isim
              write(facc(nc_f_stem+22:nc_f_stem+24),'(i3.3)') isite
c end edit
 
              call trim_c(facc, nc_facc)         
              open (ioAcc,file=facc(1:nc_facc),status='unknown')
      
              call writeAcc(ioAcc, totalWave, npw2TotalWave, pga,
     :                 isim,isite,HypoDistance,d_cd2f, fpar)
              close(ioAcc)
              
              fhusid = ' '
c edwards edit 171016
c              fhusid = f_stem(1:nc_f_stem)//'_husid_s000.out'
c              write(fhusid(nc_f_stem+9:nc_f_stem+11),'(i3.3)') isite
c add hypo and sim #
              fhusid = f_stem(1:nc_f_stem)//
     :'_husid_hypo000_sim000_s000.out'
              write(fhusid(nc_f_stem+12:nc_f_stem+14),'(i3.3)') ihypo
              write(fhusid(nc_f_stem+19:nc_f_stem+21),'(i3.3)') isim
              write(fhusid(nc_f_stem+24:nc_f_stem+26),'(i3.3)') isite
c end edit

              call trim_c(fhusid, nc_fhusid)         
              open (ioHus,file=fhusid(1:nc_fhusid),status='unknown')
              
              dur_75_05_sim1 =  real(j75-j05)*dt 
            
              arias_sim1 = (pi/(2.0*981.0))*husid_end
              
              call writeHUS(ioHus, husid, npw2TotalWave, 
     :                 arias_sim1, dur_75_05_sim1, 
     :                 isim,isite,HypoDistance,d_cd2f,fpar)     
              close(ioHUS)
c edwards edit 171016             
c            end if
c end edit

          end if

            dur=real(npw2TotalWave-1)*dt-5.0
            psa=0.0
            call computeResponse(totalWave,dt,dur,damp,nfreq,freq,psa)
            do ifreq = 1, nfreq
              call sum4AverageY(accum_psa(ifreq),psa(ifreq), 
     :                              iflagpsa_avg_over_sims)  ! dmb
            end do
            FA=0.0
            call computeFACCN(totalWave,npw2TotalWave,dt,nfreq,freq,FA)
            call sum4AverageFA(accum_fas,FA,nfreq,nnFA,
     :                                 iflagfas_avg)
     
!            write(*,"(' *** finished with trial',i4)") isim
            
            deallocate(totalWave)
            deallocate(husid)
         
          END DO loop over nsims for a hypocenter ! End loop over nsims for a given hypocenter

          call computeAverageY(accum_pga, averagePGA, 
     :                         nsims, iflagpsa_avg_over_sims)
          alogPGAsum = alogPGAsum + alog10(averagePGA)  ! In a future revision, allow for a different average over hypocenters
          call computeAverageY(accum_pgv, averagePGV, 
     :                         nsims, iflagpsa_avg_over_sims)
          alogPGVsum = alogPGVsum + alog10(averagePGV)


          DO ifreq = 1, nfreq
            call computeAverageY(accum_psa(ifreq), averagePSA(ifreq), 
     :                         nsims, iflagpsa_avg_over_sims)
            
            if (iflagfas_avg == 2 .and. accum_fas(ifreq) == 0.0) then   
              averageFA(ifreq) = -9.99
            else
              if (iflagfas_avg == 1) then  ! arithmetic
                averageFA(ifreq) = accum_fas(ifreq)/
     :                               real(nnFA(ifreq))
              else if (iflagfas_avg == 2) then  ! geometric
                averageFA(ifreq) = 10.0**(accum_fas(ifreq)/
     :                               real(nnFA(ifreq)))
              else if (iflagfas_avg == 3) then  ! rms
                averageFA(ifreq) = sqrt(accum_fas(ifreq))
              else
                print *,
     :         ' ERROR, iflagfas_avg not 1, 2, or 3; = ', iflagfas_avg
                stop
              end if
            end if
         
            if (averagePSA(ifreq) > 0.0) then
              alogPSAsum(ifreq) = alogPSAsum(ifreq) + 
     :                              alog10(averagePSA(ifreq))
              nPSAsum(ifreq) = nPSAsum(ifreq) + 1
            end if
             
            if (averageFA(ifreq) > 0.0) then
              alogFASsum(ifreq) = alogFASsum(ifreq) + 
     :                              alog10(averageFA(ifreq))
              nFASsum(ifreq) = nFASsum(ifreq) + 1
            end if
             
          END DO  ! loop over nfreq
 
          dur_75_05 = 10.0**(accum_dur_75_05/real(nsims))
          arias     = 10.0**(accum_arias/real(nsims))
          
          write(nu_dur_75_05,'(4x,i3,
     :                         2x,i3,
     :                         1x,i3, 1x,i3, 
     :                         1x,f9.2, 
     :                         1p, 1x,e10.3)') 
     :                         isite,
     :                         ihypo,
     :                         i0, j0,
     :                         dur_75_05,
     :                         arias


        END DO loop over hypocenters !loop over random hypocenters
       
! Compute average of average FAS and PSA


        avgavgPGA = 10.0**(alogPGAsum/real(n_hypocenters))
        avgavgPGV = 10.0**(alogPGVsum/real(n_hypocenters))

        DO ifreq = 1, nfreq
       
         avgavgPSA(ifreq) = 10.0**(alogPSAsum(ifreq)/
     :                        real(nPSAsum(ifreq)))
c edwards 171016 edit - avoid /0 where FAS not defined (very short periods)
c         avgavgFAS(ifreq) = 10.0**(alogFASsum(ifreq)/
c     :                        real(nFASsum(ifreq)))
         if (nFASsum(ifreq).gt.0) then
           avgavgFAS(ifreq) = 10.0**(alogFASsum(ifreq)/
     :                        real(nFASsum(ifreq)))
         else 
           avgavgFAS(ifreq) = NaN
        end if
c end edit
        END DO ! end loop over nfreq
       
        nfreq4intrp = nfreq
        do i = 1, nfreq4intrp
          freq4intrp(i) = freq(i)
          psa4intrp(i) = avgavgPSA(i)
        end do
        
        do i = 1, nfout
        
          if (freq_out(i) == -1.0) then ! pgv
            avgavgPSA_out(i) = avgavgPGV
          else if (freq_out(i) >= 99.0) then ! pga
            avgavgPSA_out(i) = avgavgPGA
          else
            avgavgPSA_out(i) = 
     :       yintrf( freq_out(i), freq4intrp, psa4intrp, nfreq4intrp)
          end if
          
        end do
 
        if (write_site_files) then

          fpsa = ' '
          fpsa = f_stem(1:nc_f_stem)//'_psa_fs_s000.out'
          write(fpsa(nc_f_stem+10:nc_f_stem+12),'(i3.3)') isite
          call trim_c(fpsa, nc_fpsa)         
          open (ioPsa,file=fpsa(1:nc_fpsa),status='unknown')

        call writePSA_FA(ioPSA, freq,avgavgPGA,avgavgPGV,
     :       minPGA, maxPGA, 
     :       avgavgPSA,nfreq,nsims,damp,amag,
     :       siteLat,siteLon,h,r,fi2,d_cd2f,avgavgFAS,
     :       isite, HypoDistance,fPar, nl, nw)
          close(ioPSA)
         
        end if
c BE 02/12/16 : 10(1x,e10.3, 1x,e10.3))' --> 50(1x,e10.3, 1x,e10.3))'
        write(nu_dist_psa, '(
     :    1x,i5,
     :    4x,f9.3, 4x,f9.3, 14x,i1,
     :    8x,f9.3, 8x,f9.3, 
     :    1x,f8.2, 1x,f8.2,  
     :    1p,50(1x,e10.3, 1x,e10.3))') 
     :    isite, 
     :    SiteLat,  SiteLon, isitecoordflag, 
     :    site_lat_degrees,  site_lon_degrees,
     :     d_jb, d_cd2f,
     :     (freq_out(i), avgavgPSA_out(i), i= 1, nfout)

       
        call DATE_AND_TIME( datx, time_stop )
! Date is returned as 'CCYYMMDD'
! Time is returned as 'hhmmss.sss'
        call time_diff(time_start, time_stop, time_elapsed)
        write(nu_Write_Params,'(1x,a, 1x,i2, a,f6.2)')  
     :       ' For site ', isite, ' elapsed time (sec) = ', time_elapsed

      END DO  loop over number of sites ! End loop over sites
!___________________________________________________________________________________________________
!///////////////////////////////////////////|\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!---------------------------------------------------------------------------------------------------

      close(nu_Write_Params) 
      
      close(nu_h)
      close(nu_dur_75_05)
      close(nu_dist_psa)

 
      STOP
      END
!ccccccccccccccccccccccc   End of Main Program    cccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!------------------- SUBPROGRAMS -----------------------------------------

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real function amplf (f,afreq,amp,namp)
!       computes amplification value at frequency f by linear
!       interpolation between its adjoining values
        dimension amp(*),afreq(*)

        if (f <= afreq(1)) then
            amplf=amp(1)
            return
        endif

        if (f >= afreq(namp)) then
            amplf=amp(namp)
            return
        endif

        do i=2,namp
            if (f <= afreq(i)) then
                d=(amp(i)-amp(i-1))/(afreq(i)-afreq(i-1))*(f-afreq(i-1))
                amplf=amp(i-1)+d
                return
            endif
        enddo

        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine avg(n,s,avamp2)
!     calculates average squared amplitude spectrum
      complex s
      dimension s(*)
      sum=0.
      do j=1,n
          amp=cabs(s(j))
          sum=sum+amp*amp
      enddo
      avamp2=sum/real(n)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sum4AverageFA(accum_fas,FA,nfreq,nn,
     :                                  iflagfas_avg)
!     calculate average Fourier spectrum    ! revised by dmb
! Dates: 05/15/15 - Removed nsims from argument list
      dimension accum_fas(*),FA(*), nn(500)
      DO i=1,nfreq
        if(FA(i) > 0)then
          nn(i)=nn(i)+1
          if      (iflagfas_avg == 1) then ! arithmetic
                accum_fas(i)=
     :              accum_fas(i) + FA(i)
          else if (iflagfas_avg == 2) then ! geometric
                accum_fas(i)=
     :              accum_fas(i) + alog10(FA(i))
          else if (iflagfas_avg == 3) then ! rms
            accum_fas(i)=
     :         (accum_fas(i)*real(nn(i)-1) + FA(i)**2.0 )/real(nn(i))
!          else if (iflagfas_avg == 3) then ! rms
!                accum_fas(i)=
!     :              accum_fas(i) + FA(i)**2.0
          else
            print *,
     :         ' ERROR, iflagfas_avg not 1, 2, or 3; = ', iflagfas_avg
            stop
          end if     
        endif
      ENDDO
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sum4AverageY(accum_y, y, iflagpsa_avg)
! Calculate accumulated sum for average response spectrum   
! Dates: 11/21/11 - Written by DMB
      real accum_y, y
      integer iflagpsa_avg
      if          (iflagpsa_avg == 1) then ! arithmetic
          accum_y=
     :              accum_y + y
          else if (iflagpsa_avg == 2) then ! geometric
                accum_y=
     :              accum_y + alog10(y)
          else if (iflagpsa_avg == 3) then ! rms
                accum_y=
     :              accum_y + y**2.0
          else
            print *,
     :         ' ERROR, iflagpsa_avg not 1, 2, or 3; = ', iflagpsa_avg
            stop
      end if     
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeAverageY(accum_y, averageY, ny, iflagpsa_avg)
! Compute average of y over ny   
! Dates: 11/21/11 - Written by DMB
      real accum_y, averageY
      integer iflagpsa_avg
      if (iflagpsa_avg == 2 .and. accum_y == 0.0) then
        averageY = -9.99
      else
        if      (iflagpsa_avg == 1) then  ! arithmetic
          averageY = accum_y/real(ny)
        else if (iflagpsa_avg == 2) then  ! geometric
          averageY = 10.0**(accum_y/real(ny))
        else if (iflagpsa_avg == 3) then  ! rms
          averageY = sqrt(accum_y/real(ny))
        else
          print *,
     :         ' ERROR, iflagpsa_avg not 1, 2, or 3; = ', iflagpsa_avg
          stop
        end if
      end if
       
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeFACCN(totalWave,npts,dt,nfreq,freq, FA)
!
!       Compute Fourier acceleration spectrum for final timeseries.
      dimension totalWave(*), freq(500), FA(500)
      complex fft(:)
      real spectrum(:)
      allocatable :: fft, spectrum
      
      allocate(fft(npts), spectrum(npts))
      
      fft = 0.0

!     Copy totalWave into new array because we need to save original for later.
      do j=1, npts
          spectrum(j)=totalWave(j)
      enddo

      call makeItComplex(spectrum,fft,npts)

!     Take padded time series into freq. domain
      call fork(npts,fft,-1.)
      ncall = npts/2 +1
      df = 1./(dt * real(npts))

!     Scale spectrum and put value back into spectrum array.
!      fft(ncall) = cmplx(real(fft(ncall)),0.)
      fft(ncall) = 0.0  ! ncall = Nyquist, set value to 0.0
      
      do jf = 1, ncall
        spectrum(jf) = dt* sqrt(real(npts)) * cabs(fft(jf))
      enddo

!     Now sample the spectrum into the same freq. bins used for PSA output.
      call sample(spectrum, ncall, df, nfreq, freq, FA)

!     Array FA now contains the nfreq values of sampled faccn.
!     TotalWave input array is unchanged.

      deallocate(fft, spectrum)
      
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeResponse(a,dt,dur,damp, Nfreq,freq,psa)

!     calculates response spectrum from acceleration time history
! Dates: 11/23/08 - Some changes by DMB (compute pi, use rdcalcdp)
!        12/06/09 - Changed "beta" to "damp_fraction" ("beta" is used elsewhere for
!                   shear-wave velocity)
!        11/23/11 - Call rscalc_interp_acc rather than rdcalcdp

      dimension Freq(*),psa(*),a(*)
      real pi,damp,Maxtim
!     damp is % critical damping. freq2 is maximum freq. to consider.

      pi = 4.0*atan(1.0)
!     Convert % damping to fraction
      damp_fraction=damp*0.01

!     Add 5 sec. to dur for shake-down
      maxtim=dur+5.0
      n=ifix(dur/dt)
      ntime=ifix(maxtim/dt)
      do j=n+1,ntime
          a(j)=0.
      enddo

!      Call new response spectrum routine for each desired freq
      do k=1,nfreq
          omega=2.*pi*freq(k)
          call rscalc_interp_acc(a, ntime, omega, damp_fraction, dt, 
     :                           sd, rv, aa)
          PSA(k)=sd*omega*omega
      enddo
      return
      end
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computeStochasticWave(
     :             seism,
     :             subFaultMoment, subfaultDistance,
     :             F0main, f0,  dur_sub,
     :             iflagscalefactor, scalingFactor)
!     calculates accelerogram from individual subfault by using
!     band-limited white noise method (Boore, BSSA, no.6/83).
!     An w-squared spectrum is assumed

! Dates: 08/20/08 - Added scale factor flag so can choose which factor
!                   to use without recompiling the program
!                   (iflagscalefactor = 1, 2, 3, 4 = vel^2, acc^2, acc, DMB respectively).
!        08/21/08 - Add scalingfactor (h) to output.
!                   Do/do not apply the scalefactor taper, depending on the
!                   value of taper_scalefactor.
!                   Doubled dimension of s
!        12/01/08 - Removed "NumberOfActiveSubs" from argument lists
!        12/17/08 - Added f0 and risetime to argument list and removed iii, jjj, z.
!                   Also, the program used dur rather than dur_sub in ndur, etc.
!                   I changed this.  I also pass dur_sub through the argument list
!                   rather than compute it here.
!        11/25/11 - The previous code used the fmax either as fmax or as kappa, which is
!                   very poor coding; I use one or the other or both, which will decrease
!                   the chance of mistakes.
!        05/15/15 - Removed nw, nl, risetime from argument list

      dimension seism(*)
      complex s(:)
      allocatable :: s
      include 'common_blocks.fi'
      
!      common /params/
!     :   rho,beta,prtitn,rtp,fs,
!     :   iKapa,fmax,
!     :   fr1, qr1, s1q, ft1, ft2, fr2, qr2, s2q, c_q, qt1, stq, 
!     :   iwind,n_crustal_amps,n_site_amps,totalMoment,
!     :   c4taper, fc4taper,
!     :   nsubs, flocut, nslope, pi, twopi
!      common /par/ iFFTsub, dt, npadl, npadt
!      common /seed/ iseed  ! Dave added this statement

!     specify point-source constants
        tmax=real(iFFTsub)*dt
        df=1.0/tmax
!        taper=0.02
        
!        prtitn=1.0/sqrt(2.0)
!        rtp=0.55
!        fs=2.0
 
!     define spectral constant, for subfaultDistance=1 km
        const=prtitn*rtp*fs*(1.e-20)/(4.*pi*rho*beta**3.)
        nnyq=iFFTsub/2+1
 
      if (dur_sub < 0.0) then
        print *,'STOPPING!!! Negative duration. Check duration model'
        stop
      end if
      
      ndur=dur_sub/dt
      ntaper=int(taper*ndur)
      nstop=ndur+ntaper+ntaper
      if(nstop>iFFTsub) then
          write(*,*) 'PROGRAM STOPS: subsource duration too long.'
          write(*,*) 'Options to fix the problem:'
          write(*,*) '1) Increase iFFTsub or dt, or'
          write(*,*) '2) Reduce distance to observation point, or'
          write(*,*) '3) Change parameters of duration model'
          stop
      endif

      allocate (s(iFFTsub))  ! iFFTsub is passed through the common block "par"

      do i=1,iFFTsub
          seism(i)=0.0
          s(i)=cmplx(0.0,0.0)
          seism(i)=0.0
      enddo

! generate Gaussian white noise with zero mean, unit variance. Noise is
! tapered.

      do i=1,nstop
          if (iwind == 0) then
            call window(i,1,nstop,ntaper,wind)
          else if (iwind == 1) then
            call wind2 (i,1,nstop,ntaper,dur_sub,0.2,0.2,wind)
          else
            print *,'STOPPING: invalid value of iwind (= ', iwind, ')'
            stop
          end if
          s(i+npadl)=wind*cmplx(gasdev(iseed),0.)
      enddo

! transform to frequency domain
      call fork (iFFTsub,s,-1.)

! find scaling factor, H
      averageMoment=totalMoment/real(nsubs)
      if (iflagscalefactor == 3) then
        ScalingFactor=sqrt(real(nsubs))*(F0main/f0)**2 !! David Boore
      else      
        sum1=0
        sum2=0
        do i=1,nnyq
          f=(i-1)*df
          if (akappa == 0.0) then
            highc=1./sqrt(1.+(f/fmax)**8.)
          else if (fmax == 0.0) then
            highc=exp(-pi*akappa*f)
          else
            highc=1./sqrt(1.+(f/fmax)**8.) * exp(-pi*akappa*f)
          end if
          if (iflagscalefactor == 1) then  ! vel spectrum
            s1=const*totalMoment*
     :           (2.*pi*f)**1.0*(1/(1.+(f/F0main)**2.))*
     :           highc
            s2=const*averageMoment*
     :           (2.*pi*f)**1.0*(1/(1.+(f/f0)**2.))*
     :           highc
          else if (iflagscalefactor == 2) then  ! acc spectrum
            s1=const*totalMoment*
     :           (2.*pi*f)**2.0*(1/(1.+(f/F0main)**2.))*
     :           highc
            s2=const*averageMoment*
     :           (2.*pi*f)**2.0*(1/(1.+(f/f0)**2.))*
     :           highc
          else
            print *,' ERROR, iflagscalefactor = ', iflagscalefactor
            print *, '; not a legal value; QUITTING!!!'
            stop
          end if
          sum1=sum1+s1**2/real(nsubs)
          sum2=sum2+s2**2
        enddo      
        scalingFactor=sqrt(sum1/sum2)
      end if

! Compute parameters used to taper spectrum toward low frequencies
      c4taper = sqrt(real(nsubs))/scalingFactor
      fc4taper = f0/sqrt(c4taper)
      
! multiply spectrum by an w-squared spectral shape after normalizing to
! square of unit spectral amplitude
      call avg (nnyq,s,avamp2)
      avamp=sqrt(avamp2)
      do i=1,nnyq
        f=(i-1)*df
        s(i)=sourceSpectra(f,subfaultDistance,const,
     :           scalingFactor,subFaultMoment,f0) *
     :                               s(i)/avamp
        if (i/=1) s(iFFTsub+2-i)=conjg(s(i))
      enddo
      s(nnyq)=cmplx(real(s(nnyq)),0.)

!     transform back to time domain
      call fork(iFFTsub,s,1.)
      afact=sqrt(real(iFFTsub))/tmax
      do i=1,iFFTsub
          seism(i)=afact*real(s(i))
      enddo

      deallocate (s)


      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine createRandomWeight(nl,nw,slipWeights)
      include 'common_blocks.fi'

      dimension slipWeights(1000,1000)
           do i=1,nl
              do j=1,nw
!**              draw random number 0 to 1.  This was changed from
!                using =ggnqf(iseed)+1. Also impose nonzero weight.
                 slipWeights(i,j)=Ran1(iseed)
                 if (slipWeights(i,j)<=0.) slipWeights(i,j)=0.001
              enddo
           enddo
       return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dur_path_cmp(subfaultDistance, 
     :                        rpathdur, pathdur, durslope, ndur_hinges,
     :                        dur_path)
!     computes duration dur due to path for a given subsource
! Dates: 12/05/09 - duration parameters contained in arrays and 
! duration computation done as in smsim (subroutine durpath in rv_td_subs)
      real rpathdur(*), pathdur(*)  
      
      include 'common_blocks.fi'
     
      if ( subfaultDistance <= rpathdur(1) ) then
        dur_path = pathdur(1)
      else if ( subfaultDistance >= rpathdur(ndur_hinges) ) then
        dur_path = pathdur(ndur_hinges) + 
     :     (subfaultDistance - rpathdur(ndur_hinges))*durslope
     :              
      else
        call locate(rpathdur, ndur_hinges, subfaultDistance, j)
        dur_path = pathdur(j) + 
     :     (subfaultDistance - rpathdur(j))*( pathdur(j+1)-pathdur(j))
     :              /(rpathdur(j+1)-rpathdur(j))
      end if

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine findDistanceAndAzimuth (FaultLat,FaultLon,SiteLat,
     *                                 SiteLon,epi,azi,isitecoordflag)
!       calculates distance and azimuth between two
!       points using their latitudes and longitudes

! azi is in degrees

! Dates: 02/17/09 - renamed some variables to avoid possible confusion
!                   between angles in degrees and radians.

! The input and output angles are in degrees

        pi = 4.0*atan(1.0)
        d2r = pi/180.0
        
        if (isitecoordflag == 1) then  ! lat,long

          re=6371.
          alat=(FaultLat+SiteLat)/2.
          alat_radians=d2r*alat
          dlat=SiteLat-FaultLat  ! degrees
          dlon=SiteLon-FaultLon  ! degrees
          r1 = d2r*re*(SiteLat-FaultLat)
          epi= d2r*re*sqrt(dlat**2.+cos(alat_radians)**2.*dlon**2.)
          if(r1>=epi)then
             azi_radians=pi
          else
             azi_radians=acos(r1/epi)
          endif

          if (dlon<=0.) azi_radians=2.0*pi-azi_radians
          azi= azi_radians/d2r

        else if (isitecoordflag == 2) then  ! R,Az
        
          epi = sitelat
          azi = sitelon
          
        else  ! assume cartesian coordinates
        
          xn = sitelat
          xe = sitelon
          
          epi = sqrt(xn**2 + xe**2)
          
          azi = 180.0*atan2(xe, xn)/pi
          if (azi < 0.0) then
            azi = 360.0 + azi
          end if
          
        end if
          
          

        return
        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function findPeakValue(nptsTotalWave,totalWave)
!     find absolute maximum acceleration in simulated time history
      dimension totalWave(*)
        peakValue=0.
        do i=1,nptsTotalWave
          if(abs(totalWave(i))>peakValue) peakValue=abs(totalWave(i))
        enddo
        findPeakValue=peakValue
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real function findSubfaultDistance (R,h,FaultStrike,
     *                  fi2,FaultDip,dl,dw,i,j)
!       computes distance subfaultDistance from center of subfault (i,j)
!       to observation point using formula (1) (Figure 1)

! Note confusion in what dip means in the figure 1 of Beresnev and Atkinson.
! What they label as "delta" is what was used in the original version of
! this subroutine as "dip".  Their figure has delta_1 as the true 
! fault dip, probably added as a result of a review. I have
! replaced "dip" by 90_faultdip here to eliminate confusion (D. Boore, 17Feb09).

        pi = 4.0*atan(1.0)
        d2r = pi/180.0
        
        a90_faultdip_radians = d2r*(90.0-FaultDip)
        phi2_strike_radians = d2r*(fi2-FaultStrike)
        t1=R*cos(phi2_strike_radians)-(2.*i-1)*dl/2.
        t2=R*sin(phi2_strike_radians)-(2.*j-1)*dw/2.*
     :       sin(a90_faultdip_radians)
        t3=-h-(2.*j-1)*dw/2.*cos(a90_faultdip_radians)
        findSubfaultDistance=sqrt(t1**2.+t2**2.+t3**2.)
        return
        end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findSubFaultMoment(nl,nw,slipWeights,totalMoment,
     *                             weightedMoment)
      dimension slipWeights(1000,1000)
      dimension weightedMoment(1000,1000)

!       calculate moments of subfaults. First calculate the total of all weights
!       (totalSlipWeights), then compute the total moment per unit 
!       weight (elem), and use this to compute the moment for each subfault.
 
      totalSlipWeights=0.
      do i=1,nl
        do j=1,nw
          totalSlipWeights=totalSlipWeights+slipWeights(i,j)
        enddo
      enddo

      elem=totalMoment/totalSlipWeights
      do i=1,nl
        do j=1,nw
          weightedMoment(i,j)=slipWeights(i,j)*elem
        enddo
      enddo

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getAmps(nu_getamps, 
     :   f_crustal_amps,
     :     n_crustal_amps, freq_crustal_amps, amp_crustal_amps, 
     :   f_site_amps,
     :     n_site_amps,    freq_site_amps,    amp_site_amps) 

! Dates: 11/30/11 - Written by D. Boore (reading of the crustal and site amps
!                   used to be within the sourceSpectra, which was a poor
!                   choice by the people who wrote the original version of exsim).

      character f_crustal_amps*(*), f_site_amps*(*)
      real freq_crustal_amps(*), amp_crustal_amps(*), 
     :        freq_site_amps(*), amp_site_amps(*)
      integer n_crustal_amps, n_site_amps
      
      call trim_c(f_crustal_amps, nc_f_crustal_amps)
      open (nu_getamps,
     :  file=f_crustal_amps(1:nc_f_crustal_amps),status='unknown')
     
      read(nu_getamps,*) n_crustal_amps
      do i=1,n_crustal_amps
        read (nu_getamps,*) freq_crustal_amps(i), amp_crustal_amps(i)
      enddo
      close (nu_getamps)
        
      call trim_c(f_site_amps, nc_f_site_amps)
      open (nu_getamps,
     :  file=f_site_amps(1:nc_f_site_amps),status='unknown')
     
      read(nu_getamps,*) n_site_amps
      do i=1,n_site_amps
        read (nu_getamps,*) freq_site_amps(i), amp_site_amps(i)
      enddo
      close (nu_getamps)
        
       
      return
      end



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine getInputParameters(nu_read_params,
     :  rho,beta,prtitn,rtp,fs,r_ref, nsprd_segs, rlow, a_s, b_s, m_s,
     :  fr1, qr1, s1q, ft1, ft2, fr2, qr2, s2q, c_q,
     :  qt1, stq, ndur_hinges,rpathdur, 
     :  pathdur, durslope,fmax, akappa_0, dkappadmag, amagkref,
     :  akappa,flocut, nslope, iwind, taper, eps_w, eta_w, f_tb2te, 
     :  f_te_xtnd, tpadl, tpadt, dt, iseed, nsims, amag, stress,
     :  FaultLat,FaultLon, FaultStrike,FaultDip,h, fault_type,
     :  FaultLength, FaultWidth, dl, dw, stress_ref,nl, nw, nsubs, 
     :  specify_length, specify_width, vrup_beta, hyp_loc_fl, 
     :  hyp_loc_fw, n_hypocenters, i0_in, j0_in, i_rise_time,
     :  iDynamicFlag, pulsingPercent, iflagscalefactor, islipweight, 
     :  slipWeights, iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :  f_crustal_amps, f_site_amps, iflagfas_avg, 
     :  iflagpsa_avg_over_sims, iflagpsa_avg_over_hypos,
     :  write_site_files, f_stem, damp,
     :  nfreq,freq1,freq2, nfout, freq_out,
     :  isitecoordflag, move_site, numberOfSites, siteLocation)

! Dates: 11/30/08 - Move computation of random slipweights 
!                   if Islipweight == 1.0 from main program to here
!                   (a more logical place).
!        12/02/08 - Get i_rise_time, to determine what type of risetime is used 
!        11/22/11 - A complete revision of the params file.
!        05/13/15 - Some small cosmetic changes (e.g., replace ".ge." with ">=")
!        05/15/15 - Removed jColumn, iRow from argument list;
!                   reduce length of row containing the reading of
!                   f_te_xtnd (note: f_te_xtnd is not used in this version of exsim_dmb)
!        05/17/15 - Compute qt1, st here.
     
      real rlow(*), a_s(*), b_s(*), m_s(*), freq_out(*)
      real rpathdur(*), pathdur(*)
      dimension slipWeights(1000,1000), siteLocation(3000,2)
      character f_stem*(*), 
     :  f_crustal_amps*30,f_site_amps*30,aline*60, fault_type*(*),
     :  version_ctl*8, version_in*30
      character cmnts2skip(50)*80, buf_in*10, f_slip_weights*80
      character buf*200
      
      logical f_exist, specify_length, specify_width,
     :        move_site, write_site_files
      
      pi = 4.0*atan(1.0)
      d2r = pi/180.0

      version_ctl = ' '
      version_ctl = '11/22/11'
      call trim_c(version_ctl,nc_version_ctl)

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      version_in = ' '
      read(nu_Read_Params,'(a)') version_in
      call trim_c(version_in,nc_version_in)

      if (version_ctl(1:nc_version_ctl) /= 
     :    version_in(1:nc_version_in)) then
        write(*,*) 
        write(*,'(2x,a)') 
     :      'STOPPING because the parameter file has the '//
     :      'wrong version number.'
        write(*,'(2x,a)') 
     :      'It is likely that some parameters and/or the '//
     :      'order in which they'
        write(*,'(2x,a)') 
     :      'are read have changed.  Please change your '//
     :      'file to match the current' 
        write(*,'(2x,a)') 
     :      'parameter file.'
        close(nu_Read_Params)
        stop
      end if
      
!
!-------------------------------------------------------------------------
! ******* Input parameters common to SMSIM and EXSIM (in the order in which 
!         they appear in the SMSIM parameter file) *******
!-------------------------------------------------------------------------
!
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      aline = ' '
      read (nu_Read_Params,'(a)') aline               ! Title
      call trim_c(aline, nc_aline)
!      print *, ' Title: '//aline(1:nc_aline)
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) rho, beta, prtitn, rtp, fs

! gsprd: 
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params, *) r_ref
      read(nu_Read_Params, *) nsprd_segs
      do i = 1, nsprd_segs
        read(nu_Read_Params, *) rlow(i), a_s(i), b_s(i), m_s(i)
      end do
 
! Q:
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params, *) 
     :   fr1, qr1, s1q, ft1, ft2, fr2, qr2, s2q, c_q
      qt1 = qr1*(ft1/fr1)**s1q
      qt2 = qr2*(ft2/fr2)**s2q
      stq = 0.0
      if (ft1 .ne. ft2) then
        stq = alog10(qt2/qt1)/alog10(ft2/ft1)
      end if
       
! Duration:
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) ndur_hinges
      do i = 1, ndur_hinges
        read(nu_Read_Params,*) rpathdur(i), pathdur(i)
      end do
      read(nu_Read_Params,*) durslope
 
! Site diminution:
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) fmax, akappa_0, dkappadmag, amagkref

! Low-cut filter:
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) flocut, nslope

! Window params:
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) iwind, taper, eps_w, eta_w, 
     :  f_tb2te, f_te_xtnd
      
! Time series params:      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) tpadl, tpadt, dt, seed, nsims   
!      read (nu_Read_Params,*) dur_fctr, dt, tshift, seed, nsims                 
      iseed = -int(abs(seed))
! NOTE: The random number generation uses either the 
! Numerical Recipes subroutines ran1 (for uniform distributions) or gasdev 
! (for normal distributions); these routines require a negative integer seed for 
! the first call to the subroutine; for this reason the following statement is used here
! iseed = -iabs(iseed).
! This statement used to be in the main body of the program, but it needs to be before 
! the first call to either ran1 or gasdev, which occurs here if random slip weights
! are requested (the call is to createRandomWeight).
! I am reading the input parameter seed as a real number, in keeping with the SMSIM params file.
! (See 11/25/08 revision comments for more information).

!
!-------------------------------------------------------------------------
! ******* Input parameters specific to EXSIM *******
!-------------------------------------------------------------------------
!

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) amag, stress

      if (akappa_0 /= 0.0) then
        akappa = akappa_0 + dkappadmag*(amag-amagkref)
      else
        akappa = 0.0
      end if
      


      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) FaultLat, FaultLon
     
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) FaultStrike, FaultDip, h
     
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      fault_type = ' '
      read(nu_Read_Params,'(a)') fault_type
      call trim_c(fault_type, nc_fault_type)
      call upstr(fault_type(1:1))
     
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) FaultLength, FaultWidth, dl, dw, 
     :   stress_ref
      
      stress_factor = (stress_ref/stress)**(1.0/3.0)
 
      specify_length = .true.
      specify_width  = .true.


      if (faultlength == 0.0) then
        specify_length = .false.
        if(fault_type(1:1) == 'S') then
          faultlength = 10.0**(-2.57+0.62*amag)
        else if(fault_type(1:1) == 'R') then
          faultlength = 10.0**(-2.42+0.58*amag)
        else if(fault_type(1:1) == 'N') then
          faultlength = 10.0**(-1.88+0.50*amag)
        else
          faultlength = 10.0**(-2.44+0.59*amag)
        end if
        
!        print *, ' WC FL = ', faultlength
        faultlength = stress_factor * faultlength
!        print *, ' WC FL, after scaling = ', faultlength
       
      end if
      if (faultwidth == 0.0) then
        specify_width  = .false.
        if(fault_type(1:1) == 'S') then
          faultwidth = 10.0**(-0.76+0.27*amag)
        else if(fault_type(1:1) == 'R') then
          faultwidth = 10.0**(-1.61+0.41*amag)
        else if(fault_type(1:1) == 'N') then
          faultwidth = 10.0**(-1.14+0.35*amag)
        else
          faultwidth = 10.0**(-1.01+0.32*amag)
        end if
        
!        print *, ' WC FW = ', faultwidth
        faultwidth = stress_factor * faultwidth
!        print *, ' WC FW, after scaling = ', faultwidth
       
      end if
      
      nl=FaultLength/dl
      if (nl < 1) nl = 1
      nw=FaultWidth/dw
      if (nw < 1) nw = 1
 
      nsubs = nl * nw
        
      dl = FaultLength/real(nl)  ! Need to reset dl and dw because nl and nw are integers      
      dw = FaultWidth/real(nw)        
c BE 06/12/16 all 200,200 array definitions changed to 1000,1000 and the following updated:
c      if(nl > 200 .or. nw > 200) then
      if(nl > 1000 .or. nw > 1000) then
        print *,"Burp!, you exceeded the
     :         the array dimentions, (1000,1000)"
        stop
      end if
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) vrup_beta
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) hyp_loc_fl, hyp_loc_fw, n_hypocenters
      if (nl == 1 .and. nw == 1) then   ! a single subsource
        i0_in = 1
        j0_in = 1
        n_hypocenters = 1
      else if (hyp_loc_fl >= 0.0 .and. hyp_loc_fw >= 0.0 ) then ! nonrandom hypocenter 
                                                                ! convert hyp_loc to subfault index
                                                                ! NOTE on 13May15: hyp_loc_fl > 0.0 in previous version
        i0_in = int(hyp_loc_fl/dl)+1
        j0_in = int(hyp_loc_fw/dw)+1
        n_hypocenters = 1       ! not a random hypocenter, force this
      else ! not a single subsource and random
        i0_in = 0
        j0_in = 0
      end if
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) i_rise_time

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*)iDynamicFlag,pulsingPercent

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) iflagscalefactor ! (1=vel^2; 2=acc^2; 3=asymptotic acc^2 (dmb))
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) islipweight
     
      f_slip_weights = ' '
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,'(a)') f_slip_weights                ! need to read a dummy string
                                                        ! even if slip weights are not
                                                        ! specified
      call trim_c(f_slip_weights, nc_f_slip_weights)
      

      if (islipweight == -1) then   ! assign unity to all cells
        do i = 1, nl
          do j = 1, nw
            slipWeights(i,j) = 1.0
          end do
        end do
      else if (islipweight == 0) then
        call get_lun(nu_slip_weights)
        open(nu_slip_weights,file=f_slip_weights(1:nc_f_slip_weights),
     :       status='unknown')
        read(nu_slip_weights,*)((slipWeights(i,j),i=1,nl),j=1,nw)
        close(nu_slip_weights)
!        write(*,'(/a)')"       ***  read slip distribution   ***"
      else if (islipweight == 1) then
        call createRandomWeight(nl,nw,slipWeights)
      else
        print *,' islipweight = ', islipweight
        print *, 'NOT A LEGAL VALUE; STOPPING'
        stop
      endif

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak
     
!
!-------------------------------------------------------------------------
! PARAMETERS RELATED TO PATH AND SITE:
!-------------------------------------------------------------------------
!

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      f_exist = .false.
      f_crustal_amps = ' '
      read (nu_Read_Params,'(a)') f_crustal_amps
      call trim_c(f_crustal_amps,nc_f_crustal_amps)
      inquire(file=f_crustal_amps(1:nc_f_crustal_amps), exist=f_exist)
      if (.not. f_exist) then
        print *,' file '//f_crustal_amps(1:nc_f_crustal_amps)//
     :          ' does not exist; QUITTING!!!'
        stop
      end if

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      f_exist = .false.
      f_site_amps = ' '
      read (nu_Read_Params,'(a)') f_site_amps
      call trim_c(f_site_amps,nc_f_site_amps)
      inquire(file=f_site_amps(1:nc_f_site_amps), exist=f_exist)
      if (.not. f_exist) then
        print *,' file '//f_site_amps(1:nc_f_site_amps)//
     :          ' does not exist; QUITTING!!!'
        stop
      end if


!
!-------------------------------------------------------------------------
! PARAMETERS RELATED TO COMPUTATIONS OF AVERAGES:
!-------------------------------------------------------------------------
!

! Get types of averages for FAS and PSA
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) iflagfas_avg
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) iflagpsa_avg_over_sims

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) iflagpsa_avg_over_hypos

!
!-------------------------------------------------------------------------
! PARAMETERS RELATED TO THE OUTPUT:
!-------------------------------------------------------------------------
!
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      buf_in = ' '
      read (nu_Read_Params,'(a)') buf_in
      call upstr(buf_in)
      call trim_c(buf_in, nc_buf_in)
      if (buf_in(1:1) == 'Y') then
         write_site_files = .true.
      else
         write_site_files = .false.
      end if
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      f_stem = ' '
      read (nu_Read_Params,'(a)') f_stem
      call trim_c(f_stem, nc_f_stem)
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) damp

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) nfreq,freq1,freq2
 
      nfreq_max = 500

      if (iabs(nfreq) > nfreq_max) then
           print *,' '
           print *, ' nfreq = ', nfreq, 
     :              ' but it cannot exceed ', nfreq_max
           print *,' STOPPING!!!'
           stop
      endif

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)      
      read(nu_Read_Params,*) nfout
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read(nu_Read_Params,*) (freq_out(i),i=1,nfout)
      
!
!-------------------------------------------------------------------------
! PARAMETERS RELATED TO THE SITES AT WHICH MOTIONS ARE COMPUTED:
!-------------------------------------------------------------------------
!

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      read (nu_Read_Params,*) isitecoordflag
      
      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      buf_in = ' '
      read (nu_Read_Params,'(a)') buf_in
      call upstr(buf_in)
      call trim_c(buf_in, nc_buf_in)
      if (buf_in(1:1) == 'Y') then
         move_site = .true.
      else
         move_site = .false.
      end if
      
      if (move_site .and. FaultStrike /= 0.0) then
        print *,' ERROR: Cannot request move site if'//
     :          ' the fault strike /= 0.0; QUITTING'
        stop
      end if

      call skipcmnt(nu_Read_Params,cmnts2skip,nc_cmnts2skip)
      
      numberOfSites = 0
      DO
        buf = ' '
        read (nu_Read_Params,'(a)') buf
        call trim_c(buf,nc_buf)
        call upstr(buf)
        if (buf(1:4) == 'STOP') exit
        numberOfSites = numberOfSites + 1
        read (buf(1:nc_buf),*) 
     :    siteLocation(numberOfSites,1), siteLocation(numberOfSites,2)        
      END DO
 

      isitecoordflag_in = isitecoordflag
      Do i=1,numberOfSites
        if (move_site .and. isitecoordflag_in == 1 ) then
          dum = 0.0
          ! eventually might do something with this case
        else if (move_site .and. isitecoordflag_in == 2 ) then 
          isitecoordflag = 3
          r = siteLocation(i,1)
          az = siteLocation(i,2)
          xn = r*cos(d2r*az)
          xe = r*sin(d2r*az)
          siteLocation(i,1) = faultLength/2.0 + xn
          siteLocation(i,2) = xe
        else if (move_site .and. isitecoordflag_in == 3) then
          if (sitelocation(i,1) == 0.0) then  ! move to midpoint
            siteLocation(i,1) = faultLength/2.0
          end if
          if (siteLocation(i,2) == 0.0) then  ! move to end
            siteLocation(i,1) = faultLength + siteLocation(i,1)
          end if
        end if
      enddo

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ggnqf(iseed)
!     generates Gaussian white noise by approx. method of Ross (1985)
      usum=0.
      do j=1,12
          usum=usum+Ran1(iseed)
      enddo
      ggnqf=usum-6.
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine locateHypocentreRandomly(i0,j0,nl,nw)
      include 'common_blocks.fi'

!**    draw random hypocentre
      integer i0,j0,nl,nw

      i0= nint(Ran1(iseed)*nl)
      if (i0 < 1) i0=1
      if (i0 > nl) i0=nl

      j0=nint(Ran1(iseed)*nw)
      if (j0 < 1) j0=1
      if (j0 > nw) j0=nw

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makeItComplex(x,cx,noOfSegPoints)
      dimension x(*)
      complex cx(*)
         do jjj = 1, noOfSegPoints
             cx(jjj) = cmplx(x(jjj),0.)
         enddo
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine MavroPapa(totalWave,nptsTotalWave,
     *            PapaGama,PapaNu,PapaT0,amag,PapaPeak)
!c     include Mavroeidis and Papageorgiou, 2003 approach

! Dates: 05/15/15 - Remove iPapaFlag from argument list;
!                   delete papaW, cpapaW 

      include 'common_blocks.fi'
      dimension totalWave(*),PapaPulseW(70000)
      complex   cTotalWave(70000),cPapaPulseW(70000) 
      complex   cTemp(70000)
      integer nptsTotalWave
      pi=4.0*atan(1.0)
      peakValue=findPeakValue(nptsTotalWave,totalWave)
      papaA=PapaPeak
      papaTp=10**(-2.9+0.5*aMag)
      papaFp=1.0/papaTp
      Fp=papaFp
      samplingRate=1./dt
      initialnptsTotalWave=nptsTotalWave

      T0=papaT0

!cccccc step1 calculate pusle wave in the time domain
      PapaPulseW=0.0
      do i=1,nptsTotalWave
        t=(i-1)*dt
          if(t>=(T0-papaGama/2./Fp).and.t<=(T0+papaGama/2./Fp))then
            c1=papaA*pi*Fp/papaGama
            c2=sin(2*pi*Fp*(t-T0)/papaGama)*cos(2*pi*Fp*(t-T0)+papaNu)
            c3=papaGama*sin(2*pi*Fp*(t-T0)+papaNu)
            c4=1+cos(2*pi*Fp*(t-T0)/papaGama)
            PapaPulseW(i)=c1*(c2+c3*c4)
        else
            PapaPulseW(i)=0.0
        endif
      enddo

!cccccc step3-a take totalWave (stochastic) to the frequency domain
      call padding(TotalWave,nptsTotalWave)
      call makeItComplex(totalWave,cTotalWave,nptsTotalWave)
      call fork(nptsTotalWave,cTotalWave,-1.)


!cccccc step3-b take PapaPulseW to the frequency domain
      call padding(PapaPulseW,nptsTotalWave)
      call makeItComplex(PapaPulseW,cPapaPulseW,nptsTotalWave)
      call fork(nptsTotalWave,cPapaPulseW,-1.)

!cccccc step 4, keep the phase of stocastic,ccalculate abs(cTotalWave)-abs(cPapaPulseW)
!cccccc calculate real and imaginary part  of the resulus
      do i=1,nptsTotalWave
        if(real(cTotalWave(i))>=0)then !!very importat when using atan
             phase=atan(aimag(cTotalWave(i))/real(cTotalWave(i)))
        else !!very importat when using atan
             phase=atan(aimag(cTotalWave(i))/real(cTotalWave(i)))+pi
        endif
            amp=abs(abs(cTotalWave(i))-abs(cPapaPulseW(i)))
          x=amp*cos(phase)
          y=amp*sin(phase)
          cTemp(i)=cmplx(x,y)
      enddo
!cccccc take care of complex conjugates
      cTemp(nptsTotalWave/2+1)=cmplx(1.,0.)*
     *                               cTemp(nptsTotalWave/2+1)
      do i=1, nptsTotalWave/2+1
                cTemp(nptsTotalWave+2-i) = conjg(cTemp(i))
      enddo
!cccccc go back to the time domain
      call fork(nptsTotalWave,cTemp,+1.)
 
!cccccc add them up
      do i=1,nptsTotalWave
         totalWave(i)=real(cTemp(i))+PapaPulseW(i)
      enddo
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function NumberOfPulsingSubs(i0,j0,i,j,nl,nw
     *                                ,NoeffectiveSubfaults)
!**   this function determines how many subfaults are active at this time (i,j).

! Dates: 12/01/08 - DMB added integer specification for Rmin, Rmax
!                   and used iabs rather than abs.

      integer Rmin, Rmax, r,NumberOfPulsingSubs
      n=0
      RMax=max(iabs(i-i0)+1,iabs(j-j0)+1)
      RMin=Rmax-NoeffectiveSubfaults
      if(RMin<0)Rmin=0

      do ii=1,nl
         do jj=1,nw
            r=max(iabs(ii-i0)+1,iabs(jj-j0)+1)
            if(r>RMin.and.r<=RMax)n=n+1
         enddo
      enddo
      NumberOfPulsingSubs=n

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine padding(x,noOfSegPoints)
      dimension x(70000)

      do ii=1,100
        ncheck=2**(ii)
        if (ncheck>=noOfSegPoints) go to 110
        if(ncheck > 70000) then
           write(*,*)' Burp !! Exceeds array dimensions.'
           stop
        endif
      enddo
110   continue

      do j = noOfSegPoints+1, ncheck
            x(j)=0.
      enddo
      noOfSegPoints=ncheck

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sample(x,npts,df,nfreq,freq,FA)
      dimension x(*), freq(500), FA(500)

!**   takes a large FFT array and samples for nfreq points,
!       in array FA.  Uses avg (sq) FA of 5 around each freq.

 
      FA=0.0
      
      do jf = 1, nfreq
      
        itarg= ifix(freq(jf)/df)+1
        
!        f=df*real(itarg-1)        ! dmb replaced itarg with itarg-1
                                   ! 05/15/15: commented out, as f is never used
        
        if (itarg < 4) then
        
            FA(jf)=-9.999   ! dmb
            
!            FA(jf)=9.999

        else
        
!             find spectrum around itarg
          sum = 0.0            ! dmb
          do ii=itarg-2, itarg+2
              sum=0.2*x(ii)*x(ii) + sum   ! dmb
          enddo
          if(sum==0) then
            FA(jf)= -9.999        ! dmb
          else
            FA(jf) = sqrt(sum)  ! dmb
          end if
        endif                              ! dmb
        
!        if(FA(jf)==0)FA(jf)=9.999
!            FA(jf) = alog10(sqrt(FA(jf)))
!        endif

      enddo

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine site_coords_in_degrees(
     :         FaultLat, FaultLon, 
     :         SiteLat, SiteLon, isitecoordflag, 
     :         site_lat_degrees, site_lon_degrees ) 
     
        pi = 4.0*atan(1.0)
        dtor = pi/180.0
     
! If site coords not in degrees, convert      
        if (isitecoordflag == 1) then  ! degrees
          site_lat_degrees = sitelat
          site_lon_degrees = sitelon         
        else if (isitecoordflag == 2) then  ! R,Az
          r = sitelat
          az = sitelon
          xn = r*cos(dtor*az)
          xe = r*sin(dtor*az)
          call km2deg_f( xn, xe, FaultLat, FaultLon,  
     :            site_lat_degrees, site_lon_degrees )          
        else  ! assume cartesian coordinates
          xn = sitelat
          xe = sitelon
          call km2deg_f( xn, xe, FaultLat, FaultLon,  
     :            site_lat_degrees, site_lon_degrees )          
        end if          
          
          
      
        return
        end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function sourceSpectra (f,subfaultDistance,
     :           const,scalingFactor,aMoment,cornerF)
!     Computes amplitude value of w-squared spectrum at
!     distance subfaultDistance.
!     Source terms correspond to subfaultDistance=1 km

! Dates: 02/07/09 - Determine gspread and q using gsprd_f and q_f (these
!                   need to be initialized and deallocated in the main program).
!        11/30/11 - Read site and crustal amps in getAmps
!                 - Remove check of namps /= 0; if no amps are desired,
!                   the amps files should contain one frequency with amplitude = 1.
!        05/17/15 - Call q_func, which is appended to the bottom of this file.

      include 'common_blocks.fi'
      w=twopi*f

      scalingFactor_hf_taper= scalingFactor * 
     :     c4taper*(1.0+(f/cornerF)**2.0)/(1.0+(f/fc4taper)**2.0)
 
      sa=const*aMoment*scalingFactor_hf_taper*w**2.0*
     :                   (1.0/(1.0+(f/cornerF)**2.0))
     
      if (f == 0.0) gamma=0.
      if (f /= 0.0) then
          Q= q_func(f, qr1,fr1,s1q, qt1,ft1,stq, ft2, qr2,fr2,s2q ) 
c edwards edit 171016          gamma=w/(2.*Q*beta)
          gamma=w/(2.*Q*c_q)
      endif
      anelas=exp(-gamma*subfaultDistance)
      if (akappa == 0.0) then
        highc=1./sqrt(1.+(f/fmax)**8.)
      else if (fmax == 0.0) then
        highc=exp(-pi*akappa*f)
      else
        highc=1./sqrt(1.+(f/fmax)**8.) * exp(-pi*akappa*f)
      end if
     
      ga = gsprd_f(subfaultdistance)
     
      am1=amplf(f, 
     :    freq_crustal_amps, amp_crustal_amps, n_crustal_amps)
      am2=amplf(f,
     :    freq_site_amps,    amp_site_amps,    n_site_amps)
      sourceSpectra=sa*ga*highc*anelas*am1*am2 * 
     :              buttrlcf( f, flocut, nslope)
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine window (i,nstart,nstop,ntaper,wind)
!     applies cosine tapered window
!     unit amplitude assumed
!     Subroutine from Dave Boore.
      wind=0.
      if (i<nstart.or.i>nstop) return
      wind=1.
      if(i>=nstart+ntaper.and.i<=nstop-ntaper) return
      pi=3.141593
!     Value of wind goes from 0 to 1 from nstart to nstart+ntaper,
!     then from 1 to 0 from nstop-ntaper to nstop
      dum1=(nstop+nstart)/2.
      dum2=(nstop-nstart-ntaper)/2.
      wind=0.5*(1.-sin(pi*(abs(real(i)-dum1)-dum2)/real(ntaper)))
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wind2 (i,nstart,nstop,ntaper,tw,eps,eta,wind)
!     applies Saragoni and Hart (1974) window, with parameters
!     tw (dur), eps (fraction of dur to reach peak), and
!     eta (fraction of peak ampl. at tw)
!     The Saragoni and Hart window is applied between
!     nstart + ntaper, and nstop - ntaper.  Outside these
!     bounds we do a quick cosine taper down to 0
!
      pi=3.141592654
      wind=0.
      if(i<nstart.or.i>nstop) return
      wind=1.

!     Apply Sargoni and Hart window

        b=-eps*alog(eta)/(1.+eps*(alog(eps)-1.))
        c=b/(eps*tw)
        a=(2.7182818/(eps*tw))**b
        if (i>=(nstart+ntaper).and.i<=(nstop-ntaper)) then
           t=(tw/real(nstop-nstart-2*ntaper)) *
     *       (real(i-nstart-ntaper+1))
           wind=a*t**b*exp(-c*t)
           return
        else
!     Cos. taper goes from  0 to 1 from nstart to nstart+ntaper, from
!                           1 to 0 from nstop-ntaper to nstop.

        if (i<(nstart+ntaper)) then
           t1=tw/real(nstop-nstart-2*ntaper)
           wf=a*t1**b*exp(-c*t1)
           wind=abs(sin((real(i-nstart)/real(ntaper))*pi/2.))*wf
        else
           wf=a*tw**b*exp(-c*tw)
           wind=abs(sin((real(nstop-i)/real(ntaper))*pi/2.))*wf
        endif
        return
      endif
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeAcc(ioAcc, totalWave, npts, peakValue,
     :  isim, isite, HypoDistance, faultDistance, fpar)
!       writes acceleration time history into specified ascii file
      include 'common_blocks.fi'
      character fpar*(*)
      dimension totalWave(*)
      integer  npts

      write (ioAcc,"('**********************************************')")
      write (ioAcc,"('*******  Acceleration Time Series ************')")
      write (ioAcc,"('*** Site #',i4)") isite
      write (ioAcc,"('*** trial ',i4)") isim
      write (ioAcc,"('Input Parameters file =  ',a)") fPar

      write  (ioAcc,50) npts
50    format (i6,' samples')
      write  (ioAcc,30) dt
30    format ('dt: ',f6.3,' sec')
      write  (ioAcc,'(a)') 'data format: (1x,f8.3,1p,2x,e11.4)'
      write  (ioAcc,10) peakValue
10    format ('maximum absolute acceleration:',f7.2)
      write(ioAcc,'("fault Dis.(km)=",f8.2)')faultDistance
      write(ioAcc,'("Hypo  Dis.(km)=",f8.2)')HypoDistance
      write  (ioAcc,'(2x,a, 1x,a)') 'time(s)', 'acc(cm/s**2)'

      do i=1,npts
         time=(i-1)*dt
         write (ioAcc,20) time,totalWave(i)
      enddo
20    format (1x,f8.3,1p,2x,e11.4)
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writeHUS( ioHus, husid, npts, 
     :                 arias, dur_75_05, 
     :                 isim,isite,HypoDistance,faultDistance,fpar)
!       writes husid time history, arias intensity, and 75%-5% duration
!       into specified ascii file
      include 'common_blocks.fi'
      character fpar*(*)
      dimension husid(*)
      integer  npts

      write (iohus,"('**********************************************')")
      write (iohus,"('*******  Husid Time Series ************')")
      write (iohus,"('*** Site #',i4)") isite
      write (iohus,"('*** trial ',i4)") isim
      write (iohus,"('Input Parameters file =  ',a)") fPar

      write  (iohus,50) npts
50    format (i6,' samples')
      write  (iohus,30) dt
30    format ('dt: ',f6.3,' sec')
      write  (iohus,'(a)') 'data format: (1x,f8.3,1p,2x,e11.4)'
      write  (iohus,'(a,1x,1pe10.3)') ' Arias intensity = ', arias
      write  (iohus,'(a,1x,1pe10.3)') ' 75%-5% duration = ', dur_75_05
      write(iohus,'("fault Dis.(km)=",f8.2)')faultDistance
      write(iohus,'("Hypo  Dis.(km)=",f8.2)')HypoDistance
      write  (iohus,'(5x,a, 7x,a)') 'time', 'husid'

      do i=1,npts
         time=(i-1)*dt
         write (iohus,'(1x,f8.3,1p,1x,e11.4)') time,husid(i)
      enddo

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writePar(nu_Write_Params,
     :   FaultStrike,FaultDip,h, FaultLat, FaultLon, 
     :   siteLocation,numberOfSites, move_site,
     :   fault_type, FaultLength, FaultWidth,nl,nw,dl,dw, 
     :   specify_length, specify_width, stress_ref,
     :   vrup_beta,
     :   hyp_loc_fl, hyp_loc_fw, i0_in,j0_in,n_hypocenters,
     :   nsims,
     :   amag,
     :   stress, 
     :   pulsingPercent,iDynamicFlag,
     :   i_rise_time,
     :   iPapaFlag,PapaGama,PapaNu,PapaT0,PapaPeak,
     :   iflagscalefactor,   
     :   iflagfas_avg,
     :   iflagpsa_avg_over_hypos, iflagpsa_avg_over_sims,
     :   tpadl, tpadt,
     :   r_ref, nsprd_segs, rlow, a_s, b_s, m_s,
     :   rpathdur, pathdur, durslope, ndur_hinges 
     :                    )
 
!     writes modeling parameters to specified ascii file

! Dates: 12/01/08 - DMB added writing of istart, istop, nptsenvelope, maxPointsOfWave
!        01/20/10 - Changed format descriptor of dt, kappa to f8.4; removed iFFTsub
!                   from output, because it has not been defined at the time that
!                   this subprogram is called.
!        11/25/11 - Some changes as a result of changing the params file
!        05/15/15 - Remove typ--it was set but not used

      logical specify_length, specify_width, move_site
      
      character fault_type*(*)
      
      dimension siteLocation(3000,2)


      real amag, r_ref, rlow(*), a_s(*), b_s(*), m_s(*)
      real rpathdur(*), pathdur(*)
      
      include 'common_blocks.fi'
      
      write(nu_Write_Params,'("         modeling parameters     ")' )
      write(nu_Write_Params,'("Fault Strike              = ",f8.2)')
     :   FaultStrike
      write(nu_Write_Params,'("Fault dip                 = ",f8.2)')
     :   FaultDip
      write(nu_Write_Params,'("Fault depth to upper edge = ",f8.2)')h
      
      if(.not. specify_length) then
        write(nu_Write_Params,'(a,1x, f6.1)') 
     :    ' Fault length from Wells and Coppersmith for fault type '//
     :      fault_type(1:1)//', using a reference stress of ', 
     :                             stress_ref
      end if     
      write(nu_Write_Params,'("Fault Length              = ",f8.2)')
     :   FaultLength
      
      if(.not. specify_width) then
        write(nu_Write_Params,'(a,1x, f6.1)') 
     :    ' Fault width from Wells and Coppersmith for fault type '//
     :      fault_type(1:1)//', using a reference stress of ', 
     :                             stress_ref
      end if     
      write(nu_Write_Params,'("Fault Width               = ",f8.2)')
     :   FaultWidth
      
      write(nu_Write_Params,
     :  '("ratio of rupture to s-wave velocity = ",f8.2)') 
     :     vrup_beta
      
      write(nu_Write_Params,'("FaultLat                  = ",f8.2)')
     :   FaultLat
      write(nu_Write_Params,'("FaultLon                  = ",f8.2)')
     :   FaultLon

      write(nu_Write_Params,'("No.of subs along strike   = ",i8)')nl
      write(nu_Write_Params,'("No.of subs along dip      = ",i8)')nw
      write(nu_Write_Params,'("subfault length           = ",f8.2)')dl
      write(nu_Write_Params,'("subfault width            = ",f8.2)')dw
      write(nu_Write_Params,'("i_rise_time (1=orig,2=1/f0) = ",i2)')
     :                      i_rise_time
      write(nu_Write_Params,'("iseed, nsims = ",1x,i11, 1x,i4)')
     :                      iseed, nsims
      write(nu_Write_Params,
     :   '("-----------------------------------------------")')
      write(nu_Write_Params,'("input hypocenter at position    = ",
     :                                             1p,2(1x,e10.3))')
     :                           hyp_loc_fl, hyp_loc_fw
      write(nu_Write_Params,'("input hypocenter at subfault    = ",
     :   2i4)')
     :                           i0_in,j0_in
      write(nu_Write_Params,'("n_hypocenters    = ",i4)') n_hypocenters
      write(nu_Write_Params,'("Mag.                      = ",f8.2)')
     :   amag
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      write(nu_Write_Params,'("dt (sec)                  = ",f8.4)')
     :   dt
      write(nu_Write_Params,'("beta (km/s)               = ",f8.2)')
     :   beta
      write(nu_Write_Params,'("density (rho), gr/cm3     = ",f8.2)')
     :   rho
      write(nu_Write_Params,'("prtitn                    = ",f8.2)')
     :   prtitn
      write(nu_Write_Params,'("rtp                       = ",f8.2)')
     :   rtp
      write(nu_Write_Params,'("fs                        = ",f8.2)')
     :   fs
      write(nu_Write_Params,'
     :   ("pulsing Percentage        = ",f8.2)')pulsingPercent
      write(nu_Write_Params,'(a,i1)') 
     :   'iflagscalefactor (1=vel^2; 2=acc^2;  '//
     :                     '3=asymptotic acc^2 (dmb)) = ',
     :               iflagscalefactor
      write(nu_Write_Params,'("flocut, nslope            = ",
     :   1x, f6.3, 1x, i2)') 
     :               flocut, nslope
      write(nu_Write_Params,'("iflagfas_avg              = ",i1)')
     :   iflagfas_avg
      write(nu_Write_Params,'("iflagpsa_avg_over_hypos   = ",i1)')
     :                                       iflagpsa_avg_over_hypos
      write(nu_Write_Params,'("iflagpsa_avg_over_sims   = ",i1)')
     :                                       iflagpsa_avg_over_sims
      write(nu_Write_Params,'("stress parameter (bars)   = ",f8.2)')
     :   stress
      write(nu_Write_Params,'("fmax                      = ",f8.2)')
     :   fmax
      write(nu_Write_Params,'("kappa                     = ",f8.4)')
     :   akappa
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      write(nu_Write_Params,'("         Corner Frequency    ")' )
      if (iDynamicFlag==1)then
          write(nu_Write_Params,'
     :   ("Dynamic Corner Frequency Flag is  ON = ",i4)')
     :                          iDynamicFlag
      else
          write(nu_Write_Params,'
     :   ("Dynamic Corner Frequency Flag is  OFF = ",i4)')
     :                          iDynamicFlag
          write(nu_Write_Params,'("Corner Frequency is static")')
      endif
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      write(nu_Write_Params,'(a)') 
     :     ' fr1, qr1, s1q, ft1, ft2, fr2, qr2, s2q, c_q, qt1, stq  = '
      write(nu_Write_Params,*) 
     :       fr1, qr1, s1q, ft1, ft2, fr2, qr2, s2q,c_q, qt1, stq 
      write(nu_Write_Params,'(a)') 
     :  ' Path duration: ndur_hinges, rpathdur, pathdur, durslope: '
      write(nu_Write_Params,'(1x,i2)') ndur_hinges
      do i = 1, ndur_hinges
        write(nu_Write_Params,'(1x,f5.1, 1x,f6.2)') rpathdur(i),
     :    pathdur(i)
      end do
      write(nu_Write_Params,'(1x,f6.3)') durslope
      
      write(nu_Write_Params,'(a)') 
     :    ' gspread: i, nsprd_segs, r_ref, rlow, a_s, b_s, m_s'
      do i = 1, nsprd_segs
        write(nu_Write_Params,*) 
     :    i, nsprd_segs, r_ref, rlow(i), a_s(i), b_s(i), m_s(i)
      end do

      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      if (iwind == 0) write (nu_Write_Params,120)
      if (iwind == 1) write (nu_Write_Params,130)
120   format ('window applied               = tapered boxcar')
130   format ('window applied               = Saragoni-Hart')
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      call trim_c(f_crustal_amps, nc_f_crustal_amps)
      write(nu_Write_Params,'(a)') ' Crustal amps from file '//
     :   f_crustal_amps(1: nc_f_crustal_amps)
      write(nu_Write_Params,'(4x,a, 5x,a)') 'freq', 'amp'
      do i = 1, n_crustal_amps
        write(nu_Write_Params, '(1x,f7.3, 1x,f7.3)') 
     :     freq_crustal_amps(i), amp_crustal_amps(i)
      end do
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      call trim_c(f_site_amps, nc_f_site_amps)
      write(nu_Write_Params,'(a)') ' site amps from file '//
     :   f_site_amps(1: nc_f_site_amps)
      write(nu_Write_Params,'(4x,a, 5x,a)') 'freq', 'amp'
      do i = 1, n_site_amps
        write(nu_Write_Params, '(1x,f7.3, 1x,f7.3)') 
     :     freq_site_amps(i), amp_site_amps(i)
      end do
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      
      
 
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      write(nu_Write_Params,'
     :   ("         Analytical Pulse parameters     ")' )
      write(nu_Write_Params,'
     :   ("-----------------------------------------------")')
      if (iPapaFlag==1)then
          write(nu_Write_Params,'("Analytical Flag is  ON   = ",i4)')
     :   iPapaFlag
          write(nu_Write_Params,'("Analytical Gama          = ",f8.3)')
     :   PapaGama
          write(nu_Write_Params,'("Analytical Nu            = ",f8.3)')
     :   PapaNu
          write(nu_Write_Params,'("Analytical T0            = ",f8.3)')
     :   PapaT0
          write(nu_Write_Params,'("Analytical Peak          = ",f8.3)')
     :   PapaPeak
      else
          write(nu_Write_Params,'("Analytical Flag is  OFF   =  ",i4)')
     :   iPapaFlag
      endif

      write(nu_Write_Params,*)
      write(nu_Write_Params,'(a,3(1x,f9.3))')
     :   'tpadl, tpadt, dt(sec)                     :',
     :    tpadl, tpadt, dt
     
      do i = 1, numberOfSites
        write(nu_Write_Params,'(a, 1x,i3, 1x,a, 1x,f8.2, 1x, f8.2)')
     :   'For site', i, 'siteLocation coordinates 1&2 = ', 
     :                siteLocation(i,1),  siteLocation(i,2)
      end do
      if (move_site) then
        write(nu_Write_Params,'(a)')
     : 'Site may have been moved to midpoint or end of '//
     : 'surface projection of upper edge of fault' 
      end if

      

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writePSA_FA(ioPSA, freq,averagePGA,averagePGV,
     :       minPGA, maxPGA, 
     :  averagePSA,nfreq,nsims,damp,amag,
     :  siteLat,siteLon,depth,r,azimuth,faultDistance,averageFA,isite,
     :  HypoDistance,fPar, nl, nw)
!     :  HypoDistance,fPar)
!       writes simulated average PSA spectrum at nfreq frequencies
!       into specified ascii file. First column is frequency in Hz,
!        second column is PSA value in cm/s**2

! Note: azimuth is in degrees

! Dates: 12/31/08 - Write sqrt(nl*nw)*psa, sqrt(nl*nw)*fas
!        01/06/09 - Write pgv, pga as first two values
!        01/13/09 - Write sqrt(nl*nw)*psa, sqrt(nl*nw)*fas and add SD, PSV.
!                 - Reverse order of FAS and PSA
!        11/21/11 - Write minPGA, maxPGA
!        04/12/13 - Change "R             =" to "distance from fault origin to site ="
!        05/15/15 - Delete fname, which is not used, and sqrt_n, which is set but not used.

      dimension freq(*),averagePSA(*),averageFA(*)
      character fPar*(*)
      real omega, omega_sq, minPGA, maxPGA
      
      pi=4.0*atan(1.0)
      twopi = 2.0 * pi
       
      write(ioPsa,"('**********************************************')")
      write(ioPsa,'(a, 1x,i5, a)') 
     :  '*******  Average PSA and Fourier spectrum over', nsims,
     :  ' simulations  *******'
      write(ioPsa,"('*******     Site #',i4)") isite
      write(ioPsa,50) int(damp)
      write(ioPsa,20) nfreq
      amag=anint(amag*10)/10.
      write(ioPsa,'(
     : "Mag.                               = ",f8.2)')amag
      write(ioPsa,'(
     : "siteLat                            = ",f8.2)')siteLat
      write(ioPsa,'(
     : "siteLon                            = ",f8.2)')siteLon
      write(ioPsa,'(
     : "depth                              = ",f8.2)')depth
      write(ioPsa,'(
     : "distance from fault origin to site = ",f8.2)')R
      write(ioPsa,'(
     : "azimuth from fault origin to site  = ",f8.2)')azimuth
      write(ioPsa,'(
     : "fault Dis.(km)                     = ",f8.2)')faultDistance
      write(ioPsa,'(
     : "Hypo  Dis.(km)                     = ",f8.2)')HypoDistance
      write(ioPsa,
     : "('Input Parameters file =  ',a)") fPar

      write (ioPsa,'(1x,a)') 'data format: '//
     :         '(1x,f10.3, 1p, 1x,e11.4, 5x,e11.4,'//
     :         ' 2x,e11.4, 1x,e11.4)'
      write(ioPSA, '(3x,a, 5x,a,
     :                     1x,a, 1x,a,
     :                     1x,a, 3x,a,
     :                     8x,a, 8x,a)') 
     :       'freq(Hz)', 'per(s)',
     :       'FAS(cm/sec)', 'avgPSA(cm/s**2)', 
     :                    'avgPSV(cm/s)', 'avgSD(cm)',
     :                    'minY', 'maxY'
 
      freq_dum = -1.0
      write (ioPsa,'(1x,f10.3, 8x,a3, 
     :               9x,a3, 1p, 13x,a3,  
     :               2x,e11.4, 9x,a3)') 
     :                   freq_dum,  'NaN', 
     :                   'NaN', 'NaN', 
     :                   averagePGV, 'NaN'
      
      freq_dum = 999.9
      per_dum = 0.0
      write (ioPsa,'(1x,f10.3, 1x,f10.3,
     :               9x,a3, 5x,es11.4, 
     :                            10x,a3, 9x,a3,
     :               1x,es11.4, 1x,es11.4)') 
     :                   freq_dum, per_dum,
     :                    'NaN', averagePGA,
     :                             'NaN','NaN',
     :                    minPGA, maxPGA
      
      do i=1, nfreq
      
        omega = twopi*freq(i)
        omega_sq = omega**2.0
        per_dum=1.0/freq(i)
        if(averageFA(i)>0) then      ! dmb            
          write (ioPsa,'(1x,f10.3, 1x,f10.3,
     :               1p, 1x,e11.4, 5x,e11.4, 
     :               2x,e11.4, 1x,e11.4)') 
     :              freq(i), per_dum,
     :                averageFA(i), averagePSA(i),
     :                averagePSA(i)/omega, averagePSA(i)/omega_sq
     
        else
        per_dum=1.0/freq(i)
          write (ioPsa,'(1x,f10.3,  1x,f10.3,
     :               9x,a3,    1p, 5x,e11.4, 
     :               1x,e11.4, 1x,e11.4)') 
     :              freq(i), per_dum,
     :              'NaN', averagePSA(i), 
     :              averagePSA(i)/omega, averagePSA(i)/omega_sq
        endif
      enddo

      write (ioPsa,*)
      write (ioPsa,*)
20    format (i4,' samples')
50    format ('damping: ',i3,'%')
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!----------------- BEGIN Q_func -----------------------------

! I added entry points so that the program could be called using a single 
! argument, without passing the other arguments through common.  Using
! a single argument is necessary when called by sme Numerical Recipes programs.

! Use:

! Call the setup entry point:
!      dummy = q_f_setup(fr1, qr1, s1q, ft1, ft2,
!     :              fr2, qr2, s2q, c_q)

! Call as a function:
!      q_fq = q_f(fq)



! Dates: 05/17/15 - Written by D.M. Boore, based on \smsim\q_f.for
 
      function q_func(f, qr1,fr1,s1q, qt1,ft1,stq, ft2, qr2,fr2,s2q ) 

      q_func = 9999.0
      if (f .eq. 0.0) return
        
      if ( f .le. ft1) then
        q_func = qr1*(f/fr1)**s1q
      else if ( f .ge. ft2) then
        q_func = qr2*(f/fr2)**s2q
      else
        q_func = qt1*(f/ft1)**stq
      end if

      return
  
      end
!----------------- END Q_func -----------------------------

 
 

      include 'exsim_dmb_included_subprograms.for'

!      include '\forprogs\accsqint.for' 
!      include '\forprogs\dcdt.for' 
!      include '\forprogs\locate.for' 
!      include '\forprogs\locate_d.for' 
!      include '\forprogs\zbrent.for'
!      include '\forprogs\dist_3df.for'
!      include '\forprogs\km2deg_f.for'
!      include '\forprogs\yintrf.for'

!      include '\smsim\gsprd_f.for'
!      include '\smsim\gsprd_q_f.for'
!      include '\smsim\gsprd_q_avg_f.for'
!      include '\smsim\q_f.for'




  
