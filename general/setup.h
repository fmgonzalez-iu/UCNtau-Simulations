#pragma once // Make sure these aren't called multiple times

#define INITBLOCK 0 // RNG Seed

//======================================================================
// Constants
//----------------------------------------------------------------------
// General Physical Constants
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27
#define HBAR 1.054571800e-34
#define EPSILON 1.0e-9

//----------------------------------------------------------------------
// Neutron Physical Constants

#define MU_N -9.6623647e-27 // Sign from fields_nate.c
#define MASS_N 1.674927471e-27 // Mass of Neutron
#define MINU -2.390352484438862e-26 // Minimum energy potential in trap
#define TAU_N 877.7 // Lifetime of neutron (from last result)

//----------------------------------------------------------------------
// Magnetic field parameters

#define B_HOLD 0.005 //holding field strength
#define B_REM 1.35 //remnant magnet field strength
#define MAG_THICK 0.0254 //thickness of layer of PM array
#define MAG_SPACE 0.05114 //characteristic spacing of magnets
#define N_TERMS 3 //how far out to go in field ripple expansion
#define KAPPA 1000 // Overlap area between asymmetric properties

//----------------------------------------------------------------------
// Absorption parameters
#define NBORONB2O3 8.824e27
#define NBORON 1.37e29
#define ABORON -0.1e-15
#define SIGMABORON 2200*3.835e-25

#define NOXYGENB2O3 1.32e28
#define AOXYGEN 5.803e-15
#define SIGMAOXYGEN 4.232e-28

#define NCARBON 1.133e29
#define ACARBON 6.6460e-15
#define SIGMACARBON 5.551e-28

#define NZINC 2.527e28
#define AZINC 5.68e-15
#define SIGMAZINC 5.241e-28

#define NSULFUR 2.527e28
#define ASULFUR 2.847e-15
#define SIGMASULFUR 1.556e-28

//======================================================================
// RUNNING FUNCTIONS
//----------------------------------------------------------------------
// Generating 
// Only randomPointTrap should be needed
#define TRACKGENERATOR randomPointTrap 

// Tracking
// #define TRACKER daggerHitTimes
// #define TRACKER fixedEffDaggerHitTime
// #define TRACKER fixedEffDaggerHitTime_PSE
// fixedEffDaggerHitTime_reinsert is the best "Dagger Only" one
#define TRACKER fixedEffDaggerHitTime_reinsert 
// #define TRACKER cleanTime
// #define TRACKER calcLyap
// #define TRACKER noabsCleanTime
// And NoAbsCleanTrack_multidata is the best "2 detector" one.
// #define TRACKER noAbsCleanTrack_multidata
//#define TRACKER timeToBounce


// Output
#define NRECORDS 50
//#define MULTIDETSAVE false

//#define WRITER writeNoabsRes
//#define WRITER writeNoAbsCleanDagRes
//#define WRITER writeNoAbsCleanDagResZ
//#define WRITER writeBounceRes
 #define WRITER writeFixedRes
// #define WRITER writeCleanRes
// #define WRITER writeLyapRes

// Turn to 1 if we want to count steps (debug)
#define TRACKANDPRINT 1

//======================================================================
// GENERATING FUNCTIONS
//----------------------------------------------------------------------
// Energy spectrum (unoptimized)
//#define EMIN 0.5 // Minimum/maximum energies
//#define EMAX 50.0
//#define EPOW 1.0 // Energy 
//#define THETAPOW 0.0 // Forward-directedness
// #define ZETACUT 0.0
// #define BTHICK 20.0

// Energy spectrum (optimized 2016/2017)
#define EMIN 7.2092 // Minimum energy
#define EMAX 34.444 // Maximum energy
#define EPOW 1.17727
#define THETAPOW 0.275457
#define ZETACUT 0.0127203 // Damaged area of dagger
#define BTHICK 5.59909 // Boron thickness

//======================================================================
// GEOMETRY
//----------------------------------------------------------------------
// Small Cleaner size
#define SMALLCLEANERXMIN (0.115529/2 + 0.212841 - 0.6604)
#define SMALLCLEANERXMAX (0.115529/2 + 0.212841)
#define SMALLCLEANERYMIN -(0.218041/2 + 0.335121 + 0.3556)
#define SMALLCLEANERYMAX -(0.218041/2 + 0.335121)

//----------------------------------------------------------------------
// Cleaning parameters
#define CLEANINGHEIGHT 0.38 // Lowered Cleaning Height
#define RAISEDCLEANINGHEIGHT 0.43 // Raised Cleaning Height
#define ECLEAN 5.571749397933261e-27 // In J

//#define GC_EFF 1.0 // Giant cleaner efficiency (comment out if you want 100% efficiency)
//#define AC_EFF 0.5 // Small active cleaner efficiency
#define AC_REINSERT true // Set to TRUE if we want to lower AC during unload.

// Cleaner Motion (for simplicity assume both AC/GC are same)
#define CLEANERVEL 5.0

//----------------------------------------------------------------------
// Dagger Motion 2016-2017
#define DAGGERACC 6.0
#define DAGGERVEL 1.6
#define DAGGERSPR 51200/2.0

//----------------------------------------------------------------------
// Aluminum block
#define BLOCKSCALE 0.0

//======================================================================
// RUN PARAMETERS
//----------------------------------------------------------------------
// Timing parameters
#define FILLINGTIME 150 // Filling time (300 for RH)
#define FILLINGEXPTIME 70 // Exponent of filling time
#define FIRSTDIPTIME 1550 // Start of first dip
#define HOLDTIMESHORT 20
#define HOLDTIME 20

#define CLEANINGTIME 50

//----------------------------------------------------------------------
// Unload Dips (NDIPS - 1)

//holdT gets replaced by HOLDTIME inside the tracking function

// #define NDIPS 3
// #define HEIGHTS {0.49, 0.38, 0.010}
// #define ENDTIMES {holdT, holdT+300, holdT+500}

#define NDIPS 4
#define HEIGHTS {0.49, 0.380, 0.250, 0.010}
#define DIPTIMES {0.0, 40.0, 60.0, 160.0}
#define ENDTIMES {HOLDTIME, HOLDTIME+40.0, HOLDTIME+60.0,HOLDTIME+160.0} // Should be DIPTIMES + holdT

// #define NDIPS 10
// #define HEIGHTS {0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}
// #define DIPTIMES {0.0,  40.0,  80.0,  100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 310.0}
// #define ENDTIMES {HOLDTIME,  HOLDTIME+40.0,  HOLDTIME+80.0,  HOLDTIME+100.0, HOLDTIME+120.0, HOLDTIME+140.0, HOLDTIME+160.0, HOLDTIME+180.0, HOLDTIME+200.0, HOLDTIME+310.0}

// #define NDIPS 13
// #define HEIGHTS {0.49, 0.38, 0.250,   0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010}
// #define ENDTIMES {0.0, 50.0, 250.0,  250.0+holdT, 250.0+holdT+20.0, 250.0+holdT+40.0, 250.0+holdT+50.0, 250.0+holdT+60.0, 250.0+holdT+70.0, 250.0+holdT+80.0, 250.0+holdT+90.0, 250.0+holdT+100.0, 250.0+holdT+120.0}

// Deep Dagger Cleaning
// #define NDIPS 10
// #define HEIGHTS {0.49, 0.38, 0.25, 0.18, 0.14, 0.11, 0.08, 0.06, 0.04, 0.01}
// #define ENDTIMES {holdT, holdT+40.0, holdT+80, holdT+100, holdT+120, holdT+140, holdT+160, holdT+180, holdT+200.0, holdT+500.0}

//======================================================================
// SYSTEMATICS
//----------------------------------------------------------------------
// Heating
#define HEATMULT 1.0
#define SAMPDT 0.0004 // Sampling frequency of vibration field

// On Supercomputer:
#define XFNAME "/N/u/frangonz/BigRed3/UCNtau-Simulations/simulations/xvals.bin"
#define YFNAME "/N/u/frangonz/BigRed3/UCNtau-Simulations/simulations/yvals.bin"
#define ZFNAME "/N/u/frangonz/BigRed3/UCNtau-Simulations/simulations/zvals.bin"
// On Localcomputer:
//#define XFNAME "./xvals.bin"
//#define YFNAME "./yvals.bin"
//#define ZFNAME "./zvals.bin"

// Magnetic field oscillation
// #define FREQ 60
// #define FREQ 1000
// #define AMPLITUDE 0.000010
// #define AMPLITUDE 0.000005
// #define AMPLITUDE 0.000001

//----------------------------------------------------------------------
// Chaos

#define LYAPTIME 5
#define NUMSEP 100

// Lyuponev scaling factors
#define XSCALE 0.2
#define YSCALE 0.35
#define ZSCALE 1.34
#define PSCALE 1.25e-27


//======================================================================
// Out of context old values
//----------------------------------------------------------------------
// Old generating functions
// #define TRACKGENERATOR randomPointTrapEdE
// #define TRACKGENERATOR randomPointTrapOptimum
// #define TRACKGENERATOR randomPointTrapOptimumCleanable
// #define TRACKGENERATOR randomPointTrapEdECleanable
// #define TRACKGENERATOR randomPointTrapOptimumOnlyCleanable

// Old Tracking Functions


// #define B_REM 1.4 //remnant magnet field strength
// #define MAG_SPACE 0.02 //characteristic spacing of magnets
// #define MAG_SPACE 0.0508 //characteristic spacing of magnets
// #define mu -9.662364e-27/1.674927351e-27 //mu in units where m=1

// #define EMIN 7.6071
// #define EMIN 34.776124570119975
// #define EPOW 1.30448
// #define THETAPOW 0.264079 // Angle
// #define ZETACUT 0.0160763
// #define BTHICK 5.76556

// #define CLEANINGHEIGHT 0.35
// #define CLEANINGHEIGHT 0.43
// #define ECLEAN 5.077298660340679e-27

// #define FIRSTDIPTIME 1380
// #define HOLDTIME 1380
