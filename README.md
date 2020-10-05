# Adaptive Ion-Kinetic Electron-Fluid (AIKEF)

AIKEF is a plasma simulation code for space physics with two main fields of application:
* The interaction between a planetary obstacle and the impinging plasma flow. 
* Plasma turbulence.
	

## Required Software
* c++ compiler:
gcc (standard compiler for all linux distributions; easiest compiler to install, already installed on most systems)
* openmpi
* libraries (should already be included for 64-bit systems): libsiloh5, libhdf5, libsz, libz, libgsl
* Visualization tool: VisIt https://wci.llnl.gov/simulation/computer-codes/visit/


## Getting Started:
Before you start you should make sure you are using a 64-bit linux system (you can type "uname -m" into a terminal, the output should be \textit{x86\_64}).
* Makefile (to compile the code) already included in the root folder
* job file for cluster submission in the bin folder

 

## Some Basic Parameters:
There are two files which are modified by every user:
defines.h : general switches that decide if you want to include neutral species, ion production fields, dust, a dipole field ... and which terms of the field equation should be calculated
parameters.cpp : nearly all parameters of your simulation... Most important aspects are:


Set your normalization values and the properties of your upstream plasma (sections 0 and 4 of the parameter file)
Set box size and resolution:

		//! Number of Root Blocks in Szenario
		const INT32 RB_X = 4;
		const INT32 RB_Y = 4;
		const INT32 RB_Z = 4;
		
		//! Number of Root Blocks in Szenario
		const INT32 BlkNds_X = 10;
		const INT32 BlkNds_Y = 10;
		const INT32 BlkNds_Z = 10;
		
		//! Radius of Moon for Plotting etc
		const D_REAL R_Moon = 252.e+3/SI_x0;

		//! Size of simulation box
		const D_REAL LX = 200.;
		const D_REAL LY = 200.;
		const D_REAL LZ = 200.;

		//! Origin of simulationbox (= Position of Obstacle)
		const D_REAL Box_Origin[3] = {0.4001 * LX, 0.5001 * LY, 0.5001 * LZ};

		//! Radius of Obstacle
		const D_REAL R_Obstacle = 0;
		
		 //! End Run when TL_MAX is reached
		const INT32 TL_MAX = 12345;
		
		//! Numerical time step
		//! NOTE: dt should be at least 5 times smaller
		//!  than Courant Criteria suggest. 
		//! TODO MAX_LEVEL is set BELOW this parameter!
		const D_REAL  dt =  0.2* LX/(MAX_LEVEL*RB_X*(BlkNds_X-2)) / (SW_v/SI_v0);

Usually, simulations become stationary after 1 - 3 box crossings of the upstream plasma...

		const INT32  TL_OUTPUT_2D_SILO    =    400;
		const INT32  TL_OUTPUT_3D_SILO	 =     0;
		const INT32  TL_OUTPUT_3D_uniform_grid =    0;
		const INT32  TL_OUTPUT_LINEOUT    =   1;
		const INT32  TL_OUTPUT_PARTICLE_DETECTOR =  0;
		const INT32  TL_OUTPUT_TRAJECTORY = 0;
		const INT32  TL_OUTPUT_PARTICLE_TRACKS = 0;
		const INT32  TL_OUTPUT_GETMESH = 0;
		
 It may happen that your simulations stops before TL\_MAX is reached (crash or wall time on cluster). Thus, it is possible to save the simulation and restart it from that point
		
		 //! save state every xx time level (set to 0 to switch off)
		const INT32 TL_SAVE_STATE  = 1000;
		
The state is saved in the folder **State** in **bin**. Every time your start a simulation, it checks if a state file exists and if so, restores that state and continues. Note that the state file contains all particles and may thus be very large (several GB). If you run the simulation on a cluster that restricts your wall time, you may also save the state after a defined real computing time instead of just every TL\_SAVE\_STATE by using 
		
		 bool resubmit_enable = true;
		const INT32 resubmit_walltime =12*3600;

Unfortunately, an automatic resubmit is not possible on most clusters, but the state is still saved. Set further properties of your plasma (section 4 of the parameter file). Amount of particles in the simulation (apart from the number of nodes in the box, this is the most critical parameter for the simulation run time)
		
		 const INT32 optimal_MPiC[] = {400, 0, 0, 0,  0, 0};
		
The number of macroparticles in each cell is adjusted to this value. For quick and dirty tests on your desktop PC, 10 should be ok, for nice simulations with planetary bodies, it should be rather 100 and turbulence simulations may require several hundreds...

Just an example: assuming a simulation with 100^3 cells, with 100 Particles in each cell. A particle has a position, a velocity and a weight as double, thus 56 Byte. This simulation would need 100^3 * 100 * 56 = 5.6 GB	of memory just for the particles.... For simulations with planetary bodies one needs to set the properties of the body. These properties are a possible conductivity,
		
		 //! decide if use resistivity inside obstacle
		//! details are specified in CBlk_EtaProfiles.cpp
		const bool use_resistive_obstacle = false;
		
and a possible core,
		
		 //! within the core no magnetic field is advanced nor smoothed
		//! (-> it will remain the initial field forever)
		//! in % of obstacle radius
		const D_REAL obstacle_core_fraction = 0.0;
		

## Output and Visualization:
The logfile tells you that the simulation runs and has generated the first output (or it is crashed and you need the output to understand why). The usual output is writing the fields and plasma moments in the silo file format which can be opened by VisIt. This could be done either in 2D (three predefined cross-sections, SILO\_2D) or in 3D (whole simulation box, SILO\_3D).

Now you want to look at this output.
Start VisIt
Load the 2d or 3d output of your simulation.

Documentation mostly adapted from Hendrik Kriegel
