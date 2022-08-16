# ifndef __blockode__
# define __blockode__

# include "../FSTI_block_gap/block_gap.h"
# include "../FSTI_input/input.h"
# include <fstream>
# include <ctime>

class cblock_ode
{

	const input& in;				// reference to input structure
	cblock_gap block_gap;		// main gap structure

	double step_angle;				// simulation step angle
	double omega;						// pump speed
	double t;								// time
	double phi;							// shaft angle
	double hmin;						// minimum film thickness
	
	double dphi;						// step angle [rad]
	double dt;							// time step
	double T;								// time for one revolution
	double t_end;						// simulation end time
	int revcounter;					// revolution counter

	time_t overall_time;		// simulation time
	struct tm* timeinfo;
	time_t clock_start;			
	time_t clock_end;

	std::ofstream out;			// output file

	bool EHD;								// true if EHD_CB and/or EHD_VP are on
	bool EHD_CB;						// true if EHD_CB is on
	bool EHD_VP;						// true if EHD_VP is on
	bool Thermal;						// true if Thermal_CB and/or Thermal_VP are on
	bool Thermal_CB;				// true if Thermal_CB is on
	bool Thermal_VP;				// true if Thermal_CB is on
	bool dense_vtk;					// dense vtk output
	bool dense_light_vtk;		// dense ligh vtk output
	bool EHD_debug;					// write debug info on the EHD loop
	int contact_opt;				// contact option: 
													//   0 -> do nothing 
													//   1 -> saturate squeeze for points at min film
													//   2 -> use squeeze velocity correction
	void print_state();			// print the solver status
	void write_output();		// write txt output file
	void end_revolution();	// check when a revolution is completed

	// derived output
	struct
	{
		double leak;						// leakage
		double Mfr;							// friction moment
		double leak_avg;				// leakage for one revolution
		double Qin_leak;				// heat in the gap due to leakage
		double Qout_leak;				// heat out the gap due to leakage
		double phid;						// heat dissipated by viscous friction 
		double Qin_leak_avg;		// heat in the gap due to leakage (average over one rev)
		double Qout_leak_avg;		// heat out the gap due to leakage (average over one rev)
		double phid_avg;				// heat dissipated by viscous friction (average over one rev)
		double Pl;							// Total power loss
		double Plf;							// Loss due to friction
		double Pll;							// Loss due to leakage
		double FTK;							//friction in the piston cylinder interface
	} output;

public:

	cblock_ode(const input& in);		// constructor
	~cblock_ode();									// destructor
	void initialize();							// initialize the ode solver
	void run();											// run the ODE
	void testing();									// testing function

};

# endif