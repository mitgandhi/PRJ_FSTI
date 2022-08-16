# ifndef __th_load__
# define __th_load__

# include <vector>
# include <map>
# include "./htr_analysis.h"

// ------------------------------------------------------------------------- //
class th_load
{
	
	friend class htr_analysis;

protected:

	htr_analysis& fea;	

public:

	// <global dof index where the load is applied, value of the load>
	std::map<int, double> lgdof_val;

	th_load(htr_analysis& _fea) : fea(_fea) {}
	~th_load() {}

	virtual void apply() = 0;

};
// ------------------------------------------------------------------------- //
class heat_flux : public th_load	// heat flux load (neaumann th_boundary)
{
	 
	const char* set;
	std::vector<double> q;

public:

	heat_flux(htr_analysis& _fea, const char* _set, std::vector<double> _q) 
		: th_load(_fea), set(_set), q(_q) {}

	void apply();
	
};
// ------------------------------------------------------------------------- //
class convection : public th_load	// convection load (mixed th_boundary)
{
	const char* set;
	double h;
	double Tinf;

public:

	convection(htr_analysis& _fea, const char* _set, double _h, double _Tinf) 
		: th_load(_fea), set(_set), h(_h), Tinf(_Tinf) {}

	void apply();
};
// ------------------------------------------------------------------------- //

# endif