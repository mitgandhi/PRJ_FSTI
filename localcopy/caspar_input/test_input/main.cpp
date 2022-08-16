# include "../caspar_input/input.h"
# include <iostream>

using namespace std;

int main()
{


	// pushed again from marco

	input my_input;	

	cout << my_input.data.thermal.valveplate.calc_deflect_on.size() << endl;
	for(unsigned int i=0; i<my_input.data.thermal.valveplate.calc_deflect_on.size(); i++)
		cout << my_input.data.thermal.valveplate.calc_deflect_on[i] << endl;
	
	return 0;
}