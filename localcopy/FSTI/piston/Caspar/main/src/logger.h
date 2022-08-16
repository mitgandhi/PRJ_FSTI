#include <iostream>
#include <fstream>

using namespace std;

extern bool coupled_caspar_simulation;
extern bool output_piston_window;

class Logger
{
public:
	ofstream logfile;
	
	//this handles "\n", scientific, setprecision, etc..
	Logger& operator<<( ostream&(*f) (ostream&))
	{
		//should we send the output to the console?
		if(!coupled_caspar_simulation && output_piston_window)
		{
			//send the output to console
			cout << f;
		}

		//if the log file is open, also to the log file
		if(logfile.is_open())
		{
			logfile << f;
		}

		return *this;
	}

	//this handles strings, doubles, etc
	template<typename T> Logger& operator<<(T& str)
	{
		//should we send the output to the console?
		if(!coupled_caspar_simulation && output_piston_window)
		{ 
			//send the output to console
			cout << str;
		}

		//if the log file is open, also to the log file
		if(logfile.is_open())
		{
			logfile << str;
		}

		return *this;
	}
	
	void open(const string filename)
	{
		//append to the file if it exists
		logfile.open(filename.c_str(), ios_base::app);

		if(logfile.is_open())
		{
			//maybe write a header here

				/*logfile << "***************************************************************" << "\n";
				logfile << "*                                                             *" << "\n";
				logfile << "*                            /\\_/\\                            *" << "\n";										// I think he looks better with horns -- dan  cout << "*                           /''~``;                           *" << "\n";
				logfile << "*                           ( o o )                           *" << "\n";
				logfile << "*   +------------------.oooO--(_)--Oooo.------------------+   *" << "\n";
				logfile << "*   |                                                     |   *" << "\n";
				logfile << "*   |              PISTON/CYLINDER INTERFACE              |   *" << "\n";
				logfile << "*   |                  FSTI Gap Design                    |   *" << "\n";
				logfile << "*   |                Maha Internal Version                |   *" << "\n";
				logfile << "*   |                                                     |   *" << "\n";
				logfile << "*   |                    .oooO                            |   *" << "\n";
				logfile << "*   |                    (   )   Oooo.                    |   *" << "\n";
				logfile << "*   +---------------------( (----(   )--------------------+   *" << "\n";
				logfile << "*                          (_)    ) )                         *" << "\n";
				logfile << "*                                (_)                          *" << "\n";
				logfile << "*                                                             *" << "\n";
				logfile << "*                 Matteo Pelosi - Dan Mizell                  *" << "\n";
				logfile << "*              Maha Fluid Power Research Center               *" << "\n";
				logfile << "*                     Purdue University                       *" << "\n";
				logfile << "*             https://engineering.purdue.edu/Maha             *" << "\n";
				logfile << "*                                                             *" << "\n";
				logfile << "***************************************************************" << "\n";
				logfile << "\n";*/
		}
	}

	~Logger()
	{
		if(logfile.is_open())
		{
			logfile.close();
		}
	}

};

extern Logger Log;