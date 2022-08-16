#include "logger.h"
#include "time.h"
#include <Windows.h>
#include <stdio.h>
#include <direct.h>

using namespace std;

//the global input log
logger inputlog(false);

logger::logger(const bool Silent) : silent(Silent)
{
	logstream.open("./input_log.txt");
	if(!logstream.is_open())
	{
		//error("Unable to open ./input_log.txt");
	}

	logstream << "*************************************************************************" << endl;
	logstream << "*                                                                       *" << endl;
	logstream << "* FSTI Gap Design                                                       *" << endl;
	logstream << "* Version: Run \"About FSTI Gap Design\"                                *" << endl;
	logstream << "* Maha Fluid Power Research Center                                      *" << endl;
	logstream << "* Purdue University                                                     *" << endl;
	logstream << "* https://engineering.purdue.edu/Maha                                   *" << endl;
	logstream << "*                                                                       *" << endl;
	logstream << "*************************************************************************" << endl;
	logstream << endl;
	logstream << "Time: " << gettime() << endl;
	logstream << "Computer: " << getcomputername() << endl;
	logstream << "Working Directory: " << workingdir() << endl;
	logstream << endl;
	logstream << "Begin processing FSTI Inputs ..." << endl;
	logstream << endl;

}

logger::~logger(void) 
{
	if(logstream.is_open())
	{
		logstream << endl;
		logstream << "Closing FSTI Inputs ..." << endl;
		logstream.close();
	}
}

//this handles endl, scientific, setprecision, etc..
logger& logger::operator<<( ostream&(*f) (ostream&))
{
	//should we send the output to the console?
	if(!silent)
	{
		//send the output to console
		cout << f;
	}

	//if the log file is open, also to the log file
	if(logstream.is_open())
	{
		logstream << f;
	}

	return *this;
}
//get the current date / time
string logger::gettime(void)
{
	time_t rawtime;
	struct tm timeinfo;
	char buffer [80];
	time ( &rawtime );
	localtime_s (&timeinfo, &rawtime);
	strftime (buffer,80,"%x %X",&timeinfo);
	return buffer;
};
//get the working program dir
string logger::workingdir(void)
{
	char cdir[FILENAME_MAX] = {0};
	_getcwd(cdir, sizeof(cdir));

	return cdir;
};
//get the windows system computer name
string logger::getcomputername(void)
{
	unsigned long computerNameLength = 256;
	char name[256];
	GetComputerNameA(name, &computerNameLength);
	return string(name);
};
