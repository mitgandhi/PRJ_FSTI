#include "CGapLog.h"

CGapLog::CGapLog() 
{
	enable_cout = true;

	if(enable_cout)
	{
		cout << "Opening slipper_log.txt" << endl;
	}
	initialize("slipper_log.txt");
}

void CGapLog::initialize(const string LogFileName) {
	LogFile = LogFileName;
	GapLog.open(("./" + LogFileName).c_str(),ios::app);
}
CGapLog::~CGapLog() {
	GapLog << "\n\n*-------------------------------------------------------------------------------*" << std::endl;
	GapLog.close();
}
void CGapLog::message(const std::string& msg) 
{
	GapLog << msg << endl;
	
	if(enable_cout)
	{
		cout << msg << endl;
	}
}

//this handles endl, scientific, setprecision, etc..
CGapLog& CGapLog::operator<<( ostream&(*f) (ostream&))
{
	//should we send the output to the console?
	if(enable_cout)
	{
		//send the output to console
		cout << f;
	}

	//if the log file is open, also to the log file
	if(GapLog.is_open())
	{
		GapLog << f;
	}

	return *this;
}

