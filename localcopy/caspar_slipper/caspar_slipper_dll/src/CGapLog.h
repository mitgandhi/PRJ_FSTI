#include "main.h"
#pragma once

class CGapLog {
	std::ofstream GapLog;
	string LogFile;

	void initialize(const string LogFileName);

public:

	//should the output also be piped to cout?
	bool enable_cout;

	// ofstrem for log text file
	CGapLog ();
	~CGapLog ();

	//used in legacy code
	void message(const std::string& msg);

	//A logger stream 'manipulator' that indicates an error
	class error
	{
		const int err_code;
	public:
		error(const int error_code = 1) : err_code(error_code)
		{
		};

		//This will cause the process to exit when the object is distructed
		~error()
		{
			exit(err_code);
		}
	};

	//new operators
	CGapLog& operator<<( ostream&(*f) (ostream&));

	//handles the custom error manipulator
	CGapLog& operator<< (const error& e)
	{
		*this << endl << endl;
		*this << "ERROR: ";
		return *this;
	};

	template<typename T> CGapLog& operator<<(const T& str)
	{
		//should we send the output to the console?
		if(enable_cout)
		{ 
			//send the output to console
			cout << str;
		}

		//if the log file is open, also to the log file
		if(GapLog.is_open())
		{
			GapLog << str;
		}

		return *this;
	};


};
