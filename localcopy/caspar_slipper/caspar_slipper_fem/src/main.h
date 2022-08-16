#include <iostream>
#include <vector>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include "../../caspar_slipper_dll/src/CGapLog.h"

using namespace std;

#pragma once

//Logging object
extern CGapLog GapLog;

namespace CasparSlipperFEM
{
	template <class T> string n2s(const T number)
	{
		ostringstream oss (ostringstream::out);
		oss.str("");
		oss << number;
		return oss.str();
	};

	enum analysis_types { NONE, ELASTIC, THERMAL };
};
