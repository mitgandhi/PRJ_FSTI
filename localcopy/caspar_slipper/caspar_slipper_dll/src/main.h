#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#pragma once

#define WINDOWS 1
#define LINUX 2

//What operating system will this be compiled for
#define OS WINDOWS

using namespace std;

const double PI = acos(-1.0);

template <class T> string n2s(const T number)
{
	ostringstream oss (ostringstream::out);
	oss.str("");
	oss << number;
	return oss.str();
};
