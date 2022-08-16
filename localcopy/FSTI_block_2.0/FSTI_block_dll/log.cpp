# include "./log.h"

gaplog::gaplog(const char* file) 
{ 
	logfile.open(file);

	logfile 

		<< "*-----------------------------------------------------------------------------*"
		<< "\n\n"
		<< "FSTI Gap Design\n"
		<<"Cylinder Block / Valve Plate Lubrication Module\n"
		<< "\nT U Dresden\n"
		<< "Maha Fluid Power Research Center\n"
		<< "https://engineering.purdue.edu/Maha"
		<< "\nPROPRIETARY AND CONFIDENTIAL.\n"
		<< "INFORMATION CONTAINED HEREIN SHALL NOT BE DECOMPILED, DISASSEMBLED,\n"
		<< "DUPLICATED, OR DISCLOSED IN WHOLE OR IN PART FOR ANY PURPOSE.\n"
		<< "UNAUTHORIZED USE IS STRICTLY PROHIBITED.\n\n"
		<< "*-----------------------------------------------------------------------------*"
		<< "\n\n\n";

}