/*****************************************************************************/
//CASPAR SLIPPER DLL HEADER FILE

//Purdue University
//Maha Fluid Power Research Center
//https://engineering.purdue.edu/Maha

//PROPRIETARY AND CONFIDENTIAL.
//CONTAINS TRADE SECRET INFORMATION.
//INFORMATION CONTAINED HEREIN SHALL NOT BE DECOMPILED, DISASSEMBLED, 
//DUPLICATED, OR DISCLOSED IN WHOLE OR IN PART FOR ANY PURPOSE.
//UNAUTHORIZED USE IS STRICTLY PROHIBITED.
/*****************************************************************************/

//Note:
//caspar_input.h should be included before including this file

#ifdef _WIN32
#include <ppl.h>

//This is the primary hook function used by coupled caspar
__declspec(dllexport) void run_slipper_gap(Concurrency::event & local_event, Concurrency::event & global_event, const input & I);

//This function will return the current version and date/time stamp of the caspar_slipper.dll
__declspec(dllexport) string slipper_version();

//This function can be called to run a standalone Caspar Slipper simulation executable
__declspec(dllexport) int slipper_standalone(int argc, char *argv[], const input & I);

#else

//This function will return the current version and date/time stamp of the caspar_slipper.dll
string slipper_version();

//This function can be called to run a standalone Caspar Slipper simulation executable
int slipper_standalone(int argc, char *argv[]);

#endif
