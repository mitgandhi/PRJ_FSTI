// *************************************************************************** //
// CASPAR BLOCK DLL HEADER FILE

// Purdue University
// Maha Fluid Power Research Center
// https://engineering.purdue.edu/Maha

// PROPRIETARY AND CONFIDENTIAL.
// CONTAINS TRADE SECRET INFORMATION.
// INFORMATION CONTAINED HEREIN SHALL NOT BE DECOMPILED, DISASSEMBLED, 
// DUPLICATED, OR DISCLOSED IN WHOLE OR IN PART FOR ANY PURPOSE.
// UNAUTHORIZED USE IS STRICTLY PROHIBITED.
// *************************************************************************** //

# include <ppl.h>
# include "../FSTI_input/input.h"


// This function can be called to run a standalone Block Slipper simulation executable
__declspec(dllexport) int run_block_standalone(int argc, char *argv[]);

// This is the primary hook function used by coupled caspar
__declspec(dllexport) void run_block_gap(Concurrency::event & local_event, Concurrency::event & global_event, const input& in);
