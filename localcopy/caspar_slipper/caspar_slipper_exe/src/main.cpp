#include "../../caspar_input/caspar_input.h"
#include "../../caspar_slipper_dll/src/caspar_slipper_dll.h"

int main(int argc, char *argv[])
{
	//create the caspar_input object
	input gapinput;

	return slipper_standalone(argc, argv, gapinput);
}
