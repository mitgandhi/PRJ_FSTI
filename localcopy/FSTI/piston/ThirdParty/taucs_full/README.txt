*---------------------------------------------------------------------------------------------*
|			Set up visual studio 2005 with tucs library			      |
*---------------------------------------------------------------------------------------------*

1) unpack the tarball archive (http://www.tau.ac.il/~stoledo/taucs/2.2/taucs_full.tgze),
   let say in the DIR directory. 
   You will get DIR\taucs_full that contains

	src
	progs
	matlab
	external
	doc
	configurator
	config
	testscript.bat
	testscript
	configure.bat
	configure

2) run the configure.but script (just double click on it)

3) open cmd (start->run->cmd) and go to DIR\taucs_full and type nmake

4) at the end of the building process you should have in DIR\taucs_full

	bin
	build
	config
	configurator
	doc
	external
	lib
	matlab
	obj
	progs
	scr
	.lastconf
	configure.bat
	makefile
	testscript
	testscript.bat
	
5) The installation is complete. Now you have to set up Visual Studio (VS). 
   The VS setup is divided in two parts. The first concern General options of VS, and
   the second the properties of the project that will be used for create your program.

   5.a) VS general options

   Open VS, and go to Tools->Options->Projects and Solutions->VC++ Directories

   Under "Include files" add these paths

	DIR\taucs_full\build\win32
	DIR\taucs_full\src

   Under "Library files" add these path

	DIR\taucs_full\external\lib\win32
	DIR\taucs_full\lib\win32

   5.b) Project options

   After the creation of the project right click on the project name in the solution
   explorer bar. A menu will appear; click on "Properties". A windows will be displayed.
   Make sure to select "All configurations" under Configuration.
   Go to Linker -> Input. In the "Additional dependencies" field, add the following entries

	libtaucs.lib
	blas_win32.lib
	lapack_win32.lib
	libatlas.lib
	libcblas.lib
	libf77blas.lib
	liblapack.lib
	libmetis.lib
	vcf2c.lib

   NOTE: the first library is the Taucs library, VS finds it in DIR\taucs_full\lib\win32 
         (see 5.a). The other libraries are external libraries that Taucs uses, VS finds 
         them in DIR\taucs_full\external\lib\win32. This libraries are precompiled and
         are part of the Taucs full distribution. You cuold compile them by your own...

   Go to C\C++ -> Code Generation. (Note that if in your project is not included any 
   C++ file) yuo will not see this menu!). Because Taucs is compiled using the 
   Multi-threaded options, yuo have to say to use MD option.
   Select Configuration -> Debug and set "runtime library" to "Multi-threaded Debug (/MTd)"
   Select Configuration -> Release and set "runtime library" to "Multi-threaded (/MT)"

6) If you want to use a the C++ compiler, so work with cpp files, you have to include the
   taucs header with this syntax;

	extern "C" {
	#include <taucs.h>
	}


That's it!

Now you shold be able to use the library in yuor own programs with VS!
Thank's to Marco Zecchi for this instructions
         




