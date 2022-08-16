# ifndef __gaplog__
# define __gaplog__

# include <iostream>
# include <fstream>

class gaplog
{

	std::ofstream logfile;
   
public:

	gaplog(const char* file);
	
	~gaplog() { logfile.close(); }
	
	// function that takes a custom stream, and returns it
  typedef gaplog& (*log_manipulator)(gaplog&);

  // take in a function with the custom signature
  gaplog& operator<<(log_manipulator manip)
  {
    // call the function, and return it's value
    return manip(*this);
  }

	// define the custom endl for log
  static gaplog& endl(gaplog& stream)
  {
    std::cout << std::endl;
    stream.logfile << std::endl;
    return stream;
  }

	// operator << overload
  template <typename T>
  gaplog& operator<<(const T& x)
  {
		logfile.copyfmt(std::cout);
		logfile << x;
    std::cout << x;
		return *this;
  }
	
};


# endif