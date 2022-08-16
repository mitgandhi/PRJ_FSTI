#	ifndef __gapvectors__
#	define __gapvectors__

#	include <cmath>
#	include <iostream>

// ---------------- block gap vector, cartesian components ----------------- //
class ca_vector
{

	double coords[3];
	double pi() const { return 3.14159265358979323846; }

public:

	// default constructor
	ca_vector();
	// construct from components
	ca_vector(const double x, const double y, const double z = 0.0);
	ca_vector(const ca_vector& v);
	// return x component
	const double x() const;
	// return y component
	const double y() const;
	// return z component
	const double z() const;
	// magnitude
	const double r() const;
	// angle
	const double theta() const;
	
	// ------------------------- operator overloading ------------------------ //

	// = operator
	ca_vector& operator=(const double val);
	ca_vector& operator=(const ca_vector& v);
	// assign component
	double& operator[](int idx);
	// change sign
	friend const ca_vector operator-(const ca_vector& v);
	// multiply with scalar
	friend const ca_vector operator*(const ca_vector& v, const double s);
	friend const ca_vector operator*(const double s, const ca_vector& v);
	// divide by a scalar
	friend const ca_vector operator/(const ca_vector& v, const double s);
	// sum
	friend const ca_vector operator+(const ca_vector& left, const ca_vector& right);
	// subctration
	friend const ca_vector operator-(const ca_vector& left, const ca_vector& right);
	// dot product
	friend double operator*(const ca_vector& left, const ca_vector& right);
	// cross product
	friend const ca_vector cross(const ca_vector& left, const ca_vector& right);
	// comparison
	friend bool operator==(const ca_vector& left, const ca_vector& right);
	friend bool operator!=(const ca_vector& left, const ca_vector& right);
	// print	
	friend std::ostream& operator<<(std::ostream& os, const ca_vector& v);

};

inline ca_vector::ca_vector() 
{
	coords[0] = coords[1] = coords[2] = 0.0;
}

inline ca_vector::ca_vector(const double x, const double y, const double z) 
{
	coords[0] = x, coords[1] = y, coords[2] = z;
}

inline ca_vector::ca_vector(const ca_vector& v) 
{
	coords[0] = v.x(), coords[1] = v.y(), coords[2] = v.z();
}

inline const double ca_vector::x() const 
{ 
	return coords[0]; 
}

inline const double ca_vector::y() const 
{ 
	return coords[1]; 
}

inline const double ca_vector::z() const 
{ 
	return coords[2]; 
}

inline const double ca_vector::r() const 
{ 
	return sqrt( pow(coords[0],2.0) + pow(coords[1],2.0) + pow(coords[2],2.0) ); 
}

inline const double ca_vector::theta() const 
{ 
	return (atan2(coords[0], coords[1]) < 0) ? 2.0*pi() + atan2(coords[0], coords[1]) : atan2(coords[0], coords[1]);
}

inline ca_vector& ca_vector::operator=(const double val) 
{
	coords[0] = coords[1] = coords[2] = val;
	return *this;
}

inline ca_vector& ca_vector::operator=(const ca_vector& v) 
{
	coords[0] = v.x(), coords[1] = v.y(), coords[2] = v.z();
	return *this;
}

inline double& ca_vector::operator[](int idx) 
{
	# ifdef _DEBUG
	if(idx < 0 || idx > 2)
		std::cout << "ca_vector::operator[] out of range!" << std::endl;
	# endif
	
	return coords[idx];	
}

inline const ca_vector operator-(const ca_vector& v) 
{
	return ca_vector(-v.x(), -v.y(), -v.z());
}

inline const ca_vector operator*(const ca_vector& v, const double s) 
{
	return ca_vector(s*v.x(), s*v.y(), s*v.z());
}

inline const ca_vector operator*(const double s, const ca_vector& v) 
{
	return ca_vector(s*v.x(), s*v.y(), s*v.z());
}

inline double operator*(const ca_vector& left, const ca_vector& right) 
{
	return (left.x()*right.x() + left.y()*right.y() + left.z()*right.z());
}

inline const ca_vector operator/(const ca_vector& v, const double s) 
{
	return ca_vector(v.x()/s, v.y()/s, v.z()/s);
}

inline const ca_vector cross(const ca_vector& left, const ca_vector& right) 
{
	return ca_vector 
	(
		left.y()*right.z() - left.z()*right.y(),
		left.z()*right.x() - left.x()*right.z(),
		left.x()*right.y() - left.y()*right.x()
	);
}

inline const ca_vector operator+(const ca_vector& left, const ca_vector& right) 
{
	return ca_vector(left.x() + right.x(), left.y() + right.y(), left.z() + right.z());
}

inline const ca_vector operator-(const ca_vector& left, const ca_vector& right) 
{
	return ca_vector(left.x() - right.x(), left.y() - right.y(), left.z() - right.z());
}

inline  bool operator==(const ca_vector& left, const ca_vector& right) 
{
	return (left.x() == right.x() && left.y() == right.y() && left.z() == right.z());
}

inline bool operator!=(const ca_vector& left, const ca_vector& right) 
{
	return (left.x() != right.x() || left.y() != right.y() || left.z() != right.z());
}

inline std::ostream& operator<<(std::ostream& os, const ca_vector& v) 
{
	os << "(" << v.x() << " , " << v.y() << " , " << v.z() << ")";
	return os;
}

// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //
// ------------------------------------------------------------------------- //

// -------------- block gap vector, cylinderical components ---------------- //

class cy_vector
{

	double coords[3];
	double pi() const { return 3.14159265358979323846; }

public:

	// default constructor
	cy_vector();
	// construct from components
	cy_vector(const double x, const double y, const double z = 0.0);
	cy_vector(const cy_vector& v);
	// return x component
	const double x() const;
	// return y component
	const double y() const;
	// return z component
	const double z() const;
	// radius
	const double r() const;
	// angle
	const double theta() const;
	
	// ------------------------- operator overloading ------------------------ //

	// = operator
	cy_vector& operator=(const double val);
	cy_vector& operator=(const cy_vector& v);
	// assign component
	double& operator[](int idx);
	// change sign
	friend const cy_vector operator-(const cy_vector& v);
	// multiply with scalar
	friend const cy_vector operator*(const cy_vector& v, const double s);
	friend const cy_vector operator*(const double s, const cy_vector& v);
	// divide by a scalar
	friend const cy_vector operator/(const cy_vector& v, const double s);
	// sum
	friend const cy_vector operator+(const cy_vector& left, const cy_vector& right);
	// subctration
	friend const cy_vector operator-(const cy_vector& left, const cy_vector& right);
	// dot product
	friend double operator*(const cy_vector& left, const cy_vector& right);
	// comparison
	friend bool operator==(const cy_vector& left, const cy_vector& right);
	friend bool operator!=(const cy_vector& left, const cy_vector& right);
	// print	
	friend std::ostream& operator<<(std::ostream& os, const cy_vector& v);

};

inline cy_vector::cy_vector() 
{
	coords[0] = coords[1] = coords[2] = 0.0;
}

inline cy_vector::cy_vector(const double r, const double theta, const double z) 
{
	coords[0] = r, coords[1] = theta, coords[2] = z;
}

inline cy_vector::cy_vector(const cy_vector& v) 
{
	coords[0] = v.r(), coords[1] = v.theta(), coords[2] = v.z();
}

inline const double cy_vector::r() const 
{ 
	return coords[0];
}

inline const double cy_vector::theta() const 
{ 
	return coords[1];
}

inline const double cy_vector::z() const 
{ 
	return coords[2]; 
}

inline const double cy_vector::x() const 
{ 
	return coords[0]*sin(coords[1]); 
}

inline const double cy_vector::y() const 
{ 
	return coords[0]*cos(coords[1]); 
}

inline cy_vector& cy_vector::operator=(const double val) 
{
	coords[0] = coords[1] = coords[2] = val;
	return *this;
}

inline cy_vector& cy_vector::operator=(const cy_vector& v) 
{
	coords[0] = v.r(), coords[1] = v.theta(), coords[2] = v.z();
	return *this;
}

inline double& cy_vector::operator[](int idx) 
{
	# ifdef _DEBUG
	if(idx < 0 || idx > 2)
		std::cout << "cy_vector::operator[] out of range!" << std::endl;
	# endif
	
	return coords[idx];	
}

inline const cy_vector operator-(const cy_vector& v) 
{
	return cy_vector(-v.r(), -v.theta(), -v.z());
}

inline const cy_vector operator*(const cy_vector& v, const double s) 
{
	return cy_vector(s*v.r(), s*v.theta(), s*v.z());
}

inline const cy_vector operator*(const double s, const cy_vector& v) 
{
	return cy_vector(s*v.r(), s*v.theta(), s*v.z());
}

inline double operator*(const cy_vector& left, const cy_vector& right) 
{
	return (left.r()*right.r() + left.theta()*right.theta() + left.z()*right.z());
}

inline const cy_vector operator/(const cy_vector& v, const double s) 
{
	return cy_vector(v.r()/s, v.theta()/s, v.z()/s);
}

inline const cy_vector operator+(const cy_vector& left, const cy_vector& right) 
{
	return cy_vector(left.r() + right.r(), left.theta() + right.theta(), left.z() + right.z());
}

inline const cy_vector operator-(const cy_vector& left, const cy_vector& right) 
{
	return cy_vector(left.r() - right.r(), left.theta() - right.theta(), left.z() - right.z());
}

inline bool operator==(const cy_vector& left, const cy_vector& right) 
{
	return (left.r() == right.r() && left.theta() == right.theta() && left.z() == right.z());
}

inline bool operator!=(const cy_vector& left, const cy_vector& right) 
{
	return (left.r() != right.r() || left.theta() != right.theta() || left.z() != right.z());
}

inline std::ostream& operator<<(std::ostream& os, const cy_vector& v) 
{
	os << "(" << v.r() << " , " << v.theta() << " , " << v.z() << ")";
	return os;
}

#	endif