#	ifndef __point__
#	define __point__

#	include <cmath>
#	include <iostream>
 
class point 
{

	double p[3];
	
	double pi() const { return 3.14159265358979323846; }
	
public:

	// default constructor
	point();
	// construct from components
	point(const double x, const double y, const double z = 0.0);
	point(const point& P);
	// return x component
	const double x() const;
	// return y component
	const double y() const;
	// return z component
	const double z() const;
	// magnitude
	const double mag() const;
	// angle
	const double ang() const;
	// distance respect another point
	const double distance (const point& pt) const;
	// squared distance respect another point
	const double distance2 (const point& pt) const;

	// --------------------------------- operator overloading -------------------------------- //

	point& operator=(const double val);
	point& operator=(const point& pr);
	// assign component
	double& operator[](int idx);
	// change sign
	friend const point operator-(const point& p);
	// multiply with scalar
	friend const point operator*(const point& p, const double s);
	friend const point operator*(const double s, const point& p);
	// divide by a scalar
	friend const point operator/(const point& p, const double s);
	// sum
	friend const point operator+(const point& left, const point& right);
	// subctration
	friend const point operator-(const point& left, const point& right);
	// dot product
	friend double operator*(const point& left, const point& right);
	// cross product
	friend const point operator^(const point& left, const point& right); 
	// comparison
	friend bool operator==(const point& left, const point& right);
	friend bool operator!=(const point& left, const point& right);
	// print	
	friend std::ostream& operator<<(std::ostream& os, const point& p);

};

inline point::point() {
	p[0] = p[1] = p[2] = 0.0;
}

inline point::point(const double x, const double y, const double z) {
	p[0] = x, p[1] = y, p[2] = z;
}

inline point::point(const point& P) {
	p[0] = P.x(), p[1] = P.y(), p[2] = P.z();
}

inline const double point::x() const { 
	return p[0]; 
}

inline const double point::y() const { 
	return p[1]; 
}

inline const double point::z() const { 
	return p[2]; 
}

inline const double point::mag() const { 
	return sqrt( pow(p[0],2.0) + pow(p[1],2.0) + pow(p[2],2.0) ); 
}

inline const double point::ang() const { 
	return (atan2(p[1],p[0]) < 0) ? 2.0*pi() + atan2(p[1],p[0]) : atan2(p[1],p[0]);
}

inline const double point::distance (const point& pt) const {
	return sqrt( pow(p[0] - pt.x(),2.0) +  pow(p[1] - pt.y(),2.0) +  pow(p[2] - pt.z(),2.0) );
}

inline const double point::distance2 (const point& pt) const {
	return ( pow(p[0] - pt.x(),2.0) +  pow(p[1] - pt.y(),2.0) +  pow(p[2] - pt.z(),2.0) );
}

inline point& point::operator=(const double val) {
	p[0] = p[1] = p[2] = val;
	return *this;
}
inline point& point::operator=(const point& pr) {
	p[0] = pr.x(), p[1] = pr.y(), p[2] = pr.z();
	return *this;
}

inline double& point::operator[](int idx) 
{
	# ifdef _DEBUG
	if(idx < 0 || idx > 2)
		std::cout << "point::operator[] out of range!" << std::endl;
	# endif
	
	return p[idx];	
}

inline const point operator-(const point& p) {
	return point(-p.x(), -p.y(), -p.z());
}

inline const point operator*(const point& p, const double s) {
	return point(s*p.x(), s*p.y(), s*p.z());
}

inline const point operator*(const double s, const point& p) {
	return point(s*p.x(), s*p.y(), s*p.z());
}

inline double operator*(const point& left, const point& right) {
	return (left.x()*right.x() + left.y()*right.y() + left.z()*right.z());
}

inline const point operator/(const point& p, const double s) {
	return point(p.x()/s, p.y()/s, p.z()/s);
}

inline const point operator+(const point& left, const point& right) {
	return point(left.x() + right.x(), left.y() + right.y(), left.z() + right.z());
}

inline const point operator-(const point& left, const point& right) {
	return point(left.x() - right.x(), left.y() - right.y(), left.z() - right.z());
}

inline const point operator^(const point& left, const point& right) {
	return point (
				left.y()*right.z() - left.z()*right.y(),
				left.z()*right.x() - left.x()*right.z(),
				left.x()*right.y() - left.y()*right.x()
			);
}

inline  bool operator==(const point& left, const point& right) {
	return (left.x() == right.x() && left.y() == right.y() && left.z() == right.z());
}

inline bool operator!=(const point& left, const point& right) {
	return (left.x() != right.x() || left.y() != right.y() || left.z() != right.z());
}

inline std::ostream& operator<<(std::ostream& os, const point& p) {
	os << "(" << p.x() << " , " << p.y() << " , " << p.z() << ")";
	return os;
}

#	endif
