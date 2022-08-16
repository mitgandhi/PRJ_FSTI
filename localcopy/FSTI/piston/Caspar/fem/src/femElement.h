class femElement {
	
protected:

	double E;							// Young module
	double nu;							// Poisson ratio

	// nodes coordinate
	double* x;
	double* y;
	double* z;
	// nodes index
	int* conn;

public:

	int* globalNumbering;				// global numbering for the element
	double** stiffMat;					// element stiffness matrix

	femElement();						// constructr
	
	// the destructor MUST BE VIRTUAL!!!!!!!
	// in this way, with a pointer to the base femElement object, the proper destructor
	// for each derived type will be called
	virtual ~femElement();				
	
	virtual void calcStiff();				// calculate the element stiffness matrix
	
};

class brick : public femElement
{
	int		intOp;					// integralOperator
	double	D[6][6];				// constitution matrix
	double	B[6][24];				// strain-displacement matrix
	double	dndr[8][3];				// derivative of shape functions respect natural reference system
	double	dndx[8][3];				// derivative of shape functions respect global reference system
	double	jacob[3][3];			// Jacobian matrix 
	double	jacobInv[3][3];			// Jacobian inverse matrix
	double	jacobDet;				// determinant of the Jacobian matrix
	double	jacobDetW;				// [J]w1w2w3 term
	double	v;						// element volume
	double	r1, r2, r3;				// integration points	
	
	// member functions
	void calcJacob();				// calculate the element Jacobian matrix
	void calcB();					// calculate the element strain-displacement matrix
	void calcD();					// calculate the element constitution matrix
	void calcBDB();					// calculate the [B]'[D][B] product

public:

	// functions
	brick(int intOp, double E, double nu, double* x, double* y, double*z, int* conn);	// constructor
	~brick();
	// calculate the element stiffness matrix
	void calcStiff();				
};

class tetra : public femElement
{

	double D[6][6];				// constitution matrix
	double B[6][12];			// strain-displacement matrix
	double jacobDet;			// determinant of the Jacobian matrix
	double v;					// element volume

	void calcD();				// calculate the costitutive matrix
	void calcB();				// calculate the element strain-displacement matrix

public:
	
	// constructor
	tetra(double E, double nu, double* x, double* y, double*z, int* conn);
	~tetra();
	// calculate the element stiffness matrix
	void calcStiff();		

};
