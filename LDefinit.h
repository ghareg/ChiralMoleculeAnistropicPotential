#ifndef LDEFINIT_H_
#define LDEFINIT_H_
#include <Eigen/Dense>
#include <cmath>
#include "constants.h"
typedef std::complex<double> Complex;
const Complex I(0.0, 1.0);
const double pi14 = std::pow(Pi, 0.25);
using Eigen::Matrix2cd;

struct PauliMatrix
{
	Matrix2cd pal0;
	Matrix2cd pal1;
	Matrix2cd pal2;
	Matrix2cd pal3;
};

struct LValStruct
{
	Matrix2cd L1;
	Matrix2cd L2;
};

struct GValStruct
{
	Matrix2cd G1;
	Matrix2cd G2;
	Matrix2cd G3;
	Matrix2cd G4;
};

struct PGValStruct
{
	Matrix2cd PG1;
	Matrix2cd PG2;
	Matrix2cd PG3;
	Matrix2cd PG4;
	Matrix2cd PG5;
};

struct Param
{

	double k;
	double alpha;
	double beta;
	double theta;
	double tau;
	double kSinal;
	double kCosal;
	double kSinth;
	double kCosth;
	double sinBeta;
	double cosBeta;
	double sinTau;
	double cosTau;
	double alphax;
	double alphay;

	void updateValues(double kc, double alphac, double betac, double thetac, double tauc, double alphaxc, double alphayc);
};

void initPauliMatrix(PauliMatrix& pal);
void updateLVal(LValStruct& LVal, int n1, int n2, const Param& pm, const PauliMatrix& pal, const double* fact, const double* HermiteMatx, const double* HermiteMaty);
#endif
