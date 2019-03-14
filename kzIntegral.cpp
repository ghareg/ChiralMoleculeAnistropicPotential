#include "kzIntegral.h"
#include "rIntegral.h"

Complex Iz1(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		return I * (std::exp(I * (k1 - k2) * a) - std::exp(I * (k1 - k2) * b)) / (k1 - k2);
	}
	else {
		return (b - a);
	}
}

Complex Iz2(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		return ((1.0 - I * (k1 - k2) * b) * std::exp(I * (k1 - k2) * b) - (1.0 - I * (k1 - k2) * a) * std::exp(I * (k1 - k2) * a)) / 
			((k1 - k2) * (k1 - k2));
	}
	else {
		return 0.5 * (b * b - a * a);
	}
}

Complex Iz3(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		Complex mult1 = 2.0 * I + 2.0 * (k1 - k2) * b - I * (k1 -k2) * (k1 - k2) * b * b;
		Complex mult2 = 2.0 * I + 2.0 * (k1 - k2) * a - I * (k1 -k2) * (k1 - k2) * a * a;
		return (mult1 * std::exp(I * (k1 - k2) * b) - mult2 * std::exp(I * (k1 - k2) * a)) / 
			((k1 - k2) * (k1 - k2) * (k1 - k2));
	}
	else {
		return (b * b * b - a * a * a) / 3.0;
	}
}

Complex Iz4(Complex k1, Complex k2, double a, double b) 
{
	if (std::abs(k1 - k2) > thresh) {
		Complex mult1 = -6.0 + 6.0 * I * (k1 - k2) * b + 3.0 * (k1 -k2) * (k1 - k2) * b * b - I * std::pow((k1 - k2) * b, 3);
		Complex mult2 = -6.0 + 6.0 * I * (k1 - k2) * a + 3.0 * (k1 -k2) * (k1 - k2) * a * a - I * std::pow((k1 - k2) * a, 3);
		return (mult1 * std::exp(I * (k1 - k2) * b) - mult2 * std::exp(I * (k1 - k2) * a)) / std::pow((k1 - k2), 4);
	}
	else {
		return (std::pow(b, 4) - std::pow(a, 4)) / 4.0;
	}
}

void updateGVal(GValStruct& GVal, const LValStruct& LVal, double kCosal, Complex kn12)
{
	if (std::abs(kCosal - kn12) < thresh) {
		GVal.G1 = (I / (8.0 * kn12 * kn12 * kn12)) * LVal.L1 + (1.0 / (4.0 * kn12 * kn12)) * LVal.L2;
		GVal.G2 = (-std::exp(2.0 * I * kn12 * Lz) * (I + 2.0 * kn12 * Lz) / (8.0 * kn12 * kn12 * kn12)) * LVal.L1 -
			(std::exp(2.0 * I * kn12 * Lz) / (4.0 * kn12 * kn12)) * LVal.L2;
		GVal.G3 = (1.0 / (4.0 * kn12 * kn12)) * LVal.L1 - (I / (2.0 * kn12)) * LVal.L2;
		GVal.G4 = (-I / (4.0 * kn12)) * LVal.L1;
	}
	else if (std::abs(kCosal + kn12) < thresh) {
		GVal.G1 = (I / (8.0 * kn12 * kn12 * kn12)) * LVal.L1 - (1.0 / (4.0 * kn12 * kn12)) * LVal.L2;
		GVal.G2 = (-I * (1.0 + 2.0 * kn12 * kn12 * Lz * Lz) / (8.0 * kn12 * kn12 * kn12)) * LVal.L1 +
			((1.0 - 2.0 * Lz * kn12 * I) / (4.0 * kn12 * kn12)) * LVal.L2;
		GVal.G3 = ((1.0 + 2.0 * kn12 * Lz * I) / (4.0 * kn12 * kn12)) * LVal.L1 + (I / (2.0 * kn12)) * LVal.L2;
		GVal.G4 = (-I / (4.0 * kn12)) * LVal.L1;
	}
	else {
		GVal.G1 = (I / (2.0 * kn12 * (kCosal - kn12) * (kCosal - kn12))) * LVal.L1 + 
			(1.0 / (2.0 * kn12 * (kCosal - kn12))) * LVal.L2;
		GVal.G2 = (-2.0 * I * kCosal / ((kCosal * kCosal - kn12 * kn12) * (kCosal * kCosal - kn12 * kn12))) * LVal.L1 -
		   (1.0 / (kCosal * kCosal - kn12 * kn12)) * LVal.L2;
		GVal.G3 = (-std::exp(I * (kCosal + kn12) * Lz) * (I + (kCosal + kn12) * Lz) / (2.0 * kn12 * (kCosal + kn12) * (kCosal + kn12))) * LVal.L1 - 
			(std::exp(I * (kCosal + kn12) * Lz) / (2.0 * kn12 * (kCosal + kn12))) * LVal.L2;
		GVal.G4 = (-1.0 / (kCosal * kCosal - kn12 * kn12)) * LVal.L1;
	}
}

void updatePGVal(PGValStruct& PGVal, const GValStruct& GVal, double kCosth, double kCosal, Complex kn12)
{
	if (std::abs(kCosal - kn12) < thresh) {
		Complex Iz1kn12th = Iz1(kn12, kCosth, 0, Lz);
		Complex Iz1mkn12th = Iz1(-kn12, kCosth, 0, Lz);
		Complex Iz2kn12th = Iz2(kn12, kCosth, 0, Lz);
		Complex Iz2mkn12th = Iz2(-kn12, kCosth, 0, Lz);
		Complex Iz3kn12th = Iz3(kn12, kCosth, 0, Lz);
		Complex Iz4kn12th = Iz4(kn12, kCosth, 0, Lz);
		PGVal.PG1 = GVal.G1 * Iz1kn12th + GVal.G2 * Iz1mkn12th + 
			GVal.G3 * Iz2kn12th + GVal.G4 * Iz3kn12th;
		PGVal.PG2 = GVal.G1 * Iz2kn12th + GVal.G2 * Iz2mkn12th +
			GVal.G3 * Iz3kn12th + GVal.G4 * Iz4kn12th;
		PGVal.PG3 = I * kn12 * Iz1kn12th * GVal.G1 - I * kn12 * Iz1mkn12th * GVal.G2 + 
			(Iz1kn12th + I * kn12 * Iz2kn12th) * GVal.G3 + (2.0 * Iz2kn12th + I * kn12 * Iz3kn12th) * GVal.G4;
		PGVal.PG4 = -I * kCosth * PGVal.PG1 + PGVal.PG3;
		PGVal.PG5 = -I * kCosth * PGVal.PG2 + I * kn12 * Iz2kn12th * GVal.G1 - I * kn12 * Iz2mkn12th * GVal.G2
			+ (Iz2kn12th + I * kn12 * Iz3kn12th) * GVal.G3 + (2.0 * Iz3kn12th + I * kn12 * Iz4kn12th) * GVal.G4;
	}
	else if (std::abs(kCosal + kn12) < thresh) {
		Complex Iz1kn12th = Iz1(kn12, kCosth, 0, Lz);
		Complex Iz1mkn12th = Iz1(-kn12, kCosth, 0, Lz);
		Complex Iz2kn12th = Iz2(kn12, kCosth, 0, Lz);
		Complex Iz2mkn12th = Iz2(-kn12, kCosth, 0, Lz);
		Complex Iz3mkn12th = Iz3(-kn12, kCosth, 0, Lz);
		Complex Iz4mkn12th = Iz4(-kn12, kCosth, 0, Lz);
		PGVal.PG1 = GVal.G1 * Iz1kn12th + GVal.G2 * Iz1mkn12th + 
			GVal.G3 * Iz2mkn12th + GVal.G4 * Iz3mkn12th;
		PGVal.PG2 = GVal.G1 * Iz2kn12th + GVal.G2 * Iz2mkn12th +
			GVal.G3 * Iz3mkn12th + GVal.G4 * Iz4mkn12th;
		PGVal.PG3 = I * kn12 * Iz1kn12th * GVal.G1 - I * kn12 * Iz1mkn12th * GVal.G2 + 
			(Iz1mkn12th - I * kn12 * Iz2mkn12th) * GVal.G3 + (2.0 * Iz2mkn12th - I * kn12 * Iz3mkn12th) * GVal.G4;
		PGVal.PG4 = -I * kCosth * PGVal.PG1 + PGVal.PG3;
		PGVal.PG5 = -I * kCosth * PGVal.PG2 + I * kn12 * Iz2kn12th * GVal.G1 - I * kn12 * Iz2mkn12th * GVal.G2
			+ (Iz2mkn12th - I * kn12 * Iz3mkn12th) * GVal.G3 + (2.0 * Iz3mkn12th - I * kn12 * Iz4mkn12th) * GVal.G4;
	
	}
	else {
		Complex Iz1kn12th = Iz1(kn12, kCosth, 0, Lz);
		Complex Iz1alth = Iz1(kCosal, kCosth, 0, Lz);
		Complex Iz1mkn12th = Iz1(-kn12, kCosth, 0, Lz);
		Complex Iz2alth = Iz2(kCosal, kCosth, 0, Lz);
		Complex Iz2kn12th = Iz2(kn12, kCosth, 0, Lz);
		Complex Iz2mkn12th = Iz2(-kn12, kCosth, 0, Lz);
		Complex Iz3alth = Iz3(kCosal, kCosth, 0, Lz);
		PGVal.PG1 = GVal.G1 * Iz1kn12th + GVal.G2 * Iz1alth +
			GVal.G3 * Iz1mkn12th + GVal.G4 * Iz2alth;
		PGVal.PG2 = GVal.G1 * Iz2kn12th + GVal.G2 * Iz2alth +
			GVal.G3 * Iz2mkn12th + GVal.G4 * Iz3alth;
		PGVal.PG3 = I * kn12 * Iz1kn12th * GVal.G1 + I * kCosal * Iz1alth * GVal.G2 -
			I * kn12 * Iz1mkn12th * GVal.G3 + (Iz1alth + I * kCosal * Iz2alth) * GVal.G4;
		PGVal.PG4 = -I * kCosth * PGVal.PG1 + PGVal.PG3;
		PGVal.PG5 = -I * kCosth * PGVal.PG2 + I * kn12 * Iz2kn12th * GVal.G1 + I * kCosal * Iz2alth * GVal.G2
			- I * kn12 * Iz2mkn12th * GVal.G3 + (Iz2alth + I * kCosal * Iz3alth) * GVal.G4;
	}
}
