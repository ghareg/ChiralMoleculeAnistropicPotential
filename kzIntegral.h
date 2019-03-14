#ifndef KZINTEGRAL_H_
#define KZINTEGRAL_H_

#include "LDefinit.h"

Complex Iz1(Complex k1, Complex k2, double a, double b);
Complex Iz2(Complex k1, Complex k2, double a, double b);
Complex Iz3(Complex k1, Complex k2, double a, double b);
Complex Iz4(Complex k1, Complex k2, double a, double b);
void updateGVal(GValStruct& GVal, const LValStruct& LVal, double kCosal, Complex kn12);
void updatePGVal(PGValStruct& PGVal, const GValStruct& GVal, double kCosth, double kCosal, Complex kn12);

#endif
