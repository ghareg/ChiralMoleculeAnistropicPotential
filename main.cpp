#include "kzIntegral.h"
#include "rIntegral.h"
#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_hermite.h>
#include <thread>
void calcFullPol(double* Pol, const double* fact, const PauliMatrix& pal);
void calcPol(double* Pol, const double* fact, const PauliMatrix& pal, int rs, int re);
void generateHermite(int nmax, double arg, double* HMat);

int main()
{
	PauliMatrix pal;
	initPauliMatrix(pal);

	double* Pol = new double[ND * ND]; 
	double fact[nmax + 1];
	for (int k = 0; k < nmax + 1; ++k) {
		if (k == 0) {
			fact[k] = 1.0;
		}
		else {
			fact[k] = k * fact[k - 1];
		}
	}

	calcFullPol(Pol, fact, pal);
	double E = 0.001;
	double theta = 0.0;
	std::cout << std::setprecision(4);
	for (int ind1 = 0; ind1 < ND; ++ind1) {
		E = 0.001;
		for (int ind2 = 0; ind2 < ND; ++ind2) {
			std::cout << theta << "\t" << E << "\t" << Pol[ind2 * ND + ind1] << std::endl;
			E += EMax / ND;
		}
		std::cout << std::endl;
		theta += Pi / ND;
	}
	delete[] Pol;
}

void calcFullPol(double* Pol, const double* fact, const PauliMatrix& pal)
{
	double* PolMat[nCore];
	int div = ND / nCore;
	int size1 = div * ND;
	int size2 = (ND - (nCore - 1) * div) * ND;
	for (int i = 0; i < nCore - 1; ++i) {
		PolMat[i] = new double[size1];
	}
	PolMat[nCore - 1] = new double[size2];
	std::thread t[nCore - 1];
	for (int i = 0; i < nCore - 1; ++i) {
		t[i] = std::thread(calcPol, PolMat[i], fact, std::cref(pal), i * div, (i + 1) * div);
	}
						
	calcPol(PolMat[nCore - 1], fact, pal, (nCore - 1) * div, ND);

	for (int i = 0; i < nCore - 1; ++i) {
		t[i].join();
	}

	int ind = 0;
	for (int i = 0; i < nCore - 1; ++i) {
		for (int j = 0; j < size1; ++j) {
			Pol[ind++] = PolMat[i][j];
		}
	}
	for (int j = 0; j < size2; ++j) {
		Pol[ind++] = PolMat[nCore - 1][j];
	}
	for (int i = 0; i < nCore; ++i) {
		delete[] PolMat[i];
	}
}

void calcPol(double* Pol, const double* fact, const PauliMatrix& pal, int rs, int re)
{
	double E = EnMult * (rs * EMax / ND + 0.001);
	double k = sqrt(E);
	double alpha = 0.8;
	double beta = 0.0;
	double theta = 0;
	double tau = 0;
	Complex kn12 = 0.0;
	double qn12 = 0.0;
	double alphaxSq = 0.0;
	double alphaySq = 0.0;
	Matrix2cd ft;
	Matrix2cd f0;
	Matrix2cd f1;
	Matrix2cd f2;
	Matrix2cd f3;
	Matrix2cd f4;
	Matrix2cd rho = 0.5 * pal.pal0;
	
	Complex I10xt;
	Complex I11xt;
	Complex I12xt;
	Complex I13xt;
	Complex I10yt;
	Complex I11yt;
	Complex I12yt;
	Complex I13yt;

	Complex I10m1xt;
	Complex I11m1xt;
	Complex I12m1xt;
	Complex I10m1yt;
	Complex I11m1yt;
	Complex I12m1yt;
	
	Complex multx;
	Complex multy;
	LValStruct LVal;
	GValStruct GVal;
	PGValStruct PGVal;

	double tmultx = 0.0;
	double tmulty = 0.0;

	double HermiteMatx1[nmax + 4];
	double HermiteMatx2[nmax + 4];
	double HermiteMaty1[nmax + 4];
	double HermiteMaty2[nmax + 4];
	Param pm;
	double Polc = 0.0;

	for (int ind1 = 0; ind1 < (re - rs); ++ind1) {
		k = sqrt(E);
		theta = 0.0;
		for (int ind2 = 0; ind2 < ND; ++ind2) {
			tau = 0.0;
			Polc = 0.0;	
			for (int ind3 = 0; ind3 < ND; ++ind3) {
				beta = 0.0;
				for (int ind4 = 0; ind4 < ND; ++ind4) {
					pm.updateValues(k, alpha, beta, theta, tau, alphaxC, alphayC);
					f0 << 0.0, 0.0, 0.0, 0.0;
					f1 << 0.0, 0.0, 0.0, 0.0;
					f2 << 0.0, 0.0, 0.0, 0.0;
					f3 << 0.0, 0.0, 0.0, 0.0;
					f4 << 0.0, 0.0, 0.0, 0.0;
					generateHermite(nmax + 4, pm.kSinal * pm.cosBeta / pm.alphax, HermiteMatx1);
					generateHermite(nmax + 4, pm.kSinal * pm.sinBeta / pm.alphay, HermiteMaty1);
				//	gsl_sf_hermite_phys_array(nmax + 4, pm.kSinal * pm.cosBeta / pm.alphax, HermiteMatx1);
				//	gsl_sf_hermite_phys_array(nmax + 4, pm.kSinal * pm.sinBeta / pm.alphay, HermiteMaty1);
					
					generateHermite(nmax + 4, -pm.kSinth * pm.cosTau / pm.alphax, HermiteMatx2);
					generateHermite(nmax + 4, -pm.kSinth * pm.sinTau / pm.alphay, HermiteMaty2);
				//	gsl_sf_hermite_phys_array(nmax + 4, -pm.kSinth * pm.cosTau / pm.alphax, HermiteMatx2);
				//	gsl_sf_hermite_phys_array(nmax + 4, -pm.kSinth * pm.sinTau / pm.alphay, HermiteMaty2);
					for (int nc1 = 0; nc1 <= nmax; ++nc1) {
						for (int nc2 = 0; nc2 <= nmax; ++nc2) {	
							alphaxSq = pm.alphax * pm.alphax;
							alphaySq = pm.alphay * pm.alphay;	
							qn12 = alphaxSq * (2.0 * nc1 + 1) + alphaySq * (2.0 * nc2 + 1);
							kn12 = std::sqrt(Complex(k * k - qn12, 0.0));
							updateLVal(LVal, nc1, nc2, pm, pal, fact, HermiteMatx1, HermiteMaty1);
							updateGVal(GVal, LVal, pm.kCosal, kn12);
							updatePGVal(PGVal, GVal, pm.kCosth, pm.kCosal, kn12);
							
							multx = (pi14 * std::pow(I , nc1) / std::sqrt(pm.alphax * std::pow(2, nc1 - 1) * fact[nc1])) *
								std::exp(-0.5 * pm.kSinth * pm.kSinth * pm.cosTau * pm.cosTau / (alphaxSq));
							multy = (pi14 * std::pow(I , nc2) / std::sqrt(pm.alphay * std::pow(2, nc2 - 1) * fact[nc2])) *
								std::exp(-0.5 * pm.kSinth * pm.kSinth * pm.sinTau * pm.sinTau / (alphaySq));
							
							I10xt = multx * I10(nc1, HermiteMatx2, pm.alphax);
							I11xt = multx * I11(nc1, HermiteMatx2, pm.alphax);
							I12xt = multx * I12(nc1, HermiteMatx2, pm.alphax);
							I13xt = multx * I13(nc1, HermiteMatx2, pm.alphax);
							I10yt = multy * I10(nc2, HermiteMaty2, pm.alphay);
							I11yt = multy * I11(nc2, HermiteMaty2, pm.alphay);
							I12yt = multy * I12(nc2, HermiteMaty2, pm.alphay);
							I13yt = multy * I13(nc2, HermiteMaty2, pm.alphay);

							if (nc1 > 0) {
								multx = (pi14 * std::pow(I , nc1 - 1) / std::sqrt(pm.alphax * std::pow(2, nc1 - 2) * fact[nc1 - 1])) *
									std::exp(-0.5 * pm.kSinth * pm.kSinth * pm.cosTau * pm.cosTau / (alphaxSq));
								I10m1xt = multx * I10(nc1 - 1, HermiteMatx2, pm.alphax);
								I11m1xt = multx * I11(nc1 - 1, HermiteMatx2, pm.alphax);
								I12m1xt = multx * I12(nc1 - 1, HermiteMatx2, pm.alphax);
							}
							else {
								I10m1xt = 0.0;
								I11m1xt = 0.0;
								I12m1xt = 0.0;
							}
							
							if (nc2 > 0) {
								multy = (pi14 * std::pow(I , nc2 - 1) / std::sqrt(pm.alphay * std::pow(2, nc2 - 2) * fact[nc2 - 1])) *
									std::exp(-0.5 * pm.kSinth * pm.kSinth * pm.sinTau * pm.sinTau / (alphaySq));
								I10m1yt = multy * I10(nc2 - 1, HermiteMaty2, pm.alphay);
								I11m1yt = multy * I11(nc2 - 1, HermiteMaty2, pm.alphay);
								I12m1yt = multy * I12(nc2 - 1, HermiteMaty2, pm.alphay);
							}
							else {
								I10m1yt = 0.0;
								I11m1yt = 0.0;
								I12m1yt = 0.0;
							}

							tmultx = std::sqrt(2 * nc1) * pm.alphax;
							tmulty = std::sqrt(2 * nc2) * pm.alphay;
							
							f0 += Iz2(pm.kCosal, pm.kCosth, 0, Lz) * I10xt * I10yt * LVal.L1 +
								Iz1(pm.kCosal, pm.kCosth, 0, Lz) * I10xt * I10yt * LVal.L2;
							f1 += F1 * I10xt * I10yt * PGVal.PG2 + (F2 * I11xt * I10yt + F3 * I10xt * I11yt) * PGVal.PG1;
							f2 += (2.0 * alphaySq * alphaySq * I11yt + F3 * I10yt) * I10xt * PGVal.PG3 + 
								(alphaxSq * alphaxSq * I12xt * (-alphaySq * I11yt + tmulty * I10m1yt) + 
								 alphaySq * alphaySq * I10xt * (-alphaySq * I13yt + tmulty * I12m1yt) + 
								 F2 * I11xt * (-alphaySq * I11yt + tmulty * I10m1yt) +
								 F3 * I10xt * (-alphaySq * I12yt + tmulty * I11m1yt)) * PGVal.PG4 +
								F1 * I10xt * (-alphaySq * I11yt + tmulty * I10m1yt) * PGVal.PG5;
							f3 += (2.0 * alphaxSq * alphaxSq * I11xt + F2 * I10xt) * I10yt * PGVal.PG3 + 
								(alphaxSq * alphaxSq * I10yt * (-alphaxSq * I13xt + tmultx * I12m1xt) + 
								 alphaySq * alphaySq * I12yt * (-alphaxSq * I11xt + tmultx * I10m1xt) + 
								 F2 * I10yt * (-alphaxSq * I12xt + tmultx * I11m1xt) +
								 F3 * I11yt * (-alphaxSq * I11xt + tmultx * I10m1xt)) * PGVal.PG4 +
								F1 * I10yt * (-alphaxSq * I11xt + tmultx * I10m1xt) * PGVal.PG5;
							f4 += ((2.0 * alphaxSq * alphaxSq * I11xt + F2 * I10xt) * (-alphaySq * I11yt + tmulty * I10m1yt) -
									(2.0 * alphaySq * alphaySq * I11yt + F3 * I10yt) * (-alphaxSq * I11xt + tmultx * I10m1xt)) 
								* PGVal.PG1;
						}
					}
					f0 = (-1.0  / (4 * Pi)) * f0;
					f1 = (-1.0 / (4 * Pi)) * f1;
					f2 = (alphaSOC / (4 * Pi)) * I * (pal.pal1 * f2);
					f3 = (-alphaSOC / (4 * Pi)) * I * pal.pal2 * f3;
					f4 = (alphaSOC / (4 * Pi)) * I * pal.pal3 * f4;
					ft = f0 + f1 + f2 + f3 + f4;
					ft = ft * rho * ft.adjoint();
					Polc += ((pal.pal3 * ft).trace() / ft.trace()).real();
					beta += 2 * Pi / ND;
				}
				tau += 2 * Pi / ND;
			}
			Pol[ind1 * ND + ind2] = Polc / (1.0 * ND * ND);
			theta += Pi / ND;
		}
		E += EnMult * EMax / ND;
	}
}

void generateHermite(int nmax, double arg, double* HMat)
{
	if (nmax > 0) {
		HMat[0] = 1.0;
	}
	if (nmax > 1) {
		HMat[1] = 2.0 * arg;
	}
	for (int i = 2; i < nmax; ++i) {
		HMat[i] = 2.0 * arg * HMat[i - 1] - 2.0 * (i - 1) * HMat[i - 2];
	}
}
