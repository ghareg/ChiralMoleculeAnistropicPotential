#include "kzIntegral.h"
#include "rIntegral.h"
#include <iostream>
#include <iomanip>
//#include <gsl/gsl_sf_hermite.h>
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
	double alphaxSq = 0.0;
	double alphaySq = 0.0;
	Matrix2cd ft;
	Matrix2cd f0;
	Matrix2cd rho = 0.5 * pal.pal0;
	
	Complex I10xt;
	Complex I10yt;

	Complex multx;
	Complex multy;
	LValStruct LVal;


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
			//	Polc = 0.0;	
				beta = 0.0;
				for (int ind4 = 0; ind4 < ND; ++ind4) {
					pm.updateValues(k, alpha, beta, theta, tau, alphaxC, alphayC);
					f0 << 0.0, 0.0, 0.0, 0.0;
					generateHermite(nmax + 4, pm.kSinal * pm.cosBeta / pm.alphax, HermiteMatx1);
					generateHermite(nmax + 4, pm.kSinal * pm.sinBeta / pm.alphay, HermiteMaty1);
					
					generateHermite(nmax + 4, -pm.kSinth * pm.cosTau / pm.alphax, HermiteMatx2);
					generateHermite(nmax + 4, -pm.kSinth * pm.sinTau / pm.alphay, HermiteMaty2);
					for (int nc1 = 0; nc1 <= nmax; ++nc1) {
						for (int nc2 = 0; nc2 <= nmax; ++nc2) {	
							alphaxSq = pm.alphax * pm.alphax;
							alphaySq = pm.alphay * pm.alphay;	
							updateLVal(LVal, nc1, nc2, pm, pal, fact, HermiteMatx1, HermiteMaty1);
							
							multx = (pi14 * std::pow(I , nc1) / std::sqrt(pm.alphax * std::pow(2, nc1 - 1) * fact[nc1])) *
								std::exp(-0.5 * pm.kSinth * pm.kSinth * pm.cosTau * pm.cosTau / (alphaxSq));
							multy = (pi14 * std::pow(I , nc2) / std::sqrt(pm.alphay * std::pow(2, nc2 - 1) * fact[nc2])) *
								std::exp(-0.5 * pm.kSinth * pm.kSinth * pm.sinTau * pm.sinTau / (alphaySq));
							
							I10xt = multx * I10(nc1, HermiteMatx2, pm.alphax);
							I10yt = multy * I10(nc2, HermiteMaty2, pm.alphay);
		
							f0 += Iz2(pm.kCosal, pm.kCosth, 0, Lz) * I10xt * I10yt * LVal.L1 +
								Iz1(pm.kCosal, pm.kCosth, 0, Lz) * I10xt * I10yt * LVal.L2;
						}
					}
					f0 = (-1.0  / (4 * Pi)) * f0;
					ft = f0;
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
