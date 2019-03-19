#include "rIntegral.h"
#include "LDefinit.h"

void initPauliMatrix(PauliMatrix& pal)
{
	pal.pal0 << 1.0, 0.0, 0.0, 1.0;
	pal.pal1 << 0.0, 1.0, 1.0, 0.0;
	pal.pal2 << 0.0, -I, I, 0.0;
	pal.pal3 << 1.0, 0.0, 0.0, -1.0;
}

void updateLVal(LValStruct& LVal, int n1, int n2, const Param& pm, const PauliMatrix& pal, const double* fact, const double* HermiteMatx, const double* HermiteMaty)
{
	Complex multx;
	Complex multy;
	multx = (pi14 * std::pow(I , n1) / std::sqrt(pm.alphax * std::pow(2, n1 - 1) * fact[n1])) * 
		std::exp(-0.5 * pm.kSinal * pm.kSinal * pm.cosBeta * pm.cosBeta / (pm.alphax * pm.alphax)); 
	multy = (pi14 * std::pow(I , n2) / std::sqrt(pm.alphay * std::pow(2, n2 - 1) * fact[n2])) * 
		std::exp(-0.5 * pm.kSinal * pm.kSinal * pm.sinBeta * pm.sinBeta / (pm.alphay * pm.alphay)); 
	
	Complex I10x = multx * I10(n1, HermiteMatx, pm.alphax);
	Complex I10y = multy * I10(n2, HermiteMaty, pm.alphay);
	Complex I11x = multx * I11(n1, HermiteMatx, pm.alphax);
	Complex I11y = multy * I11(n2, HermiteMaty, pm.alphay);

	Complex I20x = multx * I20(n1, pm.kSinal * pm.cosBeta, pm.alphax);
	Complex I20y = multy * I20(n2, pm.kSinal * pm.sinBeta, pm.alphay);

	double alphax4 = std::pow(pm.alphax, 4);
	double alphay4 = std::pow(pm.alphay, 4);

	LVal.L1 = F1 * I10x * I10y * pal.pal0;
	LVal.L2 = I * kappa * I20x * I20y * pal.pal0 + (F2 * I11x * I10y + F3 * I10x * I11y) * pal.pal0 + alphaSOC1 * I10x * (2.0 * alphay4 * pm.kCosal * I11y + (F3 * pm.kCosal - F1 * pm.kSinal * pm.sinBeta) * I10y)  * pal.pal1 -
		alphaSOC2 * I10y * (2.0 * alphax4 * pm.kCosal * I11x + (F2 * pm.kCosal - F1 * pm.kSinal * pm.cosBeta) * I10x) * pal.pal2 + 
		alphaSOC3 * (2.0 * pm.kSinal * (alphax4 * pm.sinBeta * I11x * I10y - alphay4 * pm.cosBeta * I10x * I11y) + (F2 * pm.kSinal * pm.sinBeta - F3 * pm.kSinal * pm.cosBeta) * I10x * I10y) * pal.pal3;
}

void Param::updateValues(double kc, double alphac, double betac, double thetac, double tauc, double alphaxc, double alphayc)
{
	k = kc;
	alpha = alphac;
	beta = betac;
	theta = thetac;
	tau = tauc;
	kSinal = k * std::sin(alpha);
	kCosal = k * std::cos(alpha);
	kSinth = k * std::sin(theta);
	kCosth = k * std::cos(theta);
	sinBeta = std::sin(beta);
	cosBeta = std::cos(beta);
	sinTau = std::sin(tau);
	cosTau = std::cos(tau);
	alphax = alphaxc;
	alphay = alphayc;
}
