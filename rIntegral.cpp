#include "rIntegral.h"

double powm1(int n)
{
	if (n % 2) {
		return -1.0;
	}
	return 1.0;
}

Complex I10(int n, const double* HermiteMat, double alpha)
{
	return HermiteMat[n];
}

Complex I11(int n, const double* HermiteMat, double alpha)
{
	if (n > 0) {
		return ((0.5 * HermiteMat[n + 1] - n * HermiteMat[n - 1]) / alpha) * I;
	}
	else {
		return (0.5 * HermiteMat[n + 1] / alpha) * I;
	}
	return 0.0;
}

Complex I12(int n, const double* HermiteMat, double alpha)
{
	if (n > 1) {
		return -((0.25 * HermiteMat[n + 2] - (n + 0.5) * HermiteMat[n] + n * (n - 1) * HermiteMat[n - 2]) / (alpha * alpha));
	}
	else {
		return -((0.25 * HermiteMat[n + 2] - (n + 0.5) * HermiteMat[n]) / (alpha * alpha));
	}
	return 0.0;
}

Complex I13(int n, const double* HermiteMat, double alpha)
{
	if (n > 2) {
		return -((0.125 * HermiteMat[n + 3] - 0.75 * (n + 1) * HermiteMat[n + 1] + 1.5 * n * n * HermiteMat[n - 1] - n * (n - 1) * (n - 2) * HermiteMat[n - 3]) / (alpha * alpha * alpha)) * I;
	}
	else if (n > 0) {
		return -((0.125 * HermiteMat[n + 3] - 0.75 * (n + 1) * HermiteMat[n + 1] + 1.5 * n * n * HermiteMat[n - 1]) / (alpha * alpha * alpha)) * I;
	}
	else {
		return -((0.125 * HermiteMat[n + 3] - 0.75 * (n + 1) * HermiteMat[n + 1]) / (alpha * alpha * alpha)) * I;
	}

	return 0.0;
}

Complex I20(int n, double x, double alpha)
{
	return std::pow(x / alpha, n) / sqrt(2);
}
