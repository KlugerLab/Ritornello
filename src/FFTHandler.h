/*
 * FFTHandler.h
 *
 *  Created on: Aug 12, 2014
 *      Author: Kelly Patrick Stanton
 */

#ifndef FFTHANDLER_H_
#define FFTHANDLER_H_
#include <fftw3.h>
#include <complex>
#include <vector>

#define LOWPASSFREQ 0.04
using namespace std;

class FFTHandler {
public:
	static void init(long n);
	static void destroy();
	static void forward(double* inout);
	static void reverse(double* inout);
	static double* newArray();
	static void freeArray(double* argArray);
	static double* copy(double* in);
	static void complexInPlaceDivide(fftw_complex* inout, fftw_complex* in2);
	static void complexInPlaceMultiply(fftw_complex* inout, fftw_complex* in2);
	static void complexInPlaceConjugate(fftw_complex* inout);
	static void inPlaceConvolve(double* inout, double* in2);
	static void inPlaceDeconvolve(double* inout, double* in2);
	static double* getGaussianKernel(double var);
	//double* convolve(const double* f, const double* g) const;
	//double* deconvolve(const double* f, const double* g) const;
private:

	static fftw_plan fftForwardPlan;
	static fftw_plan fftReversePlan;
	static long leased;
	static long n;
	static vector<double*> memoryPool;
};

#endif /* FFTHANDLER_H_ */
