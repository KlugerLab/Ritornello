/*
 * FFTHandler.cpp
 *
 *  Created on: Aug 12, 2014
 *      Author: Kelly Patrick Stanton
 */

#include "FFTHandler.h"
#include <math.h>
#include <omp.h>

fftw_plan FFTHandler::fftForwardPlan = NULL;
fftw_plan FFTHandler::fftReversePlan = NULL;
long FFTHandler::n = 0;
long FFTHandler::leased=0;
vector<double*> FFTHandler::memoryPool;

void FFTHandler::init(long arg_n){
//#ifndef DEBUG
	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());
//#endif
	n = arg_n;
	leased=0;
	fprintf(stderr, "Initializing fft plan of size [%ld]\n", n);
	memoryPool.resize(1);
	memoryPool[0] = (double*)fftw_malloc(sizeof(double)*2*(n/2+1));
	fftForwardPlan = fftw_plan_dft_r2c_1d(n, memoryPool[0], (fftw_complex*)memoryPool[0],FFTW_MEASURE);
	fftReversePlan = fftw_plan_dft_c2r_1d(n, (fftw_complex*)memoryPool[0], memoryPool[0],FFTW_MEASURE);
}
void FFTHandler::destroy(){
	if(fftForwardPlan!=NULL){
		fftw_destroy_plan(fftForwardPlan);
		fftForwardPlan=NULL;
	}
	if(fftReversePlan!=NULL){
		fftw_destroy_plan(fftReversePlan);
		fftReversePlan=NULL;
	}

	for(long ii =0; ii < (long)memoryPool.size();++ii)
		fftw_free(memoryPool[ii]);
	memoryPool.resize(0);
}

double* FFTHandler::newArray(){
	++leased;
	if(leased>(long)memoryPool.size()){
		memoryPool.push_back((double*)fftw_malloc(sizeof(double)*2*(n/2+1)));
	}
	if(leased> 11)
		fprintf(stderr, "Requested more than 11 genome sized arrays\n");
	for(long ii =0; ii < n; ++ii)
		memoryPool[leased-1][ii]=0;
	return memoryPool[leased-1];
}
void FFTHandler::freeArray(double* argArray){
	long ii=0;
	for(ii =0; ii < leased; ++ii){
		if(memoryPool[ii]==argArray)
			break;
	}
	if(ii>=leased){
		fprintf(stderr, "Error freeing fftarray\n");
	}
	//swap the array with the end of the pool
	double* tmp = memoryPool[ii];
	memoryPool[ii] = memoryPool[leased-1];
	memoryPool[leased-1] = tmp;
	--leased;
}
void FFTHandler::forward(double* inout){
	fftw_execute_dft_r2c(fftForwardPlan,inout, (fftw_complex*)inout);
}
void FFTHandler::reverse(double* inout){
	fftw_execute_dft_c2r(fftReversePlan,(fftw_complex*)inout, inout);
	//normalize the reverse transform
	for(long ii =0; ii < n; ++ii) inout[ii]/=n;
}
double* FFTHandler::copy(double* in){
	double* ret = newArray();
	for(int ii = 0; ii < 2*(n/2+1); ++ii){
		ret[ii]=in[ii];
	}
	return ret;
}
void FFTHandler::complexInPlaceDivide(fftw_complex* inout, fftw_complex* in2){
	for(long ii=0;ii<(n/2+1);++ii) {
		if(ii<(n/2+1)*LOWPASSFREQ){
			complex<double> tmp=complex<double>(inout[ii][0],inout[ii][1])
					/complex<double>(in2[ii][0],in2[ii][1]);
			inout[ii][0]=tmp.real();
			inout[ii][1]=tmp.imag();
		}
		else{
			inout[ii][0]=0;
			inout[ii][1]=0;
		}
	}
}
void FFTHandler::complexInPlaceMultiply(fftw_complex* inout, fftw_complex* in2){
	for(long ii=0;ii<(n/2+1);++ii) {
		complex<double> tmp=complex<double>(inout[ii][0],inout[ii][1])
				*complex<double>(in2[ii][0],in2[ii][1]);
		inout[ii][0]=tmp.real();
		inout[ii][1]=tmp.imag();
	}
}
void FFTHandler::inPlaceConvolve(double* inout, double* in2){
	forward(inout);
	forward(in2);
	complexInPlaceMultiply((fftw_complex*)inout, (fftw_complex*)in2);
	reverse(inout);
	reverse(in2);
}
void FFTHandler::inPlaceDeconvolve(double* inout, double* in2){
	forward(inout);
	forward(in2);
	complexInPlaceDivide((fftw_complex*)inout, (fftw_complex*)in2);
	reverse(inout);
	reverse(in2);
}
void FFTHandler::complexInPlaceConjugate(fftw_complex* inout){
	for(long ii=0;ii<(n/2+1);++ii) {
		inout[ii][1]*=-1;
	}
}
double* FFTHandler::getGaussianKernel(double var){
	//create Guassian kernel smoother
	double* Gaussian = FFTHandler::newArray();
	for(long ii = 0; ii < n; ++ii){
		Gaussian[ii] = 1/sqrt(var*2*M_PI)*exp(-pow(ii,2.0)/(2*var))
				+1/sqrt(var*2*M_PI)*exp(-pow(ii-(n),2.0)/(2*var));
	}
	return Gaussian;
}
