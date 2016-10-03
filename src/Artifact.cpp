/*
 * Artifact.cpp
 *
 *  Created on: May 5, 2015
 *      Author: Kelly
 */

#include "Artifact.h"
#include "Peaks.h"
#include <math.h>
#include "IOhandler.h"
#include <algorithm>
#include "FFTHandler.h"

#define FOURIERLENGTH 1009
Artifact::Artifact() {
	// TODO Auto-generated constructor stub

}

Artifact::~Artifact() {
	// TODO Auto-generated destructor stub
}

void Artifact::init(){
	FFTHandler::init(2*FOURIERLENGTH);
}
void Artifact::destroy(){
	FFTHandler::destroy();
}
double Artifact::test(long windowSize, long readLength, double* fir, Peak candidate, double testRatio, long halfLength){
	//Reserve storage for positive and negative strands
	double* ps = FFTHandler::newArray();
	double* ms = FFTHandler::newArray();

	//copy and shift by read length
	for(int ii = 0; ii < FOURIERLENGTH; ++ii){
		//TODO revisit this.  It may be better to just enforce larger windows
		long pind = 2*windowSize-FOURIERLENGTH/2 + ii - readLength/2;
		long mind = 2*windowSize-FOURIERLENGTH/2 + ii + readLength/2;
		if(pind >= 0 && pind < 4*windowSize)
			ps[ii]=candidate.pstrandExtended[pind];
		else
			ps[ii]=0;
		if(mind >= 0 && mind < 4*windowSize)
			ms[ii]=candidate.mstrandExtended[mind];
		else
			ms[ii]=0;
	}
	for(int ii = FOURIERLENGTH; ii < 2*FOURIERLENGTH; ++ii){
		ps[ii]=0;
		ms[ii]=0;
	}

	//calculate the fourier transform of the cross correlation
	FFTHandler::forward(ps);
	FFTHandler::complexInPlaceConjugate((fftw_complex*)ps);
	FFTHandler::forward(ms);
	double features[200];
	for(long ii=0;ii<100;++ii) {
		//Do the complex product of the minus strand and the reversed positive strand in the frequency domain
		complex<double> tmp=complex<double>(((fftw_complex*)ms)[ii][0],((fftw_complex*)ms)[ii][1])
				*complex<double>(((fftw_complex*)ps)[ii][0],((fftw_complex*)ps)[ii][1]);
		features[ii]=tmp.real();
		features[ii+100]=tmp.imag();
	}
	FFTHandler::freeArray(ps);
	FFTHandler::freeArray(ms);
	//Inner product the features with the coefficients and add the offset
	double activation = b;
	for(long ii=0;ii<200;++ii) {
		activation+=features[ii]*W[ii];
	}

	//apply logit function
	double logitScore = 1/(1+exp(-activation));
	// return sigmoid activation
	return logitScore;
}
