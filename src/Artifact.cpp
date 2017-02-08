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

#define FOURIERLENGTH 512

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
//returns 0 if this is an artifact, 1 otherwise
double Artifact::test(long windowSize, long readLength, double* fir, Peak candidate, double testRatio, long halfLength){
	//Reserve storage for positive and negative strands
	double* ps = FFTHandler::newArray();
	double* ms = FFTHandler::newArray();

	//copy and shift by read length
	for(int ii = 0; ii < FOURIERLENGTH+readLength; ++ii){
		//TODO revisit this.  It may be better to just enforce larger windows
		long pind = 2*windowSize-FOURIERLENGTH/2 + ii;//- readLength;
		long mind = 2*windowSize-FOURIERLENGTH/2 + ii;
		if(pind >= 0 && pind < 4*windowSize)
			ps[ii]=candidate.pstrandExtended[pind];
		else
			ps[ii]=0;
		if(mind >= 0 && mind < 4*windowSize)
			ms[ii]=candidate.mstrandExtended[mind];
		else
			ms[ii]=0;
	}
	for(int ii = FOURIERLENGTH+readLength; ii < 2*FOURIERLENGTH; ++ii){
		ps[ii]=0;
		ms[ii]=0;
	}

	//IOhandler::printDoubleArrayToFile(ps,2*FOURIERLENGTH,"ps.txt");
	//IOhandler::printDoubleArrayToFile(ms,2*FOURIERLENGTH,"ms.txt");

	//calculate the cross correlation
	FFTHandler::forward(ps);
	FFTHandler::complexInPlaceConjugate((fftw_complex*)ps);
	FFTHandler::forward(ms);
	FFTHandler::complexInPlaceMultiply((fftw_complex*)ms,(fftw_complex*)ps);
	FFTHandler::reverse(ms);

	IOhandler::printDoubleArrayToFile(ms,2*FOURIERLENGTH,"cc.txt");

	//make ms sum to 1
	double sum = 0;
	for(long ii=0;ii<2*FOURIERLENGTH;++ii)
		sum+=ms[ii];
	for(long ii=0;ii<2*FOURIERLENGTH;++ii)
		ms[ii]/=sum;

	IOhandler::printDoubleArrayToFile(ms,2*FOURIERLENGTH,"normcc.txt");

	double score1 = 0;
	for(long ii=0;ii<readLength;++ii)
			if(ms[ii]>score1)
				score1 = ms[ii];
	double score2 = 0;
	for(long ii=readLength+10;ii<min(windowSize,(long)2*FOURIERLENGTH);++ii)
			if(ms[ii]>score2)
				score2 = ms[ii];
	score2+=1;

	//copy and shift by read length
	for(int ii = 0; ii < FOURIERLENGTH+readLength; ++ii){
		//TODO revisit this.  It may be better to just enforce larger windows
		long pind = 2*windowSize-FOURIERLENGTH/2 + ii;//- readLength;
		long mind = 2*windowSize-FOURIERLENGTH/2 + ii;
		if(pind >= 0 && pind < 4*windowSize)
			ps[ii]=candidate.pstrandExtended[pind]>0 ? 1:0;
		else
			ps[ii]=0;
		if(mind >= 0 && mind < 4*windowSize)
			ms[ii]=candidate.mstrandExtended[mind]>0 ? 1:0;
		else
			ms[ii]=0;
	}
	for(int ii = FOURIERLENGTH+readLength; ii < 2*FOURIERLENGTH; ++ii){
		ps[ii]=0;
		ms[ii]=0;
	}

	IOhandler::printDoubleArrayToFile(ps,2*FOURIERLENGTH,"ps.txt");
	IOhandler::printDoubleArrayToFile(ms,2*FOURIERLENGTH,"ms.txt");

	//calculate the cross correlation
	FFTHandler::forward(ps);
	FFTHandler::complexInPlaceConjugate((fftw_complex*)ps);
	FFTHandler::forward(ms);
	FFTHandler::complexInPlaceMultiply((fftw_complex*)ms,(fftw_complex*)ps);
	FFTHandler::reverse(ms);

	IOhandler::printDoubleArrayToFile(ms,2*FOURIERLENGTH,"bincc.txt");
	double score3 = 0;
	for(long ii=0;ii<readLength;++ii)
			if(ms[ii]>score3)
				score3 = ms[ii];
	double score4 = 0;
	for(long ii=readLength+10;ii<min(windowSize,(long)2*FOURIERLENGTH);++ii)
			if(ms[ii]>score4)
				score4 = ms[ii];
	score4+=1;
	double feature1 = score1/score2;
	double feature2 = score3/score4;
	double activation = 28.70394-637.37217*feature1-26.13757*feature2;
	/*
	//go back to Fourier space
	FFTHandler::forward(ms);

	//IOhandler::printDoubleArrayToFile(ms,2*FOURIERLENGTH,"ccfft.txt");
	//Get the features
	int nfreqs = 10;
	double features[2*nfreqs];
	for(long ii=0;ii<nfreqs;++ii) {
		features[ii]=ms[2*ii];
		features[ii+nfreqs]=ms[2*ii+1];
	}

	//IOhandler::printDoubleArrayToFile(features,100,"features.txt");

	//Inner product the features with the coefficients and add the offset
	double activation = b;
	for(long ii=0;ii<2*nfreqs;++ii) {
		activation+=features[ii]*W[ii];
	}


	double rlscore = 0;
	for(int ii = 0; ii < readLength; ++ii)
		rlscore+=ms[ii];
	double flscore = 0;
	for(int ii = 2*readLength; ii <500; ++ii)
		flscore+=ms[ii];
	//double total = rlscore + flscore;
	//for(int ii = 2*readLength; ii < 500; ++ii)
		//total+=ms[ii];
	//rlscore/=total;
	//flscore/=total;
	double activation = -10.81017+-14.23095*rlscore+27.87288 *flscore;


	//features are the sum of each bin
	long binsize = 16;
	long numBins = 2*FOURIERLENGTH/binsize;
	double features[numBins];
	for(int ii = 0; ii < numBins; ++ii)
		features[ii]=0;
	for(int ii = 0; ii < 2*FOURIERLENGTH; ++ii){
		features[ii/binsize]+=ms[ii];
	}
	//Inner product the features with the coefficients and add the offset
	double activation = b;
	for(long ii=0;ii<numBins;++ii) {
		activation+=features[ii]*W[ii];
	}
	IOhandler::printDoubleArrayToFile(features,32,"features.txt");
	 */

	FFTHandler::freeArray(ps);
	FFTHandler::freeArray(ms);

	//apply logit function
	double logitScore = 1/(1+exp(-activation));
	// return sigmoid activation
	return logitScore;
	//return 1;
}
