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

Artifact::Artifact() {
	// TODO Auto-generated constructor stub

}

Artifact::~Artifact() {
	// TODO Auto-generated destructor stub
}

double Artifact::test(long windowSize, long readLength, double* fir, Peak candidate, double testRatio, long halfLength){
	int ccLength = 0.75*windowSize;
	long smoothingbandwidth = 2;
	long coverageSize = 2*windowSize;
	double cc[ccLength];
	int rlstop = readLength+10;
	int flstart = 2*readLength;

	//if the fragment length is too short to compare against for artifacts
	if(2*readLength + 10 >= ccLength){
		//use fallback classifier
		rlstop = readLength;
		flstart = readLength+10;
	}

	//Calculate cross correlation nonbinarized nonsmoothed
	for(int ii = 0; ii < ccLength; ++ii){
		cc[ii]=0;
		for(int jj = 0; jj < coverageSize-ii; ++jj){
			cc[ii]+=candidate.pstrand[jj]*candidate.mstrand[jj+ii];
		}
	}
	//Score1 is the read length score
	double score1 = 0;
	for(long ii=0;ii<rlstop;++ii)
			if(cc[ii]>score1)
				score1 = cc[ii];
	//Score2 is the fragment length score
	double score2 = 0;
	for(long ii=flstart;ii<ccLength;++ii)
			if(cc[ii]>score2)
				score2 = cc[ii];
	//Add 1 to prevent divide by zero
	score2+=1;

	//get binarized smoothed cross correlation
	double smoothedBinPstrand[coverageSize];
	double smoothedBinMstrand[coverageSize];
	//run sum Smooth and keep only unique reads (After smoothing)
	for(long ii = 0; ii < coverageSize; ++ii){
		smoothedBinPstrand[ii]=0;
		smoothedBinMstrand[ii]=0;
		//smooth by 5 nucleotides
		for(long jj = max((long)0,ii-smoothingbandwidth); jj <= min(ii+smoothingbandwidth,coverageSize); ++jj){
			if(candidate.pstrand[jj]!=0)
				smoothedBinPstrand[ii]=1;
			if(candidate.mstrand[jj]!=0)
				smoothedBinMstrand[ii]=1;
		}
	}
	//Calculate cross correlation binarized smooth
	for(int ii = 0; ii < ccLength; ++ii){
		cc[ii]=0;
		for(int jj = 0; jj < coverageSize-ii; ++jj){
			cc[ii]+=smoothedBinPstrand[jj]*smoothedBinMstrand[jj+ii];
		}
	}

	//Score3 is the read length score
	double score3 = 0;
	for(long ii=0;ii<rlstop;++ii)
			if(cc[ii]>score3)
				score3 = cc[ii];
	//Score4 is the fragment length score
	double score4 = 0;
	for(long ii=flstart;ii<ccLength;++ii)
			if(cc[ii]>score4)
				score4 = cc[ii];
	//Add 1 to prevent divide by zero
	score4+=1;

	double feature1 = score1/score2;
	double feature2 = score3/score4;
	double activation = 17.582137+-2.678528*feature1+-13.505810*feature2;
	//if the fragment length is too short to compare against for artifacts
	if(2*readLength + 10 >= ccLength){
		//use the fallback classifier
		activation = 48.91067+-19.45451*feature1-29.32183*feature2;
	}
	//apply logit function
	double logitScore = 1/(1+exp(-activation));
	// return sigmoid activation
	return logitScore;
}

