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

struct coverage{
	long startPos;
	long endPos;
	vector<double> pstrand;
	vector<double> mstrand;
	vector<Peak*> touched;
};
void getExtendedCoverage(coverage* extendedCoverage, Peak* candidate, long windowSize, long depth=0){
	const static long maxRecursionDepth =10;
	//mark that we have been here
	extendedCoverage->touched.push_back(candidate);
	if(depth>maxRecursionDepth)
		return;
	long candidateStartPos = candidate->pos - windowSize;
	long candidateEndPos = candidate->pos + windowSize;
	vector<double> candidatePstrand(candidate->pstrand, candidate->pstrand+ 2*windowSize);
	vector<double> candidateMstrand(candidate->mstrand, candidate->mstrand+ 2*windowSize);
	//for the first peak
	if(depth==0){
		//initialize extended coverage to the peaks start base pair
		extendedCoverage->startPos = candidateStartPos;
		extendedCoverage->endPos = candidateStartPos;
	}
	// if this peaks coverage starts before our extended coverage
	if (candidateStartPos < extendedCoverage->startPos){
		//get the length to insert
		long insertLength = extendedCoverage->startPos-candidateStartPos;
		//insert the pstrand
		extendedCoverage->pstrand.insert(extendedCoverage->pstrand.begin(),candidatePstrand.begin(),candidatePstrand.begin()+insertLength);
		//insert the mstrand
		extendedCoverage->mstrand.insert(extendedCoverage->mstrand.begin(),candidateMstrand.begin(),candidateMstrand.begin()+insertLength);
		//update the start position
		extendedCoverage->startPos = candidateStartPos;
	}
	// if this neighbors coverage ends after our extended coverage
	if(candidateEndPos > extendedCoverage->endPos){
		//get the length to insert
		long insertLength = candidateEndPos - extendedCoverage->endPos;
		//insert the pstrand
		extendedCoverage->pstrand.insert(extendedCoverage->pstrand.end(),candidatePstrand.end()-insertLength,candidatePstrand.end());
		//insert the mstrand
		extendedCoverage->mstrand.insert(extendedCoverage->mstrand.end(),candidateMstrand.end()-insertLength,candidateMstrand.end());
		//update the end position
		extendedCoverage->endPos = candidateEndPos;
	}
	//recurse over neighbor peaks that we havent touched yet
	for(long ii = 0; ii < (long)(candidate->localPeakPtr.size()); ++ii){
		//check if we've been here before to prevent cycling
		if(std::find(
				extendedCoverage->touched.begin(),
				extendedCoverage->touched.end(),
				candidate->localPeakPtr[ii]) != extendedCoverage->touched.end())
			continue;
		//get coverage for this node
		getExtendedCoverage(extendedCoverage,candidate->localPeakPtr[ii], windowSize, depth+1);
	}
}
double Artifact::test(long windowSize, long readLength, double* fir, Peak candidate, double testRatio, long halfLength){
	coverage extendedCoverage;
	getExtendedCoverage(&extendedCoverage,&candidate, windowSize);
	IOhandler::printDoubleArrayToFile(&extendedCoverage.pstrand[0],extendedCoverage.pstrand.size(),"extendedCoveragePstrand.txt");
	IOhandler::printDoubleArrayToFile(&extendedCoverage.mstrand[0],extendedCoverage.pstrand.size(),"extendedCoveragemstrand.txt");

	long coverageSize = extendedCoverage.endPos - extendedCoverage.startPos;
	double smoothedBinPstrand[coverageSize];
	double smoothedBinMstrand[coverageSize];
	double smoothedPstrand[coverageSize];
	double smoothedMstrand[coverageSize];
	long smoothingbandwidth = 2;

	//run sum Smooth
	for(long ii = 0; ii < coverageSize; ++ii){
		smoothedPstrand[ii]=0;
		smoothedMstrand[ii]=0;
		//smooth by 5 nucleotides
		for(long jj = max((long)0,ii-smoothingbandwidth); jj <= min(ii+smoothingbandwidth,coverageSize); ++jj){
			smoothedPstrand[ii]+=extendedCoverage.pstrand[jj];
			smoothedMstrand[ii]+=extendedCoverage.mstrand[jj];
		}
	}

	//run sum Smooth and keep only unique reads (After smoothing)
	for(long ii = 0; ii < coverageSize; ++ii){
		smoothedBinPstrand[ii]=0;
		smoothedBinMstrand[ii]=0;
		//smooth by 5 nucleotides
		for(long jj = max((long)0,ii-smoothingbandwidth); jj <= min(ii+smoothingbandwidth,coverageSize); ++jj){
			if(extendedCoverage.pstrand[jj]!=0)
				smoothedBinPstrand[ii]=1;
			if(extendedCoverage.mstrand[jj]!=0)
				smoothedBinMstrand[ii]=1;
		}
	}

	//Calculate cross correlation binarized smooth
	int ccLength = 0.75*windowSize;
	double cc[ccLength];
	for(int ii = 0; ii < ccLength; ++ii){
		cc[ii]=0;
		for(int jj = 0; jj < coverageSize-ii; ++jj){
			cc[ii]+=smoothedBinPstrand[jj]*smoothedBinMstrand[jj+ii];
		}
	}
	IOhandler::printDoubleArrayToFile(cc,ccLength,"ccBinSmooth.txt");


	double artifactScore = 0;
	//long peakCutoff = max(readLength+10, halfLength);
	long peakCutoff = readLength+10;
	double peakScore = 0;
	for(int ii = 0; ii < ccLength; ++ii){
		if(ii < readLength && cc[ii] > artifactScore)
			artifactScore = cc[ii];
		else if(ii >= peakCutoff && cc[ii] > peakScore)
			peakScore = cc[ii];
	}
	//if it passes the binarized smooth cc then do the non binarized not smoothed version

	//Calculate cross correlation nonbinarized nonsmoothed
	for(int ii = 0; ii < ccLength; ++ii){
		cc[ii]=0;
		for(int jj = 0; jj < coverageSize-ii; ++jj){
			cc[ii]+=extendedCoverage.pstrand[jj]*extendedCoverage.mstrand[jj+ii];
		}
	}
	IOhandler::printDoubleArrayToFile(cc,ccLength,"ccPlain.txt");
	//nonbin nonsmoothed
	double artifactScoreNonBin = 0;
	double peakScoreNonBin = 0;
	for(int ii = 0; ii < ccLength; ++ii){
		if(ii < readLength && cc[ii] > artifactScoreNonBin)
			artifactScoreNonBin = cc[ii];
		else if(ii >= peakCutoff && cc[ii] > peakScoreNonBin)
			peakScoreNonBin = cc[ii];
	}

	//double logitScore = -0.3208*artifactScore+0.2985*peakScore+-0.1135*artifactScoreNonBin+0.1200*peakScoreNonBin;
	double logitScore = 17.871*peakScore/artifactScore+6.035*peakScoreNonBin/artifactScoreNonBin + -25.179;
	logitScore = 1/(1+exp(-logitScore));
	if(logitScore < 0.5){
		return logitScore;
	}

	double multiArtifactScore = multiArtifactTest(windowSize, readLength, fir, candidate, testRatio, halfLength);
	return logitScore*multiArtifactScore;
}

//returns 0 if this is an artifact, 1 otherwise
double Artifact::multiArtifactTest(long windowSize, long readLength, double* fir, Peak candidate, double testRatio, long halfLength){
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

	//set up filters
	double* psFilter = FFTHandler::newArray();
	double* msFilter = FFTHandler::newArray();
	double* nullFilter = FFTHandler::newArray();
	for(int ii = 0; ii < 2*FOURIERLENGTH; ++ii){
		psFilter[ii]=0;
		psFilter[ii]=0;
		nullFilter[ii]=0;
		if(ii >= FOURIERLENGTH && ii < (FOURIERLENGTH+readLength))
			psFilter[ii]=1/sqrt(2*readLength);
		if(ii >= (FOURIERLENGTH-readLength) && ii < FOURIERLENGTH)
			msFilter[ii]=1/sqrt(2*readLength);
		if(ii >= (FOURIERLENGTH-2*readLength) && ii < (FOURIERLENGTH+2*readLength))
			nullFilter[ii]=1/sqrt(4*readLength);
	}

	//calculate the cross correlation
	FFTHandler::inPlaceConvolve(psFilter, ps);
	FFTHandler::inPlaceConvolve(msFilter, ms);
	FFTHandler::inPlaceConvolve(ps, nullFilter);
	FFTHandler::inPlaceConvolve(ms, nullFilter);

	double ratio[2*FOURIERLENGTH];
	for(long ii=0;ii<2*FOURIERLENGTH;++ii)
		ratio[ii]= (psFilter[ii]+msFilter[ii]+1)/(ps[ii]+ms[ii]+1);

	FFTHandler::freeArray(ps);
	FFTHandler::freeArray(ms);
	FFTHandler::freeArray(psFilter);
	FFTHandler::freeArray(msFilter);
	FFTHandler::freeArray(nullFilter);

	double feature = 0;
	for(long ii=0;ii<2*FOURIERLENGTH;++ii)
			if(ratio[ii]>feature)
				feature = ratio[ii];

	double activation = 2.324354-27.916630*log(feature);

	//apply logit function
	double logitScore = 1/(1+exp(-activation));
	// return sigmoid activation
	return logitScore;
}
