/*
 * Peaks.cpp
 *
 *  Created on: May 5, 2015
 *      Author: Kelly
 */

#include "Peaks.h"
#include "BufferedDepthGraphReader.h"
#include <float.h>
#include <math.h>
#include "IOhandler.h"
#include "FFTHandler.h"

Peaks::Peaks() {
}

Peaks::~Peaks() {
	for(long ii =0; ii < (long)data.size();++ii){
		delete data[ii].pstrandExtended;
		delete data[ii].mstrandExtended;
	}
}

void Peaks::clear(){
	for(long ii =0; ii < (long)data.size();++ii){
		delete[] data[ii].pstrandExtended;
		delete[] data[ii].mstrandExtended;
	}
	data.clear();
}

double getScore(int windowSize,double* pstrand, double* mstrand, double* fir){
	//Calculate peak projection
	double pscore = 0;
	double mscore = 0;
	for(int ii =0; ii < 2*windowSize; ++ii){
		if(ii < windowSize)
			pscore += pstrand[ii]*fir[windowSize-1-ii];
		else
			mscore += mstrand[ii]*fir[ii-windowSize];
	}
	return pscore+mscore;
}
double getScoreDerivative(int windowSize,double* pstrand, double* mstrand, double* derivativeFIRP, double* derivativeFIRM){
	//Calculate peak projection
	double pscore = 0;
	double mscore = 0;
	for(int ii =0; ii < 2*windowSize; ++ii){
		pscore += pstrand[ii]*derivativeFIRP[ii];
		mscore += mstrand[ii]*derivativeFIRM[ii];
	}
	return pscore+mscore;
}
int getHalfLength(int windowSize, double* fir){
	//Calculate running max smoothing window as 50% of filter
	double firSum = 0;
	for(int ii =0; ii < windowSize; ii++)
		firSum+=fir[ii];
	double cumsum = 0;
	int halfLength;
	for(halfLength =0; halfLength < windowSize; halfLength++){
		cumsum+=fir[halfLength]/firSum;
		if(cumsum > 0.5) break;
	}
	return halfLength;
}
void Peaks::addPeak(long chr, long position, double* pstrand, double* mstrand){
	Peak currentResult;
	currentResult.chr = chr;
	//The center or 0 point is mfl from the start position of the window
	currentResult.pos = position+2*windowSize;
	currentResult.lpValue=0;
	currentResult.betaAlt=0;
	currentResult.isArtifact=false;
	currentResult.pstrandExtended=new unsigned short[4*windowSize];
	currentResult.mstrandExtended=new unsigned short[4*windowSize];
	currentResult.pstrand=&(currentResult.pstrandExtended[windowSize]);
	currentResult.mstrand=&(currentResult.mstrandExtended[windowSize]);
	for(int ii = 0;  ii < 4*windowSize; ++ii){
		currentResult.pstrandExtended[ii]=pstrand[ii];
		currentResult.mstrandExtended[ii]=mstrand[ii];
	}
	data.push_back(currentResult);

}

double GaussianDerivative(int ii, double sd, double center){
	return (ii-center)*exp(-pow((double)ii-center,2)/(2*pow(sd,2)));
}

double* createGuassianDerivativeFIRP(int windowSize, double* fir, double sd){
	FFTHandler::init(2*windowSize);
	double* gd = FFTHandler::newArray();
	double* gdFIR = FFTHandler::newArray();
	for(int ii = 0; ii < 2*windowSize; ++ii){
		if(ii < windowSize){
			gd[ii]=GaussianDerivative(ii, sd, 0);
			gdFIR[ii]=fir[windowSize-1-ii];
		}else{
			gd[ii]=GaussianDerivative(ii, sd, 2*windowSize);
			gdFIR[ii]=0;
		}
	}
	FFTHandler::inPlaceConvolve(gdFIR, gd);
	double* derivativeFIR = new double[2*windowSize];
	for(int ii = 0; ii < 2*windowSize; ++ii)
		derivativeFIR[ii]= gdFIR[ii];
	FFTHandler::freeArray(gd);
	FFTHandler::freeArray(gdFIR);
	FFTHandler::destroy();
	return derivativeFIR;
}
double* createGuassianDerivativeFIRM(int windowSize, double* fir, double sd){
	FFTHandler::init(2*windowSize);
	double* gd = FFTHandler::newArray();
	double* gdFIR = FFTHandler::newArray();
	for(int ii = 0; ii < 2*windowSize; ++ii){
		if(ii < windowSize){
			gd[ii]=GaussianDerivative(ii, sd, 0);
			gdFIR[ii]=0;
		}else{
			gd[ii]=GaussianDerivative(ii, sd, 2*windowSize);
			gdFIR[ii]=fir[ii-windowSize];
		}
	}
	FFTHandler::inPlaceConvolve(gdFIR, gd);
	double* derivativeFIR = new double[2*windowSize];
	for(int ii = 0; ii < 2*windowSize; ++ii)
		derivativeFIR[ii]= gdFIR[ii];
	FFTHandler::freeArray(gd);
	FFTHandler::freeArray(gdFIR);
	FFTHandler::destroy();
	return derivativeFIR;
}

void Peaks::find(char* bamFileName, int argWindowSize, double* fir, double peakThreshold, long minReadsThreshold, bool correctPCR){
	windowSize=argWindowSize;

	fprintf(stderr, "Looking for candidate peaks\n");

	BufferedDepthGraphReader* bgr;
	bgr = new BufferedDepthGraphReader(4*windowSize);
	bgr->init();

	//set up Gaussian derivative filter
	//double sd = ((double)getHalfLength(windowSize, fir))/4.0;
	double sd = 10;

	//fprintf(stderr, "Prepping filter derivatives\n");
	double* derivativeFIRP =createGuassianDerivativeFIRP(windowSize, fir, sd);
	double* derivativeFIRM =createGuassianDerivativeFIRM(windowSize, fir, sd);
	//IOhandler::printDoubleArrayToFile(derivativeFIRP,2*windowSize,"gdp.txt");
	//IOhandler::printDoubleArrayToFile(derivativeFIRM,2*windowSize,"gdm.txt");

	//current maximum values
	double previousFilterDerivative = -1;
	double currentFilterDerivative = -1;
	// for reporting
	double previousChromosome = -1;

	positionsScanned=0;
	while(bgr->next()){
/*
#ifdef DEBUG

		if(bgr->currentChromosome<0 || bgr->getPos()+2*windowSize<0)
			continue;
		if(bgr->currentChromosome>1 || bgr->getPos()+2*windowSize>1000000000)
			break;

#endif
*/
		//first apply read count threshold
		if(bgr->getPstrandReads()<minReadsThreshold/2.0 || bgr->getMstrandReads()<minReadsThreshold/2.0){
			previousFilterDerivative = -1;
			continue;
		}
		//If the window exceeds the minimum read count threshold increment the number of windows scanned
		++positionsScanned;
		//report progress
		if(bgr->currentChromosome!=previousChromosome && bgr->currentChromosome>=0 && bgr->currentChromosome < (int)bgr->chromosomeNames.size()){
			//fprintf(stderr, "Searching chromosome [%s]\n",bgr->chromosomeNames[bgr->currentChromosome].c_str());
		}
		previousChromosome = bgr->currentChromosome;
		//get the filter derivative
		currentFilterDerivative = getScoreDerivative(windowSize,&bgr->getPstrand()[windowSize],&bgr->getMstrand()[windowSize], derivativeFIRP, derivativeFIRM);
		//if we found a local maximum
		if(previousFilterDerivative > 0 && currentFilterDerivative <=0){
			//if the score is high then we are at a good maxima
			if(getScore(windowSize,&bgr->getPstrand()[windowSize],&bgr->getMstrand()[windowSize], fir)> ((double)2*minReadsThreshold)/(2*windowSize))
				addPeak(bgr->currentChromosome, bgr->getPos(), bgr->getPstrand(), bgr->getMstrand());
		}
		//store the derivative
		previousFilterDerivative = currentFilterDerivative;

		//next apply filter value threshold
		//if(getScore(windowSize,bgr.pstrand,bgr.mstrand, fir)< peakThreshold)
		//	continue;
	}
	bgr->close();
	delete bgr;
	delete[] derivativeFIRP;
	delete[] derivativeFIRM;
	findPeakNeighbors();
	fprintf(stderr, "[%ld] candidate peaks found from local maxima after scanning [%ld] positions\n",data.size(), positionsScanned);
}

void Peaks::findPeakNeighbors(){
	for(long ii =0; ii < (long)data.size(); ++ii){
		//add the current peak
		data[ii].localPeakPos.push_back(windowSize);
		data[ii].localPeakPtr.push_back(&(data[ii]));
		//add peaks before it within the window size
		for(long jj=ii-1; jj>=0 && data[jj].pos > data[ii].pos - windowSize && data[jj].chr==data[ii].chr; --jj){
			data[ii].localPeakPos.push_back(data[jj].pos-data[ii].pos+windowSize);
			data[ii].localPeakPtr.push_back(&(data[jj]));
		}
		//add peaks after it within the window size
		for(long jj=ii+1; jj<(long)data.size() && data[jj].pos < data[ii].pos + windowSize && data[jj].chr==data[ii].chr; ++jj){
			data[ii].localPeakPos.push_back(data[jj].pos-data[ii].pos+windowSize);
			data[ii].localPeakPtr.push_back(&(data[jj]));
		}
	}
}


//	//current maximum values
//	double currentMaxValue = -1;
//	long currentMaxPosition = -1;
//	long currentMaxChr = -1;
//	double currentMaxPstrand[2*windowSize];
//	double currentMaxMstrand[2*windowSize];
//
//	while(bgr.next()){
//#ifdef DEBUG
//		//FIXME debugging remove later
//		//if(bgr.currentChromosome!=0 ||bgr.getPos() +windowSize < 183250470 -1000 || bgr.getPos() +windowSize > 183250470 +1000)
//			//continue;
//		if(bgr.currentChromosome>0 || bgr.getPos()>97201154 +1000)
//		//if(bgr.currentChromosome!=8 || bgr.getPos()+windowSize <95640326-100 || bgr.getPos()+windowSize  > 95640326+100)
//		//	continue;
//#endif
//		if(bgr._pstrandReads+bgr._mstrandReads<minReadsThreshold)
//			continue;
//		//if we move out of the window analyze its local maxima
//		if(bgr.getPos()>= currentMaxPosition+runMaxWindow ||
//				bgr.currentChromosome!=currentMaxChr)
//		{
//			if(bgr.currentChromosome!=currentMaxChr && bgr.currentChromosome>=0 && bgr.currentChromosome < (int)bgr.chromosomeNames.size()){
//				fprintf(stderr, "Processing chromosome [%s]\n",bgr.chromosomeNames[bgr.currentChromosome].c_str());
//			}
//			//if the filtered signal exceeds the threshold record it
//			if(currentMaxValue > peakThreshold)
//			{
//				Peak currentResult;
//				currentResult.chr = currentMaxChr;
//				//The center or 0 point is mfl from the start position of the window
//				currentResult.pos = currentMaxPosition+windowSize;
//				currentResult.lpValue=0;
//				currentResult.betaAlt=0;
//				currentResult.isArtifact=false;
//				currentResult.pstrand=new double[2*windowSize];
//				currentResult.mstrand=new double[2*windowSize];
//				for(int ii = 0;  ii < 2*windowSize; ++ii){
//					currentResult.pstrand[ii]=currentMaxPstrand[ii];
//					currentResult.mstrand[ii]=currentMaxMstrand[ii];
//				}
//				data.push_back(currentResult);
//			}
//			//set local max to zero
//			currentMaxValue = -1;
//			currentMaxPosition=bgr.getPos();
//			currentMaxChr = bgr.currentChromosome;
//		}
//
//		//Calculate peak projection
//		double pscore = 0;
//		double mscore = 0;
//		for(int ii =0; ii < 2*windowSize; ++ii){
//			if(ii < windowSize)
//				pscore += bgr.pstrand[ii]*fir[windowSize-1-ii];
//			else
//				mscore += bgr.mstrand[ii]*fir[ii-windowSize];
//		}
//
//		//if it is higher record it
//		if(pscore+mscore > currentMaxValue &&
//				bgr.getPos()< currentMaxPosition+2*windowSize &&
//				bgr.currentChromosome==currentMaxChr){
//			currentMaxValue = pscore+mscore;
//			currentMaxPosition = bgr.getPos();
//			currentMaxChr = bgr.currentChromosome;
//			for(long ii = 0; ii < 2*windowSize;++ii){
//				currentMaxPstrand[ii]=bgr.pstrand[ii];
//				currentMaxMstrand[ii]=bgr.mstrand[ii];
//			}
//		}
//	}
//	bgr.close();
