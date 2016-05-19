/*
 * Filter.cpp
 *
 *  Created on: Jun 15, 2014
 *      Author: Kelly Patrick Stanton
 */

#include "FLD.h"
#include "BufferedGenomeReader.h"
#include "FFTHandler.h"
#include <math.h>
#include "IOhandler.h"
#include <boost/math/distributions/chi_squared.hpp>
#include "PCRCorrectGenomeReader.h"

namespace FLD{
double* calculate(char* bamFileName, int windowSize, bool correctPCR){
	//test to see if it is paired end data
	bool isPairedEnd=false;
	samfile_t* samFile = samopen(bamFileName, "rb", 0);
	bam1_t* alignment;
	alignment = new bam1_t();
	int bytesRead = samread(samFile,alignment);
	if(alignment->core.tid==-1 || bytesRead <=0)
	{
		fprintf(stderr, "Failed to read any mapping alignment%s\n", bamFileName);
		exit(EXIT_FAILURE);
	}
	if(alignment->core.isize > 0)
		isPairedEnd=true;
	delete alignment;
	samclose(samFile);
	//Get clean FLD
	if(isPairedEnd)
		return calculatePairedEnd(bamFileName, windowSize);
	else
		return calculateSingleEnd(bamFileName, windowSize, correctPCR);
}
double* calculatePairedEnd(char* bamFileName, int windowSize){
	double* fld = new double[windowSize];
	for(int ii = 0; ii <windowSize; ++ii )
		fld[ii]=0;


	//open sam file for reading
	samfile_t* samFile = samopen(bamFileName, "rb", 0);
	bam1_t* alignment;
	alignment = new bam1_t();

	int bytesRead = samread(samFile,alignment);
	//tally up the fragment lengths
	while(bytesRead>0){
		int fragmentLength = alignment->core.isize;
		if( fragmentLength < windowSize && fragmentLength > 0)
			++fld[fragmentLength];
		bytesRead = samread(samFile,alignment);
	}
	delete alignment;
	samclose(samFile);

	//normalize to sum to 1
	double sum = 0;
	for(int ii = 0; ii <windowSize; ++ii )
		sum+=fld[ii];
	for(int ii = 0; ii <windowSize; ++ii)
		fld[ii]/=sum;
	return fld;
}
double* calculateSingleEnd(char* bamFileName, int windowSize, bool correctPCR){

	boost::math::chi_squared_distribution<double> chi2(4*windowSize-1);
	double threshold=0.5;

	FFTHandler::init(4*windowSize);
	double* ac = FFTHandler::newArray();
	double* cc = FFTHandler::newArray();
	for(int ii = 0; ii < 4*windowSize; ++ii){
		cc[ii]=0;
		ac[ii]=0;
	}
	BufferedGenomeReader* bgr;
	if(correctPCR)
		bgr = new PCRCorrectGenomeReader(4*windowSize);
	else
		bgr = new BufferedGenomeReader(4*windowSize);
	bgr->init(bamFileName);
	int iterations = 0;
#ifdef DEBUG
	while(bgr->next()){
	//while(bgr.next() && iterations < 100000000){
#endif
#ifndef DEBUG
	while(bgr->next()){
#endif
		if(bgr->getPstrandReads()==0 || bgr->getMstrandReads()==0)
			continue;

		if(bgr->getPstrand()[2*windowSize]!=0){
			//calcualate mean coverage
			double pmean = bgr->getPstrandReads()/(double)(4*windowSize);
			//Test for uniform using Chi2 goodness of fit
			double X2=0;
			for(int ii =0; ii < 4*windowSize; ++ii){
				X2+=pow(bgr->getPstrand()[ii]-pmean,2)/pmean;
			}
			double chiX2 = 1-cdf(chi2,X2);
			if(chiX2>threshold){
				//add the autocorrelation and cross-correlation assuming pstrand is independent
				for(int ii = 0; ii < 4*windowSize; ++ii){
					cc[ii]+=bgr->getMstrand()[ii]*bgr->getPstrand()[2*windowSize];
					ac[ii]+=bgr->getPstrand()[ii]*bgr->getPstrand()[2*windowSize];
				}
			}

		}
		if(bgr->getMstrand()[2*windowSize]!=0){
			//calcualate mean coverage
			double mmean = bgr->getMstrandReads()/(double)(4*windowSize);
			//Test for uniform using Chi2 goodness of fit
			double X2=0;
			for(int ii =0; ii < 4*windowSize; ++ii){
				X2+=pow(bgr->getMstrand()[ii]-mmean,2)/mmean;
			}
			double chiX2 = 1-cdf(chi2,X2);
			if(chiX2>threshold){
				//add the autocorrelation and cross-correlation assuming mstrand is independent
				for(int ii = 0; ii < 4*windowSize; ++ii){
					cc[ii]+=bgr->getPstrand()[4*windowSize-ii-1]*bgr->getMstrand()[2*windowSize];
					ac[ii]+=bgr->getMstrand()[ii]*bgr->getMstrand()[2*windowSize];
				}
			}

		}
		++iterations;
	}

	IOhandler::printDoubleArrayToFile(cc,4*windowSize,"cc.txt");
	IOhandler::printDoubleArrayToFile(ac,4*windowSize,"ac.txt");

	//deconvolve autocorrelation from cross correlation
	FFTHandler::inPlaceDeconvolve(cc,ac);
	//copy result to fld vector
	double* fld = new double[windowSize];
	for(int ii = 0; ii <windowSize; ++ii )
		fld[ii]=cc[ii];
	//free FFT memory
	FFTHandler::freeArray(cc);
	FFTHandler::freeArray(ac);
	FFTHandler::destroy();

	bgr->close();
	delete bgr;

//	//find the max of areas that are supposed to be zero
//	double maxOffset = 0;
//	for(int ii =0; ii <windowSize; ++ii )
//			if(fld[ii]>maxOffset && (ii < bgr.readLength || ii >= windowSize-50))
//				maxOffset = fld[ii];
//
//	//Subtract y offset and set negative areas to 0
//	for(int ii =0; ii <windowSize; ++ii )
//		fld[ii]=max(fld[ii]-maxOffset,0.0);
//
//	//set area below the read length to 0
//	for(int ii =0; ii < bgr.readLength; ++ii )
//		fld[ii]=0;

	//subtract off y offset
	double meanNoise = 0;
	for(int ii =4*windowSize/6; ii <windowSize; ++ii )
		meanNoise+=fld[ii];
	meanNoise/=windowSize-4*windowSize/6;
	for(int ii = 0; ii <windowSize; ++ii)
		fld[ii]-=meanNoise;

	//set area below the read length to 0
	for(int ii =0; ii < bgr->readLength; ++ii )
		fld[ii]=0;

	//interpolate to 2x fragment length
	for(int ii =bgr->readLength; ii < 2*bgr->readLength; ++ii )
		fld[ii]=(ii-bgr->readLength)*fld[2*bgr->readLength]/bgr->readLength;
	//normalize to sum to 1
	double sum = 0;
	for(int ii = 0; ii <windowSize; ++ii )
		sum+=fld[ii];
	for(int ii = 0; ii <windowSize; ++ii)
		fld[ii]/=sum;

	//return the FLD
	return fld;
}
}
