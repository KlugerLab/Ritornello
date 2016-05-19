/*
 * PCRCorrectGenomeReader.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: Kelly P. Stanton
 */

#include "PCRCorrectGenomeReader.h"
#include <math.h>

long PCRCorrectGenomeReader::bandwidth=30;
PCRCorrectGenomeReader::PCRCorrectGenomeReader(long argWindowSize)
:BufferedGenomeReader(argWindowSize+2*bandwidth){
	_pstrandReadsCorrected = 0;
	_mstrandReadsCorrected = 0;
}
PCRCorrectGenomeReader::~PCRCorrectGenomeReader() {
}
void PCRCorrectGenomeReader::init(const char* samFileName){
	BufferedGenomeReader::init(samFileName);
	correctedPstrandBuffer = new double[2*windowSize-2*bandwidth];
	correctedMstrandBuffer = new double[2*windowSize-2*bandwidth];
	correctedPstrand=correctedPstrandBuffer;
	correctedMstrand=correctedMstrandBuffer;
	for(long ii = 0; ii < windowSize-bandwidth; ++ii){
		correctedPstrand[ii]=-1;
		correctedMstrand[ii]=-1;
	}
}
void PCRCorrectGenomeReader::onChrChange(){
	_pstrandReadsCorrected = 0;
	_mstrandReadsCorrected = 0;
	for(long ii = 0; ii < windowSize-bandwidth; ++ii){
		correctedPstrand[ii]=-1;
		correctedMstrand[ii]=-1;
	}
}
void PCRCorrectGenomeReader::incrementBuffer(){
	BufferedGenomeReader::incrementBuffer();
	//subtract off first element from read counts
	_pstrandReadsCorrected-=correctedPstrand[0];
	_mstrandReadsCorrected-=correctedMstrand[0];

	//reset pstrand and mstrands back to the start of the buffers if they exceed the size
	if(correctedPstrand==&correctedPstrandBuffer[windowSize-bandwidth]){
		correctedPstrand = correctedPstrandBuffer;
		correctedMstrand = correctedMstrandBuffer;
		for(long ii = 0; ii < windowSize -1; ++ii){
			correctedPstrand[ii]=correctedPstrandBuffer[windowSize-bandwidth+ii+1];
			correctedMstrand[ii]=correctedMstrandBuffer[windowSize-bandwidth+ii+1];
		}
	}
	//else increment their pointers to the next position of the buffers
	else{
		++correctedPstrand;
		++correctedMstrand;
	}

	//set the new position to -1 (untouched
	correctedPstrand[windowSize-bandwidth-1]=-1;
	correctedMstrand[windowSize-bandwidth-1]=-1;
}
long PCRCorrectGenomeReader::getPos(){
	return BufferedGenomeReader::getPos()+bandwidth/2;
}
double* PCRCorrectGenomeReader::getPstrand(){
	return correctedPstrand;
}
double* PCRCorrectGenomeReader::getMstrand(){
	return correctedMstrand;
}

long PCRCorrectGenomeReader::getPstrandReads(){
	return _pstrandReadsCorrected;
}
long PCRCorrectGenomeReader::getMstrandReads(){
	return _mstrandReadsCorrected;
}

int PCRCorrectGenomeReader::next(){
	int ret = BufferedGenomeReader::next();
	//fill in any reads that should be in this window
	for(int ii =0; ii < windowSize-bandwidth; ++ii){
		if(correctedPstrand[ii]==-1){
			correctedPstrand[ii]=correctedPReadsAt(ii+bandwidth/2);
			_pstrandReadsCorrected+=correctedPstrand[ii];
			correctedMstrand[ii]=correctedMReadsAt(ii+bandwidth/2);
			_mstrandReadsCorrected+=correctedMstrand[ii];
		}
	}
	return ret;
}

//log Poisson distribution pdf
inline double lPois(long k, double mu){
	return k*log(mu) -mu -lgamma(k+1);
}
long PCRCorrectGenomeReader::correctedPReadsAt(long pos){
	//get local mean
	double mean = 0;
	for(int ii = pos-bandwidth/2;ii < pos+bandwidth/2;++ii){
		mean += pstrand[ii];
	}
	mean/=bandwidth;
	if(lPois(pstrand[pos],mean)<log(0.01)){
		return (int)ceil(mean);
	}
	else{
		return pstrand[pos];
	}
}
long PCRCorrectGenomeReader::correctedMReadsAt(long pos){
	//get local mean
	double mean = 0;
	for(int ii = pos-bandwidth/2;ii < pos+bandwidth/2;++ii){
		mean += mstrand[ii];
	}
	mean/=bandwidth;
	if(lPois(mstrand[pos],mean)<log(0.01)){
		return (int)ceil(mean);
	}
	else{
		return mstrand[pos];
	}
}
