/*
 * BufferedDepthGraphReader.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: kps3
 */

#include "BufferedDepthGraphReader.h"

BufferedDepthGraphReader::BufferedDepthGraphReader(long argWindowSize) {
	windowSize = argWindowSize;
	pstrandBuffer = 0;
	mstrandBuffer = 0;
	pstrand=0;
	mstrand=0;
	_pstrandReads = 0;
	_mstrandReads = 0;
	readLength = 0;
}

BufferedDepthGraphReader::~BufferedDepthGraphReader() {
	fprintf(stderr, "cleaning up buffered reader\n");
	delete[] pstrandBuffer;
	delete[] mstrandBuffer;
}

void BufferedDepthGraphReader::close(){
	dg.close();
}

void BufferedDepthGraphReader::init(const string& outPrefix){

	//Read coverage from depth graph
	dg.open(outPrefix);
	currentChromosome= dg.currentChromosome;
	currentStartPosition = -1;

	pstrandBuffer = new double[2*windowSize];
	mstrandBuffer = new double[2*windowSize];
	pstrand=pstrandBuffer;
	mstrand=mstrandBuffer;
	for(long ii = 0; ii < windowSize; ++ii){
		pstrand[ii]=0;
		mstrand[ii]=0;
	}
}

void BufferedDepthGraphReader::incrementBuffer(){
	//subtract off first element from read counts
	_pstrandReads-=pstrand[0];
	_mstrandReads-=mstrand[0];

	//reset pstrand and mstrands back to the start of the buffers if they exceed the size
	if(pstrand==&pstrandBuffer[windowSize]){
		pstrand = pstrandBuffer;
		mstrand = mstrandBuffer;
		for(long ii = 0; ii < windowSize -1; ++ii){
			pstrand[ii]=pstrandBuffer[windowSize+ii+1];
			mstrand[ii]=mstrandBuffer[windowSize+ii+1];
		}
	}
	//else increment their pointers to the next position of the buffers
	else{
		++pstrand;
		++mstrand;
	}

	//set the new position to 0
	pstrand[windowSize-1]=0;
	mstrand[windowSize-1]=0;

	//Shift the window
	++currentStartPosition;
}
void BufferedDepthGraphReader::onChrChange(){
}
int BufferedDepthGraphReader::next(){

	//if the strands are both empty
	if(_pstrandReads ==0 && _mstrandReads==0){
		//if we have no further reads we are done
		if(!dg.hasNext()) return 0;

		//otherwise skip to the next reads chromosome
		currentChromosome = dg.currentChromosome;
		//and a windowSize before its position
		currentStartPosition = max((long)0,dg.currentPosition-windowSize+1);
		onChrChange();
	}
	// else increment the buffer
	else{
		incrementBuffer();
	}

	//fill in any reads that should be in this window
	while(dg.hasNext()
			&& dg.currentChromosome==currentChromosome
			&& dg.currentPosition < currentStartPosition+windowSize
			){
			pstrand[dg.currentPosition-currentStartPosition]=dg.pReads;
			_pstrandReads+=dg.pReads;
			mstrand[dg.currentPosition-currentStartPosition]=dg.mReads;
			_mstrandReads+=dg.mReads;
			dg.next();
	}
	return 1;
}
long BufferedDepthGraphReader::getPos(){
	return currentStartPosition;
}
double* BufferedDepthGraphReader::getPstrand(){
	return pstrand;
}
double* BufferedDepthGraphReader::getMstrand(){
	return mstrand;
}

long BufferedDepthGraphReader::getPstrandReads(){
	return _pstrandReads;
}
long BufferedDepthGraphReader::getMstrandReads(){
	return _mstrandReads;
}
