/*
 * BufferedGenomeReader.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: kps3
 */

#include "BufferedGenomeReader.h"

BufferedGenomeReader::BufferedGenomeReader(long argWindowSize) {
	windowSize = argWindowSize;
	pstrandBuffer = 0;
	mstrandBuffer = 0;
	pstrand=0;
	mstrand=0;
	_pstrandReads = 0;
	_mstrandReads = 0;
}

BufferedGenomeReader::~BufferedGenomeReader() {
	fprintf(stderr, "cleaning up buffered reader\n");
	delete[] pstrandBuffer;
	delete[] mstrandBuffer;
}

void BufferedGenomeReader::close(){
	samstream.close();
}

void BufferedGenomeReader::init(const char* samFileName){
	samstream.init(samFileName);
	samstream.next();

	genomeLength = samstream.genomeLength;
	//Read the header
	chromosomeNames = samstream.chromosomeNames;
	readLength = samstream.readLength;
	currentChromosome= samstream.currentChromosome;
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

void BufferedGenomeReader::incrementBuffer(){
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
void BufferedGenomeReader::onChrChange(){
}
int BufferedGenomeReader::next(){

	//if the strands are both empty
	if(_pstrandReads ==0 && _mstrandReads==0){
		//if we have no further reads we are done
		if(!samstream.hasNext()) return 0;

		//otherwise skip to the next reads chromosome
		currentChromosome = samstream.currentChromosome;
		//and a windowSize before its position
		currentStartPosition = max((long)0,samstream.currentPosition-windowSize+1);
		onChrChange();
	}
	// else increment the buffer
	else{
		incrementBuffer();
	}

	//fill in any reads that should be in this window
	while(samstream.hasNext()
			&& samstream.currentChromosome==currentChromosome
			&& samstream.currentPosition < currentStartPosition+windowSize
			){
			pstrand[samstream.currentPosition-currentStartPosition]=samstream.pReads;
			_pstrandReads+=samstream.pReads;
			mstrand[samstream.currentPosition-currentStartPosition]=samstream.mReads;
			_mstrandReads+=samstream.mReads;
			samstream.next();
	}
	return 1;
}
long BufferedGenomeReader::getPos(){
	return currentStartPosition;
}
double* BufferedGenomeReader::getPstrand(){
	return pstrand;
}
double* BufferedGenomeReader::getMstrand(){
	return mstrand;
}

long BufferedGenomeReader::getPstrandReads(){
	return _pstrandReads;
}
long BufferedGenomeReader::getMstrandReads(){
	return _mstrandReads;
}
