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
}

BufferedDepthGraphReader::~BufferedDepthGraphReader() {
	fprintf(stderr, "cleaning up buffered reader\n");
	delete[] pstrandBuffer;
	delete[] mstrandBuffer;
}

void BufferedDepthGraphReader::close(){
	DepthGraphFile.close();
}

void BufferedDepthGraphReader::init(const string& outPrefix){

	//read key file info
	ifstream DepthGraphKeyFile.open((outPrefix+".RitorDepthGraphKey").c_str(),ios::in);
	readLength<<DepthGraphKeyFile;
	genomeLength<<DepthGraphKeyFile;
	//Read the header
	chromosomeNames.clear();
	while(!DepthGraphKeyFile.eof()){
		string line<<DepthGraphKeyFile;
		chromosomeNames.push_back(line);
	}

	//Read coverage from depth graph
	DepthGraphFile.open((outPrefix+".RitorDepthGraph").c_str(),ios::in);
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
