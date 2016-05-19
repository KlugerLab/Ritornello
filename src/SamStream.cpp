/*
 * SamStream.cpp
 *
 *  Created on: Apr 30, 2015
 *      Author: kps3
 */

#include "SamStream.h"

SamStream::SamStream() {

}

SamStream::~SamStream() {

}

void SamStream::init(const char* samFileName) {
	//Open sam file
	samFile = samopen(samFileName, "rb", 0);
	if (samFile == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", samFileName);
		exit(EXIT_FAILURE);
	}

	genomeLength=0;
	//Read the header
	chromosomeNames.resize(samFile->header->n_targets);
	for(int ii=0;ii < samFile->header->n_targets; ++ii){
		chromosomeNames[ii] = samFile->header->target_name[ii];
		genomeLength +=samFile->header->target_len[ii];
	}


	//check that we can read from it.  Get read length
	alignment = new bam1_t();
	bytesRead = samread(samFile,alignment);
	if(alignment->core.tid==-1 || bytesRead <=0)
	{
		fprintf(stderr, "Failed to read any mapping alignment%s\n", samFileName);
		exit(EXIT_FAILURE);
	}
	readLength = alignment->core.l_qseq;
	currentChromosome=alignment->core.tid;
	currentPosition = -1;

	pBigBuffer = new double[2*readLength];
	mBigBuffer = new double[2*readLength];
	for(long ii = 0; ii < readLength; ++ii){
		pBigBuffer[ii]=0;
		mBigBuffer[ii]=0;
	}
	pBuffer = pBigBuffer;
	mBuffer = mBigBuffer;
	pBufferReads=0;
	mBufferReads=0;
}

void SamStream::close(){
	delete pBigBuffer;
	delete mBigBuffer;
	delete alignment;
	fprintf(stderr, "closing sam file\n");
	samclose(samFile);
}

void SamStream::lookaheadFill(){
	//read all new alignments a read length ahead
	bool forwardStrand = (alignment->core.flag&BAM_FREVERSE) == 0;
	long position = forwardStrand ? alignment->core.pos : alignment->core.pos+readLength-1;
	//look ahead the read length
	while(alignment->core.tid==currentChromosome && position <= currentPosition+readLength-1 &&
			bytesRead>0 && alignment->core.tid!=-1){
		if(forwardStrand){
			++pBuffer[position-currentPosition];
			++pBufferReads;
		}
		else{
			++mBuffer[position-currentPosition];
			++mBufferReads;
		}
		bytesRead = samread(samFile,alignment);
		forwardStrand = (alignment->core.flag&BAM_FREVERSE) == 0;
		position = forwardStrand ? alignment->core.pos : alignment->core.pos+readLength-1;
	}

}
void SamStream::incrementBuffer(){
	//subtract first element reads
	pBufferReads-=pBuffer[0];
	mBufferReads-=mBuffer[0];

	//reset pBuffer and mBuffer back to the start of the buffers if they exceed the size
	if(pBuffer==&pBigBuffer[readLength]){
		pBuffer = pBigBuffer;
		mBuffer = mBigBuffer;
		for(long ii = 0; ii < readLength -1; ++ii){
			pBuffer[ii]=pBigBuffer[readLength+ii+1];
			mBuffer[ii]=mBigBuffer[readLength+ii+1];
		}
	}
	//else increment their pointers to the next position of the big buffers
	else{
		++pBuffer;
		++mBuffer;
	}

	//set the last element to 0 since it will be new
	pBuffer[readLength-1]=0;
	mBuffer[readLength-1]=0;

	//increment the position
	++currentPosition;
}

int SamStream::next(){
	//increment them
	incrementBuffer();
	//and lookahead fill them
	lookaheadFill();
	while(pBuffer[0]==0 && mBuffer[0]==0){
		//if the buffers are both empty
		if(pBufferReads ==0 && mBufferReads==0){
			//if we have no further reads we are done
			if(bytesRead<=0 || alignment->core.tid==-1) return 0;

			//otherwise skip to the next reads chromosome
			currentChromosome = alignment->core.tid;
			//and a readlength before its position
			bool forwardStrand = (alignment->core.flag&BAM_FREVERSE) == 0;
			long position = forwardStrand ? alignment->core.pos : alignment->core.pos+readLength-1;
			currentPosition = position - readLength+1;
			//and Fill the buffer with the new reads
			lookaheadFill();
		}
		//else increment
		else{
			//increment them
			incrementBuffer();
			//and lookahead fill them
			lookaheadFill();
		}
	}

	//Set pReads and mReads from the lookahead buffers
	pReads = pBuffer[0];
	mReads = mBuffer[0];
	return 1;
}

bool SamStream::hasNext(){
	//if the buffers are empty and there are no more reads from the sam file
	if(pBufferReads ==0 && mBufferReads==0 &&
			(bytesRead<=0 || alignment->core.tid==-1)){
		return false;
	}
	return true;
}
