/*
 * SamStream.h
 *
 *  Created on: Apr 30, 2015
 *      Author: kps3
 */

#ifndef SAMSTREAM_H_
#define SAMSTREAM_H_

#include <vector>
#include <string>
#include <sam.h>
using namespace std;
class SamStream {
public:
	SamStream();
	virtual ~SamStream();
	void close();
	void init(const char* samFileName);
	int next();
	bool hasNext();
	long readLength;
	vector<string> chromosomeNames;
	long genomeLength;
	long currentChromosome;
	long currentPosition;
	double pReads;
	double mReads;
private:
	void incrementBuffer();
	void lookaheadFill();
	bam1_t* alignment;
	int bytesRead;
	samfile_t* samFile;
	double* pBuffer;
	double* pBigBuffer;
	long pBufferReads;
	double* mBuffer;
	double* mBigBuffer;
	long mBufferReads;
};

#endif /* SAMSTREAM_H_ */
