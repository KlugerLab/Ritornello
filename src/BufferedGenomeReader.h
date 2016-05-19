/*
 * BufferedGenomeReader.h
 *
 *  Created on: Feb 11, 2015
 *      Author: kps3
 */

#ifndef BUFFEREDGENOMEREADER_H_
#define BUFFEREDGENOMEREADER_H_

#include <vector>
#include <string>
#include <queue>
#include <sam.h>
#include "SamStream.h"
using namespace std;
class BufferedGenomeReader {
public:
	BufferedGenomeReader(long argWindowSize);
	virtual ~BufferedGenomeReader();
	virtual void close();
	virtual void init(const char* samFileName);
	virtual int next();
	virtual long getPos();
	virtual double* getPstrand();
	virtual double* getMstrand();
	long readLength;
	long currentChromosome;
	long genomeLength;
	virtual long getPstrandReads();
	virtual long getMstrandReads();
	virtual void onChrChange();

	vector<string> chromosomeNames;
protected:
	virtual void incrementBuffer();
	long windowSize;
	long currentStartPosition;
	double* pstrandBuffer;
	double* mstrandBuffer;
	SamStream samstream;
	bool additionalReads;
	double* pstrand;
	double* mstrand;
	long _pstrandReads;
	long _mstrandReads;
};

#endif /* BUFFEREDGENOMEREADER_H_ */
