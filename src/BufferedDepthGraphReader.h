/*
 * BufferedDepthGraphReader.h
 *
 *  Created on: Feb 11, 2015
 *      Author: kps3
 */

#ifndef BUFFEREDDEPTHGRAPHREADER_H_
#define BUFFEREDDEPTHGRAPHREADER_H_

#include <vector>
#include <string>
#include <queue>
#include <fstream>
using namespace std;
class BufferedDepthGraphReader {
public:
	BufferedDepthGraphReader(long argWindowSize);
	virtual ~BufferedDepthGraphReader();
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
	ifstream DepthGraphFile;
	bool additionalReads;
	double* pstrand;
	double* mstrand;
	long _pstrandReads;
	long _mstrandReads;
};

#endif /* BUFFEREDDEPTHGRAPHREADER_H_ */
