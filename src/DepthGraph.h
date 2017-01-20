/*
 * SamToDepthGraph.h
 *
 *  Created on: Jan 17, 2017
 *      Author: Kelly P. Stanton
 */

#ifndef DEPTHGRAPH_H_
#define DEPTHGRAPH_H_

#include <string>
#include <fstream>
#include <vector>
using namespace std;

//position with coverage
struct position{
	unsigned char chr;
	unsigned short preads;
	unsigned short mreads;
	long pos;
};

class DepthGraph {
public:
	DepthGraph();
	virtual ~DepthGraph();
	void Sam2DepthGraph(const string& bamFileName, const string& outPrefix);
	void PCRcorrect(const string& outPrefix,int bandwidth);
	void open(const string& outPrefix);
	void close();
	int next();
	bool hasNext();

	long readLength;
	vector<string> chromosomeNames;
	long genomeLength;
	long currentChromosome;
	long currentPosition;
	short pReads;
	short mReads;

private:
	ifstream depthFile;
	position positionBuffer;
};


#endif /* DEPTHGRAPH_H_ */
