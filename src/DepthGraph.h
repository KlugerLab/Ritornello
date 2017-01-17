/*
 * SamToDepthGraph.h
 *
 *  Created on: Jan 17, 2017
 *      Author: Kelly P. Stanton
 */

#ifndef DEPTHGRAPH_H_
#define DEPTHGRAPH_H_

#include <string>
using namespace std;

class DepthGraph {
public:
	DepthGraph();
	virtual ~DepthGraph();
	void Sam2DepthGraph(const string& bamFileName, const string& outPrefix);
	void PCRcorrect(const string& outPrefix,int bandwidth);

};


#endif /* DEPTHGRAPH_H_ */
