/*
 * SamToDepthGraph.cpp
 *
 *  Created on: Jan 17, 2017
 *      Author: Kelly P. Stanton
 */

#include "DepthGraph.h"
#include "BufferedGenomeReader.h"
#include "SamStream.h"
#include <fstream>
#include <deque>
#include <math.h>

DepthGraph::DepthGraph() {
	// TODO Auto-generated constructor stub

}

DepthGraph::~DepthGraph() {
	// TODO Auto-generated destructor stub
}

//position with coverage
struct position{
	unsigned char chr;
	unsigned short preads;
	unsigned short mreads;
	long pos;
};
void DepthGraph::Sam2DepthGraph(const string& bamFileName, const string& outPrefix){
	fprintf(stderr, "Converting sam file to depth graph\n");
	//open the sam file
	SamStream samstream;
	samstream.init(bamFileName.c_str());
	samstream.next();

	//write Key file
	ofstream keyFile((outPrefix+".RitorDepthGraphKey").c_str(),ios::out);
	keyFile<<samstream.readLength<<endl;
	for(unsigned int ii = 0; ii < samstream.chromosomeNames.size();++ii){
		keyFile<<samstream.chromosomeNames[ii]<<endl;samstream.chromosomeNames[ii];
	}
	keyFile.close();

	//open the depth file to write to
    ofstream depthFile((outPrefix+".RitorDepthGraph").c_str(), ios::out | ios::binary);

    //output reads on each strand at each position to the depth file
    position buffer;
	while(samstream.hasNext()){
		buffer.chr=samstream.currentChromosome;
		buffer.pos=samstream.currentPosition;
		buffer.preads=samstream.pReads;
		buffer.mreads=samstream.mReads;
		samstream.next();
		depthFile.write((char*)(&buffer), sizeof(buffer));
	}
	//close both files
	samstream.close();
	depthFile.close();
	fprintf(stderr, "Finished converting sam file to depth graph\n");
}

//log Poisson distribution pdf
inline double lPois(long k, double mu){
	return k*log(mu) -mu -lgamma(k+1);
}
//PCR correct the depth file and output a new depth file
void DepthGraph::PCRcorrect(const string& outPrefix, int bandwidth){
	fprintf(stderr, "PCR correcting depth graph\n");
	deque<position> window;
	unsigned int ind=0;
	position positionBuffer;
	position currentPositionPCRcorrected;
	short windowPreads=0;
	short windowMreads=0;

	//open file
	ifstream depthFile ((outPrefix+".RitorDepthGraph").c_str(), ios::in | ios::binary);
	if (!depthFile) {
		// An error occurred!
		throw "Could not open depth file";
	}

	//open the PCR corrected depth file to write to
    ofstream PCRcorrectDepthFile ((outPrefix+"-PCRCorrect.RitorDepthGraph").c_str(), ios::out | ios::binary);

	//seed a new position to the buffer
	depthFile.read((char*)(&positionBuffer), sizeof(positionBuffer));

	//while we still have positions
	while(!depthFile.eof() || ind!=window.size()) {

		//if we are at end of the queue get another entry
		if(ind == window.size()){
			//push a new position onto the back
			window.push_back(positionBuffer);
			//add its reads
			windowPreads+=positionBuffer.preads;
			windowMreads+=positionBuffer.mreads;
			depthFile.read((char*)(&positionBuffer), sizeof(positionBuffer));
		}

		//Fill up the positions after the current
		while(!depthFile.eof()
				&& positionBuffer.pos < window.at(ind).pos+bandwidth/2+1
				&& positionBuffer.chr ==window.at(ind).chr){
			//add it to the positions after queue
			window.push_back(positionBuffer);
			//add its reads
			windowPreads+=positionBuffer.preads;
			windowMreads+=positionBuffer.mreads;
			depthFile.read((char*)(&positionBuffer), sizeof(positionBuffer));
		}
		//remove the previous positions now too far
		while(window.front().pos < window.at(ind).pos-bandwidth/2
				|| window.front().chr != window.at(ind).chr){
			windowPreads-=window.front().preads;
			windowMreads-=window.front().mreads;
			window.pop_front();
			--ind;
		}

		//PCR correct the current position
		currentPositionPCRcorrected = window.at(ind);
		//if there is more than 1 read on the positive strand PCR correct the current position
		if(currentPositionPCRcorrected.preads >1){
			double mean = ((double)(windowPreads-currentPositionPCRcorrected.preads))/bandwidth;
			if(currentPositionPCRcorrected.preads > mean && lPois(currentPositionPCRcorrected.preads,mean)<log(0.01))
				currentPositionPCRcorrected.preads = max((int)ceil(mean),1);
		}
		//if there is more than 1 read on the positive strand PCR correct the current position
		if(currentPositionPCRcorrected.mreads >1){
			double mean = ((double)(windowMreads-currentPositionPCRcorrected.mreads))/bandwidth;
			if(currentPositionPCRcorrected.mreads > mean && lPois(currentPositionPCRcorrected.mreads,mean)<log(0.01))
				currentPositionPCRcorrected.mreads = max((int)ceil(mean),1);
		}

		//write the PCR corrected position
		PCRcorrectDepthFile.write((char*)(&currentPositionPCRcorrected), sizeof(currentPositionPCRcorrected));

		//increment position
		++ind;
	}
	depthFile.close();
	PCRcorrectDepthFile.close();
	fprintf(stderr, "Finished PCR correcting depth graph\n");
}

