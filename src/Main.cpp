// Author      : Kelly Patrick Stanton
// Version     :
// Copyright   : Copyright Kluger Lab
// Description : Ritornello
//============================================================================
#include "Config.h"
#include "Ritornello.h"
#include "SamStream.h"
#include "BufferedGenomeReader.h"
#include "PCRCorrectGenomeReader.h"
#include "BufferedDepthGraphReader.h"
#include "DepthGraph.h"
#include <string>
#include <fstream>
using namespace std;
/**
 * Main entry point for the program
 *
 * Ritornello starts here.  It first calculates a suitable fragment length distribution using artifact masking
 * and Fourier deconvolution.  Next it identifies peaks according to a filter constructed form the predicted
 * binding pattern based on the measured fragment length distribution (FLD).  Finally it calls peaks and tests
 * each one to determine if it is actually a read length artifact instead.
 *
 * @param argc number of command line parameters
 * @param argv array of characters pointers.  Each string is a command line argument.
 * @return
 */
int main(int argc, char *argv[]){
	//get parameters from the command line
	Config parms;
	parms.parseCmdLine(argc,argv);
	if(parms.getDebugFolder()){
		IOhandler::outputFolder=parms.getDebugFolder();
	}

	DepthGraph::outPrefix =parms.getOutputPrefix();
	DepthGraph depthgraph;
	//Check if we need to convert the sam to a Depth Graph
	ifstream f1((string(parms.getOutputPrefix())+".RitorDepthGraph").c_str(),ios::in);
	ifstream f2((string(parms.getOutputPrefix())+".RitorDepthGraphKey").c_str(),ios::in);
	if(!f1.good()||!f2.good())
		depthgraph.Sam2DepthGraph(parms.getBamFileName());
	ifstream f3((string(parms.getOutputPrefix())+"-PCRCorrect.RitorDepthGraph").c_str(),ios::in);
	if(!f3.good())
		depthgraph.PCRcorrect(20);
	f1.close();
	f2.close();
	f3.close();



	/*
	//create a buffered reader to get the chromosome names
	//BufferedGenomeReader bgr(1024);
	PCRCorrectGenomeReader bgr(1024);
	bgr.init(parms.getBamFileName());
	BufferedGenomeReader bdr(1044);
	bdr.init(parms.getBamFileName());
	//BufferedDepthGraphReader bdr(1024);
	//bdr.init(parms.getOutputPrefix());
	//for (int ii = 0; ii < 30; ++ii)
	//	bgr.next();

	int count = 0;
	while(bdr.next()){
		bgr.next();
		double* bgrpstrand = bgr.getPstrand();
		double* bdrpstrand = bdr.getPstrand();
		double* bgrmstrand = bgr.getMstrand();
		double* bdrmstrand = bdr.getMstrand();
		for(int ii =0; ii < 1024; ++ii){
			if(bgrpstrand[ii] != bdrpstrand[ii] || bgrmstrand[ii] != bdrmstrand[ii])
				fprintf(stderr,"strands not equal");
		}
		++count;
	}
	*/

	//initialize program
	Ritornello ritornello(parms);

	//calculate fragment length distribution
	if(parms.getFLDFile()==0 && parms.getFilterFile()==0)
		ritornello.calculateFLD();
	else
		ritornello.readFLDFromFile();
	ritornello.calculateMaxFragmentLength();

	//Get average read coverage
	ritornello.calcWindowReadCountDist();

	//calculate filter shape
	if(parms.getFilterFile()==0){
		ritornello.estimateTrainingPeaks();
		ritornello.calculateFIR();
	}
	else
		ritornello.readFIRFromFile();

	if(parms.getTRAINFile()==0){
		ritornello.findPeaks(1.0);
		ritornello.savePeaks();
	}
	else
		ritornello.readTrainingPeaksFromFile();

	//test for peaks
	ritornello.testPeakProjections(1.0);
	//write peaks to file
	ritornello.writeResults();

    fprintf(stderr, "Finished!!!\n");
    exit(EXIT_SUCCESS);
}
