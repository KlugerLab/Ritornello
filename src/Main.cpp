// Author      : Kelly Patrick Stanton
// Version     :
// Copyright   : Copyright Kluger Lab
// Description : Ritornello
//============================================================================
#include "Config.h"
#include "Ritornello.h"
#include "SamStream.h"
#include "BufferedGenomeReader.h"
#include "DepthGraph.h"
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

	DepthGraph depthgraph;
	depthgraph.Sam2DepthGraph(parms.getBamFileName(),parms.getOutputPrefix());
	depthgraph.PCRcorrect(parms.getOutputPrefix(),20);

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
