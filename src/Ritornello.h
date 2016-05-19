/*
 * Ritornello.h
 *
 *  Created on: Dec 19, 2014
 *      Author: kps3
 */

#ifndef RITORNELLO_H_
#define RITORNELLO_H_
#include <math.h>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include "Config.h"
#include "FFTHandler.h"
#include "IOhandler.h"
#include "LikelihoodRatioTest.h"
#include "Peaks.h"

class Ritornello {
public:
	Ritornello(Config argParms);
	virtual ~Ritornello();
	void readFullGenomeCoverage();
	void calcFLD();
	void calculateFLD();
	void calculateFIR();
	void calculatePeakFilter();
	void testPeakProjections(double artifactTestRatio);
	void writeResults();
	void readFIRFromFile();
	void readFLDFromFile();
	void calculateMaxFragmentLength();
	void removeArtifacts();
	void estimateTrainingPeaks();
	void readTrainingPeaksFromFile();
	void calcWindowReadCountDist();
	void FDRcorrect();
	void findPeaks(double peakThresholdRatio);
private:
	long _getPeakNumber(double* array);
	Config parms;
	Peaks peaks;
	double* pstrand;
	double* mstrand;
	double* fld;
	double* fir;
	long printLength;
	double* peakFilter;
	double* candidatePeaks;
	double* rlArtifacts;
	double* likelihoodRatioTestScores;
	double* likelihoodPeakCounts;
	long readLength;
	long genomeLength;
	vector<double> windowReadCountDist;
	long medianWindowReadCount;
};
#endif /* RITORNELLO_H_ */
