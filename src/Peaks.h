/*
 * Peaks.h
 *
 *  Created on: May 5, 2015
 *      Author: Kelly
 */

#ifndef SRC_PEAKS_H_
#define SRC_PEAKS_H_
#include <vector>

using namespace std;
struct Peak{
	int chr;
	int pos;
	double lpValue;
	double lpValueOpenChromatin;
	double lqValue;
	double betaAlt;
	bool isArtifact;
	double artifactScore;
	unsigned short* pstrandExtended;
	unsigned short* mstrandExtended;
	unsigned short* pstrand;
	unsigned short* mstrand;
	vector<long> localPeakPos;
	vector<Peak*> localPeakPtr;
};

class Peaks {
public:
	Peaks();
	virtual ~Peaks();
	void clear();
	void find(char* bamFileName, int argWindowSize, double* fir, double peakThreshold, long minReadsThreshold, bool correctPCR);
	void findPeakNeighbors();
	void addPeak(long chr, long position, double* pstrand, double* mstrand);
	vector<Peak> data;
	long windowSize;
	long positionsScanned;
};

#endif /* SRC_PEAKS_H_ */
