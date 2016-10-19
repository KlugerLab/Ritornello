/*
 * IOhandler.h
 *
 *  Created on: Jun 19, 2014
 *      Author: Kelly Patrick Stanton
 */

#ifndef IOHANDLER_H_
#define IOHANDLER_H_

#include <string>
#include <vector>
#include "Peaks.h"
using namespace std;

class IOhandler {
public:
	IOhandler();
	virtual ~IOhandler();
	static string outputFolder;
	static vector<double> readDoubleArrayFromFile(const char* fileName );
	static void printDoubleArrayToFile(double* array, int length, const char* fileName );
	static void debugPeaks(long windowSize, const Peaks& peaks);
	static vector<Peak> readPeaks(long windowSize, const char* fileName);
};

#endif /* IOHANDLER_H_ */
