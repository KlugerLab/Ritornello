/*
 * Artifact.h
 *
 *  Created on: May 5, 2015
 *      Author: Kelly
 */

#ifndef SRC_ARTIFACT_H_
#define SRC_ARTIFACT_H_
#include "Peaks.h"

class Artifact {
public:
	Artifact();
	virtual ~Artifact();
	static double test(long windowSize, long readLength, double* fir, Peak cadidate, double testRatio, long halfLength);
};

#endif /* SRC_ARTIFACT_H_ */
