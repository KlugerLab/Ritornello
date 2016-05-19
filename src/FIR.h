/*
 * FIR.h
 *
 *  Created on: Apr 10, 2015
 *      Author: kps3
 */

#ifndef FIR_H_
#define FIR_H_
#include "Peaks.h"

namespace FIR {
double* calculateFilter(double parm, long windowSize, double* fld);
double calculateAlpha(const Peaks& peaks, long windowSize, double* fld, int numThreads);
double nextTestParmFract(long ii, double* peakScores,double*peakParms);
}
#endif /* FIR_H_ */
