/*
 * CleanFLD.h
 *
 *  Created on: Sep 22, 2014
 *      Author: Kelly Patrick Stanton
 */

#ifndef FLD_H_
#define FLD_H_

namespace FLD{
double* calculate(char* bamFileName, int windowSize, bool correctPCR);
double* calculatePairedEnd(char* bamFileName, int windowSize);
double* calculateSingleEnd(char* bamFileName, int windowSize, bool correctPCR);
}
#endif /* FLD_H_ */
