/*
 * FIR.cpp
 *
 *  Created on: Apr 10, 2015
 *      Author: kps3
 */

#include "FIR.h"
#include "IOhandler.h"
#include "float.h"
#include "FFTHandler.h"
#include <fstream>
#include <string>
#include "LikelihoodRatioTest.h"
#include <math.h>
#include <omp.h>

namespace FIR {

double pK(double parm, double k)
{
	if(k >= 1 || k <= 0) return 0;
	return pow(k, parm-1)*pow(1-k, parm-1);
}

double* calculateFilter(double parm, long windowSize, double* fld){
	double* filterp = new double[windowSize];
	double sum = 0;
	for(int ii = 0; ii <windowSize; ++ii)
	{
		filterp[ii] = 0;
		for(int jj = 0; jj < windowSize; ++jj)
		{
			filterp[ii] += fld[jj+1]*pK(parm,((double)(ii+1))/(jj+1))/abs((jj+1));
		}
		sum+=filterp[ii];
		//sum+=pow(filterp[ii],2.0);
	}
	//sum=sqrt(sum);
	for(int ii = 0; ii<windowSize;++ii){
		filterp[ii]/=sum;
	}
	for(long ii=windowSize; ii < windowSize; ++ii){
			filterp[ii]=0;
	}
	return filterp;
}

double calculateAlpha(const Peaks& peaks, long windowSize, double* fld, int numThreads){
	fprintf(stderr, "Parameterizing alpha\n");
	int numTestPeaks=min(200, (int)peaks.data.size());
	int numNonConvergentTests = 0;

	double bestAlpha=-1;
	double bestAlphall=-DBL_MAX;
	//try a bunch of alpha's
	for(int ii = 9; ii < 31; ++ii){
		double alpha = ii/10.0 +0.1;
		double alphall = 0;
		//create an array of tests 1 for each thread
		LikelihoodRatioTest* likelihoodRatioTest[numThreads];
		//call constructor on each
		for(int ii = 0; ii < numThreads; ++ii){
			likelihoodRatioTest[ii] = new LikelihoodRatioTest(windowSize,calculateFilter(alpha, windowSize, fld), fld, false);
		}
		//for each peak
		omp_lock_t writelock;
		omp_init_lock(&writelock);
		#pragma omp parallel for reduction(+:alphall)
		for(int jj = 0; jj < numTestPeaks; ++jj){
			int tid = omp_get_thread_num();
			int ret = likelihoodRatioTest[tid]->test(peaks.data[jj].pstrand,peaks.data[jj].mstrand,
					peaks.data[jj].pstrandExtended,peaks.data[jj].mstrandExtended,
					peaks.data[jj].localPeakPos.size(),&peaks.data[jj].localPeakPos[0]);
			if(ret!=0){
				omp_set_lock(&writelock);
				++numNonConvergentTests;
				omp_unset_lock(&writelock);
				//FIXME this should not be zero in the case of non-convergence
				alphall=0;
			}
			else
				alphall += likelihoodRatioTest[tid]->llAlt;
		}
		//delete all test objects
		for(int ii = 0; ii < numThreads; ++ii){
			delete likelihoodRatioTest[ii];
		}
		fprintf(stderr, "Alpha parameter [%f] gave likelihood [%f]\n", alpha, alphall);
		if(alphall>bestAlphall){
			bestAlphall = alphall;
			bestAlpha=alpha;
		}
	}

	fprintf(stderr, "Best alpha parameter is [%f]\n", bestAlpha);
	fprintf(stderr, "Non-convergent likelihood tests[%d]\n", numNonConvergentTests);
	return bestAlpha;
}

} /* namespace FIR */
