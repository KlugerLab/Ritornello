/*
 * LikelihoodRatioTest.h
 *
 *  Created on: Jan 16, 2015
 *      Author: kps3
 */

#ifndef LikelihoodRatioTest_H_
#define LikelihoodRatioTest_H_
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace UBL = boost::numeric::ublas;

class LikelihoodRatioTest{
public:
	LikelihoodRatioTest();
	LikelihoodRatioTest(long argMaxFragmentLength,double* argPeakFilter, double* fld,bool argModelOpenChromatinEffect);
	virtual ~LikelihoodRatioTest();
	void calcDesignMatrices(int numLocalPeaks, const long* localPeaks, double* pstrand, double* mstrand);
	void DesignMatricesAddOpenChromatinEffect(int numLocalPeaks, const long* localPeaks, double* pstrandExtended, double* mstrandExtended);
	double ll(double* strand, const UBL::matrix<double>& filterMat, UBL::vector<double> theta);
	UBL::vector<double> gradient(double* strand, const UBL::matrix<double>& filterMat, UBL::vector<double> theta);
	UBL::matrix<double> hessian(double* strand, const UBL::matrix<double>& filterMat,const UBL::matrix<double>& tfilterMat, UBL::vector<double> theta);
	int maximizeLL(double* strand, const UBL::matrix<double>& filterMat,const UBL::matrix<double>& tfilterMat,UBL::vector<double>& theta, double& maxll);
	int test(double* pstrand, double* mstrand,double* pstrandExtended, double* mstrandExtended, const long numLocalPeaks, const long* localPeaks);
	inline void constrain(UBL::vector<double>& theta, double reads);
	double peakReads;
	double lpValue;
	double lpValueOpenChromatin;
	double llAlt;
private:
	long maxFragmentLength;
	double* peakFilter;
	double* fld;
	bool modelOpenChromatinEffect;
	UBL::vector<double> uFilter;
	UBL::matrix<double> pFilterMat;
	UBL::matrix<double> mFilterMat;
	UBL::matrix<double> tpFilterMat;
	UBL::matrix<double> tmFilterMat;
};

#endif /* LikelihoodRatioTest_H_ */
