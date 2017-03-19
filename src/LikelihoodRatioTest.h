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
	void calcDesignMatrices(int numLocalPeaks, const long* localPeaks, unsigned short* pstrand, unsigned short* mstrand);
	void DesignMatricesAddOpenChromatinEffect(int numLocalPeaks, const long* localPeaks, unsigned short* pstrandExtended, unsigned short* mstrandExtended);
	double ll(unsigned short* strand, const UBL::matrix<double>& filterMat, UBL::vector<double> theta);
	UBL::vector<double> gradient(unsigned short* strand, const UBL::matrix<double>& filterMat, UBL::vector<double> theta);
	UBL::matrix<double> hessian(unsigned short* strand, const UBL::matrix<double>& filterMat,const UBL::matrix<double>& tfilterMat, UBL::vector<double> theta);
	int maximizeLL(unsigned short* strand, const UBL::matrix<double>& filterMat,const UBL::matrix<double>& tfilterMat,UBL::vector<double>& theta, double& maxll);
	int test(unsigned short* pstrand, unsigned short* mstrand,unsigned short* pstrandExtended, unsigned short* mstrandExtended, const long numLocalPeaks, const long* localPeaks);
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
