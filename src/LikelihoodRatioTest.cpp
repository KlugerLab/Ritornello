/*
 * LikelihoodRatioTest.cpp
 *
 *  Created on: Jan 16, 2015
 *      Author: kps3
 */
#include "LikelihoodRatioTest.h"
#include <math.h>
#include <stdlib.h>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
//#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <float.h>
#include "IOhandler.h"
#include <algorithm>
#include <fstream>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/distributions/normal.hpp>

using namespace boost::math;
using namespace boost::numeric::ublas;
using boost::numeric::ublas::matrix_range;

/*
#include <iostream>

using namespace std;

// C++ calls functions in a different way, so you need to change specify that this is a C/FORTRAN function call

extern "C"{
// FORTRAN adds _ after all the function names
// and all variables are called by reference
double ddot_( const int *N, const double *a, const int *inca, const double *b, const int *incb );
}

double ddot( int N, double *a, int inca, double *b, int incb ){
  return ddot_( &N, a, &inca, b, &incb );
};
*/


LikelihoodRatioTest::LikelihoodRatioTest(){
	maxFragmentLength = 0;
	peakFilter = 0;
	fld=0;
	peakReads = 0;
	lpValue = 0;
	lpValueOpenChromatin = 0;
	llAlt = NAN;
	modelOpenChromatinEffect=true;
}
LikelihoodRatioTest::LikelihoodRatioTest(long argMaxFragmentLength,double* argPeakFilter, double* argFld, bool argModelOpenChromatinEffect){
	maxFragmentLength = argMaxFragmentLength;
	peakFilter = new double[maxFragmentLength];
	fld = new double[maxFragmentLength];
	uFilter.resize(2*maxFragmentLength);
	peakReads = 0;
	lpValue = 0;
	lpValueOpenChromatin = 0;
	llAlt = NAN;
	modelOpenChromatinEffect=argModelOpenChromatinEffect;

	double pSum = 0;
	double uSum = 0;
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		uFilter[ii] = 1;
		uSum += uFilter[ii];
	}
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii)
		uFilter[ii]/=uSum;

	for(long ii = 0; ii < maxFragmentLength; ++ii){
		peakFilter[ii]=argPeakFilter[ii];
		pSum += peakFilter[ii];
	}
	for(long ii = 0; ii < maxFragmentLength; ++ii)
		peakFilter[ii] /= pSum;

	for(long ii = 0; ii < maxFragmentLength; ++ii)
		fld[ii] = argFld[ii];
}
LikelihoodRatioTest::~LikelihoodRatioTest() {
	delete[] peakFilter;
	delete[] fld;
}

void LikelihoodRatioTest::calcDesignMatrices(int numLocalPeaks, const long* localPeaks, double* pstrandExtended, double* mstrandExtended){
	int n = 2*maxFragmentLength;
	//set up model
	pFilterMat.resize(n,numLocalPeaks+1);
	mFilterMat.resize(n,numLocalPeaks+1);

	//add uniform filter
	for(int ii = 0; ii <n;++ii){
		pFilterMat(ii,0) = uFilter[ii];
		mFilterMat(ii,0) = uFilter[ii];
	}

	//add peak filters
	for(int jj = 1; jj <numLocalPeaks+1;++jj){
		long offset = localPeaks[jj-1];
		for(int ii = 0; ii <n;++ii){
			int pind = -ii+offset;
			pFilterMat(ii,jj) = (pind>=0&&pind<maxFragmentLength)?peakFilter[pind]:0;
			int mind = ii-offset;
			mFilterMat(ii,jj) = (mind>=0&&mind<maxFragmentLength)?peakFilter[mind]:0;
		}
	}

	if(modelOpenChromatinEffect){
		DesignMatricesAddOpenChromatinEffect(numLocalPeaks,localPeaks,pstrandExtended,mstrandExtended);
	}

	tpFilterMat = UBL::trans(pFilterMat);
	tmFilterMat = UBL::trans(mFilterMat);
}

void LikelihoodRatioTest::DesignMatricesAddOpenChromatinEffect(int numLocalPeaks, const long* localPeaks, double* pstrandExtended, double* mstrandExtended){
	int n = 2*maxFragmentLength;
	//set up model
	pFilterMat.resize(n,numLocalPeaks+2);
	mFilterMat.resize(n,numLocalPeaks+2);

	// add open chromatin effect filter
	boost::math::normal norm(0,50);
	double pg[4*maxFragmentLength];
	double mg[4*maxFragmentLength];

	/*
	 * //convolve the positive and negative strands against the fld for the predicted null coverage of the opposing strand
	double pu[4*maxFragmentLength];
	double mu[4*maxFragmentLength];

	for(int ii = 0; ii < 4*maxFragmentLength; ++ii){
		pu[ii]=0;
		mu[ii]=0;
		for(int jj = 0; jj < 4*maxFragmentLength; ++jj){
			int pind = -ii+jj;
			if(pind > 0 && pind < maxFragmentLength)
				pu[ii]+=mstrandExtended[jj]* fld[pind];
			int mind = ii-jj;
			if(mind > 0 && mind < maxFragmentLength)
				mu[ii]+=pstrandExtended[jj]* fld[mind];
		}

	}
	 */
	//convolve with guassian
	for(int ii = maxFragmentLength; ii < 4*maxFragmentLength; ++ii){
		pg[ii]=0;
		mg[ii]=0;

		for(int jj = 0; jj < 4*maxFragmentLength; ++jj){
			double gaus = pdf(norm,ii-jj);
			pg[ii]+=pstrandExtended[jj]* gaus;
			mg[ii]+=mstrandExtended[jj]* gaus;
		}
	}
	//normalize
	double pSum = 0;
	double mSum = 0;
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		pSum += pg[ii+maxFragmentLength];
		mSum += mg[ii+maxFragmentLength];
	}
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		pg[ii+maxFragmentLength]/=pSum;
		mg[ii+maxFragmentLength]/=mSum;
	}
	//IOhandler::printDoubleArrayToFile(&pg[maxFragmentLength],2*maxFragmentLength,"pg.txt");
	//IOhandler::printDoubleArrayToFile(&mg[maxFragmentLength],2*maxFragmentLength,"mg.txt");
	//add open chr effect filter
	for(int ii = 0; ii <n;++ii){
		pFilterMat(ii,pFilterMat.size2()-1) = pg[ii+maxFragmentLength];
		mFilterMat(ii,mFilterMat.size2()-1) = mg[ii+maxFragmentLength];
	}
}

//log Poisson distribution pdf
inline double lPois(long k, double mu){
	return k*log(mu) -mu -lgamma(k+1);
}
//first partial with respect to m
inline double dlPoisdMu(long k, double mu){
	return k/mu-1.0;
}
//second partial with respect to m
inline double ddlPoisddMu(long k, double mu){
	return -k/(mu*mu);
}

inline void LikelihoodRatioTest::constrain(UBL::vector<double>& beta, double reads){
	static const double minValB = 0.01;
	static const double maxValB = 2*reads;

	//beta0
	if(beta[0] < minValB)
		beta[0]=minValB;
	if(beta[0] > maxValB)
		beta[0]=maxValB;
	//betaN
	for(unsigned int ii = 1; ii < beta.size(); ++ii){
		if(beta[ii] < 0)
			beta[ii]=0;
		if(beta[ii] > 2*reads)
			beta[ii]=2*reads;
	}
}
/*
int testBlas( int argc, char** argv ){
  // you can define the arrays in one of two ways
  // on the heap
  double *a = new double[3];
  a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
  // on the stack
  double b[3] = { 4.0, 5.0, 6.0 };

  double dot_product = ddot( 3, a, 1, b, 1 );
  cout <<" The dot product is: " <<  dot_product << endl;

  return 0;
};
*/
inline int ginv(UBL::matrix<double>& mat){
	UBL::matrix<double> argMat=mat;
	// perform LU-factorization
	typedef permutation_matrix<std::size_t> pmatrix;
	pmatrix pm(argMat.size1());
	int res = lu_factorize(argMat,pm);
	if( res != 0 )
		return -1;
	// create identity matrix of "inverse"
	mat.assign(UBL::identity_matrix<double>(argMat.size1()));
	// backsubstitute to get the inverse
	swap_rows (pm, mat);
	UBL::matrix<double> cm1 (mat);
	inplace_solve (argMat, mat, unit_lower_tag ());
	bool check1 = UBL::detail::expression_type_check (prod (triangular_adaptor<UBL::matrix<double>, unit_lower> (argMat), mat), cm1);
	UBL::matrix<double> cm2 (mat);
	inplace_solve (argMat, mat, upper_tag ());
	bool check2 = UBL::detail::expression_type_check (prod (triangular_adaptor<UBL::matrix<double>, upper> (argMat), mat), cm2);

	if(check1 && check2)
		return 0;
	else
		return -1;
}
double LikelihoodRatioTest::ll(double* strand, const UBL::matrix<double>& filterMat, UBL::vector<double> beta){
	UBL::vector<double> lambda = prod(filterMat,beta);
	double ll = 0;
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		if(lambda[ii] <= 0) lambda[ii]=DBL_MIN;
		ll+= lPois(strand[ii],lambda[ii]);
	}
	return ll;
}
UBL::vector<double> LikelihoodRatioTest::gradient(double* strand, const UBL::matrix<double>& filterMat, UBL::vector<double> beta){
	UBL::vector<double> lambda = prod<UBL::vector<double> >(filterMat,beta);
	UBL::vector<double> grad(2*maxFragmentLength,0.0);
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		grad[ii] = dlPoisdMu(strand[ii], lambda[ii]);
	}
	grad = prod(grad,filterMat);
	return grad;
}
UBL::matrix<double> LikelihoodRatioTest::hessian(double* strand, const UBL::matrix<double>& filterMat, const UBL::matrix<double>& tfilterMat, UBL::vector<double> beta){
	//setup up hess for return
	UBL::matrix<double> hess(beta.size(),beta.size(),0.0);
	UBL::vector<double> lambda = prod<UBL::vector<double> >(filterMat,beta);
	UBL::vector<double> ddm(2*maxFragmentLength);
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		ddm(ii) = ddlPoisddMu(strand[ii], lambda[ii]);
	}
	//do row wise product to avoid large diagonal matrix
	UBL::matrix<double> tmp = filterMat;
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		matrix_row<matrix<double> > row(tmp, ii);
		row*=ddm[ii];
	}
	hess=prod(tfilterMat,tmp);
	return hess;
}
int LikelihoodRatioTest::maximizeLL(double* strand, const UBL::matrix<double>& filterMat, const UBL::matrix<double>& tfilterMat, UBL::vector<double>& beta, double& maxll){
	int convergedSteps=0;
	double reads = 0;
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii)
		reads+=strand[ii];

	//constrain theta
	constrain(beta,reads);
	maxll = ll(strand, filterMat, beta);
	double error=0.01;
	//parameterize the  model
	for(int ii = 0; ii < 100000; ++ii){
		//initialize the step
		double alpha = 1;
		//calculate the gradient
		UBL::vector<double> grad = gradient(strand,filterMat,beta);
		//calculate the step
		UBL::vector<double> step;
		//Dogleg optimization uses gradient descent for the first couple of iterations then switch to newton
		if(ii <3){
			step = -grad;
		}
		else{
			//calculate the hessian
			UBL::matrix<double> hess = hessian(strand,filterMat,tfilterMat,beta);
			//invert the hessian
			if(ginv(hess)==-1){
				//fprintf(stderr, "Non-invertible hessian switching to gradient descent\n");
				step=-grad;
			}
			else
				step = prod(hess,grad);
		}
		double newLL = -DBL_MAX;
		UBL::vector<double> newBeta;
		//calculate new alpha and parameters according to half stepping
		do {
			//get new model parameters
			newBeta = beta - alpha*step;
			//if the beta params havn't converged hold R still
			if(convergedSteps<1)
				newBeta[newBeta.size()-1] = beta[beta.size()-1];
			constrain(newBeta,reads);
			newLL = ll(strand, filterMat, newBeta);
			alpha/=2;
		} while(newLL<maxll);

		//check for convergence
		double stepError = abs(newLL - maxll);
		if(stepError<error){
			++convergedSteps;
		}else
			convergedSteps=0;
		//update the parameters
		beta = newBeta;
		maxll = newLL;
		//stop if we converged
		if(convergedSteps>2)
			return 0;
	}
	return -1;
}

int LikelihoodRatioTest::test(double* pstrand, double* mstrand,double* pstrandExtended, double* mstrandExtended, const long numLocalPeaks, const long* localPeaks){
	int converged = 0;
	//count reads
	double preads = 0;
	double mreads = 0;
	for(long ii = 0; ii < 2*maxFragmentLength; ++ii){
		preads+=pstrand[ii];
		mreads+=mstrand[ii];
	}
	//IOhandler::printDoubleArrayToFile(pstrand,2*maxFragmentLength,"pstrandTest.txt");
	//IOhandler::printDoubleArrayToFile(mstrand,2*maxFragmentLength,"mstrandTest.txt");

	//Set up the null model
	long numNullLocalPeaks = numLocalPeaks-1;
	long nullLocalPeaks[numNullLocalPeaks];
	//we want to exclude the center peak that we are testing for the null model which is the first
	//peak in numLocalPeaks
	for(int ii = 0; ii < numNullLocalPeaks; ++ii)
		nullLocalPeaks[ii] = localPeaks[ii+1];
	//calculate design matrices for the null model
	calcDesignMatrices(numNullLocalPeaks,nullLocalPeaks, pstrandExtended, mstrandExtended);

	//init beta0 to read counts and beta1 to no reads
	UBL::vector<double> pBetaNull(pFilterMat.size2(),0.0);
	UBL::vector<double> mBetaNull(mFilterMat.size2(),0.0);
	pBetaNull[0]=preads;
	mBetaNull[0]=mreads;

	//optimize the null model likelihood
	double llNullp = 0;
	double llNullm = 0;
	converged += maximizeLL(pstrand,pFilterMat,tpFilterMat, pBetaNull,llNullp);
	converged += maximizeLL(mstrand,mFilterMat,tmFilterMat, mBetaNull,llNullm);
	double llNull = llNullp + llNullm;

	//calculate design matrices for the alt model
	calcDesignMatrices(numLocalPeaks,localPeaks, pstrandExtended, mstrandExtended);
	//init beta0 to read counts and beta1 to no reads
	UBL::vector<double> pBetaAlt(pFilterMat.size2(),0.0);
	UBL::vector<double> mBetaAlt(mFilterMat.size2(),0.0);
	//initialize the parameters to what we learned for the null model
	pBetaAlt[0]=pBetaNull[0];
	mBetaAlt[0]=mBetaNull[0];
	for(unsigned int ii = 2; ii < pBetaAlt.size(); ++ii){
		pBetaAlt[ii]=pBetaNull[ii-1];
		mBetaAlt[ii]=mBetaNull[ii-1];
	}

	//optimize the alternative model likelihood
	double llAltp = 0;
	double llAltm = 0;
	converged += maximizeLL(pstrand,pFilterMat,tpFilterMat, pBetaAlt,llAltp);
	converged += maximizeLL(mstrand,mFilterMat,tmFilterMat, mBetaAlt,llAltm);
	llAlt = llAltp + llAltm;

	peakReads= pBetaAlt[1]+mBetaAlt[1];
	//check to see if we converged (do we have NaNs)
	if(llNull==llNull && llAlt==llAlt && !std::isinf(llNull) && !std::isinf(llAlt)){
		double lrtStatistic = -2*llNull + 2* llAlt;
		lpValue = lrtStatistic/2.0/log(10);
	}
	else{
		lpValue=0;
	}
/*
	//Set up the open chromatin model
	//calculate design matrices for the null model
	calcDesignMatricesOpenChromatin(numNullLocalPeaks,nullLocalPeaks, pstrandExtended, mstrandExtended);

	//init betaOpenChromatin to beta0
	//init beta0 to read counts and beta1 to no reads
	UBL::vector<double> pBetaOpenChromatin(numNullLocalPeaks+2,0.0);
	UBL::vector<double> mBetaOpenChromatin(numNullLocalPeaks+2,0.0);
	//initialize the parameters to what we learned for the null model
	pBetaOpenChromatin[0]=pBetaNull[0];
	mBetaOpenChromatin[0]=mBetaNull[0];
	for(int ii = 1; ii < numNullLocalPeaks+1; ++ii){
		pBetaOpenChromatin[ii]=pBetaNull[ii];
		mBetaOpenChromatin[ii]=mBetaNull[ii];
	}

	//optimize the OpenChromatin model likelihood
	double llOpenChromatinp = 0;
	double llOpenChromatinm = 0;
	converged += maximizeLL(pstrand,pFilterMat,tpFilterMat, pBetaOpenChromatin,llOpenChromatinp);
	converged += maximizeLL(mstrand,mFilterMat,tmFilterMat, mBetaOpenChromatin,llOpenChromatinm);
	double llOpenChromatin = llOpenChromatinp + llOpenChromatinm;

	if(llNull==llNull && llOpenChromatin==llOpenChromatin && !std::isinf(llNull) && !std::isinf(llOpenChromatin)){
		double lrtStatistic = -2*llNull + 2* llOpenChromatin;
		lpValueOpenChromatin = lrtStatistic/2.0/log(10);
	}
	else{
		lpValueOpenChromatin=0;
	}
	*/
	return converged;
}

