/*
 * Ritornello.cpp
 *
 *  Created on: Dec 19, 2014
 *      Author: kps3
 */

#include "Ritornello.h"
#include "BufferedDepthGraphReader.h"
#include "FLD.h"
#include "FIR.h"
#include "Artifact.h"
#include <fstream>
#include <float.h>
#include <omp.h>

Ritornello::Ritornello(Config argParms) {
	parms=argParms;
	pstrand=NULL;
	mstrand=NULL;
	fld=NULL;
	printLength=0;
	peakFilter=NULL;
	candidatePeaks=NULL;
	rlArtifacts=NULL;
	omp_set_num_threads(argParms.getNumThreads());
}

Ritornello::~Ritornello() {
}

void Ritornello::calculateMaxFragmentLength(){
	double sum = 0;
	//if MFL is not too far set it
	for(long ii =0; ii < 0.75*parms.getMaxFragmentLength(); ++ii)
	{
		sum+= fld[ii];
		if(sum >= 0.75) {
			parms.setMaxFragmentLength(ii/0.5);
			break;
		}
	}
	fprintf(stderr, "Estimated maximum fragment length [%ld bp]\n", parms.getMaxFragmentLength());
}

void Ritornello::readFIRFromFile(){
	vector<double> peakFilterVector = IOhandler::readDoubleArrayFromFile(parms.getFilterFile());
	parms.setMaxFragmentLength(peakFilterVector.size());
	fir = new double[parms.getMaxFragmentLength()];
	for(int ii=0;ii < parms.getMaxFragmentLength(); ++ii)
		fir[ii]=peakFilterVector[ii];
}

void Ritornello::readFLDFromFile(){
	vector<double> peakFLDVector = IOhandler::readDoubleArrayFromFile(parms.getFLDFile());
	parms.setMaxFragmentLength(peakFLDVector.size());
	fld = new double[parms.getMaxFragmentLength()];
	for(int ii=0;ii < parms.getMaxFragmentLength(); ++ii)
		fld[ii]=peakFLDVector[ii];
}
void Ritornello::findPeaks(double peakThresholdRatio){
	peaks.clear();
	long minReadTheshold = max(2*medianWindowReadCount,(long)parms.getMinReadThreshold());
	//long minReadTheshold = (long)parms.getMinReadThreshold();
	peaks.find(parms.getBamFileName(),parms.getMaxFragmentLength(),fir,parms.getPeakThreshold(),peakThresholdRatio*minReadTheshold,
			parms.getCorrectPCR());
}

void Ritornello::calculateFLD(){
	//Get clean FLD
	fld = FLD::calculate(parms.getBamFileName(), parms.getMaxFragmentLength(), parms.getCorrectPCR());
	IOhandler::printDoubleArrayToFile(fld,parms.getMaxFragmentLength(),"fld.txt");
}
void Ritornello::calculateFIR(){
	//Learn FIR alpha from the training peaks
	double alpha = FIR::calculateAlpha(peaks, parms.getMaxFragmentLength(), fld, parms.getNumThreads());
	//clear training peaks
	peaks.clear();
	fprintf(stderr, "K ~ Beta(alpha,alpha), alpha parameter=[%f]\n",alpha);
	//construct alpha
	fir = FIR::calculateFilter(alpha, parms.getMaxFragmentLength(), fld);
	IOhandler::printDoubleArrayToFile(fir,parms.getMaxFragmentLength(),"fir.txt");
}
bool comparepValue(Peak i,Peak j) { return (i.lpValue>j.lpValue); }
bool compareBetaAlt(Peak i,Peak j) { return (i.betaAlt>j.betaAlt); }
void Ritornello::estimateTrainingPeaks(){
	//guess alpha to be 1
	fir = FIR::calculateFilter(1.0, parms.getMaxFragmentLength(), fld);

	findPeaks(2.0);

	//calculate peaks based on this guess
	testPeakProjections(1.0);
	delete fir;

	//remove artifacts
	removeArtifacts();
	fprintf(stderr, "Sorting by effect size\n");
	std::sort (peaks.data.begin(), peaks.data.end(), compareBetaAlt);
	fprintf(stderr, "Sorting by pvalue\n");
	std::sort (peaks.data.begin(), peaks.data.end(), comparepValue);
}

void Ritornello::savePeaks(){
	//print out peaks
	IOhandler::debugPeaks(parms.getMaxFragmentLength(), peaks);
}

void Ritornello::readTrainingPeaksFromFile(){
	peaks.data = IOhandler::readPeaks(parms.getMaxFragmentLength(), parms.getTRAINFile());
	peaks.windowSize = parms.getMaxFragmentLength();
}


void Ritornello::testPeakProjections(double artifactTestRatio){
	BufferedDepthGraphReader bgr(1);
	bgr.init();
	readLength = bgr.readLength;
	genomeLength = bgr.genomeLength;
	bgr.close();

	int numNonConvergentTests=0;
	//Calculate artifact test window as 50% of filter
	double firSum = 0;
	for(int ii =0; ii < parms.getMaxFragmentLength(); ii++)
		firSum+=fir[ii];
	double cumsum = 0;
	long halfLength;
	for(halfLength =0; halfLength < parms.getMaxFragmentLength(); halfLength++){
		cumsum+=fir[halfLength]/firSum;
		if(cumsum > 0.5) break;
	}
	//create an array of tests 1 for each thread
	LikelihoodRatioTest* likelihoodRatioTest[parms.getNumThreads()];
	//call constructor on each
	for(int ii = 0; ii < parms.getNumThreads(); ++ii){
		likelihoodRatioTest[ii] = new LikelihoodRatioTest(parms.getMaxFragmentLength(),fir, fld,parms.getModelOpenChromatinEffect());
	}

	Artifact::init();

	omp_lock_t writelock;
	omp_init_lock(&writelock);

	fprintf(stderr, "Testing peaks for artifacts and significance\n");
	fprintf(stderr, "|--------------------------------------------------------------------------------------------------|\n");
	long progress = 0;
	unsigned int percent = 0;
#ifndef DEBUG
	#pragma omp parallel for
	for(long ii =0; ii < (long)peaks.data.size(); ++ii){
		int tid = omp_get_thread_num();
#endif
#ifdef DEBUG
	for(long ii =0; ii < (long)peaks.data.size(); ++ii){
		int tid = 0;
#endif
		//print out progress
		omp_set_lock(&writelock);
		++progress;
		if(progress*100/peaks.data.size() > percent){
			percent = progress*100/peaks.data.size();
			fprintf(stderr, "*");
		}
		//The artifact test also requires write lock due to the FFTHandler
		peaks.data[ii].artifactScore = Artifact::test(parms.getMaxFragmentLength(),
				readLength, fir, peaks.data[ii],artifactTestRatio, halfLength);
		omp_unset_lock(&writelock);
		if(peaks.data[ii].artifactScore<0.5)
			peaks.data[ii].isArtifact = true;
		else
			peaks.data[ii].isArtifact = false;
		if(!peaks.data[ii].isArtifact){
			//do likelihood ratio test
			int ret =likelihoodRatioTest[tid]->test(peaks.data[ii].pstrand,peaks.data[ii].mstrand,
					peaks.data[ii].pstrandExtended,peaks.data[ii].mstrandExtended,
					peaks.data[ii].localPeakPos.size(),&peaks.data[ii].localPeakPos[0]);
			if(ret!=0){
				omp_set_lock(&writelock);
				++numNonConvergentTests;
				omp_unset_lock(&writelock);
			}
			peaks.data[ii].lpValue = likelihoodRatioTest[tid]->lpValue;
			peaks.data[ii].lpValueOpenChromatin = likelihoodRatioTest[tid]->lpValueOpenChromatin;
			peaks.data[ii].betaAlt = likelihoodRatioTest[tid]->peakReads;
		}
	}
	fprintf(stderr, "\n");
	//delete all test objects
	for(int ii = 0; ii < parms.getNumThreads(); ++ii){
		delete likelihoodRatioTest[ii];
	}
	Artifact::destroy();

	//count artifacts and candidate peaks
	long numArtifacts = 0;
	long numCandidatePeaks = 0;
	for(long ii =0; ii < (long)peaks.data.size(); ++ii)
		if(peaks.data[ii].isArtifact)
			++numArtifacts;
		else
			++numCandidatePeaks;
	fprintf(stderr, "Non-convergent likelihood tests[%d]\n", numNonConvergentTests);
	fprintf(stderr, "[%ld] read length Artifacts found\n",numArtifacts);
	fprintf(stderr, "Finished analyzing all Chromosomes found [%d] candidate peaks\n",(int)numCandidatePeaks);
	FDRcorrect();
}

void Ritornello::FDRcorrect(){

	fprintf(stderr, "Sorting by log ratio test p-value\n");
	std::sort (peaks.data.begin(), peaks.data.end(), comparepValue);

	//Benjamini hochberg correction
	fprintf(stderr, "Applying Benjamini Hochberg\n");
	int rank = 1;
	for(long ii =0; ii < (long)peaks.data.size();++ii){
		if(peaks.data[ii].isArtifact)
			continue;
			//peaks.data[ii].lqValue = peaks.data[ii].lpValue - log10((long)peaks.data.size()/rank);
		peaks.data[ii].lqValue = peaks.data[ii].lpValue - log10((double)peaks.positionsScanned/(double)rank);
		++rank;
	}


	//Yekutieli and Benjamini (1999) correction
	for(long ii =0; ii < (long)peaks.data.size();++ii){
		if(peaks.data[ii].isArtifact)
			continue;

		double bestlQVal = 0;
		for(long jj =ii; jj < (long)peaks.data.size();++jj){
			if(peaks.data[jj].isArtifact)
				continue;

			if(peaks.data[jj].lqValue > bestlQVal)
				bestlQVal = peaks.data[jj].lqValue;

		}
		peaks.data[ii].lqValue = bestlQVal;
	}

	long numPeaks =0;
	for(long ii =0; ii < (long)peaks.data.size();++ii){
		if(peaks.data[ii].lqValue >= parms.getLogQValSignifThreshold()
			&& peaks.data[ii].betaAlt >=parms.getSignalValueThreshold()
			&& !peaks.data[ii].isArtifact)
			++numPeaks;
	}

	fprintf(stderr, "[%ld] Peaks are significant after multiple hypothesis correction\n",numPeaks);

}

bool predicate(Peak peak){
	return peak.isArtifact;
}
void Ritornello::removeArtifacts(){
	for(vector<Peak>::iterator iter = peaks.data.begin(); iter != peaks.data.end(); ++iter)
		if(iter->isArtifact){
			delete iter->pstrandExtended;
			delete iter->mstrandExtended;
		}
	peaks.data.erase(std::remove_if(peaks.data.begin(), peaks.data.end(), predicate), peaks.data.end());
}
void Ritornello::writeResults(){
	string filename = (string(parms.getOutputPrefix())+"-peakSummary.narrowPeak");
	fprintf(stderr, "Writing results to [%s]\n",filename.c_str());

	//create a buffered reader to get the chromosome names
	BufferedDepthGraphReader bgr(2*parms.getMaxFragmentLength());
	bgr.init();
	//open the output file

	ofstream outputFile(filename.c_str());
	//for each peak
	for(long ii = 0; ii < (long)peaks.data.size(); ++ii){
		//if it is significant and not an artifact
		if(peaks.data[ii].lqValue>=parms.getLogQValSignifThreshold()
			&& peaks.data[ii].betaAlt >=parms.getSignalValueThreshold()
				&& !peaks.data[ii].isArtifact)
			//print it
			outputFile<<bgr.chromosomeNames[peaks.data[ii].chr]<<"\t"<<
				peaks.data[ii].pos<<"\t"<<
				peaks.data[ii].pos+1<<"\t"<<
				"."<<"\t"<<
				peaks.data[ii].betaAlt<<"\t"<<
				"."<<"\t"<<
				peaks.data[ii].betaAlt<<"\t"<<
				peaks.data[ii].lpValue<<"\t"<<
				peaks.data[ii].lqValue<<"\t"<<
				"0"<<
				//"\t"<<peaks.data[ii].lpValueOpenChromatin<<"\t"<<
				//peaks.data[ii].isArtifact<<"\t"<<
				endl;
	}
	outputFile.close();
	//open the output file
	filename = (string(parms.getOutputPrefix())+"-artifactSummary.bed");
	fprintf(stderr, "Writing artifacts to [%s]\n",filename.c_str());
	ofstream artifactOutputFile(filename.c_str());
	//for each peak
	for(long ii = 0; ii < (long)peaks.data.size(); ++ii){
		//if it is an artifact
		if(peaks.data[ii].isArtifact)
			//print it
			artifactOutputFile<<bgr.chromosomeNames[peaks.data[ii].chr]<<"\t"<<
				peaks.data[ii].pos<<"\t"<<
				peaks.data[ii].pos+1<<"\t"<<
				peaks.data[ii].artifactScore<<"\t"<<
				endl;
	}
	artifactOutputFile.close();
}

void Ritornello::calcWindowReadCountDist(){
#ifdef DEBUG
	medianWindowReadCount=10;
	fprintf(stderr, "Set median non-zero reads per window=[%ld]\n",medianWindowReadCount);
	return;
#endif
	fprintf(stderr, "Calculating median non-zero reads per window\n");
	long windowSize = parms.getMaxFragmentLength();
	BufferedDepthGraphReader* bgr;
	bgr = new BufferedDepthGraphReader(2*windowSize);
	bgr->init();
	while(bgr->next()){
		//calcualate mean coverage
		double numReads = bgr->getPstrandReads() + bgr->getMstrandReads();
		if(numReads >= windowReadCountDist.size()){
			long oldSize = windowReadCountDist.size();
			windowReadCountDist.resize(numReads+1);
			for(long ii = oldSize; ii < (long)windowReadCountDist.size();++ii){
				windowReadCountDist[ii]=0;
			}
		}
		++windowReadCountDist[numReads];
	}
	bgr->close();
	delete bgr;
	IOhandler::printDoubleArrayToFile(&windowReadCountDist[0],windowReadCountDist.size(),"windowReadCountDist.txt");
	//get 50 percentile for non zero read coverage
	double sum = 0;
	for(long ii = 1; ii < (long)windowReadCountDist.size();++ii){
		sum+=windowReadCountDist[ii];
	}
	double cumsum =0;
	long ii;
	for(ii = 1; ii < (long)windowReadCountDist.size();++ii){
		cumsum+=windowReadCountDist[ii];
		if(2*cumsum >= sum)
			break;
	}
	medianWindowReadCount = ii;
	fprintf(stderr, "Calculated median non-zero reads per window=[%ld]\n",medianWindowReadCount);

}


