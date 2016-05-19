/*
 * Config.h
 *
 *  Created on: Dec 18, 2014
 *      Author: Kelly Stanton
 */
#ifndef CONFIG_H_
#define CONFIG_H_
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
using namespace std;
class Config {
public:
	Config();
	virtual ~Config();
	int parseCmdLine(int argc, char* argv[]);
	//Accessors
	char* getBamFileName();
	char* getOutputPrefix();
	double getPeakThreshold();
	long getMaxFragmentLength();
	char* getDebugFolder();
	char* getFilterFile();
	char* getFLDFile();
	char* getTRAINFile();
	bool getModelOpenChromatinEffect();
	bool getCorrectPCR();
	int getNumThreads();
	double getMinReadThreshold();
	double getSignalValueThreshold();
	double getLogQValSignifThreshold();
	void setBamFileName(char* argBamFileName);
	void setPeakThreshold(double argPeakThreshold);
	void setMaxFragmentLength(long argMaxFragmentLength);
	void setDebugFolder(char* argDebugFolder);
	void setFilterFile(char* argFilterFile);
	void setFLDFile(char* argFLDFile);
	void setTRAINFile(char* argTRAINFile);
	void setNumThreads(int argNumThreads);
private:
	//internal function to get the value of a command line option
	char* _getOption(char ** begin, char ** end, const std::string & option);
	//internal function to get if a command line option is present
	bool _optionExists(char** begin, char** end, const std::string& option);
	//The filename for ChIP-seq bam file.  This is the main input to the program
	char* bamFileName;
	//The ouput file for ChIP-seq narrowPeak file.  This is the main output of the program
	char* outputPrefix;
	//Cutoff used for identifying local maxima.
	double peakThreshold;
	//The maximum fragment length.  This gets updated as the program runs.
	long maxFragmentLength;
	//Location to print out debug files
	char* debugFolder;
	//File containing the filter shape if one is provided
	char* filterFile;
	//File containing the fragment length distribution
	char* FLDFile;
	//File containing the training peak information
	char* TRAINFile;
	//number of threads to use for parallelization
	int numThreads;
	//log q-value threshold below which peaks are no longer reported
	double logQValSignifThreshold;
	//minimum read count threshold
	double minReadThreshold;
	//signal value threshold
	double signalValueThreshold;
	//model open chromatin
	bool modelOpenChromatinEffect;
	//flag to turn on PCR correction
	bool correctPCR;
};

#endif /* CONFIG_H_ */
