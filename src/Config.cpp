/*
 * Config.cpp
 *
 *  Created on: Dec 18, 2014
 *      Author: Kelly Stanton
 */
#include "Config.h"
#include <omp.h>
Config::Config() {
	//input bam file to run ritornello on
	bamFileName = NULL;
	//Threshold for peak filter
	peakThreshold=0.5;
	//maximum fragment length.  Note that this will change as the algorithm runs
	maxFragmentLength=1024;
	//Folder to output debugging files
	debugFolder =NULL;
	// Get the number of processors in this system
	numThreads = omp_get_num_procs();
	//file which contains the filter shape is provided
	filterFile=NULL;
	//file with peaks and their coverage
	TRAINFile=NULL;
	//Output file
	outputPrefix = NULL;
	//log q-value threshold below which peaks are no longer reported
	logQValSignifThreshold = 2.0;
	//minimum read count threshold
	minReadThreshold = 20;
	//signal value threshold
	signalValueThreshold = 40;
	//model open chromatin
	modelOpenChromatinEffect = false;
	//Flag to turn on PCR correction
	correctPCR = false;
	//Flag to turn on PCR correction
	handleArtifacts = true;
}
Config::~Config() {
}
void printVersion(){
	const char* versionInfo =
	"Ritornello version 1.0.0"
	;
	fprintf(stderr,"%s\n",versionInfo);
}
void printHelp(){
	const char* helpInfo = "Please refer to the documentation in the README.md or at https://github.com/KlugerLab/Ritornello";
	fprintf(stderr,"%s\n",helpInfo);
}
/**
 * Parses the command line setting the global parameters of the program
 *
 * @param begin address of the first command line argument
 * @param end address of the last command line argument
 */
int Config::parseCmdLine(int argc, char* argv[]){
	if(_optionExists(argv, argv + argc, "--version")){
		printVersion();
		exit(0);
	}
	if(_optionExists(argv, argv + argc, "--help")){
		printHelp();
		exit(0);
	}

	bamFileName = _getOption(argv, argv + argc, "-f");
	if(bamFileName==NULL){
		 fprintf(stderr, "-f<ChIP.bam> is required.  in.bam is a sorted bam file\n");
		 exit(EXIT_FAILURE);
	}
	outputPrefix = _getOption(argv, argv + argc, "-o");
	if(outputPrefix ==NULL){
		outputPrefix = bamFileName;
	} else {
		fprintf(stderr, "Output prefix set to [%s]\n", outputPrefix);
	}
	/** deprecated
	char * MAX_FRAGMENT_LENGTH_STR = _getOption(argv, argv + argc, "-m");
	if(MAX_FRAGMENT_LENGTH_STR){
		setMaxFragmentLength(atoi(MAX_FRAGMENT_LENGTH_STR));
	}
	*/
	debugFolder = _getOption(argv, argv + argc, "--debug-folder");
	if(debugFolder!=NULL){
		fprintf(stderr, "Debug folder set to [%s]\n", debugFolder);
	}
	filterFile = _getOption(argv, argv + argc, "--filter-file");
	if(filterFile!=NULL){
		fprintf(stderr, "Reading filter from file [%s]\n", filterFile);
	}
	FLDFile = _getOption(argv, argv + argc, "--FLD-file");
	if(FLDFile!=NULL){
		fprintf(stderr, "Reading FLD from file [%s]\n", FLDFile);
	}
	TRAINFile = _getOption(argv, argv + argc, "--TRAIN-file");
	if(TRAINFile!=NULL){
		fprintf(stderr, "Reading training peaks from from file [%s]\n", TRAINFile);
	}
	char* NUM_THREADS_STR = _getOption(argv, argv + argc, "-p");
	if(NUM_THREADS_STR){
		numThreads = atoi(NUM_THREADS_STR);
		fprintf(stderr, "Set number to threads to [%d]\n", numThreads);
	}
	char* LQVAL_SIGNIF_STR = _getOption(argv, argv + argc, "-q");
	if(LQVAL_SIGNIF_STR){
		logQValSignifThreshold = atof(LQVAL_SIGNIF_STR);
		fprintf(stderr, "Setting the log q-Value significance threshold to [%f]\n", logQValSignifThreshold);
	}
	char* SIG_VAL_STR = _getOption(argv, argv + argc, "-s");
	if(SIG_VAL_STR){
		signalValueThreshold = atof(SIG_VAL_STR);
		fprintf(stderr, "Setting the signal value threshold to [%f]\n", signalValueThreshold);
	}
	char* MIN_READ_STR = _getOption(argv, argv + argc, "-n");
	if(MIN_READ_STR){
		minReadThreshold = atof(MIN_READ_STR);
		fprintf(stderr, "Setting the minimum read and matched filter threshold to [%f]\n", minReadThreshold);
	}
	if(_optionExists(argv, argv + argc, "--OCE")){
		modelOpenChromatinEffect=true;
		fprintf(stderr, "Adding open chromatin effect to the model\n");
	}
	if(_optionExists(argv, argv + argc, "--Correct-PCR")){
		correctPCR=true;
		fprintf(stderr, "PCR correction preprocessing on\n");
	}
	if(_optionExists(argv, argv + argc, "--no-artifact-handling")){
		handleArtifacts=false;
		fprintf(stderr, "Turning off artifact handling\n");
	}
	return 0;
}
//Accessors
char* Config::getBamFileName(){
	return bamFileName;
}
char* Config::getOutputPrefix(){
	return outputPrefix;
}
double Config::getPeakThreshold(){
	return peakThreshold;
}
long Config::getMaxFragmentLength(){
	return maxFragmentLength;
}
char* Config::getDebugFolder(){
	return debugFolder;
}
char* Config::getFilterFile(){
	return filterFile;
}
char* Config::getFLDFile(){
	return FLDFile;
}
char* Config::getTRAINFile(){
	return TRAINFile;
}
int Config::getNumThreads(){
	return numThreads;
}
double Config::getLogQValSignifThreshold(){
	return logQValSignifThreshold;
}
double Config::getSignalValueThreshold(){
	return signalValueThreshold;
}
double Config::getMinReadThreshold(){
	return minReadThreshold;
}
bool Config::getModelOpenChromatinEffect(){
	return modelOpenChromatinEffect;
}
bool Config::getCorrectPCR(){
	return correctPCR;
}
bool Config::getHandleArtifacts(){
	return handleArtifacts;
}
void Config::setBamFileName(char* argBamFileName){
	if(bamFileName!=NULL) delete bamFileName;
	bamFileName=argBamFileName;
}
void Config::setPeakThreshold(double argPeakThreshold){
	peakThreshold=argPeakThreshold;
}
void Config::setMaxFragmentLength(long argMaxFragmentLength){
	maxFragmentLength=argMaxFragmentLength;
	fprintf(stderr, "Max fragment length set to [%ld]\n", maxFragmentLength);
}
void Config::setDebugFolder(char* argDebugFolder){
	if(debugFolder!=NULL) delete debugFolder;
	debugFolder=argDebugFolder;
}
void Config::setFilterFile(char* argFilterFile){
	if(filterFile!=NULL) delete filterFile;
	filterFile=argFilterFile;
}
void Config::setFLDFile(char* argFLDFile){
	if(FLDFile!=NULL) delete FLDFile;
	FLDFile=argFLDFile;
}
void Config::setTRAINFile(char* argTRAINFile){
	if(TRAINFile!=NULL) delete TRAINFile;
	TRAINFile=argTRAINFile;
}
void Config::setNumThreads(int argNumThreads){
	numThreads=argNumThreads;
}
/**
 * Get the value of a command line option
 *
 * @param begin address of the first command line argument
 * @param end address of the last command line argument
 * @param option command line option for which we want the associated argument
 * @return string containing the value following the command line option
 */
char* Config::_getOption(char ** begin, char ** end, const std::string & option){
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end){
        return *itr;
    }
    return 0;
}
/**
 * Get whether a command line option is set
 *
 * @param begin address of the first command line argument
 * @param end address of the last command line argument
 * @param option the command line option to check whether it is present
 * @return bool that is true if the option is present and false otherwise
 */
bool Config::_optionExists(char** begin, char** end, const std::string& option){
    return std::find(begin, end, option) != end;
}
