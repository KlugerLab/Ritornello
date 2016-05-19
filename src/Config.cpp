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
	const char* helpInfo =
	"Ritornello options:\n\n"
	"	--help	print this message\n\n"
	"	--version	print Ritornello's current version\n\n"
	"	-f <ChIP.bam>	ChIP.bam is a ChIP-seq sorted bam file to call peaks on. "
	"If you're bamfile is not sorted, please use the samtools sort utility.  Additionally, "
	"running the samtools index utility "
	"may be required by some visualization tools (IGV etc.) but is not required by Ritornello.\n\n"
	"	-o <OutputPrefix>	Specifies a prefix to use when reporting output.  This can be a file path. "
	"Ex. -o /home/MyUser/MyOutputPrefix would report called peaks to /home/MyUser/MyOutputPrefix-peakSummary.narrowPeak \n\n"
	"	-p <int>	maximum number of threads to use during execution\n\n"

	"	-q <decimal>	The -log10(q-value) used to threshold reported peaks.  Ex.  Specifying -q 2 (the default) will report all peaks "
	"that are more significant than q-value=0.01.  Set -q 0 when calling peaks for input to the "
	"Irreproducible Discovery Rate software.\n\n"

	"	-s <decimal>	The signal value used to threshold reported peaks.  The signal value is the effect size for the "
	"reported peak and has units of read count.  It is the sum of beta1 for the positive and negative strands around a reported peak,"
	"and can best be interpreted as the maximum likelihood number of reads due to binding at the reported peak.  Ex. specifying"
	"-s 40 (the default) will report all peaks with signal value greater than 40.  Set -s 0 when calling peaks for input to the "
	"Irreproducible Discovery Rate software.\n\n"

	"	-n <decimal>	The minimum read and matched filter threshold.  This specifies the minimum number of reads per window "
	"of size twice the maximum fragment length centered around the peak required to perform a likelihood ratio test.  Additionally "
	"the matched filter (which is also normalized to units of read counts is thresholded by this number. "
	"This threshold is mostly used to control compute time by limiting the number of tests performed.  -n 20 is the default.  Setting "
	"it too high (larger than effect size) may cause lower expressed peaks to be missed.  Setting it lower generally increases runtime "
	"and memory usage. Set -n 10 when calling peaks for input to the Irreproducible Discovery Rate software.\n\n"

	"	--OCE	Specifying the OCE option tells Ritornello to include an additional term in the likelihood ratio test "
	"to control for Open Chromatin Effects.  These are areas of high coverage which are generally uniform and also present "
	"in sonicated input DNA.  Ritornello can call spurious peaks at the bounderies of these regions where coverage changes "
	"abruptly.  --OCE is useful to avoid spurious results in highly sequenced small genomes (yeast), but may cause a loss of "
	"sensitivity and not recommended for mouse, human, etc.\n\n"

	"	--Correct-PCR	Specifying the --Correct-PCR option tells Ritornello to preprocess the read coverage and control outliers "
	"likely due to PCR amplification bias.  We recommend using this option if Ritornello calls many spikey false positives\n\n"

	"Ritornello advanced options:\n"
	"	--debug-folder\n"
	"	--filter-file\n"
	"	--FLD-file\n"
	;
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
	/** deprecated
	TRAINFile = _getOption(argv, argv + argc, "--TRAIN-file");
	if(TRAINFile!=NULL){
		fprintf(stderr, "Reading training peaks from from file [%s]\n", TRAINFile);
	}
	*/
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
