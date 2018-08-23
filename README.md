# Ritornello
Ritornello is a ChIP-seq peak calling algorithm based on signal processing that can accurately call binding events without the need to do a pair total DNA input or IgG control sample.  It has been tested for use with narrow binding events such as transcription factor ChIP-seq.

[Ritornello Preprint](http://biorxiv.org/content/early/2015/12/11/034090)

# Download precompiled binaries

Currently compiled for ubuntu x64.  Mac and Windows versions coming soon.

[Download](https://github.com/kstant0725/Ritornello/releases)

# Compiling on ubuntu:

install dependencies:

-samtools

`sudo apt-get install libbam-dev`

-FFTW

`sudo apt-get install libfftw3-dev`

-boost

`sudo apt-get install libboost-dev`

Checkout the source using:

`git clone https://github.com/KlugerLab/Ritornello.git`  

Compile:

`cd Ritornello`

`make`

The executable is then made at `bin/Ritornello`.
You can move it where ever you like and/or add it to your path, or simply run it from its bin location.

# Creating a sorted bam file:

This tutorial assumes the user starts with sequenced ChIP-seq reads in the fastq format, `MyFile.fastq` for singled end or `MyFile_1.fastq` and `MyFile_2.fastq` for paired-end.

The first step in preparing a sorted bam file, the input for Ritornello, is to map the ChIP-seq reads to the reference genome.  This can be done using any comparative genome alignment tool.  For this tutorial we will use Bowtie, which can be downloaded [here](http://bowtie-bio.sourceforge.net/index.shtml)

Also install samtools, which can be downloaded [here](http://samtools.sourceforge.net/) or on Ubuntu installed using:

`sudo apt-get install samtools`

Map the reads to the genome using the following command.  In this example hg19 (human genome version 19) is the prefix for the bowtie index (reference genome), so make sure to download the correct index (also found on the bowtie website for most organisms) and specify the prefix that is appropriate for your organism.  `MyFile.fastq` for singled end or `MyFile_1.fastq` and `MyFile_2.fastq` for paired-end are user provided fastq files containing the ChIP-seq reads from the sequencer.  For single end run:  

`bowtie -S -n 2 -k 1 -m 1 -X 1000 -I 0 --best --strata hg19 MyFile.fastq | samtools view -bS - > MyBamFile.bam`

or if you have paired end data, run:

`bowtie -S -n 2 -k 1 -m 1 -X 1000 -I 0 --best --strata -1 MyFile_1.fastq -2 MyFile_2.fastq hg19 | samtools view -bS - > MyBamFile.bam`
        
This will create a bam (alignment) file in the same directory.  See the bowtie documentation for a detailed explanation of each option.

To create a sorted bam file, call the samtools sort command on your bam file as follows:

`samtools sort MyBamFile.bam MySortedBamFile`

Indexing the bam file is useful so that it can be used with other tools such as the IGV genome browser, but not required to run ritornello.  To index the bam file, run the following:

`samtools index MySortedBamFile.bam`

This will create an index `MySortedBamFile.bam.bai` file in the same directory

# Using Ritornello:
-basic usage

`./Ritornello -f MySortedBamFile.bam `

Where `MySortedBamFile.bam` is an index/sorted bam file that can be obtained by first mapping the fastq files using an aligner (such as bowtie) and then sorting and indexing using samtools.

# Analyzing the output:

The output is printed to `MySortedBamFile.bam-peakSummary.narrowPeak` where `MySortedBamFile.bam` is replaced with the output prefix if the `-o` option is specified.  The output is a [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) format file, an extension of the bed format.

To analyze these results in R we recommend using the `rtracklayer` package.  A sample analysis for comparing overlapping TSS annotation is as follows:

	#import the rtracklayer library
	library(rtracklayer)
	
	#define the extra columns need for the narrowpeak format (in addition to the bed format 
	#columns)
	extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
		qValue = "numeric", peak = "integer")
		
	#read in the peaks called by Ritornello to a GRanges object
	RitornelloPeaks = import.bed(con="MySortedBamFile.bam-peakSummary.narrowPeak",
		extraCols = extraCols_narrowPeak)
	
	#read in the tss annotation
	TSS = import.bed(con="tss.bed")
	
	extend = function(x,n){
		start(x) = start(x)-n
		end(x) = end(x)+n
		return(x)
	}
	
	#extend the GRanges object by 100 bp upstream and downstream
	RitornelloPeaksWithin100bp =extend(RitornelloPeaks,100)
	
	#finally we can check for peaks that overlapped the TSS annotation
	RitornelloTSSPeaks =  RitornelloPeaks[countOverlaps(RitornelloPeaksWithin100bp, TSS)>0]

The full script with details on how to create all required files and run the script is provided
[here](https://github.com/KlugerLab/Ritornello/blob/master/Scripts/AnalyzeRitornelloOutput.R)
 
# Ritornello options:
`--help`	print the help message

`--version`	print Ritornello's current version

`-f <ChIP.bam>`	ChIP.bam is a ChIP-seq sorted bam file to call peaks on.  If you're bamfile is not sorted, please use the samtools sort utility.  Additionally, running the samtools index utility may be required by some visualization tools (IGV etc.) but is not required by Ritornello.

`-o <OutputPrefix>`	Specifies a prefix to use when reporting output.  This can be a file path.
Ex. `-o /home/MyUser/MyOutputPrefix` would report called peaks to `/home/MyUser/MyOutputPrefix-peakSummary.narrowPeak`

`-p <int>`	maximum number of threads to use during execution

`-q <decimal>`	The -log10(q-value) used to threshold reported peaks.  Ex.  Specifying `-q 2` (the default) will report all peaks that are more significant than q-value=0.01.  Set `-q 0` when calling peaks for input to the Irreproducible Discovery Rate software.

`-s <decimal>`	The signal value used to threshold reported peaks.  The signal value is the effect size for the reported peak and has units of read count.  It is the sum of beta1 for the positive and negative strands around a reported peak, and can best be interpreted as the maximum likelihood number of reads due to binding at the reported peak.  Ex. specifying `-s 40` (the default) will report all peaks with signal value greater than 40.  Set `-s 0` when calling peaks for input to the Irreproducible Discovery Rate software.

`-n <decimal>`	The minimum read and matched filter threshold.  This specifies the minimum number of reads per window of size twice the maximum fragment length centered around the peak required to perform a likelihood ratio test.  Additionally the matched filter (which is also normalized to units of read counts is thresholded by this number. This threshold is mostly used to control compute time by limiting the number of tests performed.  `-n 20` is the default.  Setting it too high (larger than effect size) may cause lower expressed peaks to be missed.  Setting it lower generally increases runtime and memory usage. Set `-n 10` when calling peaks for input to the Irreproducible Discovery Rate software.

`--OCE`	Specifying the OCE option tells Ritornello to include an additional term in the likelihood ratio test to more strictly control for Open Chromatin Effects.  These are areas of high coverage which are generally uniform and also present in sonicated input DNA.  Ritornello's background coverage term can usually control for most open chromatin effects, however, when coverage is extremely high, it can call spurious peaks at the bounderies of these regions where coverage changes abruptly.  `--OCE` is useful to avoid spurious results in highly sequenced small genomes (yeast), but may cause a loss of sensitivity and not recommended for mouse, human, etc.

`--no-artifact-handling`	Specifying the `--no-artifact-handling` option tells Ritornello not to try to detect and remove read length artifacts.  You can add this option if your reads are fairly long and you suspect there wont be any issues with artifacts related to mismapping of reads.

# Ritornello advanced options:

`--debug-folder`	Specify a folder to print out debug files (must end with a "/").  This is mainly used for development purposes.

`--filter-file`	Specify a filter shape to use for this run.  Filter shapes are printed to fir.txt in the debug folder when it is specified.  It is simply a vector of numbers giving the negative strand (or reverse positive strand) filter shape.  --FLD-file must be set to use this option.  Mainly used for development purposes

`--FLD-file`		Specify a fragment length distribution for use with this run.  FLD files are printed to fld.txt in the debug folder when the option is specified.  It is simply a vector giving the distibution of the fragment lengths.  Mainly used for development purposes