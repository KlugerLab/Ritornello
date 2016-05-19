# Ritornello
Ritornello is a ChIP-seq peak calling algorithm based on signal processing that can accurately call binding events without the need to do a pair total DNA input or IgG control sample.

[Ritornello Preprint](http://biorxiv.org/content/early/2015/12/11/034090)

# Download precompiled binaries

Currently compiled for ubuntu x64.  Mac and Windows versions coming soon.

[Download](https://github.com/kstant0725/Ritornello/releases)

# Compiling on ubuntu:

1.  install dependencies

-samtools

sudo apt-get install libbam-dev

-FFTW

sudo apt-get install fftw3-dev

-boost

sudo apt-get install libboost-dev

2.  cd to Ritornello directory and type make


#Usage:
-basic usage

./Ritornello -f MySortedBamFile.bam

Where MySortedBamFile.bam is an index/sorted bam file that can be obtained by first mapping the fastq files using an aligner (such as bowtie) and then sorting and indexing using samtools

-p [numProcs]  specify the maximum number of threads to use during Ritornello's execution

-o [outPrefix]  A string specifying the prefix to use when generating output.  This can be a file path.  Ex.   /home/MyUser/MyOutputPrefix

-q [-log10(q-value) cutoff] (default is 2).  For example, if you would like to threshold at q-value 0.01, specify -q 2

-s  [signal value threshold].  The minimum effect size to report a peak.  This is Beta1  (summed for both strands) from the publication.  Default is -s 40

#Ritornello options:
--help	print the help message

--version	print Ritornello's current version

-f <ChIP.bam>	ChIP.bam is a ChIP-seq sorted bam file to call peaks on.  If you're bamfile is not sorted, please use the samtools sort utility.  Additionally, running the samtools index utility may be required by some visualization tools (IGV etc.) but is not required by Ritornello.

-o <OutputPrefix>	Specifies a prefix to use when reporting output.  This can be a file path.
Ex. -o /home/MyUser/MyOutputPrefix would report called peaks to /home/MyUser/MyOutputPrefix-peakSummary.narrowPeak

-p <int>	maximum number of threads to use during execution

-q <decimal>	The -log10(q-value) used to threshold reported peaks.  Ex.  Specifying -q 2 (the default) will report all peaks that are more significant than q-value=0.01.  Set -q 0 when calling peaks for input to the Irreproducible Discovery Rate software.

-s <decimal>	The signal value used to threshold reported peaks.  The signal value is the effect size for the reported peak and has units of read count.  It is the sum of beta1 for the positive and negative strands around a reported peak, and can best be interpreted as the maximum likelihood number of reads due to binding at the reported peak.  Ex. specifying -s 40 (the default) will report all peaks with signal value greater than 40.  Set -s 0 when calling peaks for input to the Irreproducible Discovery Rate software.

-n <decimal>	The minimum read and matched filter threshold.  This specifies the minimum number of reads per window of size twice the maximum fragment length centered around the peak required to perform a likelihood ratio test.  Additionally the matched filter (which is also normalized to units of read counts is thresholded by this number. This threshold is mostly used to control compute time by limiting the number of tests performed.  -n 20 is the default.  Setting it too high (larger than effect size) may cause lower expressed peaks to be missed.  Setting it lower generally increases runtime and memory usage. Set -n 10 when calling peaks for input to the Irreproducible Discovery Rate software.

--OCE	Specifying the OCE option tells Ritornello to include an additional term in the likelihood ratio test to control for Open Chromatin Effects.  These are areas of high coverage which are generally uniform and also present in sonicated input DNA.  Ritornello can call spurious peaks at the bounderies of these regions where coverage changes abruptly.  --OCE is useful to avoid spurious results in highly sequenced small genomes (yeast), but may cause a loss of sensitivity and not recommended for mouse, human, etc.

--Correct-PCR	Specifying the --Correct-PCR option tells Ritornello to preprocess the read coverage and control outliers likely due to PCR amplification bias.  We recommend using this option if Ritornello calls many spikey false positives

#Ritornello advanced options:

--debug-folder	Specify a folder to print out debug files (must end with a "/")

--filter-file	Specify a filter shape to use for this run.  Filter shapes are printed to fir.txt in the debug folder when it is specified.  It is simply a vector of numbers giving the negative strand (or reverse positive strand) filter shape.  --FLD-file must be set to use this option

--FLD-file		Specify a fragment length distribution for use with this run.  FLD files are printed to fld.txt in the debug folder when the option is specified.  It is simply a vector giving the distibution of the fragment lengths.