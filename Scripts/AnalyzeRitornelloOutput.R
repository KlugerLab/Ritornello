#Make sure to install rtracklayer using 
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")

#import the rtracklayer library
library(rtracklayer)

#define the extra columns need for the narrowpeak format (in addition to the bed format columns)
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")

#note that if your narrowPeak Ritornello output file or your annotation bed file (tss.bed) in this case
#is not in the current working directory you must specify the full path.
#alternatively you can use the setwd() command to set your working directory in R to where your narrowPeak file is
setwd("~/git/Ritornello/Scripts")

#read in the peaks called by Ritornello to a GRanges object
RitornelloPeaks = import.bed(con="MySortedBamFile.bam-peakSummary.narrowPeak",extraCols = extraCols_narrowPeak)

#read in any annotations you would like to compare against.  Here we compare against upstream of genes (TSS) in hg19 
#to create tss.bed
#1. go to https://genome.ucsc.edu/cgi-bin/hgTables
#2. set:
#		Clade: Mammal
#		Genome: Human
#		Group: Genes and Gene Predictions
#		Track: RefSeq Genes
#		Table: refGene
#		Output format: BED - Bed Extensible Data
#3. Click "get output"
#4. On the next screen select "upstream by" and set it to 1 bases
#5. Click "get BED"
#6. Right click on the screen of text and select "save as..."
#7. Rename the file to "tss.bed" and move it to the directory with this R script and your ritornello output

#read in the tss annotation
TSS = import.bed(con="tss.bed")

#Ritornello reports only the single most likely position where the transcription factor binding was detected in
#the ChIP-seq data.  In order to look at features close by, we will extend the genomic range from 1 bp to 100 bp
#upstream and downstream of the reported binding event.
extend = function(x,n){
	start(x) = start(x)-n
	end(x) = end(x)+n
	return(x)
}

#extend the GRanges object by 100 bp upstream and downstream
RitornelloPeaksWithin100bp =extend(RitornelloPeaks,100)

#finally we can check for peaks that overlapped the TSS annotation
RitornelloTSSPeaks =  RitornelloPeaks[countOverlaps(RitornelloPeaksWithin100bp, TSS)>0]

#print the peaks
RitornelloTSSPeaks

#if you prefer not to work with the GRanges object you can always convert it to a dataframe using
RitornelloTSSPeaksDataFrame = as.data.frame(RitornelloTSSPeaks)

#print the dataframe of peaks
RitornelloTSSPeaksDataFrame
