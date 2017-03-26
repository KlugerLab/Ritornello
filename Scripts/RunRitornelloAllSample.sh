#!/bin/bash
#BSUB -J Ritornello
#BSUB -n 8
#BSUB -q shared
#BSUB -W 24:00

ii=$LSB_JOBINDEX
cd ~

FOLDER_ROOT="EncodeTest/RitornelloPeaks3-19-2017"
FileList=(
'ENCFF000YVE.sorted.bam'
'ENCFF000YTJ.sorted.bam'
'ENCFF000OPO_rep1.sorted.bam'
'ENCFF000OPR_rep2.sorted.bam'
'ENCFF000YTD.sorted.bam'
'ENCFF001NPX.sorted.bam'
'ENCFF001NPY.sorted.bam'
'ENCFF001NSH.sorted.bam'
'ENCFF001NSI.sorted.bam'
'ENCFF000YVG.sorted.bam'
'ENCFF000YMP.sorted.bam'
'ENCFF000YMM.sorted.bam'
'ENCFF000OWI.sorted.bam'
'ENCFF000VUG_rep1.sorted.bam'
'ENCFF000BNJ.sorted.bam'
'ENCFF000VUH_rep2.sorted.bam'
'ENCFF000OQE.sorted.bam'
'ENCFF000QCM_rep2.sorted.bam'
'ENCFF000OFD.sorted.bam'
'ENCFF000OUM.sorted.bam'
'ENCFF000OUH.sorted.bam'
'ENCFF000OQH.sorted.bam'
'ENCFF000BNL.sorted.bam'
'ENCFF000QCL_rep1.sorted.bam'
'ENCFF000OWG.sorted.bam'
'ENCFF000OFG.sorted.bam'
'ENCFF000OME_rep1.sorted.bam'
'ENCFF000OMH_rep2.sorted.bam'
'SRR1009108.sorted.bam'
'SRR1009109.sorted.bam'
'SRR1009111.sorted.bam'
)
fileName=${FileList[ii]}
rm -R -f $FOLDER_ROOT/${fileName}_Debug
mkdir $FOLDER_ROOT/${fileName}_Debug
Ritornello/bin/Ritornello -f EncodeTest/Data/${fileName} --debug-folder $FOLDER_ROOT/${fileName}_Debug/ 2> $FOLDER_ROOT/${fileName}_Debug/log.txt -o $FOLDER_ROOT/${fileName}

