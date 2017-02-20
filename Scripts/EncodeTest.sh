#BSUB -J Ritornello
#BSUB -n 8
#BSUB -q shared
#BSUB -W 24:00

ii=$LSB_JOBINDEX
Ritornello/bin/Ritornello -f RitornelloTest1/ENCFF000YMP.sorted.bam --debug-folder RitornelloTest$ii/ -o RitornelloTest$ii/ENCFF000YMP.sorted.bam --Correct-PCR -p 8 2> RitornelloTest$ii/StabilityLog.txt
