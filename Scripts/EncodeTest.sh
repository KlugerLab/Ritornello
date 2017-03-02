#BSUB -J Ritornello
#BSUB -n 8
#BSUB -q shared
#BSUB -W 24:00

ii=$LSB_JOBINDEX
cd ~
Ritornello/bin/Ritornello -p 1 -f RitornelloTest1/ENCFF000YMP.sorted.bam --debug-folder RitornelloTest$ii/ -o RitornelloTest$ii/ENCFF000YMP.sorted.bam -p 1 2> RitornelloTest$ii/StabilityLog.txt
