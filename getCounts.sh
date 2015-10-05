#!/bin/bash
 
#SBATCH -n 24                 # Number of cores
#SBATCH -t 6-12:05           # Runtime in D-HH:MM
#SBATCH -p serial_requeue    # Partition to submit to
#SBATCH --mem-per-cpu=4000            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o 1000GP_associations.out      # File to which STDOUT will be written
#SBATCH -e 1000GP_associations.err      # File to which STDERR will be written
#SBATCH --mail-type=END      # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dschlauch@fas.harvard.edu  # Email to which notifications will be sent

minAC=20
maxAC=50

mkdir counts_${minAC}_${maxAC}

python getMACforSamples.py 1 $minAC $maxAC &
python getMACforSamples.py 2 $minAC $maxAC &
python getMACforSamples.py 3 $minAC $maxAC &
python getMACforSamples.py 4 $minAC $maxAC &
python getMACforSamples.py 5 $minAC $maxAC &
python getMACforSamples.py 6 $minAC $maxAC &
python getMACforSamples.py 7 $minAC $maxAC &
python getMACforSamples.py 8 $minAC $maxAC &
python getMACforSamples.py 9 $minAC $maxAC &
python getMACforSamples.py 10 $minAC $maxAC &
python getMACforSamples.py 11 $minAC $maxAC &
python getMACforSamples.py 12 $minAC $maxAC &
python getMACforSamples.py 13 $minAC $maxAC &
python getMACforSamples.py 14 $minAC $maxAC &
python getMACforSamples.py 15 $minAC $maxAC &
python getMACforSamples.py 16 $minAC $maxAC &
python getMACforSamples.py 17 $minAC $maxAC &
python getMACforSamples.py 18 $minAC $maxAC &
python getMACforSamples.py 19 $minAC $maxAC &
python getMACforSamples.py 20 $minAC $maxAC &
python getMACforSamples.py 21 $minAC $maxAC &
python getMACforSamples.py 22 $minAC $maxAC

