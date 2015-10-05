#!/bin/bash
 
#SBATCH -n 24                 # Number of cores
#SBATCH -N 1                 # Ensure that all cores are on one machine
#SBATCH -t 6-12:05           # Runtime in D-HH:MM
#SBATCH -p serial_requeue    # Partition to submit to
#SBATCH --mem-per-cpu=4000            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o 1000GP_jaccard.out      # File to which STDOUT will be written
#SBATCH -e 1000GP_jaccard.err      # File to which STDERR will be written
#SBATCH --mail-type=END      # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dschlauch@fas.harvard.edu  # Email to which notifications will be sent

module load centos6/python-2.7.3
module load centos6/numpy-1.7.1_python-2.7.3
module load centos6/R-3.0.2

minAC=0
maxAC=5000

mkdir output_${minAC}_${maxAC}

python jaccard_chr_gz.py 1 $minAC $maxAC &
python jaccard_chr_gz.py 2 $minAC $maxAC &
python jaccard_chr_gz.py 3 $minAC $maxAC &
python jaccard_chr_gz.py 4 $minAC $maxAC &
python jaccard_chr_gz.py 5 $minAC $maxAC &
python jaccard_chr_gz.py 6 $minAC $maxAC &
python jaccard_chr_gz.py 7 $minAC $maxAC &
python jaccard_chr_gz.py 8 $minAC $maxAC &
python jaccard_chr_gz.py 9 $minAC $maxAC &
python jaccard_chr_gz.py 10 $minAC $maxAC &
python jaccard_chr_gz.py 11 $minAC $maxAC &
python jaccard_chr_gz.py 12 $minAC $maxAC &
python jaccard_chr_gz.py 13 $minAC $maxAC &
python jaccard_chr_gz.py 14 $minAC $maxAC &
python jaccard_chr_gz.py 15 $minAC $maxAC &
python jaccard_chr_gz.py 16 $minAC $maxAC &
python jaccard_chr_gz.py 17 $minAC $maxAC &
python jaccard_chr_gz.py 18 $minAC $maxAC &
python jaccard_chr_gz.py 19 $minAC $maxAC &
python jaccard_chr_gz.py 20 $minAC $maxAC &
python jaccard_chr_gz.py 21 $minAC $maxAC &
python jaccard_chr_gz.py 22 $minAC $maxAC

wait

Rscript jaccard_processing.R $minAC $maxAC

