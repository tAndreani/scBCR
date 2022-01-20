#!/bin/sh
#SBATCH --job-name=BCR_Reconstruction      # job name
#SBATCH --array=1-48                       # Number of Job (according to the number of splitted files in the exported variables)
#SBATCH --nodes=2                          # nodes
#SBATCH -c 2                               # cores
#SBATCH --mem=6G                           # memory
#SBATCH --time=72:00:00                    # time
#SBATCH --error=BCR_TCR.err                # error file name
#SBATCH --output=BCR_TCR.out               # output file name
#SBATCH --mail-user=Tommaso.Andreani@sanofi.com  # email
#SBATCH --mail-type=ALL                          # type notification

files="/path/input_fastq/"
output="/path/output/"

mkdir $output/$fastq;
mkdir $output/$fastq.tmp;

~/anaconda3/envs/basic/bin/BASIC.py -b ~/anaconda3/envs/basic/bin/ -PE_1 $files/${fastq}_1.fastq.gz -PE_2 $files/${fastq}_2.fastq.gz -g human -i BCR -o $output/$fastq -t $output/$fastq.tmp -n $fastq
