
#!/bin/sh
#SBATCH --job-name=BCR_Reconstruction      # job name
#SBATCH --array=1-48                       # Number of Job (according to the number of splitted files in the exported variables)
#SBATCH --nodes=2                          # nodes
#SBATCH -c 2                               # cores
#SBATCH --mem=6G                       # memory
#SBATCH --time=72:00:00                    # time
#SBATCH --error=BCR_TCR.err                # error file name
#SBATCH --output=BCR_TCR.out               # output file name
#SBATCH --mail-user=Tommaso.Andreani@sanofi.com  # email
#SBATCH --mail-type=ALL                          # type notification

#####module load bracer

module load bracer/1.0


export fastq=`sed -n "$SGE_TASK_ID"p files.txt`

##paired end
files="/path_to_files/"

start=$SECONDS

python /path/bracer assemble -p 2 -c /path/bracer.conf -s Hsap $fastq $files/ $fastq.R1.gz $fastq.R2.gz
