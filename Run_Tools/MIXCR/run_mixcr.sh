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



export fastq=`sed -n "$SLURM_ARRAY_TASK_ID"p files.txt`


files="/path/files/input"
output="/path/files/output"

start=$SECONDS

##Paired end
mixcr analyze shotgun -s hsa --starting-material rna --receptor-type bcr --contig-assembly --impute-germline-on-export --assemble '-OseparateByV=false' --export '-nFeatureImputed VDJRegion -aaFeatureImputed VDJRegion' $files/$fastq.R1.25bp.length.250K.cov.fastq.gz $files/$fastq.R2.25bp.length.250K.cov.fastq.gz $fastq; mixcr exportClones $fastq.contigs.clns $fastq.contigs.full.txt; mv $fastq* $output/
