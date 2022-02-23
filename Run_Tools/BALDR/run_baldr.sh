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

export fastq=`sed -n "$SGE_TASK_ID"p files.txt`

module load star/2.5.2b                            # load star software into your job environment
module load baldr/1.0
module load java/1.8.0u102
module load samtools/1.3.1
module load bowtie2/2.3.0



files="/site/ne/home/i0439277/Datasets_Benchmark/BCR_Baldr_dataset_SRP126429/aw1_human_plasbmablast/50bp_length/50k_coverage/"
output="/site/ne/home/i0439277/baldr/BALDER2/BALDR-master/baldr_dataset/paired_end/50bp/cov50k"

start=$SECONDS

./BALDR --paired $files/$fastq.R1.50bp.length.50K.cov.fastq.gz,$files/$fastq.R2.50bp.length.50K.cov.fastq.gz --trinity /site/ne/home/i0439277/anaconda3/bin/Trinity --adapter /site/ne/app/x86_64/trimmomatic/v0.32/adapters/NexteraPE-PE.fa --trimmomatic /site/ne/app/x86_64/trimmomatic/v0.32/trimmomatic-0.32.jar --igblastn /site/ne/app/x86_64/ncbi/igblast/v1.6.1/bin/igblastn --STAR /site/ne/app/x86_64/star/v2.5.2b/STAR --STAR_index /site/ne/home/i0439277/baldr/BALDER2/BALDR-master/STAR_GRCh38/STAR_GRCh38_index --BALDR /site/ne/home/i0439277/baldr/BALDER2/BALDR-master/ --memory 40G --threads 6

end=$SECONDS

echo "duration: $((end-start)) seconds."
