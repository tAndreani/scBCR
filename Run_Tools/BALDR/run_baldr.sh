
#!/bin/bash

#$ -S /bin/bash                                    # you can overwrite the shell to use csh or sh
#$ -N baldr_data_50bp_50k                               # name of your job
#$ -j y                                            # merge standard output and error into a file
#$ -m bea                                          # alert when the job begins, ends and abort
#$ -l mem_total=6G                                 # looking for a server with at least 48G RAM free
#$ -pe threaded 2                                  # asks for a server with at least 6 cores available and reserve 6 cores for your job to run
#$ -q c32.q                                        # asks for the job to run in the all.q queue
#$ -b y                                            # run binary on a node
#$ -cwd                                            # using current working directory
#$ -t 1-49:1

. /etc/profile.d/00-site.sh                        # load profile into your shell
. /etc/profile.d/modules.sh                        # load module profile into your shell in order to use module command as below

export fastq=`sed -n "$SGE_TASK_ID"p files_pe.txt`

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
