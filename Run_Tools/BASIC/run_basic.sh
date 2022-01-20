files="/site/ne/home/i0439277/TRUST4/Pilot1/mapping2"
output="/site/ne/home/i0439277/basic/benchmark/pilot1_TT_plus_Linda"

mkdir $output/$fastq;
mkdir $output/$fastq.tmp;

~/anaconda3/envs/basic/bin/BASIC.py -b ~/anaconda3/envs/basic/bin/ -PE_1 $files/${fastq}_1.fastq.gz -PE_2 $files/${fastq}_2.fastq.gz -g human -i BCR -o $output/$fastq -t $output/$fastq.tmp -n $fastq
