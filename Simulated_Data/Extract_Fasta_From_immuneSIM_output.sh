#### Author: Tommaso Andreani
#### What this script does: It takes in input the immuneSIM output for a particular number of chains simulated (Heavy or Light - in our case 100) and output a fasta file for each of the sequences simulated.
#### Application: each fasta file will be used as a reference for ART tool which will create Illumina Libraries for each of the sample. Art script can be found in this folder.


#Introduce the ">" typical of the fasta format
sed 's/^/>/'  simulated_sequences_with_SHM_15shm.tsv > simulated_sequences_with_SHM_h_chain_15shm_.tsv ;

#take the first two columns: the first will be the Id and the second the DNA sequence
cat simulated_sequences_with_SHM_h_chain_15shm_.tsv | cut -f 1,2 > shm15_h.txt;

#The each DNA sequence in the second column must go below the Id of the first column: the FASTA file is not created for all the 100 sequences
cat shm15_h.txt | tr '\t' '\n' | grep sequence -v > shm15_h.fasta;

#Annotation using Immcantation
AssignGenes.py igblast -s shm15_h.fasta -b /usr/local/share/igblast --loci ig --organism human;
MakeDb.py igblast -i shm15_h_igblast.fmt7 -s shm15_h.fasta -r /usr/local/share/germlines/imgt/human/vdj --log shm15_h.log --extended;

#Now it is the time to split each sequence in 1 fasta file: a folder is created and a folder name for each fasta file is created inside this folder

mkdir split_fasta
cp shm15_h.fasta split_fasta/
split -l 2 shm15_h.fasta

# each file is renamed
a=`ls x*`
for i in $a; do mv $i $i.fasta; done

#this is the file name with all the names of each fasta sequence
ls x*fasta | awk -F "." '{print $1}' > list_files.txt
