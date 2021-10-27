# second create the fastq files: Input the fasta sequence for the sintetic BCR above with 15 mutations in the variable domain. 
#in sabre: path to "Art" the tool to simulate fastq reads: art_illumina="/site/ne/app/x86_64/igsimulator/v1/src/art_bin_VanillaIceCream/./art_illumina"

# Parameters:
# 1) -p is paired end
# 2) -l length 75 
# 3) -f coverage 1 mln reads
# 4) -ir and -ir2 is error rate for R1 and R2 set at the minimum. 
# 5) -m mean size of DNA/RNA fragments for paired-end simulations
# 6) -s the standard deviation of DNA/RNA fragment size for paired-end simulations

a="list of identifiers of 100 fasta sequences simulated above"

for i in $a;
do 
    $art_illumina -i $i.sequence.fasta -p -l 75 -f 500000 -ir 0.00001 -ir2 0.00001 -m 150 -s 10 -o $i.immuneSIM_15_mutations_IGH; 
done
