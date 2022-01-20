# This bash script allows to take the input from BRACER folder with the reconstructed BCR filtered for each sample and annotate the H and L chain (Kappa or Lambda)
#  files.txt is a file with the identifier for each sample of each real dataset

start=$SECONDS
a=`cat files.txt`
mkdir test;
for i in $a;
        do cat ${i}_BCRseqs.fa | grep -A 1 'BRACER|BCR|K'  >  test/$i.light_chain_K.fasta;
done
for i in $a;
        do cat ${i}_BCRseqs.fa | grep -A 1 'BRACER|BCR|L'  >  test/$i.light_chain_L.fasta;
done
for i in $a;
        do cat ${i}_BCRseqs.fa | grep -A 1 'H' >  test/$i.heavy_chain.fasta;
done
#Run blast
for i in $a; 
do
        AssignGenes.py igblast -s test/$i.heavy_chain.fasta -b /usr/local/share/igblast --loci ig --organism human;
        AssignGenes.py igblast -s test/$i.light_chain_L.fasta -b /usr/local/share/igblast --loci ig --organism human;
        AssignGenes.py igblast -s test/$i.light_chain_K.fasta -b /usr/local/share/igblast --loci ig --organism human;
done
mkdir test/tsv_files test/fmt7_files test/fasta_files test/log_files;
echo "create the database"
#Create database
for i in $a;
do
        MakeDb.py igblast -i test/${i}.heavy_chain_igblast.fmt7 -s test/$i.heavy_chain.fasta -r /usr/local/share/germlines/imgt/human/vdj --log test/$i.IGH.log --extended;
        MakeDb.py igblast -i test/${i}.light_chain_L_igblast.fmt7 -s test/$i.light_chain_L.fasta -r /usr/local/share/germlines/imgt/human/vdj --log test/$i.IGL.log --extended;
        MakeDb.py igblast -i test/${i}.light_chain_K_igblast.fmt7 -s test/$i.light_chain_K.fasta -r /usr/local/share/germlines/imgt/human/vdj --log test/$i.IGK.log --extended;
done
mv test/*fasta test/fasta_files/;
mv test/*tsv test/tsv_files/;
mv test/*fmt7 test/fmt7_files/;
end=$SECONDS
echo "duration: $((end-start)) seconds."
echo "your job is done"
