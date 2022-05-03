#Run blast
AssignGenes.py igblast -s Heavy_Chains.fasta -b /usr/local/share/igblast --loci ig --organism human;
echo "create the database"
#Create database
MakeDb.py igblast -i Heavy_Chains_igblast.fmt7 -s Heavy_Chains.fasta -r /usr/local/share/germlines/imgt/human/vdj --log Heavy_Chains.log --extended;

#Run blast
AssignGenes.py igblast -s Light_Chains_Kappa.fasta -b /usr/local/share/igblast --loci ig --organism human;
echo "create the database"
#Create database
MakeDb.py igblast -i Light_Chains_Kappa_igblast.fmt7 -s Light_Chains_Kappa.fasta -r /usr/local/share/germlines/imgt/human/vdj --log Light_Chains_Kappa.log --extended;

#Run blast
AssignGenes.py igblast -s Light_Chains_Lambda.fasta -b /usr/local/share/igblast --loci ig --organism human;
echo "create the database"
#Create database
MakeDb.py igblast -i Light_Chains_Lambda_igblast.fmt7 -s Light_Chains_Lambda.fasta -r /usr/local/share/germlines/imgt/human/vdj --log Light_Chains_Lambda.log --extended;
