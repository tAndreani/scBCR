start=$SECONDS
a=`cat files.txt`
echo "lets start by parsing the files!"
mkdir test;

for i in $a;
do
	cat $i/$i.fasta | grep light_chain -A 1 >  test/$i.light_chain.fasta;
	cat $i/$i.fasta | grep light_chain -A 1 | grep constant_region_contig -v |grep "-" -v | grep -A 1 light_chain | head -2 >  test/$i.light_chain.1.fasta;
	cat $i/$i.fasta | grep light_chain -A 1 | grep constant_region_contig -v |grep "-" -v | grep -A 1 light_chain | tail -2 >  test/$i.light_chain.2.fasta;
done

for i in $a;
do
        cat $i/$i.fasta | grep heavy_chain -A 1 >  test/$i.heavy_chain.fasta;
	cat $i/$i.fasta | grep heavy_chain -A 1 | grep constant_region_contig -v |grep "-" -v | grep -A 1 heavy_chain | head -2 >  test/$i.heavy_chain.1.fasta;
	cat $i/$i.fasta | grep heavy_chain -A 1 | grep constant_region_contig -v |grep "-" -v | grep -A 1 heavy_chain | tail -2 >  test/$i.heavy_chain.2.fasta;
done

#Run blast
echo "run blast!"
cd test/;
for i in $a;
do
        AssignGenes.py igblast -s $i.heavy_chain.fasta -b /usr/local/share/igblast --loci ig --organism human;
	AssignGenes.py igblast -s $i.heavy_chain.1.fasta -b /usr/local/share/igblast --loci ig --organism human;
        AssignGenes.py igblast -s $i.heavy_chain.2.fasta -b /usr/local/share/igblast --loci ig --organism human;

done

for i in $a;
do
        AssignGenes.py igblast -s $i.light_chain.fasta -b /usr/local/share/igblast --loci ig --organism human;
        AssignGenes.py igblast -s $i.light_chain.1.fasta -b /usr/local/share/igblast --loci ig --organism human;
        AssignGenes.py igblast -s $i.light_chain.2.fasta -b /usr/local/share/igblast --loci ig --organism human;

done

echo "delete empy files"
find -type f -empty -delete

echo "make database"
#Create database
for i in $a;
do 
        MakeDb.py igblast -i ${i}.light_chain_igblast.fmt7 -s $i.light_chain.fasta -r /usr/local/share/germlines/imgt/human/vdj --log $i.basic_IGKL.log --extended;
        MakeDb.py igblast -i ${i}.heavy_chain_igblast.fmt7 -s $i.heavy_chain.fasta -r /usr/local/share/germlines/imgt/human/vdj --log $i.basic_IGH.log --extended
	MakeDb.py igblast -i ${i}.light_chain.1_igblast.fmt7 -s $i.light_chain.1.fasta -r /usr/local/share/germlines/imgt/human/vdj --log $i.basic_IGKL.1.log --extended;
	MakeDb.py igblast -i ${i}.light_chain.2_igblast.fmt7 -s $i.light_chain.2.fasta -r /usr/local/share/germlines/imgt/human/vdj --log $i.basic_IGKL.2.log --extended;
        MakeDb.py igblast -i ${i}.heavy_chain.1_igblast.fmt7 -s $i.heavy_chain.1.fasta -r /usr/local/share/germlines/imgt/human/vdj --log $i.basic_IGH.1.log --extended;
        MakeDb.py igblast -i ${i}.heavy_chain.2_igblast.fmt7 -s $i.heavy_chain.2.fasta -r /usr/local/share/germlines/imgt/human/vdj --log $i.basic_IGH.2.log --extended;
done

echo "order the folders"
mkdir log_files;
mkdir tsv_files;
mkdir fmt7_files;
mkdir fasta_files;

mv -f *fasta fasta_files/;
mv -f *log log_files/;
mv -f *fmt7 fmt7_files/;
mv -f *tsv tsv_files/

end=$SECONDS
echo "duration: $((end-start)) seconds."
echo "your job is done"
