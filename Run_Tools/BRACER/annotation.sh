#The script takes in input the chains annotated by immmcantation from the script "immcantation.sh" and extract gene names and productivity 

ls *tsv | awk -F "." '{print $1}' > tmp;
a=`cat tmp`;
mkdir IGH IGK IGL;
for i in $a; do echo $i > $i.name; done;
for i in $a; do cat $i.heavy_chain_igblast_db-pass.tsv | head -2 | cut -f 4,5,6,7 | grep productive -v > IGH/$i.info.hc ; done;
for i in $a; do cat $i.light_chain_K_igblast_db-pass.tsv | head -2 |cut -f 4,5,7 | grep productive -v > IGK/$i.info.lck ; done;
for i in $a; do cat $i.light_chain_L_igblast_db-pass.tsv | head -2 |cut -f 4,5,7 | grep productive -v > IGL/$i.info.lcl ; done;
for i in $a; do paste $i.name IGH/$i.info.hc > IGH/$i.heavy.tsv ; done;
for i in $a; do paste $i.name IGK/$i.info.lck > IGK/$i.lck.tsv ; done;
for i in $a; do paste $i.name IGL/$i.info.lcl > IGL/$i.lcl.tsv ; done;
cat IGH/*tsv > IGH.tsv
cat IGL/*tsv > IGL.tsv
cat IGK/*tsv > IGK.tsv
cat IGH.tsv | awk '$2 == "T" || $2 == "F"' > IGH_parsed.tsv
cat IGK.tsv | awk '$2 == "T" || $2 == "F"' > IGK_parsed.tsv
cat IGL.tsv | awk '$2 == "T" || $2 == "F"' > IGL_parsed.tsv
echo "done"
