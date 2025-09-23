cat $1 | awk '/ATOM/ && $3 == "CA"  {print $4}' | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g' > native.seq
python RandSeq.py 300 native.seq > peptide
cat $1 | awk '/ATOM/ && $3 == "CA"  {print $6}' > idx
cat $1 | awk '/ATOM/ && $3 == "CA"  {print $5}' > cid
var=$(tail -1 idx)
var2=$(tail -1 cid)
var3=$(head -1 idx)
var4=$(head -1 cid)
var5=$(cat idx | wc -l)
sed -i "s/92A/${var}${var2}/" test.xml
sed -i "s/92A/${var}${var2}/" native.xml
sed -i "s/1A/${var3}${var4}/" test.xml
sed -i "s/1A/${var3}${var4}/" native.xml
sed -i "s/\"1753\"/\"${var5}\"/" native.xml
sed -i "s/\"1753\"/\"${var5}\"/" test.xml
cp $1 3gso.pdb
#sbatch job.pbs
