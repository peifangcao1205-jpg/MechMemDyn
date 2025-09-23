/mnt/storage/Previous_Files/Transfer/nots_software/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol  native.xml -s 3gso.pdb -overwrite -ignore_zero_occupancy false -pH -pH_mode true -pH:value_pH 7 > native.log
cat native.log | grep "ResResE" > temp
mv temp native.log
let "j=0"
while read line; do
	let "j=j+1"
	sed -i "s/GKRSNTTGK/$line/" test.xml
	#/work/pw8/mc70/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol  test.xml -s 3gso.pdb
	/mnt/storage/Previous_Files/Transfer/nots_software/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol test.xml -s 3gso.pdb -overwrite -ignore_zero_occupancy false -pH -pH_mode true -pH:value_pH 7 > $j.log
	cat $j.log | grep "ResResE" > temp
	mv temp $j.log
	mv 3gso_0001.pdb $j.pdb
	sed -i "s/$line/GKRSNTTGK/" test.xml
done < peptide
python /home/ps/Scripts/python_tools/frustration/Rosetta/Frust_Post_v2_MultiChain_ManyBody_MultiLig_Gr_ContactList_7Dec.py 276 -2.5 0.5 300 9 0 Function1 -2.5 0.5
