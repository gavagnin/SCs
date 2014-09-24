#!/bin/bash


IC=sph_4
config=(sphere)
extent=(5)
echo "************************************"
echo "* Mergefiles_multi.py is running...*"
echo "************************************"
echo "Input params:"
echo "Files to be merged are in folder :" ${IC}
echo "Initial configuration: " ${config}
for Z in Z01 Z1
do 
    echo "Metallicity: " ${Z}
    for run in A B C D E F G H I L 
    do
	echo "Run: " ${run}
	out_file=(/Users/egavagnin/Starlab_M/YSC/${IC}/Merged/${Z}/merged_${Z}_${IC}_${extent}pc_run${run}.txt.gz)
	dir=(/Users/egavagnin/Starlab_M/YSC/${IC}/IC/${Z}/)

	cd /Users/egavagnin/Starlab_M/YSC/${IC}/IC/${Z}/
	input_files=(*run${run}*)
	echo "Number of files to merge: "${#input_files[@]}
	echo "Output merged file is: " ${out_file}
	python /Users/egavagnin/Starlab_M/my_scripts/mergefiles_multi.py ${#input_files[@]} $dir${input_files[@]} ${out_file} $config $extent
    done
done
