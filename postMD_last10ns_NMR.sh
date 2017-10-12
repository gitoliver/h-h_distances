#!/bin/bash
for folder in `ls -d Docked*`
do 
    cd $folder
     mkdir Docked_last10ps_NMR
     cd Docked_last10ps_NMR
      rm *.out
      >$folder.out
      cpptraj -i ../../ptraj_last1ns_split.in
      for file in `ls ss*`
      do
          echo "Doing $folder/$file"
          grep "OME\|6VA\|0SA" $file > ligand.pdb
          grep -v "OME\|6VA\|0SA" $file > receptor.pdb
          h-h.exe receptor.pdb ligand.pdb > tmp
     
          grep "total" tmp | sort -k5 > NMR_Results.txt

          # Now sort via the experimental order and get a correlation value
          > sorted_NMR
          grep "6VA:H1 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "6VA:H3 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "6VA:H5 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "6VA:NAc " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "OME:NAc " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "6VA:H4 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "6VA:H6 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "6VA:H2 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H3E " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H5 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H9R " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H8H9 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H3A " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H7 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:H4H6 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
          grep "0SA:NAc " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
   
          cut -f2 /home/oliver/Programs/Cplusplus/h-h_distances/experimental_target.txt | sort -n > column2.txt
          r=`/home/oliver/Dropbox/1.glylib_scripts/zg.Correlation_Calc/correl.exe sorted_NMR column2.txt`
          echo "$r" >> $folder.out
      done
      cut -d = -f2 $folder.out > mean.out
      mean=`mean.exe mean.out norm`
      echo $mean > mean.out
     cd ../
    cd ../
done
