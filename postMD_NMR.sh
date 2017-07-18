#!/bin/bash
for folder in `ls -d Docked*`
do 
    cd $folder
     echo "Doing $folder"
     grep "OME\|6VA\|0SA" md.pdb > ligand.pdb
     grep -v "OME\|6VA\|0SA" md.pdb > receptor.pdb
     h-h.exe receptor.pdb ligand.pdb > tmp
     
     grep "total" tmp | sort -k5 > NMR_Results.txt

     # Now sort via the experimental order and get a correlation value
     > sorted_NMR
     grep "6VA:H1 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "6VA:H3 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "6VA:H5 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "OME:NAc " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "6VA:H4 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "6VA:H6 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "6VA:H2 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H3E " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H5 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H9 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H8H9 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H3A " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H7 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:H4H6 " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
     grep "0SA:NAc " NMR_Results.txt | cut -d \  -f5 >> sorted_NMR
    cd ../
done
