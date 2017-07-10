#!/bin/bash
    for folder in `ls -d Docked*`
    do 
        cd $folder
         echo "Doing $folder"
         grep "OME\|6VA\|0SA" md.pdb > ligand.pdb
         grep -v "OME\|6VA\|0SA" md.pdb > receptor.pdb
         h-h.exe receptor.pdb ligand.pdb
         
         grep "total" tmp | sort -k5 > NMR_Results.txt
        cd ../
    done
