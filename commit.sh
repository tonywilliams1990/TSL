#!/bin/sh
# Script for commiting TSL to GitHub
echo "Commit script"

make clean

# Move unnecessary files
mv ./figs ~/Desktop/figs_temp
mv ./DATA ~/Desktop/DATA_temp
mv ./petsc-3.9.1 ~/Desktop/petsc-3.9.1
mv ./slepc-3.9.1 ~/Desktop/slepc-3.9.1
mkdir ./DATA
cp ~/Desktop/DATA_temp/.gitkeep ./DATA/.gitkeep
rm -rf .sconf_temp

# Setup commit message
if [ $# -eq 0 ] ;
then
  message="empty"
else
  message=$1
fi

echo $message

# GitHub stuff
git add .
git commit -m $message
git push origin master
#git push origin Ridge
#git push origin Rayleigh

# Move the files back
mv ~/Desktop/figs_temp ./figs
mv ~/Desktop/petsc-3.9.1 ./petsc-3.9.1
mv ~/Desktop/slepc-3.9.1 ./slepc-3.9.1
rm -r DATA
mv ~/Desktop/DATA_temp ./DATA

exit 0
