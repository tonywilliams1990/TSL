#!/bin/sh
# Script for commiting TSL to GitHub
echo "Commit script"

make clean

# Move unnecessary files
mv ./figs ~/Desktop/figs_temp
mv ./DATA ~/Desktop/DATA_temp
mkdir ./DATA
cp ~/Desktop/DATA_temp/.gitkeep ./DATA/.gitkeep

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
#git push origin master
git push origin Ridge

# Move the files back
mv ~/Desktop/figs_temp ./figs
rm -r DATA
mv ~/Desktop/DATA_temp ./DATA

exit 0
