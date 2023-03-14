#!/bin/bash

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Define the name for the output tar ball
OUTPUTTAR="trackPairEfficiencyAnalysis.tar.gz"

# If the tar ball already exists, remove it
if [ -f $OUTPUTTAR ]; then 
rm $OUTPUTTAR 
fi

# Replace the placeholder string in the main analysis file by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' trackPairEfficiencyAnalysis.cxx

# Make sure there are no object files going to the tar
make clean

# Create the new tar ball
tar -cvzf $OUTPUTTAR Makefile trackPairEfficiencyAnalysis.cxx src trackCorrectionTables

# Put placeholder string back to the main analysis file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' trackPairEfficiencyAnalysis.cxx
