#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage of the script:"
  echo "$0 [inputFile] [outputFile]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  exit
fi

INPUT=$1    # Name of the input file
OUTPUT=$2   # Name of the output file

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/projectTrackPairEfficiencyHistograms.C

# Project event information and jet histograms
root -l -b -q 'plotting/projectTrackPairEfficiencyHistograms.C("'${INPUT}'","'${OUTPUT}'",3)'

# Project track histograms
root -l -b -q 'plotting/projectTrackPairEfficiencyHistograms.C("'${INPUT}'","'${OUTPUT}'",4)'

# Project uncorrected track histograms
root -l -b -q 'plotting/projectTrackPairEfficiencyHistograms.C("'${INPUT}'","'${OUTPUT}'",8)'

# Project generator level particle histograms
root -l -b -q 'plotting/projectTrackPairEfficiencyHistograms.C("'${INPUT}'","'${OUTPUT}'",16)'

# Project the track pair histograms
root -l -b -q 'plotting/projectTrackPairEfficiencyHistograms.C("'${INPUT}'","'${OUTPUT}'",32)'

# Project the generator level particle pair histograms
root -l -b -q 'plotting/projectTrackPairEfficiencyHistograms.C("'${INPUT}'","'${OUTPUT}'",64)'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/projectTrackPairEfficiencyHistograms.C
