/*
 * Implementation of the TrackPairEfficiencyCard class
 */

// Own includes
#include "TrackPairEfficiencyCard.h"

/*
 * Contructor with input file
 *
 *  TFile *inFile = Input file
 */
TrackPairEfficiencyCard::TrackPairEfficiencyCard(TFile *inFile):
  fInputFile(inFile),
  fCardDirectory("JCard"),
  fDataType(-1),
  fDataTypeString(""),
  fAlternativeDataTypeString(""),
  fGitHash(0),
  fProjectionGitHash(0)
{

  // Read the vectors from the input file
  fInputFile->cd(fCardDirectory.Data());
  ReadVectors();
  fDataType = (*fCardEntries[kDataType])[1];
  
  FindDataTypeString();
  
}

/*
 * Destructor
 */
TrackPairEfficiencyCard::~TrackPairEfficiencyCard(){

}

/*
 * Reader for all the vectors from the input file
 */
void TrackPairEfficiencyCard::ReadVectors(){
  
  // Read the git hash
  fGitHash = (TObjString*) gDirectory->Get("GitHash");
  fProjectionGitHash = (TObjString*) gDirectory->Get("ProjectionGitHash");

  // Read the TVectorT<float>:s
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    fCardEntries[iEntry] = (TVectorT<float>*) gDirectory->Get(fCardEntryNames[iEntry]);
  }
  
  // Read the file names
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    fFileNames[iFileName] = (TObjString*) gDirectory->Get(fFileNameSaveName[iFileName]);
  }

}

/*
 * Construct a data type string based on information on the card
 */
void TrackPairEfficiencyCard::FindDataTypeString(){
  
  // Define the different data types corresponding to certain indices
  TString dataTypes[5] = {"pp","PbPb","pp MC","PbPb MC","localTest"};
  TString alternativeDataTypes[5] = {"pp","PbPb","Pythia8","Pythia+Hydjet","localTest"};
  if(fDataType < 0 || fDataType > 4){
    fDataTypeString = "Unknown";
    fAlternativeDataTypeString = "Unknown";
    return;
  }
  
  // Remember the constructed data type string
  fDataTypeString = dataTypes[fDataType];
  fAlternativeDataTypeString = alternativeDataTypes[fDataType];
}

/*
 *  Getter for data type string
 */
TString TrackPairEfficiencyCard::GetDataType() const{
  return fDataTypeString;
}

/*
 *  Getter for alternative data type string
 */
TString TrackPairEfficiencyCard::GetAlternativeDataType() const{
  return fAlternativeDataTypeString;
}

/*
 *  Getter for jet type
 */
int TrackPairEfficiencyCard::GetJetType() const{
  return (*fCardEntries[kJetType])[1];
}

// Getter for the minimum jet pT cut
double TrackPairEfficiencyCard::GetJetPtCut() const{
  return (*fCardEntries[kMinPtCut])[1];
}

/*
 * Get the number of bins for internal index
 * If no vector is found in the index, return 0.
 */
int TrackPairEfficiencyCard::GetNBins(const int index) const{
  if(fCardEntries[index]) return fCardEntries[index]->GetNoElements()-1;
  return 0;
}

// Get the number of centrality bins
int TrackPairEfficiencyCard::GetNCentralityBins() const{
  return GetNBins(kCentralityBinEdges);
}

// Get the number of track pT bins
int TrackPairEfficiencyCard::GetNTrackPtBins() const{
  return GetNBins(kTrackPtBinEdges);
}

// Get the number of track pT bins in track pair histograms
int TrackPairEfficiencyCard::GetNTrackPairPtBins() const{
  return GetNBins(kTrackPairPtBinEdges);
}

/*
 * Get a bin index based on a given value.
 * If value is out of bound, return -1
 */
int TrackPairEfficiencyCard::GetBinIndex(const int index, const double value) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);
  
  // If the number given is smaller than the first value, it is out of bounds. Return -1 in this case.
  if(value < (*fCardEntries[index])[1]) return -1;
  
  // Find the bin in which the given value is and return it
  for(int iBin = 2; iBin <= nBins+1; iBin++){
    if(value < (*fCardEntries[index])[iBin]) return iBin-2;
  }
  
  // If a value was not found, it must be larger than anything in the array. Return -1 in this case.
  return -1;
}

// Get the bin index for a given centrality value
int TrackPairEfficiencyCard::GetBinIndexCentrality(const double value) const{
  return GetBinIndex(kCentralityBinEdges,value);
}

// Get the bin index for a given track pT value
int TrackPairEfficiencyCard::GetBinIndexTrackPt(const double value) const{
  return GetBinIndex(kTrackPtBinEdges,value);
}

// Get the bin index for a given track pT value in track pair histograms
int TrackPairEfficiencyCard::GetBinIndexTrackPairPt(const double value) const{
  return GetBinIndex(kTrackPairPtBinEdges,value);
}

// Get the low border of i:th bin from internal index
double TrackPairEfficiencyCard::GetLowBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin > GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+1];
  return -1;
}

// Get the low border of i:th centrality bin
double TrackPairEfficiencyCard::GetLowBinBorderCentrality(const int iBin) const{
  return GetLowBinBorder(kCentralityBinEdges,iBin);
}

// Get the low border of i:th track pT bin
double TrackPairEfficiencyCard::GetLowBinBorderTrackPt(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdges,iBin);
}

// Get the low border of i:th track pT bin in track pair histograms
double TrackPairEfficiencyCard::GetLowBinBorderTrackPairPt(const int iBin) const{
  return GetLowBinBorder(kTrackPairPtBinEdges,iBin);
}

// Get the high border of i:th bin from internal index
double TrackPairEfficiencyCard::GetHighBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+2];
  return -1;
}

// Get the high border of i:th centrality bin
double TrackPairEfficiencyCard::GetHighBinBorderCentrality(const int iBin) const{
  return GetHighBinBorder(kCentralityBinEdges,iBin);
}

// Get the high border of i:th track pT bin
double TrackPairEfficiencyCard::GetHighBinBorderTrackPt(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdges,iBin);
}

// Get the high border of i:th track pT bin in track pair histograms
double TrackPairEfficiencyCard::GetHighBinBorderTrackPairPt(const int iBin) const{
  return GetHighBinBorder(kTrackPairPtBinEdges,iBin);
}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  float entryContent = Content to be put into the vector
 */
void TrackPairEfficiencyCard::AddOneDimensionalVector(int entryIndex, float entryContent){
  
  // Only allow addition to postprocessing vectors
  if(entryIndex <= kPtHatBinEdges) return;
  
  // Make a new one dimensional vector to the desired index with given content
  float contents[1] = {entryContent};
  fCardEntries[entryIndex] = new TVectorT<float>(1,1,contents);

}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  int dimension = Number of entries in the given array
 *  float *contents = Content to be put into the vector
 */
void TrackPairEfficiencyCard::AddVector(int entryIndex, int dimension, double *contents){
  
  // Convert double pointer to float pointer
  float* convertedContents = new float[dimension];
  for(int i = 0; i < dimension; i++){
    convertedContents[i] = contents[i];
  }
  
  // Make a new vector to the desired index with given content
  fCardEntries[entryIndex] = new TVectorT<float>(1,dimension,convertedContents);
  
  // Delete the converted contents array
  delete[] convertedContents;
  
}

/*
 * Add file name to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added file name in file name array
 *  TString fileName = Added file name
 */
void TrackPairEfficiencyCard::AddFileName(int entryIndex, TString fileName){
  
  // Make convert the string to TObjString and add it to the file name array
  fFileNames[entryIndex] = new TObjString(fileName.Data());
  
}

/*
 * Add git hash for the projecting code used to project the histograms to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for projection
 */
void TrackPairEfficiencyCard::AddProjectionGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and add it to the file name array
  fProjectionGitHash = new TObjString(gitHash);
  
}

/*
 * Print the contents of the card to the console
 */
void TrackPairEfficiencyCard::Print() const{

  const char* gitHash = "Unknown";
  if(fGitHash != NULL) gitHash = fGitHash->String().Data();
  
  std::cout<<std::endl<<"========================= TrackPairEfficiencyCard =========================="<<std::endl;
  std::cout << Form("%25s","GitHash"); //print keyword
  std::cout << ": " << gitHash << std::endl;
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry]){
      std::cout << Form("%25s",fCardEntryNames[iEntry]); //print keyword
      std::cout << " (dim = "<<fCardEntries[iEntry]->GetNoElements() << ") ";//print size of TVector
      for(int iElement = 1; iElement <= fCardEntries[iEntry]->GetNoElements(); iElement++){
        std::cout << (*fCardEntries[iEntry])[iElement] << " ";//TVector components
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;
  
  if(fProjectionGitHash != NULL) std::cout << "Git hash for projections: " << fProjectionGitHash->String().Data() << std::endl;
  
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]){
      std::cout << "Used " << fFileNameType[iFileName] << " file: " << fFileNames[iFileName]->String().Data() << std::endl;
    }
  }
}

/*
 * Write the contents of the TrackPairEfficiencyCard to a file
 */
void TrackPairEfficiencyCard::Write(TDirectory *file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash to the file
  if(fGitHash) fGitHash->Write("GitHash");
  
  // Write all the vectors to the file. Not all of these exist in older versions of cards, thus check if exists before writing.
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }
   
  // Write the git hash used for projections to the file
  if(fProjectionGitHash) fProjectionGitHash->Write("ProjectionGitHash");
  
  // Write all the data names to the file.
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]) fFileNames[iFileName]->Write(fFileNameSaveName[iFileName]);
  }
  
  // Return back to the main directory
  file->cd("../");
}
