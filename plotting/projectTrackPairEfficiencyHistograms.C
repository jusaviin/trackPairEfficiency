#include "TrackPairEfficiencyCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TrackPairEfficiencyHistogramManager.h"

#include <bitset>

/*
 * Macro for projecting the histograms needed in the energy-energy correlator analysis from the THnSparses
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 */
void projectTrackPairEfficiencyHistograms(TString inputFileName = "veryCoolData.root", const char* outputFileName = "veryCoolData_processed.root", int histogramSelection = 127){

  // Print the file name to console
  cout << "Projecting histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // Choose which figure sets to draw
  bool loadEventInformation = false;
  bool loadJets = false;
  bool loadTracks = false;
  bool loadUncorrectedTracks = false;
  bool loadGenParticles = false;
  bool loadTrackPairHistograms = false;
  bool loadGenParticlePairHistograms = false;
  
  /*
   * Loading only selected histograms. Done with bitwise check of an integer
   *
   *  Bit 0 = Load event information histograms (to set: 1)
   *  Bit 1 = Load jet histograms (to set: 2)
   *  Bit 2 = Load track histograms (to set: 4)
   *  Bit 3 = Load uncorrected track histograms (to set: 8)
   *  Bit 4 = Load generator level particle histograms (to set: 16)
   *  Bit 5 = Load track pair histograms (to set: 32)
   *  Bit 6 = Load generator level particle pair histograms (to set: 64)
   */
  if(histogramSelection > 0){
    std::bitset<7> bitChecker(histogramSelection);
    loadEventInformation = bitChecker.test(0);
    loadJets = bitChecker.test(1);
    loadTracks = bitChecker.test(2);
    loadUncorrectedTracks = bitChecker.test(3);
    loadGenParticles = bitChecker.test(4);
    loadTrackPairHistograms = bitChecker.test(5);
    loadGenParticlePairHistograms = bitChecker.test(6);
  }
  
  // ====================================================
  //  Binning configuration for the projected histograms
  // ====================================================
  
  // Option to read all the binning information from TrackPairEfficiencyCard used to create the file
  const bool readCentralityBinsFromFile = false;
  const bool readTrackPtBinsFromFile = true;
  const bool readTrackPairPtBinsFromFile = true;
  
  // If not reading the bins from the file, manually define new bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nTrackPairPtBins = 7;
  const int nAverageEtaBins = 6;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,90};      // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,12,300};     // Bin borders for track pT
  double trackPairPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,12,300}; // Bin borders for track pT in track pair histograms
  double averageEtaBinBorders[nAverageEtaBins+1] = {-2.4, -1, -0.5, 0, 0.5, 1, 2.4};  // Bin borders for average eta slices
  
  // Projected bin range
  int firstProjectedCentralityBin = 0;
  int lastProjectedCentralityBin = nCentralityBins-1;
  
  int firstProjectedTrackPtBin = 0;
  int lastProjectedTrackPtBin = nTrackPtBins-1;
  
  int firstProjectedTrackPairPtBin = 0;
  int lastProjectedTrackPairPtBin = nTrackPairPtBins-1;
  
  int firstProjectedAverageEtaBin = 0;
  int lastProjectedAverageEtaBin = nAverageEtaBins;  // Histograms without eta selection are in the last bin

  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  TrackPairEfficiencyCard *card = new TrackPairEfficiencyCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    lastProjectedCentralityBin = 0;
    centralityBinBorders[0] = -0.5;
  }
  
  // If we change the binning, save the new binning to the card
  if(!readCentralityBinsFromFile) card->AddVector(TrackPairEfficiencyCard::kCentralityBinEdges,nCentralityBins+1,centralityBinBorders);
  if(!readTrackPtBinsFromFile) card->AddVector(TrackPairEfficiencyCard::kTrackPtBinEdges,nTrackPtBins+1,trackPtBinBorders);
  if(!readTrackPairPtBinsFromFile) card->AddVector(TrackPairEfficiencyCard::kTrackPairPtBinEdges,nTrackPairPtBins+1,trackPairPtBinBorders);
  card->AddVector(TrackPairEfficiencyCard::kAverageEtaBinEdges,nAverageEtaBins+1,averageEtaBinBorders);
  
  // Add information about the used input files to the card
  card->AddFileName(TrackPairEfficiencyCard::kInputFileName,inputFileName);
  
  // The git hash here will be replaced by the latest commit hash by projectHistogramsInSteps.sh script
  const char* gitHash = "";
  card->AddProjectionGitHash(gitHash);
  
  // ========================================== //
  //     TrackPairEfficiencyHistogramManager    //
  // ========================================== //
    
  // Create and setup a new histogram manager to project and handle the histograms
  TrackPairEfficiencyHistogramManager *histograms = new TrackPairEfficiencyHistogramManager(inputFile,card);
  
  // Set which histograms to project from the input file
  histograms->SetLoadEventInformation(loadEventInformation);
  histograms->SetLoadJetHistograms(loadJets);
  histograms->SetLoadTracks(loadTracks);
  histograms->SetLoadTracksUncorrected(loadUncorrectedTracks);
  histograms->SetLoadGenParticles(loadGenParticles);
  histograms->SetLoadTrackPairs(loadTrackPairHistograms);
  histograms->SetLoadGenParticlePairs(loadGenParticlePairHistograms);
  
  histograms->SetLoad2DHistograms(true);

  // Set the binning information
  histograms->SetCentralityBins(readCentralityBinsFromFile,nCentralityBins,centralityBinBorders,true);
  if(!readCentralityBinsFromFile) histograms->SetCentralityBinRange(firstProjectedCentralityBin,lastProjectedCentralityBin);
  histograms->SetTrackPtBins(readTrackPtBinsFromFile,nTrackPtBins,trackPtBinBorders,true);
  if(!readTrackPtBinsFromFile) histograms->SetTrackPtBinRange(firstProjectedTrackPtBin,lastProjectedTrackPtBin);
  histograms->SetTrackPairPtBins(readTrackPairPtBinsFromFile,nTrackPairPtBins,trackPairPtBinBorders,true);
  if(!readTrackPairPtBinsFromFile) histograms->SetTrackPairPtBinRange(firstProjectedTrackPairPtBin,lastProjectedTrackPairPtBin);
  histograms->SetAverageEtaBins(false,nAverageEtaBins,averageEtaBinBorders,true);
  histograms->SetAverageEtaBinRange(firstProjectedAverageEtaBin,lastProjectedAverageEtaBin);
  
  // Project the one dimensional histograms from the THnSparses
  histograms->LoadHistograms();
  
  // Save the histograms to an output file
  histograms->Write(outputFileName,fileWriteMode);
  
}
