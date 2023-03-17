#include "TrackPairEfficiencyHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TrackPairEfficiencyCard.h"
#include "JDrawer.h"
#include "../src/TrackPairEfficiencyHistograms.h"

/*
 * Macro for determining track pair efficiency from reconstructed and generator level simulations
 */
void trackPairEfficiencyPlotter(){
  
  // File containing the track pair distributions used to determine the track pair efficiency
  TString fileName = "data/trackPairEfficiency_triggerFiducialCut_processed_2023-03-16.root";
  // trackPairEfficiency_triggerFiducialCut_processed_2023-03-16.root
  // trackPairEfficiency_wholeTracker_updatePt_first500_projected_2023-03-15.root
  
  // Open the file and check that it exists
  TFile* inputFile = TFile::Open(fileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Check if we are using PbPb or pp data
  TrackPairEfficiencyCard *systemCard = new TrackPairEfficiencyCard(inputFile);
  TString collisionSystem = systemCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  const int nCentralityBins = (isPbPbData) ? systemCard->GetNCentralityBins() : 1;
  const int nTrackPtBins = systemCard->GetNTrackPairPtBins();
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Select the bin range to be drawn
  const int firstDrawCentralityBin = 0;
  const int lastDrawnCentralityBin = 0;
  
  const int firstDrawnTriggerPtBin = 0;
  const int lastDrawnTriggerPtBin = nTrackPtBins-1;
  
  const int firstDrawnAssociatedPtBin = 0;
  const int lastDrawnAssociatedPtBin = nTrackPtBins-1;
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Histogram saving
  const bool writeSmoothedHistograms = true; // Write flag for smoothed histograms
  
  // Create and setup a new histogram managers to project and handle the histograms
  TrackPairEfficiencyHistogramManager *histograms;
  histograms = new TrackPairEfficiencyHistogramManager(inputFile);
  histograms->SetCentralityBinRange(firstDrawCentralityBin, lastDrawnCentralityBin);
  histograms->SetTrackPairPtBinRange(TMath::Min(firstDrawnTriggerPtBin, firstDrawnAssociatedPtBin), TMath::Max(lastDrawnTriggerPtBin, lastDrawnAssociatedPtBin));
  histograms->SetLoadAllTrackPairs(true,true);
  histograms->LoadProcessedHistograms();
  
  // Jet pT histograms to calculate trigger turn on curve
  TH1D* hTrackDeltaR[nCentralityBins][nTrackPtBins][nTrackPtBins][TrackPairEfficiencyHistograms::knDataLevels];   // DeltaR distribution for track pairs
  TH1D* hRatioDeltaR[nCentralityBins][nTrackPtBins][nTrackPtBins]; // Reconstructed to generator level ratio of track pair DeltaR distributions
  TH1D* hRatioDeltaRSmoothed[nCentralityBins][nTrackPtBins][nTrackPtBins]; // Smoothed reconstructed to generator level ratio of track pair DeltaR distributions
  
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < nTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < nTrackPtBins; iAssociatedPt++){
        for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
          hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iDataLevel] = NULL;
        } // Data level loop
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt] = NULL;
        hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt] = NULL;
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Read the histograms from the file
  for(int iCentrality = firstDrawCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
          hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iDataLevel] = histograms->GetHistogramTrackPairDeltaR(iCentrality, iTriggerPt, iAssociatedPt, iDataLevel);
        } // Data level loop
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Calculate the generator level to reconstructed ratio to determine track pair efficiency
  double genLevelIntegral;
  double integralToMatch;
  int lowBinLimit = hTrackDeltaR[firstDrawCentralityBin][firstDrawnTriggerPtBin][firstDrawnAssociatedPtBin][0]->FindBin(0.1);
  int highBinLimit = hTrackDeltaR[firstDrawCentralityBin][firstDrawnTriggerPtBin][firstDrawnAssociatedPtBin][0]->FindBin(0.4);
  for(int iCentrality = firstDrawCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        
        // First, normalize the distributions such that the integrals match in the region 0.1-0.4
        integralToMatch = hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kReconstructed]->Integral(lowBinLimit, highBinLimit, "width");
        genLevelIntegral = hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kGeneratorLevel]->Integral(lowBinLimit, highBinLimit, "width");
        hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kGeneratorLevel]->Scale(integralToMatch/genLevelIntegral);
        
        // Then, calculate the ratio of the histograms
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt] = (TH1D*) hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kReconstructed]->Clone(Form("ratio%d%d%d", iCentrality, iTriggerPt, iAssociatedPt));
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt]->Divide(hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kGeneratorLevel]);
        
        // Finally, smooth the distribution to suppress fluctuations
        hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt] = (TH1D*) hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt]->Clone(Form("smoothedTrackPairEfficiencyCorrection_C%dT%dA%d", iCentrality, iTriggerPt, iAssociatedPt));
        hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt]->Smooth();
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Manual extra smoothing. The smoothing does generally good job, but some histogram require a manual extra touch to get the smoothing satisfactory
  // The values here fix the smoothing for file trackPairEfficiency_triggerFiducialCut_processed_2023-03-16.root
  
  if(firstDrawCentralityBin == 0){
    
    // Centrality 0-10, trigger pT 3-4, associated pT 3-4
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(5, 0.64);
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(6, 0.672);
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(7, 0.704);
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(8, 0.736);
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(9, 0.768);
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(10, 0.8);
    hRatioDeltaRSmoothed[0][3][3]->SetBinContent(11, 0.832);
    
    // Centrality 0-10, trigger pT 4-6, associated pT 2-3
    hRatioDeltaRSmoothed[0][4][2]->SetBinContent(1, 0.4);
    hRatioDeltaRSmoothed[0][4][2]->SetBinContent(2, 0.49);
    hRatioDeltaRSmoothed[0][4][2]->SetBinContent(3, 0.56);
    hRatioDeltaRSmoothed[0][4][2]->SetBinContent(4, 0.63);
    hRatioDeltaRSmoothed[0][4][2]->SetBinContent(5, 0.69);
    
    // Centrality 0-10, trigger pT 4-6, associated pT 3-4
    hRatioDeltaRSmoothed[0][4][3]->SetBinContent(2, 0.59);
    hRatioDeltaRSmoothed[0][4][3]->SetBinContent(3, 0.62);
    hRatioDeltaRSmoothed[0][4][3]->SetBinContent(4, 0.65);
    hRatioDeltaRSmoothed[0][4][3]->SetBinContent(5, 0.68);
    hRatioDeltaRSmoothed[0][4][3]->SetBinContent(6, 0.705);
    
    // Centrality 0-10, trigger pT 6-8, associated pT 2-3
    hRatioDeltaRSmoothed[0][5][2]->SetBinContent(2, 0.55);
    hRatioDeltaRSmoothed[0][5][2]->SetBinContent(3, 0.59);
    hRatioDeltaRSmoothed[0][5][2]->SetBinContent(4, 0.64);
    hRatioDeltaRSmoothed[0][5][2]->SetBinContent(5, 0.66);
    hRatioDeltaRSmoothed[0][5][2]->SetBinContent(6, 0.68);
    hRatioDeltaRSmoothed[0][5][2]->SetBinContent(7, 0.7);
    
    // Centrality 0-10, trigger pT 6-8, associated pT 3-4
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(2, 0.36);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(3, 0.41);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(4, 0.46);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(5, 0.5);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(6, 0.56);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(7, 0.62);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(8, 0.68);
    hRatioDeltaRSmoothed[0][5][3]->SetBinContent(9, 0.73);
    
    // Centrality 0-10, trigger pT 8-12, associated pT 2-3
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(2, 0.45);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(3, 0.51);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(4, 0.57);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(5, 0.63);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(6, 0.69);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(7, 0.72);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(8, 0.75);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(9, 0.78);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(10, 0.82);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(11, 0.85);
    hRatioDeltaRSmoothed[0][6][2]->SetBinContent(12, 0.88);
    
    // Centrality 0-10, trigger pT 8-12, associated pT 4-6
    hRatioDeltaRSmoothed[0][6][4]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[0][6][4]->SetBinContent(2, 0.4);
    hRatioDeltaRSmoothed[0][6][4]->SetBinContent(3, 0.5);
    hRatioDeltaRSmoothed[0][6][4]->SetBinContent(4, 0.6);
    
    // Centrality 0-10, trigger pT 12-16, associated pT 4-6
    hRatioDeltaRSmoothed[0][7][4]->SetBinContent(1, 0.2);
    hRatioDeltaRSmoothed[0][7][4]->SetBinContent(2, 0.32);
    hRatioDeltaRSmoothed[0][7][4]->SetBinContent(3, 0.4);
    hRatioDeltaRSmoothed[0][7][4]->SetBinContent(4, 0.47);
    hRatioDeltaRSmoothed[0][7][4]->SetBinContent(5, 0.54);
    
    // Centrality 0-10, trigger pT 12-16, associated pT 12-16
    hRatioDeltaRSmoothed[0][7][7]->SetBinContent(1, 0.2);
    hRatioDeltaRSmoothed[0][7][7]->SetBinContent(2, 0.29);
    hRatioDeltaRSmoothed[0][7][7]->SetBinContent(3, 0.37);
    hRatioDeltaRSmoothed[0][7][7]->SetBinContent(4, 0.43);
    hRatioDeltaRSmoothed[0][7][7]->SetBinContent(5, 0.48);
    
    // Centrality 0-10, trigger pT 16-20 associated pT 2-3
    hRatioDeltaRSmoothed[0][8][2]->SetBinContent(1, 0.35);
    hRatioDeltaRSmoothed[0][8][2]->SetBinContent(2, 0.44);
    hRatioDeltaRSmoothed[0][8][2]->SetBinContent(3, 0.52);
    hRatioDeltaRSmoothed[0][8][2]->SetBinContent(4, 0.59);
    hRatioDeltaRSmoothed[0][8][2]->SetBinContent(5, 0.65);
    
    // Centrality 0-10, trigger pT 16-20 associated pT 3-4
    hRatioDeltaRSmoothed[0][8][3]->SetBinContent(3, 0.3);
    hRatioDeltaRSmoothed[0][8][3]->SetBinContent(3, 0.4);
    hRatioDeltaRSmoothed[0][8][3]->SetBinContent(4, 0.47);
    hRatioDeltaRSmoothed[0][8][3]->SetBinContent(5, 0.54);
    hRatioDeltaRSmoothed[0][8][3]->SetBinContent(6, 0.6);
    
    // Centrality 0-10, trigger pT 16-20 associated pT 8-12
    hRatioDeltaRSmoothed[0][8][6]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[0][8][6]->SetBinContent(2, 0.37);
    hRatioDeltaRSmoothed[0][8][6]->SetBinContent(3, 0.44);
    hRatioDeltaRSmoothed[0][8][6]->SetBinContent(4, 0.5);
    
    // Centrality 0-10, trigger pT 16-20 associated pT 12-16
    hRatioDeltaRSmoothed[0][8][7]->SetBinContent(1, 0.25);
    hRatioDeltaRSmoothed[0][8][7]->SetBinContent(2, 0.32);
    hRatioDeltaRSmoothed[0][8][7]->SetBinContent(3, 0.39);
    hRatioDeltaRSmoothed[0][8][7]->SetBinContent(4, 0.45);
    
    // Centrality 0-10, trigger pT 20-300 associated pT 2-3
    hRatioDeltaRSmoothed[0][9][2]->SetBinContent(1, 0.4);
    
    // Centrality 0-10, trigger pT 20-300 associated pT 3-4
    hRatioDeltaRSmoothed[0][9][3]->SetBinContent(2, 0.38);
    hRatioDeltaRSmoothed[0][9][3]->SetBinContent(3, 0.6);
    hRatioDeltaRSmoothed[0][9][3]->SetBinContent(4, 0.62);
    hRatioDeltaRSmoothed[0][9][3]->SetBinContent(5, 0.635);
    hRatioDeltaRSmoothed[0][9][3]->SetBinContent(6, 0.65);
    hRatioDeltaRSmoothed[0][9][3]->SetBinContent(7, 0.66);
    
    // Centrality 0-10, trigger pT 20-300 associated pT 4-6
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(1, 0.2);
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(2, 0.26);
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(3, 0.32);
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(4, 0.38);
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(5, 0.44);
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(6, 0.49);
    hRatioDeltaRSmoothed[0][9][4]->SetBinContent(7, 0.54);
  }
  
  if(firstDrawCentralityBin == 1 || (firstDrawCentralityBin < 1 && lastDrawnCentralityBin >= 1)){
    
    // Centrality 10-30, trigger pT 4-6, associated pT 3-4
    hRatioDeltaRSmoothed[1][4][3]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[1][4][3]->SetBinContent(2, 0.37);
    hRatioDeltaRSmoothed[1][4][3]->SetBinContent(3, 0.42);
    hRatioDeltaRSmoothed[1][4][3]->SetBinContent(4, 0.47);
    hRatioDeltaRSmoothed[1][4][3]->SetBinContent(5, 0.51);
    
    // Centrality 10-30, trigger pT 8-12, associated pT 3-4
    hRatioDeltaRSmoothed[1][6][3]->SetBinContent(1, 0.35);
    hRatioDeltaRSmoothed[1][6][3]->SetBinContent(2, 0.52);
    
    // Centrality 10-30, trigger pT 8-12, associated pT 8-12
    hRatioDeltaRSmoothed[1][6][6]->SetBinContent(1, 0.4);
    hRatioDeltaRSmoothed[1][6][6]->SetBinContent(2, 0.42);
    hRatioDeltaRSmoothed[1][6][6]->SetBinContent(3, 0.44);
    hRatioDeltaRSmoothed[1][6][6]->SetBinContent(4, 0.46);
   
    // Centrality 10-30, trigger pT 12-16, associated pT 3-4
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(1, 0.2);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(2, 0.34);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(3, 0.45);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(4, 0.55);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(5, 0.60);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(6, 0.65);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(7, 0.70);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(8, 0.74);
    hRatioDeltaRSmoothed[1][7][3]->SetBinContent(9, 0.78);
    
    // Centrality 10-30, trigger pT 16-20, associated pT 2-3
    hRatioDeltaRSmoothed[1][8][2]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[1][8][2]->SetBinContent(2, 0.55);
    
    // Centrality 10-30, trigger pT 16-20, associated pT 3-4
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(1, 0.25);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(2, 0.32);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(3, 0.45);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(4, 0.55);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(5, 0.60);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(6, 0.65);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(7, 0.70);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(8, 0.74);
    hRatioDeltaRSmoothed[1][8][3]->SetBinContent(9, 0.78);
    
    // Centrality 10-30, trigger pT 16-20, associated pT 4-6
    hRatioDeltaRSmoothed[1][8][4]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[1][8][4]->SetBinContent(2, 0.46);
    hRatioDeltaRSmoothed[1][8][4]->SetBinContent(3, 0.52);
    hRatioDeltaRSmoothed[1][8][4]->SetBinContent(4, 0.56);
    hRatioDeltaRSmoothed[1][8][4]->SetBinContent(5, 0.6);
    
    // Centrality 10-30, trigger pT 16-20, associated pT 6-8
    hRatioDeltaRSmoothed[1][8][5]->SetBinContent(1, 0.4);
    hRatioDeltaRSmoothed[1][8][5]->SetBinContent(2, 0.43);
    hRatioDeltaRSmoothed[1][8][5]->SetBinContent(3, 0.46);
    hRatioDeltaRSmoothed[1][8][5]->SetBinContent(4, 0.52);
    hRatioDeltaRSmoothed[1][8][5]->SetBinContent(5, 0.56);
    hRatioDeltaRSmoothed[1][8][5]->SetBinContent(6, 0.6);
    
    // Centrality 10-30, trigger pT 16-20, associated pT 16-20
    hRatioDeltaRSmoothed[1][8][8]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[1][8][8]->SetBinContent(2, 0.35);
    
    // Centrality 10-30, trigger pT 20-300, associated pT 4-6
    hRatioDeltaRSmoothed[1][9][4]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[1][9][4]->SetBinContent(2, 0.38);
    hRatioDeltaRSmoothed[1][9][4]->SetBinContent(3, 0.44);
    hRatioDeltaRSmoothed[1][9][4]->SetBinContent(4, 0.52);
    hRatioDeltaRSmoothed[1][9][4]->SetBinContent(5, 0.57);
    hRatioDeltaRSmoothed[1][9][4]->SetBinContent(6, 0.62);
    
    // Centrality 10-30, trigger pT 20-300, associated pT 6-8
    hRatioDeltaRSmoothed[1][9][5]->SetBinContent(1, 0.5);
    hRatioDeltaRSmoothed[1][9][5]->SetBinContent(2, 0.55);
    
    // Centrality 10-30, trigger pT 20-300, associated pT 8-12
    hRatioDeltaRSmoothed[1][9][6]->SetBinContent(1, 0.32);
    hRatioDeltaRSmoothed[1][9][6]->SetBinContent(2, 0.42);
    hRatioDeltaRSmoothed[1][9][6]->SetBinContent(3, 0.5);
    hRatioDeltaRSmoothed[1][9][6]->SetBinContent(4, 0.56);
    
    // Centrality 10-30, trigger pT 20-300, associated pT 16-20
    hRatioDeltaRSmoothed[1][9][8]->SetBinContent(1, 0.4);
    hRatioDeltaRSmoothed[1][9][8]->SetBinContent(2, 0.45);
    hRatioDeltaRSmoothed[1][9][8]->SetBinContent(3, 0.5);
    hRatioDeltaRSmoothed[1][9][8]->SetBinContent(4, 0.55);
  }
  
  if(firstDrawCentralityBin == 2 || (firstDrawCentralityBin < 2 && lastDrawnCentralityBin >= 2)){
    
    // Centrality 30-50, trigger pT 4-6, associated pT 3-4
    hRatioDeltaRSmoothed[2][4][3]->SetBinContent(1, 0.28);
    hRatioDeltaRSmoothed[2][4][3]->SetBinContent(2, 0.5);
    hRatioDeltaRSmoothed[2][4][3]->SetBinContent(3, 0.54);
    hRatioDeltaRSmoothed[2][4][3]->SetBinContent(4, 0.58);
    hRatioDeltaRSmoothed[2][4][3]->SetBinContent(5, 0.62);
    
    // Centrality 30-50, trigger pT 4-6, associated pT 4-6
    hRatioDeltaRSmoothed[2][4][4]->SetBinContent(1, 0.2);
    hRatioDeltaRSmoothed[2][4][4]->SetBinContent(2, 0.38);
    hRatioDeltaRSmoothed[2][4][4]->SetBinContent(3, 0.44);
    hRatioDeltaRSmoothed[2][4][4]->SetBinContent(5, 0.5);
    
    // Centrality 30-50, trigger pT 6-8, associated pT 3-4
    hRatioDeltaRSmoothed[2][5][3]->SetBinContent(1, 0.2);
    
    // Centrality 30-50, trigger pT 6-8, associated pT 6-8
    hRatioDeltaRSmoothed[2][5][5]->SetBinContent(1, 0.15);
    hRatioDeltaRSmoothed[2][5][5]->SetBinContent(2, 0.4);
    
    // Centrality 30-50, trigger pT 8-12, associated pT 2-3
    hRatioDeltaRSmoothed[2][6][2]->SetBinContent(1, 0.4);
    
    // Centrality 30-50, trigger pT 8-12, associated pT 3-4
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(1, 0.25);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(2, 0.45);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(3, 0.52);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(4, 0.55);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(5, 0.58);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(6, 0.61);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(7, 0.64);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(8, 0.67);
    hRatioDeltaRSmoothed[2][6][3]->SetBinContent(9, 0.7);
    
    // Centrality 30-50, trigger pT 8-12, associated pT 4-6
    hRatioDeltaRSmoothed[2][6][4]->SetBinContent(1, 0.3);
    hRatioDeltaRSmoothed[2][6][4]->SetBinContent(2, 0.52);
    hRatioDeltaRSmoothed[2][6][4]->SetBinContent(3, 0.57);
    hRatioDeltaRSmoothed[2][6][4]->SetBinContent(4, 0.62);
    
    // Centrality 30-50, trigger pT 12-16, associated pT 2-3
    hRatioDeltaRSmoothed[2][7][2]->SetBinContent(1, 0.25);
    hRatioDeltaRSmoothed[2][7][2]->SetBinContent(2, 0.52);
    hRatioDeltaRSmoothed[2][7][2]->SetBinContent(3, 0.65);
    
    // Centrality 30-50, trigger pT 12-16, associated pT 12-16
    hRatioDeltaRSmoothed[2][7][7]->SetBinContent(1, 0.2);
    hRatioDeltaRSmoothed[2][7][7]->SetBinContent(2, 0.5);
    hRatioDeltaRSmoothed[2][7][7]->SetBinContent(3, 0.58);
    
    // Centrality 30-50, trigger pT 16-20, associated pT 2-3
    hRatioDeltaRSmoothed[2][8][2]->SetBinContent(1, 0.2);
    
    // Centrality 30-50, trigger pT 16-20, associated pT 3-4
    hRatioDeltaRSmoothed[2][8][3]->SetBinContent(1, 0.2);
    
    // Centrality 30-50, trigger pT 16-20, associated pT 4-6
    hRatioDeltaRSmoothed[2][8][4]->SetBinContent(2, 0.55);
    hRatioDeltaRSmoothed[2][8][4]->SetBinContent(3, 0.6);
    hRatioDeltaRSmoothed[2][8][4]->SetBinContent(4, 0.63);
    hRatioDeltaRSmoothed[2][8][4]->SetBinContent(5, 0.66);
    hRatioDeltaRSmoothed[2][8][4]->SetBinContent(6, 0.69);
    
    // Centrality 30-50, trigger pT 16-20, associated pT 8-12
    hRatioDeltaRSmoothed[2][8][6]->SetBinContent(1, 0.42);
    hRatioDeltaRSmoothed[2][8][6]->SetBinContent(2, 0.47);
    hRatioDeltaRSmoothed[2][8][6]->SetBinContent(3, 0.52);
    hRatioDeltaRSmoothed[2][8][6]->SetBinContent(4, 0.57);
    
    // Centrality 30-50, trigger pT 16-20, associated pT 16-20
    hRatioDeltaRSmoothed[2][8][8]->SetBinContent(1, 0.38);
    hRatioDeltaRSmoothed[2][8][8]->SetBinContent(2, 0.45);
    
    // Centrality 30-50, trigger pT 20-300, associated pT 8-12
    hRatioDeltaRSmoothed[2][9][6]->SetBinContent(1, 0.52);
    hRatioDeltaRSmoothed[2][9][6]->SetBinContent(2, 0.62);
    hRatioDeltaRSmoothed[2][9][6]->SetBinContent(3, 0.65);
    hRatioDeltaRSmoothed[2][9][6]->SetBinContent(4, 0.68);
  }
  
  if(firstDrawCentralityBin == 3 || (firstDrawCentralityBin < 3 && lastDrawnCentralityBin >= 3)){
   
    // Centrality 50-90, trigger pT 2-3, associated pT 2-3
    hRatioDeltaRSmoothed[3][2][2]->SetBinContent(1, 0.18);
    
    // Centrality 50-90, trigger pT 3-4, associated pT 2-3
    hRatioDeltaRSmoothed[3][3][2]->SetBinContent(1, 0.18);
    
    // Centrality 50-90, trigger pT 4-6, associated pT 4-6
    hRatioDeltaRSmoothed[3][4][4]->SetBinContent(1, 0.16);
    
    // Centrality 50-90, trigger pT 6-8, associated pT 4-6
    hRatioDeltaRSmoothed[3][5][4]->SetBinContent(1, 0.3);
    
    // Centrality 50-90, trigger pT 6-8, associated pT 6-8
    hRatioDeltaRSmoothed[3][5][5]->SetBinContent(1, 0.18);
    
    // Centrality 50-90, trigger pT 8-12, associated pT 2-3
    hRatioDeltaRSmoothed[3][6][2]->SetBinContent(1, 0.2);
    
    // Centrality 50-90, trigger pT 8-12, associated pT 3-4
    hRatioDeltaRSmoothed[3][6][3]->SetBinContent(1, 0.2);
    
    // Centrality 50-90, trigger pT 8-12, associated pT 4-6
    hRatioDeltaRSmoothed[3][6][4]->SetBinContent(1, 0.35);
    
    // Centrality 50-90, trigger pT 8-12, associated pT 6-8
    hRatioDeltaRSmoothed[3][6][5]->SetBinContent(1, 0.25);
    
    // Centrality 50-90, trigger pT 12-16, associated pT 2-3
    hRatioDeltaRSmoothed[3][7][2]->SetBinContent(1, 0.42);
    
    // Centrality 50-90, trigger pT 12-16, associated pT 8-12
    hRatioDeltaRSmoothed[3][7][6]->SetBinContent(1, 0.48);
    hRatioDeltaRSmoothed[3][7][6]->SetBinContent(2, 0.54);
    hRatioDeltaRSmoothed[3][7][6]->SetBinContent(3, 0.6);
    hRatioDeltaRSmoothed[3][7][6]->SetBinContent(4, 0.66);
    
    // Centrality 50-90, trigger pT 12-16, associated pT 12-16
    hRatioDeltaRSmoothed[3][7][7]->SetBinContent(1, 0.45);
    
    // Centrality 50-90, trigger pT 16-20, associated pT 2-3
    hRatioDeltaRSmoothed[3][8][2]->SetBinContent(1, 0.18);
    
    // Centrality 50-90, trigger pT 16-20, associated pT 4-6
    hRatioDeltaRSmoothed[3][8][4]->SetBinContent(1, 0.34);
    hRatioDeltaRSmoothed[3][8][4]->SetBinContent(2, 0.6);
    hRatioDeltaRSmoothed[3][8][4]->SetBinContent(3, 0.63);
    hRatioDeltaRSmoothed[3][8][4]->SetBinContent(4, 0.66);
    hRatioDeltaRSmoothed[3][8][4]->SetBinContent(5, 0.68);
    hRatioDeltaRSmoothed[3][8][4]->SetBinContent(6, 0.7);
    
    // Centrality 50-90, trigger pT 16-20, associated pT 6-8
    hRatioDeltaRSmoothed[3][8][5]->SetBinContent(1, 0.38);
    hRatioDeltaRSmoothed[3][8][5]->SetBinContent(2, 0.6);
    hRatioDeltaRSmoothed[3][8][5]->SetBinContent(3, 0.63);
    hRatioDeltaRSmoothed[3][8][5]->SetBinContent(4, 0.66);
    
    // Centrality 50-90, trigger pT 16-20, associated pT 8-12
    hRatioDeltaRSmoothed[3][8][6]->SetBinContent(1, 0.28);
    
    // Centrality 50-90, trigger pT 16-20, associated pT 12-16
    hRatioDeltaRSmoothed[3][8][7]->SetBinContent(1, 0.4);
    
    // Centrality 50-90, trigger pT 16-20, associated pT 16-20
    hRatioDeltaRSmoothed[3][8][8]->SetBinContent(1, 0.4);
    hRatioDeltaRSmoothed[3][8][8]->SetBinContent(2, 0.43);
    hRatioDeltaRSmoothed[3][8][8]->SetBinContent(3, 0.46);
    
    // Centrality 50-90, trigger pT 20-300, associated pT 4-6
    hRatioDeltaRSmoothed[3][9][4]->SetBinContent(1, 0.44);
    hRatioDeltaRSmoothed[3][9][4]->SetBinContent(2, 0.63);
    
    // Centrality 50-90, trigger pT 20-300, associated pT 8-12
    hRatioDeltaRSmoothed[3][9][6]->SetBinContent(1, 0.49);
    hRatioDeltaRSmoothed[3][9][6]->SetBinContent(2, 0.64);
    hRatioDeltaRSmoothed[3][9][6]->SetBinContent(3, 0.66);
    hRatioDeltaRSmoothed[3][9][6]->SetBinContent(4, 0.68);
    
    // Centrality 50-90, trigger pT 20-300, associated pT 12-16
    hRatioDeltaRSmoothed[3][9][7]->SetBinContent(1, 0.5);
    hRatioDeltaRSmoothed[3][9][7]->SetBinContent(2, 0.58);
    hRatioDeltaRSmoothed[3][9][7]->SetBinContent(3, 0.62);
    hRatioDeltaRSmoothed[3][9][7]->SetBinContent(4, 0.65);
    hRatioDeltaRSmoothed[3][9][7]->SetBinContent(5, 0.68);
    hRatioDeltaRSmoothed[3][9][7]->SetBinContent(6, 0.7);
    
    // Centrality 50-90, trigger pT 20-300, associated pT 20-300
    hRatioDeltaRSmoothed[3][9][9]->SetBinContent(1, 0.45);
    hRatioDeltaRSmoothed[3][9][9]->SetBinContent(2, 0.5);
    hRatioDeltaRSmoothed[3][9][9]->SetBinContent(3, 0.55);
    hRatioDeltaRSmoothed[3][9][9]->SetBinContent(4, 0.6);
    
  }
  
  // ===============================================
  //        Draw the trigger turn-on curves
  // ===============================================
  
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  //TLine *oneLine = new TLine(0,1,500,1);
  //oneLine->SetLineStyle(2);
  //oneLine->SetLineColor(kBlack);
  //TLine *oneTwentyLine = new TLine(120,0,120,1.3);
  //oneTwentyLine->SetLineStyle(2);
  //oneTwentyLine->SetLineColor(kRed);
  TString centralityString;
  TString compactCentralityString;
  TString triggerPtString;
  TString compactTriggerPtString;
  TString associatedPtString;
  TString compactAssociatedPtString;
  
  // Selection of colors and styles
  int color[] = {kBlack,kRed,kBlue,kGreen+4};
  int markerStyle[] = {kOpenSquare, kOpenCircle, kOpenDiamond, kOpenCross};
  
  // Draw the track pair efficiency plots
  for(int iCentrality = firstDrawCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    if(isPbPbData){
      centralityString = Form("Centrality: %.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C%.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "pp";
      compactCentralityString = "_pp";
    }
    
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      
      triggerPtString = Form("%.1f < p_{T,t} < %.1f GeV", systemCard->GetLowBinBorderTrackPairPt(iTriggerPt), systemCard->GetHighBinBorderTrackPairPt(iTriggerPt));
      compactTriggerPtString = Form("_T%.0f-%.0f", systemCard->GetLowBinBorderTrackPairPt(iTriggerPt), systemCard->GetHighBinBorderTrackPairPt(iTriggerPt));
      
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        
        associatedPtString = Form("%.1f < p_{T,a} < %.1f GeV", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
        compactAssociatedPtString = Form("_T%.0f-%.0f", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
        
        // Create a legend for the figure
        legend = new TLegend(0.2,0.73,0.5,0.88);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        // Draw the track pair efficiency
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt]->GetYaxis()->SetRangeUser(0,2);
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt]->GetXaxis()->SetRangeUser(0,0.06);
        drawer->DrawHistogram(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt], "#DeltaR", "Track pair efficiency", " ");
        
        // Draw the smoothed histogram to the same figure
        hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt]->SetLineColor(kRed);
        hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt]->Draw("same");
        
        // Add information to legend
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, triggerPtString.Data(), "");
        legend->AddEntry((TObject*)0, associatedPtString.Data(), "");
        
        // Draw the legend
        legend->Draw();
        
        // Draw a lines to one and 120 GeV
        //oneLine->Draw("same");
        //oneTwentyLine->Draw("same");
        
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/trackPairEfficiency%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactTriggerPtString.Data(), compactAssociatedPtString.Data(), figureFormat));
        }
        
      } // Associated particle pT loop
    } // Trigger particle pT loop
  } // Centrality loop
  
  // Option to write the smoothed histograms to a file to be used as correction
  if(writeSmoothedHistograms){
    TFile* outputFile = TFile::Open("trackPairEfficiencyCorrectionTable_PbPb2018.root","RECREATE");
    
    for(int iCentrality = firstDrawCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt]->Write();
        } // Associated pT loo
      } // Trigger pT loop
    } // Centrality loop
  }
  
}
