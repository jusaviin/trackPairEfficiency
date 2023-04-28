#include "TrackPairEfficiencyHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TrackPairEfficiencyCard.h"
#include "JDrawer.h"
#include "../src/TrackPairEfficiencyHistograms.h"

/*
 * Macro for determining track pair efficiency from reconstructed and generator level simulations
 */
void trackPairEfficiencyPlotter(){
  
  // File containing the track pair distributions used to determine the track pair efficiency
  TString fileName = "data/trackPairEfficiency_triggerFiducialCut_newBins_fixedCentrality_processed_2023-04-10.root";
  // trackPairEfficiencyPp_triggerEta1p6_newBins_processed_2023-04-10.root
  // trackPairEfficiency_triggerFiducialCut_newBins_fixedCentrality_processed_2023-04-10.root
  
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
  const int nAverageEtaBins = systemCard->GetNAverageEtaBins();
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Select the bin range to be drawn
  const int firstDrawnCentralityBin = 0;
  const int lastDrawnCentralityBin = nCentralityBins-1;
  
  const int firstDrawnTriggerPtBin = 0;
  const int lastDrawnTriggerPtBin = nTrackPtBins-1;
  
  const int firstDrawnAssociatedPtBin = 0;
  const int lastDrawnAssociatedPtBin = nTrackPtBins-1;
  
  int firstDrawnAverageEtaBin = nAverageEtaBins;
  int lastDrawnAverageEtaBin = nAverageEtaBins;
  
  // For eta comparison drawing, can do manual selection of which bins to draw
  const bool useSelectedBins = false;
  const int nSelectedAverageEtaBins = 3;
  const int selectedAverageEtaBins[nSelectedAverageEtaBins] = {0, nAverageEtaBins-1, nAverageEtaBins};
  if(useSelectedBins){
    firstDrawnAverageEtaBin = selectedAverageEtaBins[0];
    lastDrawnAverageEtaBin = selectedAverageEtaBins[nSelectedAverageEtaBins-1];
  }
  
  // Histogram drawing
  const bool drawDistributionComparison = false; // Draw the comparison of the raw distributions
  const bool drawSmoothedComparison = false;     // Draw the figures comparing smoothed and non-smoothed distributions in each bin
  const bool drawPtBinComparison = false;        // Draw a comparison of different pT bins to find pT trends
  const bool drawEtaRegionComparison = false;    // Draw the figures comparing the different average eta selections for the same pT and centrality bins
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Histogram saving
  const bool writeSmoothedHistograms = true; // Write flag for smoothed histograms
  
  // Create and setup a new histogram managers to project and handle the histograms
  TrackPairEfficiencyHistogramManager *histograms;
  histograms = new TrackPairEfficiencyHistogramManager(inputFile);
  histograms->SetCentralityBinRange(firstDrawnCentralityBin, lastDrawnCentralityBin);
  histograms->SetTrackPairPtBinRange(TMath::Min(firstDrawnTriggerPtBin, firstDrawnAssociatedPtBin), TMath::Max(lastDrawnTriggerPtBin, lastDrawnAssociatedPtBin));
  histograms->SetAverageEtaBinRange(firstDrawnAverageEtaBin, lastDrawnAverageEtaBin);
  histograms->SetLoadAllTrackPairs(true,true);
  histograms->LoadProcessedHistograms();
  
  // Jet pT histograms to calculate trigger turn on curve
  TH1D* hTrackDeltaR[nCentralityBins][nTrackPtBins][nTrackPtBins][nAverageEtaBins+1][TrackPairEfficiencyHistograms::knDataLevels];   // DeltaR distribution for track pairs
  TH1D* hRatioDeltaR[nCentralityBins][nTrackPtBins][nTrackPtBins][nAverageEtaBins+1]; // Reconstructed to generator level ratio of track pair DeltaR distributions
  TH1D* hRatioDeltaRSmoothed[nCentralityBins][nTrackPtBins][nTrackPtBins][nAverageEtaBins+1]; // Smoothed reconstructed to generator level ratio of track pair DeltaR distributions
  
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < nTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < nTrackPtBins; iAssociatedPt++){
        for(int iAverageEta = 0; iAverageEta < nAverageEtaBins+1; iAverageEta++){
          for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
            hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel] = NULL;
          } // Data level loop
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta] = NULL;
          hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta] = NULL;
        }
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Read the histograms from the file
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
          for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
            hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel] = histograms->GetHistogramTrackPairDeltaR(iCentrality, iTriggerPt, iAssociatedPt, iAverageEta, iDataLevel);
          } // Data level loop
        } // Average pair eta loop
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Calculate the generator level to reconstructed ratio to determine track pair efficiency
  double genLevelIntegral;
  double integralToMatch;
  TString histogramNamer;
  int lowBinLimit = hTrackDeltaR[firstDrawnCentralityBin][firstDrawnTriggerPtBin][firstDrawnAssociatedPtBin][firstDrawnAverageEtaBin][0]->FindBin(0.1);
  int highBinLimit = hTrackDeltaR[firstDrawnCentralityBin][firstDrawnTriggerPtBin][firstDrawnAssociatedPtBin][firstDrawnAverageEtaBin][0]->FindBin(0.4);
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
          
          // First, normalize the distributions such that the integrals match in the region 0.1-0.4
          integralToMatch = hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kReconstructed]->Integral(lowBinLimit, highBinLimit, "width");
          genLevelIntegral = hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kGeneratorLevel]->Integral(lowBinLimit, highBinLimit, "width");
          hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kGeneratorLevel]->Scale(integralToMatch/genLevelIntegral);
          
          // Then, calculate the ratio of the histograms
          if(iAverageEta == nAverageEtaBins){
            histogramNamer = Form("trackPairEfficiencyCorrection_C%dT%dA%d", iCentrality, iTriggerPt, iAssociatedPt);
          } else {
            histogramNamer = Form("trackPairEfficiencyCorrection_C%dT%dA%dE%d", iCentrality, iTriggerPt, iAssociatedPt, iAverageEta);
          }
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta] = (TH1D*) hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kReconstructed]->Clone(histogramNamer.Data());
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Divide(hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kGeneratorLevel]);
          
          // Finally, smooth the distribution to suppress fluctuations
          if(iAverageEta == nAverageEtaBins){
            histogramNamer = Form("smoothedTrackPairEfficiencyCorrection_C%dT%dA%d", iCentrality, iTriggerPt, iAssociatedPt);
          } else {
            histogramNamer = Form("smoothedTrackPairEfficiencyCorrection_C%dT%dA%dE%d", iCentrality, iTriggerPt, iAssociatedPt, iAverageEta);
          }
          hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta] = (TH1D*) hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Clone(histogramNamer.Data());
          hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Smooth();
        } // Average pair eta loop
      } // Associated pT loop
    } // Trigger pT loop
  } // Centrality loop
  
  // Manual extra smoothing. The smoothing does generally good job, but some histogram require a manual extra touch to get the smoothing satisfactory
  // The values here fix the smoothing for file trackPairEfficiency_triggerFiducialCut_processed_2023-03-16.root
  
  /*if(firstDrawnCentralityBin == 0){
    
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
  
  if(firstDrawnCentralityBin == 1 || (firstDrawnCentralityBin < 1 && lastDrawnCentralityBin >= 1)){
    
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
  
  if(firstDrawnCentralityBin == 2 || (firstDrawnCentralityBin < 2 && lastDrawnCentralityBin >= 2)){
    
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
  
  if(firstDrawnCentralityBin == 3 || (firstDrawnCentralityBin < 3 && lastDrawnCentralityBin >= 3)){
   
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
    
  }*/
  
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
  TString averageEtaString;
  TString compactAverageEtaString;
  double legendY1;
  
  // Selection of colors and styles
  int color[] = {kBlack,kRed,kBlue,kGreen+2,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  int markerStyle[] = {kOpenSquare, kOpenCircle, kOpenDiamond, kOpenCross, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenCrossX, kOpenDoubleDiamond, kOpenFourTrianglesPlus};
  
  // Draw the track pair efficiency histogram comparison with the smoothed histograms
  if(drawDistributionComparison){
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
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
          
          for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
            
            if(iAverageEta < nAverageEtaBins){
              averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
              compactAverageEtaString = Form("_E%.1f-%.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
              compactAverageEtaString.ReplaceAll(".","p");
              legendY1 = 0.58;
            } else {
              averageEtaString = "";
              compactAverageEtaString = "";
              legendY1 = 0.63;
            }
            
            // Create a legend for the figure
            legend = new TLegend(0.2,legendY1,0.5,0.88);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            
            // Draw the track pair efficiency
            drawer->DrawHistogram(hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kReconstructed], "#DeltaR", "Track pair yield", " ");
            
            // Draw the smoothed histogram to the same figure
            hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kGeneratorLevel]->SetLineColor(kRed);
            hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kGeneratorLevel]->Draw("same");
            
            // Add information to legend
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, triggerPtString.Data(), "");
            legend->AddEntry((TObject*)0, associatedPtString.Data(), "");
            if(iAverageEta < nAverageEtaBins) legend->AddEntry((TObject*)0, averageEtaString.Data(), "");
            legend->AddEntry(hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kReconstructed], "Reconstructed tracks", "l");
            legend->AddEntry(hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kGeneratorLevel], "Generator level particles", "l");
            
            // Draw the legend
            legend->Draw();
            
            // Draw a lines to one and 120 GeV
            //oneLine->Draw("same");
            //oneTwentyLine->Draw("same");
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/trackPairDistribution%s%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactTriggerPtString.Data(), compactAssociatedPtString.Data(), compactAverageEtaString.Data(), figureFormat));
            }
            
          } // Average pair eta loop
        } // Associated particle pT loop
      } // Trigger particle pT loop
    } // Centrality loop
  } // Draw the smoothed comparison
  
  // Draw the track pair efficiency histogram comparison with the smoothed histograms
  if(drawSmoothedComparison){
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
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
          compactAssociatedPtString = Form("_A%.0f-%.0f", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
          
          for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
            
            if(iAverageEta < nAverageEtaBins){
              averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
              compactAverageEtaString = Form("_E%.1f-%.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
              compactAverageEtaString.ReplaceAll(".","p");
              legendY1 = 0.57;
            } else {
              averageEtaString = "";
              compactAverageEtaString = "";
              legendY1 = 0.62;
            }
            
            // Create a legend for the figure
            legend = new TLegend(0.2,legendY1,0.5,0.92);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            
            // Draw the track pair efficiency
            hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->GetYaxis()->SetRangeUser(0,2);
            hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->GetXaxis()->SetRangeUser(0,0.06);
            drawer->DrawHistogram(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta], "#DeltaR", "Track pair efficiency", " ");
            
            // Draw the smoothed histogram to the same figure
            hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->SetLineColor(kRed);
            hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Draw("same");
            
            // Add information to legend
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, triggerPtString.Data(), "");
            legend->AddEntry((TObject*)0, associatedPtString.Data(), "");
            if(iAverageEta < nAverageEtaBins) legend->AddEntry((TObject*)0, averageEtaString.Data(), "");
            legend->AddEntry(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta], "Yield ratio", "l");
            legend->AddEntry(hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta], "Smoothed ratio", "l");
            
            // Draw the legend
            legend->Draw();
            
            // Draw a lines to one and 120 GeV
            //oneLine->Draw("same");
            //oneTwentyLine->Draw("same");
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/trackPairEfficiency%s%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactTriggerPtString.Data(), compactAssociatedPtString.Data(), compactAverageEtaString.Data(), figureFormat));
            }
            
          } // Average pair eta loop
        } // Associated particle pT loop
      } // Trigger particle pT loop
    } // Centrality loop
  } // Draw the smoothed comparison
  
  // Draw the track pair efficiency histogram comparison with the smoothed histograms
  if(drawPtBinComparison){
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      if(isPbPbData){
        centralityString = Form("Centrality: %.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", systemCard->GetLowBinBorderCentrality(iCentrality), systemCard->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "pp";
        compactCentralityString = "_pp";
      }
      
      for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
        
        if(iAverageEta < nAverageEtaBins){
          averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
          compactAverageEtaString = Form("_E%.1f-%.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
          compactAverageEtaString.ReplaceAll(".","p");
          legendY1 = 0.57;
        } else {
          averageEtaString = "";
          compactAverageEtaString = "";
          legendY1 = 0.62;
        }
        
        // Create a legend for the figure
        legend = new TLegend(0.2,legendY1,0.5,0.92);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        // Draw the track pair efficiency in the first pT bin
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->GetYaxis()->SetRangeUser(0,2);
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->GetXaxis()->SetRangeUser(0,0.06);
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->SetLineColor(color[0]);
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->SetMarkerColor(color[0]);
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->SetMarkerStyle(markerStyle[0]);
        drawer->DrawHistogram(hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta], "#DeltaR", "Track pair efficiency", " ");
        
        // Draw all the other track pT bins to the same figure
        for(int iTrackPt = firstDrawnTriggerPtBin+1; iTrackPt <= lastDrawnTriggerPtBin; iTrackPt++){
          
          hRatioDeltaR[iCentrality][iTrackPt][iTrackPt][iAverageEta]->SetLineColor(color[iTrackPt-firstDrawnTriggerPtBin]);
          hRatioDeltaR[iCentrality][iTrackPt][iTrackPt][iAverageEta]->SetMarkerColor(color[iTrackPt-firstDrawnTriggerPtBin]);
          hRatioDeltaR[iCentrality][iTrackPt][iTrackPt][iAverageEta]->SetMarkerStyle(markerStyle[iTrackPt-firstDrawnTriggerPtBin]);
          hRatioDeltaR[iCentrality][iTrackPt][iTrackPt][iAverageEta]->Draw("same");
          
        }
        
        // Add information to the legend
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        if(iAverageEta < nAverageEtaBins) legend->AddEntry((TObject*)0, averageEtaString.Data(), "");
        
        for(int iTrackPt = firstDrawnTriggerPtBin; iTrackPt <= lastDrawnTriggerPtBin; iTrackPt++){
          
          triggerPtString = Form("%.1f < p_{T,t&a} < %.1f GeV", systemCard->GetLowBinBorderTrackPairPt(iTrackPt), systemCard->GetHighBinBorderTrackPairPt(iTrackPt));
          legend->AddEntry(hRatioDeltaR[iCentrality][iTrackPt][iTrackPt][iAverageEta], triggerPtString.Data(), "p");
          
        }
        
        // Draw the legend
        legend->Draw();
        
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/trackPairEfficiencyPtTrend%s%s%s.%s", saveComment, compactCentralityString.Data(), compactAverageEtaString.Data(), figureFormat));
        }
        
      } // Average pair eta loop
    } // Centrality loop
  } // Draw the smoothed comparison
  
  // Draw the track pair efficiency histogram average eta region comparison
  if(drawEtaRegionComparison){
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
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
          
          if(useSelectedBins){
            legendY1 = 0.88 - 0.12 - nSelectedAverageEtaBins*0.04;
          } else {
            legendY1 = 0.88 - 0.12 - (lastDrawnAverageEtaBin - firstDrawnAverageEtaBin + 1)*0.04;
          }
          
          // Create a legend for the figure
          legend = new TLegend(0.2,legendY1,0.5,0.88);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          // Draw the first average pair eta bin
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->SetLineColor(color[0]);
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->SetMarkerColor(color[0]);
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->SetMarkerStyle(markerStyle[0]);
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->GetYaxis()->SetRangeUser(0,2);
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->GetXaxis()->SetRangeUser(0,0.06);
          drawer->DrawHistogram(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin], "#DeltaR", "Track pair efficiency", " ");
          
          // Draw the other average eta bins to the same figure
          if(useSelectedBins){
            for(int iSelected = 1; iSelected < nSelectedAverageEtaBins; iSelected++){
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][selectedAverageEtaBins[iSelected]]->SetLineColor(color[iSelected]);
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][selectedAverageEtaBins[iSelected]]->SetMarkerColor(color[iSelected]);
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][selectedAverageEtaBins[iSelected]]->SetMarkerStyle(markerStyle[iSelected]);
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][selectedAverageEtaBins[iSelected]]->Draw("same");
            }
          } else {
            for(int iAverageEta = firstDrawnAverageEtaBin+1; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->SetLineColor(color[iAverageEta-firstDrawnAverageEtaBin]);
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->SetMarkerColor(color[iAverageEta-firstDrawnAverageEtaBin]);
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->SetMarkerStyle(markerStyle[iAverageEta-firstDrawnAverageEtaBin]);
              hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Draw("same");
            }
          }
          
          // Add information to legend
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, triggerPtString.Data(), "");
          legend->AddEntry((TObject*)0, associatedPtString.Data(), "");
          
          if(useSelectedBins){
            for(int iSelected = 0; iSelected < nSelectedAverageEtaBins; iSelected++){
              if(selectedAverageEtaBins[iSelected] < nAverageEtaBins){
                averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(selectedAverageEtaBins[iSelected]), systemCard->GetHighBinBorderAverageEta(selectedAverageEtaBins[iSelected]));
              } else {
                averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(0), systemCard->GetHighBinBorderAverageEta(nAverageEtaBins-1));
              }
              legend->AddEntry(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][selectedAverageEtaBins[iSelected]], averageEtaString.Data(), "p");
            }
          } else {
            for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
              if(iAverageEta < nAverageEtaBins){
                averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(iAverageEta), systemCard->GetHighBinBorderAverageEta(iAverageEta));
              } else {
                averageEtaString = Form("%.1f < #LT#eta_{pair}#GT < %.1f", systemCard->GetLowBinBorderAverageEta(0), systemCard->GetHighBinBorderAverageEta(nAverageEtaBins-1));
              }
              legend->AddEntry(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta], averageEtaString.Data(), "p");
            }
          }
          
          // Draw the legend
          legend->Draw();
          
          // Draw a lines to one and 120 GeV
          //oneLine->Draw("same");
          //oneTwentyLine->Draw("same");
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/trackPairEfficiencyEtaComparison%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactTriggerPtString.Data(), compactAssociatedPtString.Data(), figureFormat));
          }
          
        } // Associated particle pT loop
      } // Trigger particle pT loop
    } // Centrality loop
  } // Draw the smoothed comparison
  
  // Option to write the smoothed histograms to a file to be used as correction
  if(writeSmoothedHistograms){
    TFile* outputFile = TFile::Open("trackPairEfficiencyCorrectionTableNewBinning_PbPb2018_2023-04-17.root","RECREATE");
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
            hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Write();
            hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Write();
          }
        } // Associated pT loo
      } // Trigger pT loop
    } // Centrality loop
    
    // Write the card to the file
    systemCard->Write(outputFile);
  }
  
}
