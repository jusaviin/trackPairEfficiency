#include "TrackPairEfficiencyHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TrackPairEfficiencyCard.h"
#include "JDrawer.h"
#include "../src/TrackPairEfficiencyHistograms.h"

/*
 * Macro for determining track pair efficiency from reconstructed and generator level simulations
 */
void trackPairEfficiencyPlotter(){
  
  // File containing the track pair distributions used to determine the track pair efficiency
  TString fileName = "data/trackPairEfficiencyPp_triggerEta1p6_32DeltaRBins_processed_2023-07-11.root";
  // trackPairEfficiencyPp_triggerEta1p6_newBins_processed_2023-04-10.root
  // trackPairEfficiency_triggerFiducialCut_newBins_fixedCentrality_processed_2023-04-10.root
  // trackPairEfficiencyPp_triggerEta1p6_32DeltaRBins_processed_2023-05-16.root
  // trackPairEfficiencyPp_triggerEta1p6_32DeltaRBins_combineHighestTrackPtBins_processed_2023-05-16.root
  // trackPairEfficiency_triggerFiducialCut_32DeltaRBins_processed_2023-05-16.root
  // trackPairEfficiency_triggerFiducialCut_32DeltaRBins_combineHighestTrackPtBins_processed_2023-05-16.root
  // trackPairEfficiency_triggerFiducialCut_32DeltaRBins_processed_2023-07-11.root
  // trackPairEfficiencyPp_triggerEta1p6_32DeltaRBins_processed_2023-07-11.root
  
  // Open the file and check that it exists
  TFile* inputFile = TFile::Open(fileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Check if we are using PbPb or pp data
  TrackPairEfficiencyCard* systemCard = new TrackPairEfficiencyCard(inputFile);
  TString collisionSystem = systemCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  const int nCentralityBins = (isPbPbData) ? systemCard->GetNCentralityBins() : 1;
  const int nTrackPtBins = systemCard->GetNTrackPairPtBins();
  const int nJetPtBins = systemCard->GetNJetPtBins();
  const int nAverageEtaBins = systemCard->GetNAverageEtaBins();
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Select the bin range to be drawn
  const int firstDrawnCentralityBin = 0;
  const int lastDrawnCentralityBin = nCentralityBins-1;
  
  const int firstDrawnTriggerPtBin = 5;
  const int lastDrawnTriggerPtBin = 5;
  
  const int firstDrawnAssociatedPtBin = 3;
  const int lastDrawnAssociatedPtBin = 3;

  const int firstDrawnJetPtBin = 0;
  const int lastDrawnJetPtBin = nJetPtBins-1;
  
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
  const bool drawDistributionComparison = false;       // Draw the comparison of the raw distributions
  const bool drawSmoothedComparison = true;           // Draw the figures comparing smoothed and non-smoothed distributions in each bin
  const bool drawPtBinComparison = false;              // Draw a comparison of different pT bins to find pT trends
  const bool drawEtaRegionComparison = false;          // Draw the figures comparing the different average eta selections for the same pT and centrality bins
  const bool drawJetPtComparison = false;              // Draw the comparison of different jet pT selections to see if there is a trend
  const bool drawSmoothedComparisonCloseToJets = false; // Draw the figures comparing smoothed and non-smoothed distributions close to jets
  
  // Axis zoom for the drawn histgrams
  const double minXzoom = 0.006;
  const double maxXzoom = 0.39;
  const double minRatioZoom = 0.4;
  const double maxRatioZoom = 1.6;

  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_pythia";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Histogram saving
  const bool writeCorrectionToFile = false; // Write flag for the correction histograms
  
  // Create and setup a new histogram managers to project and handle the histograms
  TrackPairEfficiencyHistogramManager* histograms;
  histograms = new TrackPairEfficiencyHistogramManager(inputFile);
  histograms->SetCentralityBinRange(firstDrawnCentralityBin, lastDrawnCentralityBin);
  histograms->SetTrackPairPtBinRange(TMath::Min(firstDrawnTriggerPtBin, firstDrawnAssociatedPtBin), TMath::Max(lastDrawnTriggerPtBin, lastDrawnAssociatedPtBin));
  histograms->SetJetPtBinRange(firstDrawnJetPtBin, lastDrawnJetPtBin);
  histograms->SetAverageEtaBinRange(firstDrawnAverageEtaBin, lastDrawnAverageEtaBin);
  histograms->SetLoadAllTrackPairs(true,true);
  histograms->SetLoadAllTrackPairsCloseToJets(true,true);
  histograms->LoadProcessedHistograms();
  
  // Track pair histograms to calculate the pair efficiency
  TH1D* hTrackDeltaR[nCentralityBins][nTrackPtBins][nTrackPtBins][nAverageEtaBins+1][TrackPairEfficiencyHistograms::knDataLevels];   // DeltaR distribution for track pairs
  TH1D* hTrackDeltaRCloseToJet[TrackPairEfficiencyHistograms::knDataLevels][TrackPairEfficiencyHistograms::knDataLevels][nCentralityBins][nTrackPtBins][nTrackPtBins][nJetPtBins+1];   // DeltaR distribution for track pairs close to jets
  TH1D* hRatioDeltaR[nCentralityBins][nTrackPtBins][nTrackPtBins][nAverageEtaBins+1]; // Reconstructed to generator level ratio of track pair DeltaR distributions
  TH1D* hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::knDataLevels][nCentralityBins][nTrackPtBins][nTrackPtBins][nJetPtBins+1]; // Reconstructed to generator level ratio of track pair DeltaR distributions close to jets
  TH1D* hRatioDeltaRSmoothed[nCentralityBins][nTrackPtBins][nTrackPtBins][nAverageEtaBins+1]; // Smoothed reconstructed to generator level ratio of track pair DeltaR distributions
  TH1D* hRatioDeltaRCloseToJetSmoothed[TrackPairEfficiencyHistograms::knDataLevels][nCentralityBins][nTrackPtBins][nTrackPtBins][nJetPtBins+1]; // Reconstructed to generator level ratio of track pair DeltaR distributions close to jets
  
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < nTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < nTrackPtBins; iAssociatedPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBins+1; iJetPt++){
          for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
            for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
              hTrackDeltaRCloseToJet[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = NULL;
            } // Track data level loop
            hRatioDeltaRCloseToJet[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = NULL;
            hRatioDeltaRCloseToJetSmoothed[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = NULL;
          } // Data level loop for jets
        } // Jet pT loop

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
        for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
          for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
            for(int iJetPt = firstDrawnJetPtBin; iJetPt <= lastDrawnJetPtBin; iJetPt++){
              hTrackDeltaRCloseToJet[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = histograms->GetHistogramTrackPairDeltaRCloseToJets(iDataLevelJets, iDataLevelTracks, iCentrality, iTriggerPt, iAssociatedPt, iJetPt);
            } // Jet pT loop
          } // Data level loop for jets
        } // Data level loop for tracks
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

        // Histograms in jet pT bins
        for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
          for(int iJetPt = firstDrawnJetPtBin; iJetPt <= lastDrawnJetPtBin; iJetPt++){

            // First, normalize the distributions such that the integrals match in the region 0.1-0.4
          integralToMatch = hTrackDeltaRCloseToJet[iDataLevelJets][TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Integral(lowBinLimit, highBinLimit, "width");
          genLevelIntegral = hTrackDeltaRCloseToJet[iDataLevelJets][TrackPairEfficiencyHistograms::kGeneratorLevel][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Integral(lowBinLimit, highBinLimit, "width");
          hTrackDeltaRCloseToJet[iDataLevelJets][TrackPairEfficiencyHistograms::kGeneratorLevel][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Scale(integralToMatch/genLevelIntegral);

          // Then, calculate the ratio of the histograms
          if(iJetPt == nJetPtBins){
            histogramNamer = Form("trackPairEfficiencyCorrectionCloseTo%sJet_C%dT%dA%d", histograms->GetDataLevelName(iDataLevelJets), iCentrality, iTriggerPt, iAssociatedPt);
          } else {
            histogramNamer = Form("trackPairEfficiencyCorrectionCloseTo%sJet_C%dT%dA%dJ%d", histograms->GetDataLevelName(iDataLevelJets), iCentrality, iTriggerPt, iAssociatedPt, iJetPt);
          }

          hRatioDeltaRCloseToJet[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = (TH1D*) hTrackDeltaRCloseToJet[iDataLevelJets][TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Clone(histogramNamer.Data());
          hRatioDeltaRCloseToJet[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Divide(hTrackDeltaRCloseToJet[iDataLevelJets][TrackPairEfficiencyHistograms::kGeneratorLevel][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]);

          // Use the smoothing also for close-to-jet distributions to suppress fluctuations in the correction
          // Then, calculate the ratio of the histograms
          if(iJetPt == nJetPtBins){
            histogramNamer = Form("smoothedTrackPairEfficiencyCorrectionCloseTo%sJet_C%dT%dA%d", histograms->GetDataLevelName(iDataLevelJets), iCentrality, iTriggerPt, iAssociatedPt);
          } else {
            histogramNamer = Form("smoothedTrackPairEfficiencyCorrectionCloseTo%sJet_C%dT%dA%dJ%d", histograms->GetDataLevelName(iDataLevelJets), iCentrality, iTriggerPt, iAssociatedPt, iJetPt);
          }
          hRatioDeltaRCloseToJetSmoothed[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = (TH1D*) hRatioDeltaRCloseToJet[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Clone(histogramNamer.Data());
          hRatioDeltaRCloseToJetSmoothed[iDataLevelJets][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Smooth();

          } // Jet pT loop
        } // Data level loop for jets
      } // Associated pT loop
    } // Trigger pT loop
  } // Centrality loop
  
  // ===============================================
  //        Draw the trigger turn-on curves
  // ===============================================
  
  JDrawer* drawer = new JDrawer();
  TLegend* legend;
  TLine* oneLine = new TLine(minXzoom,1,maxXzoom,1);
  oneLine->SetLineStyle(2);
  oneLine->SetLineColor(kBlack);

  TString centralityString;
  TString compactCentralityString;
  TString triggerPtString;
  TString compactTriggerPtString;
  TString associatedPtString;
  TString compactAssociatedPtString;
  TString averageEtaString;
  TString compactAverageEtaString;
  TString jetPtString;
  TString compactJetPtString;
  double legendY1;
  
  // Selection of colors and styles
  int color[] = {kBlack,kRed,kBlue,kGreen+2,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  int markerStyle[] = {kOpenSquare, kOpenCircle, kOpenDiamond, kOpenCross, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenCrossX, kOpenDoubleDiamond, kOpenFourTrianglesPlus};

  // Use logarithmic x-axis for the figures
  drawer->SetLogX(true);
  
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
        
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= TMath::Min(iTriggerPt, lastDrawnAssociatedPtBin); iAssociatedPt++){
          
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
            hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][TrackPairEfficiencyHistograms::kReconstructed]->GetXaxis()->SetRangeUser(minXzoom,maxXzoom);
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
        
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= TMath::Min(iTriggerPt, lastDrawnAssociatedPtBin); iAssociatedPt++){
          
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
            hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->GetYaxis()->SetRangeUser(minRatioZoom,maxRatioZoom);
            hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->GetXaxis()->SetRangeUser(minXzoom,maxXzoom);
            drawer->DrawHistogram(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta], "#DeltaR", "Reco yield / Gen yield", " ");
            
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
            
            // Draw a lines to one
            oneLine->Draw("same");
            
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
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->GetYaxis()->SetRangeUser(minRatioZoom,maxRatioZoom);
        hRatioDeltaR[iCentrality][firstDrawnTriggerPtBin][firstDrawnTriggerPtBin][iAverageEta]->GetXaxis()->SetRangeUser(minXzoom,maxXzoom);
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
        
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= TMath::Min(iTriggerPt, lastDrawnAssociatedPtBin); iAssociatedPt++){
          
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
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->GetYaxis()->SetRangeUser(minRatioZoom,maxRatioZoom);
          hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][firstDrawnAverageEtaBin]->GetXaxis()->SetRangeUser(minXzoom,maxXzoom);
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
  } // Draw the eta region comparison

  // Draw the track pair efficiency histogram jet pT selection comparison
  if(drawJetPtComparison){
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
        
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= TMath::Min(iTriggerPt, lastDrawnAssociatedPtBin); iAssociatedPt++){
          
          associatedPtString = Form("%.1f < p_{T,a} < %.1f GeV", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
          compactAssociatedPtString = Form("_T%.0f-%.0f", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
          
          legendY1 = 0.88 - 0.12 - (lastDrawnJetPtBin - firstDrawnJetPtBin + 1)*0.04;
          
          // Create a legend for the figure
          legend = new TLegend(0.2,legendY1,0.5,0.88);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          // Draw the first jet pT bin
          hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][firstDrawnJetPtBin]->SetLineColor(color[0]);
          hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][firstDrawnJetPtBin]->SetMarkerColor(color[0]);
          hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][firstDrawnJetPtBin]->SetMarkerStyle(markerStyle[0]);
          hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][firstDrawnJetPtBin]->GetYaxis()->SetRangeUser(minRatioZoom,maxRatioZoom);
          hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][firstDrawnJetPtBin]->GetXaxis()->SetRangeUser(minXzoom,maxXzoom);
          drawer->DrawHistogram(hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][firstDrawnJetPtBin], "#DeltaR", "Track pair efficiency", " ");
          
          // Draw the other jet pT bins to the same figure
          for(int iJetPt = firstDrawnJetPtBin+1; iJetPt <= lastDrawnJetPtBin; iJetPt++){
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->SetLineColor(color[iJetPt-firstDrawnJetPtBin]);
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->SetMarkerColor(color[iJetPt-firstDrawnJetPtBin]);
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->SetMarkerStyle(markerStyle[iJetPt-firstDrawnJetPtBin]);
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Draw("same");
          }
          
          
          // Add information to legend
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, triggerPtString.Data(), "");
          legend->AddEntry((TObject*)0, associatedPtString.Data(), "");
          
          for(int iJetPt = firstDrawnJetPtBin; iJetPt <= lastDrawnJetPtBin; iJetPt++){
            if(iJetPt < nJetPtBins){
                jetPtString = Form("%.0f < jet p_{T} < %.0f", systemCard->GetLowBinBorderJetPt(iJetPt), systemCard->GetHighBinBorderJetPt(iJetPt));
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.1f", systemCard->GetLowBinBorderJetPt(0), systemCard->GetHighBinBorderJetPt(nJetPtBins-1));
            }
            legend->AddEntry(hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt], jetPtString.Data(), "p");
          }
          
          // Draw the legend
          legend->Draw();
          
          // Draw a lines to one and 120 GeV
          //oneLine->Draw("same");
          //oneTwentyLine->Draw("same");
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/trackPairEfficiencyJetPtComparison%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactTriggerPtString.Data(), compactAssociatedPtString.Data(), figureFormat));
          }
          
        } // Associated particle pT loop
      } // Trigger particle pT loop
    } // Centrality loop
  } // Draw the eta region comparison

  // Draw the track pair efficiency histogram comparison with the smoothed histograms
  if(drawSmoothedComparisonCloseToJets){
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
        
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= TMath::Min(iTriggerPt, lastDrawnAssociatedPtBin); iAssociatedPt++){
          
          associatedPtString = Form("%.1f < p_{T,a} < %.1f GeV", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
          compactAssociatedPtString = Form("_A%.0f-%.0f", systemCard->GetLowBinBorderTrackPairPt(iAssociatedPt), systemCard->GetHighBinBorderTrackPairPt(iAssociatedPt));
          
          for(int iJetPt = firstDrawnJetPtBin; iJetPt <= lastDrawnJetPtBin; iJetPt++){
            
            if(iJetPt < nJetPtBins){
              jetPtString = Form("%.0f < jet p_{T} < %.0f", systemCard->GetLowBinBorderJetPt(iJetPt), systemCard->GetHighBinBorderJetPt(iJetPt));
              compactJetPtString = Form("_J%.0f-%.0f", systemCard->GetLowBinBorderJetPt(iJetPt), systemCard->GetHighBinBorderJetPt(iJetPt));
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.1f", systemCard->GetLowBinBorderJetPt(0), systemCard->GetHighBinBorderJetPt(nJetPtBins-1));
              compactJetPtString = Form("_J%.0f-%.0f", systemCard->GetLowBinBorderJetPt(0), systemCard->GetHighBinBorderJetPt(nJetPtBins-1));
            }
            
            // Create a legend for the figure
            legend = new TLegend(0.2,0.58,0.5,0.92);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            
            // Draw the track pair efficiency
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->GetYaxis()->SetRangeUser(minRatioZoom,maxRatioZoom);
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->GetXaxis()->SetRangeUser(minXzoom,maxXzoom);
            drawer->DrawHistogram(hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt], "#DeltaR", "Track pair efficiency", " ");
            
            // Draw the smoothed histogram to the same figure
            hRatioDeltaRCloseToJetSmoothed[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->SetLineColor(kRed);
            hRatioDeltaRCloseToJetSmoothed[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Draw("same");
            
            // Add information to legend
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, jetPtString.Data(), "");
            legend->AddEntry((TObject*)0, triggerPtString.Data(), "");
            legend->AddEntry((TObject*)0, associatedPtString.Data(), "");
            legend->AddEntry(hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt], "Yield ratio", "l");
            legend->AddEntry(hRatioDeltaRCloseToJetSmoothed[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt], "Smoothed ratio", "l");
            
            // Draw the legend
            legend->Draw();
            
            // Draw a lines to one
            oneLine->Draw("same");
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/trackPairEfficiencyCloseToJets%s%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactTriggerPtString.Data(), compactAssociatedPtString.Data(), compactJetPtString.Data(), figureFormat));
            }
            
          } // Jet pT loop
        } // Associated particle pT loop
      } // Trigger particle pT loop
    } // Centrality loop
  } // Draw the smoothed comparison
  
  // Option to write the smoothed histograms to a file to be used as correction
  if(writeCorrectionToFile){
    TFile* outputFile = TFile::Open("trackPairEfficiencyCorrection_PbPb2018_32DeltaRBins_noSmoothing_2023-07-13.root","RECREATE");
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
        for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          for(int iAverageEta = firstDrawnAverageEtaBin; iAverageEta <= lastDrawnAverageEtaBin; iAverageEta++){
            hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Write();
            //hRatioDeltaRSmoothed[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta]->Write();
          } // Average pair eta loop
          for(int iJetPt = firstDrawnJetPtBin; iJetPt <= lastDrawnJetPtBin; iJetPt++){
            hRatioDeltaRCloseToJet[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Write();
            //hRatioDeltaRCloseToJetSmoothed[TrackPairEfficiencyHistograms::kReconstructed][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Write();
          }
        } // Associated pT loop
      } // Trigger pT loop
    } // Centrality loop
    
    // Write the card to the file
    systemCard->Write(outputFile);

    // Close the output file
    outputFile->Close();
  }
  
}
