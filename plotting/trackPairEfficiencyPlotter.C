#include "TrackPairEfficiencyHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "TrackPairEfficiencyCard.h"
#include "JDrawer.h"
#include "../src/TrackPairEfficiencyHistograms.h"

/*
 * Macro for determining track pair efficiency from reconstructed and generator level simulations
 */
void trackPairEfficiencyPlotter(){
  
  // File containing the track pair distributions used to determine the track pair efficiency
  TString fileName = "veryCoolData_processed.root";
  
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
  const int lastDrawnTriggerPtBin = 3;
  
  const int firstDrawnAssociatedPtBin = 0;
  const int lastDrawnAssociatedPtBin = 3;
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
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
  
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < nTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < nTrackPtBins; iAssociatedPt++){
        for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
          hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iDataLevel] = NULL;
        } // Data level loop
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt] = NULL;
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Read the histograms from the file
  for(int iCentrality = firstDrawCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= lastDrawnAssociatedPtBin; iAssociatedPt++){
        for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
          hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iDataLevel] = histograms->GetHistogramTrackPairDeltaR(iCentrality, iTriggerPt, iAssociatedPt, iDataLevel);
        } // Data level loop
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  // Calculate the generator level to reconstructed ratio to determine track pair efficiency
  for(int iCentrality = firstDrawCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTriggerPt = firstDrawnTriggerPtBin; iTriggerPt <= lastDrawnTriggerPtBin; iTriggerPt++){
      for(int iAssociatedPt = firstDrawnAssociatedPtBin; iAssociatedPt <= lastDrawnAssociatedPtBin; iAssociatedPt++){
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt] = (TH1D*) hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kReconstructed]->Clone(Form("ratio%d%d%d", iCentrality, iTriggerPt, iAssociatedPt));
        hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt]->Divide(hTrackDeltaR[iCentrality][iTriggerPt][iAssociatedPt][TrackPairEfficiencyHistograms::kGeneratorLevel]);
      } // Associated pT loo
    } // Trigger pT loop
  } // Centrality loop
  
  
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
        drawer->DrawHistogram(hRatioDeltaR[iCentrality][iTriggerPt][iAssociatedPt], "#DeltaR", "Track pair efficiency", " ");
        
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
}
