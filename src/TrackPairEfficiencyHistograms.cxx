// Histograms needed in the track pair efficiency analysis

// C++ includes
#include <assert.h>

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "TrackPairEfficiencyHistograms.h"

/*
 * Default constructor
 */
TrackPairEfficiencyHistograms::TrackPairEfficiencyHistograms() :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhTrackCuts(0),
  fhGenParticleSelections(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhGenParticle(0),
  fhInclusiveJet(0),
  fhTrackPairs(0),
  fhGenParticlePairs(0),
  fCard(0)
{
  // Default constructor
  
}

/*
 * Custom constructor
 */
TrackPairEfficiencyHistograms::TrackPairEfficiencyHistograms(ConfigurationCard *newCard) :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhTrackCuts(0),
  fhGenParticleSelections(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhGenParticle(0),
  fhInclusiveJet(0),
  fhTrackPairs(0),
  fhGenParticlePairs(0),
  fCard(newCard)
{
  // Custom constructor

}

/*
 * Copy constructor
 */
TrackPairEfficiencyHistograms::TrackPairEfficiencyHistograms(const TrackPairEfficiencyHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
  fhPtHat(in.fhPtHat),
  fhPtHatWeighted(in.fhPtHatWeighted),
  fhTrackCuts(in.fhTrackCuts),
  fhGenParticleSelections(in.fhGenParticleSelections),
  fhTrack(in.fhTrack),
  fhTrackUncorrected(in.fhTrackUncorrected),
  fhGenParticle(in.fhGenParticle),
  fhInclusiveJet(in.fhInclusiveJet),
  fhTrackPairs(in.fhTrackPairs),
  fhGenParticlePairs(in.fhGenParticlePairs),
  fCard(in.fCard)
{
  // Copy constructor
  
}

/*
 * Assingment operator
 */
TrackPairEfficiencyHistograms& TrackPairEfficiencyHistograms::operator=(const TrackPairEfficiencyHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhVertexZWeighted = in.fhVertexZWeighted;
  fhEvents = in.fhEvents;
  fhCentrality = in.fhCentrality;
  fhCentralityWeighted = in.fhCentralityWeighted;
  fhPtHat = in.fhPtHat;
  fhPtHatWeighted = in.fhPtHatWeighted;
  fhTrackCuts = in.fhTrackCuts;
  fhGenParticleSelections = in.fhGenParticleSelections;
  fhTrack = in.fhTrack;
  fhTrackUncorrected = in.fhTrackUncorrected;
  fhGenParticle = in.fhGenParticle;
  fhInclusiveJet = in.fhInclusiveJet;
  fhTrackPairs = in.fhTrackPairs;
  fhGenParticlePairs = in.fhGenParticlePairs;
  fCard = in.fCard;
  
  return *this;
}

/*
 * Destructor
 */
TrackPairEfficiencyHistograms::~TrackPairEfficiencyHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhVertexZWeighted;
  delete fhEvents;
  delete fhCentrality;
  delete fhCentralityWeighted;
  delete fhPtHat;
  delete fhPtHatWeighted;
  delete fhTrackCuts;
  delete fhGenParticleSelections;
  delete fhTrack;
  delete fhTrackUncorrected;
  delete fhGenParticle;
  delete fhInclusiveJet;
  delete fhTrackPairs;
  delete fhGenParticlePairs;
}

/*
 * Set the configuration card used for the histogram class
 */
void TrackPairEfficiencyHistograms::SetCard(ConfigurationCard *newCard){
  fCard = newCard;
}

/*
 * Create the necessary histograms
 */
void TrackPairEfficiencyHistograms::CreateHistograms(){
  
  // ======== Common binning information for histograms =========
  
  // Centrality
  const Double_t minCentrality = -0.75;   // Minimum centrality bin, is negative since hiBin is -1 for pp
  const Double_t maxCentrality = 100.25;  // Maximum centrality bin
  const Int_t nCentralityBins = 202;      // Number of centrality bins
  
  //Track pT
  const Double_t minPtTrack = 0;   // Minimum track pT for track histograms
  const Double_t maxPtTrack = 50;  // Maximum track pT for track histograms
  const Int_t nPtBinsTrack = 500;  // Number of track pT bins for track histograms
  
  // Jet pT
  const Double_t minPtJet = 0;     // Minimum jet pT
  const Double_t maxPtJet = 500;   // Maximum jet pT
  const Int_t nPtBinsJet = 100;    // Number of jet pT bins
  
  // Phi
  const Double_t minPhi = -TMath::Pi();  // Minimum phi
  const Double_t maxPhi = TMath::Pi();   // Maximum phi
  const Int_t nPhiBins = 64;             // Number of phi bins
  
  // Eta
  const Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  const Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  const Int_t nEtaBins = 50;       // Number of eta bins
  
  // Vertex z-position
  const Double_t minVz = -20;   // Minimum vz
  const Double_t maxVz = 20;    // Maximum vz
  const Int_t nVzBins = 80;     // Number of vz bins
  
  // pT hat
  const Double_t minPtHat = 0;     // Minimum pT hat
  const Double_t maxPtHat = 460;   // Maximum pT hat
  const Int_t nFinePtHatBins = 230; // Number of fine pT hat bins
  
  // Reconstructed/generator level
  const Double_t minDataLevel = -0.5;
  const Double_t maxDataLevel = knDataLevels-0.5;
  const Int_t nDataLevelBins = knDataLevels;
  
  // Centrality bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideCentralityBins = fCard->GetNBin("CentralityBinEdges");
  Double_t wideCentralityBins[nWideCentralityBins+1];
  for(Int_t iCentrality = 0; iCentrality < nWideCentralityBins+1; iCentrality++){
    wideCentralityBins[iCentrality] = fCard->Get("CentralityBinEdges",iCentrality);
  }
  
  // Wide track pT bins for the track pair histograms
  const Int_t nWideTrackPtBins = fCard->GetNBin("TrackPairPtBinEdges");
  Double_t wideTrackPtBins[nWideTrackPtBins+1];
  for(Int_t iTrackPt = 0; iTrackPt < nWideTrackPtBins+1; iTrackPt++){
    wideTrackPtBins[iTrackPt] = fCard->Get("TrackPairPtBinEdges",iTrackPt);
  }
  const Double_t minWideTrackPt = wideTrackPtBins[0];
  const Double_t maxWideTrackPt = wideTrackPtBins[nWideTrackPtBins];
  
  // Bins for the pT hat histogram
  const Int_t nPtHatBins = fCard->GetNBin("PtHatBinEdges");
  Double_t ptHatBins[nPtHatBins+1];
  for(Int_t iPtHat = 0; iPtHat < nPtHatBins+1; iPtHat++){
    ptHatBins[iPtHat] = fCard->Get("PtHatBinEdges",iPtHat);
  }
  
  // Logarithmic deltaR binning for energy-energy correlator histograms
  const Int_t nDeltaRBinsEEC = 60;
  const Double_t minDeltaREEC = 0;
  const Double_t maxDeltaREEC = 0.8;
  const Double_t binnerShift = 0.01;
  const Double_t deltaRlogBinWidth = (TMath::Log(maxDeltaREEC+binnerShift) - TMath::Log(minDeltaREEC+binnerShift)) / nDeltaRBinsEEC;
  Double_t deltaRBinsEEC[nDeltaRBinsEEC+1];
  for(int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++){
    deltaRBinsEEC[iDeltaR] = (minDeltaREEC+binnerShift)*TMath::Exp(iDeltaR*deltaRlogBinWidth)-binnerShift;
  }
  
  const Int_t nAxesJet = 5;
  Int_t nBinsJet[nAxesJet];
  Double_t lowBinBorderJet[nAxesJet];
  Double_t highBinBorderJet[nAxesJet];
  
  const Int_t nAxesTrack = 4;
  Int_t nBinsTrack[nAxesTrack];
  Double_t lowBinBorderTrack[nAxesTrack];
  Double_t highBinBorderTrack[nAxesTrack];
  
  const Int_t nAxesTrackPair = 6;
  Int_t nBinsTrackPair[nAxesTrackPair];
  Double_t lowBinBorderTrackPair[nAxesTrackPair];
  Double_t highBinBorderTrackPair[nAxesTrackPair];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1F("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhVertexZWeighted = new TH1F("vertexZweighted","vertexZweighted",nVzBins,minVz,maxVz); fhVertexZWeighted->Sumw2();
  fhEvents = new TH1F("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhCentrality = new TH1F("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  fhCentralityWeighted = new TH1F("centralityWeighted","centralityWeighted",nCentralityBins,minCentrality,maxCentrality); fhCentralityWeighted->Sumw2();
  fhPtHat = new TH1F("pthat","pthat",nPtHatBins,ptHatBins); fhPtHat->Sumw2();
  fhPtHatWeighted = new TH1F("pthatWeighted","pthatWeighted",nFinePtHatBins,minPtHat,maxPtHat); fhPtHatWeighted->Sumw2();
  fhTrackCuts = new TH1F("trackCuts","trackCuts",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCuts->Sumw2();
  fhGenParticleSelections = new TH1F("genParticleSelections","genParticleSelections",knGenParticleCuts,-0.5,knGenParticleCuts-0.5); fhGenParticleSelections->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // For the track cut histogram, label each bin corresponding to a track cut
  for(Int_t i = 0; i < knTrackCuts; i++){
    fhTrackCuts->GetXaxis()->SetBinLabel(i+1,kTrackCutStrings[i]);
  }
  
  // For the generator level particle selection histogram, label each bin corresponding to a particle selection
  for(Int_t i = 0; i < knGenParticleCuts; i++){
    fhGenParticleSelections->GetXaxis()->SetBinLabel(i+1,kGenParticleSelectionStrings[i]);
  }
  
  // If we are using PF jets, change the axis label for that
  if(fCard->Get("JetType") == 1) fhEvents->GetXaxis()->SetBinLabel(kCaloJet+1,"PFJet");
  
  // ======== THnSparses for tracks and uncorrected tracks ========
  
  // Axis 0 for the track histogram: track pT
  nBinsTrack[0] = nPtBinsTrack;         // nBins for track pT
  lowBinBorderTrack[0] = minPtTrack;    // low bin border for track pT
  highBinBorderTrack[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track histogram: track phi
  nBinsTrack[1] = nPhiBins;         // nBins for track phi
  lowBinBorderTrack[1] = minPhi;    // low bin border for track phi
  highBinBorderTrack[1] = maxPhi;   // high bin border for track phi
  
  // Axis 2 for the track histogram: track eta
  nBinsTrack[2] = nEtaBins;         // nBins for track eta
  lowBinBorderTrack[2] = minEta;    // low bin border for track eta
  highBinBorderTrack[2] = maxEta;   // high bin border for track eta
  
  // Axis 3 for the track histogram: centrality
  nBinsTrack[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderTrack[3] = minCentrality;  // low bin border for centrality
  highBinBorderTrack[3] = maxCentrality; // high bin border for centrality
  
  // Create the histograms for tracks and uncorrected tracks using the above binning information
  fhTrack = new THnSparseF("track","track",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrack->Sumw2();
  fhTrackUncorrected = new THnSparseF("trackUncorrected","trackUncorrected",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackUncorrected->Sumw2();
  fhGenParticle = new THnSparseF("genParticle","genParticle",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhGenParticle->Sumw2();

  // Set custom centrality bins for histograms
  fhTrack->SetBinEdges(3,wideCentralityBins);
  fhTrackUncorrected->SetBinEdges(3,wideCentralityBins);
  fhGenParticle->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the jet histogram: jet pT
  nBinsJet[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorderJet[0] = minPtJet;    // low bin border for any jet pT
  highBinBorderJet[0] = maxPtJet;   // high bin border for any jet pT
  
  // Axis 1 for the jet histogram: jet phi
  nBinsJet[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorderJet[1] = minPhi;   // low bin border for any jet phi
  highBinBorderJet[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the jet histogram: jet eta
  nBinsJet[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorderJet[2] = minEta;   // low bin border for any jet eta
  highBinBorderJet[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the jet histogram: centrality
  nBinsJet[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderJet[3] = minCentrality;  // low bin border for centrality
  highBinBorderJet[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: data level (reconstructed / generator level)
  nBinsJet[4] = nDataLevelBins;       // nBins different data levels
  lowBinBorderJet[4] = minDataLevel;  // low bin border for data levels
  highBinBorderJet[4] = maxDataLevel; // high bin border for data levels
  
  // Create the histogram for all jets using the above binning information
  fhInclusiveJet = new THnSparseF("inclusiveJet","inclusiveJet",nAxesJet,nBinsJet,lowBinBorderJet,highBinBorderJet); fhInclusiveJet->Sumw2();

  // Set custom centrality bins for histograms
  fhInclusiveJet->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for track pairs ========
  
  // Axis 0 for the track pair histogram: deltaR between the two tracks
  nBinsTrackPair[0] = nDeltaRBinsEEC;         // nBins for deltaR between the tracks
  lowBinBorderTrackPair[0] = minDeltaREEC;    // low bin border for deltaR between the tracks
  highBinBorderTrackPair[0] = maxDeltaREEC;   // high bin border for deltaR between the tracks
  
  // Axis 1 for the track pair histogram: track pT for the first track
  nBinsTrackPair[1] = nWideTrackPtBins;         // nBins for wide track pT
  lowBinBorderTrackPair[1] = minWideTrackPt;    // low bin border for wide track pT
  highBinBorderTrackPair[1] = maxWideTrackPt;   // high bin border for wide track pT
  
  // Axis 2 for the track pair histogram: track pT for the second track
  nBinsTrackPair[2] = nWideTrackPtBins;         // nBins for wide track pT
  lowBinBorderTrackPair[2] = minWideTrackPt;    // low bin border for wide track pT
  highBinBorderTrackPair[2] = maxWideTrackPt;   // high bin border for wide track pT
  
  // Axis 3 for the track histogram: average pair phi position
  nBinsTrackPair[3] = nPhiBins;         // nBins for track eta
  lowBinBorderTrackPair[3] = minPhi;    // low bin border for track eta
  highBinBorderTrackPair[3] = maxPhi;   // high bin border for track eta
  
  // Axis 4 for the track histogram: average pair eta position
  nBinsTrackPair[4] = nEtaBins;         // nBins for track eta
  lowBinBorderTrackPair[4] = minEta;    // low bin border for track eta
  highBinBorderTrackPair[4] = maxEta;   // high bin border for track eta
  
  // Axis 5 for the track histogram: centrality
  nBinsTrackPair[5] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderTrackPair[5] = minCentrality;  // low bin border for centrality
  highBinBorderTrackPair[5] = maxCentrality; // high bin border for centrality
  
  // Create the histograms for tracks and uncorrected tracks using the above binning information
  fhTrackPairs = new THnSparseF("trackPairs","trackPairs",nAxesTrackPair,nBinsTrackPair,lowBinBorderTrackPair,highBinBorderTrackPair); fhTrackPairs->Sumw2();
  fhGenParticlePairs = new THnSparseF("genParticlePairs","genParticlePairs",nAxesTrackPair,nBinsTrackPair,lowBinBorderTrackPair,highBinBorderTrackPair); fhGenParticlePairs->Sumw2();

  // Set custom deltaR bins for histograms
  fhTrackPairs->SetBinEdges(0,deltaRBinsEEC);
  fhGenParticlePairs->SetBinEdges(0,deltaRBinsEEC);
  
  // Set custom track pT bins for histograms
  fhTrackPairs->SetBinEdges(1,wideTrackPtBins);
  fhGenParticlePairs->SetBinEdges(1,wideTrackPtBins);
  fhTrackPairs->SetBinEdges(2,wideTrackPtBins);
  fhGenParticlePairs->SetBinEdges(2,wideTrackPtBins);
  
  // Set custom centrality bins for histograms
  fhTrackPairs->SetBinEdges(5,wideCentralityBins);
  fhGenParticlePairs->SetBinEdges(5,wideCentralityBins);
}

/*
 * Write the histograms to file
 */
void TrackPairEfficiencyHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhVertexZWeighted->Write();
  fhEvents->Write();
  fhCentrality->Write();
  fhCentralityWeighted->Write();
  fhPtHat->Write();
  fhPtHatWeighted->Write();
  fhTrackCuts->Write();
  fhGenParticleSelections->Write();
  fhTrack->Write();
  fhTrackUncorrected->Write();
  fhGenParticle->Write();
  fhInclusiveJet->Write();
  fhTrackPairs->Write();
  fhGenParticlePairs->Write();
  
}

/*
 * Write the histograms to a given file
 */
void TrackPairEfficiencyHistograms::Write(TString outputFileName) const{
  
  // Define the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
}


