// Class for histograms needed in the trigger analysis

#ifndef TRACKPAIREFFICIENCYHISTOGRAMS_H
#define TRACKPAIREFFICIENCYHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class TrackPairEfficiencyHistograms{
  
public:
  
  // Enumeration for event types to event histogram and track cuts for track cut histogram
  enum enumEventTypes {kAll, kPrimaryVertex, kHfCoincidence, kClusterCompatibility, kBeamScraping, kCaloJet, kVzCut, knEventTypes};
  enum enumTrackCuts {kAllTracks, kPtCuts, kEtaCut, kHighPurity, kPtError, kVertexDistance, kCaloSignal, kReconstructionQuality, knTrackCuts};
  enum enumGenParticleSelection {kAllGenParticles, kMcCharge, kMcSube, kGenParticlePtCut, kGenParticleEtaCut, knGenParticleCuts};
  enum enumDataLevel {kReconstructed, kGeneratorLevel, knDataLevels};
  
  // Constructors and destructor
  TrackPairEfficiencyHistograms(); // Default constructor
  TrackPairEfficiencyHistograms(ConfigurationCard *newCard); // Custom constructor
  TrackPairEfficiencyHistograms(const TrackPairEfficiencyHistograms& in); // Copy constructor
  virtual ~TrackPairEfficiencyHistograms(); // Destructor
  TrackPairEfficiencyHistograms& operator=(const TrackPairEfficiencyHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                      // Create all histograms
  void Write() const;                           // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;     // Write the histograms to a file
  void SetCard(ConfigurationCard *newCard);     // Set a new configuration card for the histogram class
  
  // Histograms defined public to allow easier access to them. Should not be abused
  // Notation in comments: l = leading jet, s = subleading jet, inc - inclusive jet, uc = uncorrected, ptw = pT weighted
  TH1F *fhVertexZ;                 // Vertex z-position
  TH1F *fhVertexZWeighted;         // Weighted vertex z-position (only meaningfull for MC)
  TH1F *fhEvents;                  // Number of events. For binning see enumEventTypes.
  TH1F *fhCentrality;              // Centrality information. -0.5 for pp or PYTHIA.
  TH1F *fhCentralityWeighted;      // Weighted centrality distribution (only meaningful for MC)
  TH1F *fhPtHat;                   // pT hat for MC events (only meaningful for MC)
  TH1F *fhPtHatWeighted;           // Weighted pT hat distribution
  TH1F *fhTrackCuts;               // Number of tracks passing cuts. For binning see enumTrackCuts.
  TH1F *fhGenParticleSelections;   // Number of generator level particles passing selections. For binning see enumGenParticleSelection.
  THnSparseF *fhTrack;             // Track histogram. Axes: [pT][phi][eta][cent]
  THnSparseF *fhTrackUncorrected;  // Track histogram for uncorrected tracks. Axes: [uc pT][uc phi][uc eta][cent]
  THnSparseF *fhGenParticle;       // Generator level particle histogram. Axes: [pT][phi][eta][cent]
  THnSparseF *fhInclusiveJet;      // Inclusive jet information. Axes: [jet pT][jet phi][jet eta][cent][reco/gen][trigger]
  THnSparseF *fhTrackPairs;        // Track pair histogram
  THnSparseF *fhGenParticlePairs;  // Generator level particle pair histogram
  
private:
  
  ConfigurationCard *fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "HfCoin2Th4", "ClustCompt", "BeamScrape", "CaloJet", "v_{z} cut"}; // Strings corresponding to event types
  const TString kTrackCutStrings[knTrackCuts] = {"All", "p_{T} cut", "#eta cut", "HighPurity", "p_{T} error", "vertexDist", "caloSignal", "RecoQuality"}; // String corresponding to track cuts
  const TString kGenParticleSelectionStrings[knTrackCuts] = {"All", "MC Charge", "MC sube", "p_{T} cut", "#eta cut"}; // String corresponding to generator level particle selections
  
};

#endif
