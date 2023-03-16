// Class for the main analysis algorithms for track pair efficiency analysis

#ifndef TRACKPAIREFFICIENCYANALYZER_H
#define TRACKPAIREFFICIENCYANALYZER_H

// C++ includes
#include <vector>
#include <bitset>
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>      // For returning several arguments in a transparent manner
#include <fstream>
#include <string>
#include <tuple>

// Root includes
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>

// Own includes
#include "ConfigurationCard.h"
#include "TrackPairEfficiencyHistograms.h"
#include "ForestReader.h"
#include "trackingEfficiency2018PbPb.h"
#include "trackingEfficiency2017pp.h"
#include "TrackingEfficiencyInterface.h"

class TrackPairEfficiencyAnalyzer{
  
private:
  
  enum enumSubeventCuts{kSubeventZero,kSubeventNonZero,kSubeventAny,knSubeventCuts}; // Cuts for subevent index
  enum enumTupleDecoder{kTrackPt, kTrackEta, kTrackPhi, kTrackEfficiency}; // Components of the n-tuple in track vector
  
public:
  
  // Constructors and destructor
  TrackPairEfficiencyAnalyzer(); // Default constructor
  TrackPairEfficiencyAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard); // Custom constructor
  TrackPairEfficiencyAnalyzer(const TrackPairEfficiencyAnalyzer& in); // Copy constructor
  virtual ~TrackPairEfficiencyAnalyzer(); // Destructor
  TrackPairEfficiencyAnalyzer& operator=(const TrackPairEfficiencyAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                     // Run the dijet analysis
  TrackPairEfficiencyHistograms* GetHistograms() const;   // Getter for histograms
  
private:
  
  // Private methods
  void ReadConfigurationFromCard(); // Read all the configuration from the input card
  
  Bool_t PassEventCuts(ForestReader *eventReader); // Check if the event passes the event cuts
  Double_t GetVzWeight(const Double_t vz) const;  // Get the proper vz weighting depending on analyzed system
  Double_t GetCentralityWeight(const Int_t hiBin) const; // Get the proper centrality weighting depending on analyzed system
  Double_t GetJetPtWeight(const Double_t jetPt) const; // Get the proper jet pT weighting for 2017 and 2018 MC
  
  Bool_t PassGenParticleSelection(ForestReader *trackReader, const Int_t iTrack, TH1F *trackCutHistogram, const Bool_t bypassFill);
  Bool_t PassTrackCuts(ForestReader *trackReader, const Int_t iTrack, TH1F *trackCutHistogram, const Bool_t bypassFill);
  Bool_t PassSubeventCut(const Int_t subeventIndex) const;  // Check if the track passes the set subevent cut
  
  Double_t GetTrackEfficiencyCorrection(const Int_t iTrack); // Get the track efficiency correction for a given track
  Double_t  GetTrackEfficiencyCorrection(const Float_t trackPt, const Float_t trackEta, const Int_t hiBin); // Get the track efficiency correction for given track and event information
  
  Double_t GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const; // Get deltaR between two objects
  Double_t GetAveragePhi(const Double_t phi1, const Double_t phi2) const; // Get an average of two phi values
  
  // Private data members
  ForestReader *fEventReader;               // Reader for objects in the event
  std::vector<TString> fFileNames;          // Vector for all the files to loop over
  ConfigurationCard *fCard;                 // Configuration card for the analysis
  TrackPairEfficiencyHistograms *fHistograms;           // Filled histograms
  TF1 *fVzWeightFunction;                   // Weighting function for vz. Needed for MC.
  TF1 *fCentralityWeightFunctionCentral;    // Weighting function for central centrality classes. Needed for MC.
  TF1 *fCentralityWeightFunctionPeripheral; // Weighting function for peripheral centrality classes. Needed for MC.
  TF1 *fPtWeightFunction;                   // Weighting function for jet pT. Needed for MC.
  TrackingEfficiencyInterface *fTrackEfficiencyCorrector2018;  // Tracking efficiency corrector for 2018 PbPb and 2017 pp data.
  
  // Analyzed data and forest types
  Int_t fDataType;                   // Analyzed data type
  Int_t fJetType;                    // Type of jets used for analysis. 0 = Calo jets, 1 = PF jets
  Bool_t fUseTrigger;                // Flag for applying the jet trigger. False = Do not use jet trigger. True = Use jet trigger
  Int_t fDebugLevel;                 // Amount of debug messages printed to console
  
  // Weights for filling the MC histograms
  Double_t fVzWeight;                // Weight for vz in MC
  Double_t fCentralityWeight;        // Weight for centrality in MC
  Double_t fPtHatWeight;             // Weight for pT hat in MC
  Double_t fTotalEventWeight;        // Combined weight factor for MC
  
  // Jet selection cuts
  Int_t fJetAxis;                      // Used jet axis type. 0 = Anti-kT jet axis, 1 = Axis from leading PF candidate
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fMinimumPtHat;              // Minimum accepted pT hat value
  Double_t fMaximumPtHat;              // Maximum accepted pT hat value
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetMinimumPtCut;           // Minimum pT cut for jets
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for jets (and tracks)
  Bool_t fCutBadPhiRegion;             // Cut the phi region with bad tracker performance from the analysis
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT
  
  // Track selection cuts
  Double_t fTrackEtaCut;               // Eta cut around midrapidity
  Double_t fTriggerEtaCut;             // Stricter eta cut for the trigger particle
  Bool_t fCutBadPhiRegionTrigger;      // Do not let the trigger particle to be in the phi region with bad tracker performance
  Double_t fTrackMinPtCut;             // Minimum pT cut
  Double_t fTrackMaxPtCut;             // Maximum pT cut
  Double_t fMaxTrackPtRelativeError;   // Maximum relative error for pT
  Double_t fMaxTrackDistanceToVertex;  // Maximum distance to primary vetrex
  Double_t fCalorimeterSignalLimitPt;  // Require signal in calorimeters for track above this pT
  Double_t fHighPtEtFraction;          // For high pT tracks, minimum required Et as a fraction of track pT
  Double_t fChi2QualityCut;            // Quality cut for track reconstruction
  Double_t fMinimumTrackHits;          // Quality cut for track hits
  Int_t fSubeventCut;                  // Cut for the subevent index

};

#endif
