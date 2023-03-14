#ifndef TRACKPAIREFFICIENCYHISTOGRAMMANAGER_H
#define TRACKPAIREFFICIENCYHISTOGRAMMANAGER_H

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "TrackPairEfficiencyCard.h"
#include "../src/TrackPairEfficiencyHistograms.h"

/*
 * Class for managing the histograms produced in the track pair efficiency analysis
 */
class TrackPairEfficiencyHistogramManager {
 
public:
  
  // Dimensions for histogram arrays
  static const int kMaxCentralityBins = 5;       // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 15;         // Maximum allowed number of track pT bins
  
  // Indices for different track histogram categories
  enum enumTrackHistograms{kTrack, kUncorrectedTrack, kGenParticle, knTrackCategories};
  
private:
  
  // Naming for jet histograms
  const char* fJetHistogramName = "inclusiveJet";
  const char* fJetAxisName = "Jet";
  
  // Naming for track histograms
  const char* fTrackHistogramNames[knTrackCategories] = {"track", "trackUncorrected", "genParticle"}; // Names that different track histograms have in the input file
  const char* fTrackAxisNames[knTrackCategories] = {"Track", "Uncorrected track", "Gen particle"};    // Names attached to the figure axes
  
  // Naming for track pair histograms
  const char* fTrackPairHistogramNames[TrackPairEfficiencyHistograms::knDataLevels] = {"trackPairs", "genParticlePairs"};
  
  // Naming for data levels
  const char* fDataLevelName[TrackPairEfficiencyHistograms::knDataLevels] = {"", "GeneratorLevel"};
  
public:
  
  TrackPairEfficiencyHistogramManager();                                    // Default constructor
  TrackPairEfficiencyHistogramManager(TFile *inputFile);                    // Constructor
  TrackPairEfficiencyHistogramManager(TFile *inputFile, TrackPairEfficiencyCard *card);   // Constructor with card
  TrackPairEfficiencyHistogramManager(const TrackPairEfficiencyHistogramManager& in);     // Copy constructor
  ~TrackPairEfficiencyHistogramManager();                                   // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void Write(const char* fileName, const char* fileOption);          // Write all the loaded histograms into a file
  void LoadProcessedHistograms(); // Load processed histograms from the inputfile
  
  // Setters for binning information
  void SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true);    // Set up track pT bin indices according to provided bin borders
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for loading event information
  
  // Setters for jets
  void SetLoadJetHistograms(const bool loadOrNot);        // Setter for loading all jet histograms
  
  // Setters for tracks
  void SetLoadTracks(const bool loadOrNot);            // Setter for loading tracks
  void SetLoadTracksUncorrected(const bool loadOrNot); // Setter for loading uncorrected tracks
  void SetLoadGenParticles(const bool loadOrNot);      // Setter for loading generator level tracks
  void SetLoadAllTracks(const bool loadTracks, const bool loadUncorrected, const bool loadGenParticles); // Setter for loading all track histograms
  
  // Setters for track pairs
  void SetLoadTrackPairs(const bool loadOrNot);        // Setter for loading track pairs
  void SetLoadGenParticlePairs(const bool loadOrNot);  // Setter for loading generator level particle pairs
  void SetLoadAllTrackPairs(const bool loadTracks, const bool loadGenParticles); // Setter for loading all track pair histograms
  
  // Setter for loading additional histograms
  void SetLoad2DHistograms(const bool loadOrNot);           // Setter for loading two-dimensional histograms
  
  // Setters for ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for centrality bin range
  void SetTrackPtBinRange(const int first, const int last);    // Setter for track pT bin range
  
  // Getters for number of bins in histograms
  int GetNCentralityBins() const;  // Getter for the number of centrality bins
  int GetNTrackPtBins() const;     // Getter for the number of track pT bins
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  double GetTrackPtBinBorder(const int iTrackPt) const;        // Getter for i:th track pT bin border
  
  // Getters for histogram and axis naming
  const char* GetJetHistogramName() const; // Getter for the jet histogram name
  const char* GetJetAxisName() const;      // Getter for name suitable for x-axis in a jet histogram
  
  TString GetSystem() const;  // Getter for collision system
  
  // Getters for event information histograms
  TH1D* GetHistogramVertexZ() const;            // Getter for z-vertex histogram
  TH1D* GetHistogramVertexZWeighted() const;    // Getter for weighted z-vertex histogram
  TH1D* GetHistogramEvents() const;             // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramCentrality() const;         // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityWeighted() const; // Getter for weighted centrality histogram in all events
  
  // Getters for inclusive jet histograms
  TH1D* GetHistogramJetPt(int iCentrality, const int iDataLevel) const;     // Jet pT histograms
  TH1D* GetHistogramJetPhi(int iCentrality, const int iDataLevel) const;    // Jet phi histograms
  TH1D* GetHistogramJetEta(int iCentrality, const int iDataLevel) const;    // Jet eta histograms
  TH2D* GetHistogramJetEtaPhi(int iCentrality, const int iDataLevel) const; // 2D eta-phi histogram for jets
  
  // Getters for track histograms
  TH1D* GetHistogramTrackPt(const int iTrackType, const int iCentrality) const;                         // Track pT histograms
  TH1D* GetHistogramTrackPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const;    // Track phi histograms
  TH1D* GetHistogramTrackEta(const int iTrackType, const int iCentrality, const int iTrackPt) const;    // Track eta histograms
  TH2D* GetHistogramTrackEtaPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const; // 2D eta-phi histogram for track
  
  // Getters for track pair histograms
  TH1D* GetHistogramTrackPairDeltaR(const int iCentrality, const int iTriggerPt, const int iAssociatedPt, const int iDataLevel) const; // Track pair DeltaR histograms
  
  // Getters for the loaded centrality bins
  int GetFirstCentralityBin() const;  // Get the first loaded centrality bin
  int GetLastCentralityBin() const;   // Get the last loaded centrality bin
  
  // Getters for normalization information
  int GetNEvents() const;                      // Getter for the number of events passing the cuts
  double GetJetPtIntegral(const int iCentrality, const int iDataLevel) const; // Getter for integral over inclusive jet pT in a given centrality
  double GetJetPtIntegral(const int iCentrality, const int iDataLevel, const double minPt, const double maxPt) const; // Getter for integral over inclusive jet pT in a given pT range within a given centrality bin
  
  // Getter for the card
  TrackPairEfficiencyCard* GetCard() const;  // Getter for the JCard
  
private:
  
  // Data members
  TFile *fInputFile;                  // File from which the histograms are read
  TrackPairEfficiencyCard *fCard;                   // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                                        // Load the event information histograms
  bool fLoadJets;                                                    // Load the jet histograms
  bool fLoadTracks[knTrackCategories];                               // Load the track histograms
  bool fLoadTrackPairs[TrackPairEfficiencyHistograms::knDataLevels]; // Load the track pair histograms
  bool fLoad2DHistograms;                                            // Load also two-dimensional (eta,phi) histograms
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  int fFirstLoadedCentralityBin;  // First centrality bin that is loaded
  int fLastLoadedCentralityBin;   // Last centrality bin that is loaded
  int fFirstLoadedTrackPtBin;     // First track pT bin that is loaded
  int fLastLoadedTrackPtBin;      // Last track pT bin that is loaded
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[kMaxCentralityBins+1];    // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[kMaxCentralityBins+1]; // Centrality bin borders, from which bin indices are obtained
  int fTrackPtBinIndices[kMaxTrackPtBins+1];          // Indices for track pT bins in track eta and phi histograms
  double fTrackPtBinBorders[kMaxTrackPtBins+1];       // Track pT bin borders, from which bin indices are obtained
  int fnCentralityBins;                               // Number of centrality bins in the JCard of the data file
  int fnTrackPtBins;                                  // Number of track pT bins in the JCard of the data file
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Event information histograms
  TH1D *fhVertexZ;               // Vertex z position
  TH1D *fhVertexZWeighted;       // Weighted vertex z-position (only meaningfull for MC)
  TH1D *fhEvents;                // Number of events surviving different event cuts
  TH1D *fhCentrality;            // Centrality of all events
  TH1D *fhCentralityWeighted;    // Weighted centrality distribution in all events (only meaningful for MC)
  TH1D *fhPtHat;                 // pT hat for MC events (only meaningful for MC)
  TH1D *fhPtHatWeighted;         // Weighted pT hat distribution (only meaningful for MC)
  TH1D *fhTrackCuts;             // Number of tracks surviving different track cuts
  TH1D *fhGenParticleSelections; // Number of generator level particles surviving different selections
  
  // Histograms for inclusive jets
  TH1D *fhJetPt[kMaxCentralityBins][TrackPairEfficiencyHistograms::knDataLevels];      // Jet pT histograms
  TH1D *fhJetPhi[kMaxCentralityBins][TrackPairEfficiencyHistograms::knDataLevels];     // Jet phi histograms
  TH1D *fhJetEta[kMaxCentralityBins][TrackPairEfficiencyHistograms::knDataLevels];     // Jet eta histograms
  TH2D *fhJetEtaPhi[kMaxCentralityBins][TrackPairEfficiencyHistograms::knDataLevels];  // 2D eta-phi histogram for jets
  
  // Histograms for tracks
  TH1D *fhTrackPt[knTrackCategories][kMaxCentralityBins];                        // Track pT histograms
  TH1D *fhTrackPhi[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];    // Track phi histograms
  TH1D *fhTrackEta[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];    // Track eta histograms
  TH2D *fhTrackEtaPhi[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1]; // 2D eta-phi histograms for tracks
  
  // Histograms for track pairs
  TH1D *fhTrackPairDeltaR[kMaxCentralityBins][kMaxTrackPtBins][kMaxTrackPtBins][TrackPairEfficiencyHistograms::knDataLevels];
  
  // Private methods
  void InitializeFromCard(); // Initialize several member variables from TrackPairEfficiencyCard
  
  // Binning related methods
  void SetBinIndices(const char* histogramName, const int nBins, int *binIndices, const double *binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadJetHistograms();       // Loader for jet histograms
  void LoadTrackHistograms();     // Loader for track histograms
  void LoadTrackPairHistograms(); // Loader for track pair histograms
  
  // Generic setter for bin indice and borders
  void SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double *binBorders, const char* errorMessage, const int maxBins, const bool setIndices); // Generic bin setter
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int& first, int& last); // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Methods for histogram writing
  void WriteJetHistograms();           // Write the jet histograms to the file that is currently open
  void WriteTrackHistograms();         // Write the track histograms to the file that is currently open
  void WriteTrackPairHistograms();     // Write the track pair histograms to the file that is currently open
  
};

#endif
