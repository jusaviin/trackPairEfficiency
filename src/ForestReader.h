// Reader for jet trees from CMS data
//
//===========================================================
// ForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef FORESTREADER_H
#define FORESTREADER_H

// C++ includes
#include <iostream>
#include <assert.h>
#include <vector>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

using namespace std;

class ForestReader{
  
private:
  static const Int_t fnMaxJet = 250;        // Maximum number of jets in an event
  static const Int_t fnMaxTrack = 60000;    // Maximum number of tracks in an event
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                          // Default constructor
  ForestReader(Int_t dataType, Int_t jetType, Int_t jetAxis, Bool_t useTrigger); // Custom constructor
  ForestReader(const ForestReader& in);                    // Copy constructor
  ~ForestReader();                                 // Destructor
  ForestReader& operator=(const ForestReader& obj);        // Equal sign operator
  
  // Methods
  void GetEvent(Int_t nEvent);                 // Get the nth event in tree
  Int_t GetNEvents() const;                        // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> fileList);   // Read the forest from a file list
  void BurnForest();                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  Float_t GetEventWeight() const;     // Getter for event weight in MC
  
  // Getters for leaves in jet tree
  Int_t GetNJets() const;                     // Getter for number of jets
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  Int_t GetNGeneratorJets() const;                   // Getter for number of generator level jets
  Float_t GetGeneratorJetPt(Int_t iJet) const;       // Getter for generator level jet pT
  Float_t GetGeneratorJetPhi(Int_t iJet) const;      // Getter for generator level jet phi
  Float_t GetGeneratorJetEta(Int_t iJet) const;      // Getter for generator level jet eta
  
  // Getters for leaves in HLT tree
  Int_t GetJetFilterBit() const;                     // Getter for the jet filter bit
  
  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  Int_t GetHfCoincidenceFilterBit() const;           // Getter for hadronic forward coincidence filter bit
  Int_t GetClusterCompatibilityFilterBit() const;    // Getter for cluster compatibility filter bit
  
  // Getters for leaves in track tree
  Int_t GetNTracks() const;                                  // Getter for number of tracks
  Float_t GetTrackPt(Int_t iTrack) const;                    // Getter for track pT
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Float_t GetTrackPhi(Int_t iTrack) const;                   // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;                   // Getter for track eta
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackNormalizedChi2(Int_t iTrack) const;        // Getter for normalized track chi2 value from reconstruction fit
  Float_t GetTrackChi2(Int_t iTrack) const;                  // Getter for track chi2 value from reconstruction fit
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;       // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  
  // Getters for leaves in generator level particle tree
  Int_t GetNGenParticles() const;                            // Getter for number of generator level particles
  Float_t GetGenParticlePt(Int_t iTrack) const;              // Getter for generator level particle pT
  Float_t GetGenParticlePhi(Int_t iTrack) const;             // Getter for generator level particle phi
  Float_t GetGenParticleEta(Int_t iTrack) const;             // Getter for generator level particle eta
  Int_t GetGenParticleCharge(Int_t iTrack) const;            // Getter for generator level particle charge
  Int_t GetGenParticleSubevent(Int_t iTrack) const;          // Getter for generator level particle subevent index
  
  // Setter for data type
  void SetDataType(Int_t dataType); // Setter for data type
  
private:
  
  // Methods
  void Initialize();      // Connect the branches to the tree
    
  Int_t fDataType;        // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC
  Int_t fJetType;         // Choose the type of jets usedfor analysis. 0 = Calo jets, 1 = PF jets
  Int_t fJetAxis;         // Jet axis used for the jets. 0 = Anti-kT, 1 = WTA
  Bool_t fUseTrigger;     // Flag for applying jet trigger selection to the analysis
  Bool_t fIsMiniAOD;      // Flag for type of the forest True = MiniAOD forest, False = AOD forest
  
  // Trees in the forest
  TTree *fHeavyIonTree;    // Tree for heavy ion event information
  TTree *fJetTree;         // Tree for jet information
  TTree *fHltTree;         // Tree for high level trigger information
  TTree *fSkimTree;        // Tree for event selection information
  TTree *fTrackTree;       // Tree for reconstructed tracks
  TTree *fGenParticleTree; // Tree for generator level particles
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;                   // Branch for vertex z-position
  TBranch *fHiBinBranch;                  // Branch for centrality
  TBranch *fPtHatBranch;                  // Branch for pT hat
  TBranch *fEventWeightBranch;            // Branch for event weight
  
  // Branches for jet tree
  TBranch *fnJetsBranch;         // Branch for number of jets
  TBranch *fJetPtBranch;         // Branch for jet pT
  TBranch *fJetPhiBranch;        // Branch for jet phi
  TBranch *fJetEtaBranch;        // Branch for jet eta
  TBranch *fJetRawPtBranch;      // Branch for raw jet pT
  TBranch *fJetMaxTrackPtBranch; // Maximum pT for a track inside a jet
  
  TBranch *fnGenJetsBranch;      // Branch for number of generator level jets
  TBranch *fGenJetPtBranch;      // Branch for generator level jet pT
  TBranch *fGenJetPhiBranch;     // Branch for generator level jet phi
  TBranch *fGenJetEtaBranch;     // Branch for generator level jet eta
  
  // Branches for HLT tree
  TBranch *fJetFilterBranch;     // Branch for the jet trigger bit
  
  // Branches for skim tree
  TBranch *fPrimaryVertexBranch;           // Branch for primary vertex filter bit
  TBranch *fBeamScrapingBranch;            // Branch for beam scraping filter bit
  TBranch *fHfCoincidenceBranch;           // Branch for energy recorded in at least 3 HF calorimeter towers
  TBranch *fClusterCompatibilityBranch;    // Branch for cluster compatibility
    
  // Branches for track tree
  TBranch *fnTracksBranch;                    // Branch for number of tracks
  TBranch *fTrackPtBranch;                    // Branch for track pT
  TBranch *fTrackPtErrorBranch;               // Branch for track pT error
  TBranch *fTrackPhiBranch;                   // Branch for track phi
  TBranch *fTrackEtaBranch;                   // Branch for track eta
  TBranch *fHighPurityTrackBranch;            // Branch for high purity of the track
  TBranch *fTrackVertexDistanceZBranch;       // Branch for track distance from primary vertex in z-direction
  TBranch *fTrackVertexDistanceZErrorBranch;  // Branch for error for track distance from primary vertex in z-direction
  TBranch *fTrackVertexDistanceXYBranch;      // Branch for track distance from primary vertex in xy-direction
  TBranch *fTrackVertexDistanceXYErrorBranch; // Branch for error for track distance from primary vertex in xy-direction
  TBranch *fTrackChi2Branch;                  // Branch for track chi2 value from reconstruction fit
  TBranch *fnTrackDegreesOfFreedomBranch;     // Branch for number of degrees of freedom in reconstruction fit
  TBranch *fnHitsTrackerLayerBranch;          // Branch for number of hits in tracker layers
  TBranch *fnHitsTrackBranch;                 // Branch for number of hits for the track
  TBranch *fTrackEnergyEcalBranch;            // Branch for track energy in ECal
  TBranch *fTrackEnergyHcalBranch;            // Branch for track energy in HCal
  
  // Branches for genenerator level particle tree
  TBranch *fGenParticlePtBranch;         // Branch for generator level particle pT:s
  TBranch *fGenParticlePhiBranch;        // Branch for generator level particle phis
  TBranch *fGenParticleEtaBranch;        // Branch for generator level particle etas
  TBranch *fGenParticleChargeBranch;     // Branch for generator level particle charges
  TBranch *fGenParticleSubeventBranch;   // Branch for generator level particle subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat
  
  // Leaves for jet tree
  Int_t fnJets;          // number of jets in an event
  Int_t fnGenJets;       // Number of generator level jets in an event
  Float_t fEventWeight;  // jet weight in the MC tree
  
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
  
  Float_t fGenJetPtArray[fnMaxJet] = {0};      // pT:s of the generator level jets in an event
  Float_t fGenJetPhiArray[fnMaxJet] = {0};     // phis of the generator level jets in an event
  Float_t fGenJetEtaArray[fnMaxJet] = {0};     // etas of the generator level jets in an event
  
  // Leaves for the HLT tree
  Int_t fJetFilterBit;  // Filter bit for the jet trigger
  
  // Leaves for the skim tree
  Int_t fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  Int_t fBeamScrapingFilterBit;            // Filter bit for beam scraping
  Int_t fHfCoincidenceFilterBit;           // Filter bit for energy recorded in at least 3 HF calorimeter towers
  Int_t fClusterCompatibilityFilterBit;    // Filter bit for cluster compatibility
  
  // Leaves for the track tree regardless of forest type
  Int_t fnTracks;  // Number of tracks
  
  // Leaves for the track tree in AOD forests
  Float_t fTrackPtArray[fnMaxTrack] = {0};                    // Array for track pT:s
  Float_t fTrackPtErrorArray[fnMaxTrack] = {0};               // Array for track pT errors
  Float_t fTrackPhiArray[fnMaxTrack] = {0};                   // Array for track phis
  Float_t fTrackEtaArray[fnMaxTrack] = {0};                   // Array for track etas
  Bool_t fHighPurityTrackArray[fnMaxTrack] = {0};             // Array for the high purity of tracks
  Float_t fTrackVertexDistanceZArray[fnMaxTrack] = {0};       // Array for track distance from primary vertex in z-direction
  Float_t fTrackVertexDistanceZErrorArray[fnMaxTrack] = {0};  // Array for error for track distance from primary vertex in z-direction
  Float_t fTrackVertexDistanceXYArray[fnMaxTrack] = {0};      // Array for track distance from primary vertex in xy-direction
  Float_t fTrackVertexDistanceXYErrorArray[fnMaxTrack] = {0}; // Array for error for track distance from primary vertex in xy-direction
  Float_t fTrackChi2Array[fnMaxTrack] = {0};                  // Array for track chi2 value from reconstruction fit
  UChar_t fnTrackDegreesOfFreedomArray[fnMaxTrack] = {0};     // Array for number of degrees of freedom in reconstruction fit
  UChar_t fnHitsTrackerLayerArray[fnMaxTrack] = {0};          // Array for number of hits in tracker layers
  UChar_t fnHitsTrackArray[fnMaxTrack] = {0};                 // Array for number of hits for the track
  Float_t fTrackEnergyEcalArray[fnMaxTrack] = {0};            // Array for track energy in ECal
  Float_t fTrackEnergyHcalArray[fnMaxTrack] = {0};            // Array for track energy in HCal
  
  // Leaves for the track tree in MiniAOD forests
  vector<float> *fTrackPtVector;                    // Vector for track pT:s
  vector<float> *fTrackPtErrorVector;               // Vector for track pT errors
  vector<float> *fTrackPhiVector;                   // Vector for track phis
  vector<float> *fTrackEtaVector;                   // Vector for track etas
  vector<bool> *fHighPurityTrackVector;             // Vector for the high purity of tracks
  vector<float> *fTrackVertexDistanceZVector;       // Vector for track distance from primary vertex in z-direction
  vector<float> *fTrackVertexDistanceZErrorVector;  // Vector for error for track distance from primary vertex in z-direction
  vector<float> *fTrackVertexDistanceXYVector;      // Vector for track distance from primary vertex in xy-direction
  vector<float> *fTrackVertexDistanceXYErrorVector; // Vector for error for track distance from primary vertex in xy-direction
  vector<float> *fTrackNormalizedChi2Vector;        // Vector for normalized track chi2 value from reconstruction fit
  vector<char> *fnHitsTrackerLayerVector;           // Vector for number of hits in tracker layers
  vector<char> *fnHitsTrackVector;                  // Vector for number of hits for the track
  vector<float> *fTrackEnergyEcalVector;            // Vector for track energy in ECal
  vector<float> *fTrackEnergyHcalVector;            // Vector for track energy in HCal
  
  // Leaves for the generator level particle tree
  Int_t fnGenParticles;                     // Number of generator level particles
  vector<float> *fGenParticlePtArray;       // Array for generator level particle pT:s
  vector<float> *fGenParticlePhiArray;      // Array for generator level particle phis
  vector<float> *fGenParticleEtaArray;      // Array for generator level particle etas
  vector<int> *fGenParticleChargeArray;     // Array for generator level particle charges
  vector<int> *fGenParticleSubeventArray;   // Array for generator level particle subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
};

#endif
