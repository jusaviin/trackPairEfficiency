// Implementation for ForestReader

// Own includes
#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fDataType(0),
  fJetType(0),
  fJetAxis(0),
  fUseTrigger(false),
  fIsMiniAOD(false),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fGenParticleTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetEtaBranch(0),
  fJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fnTracksBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fGenParticlePtBranch(0),
  fGenParticlePhiBranch(0),
  fGenParticleEtaBranch(0),
  fGenParticleChargeBranch(0),
  fGenParticleSubeventBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fnJets(0),
  fnGenJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fGenJetPtArray(),
  fGenJetPhiArray(),
  fGenJetEtaArray(),
  fJetFilterBit(1),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fnTracks(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fnGenParticles(0),
  fGenParticlePtArray(0),
  fGenParticlePhiArray(0),
  fGenParticleEtaArray(0),
  fGenParticleChargeArray(0),
  fGenParticleSubeventArray(0)
{
  // Default constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC
 *   Int_t jetType: 0 = Calo jets, 1 = CSPF jets, 2 = PuPF jets, 3 = Flow jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 *   Bool_t useTrigger: True: Apply jet trigger selection. False: Do not check trigger
 */
ForestReader::ForestReader(Int_t dataType, Int_t jetType, Int_t jetAxis, Bool_t useTrigger) :
  fDataType(0),
  fJetType(jetType),
  fJetAxis(jetAxis),
  fUseTrigger(useTrigger),
  fIsMiniAOD(false),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fGenParticleTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetEtaBranch(0),
  fJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fnTracksBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fGenParticlePtBranch(0),
  fGenParticlePhiBranch(0),
  fGenParticleEtaBranch(0),
  fGenParticleChargeBranch(0),
  fGenParticleSubeventBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fnJets(0),
  fnGenJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fGenJetPtArray(),
  fGenJetPhiArray(),
  fGenJetEtaArray(),
  fJetFilterBit(1),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fnTracks(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fnGenParticles(0),
  fGenParticlePtArray(0),
  fGenParticlePhiArray(0),
  fGenParticleEtaArray(0),
  fGenParticleChargeArray(0),
  fGenParticleSubeventArray(0)
{
  // Custom constructor
  
  SetDataType(dataType);
  
  // Initialize fJetMaxTrackPtArray to -1
  for(int i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fDataType(in.fDataType),
  fJetType(in.fJetType),
  fJetAxis(in.fJetAxis),
  fUseTrigger(in.fUseTrigger),
  fIsMiniAOD(in.fIsMiniAOD),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fGenParticleTree(in.fGenParticleTree),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fPtHatBranch(in.fPtHatBranch),
  fEventWeightBranch(in.fEventWeightBranch),
  fnJetsBranch(in.fnJetsBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fnGenJetsBranch(in.fnGenJetsBranch),
  fGenJetPtBranch(in.fGenJetPtBranch),
  fGenJetPhiBranch(in.fGenJetPhiBranch),
  fGenJetEtaBranch(in.fGenJetEtaBranch),
  fJetFilterBranch(in.fJetFilterBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fBeamScrapingBranch(in.fBeamScrapingBranch),
  fHfCoincidenceBranch(in.fHfCoincidenceBranch),
  fClusterCompatibilityBranch(in.fClusterCompatibilityBranch),
  fnTracksBranch(in.fnTracksBranch),
  fTrackPtBranch(in.fTrackPtBranch),
  fTrackPtErrorBranch(in.fTrackPtErrorBranch),
  fTrackPhiBranch(in.fTrackPhiBranch),
  fTrackEtaBranch(in.fTrackEtaBranch),
  fHighPurityTrackBranch(in.fHighPurityTrackBranch),
  fTrackVertexDistanceZBranch(in.fTrackVertexDistanceZBranch),
  fTrackVertexDistanceZErrorBranch(in.fTrackVertexDistanceZErrorBranch),
  fTrackVertexDistanceXYBranch(in.fTrackVertexDistanceXYBranch),
  fTrackVertexDistanceXYErrorBranch(in.fTrackVertexDistanceXYErrorBranch),
  fTrackChi2Branch(in.fTrackChi2Branch),
  fnTrackDegreesOfFreedomBranch(in.fnTrackDegreesOfFreedomBranch),
  fnHitsTrackerLayerBranch(in.fnHitsTrackerLayerBranch),
  fnHitsTrackBranch(in.fnHitsTrackBranch),
  fTrackEnergyEcalBranch(in.fTrackEnergyEcalBranch),
  fTrackEnergyHcalBranch(in.fTrackEnergyHcalBranch),
  fGenParticlePtBranch(in.fGenParticlePtBranch),
  fGenParticlePhiBranch(in.fGenParticlePhiBranch),
  fGenParticleEtaBranch(in.fGenParticleEtaBranch),
  fGenParticleChargeBranch(in.fGenParticleChargeBranch),
  fGenParticleSubeventBranch(in.fGenParticleSubeventBranch),
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fPtHat(in.fPtHat),
  fnJets(in.fnJets),
  fnGenJets(in.fnGenJets),
  fEventWeight(in.fEventWeight),
  fJetFilterBit(in.fJetFilterBit),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fHfCoincidenceFilterBit(in.fHfCoincidenceFilterBit),
  fClusterCompatibilityFilterBit(in.fClusterCompatibilityFilterBit),
  fnTracks(in.fnTracks),
  fTrackPtVector(in.fTrackPtVector),
  fTrackPhiVector(in.fTrackPhiVector),
  fTrackEtaVector(in.fTrackEtaVector),
  fHighPurityTrackVector(in.fHighPurityTrackVector),
  fTrackVertexDistanceZVector(in.fTrackVertexDistanceZVector),
  fTrackVertexDistanceZErrorVector(in.fTrackVertexDistanceZErrorVector),
  fTrackVertexDistanceXYVector(in.fTrackVertexDistanceXYVector),
  fTrackVertexDistanceXYErrorVector(in.fTrackVertexDistanceXYErrorVector),
  fTrackNormalizedChi2Vector(in.fTrackNormalizedChi2Vector),
  fnHitsTrackerLayerVector(in.fnHitsTrackerLayerVector),
  fnHitsTrackVector(in.fnHitsTrackVector),
  fTrackEnergyEcalVector(in.fTrackEnergyEcalVector),
  fTrackEnergyHcalVector(in.fTrackEnergyHcalVector),
  fnGenParticles(in.fnGenParticles),
  fGenParticlePtArray(in.fGenParticlePtArray),
  fGenParticlePhiArray(in.fGenParticlePhiArray),
  fGenParticleEtaArray(in.fGenParticleEtaArray),
  fGenParticleChargeArray(in.fGenParticleChargeArray),
  fGenParticleSubeventArray(in.fGenParticleSubeventArray)
{
  // Copy constructor
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
    fGenJetPtArray[i] = in.fGenJetPtArray[i];
    fGenJetPhiArray[i] = in.fGenJetPhiArray[i];
    fGenJetEtaArray[i] = in.fGenJetEtaArray[i];
  }
  
  // Copy the track arrays
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
  }
}

/*
 * Assignment operator
 */
ForestReader& ForestReader::operator=(const ForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fJetType = in.fJetType;
  fJetAxis = in.fJetAxis;
  fUseTrigger = in.fUseTrigger;
  fIsMiniAOD = in.fIsMiniAOD;
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fGenParticleTree = in.fGenParticleTree;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fPtHatBranch = in.fPtHatBranch;
  fEventWeightBranch = in.fEventWeightBranch;
  fnJetsBranch = in.fnJetsBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fnGenJetsBranch = in.fnGenJetsBranch;
  fGenJetPtBranch = in.fGenJetPtBranch;
  fGenJetPhiBranch = in.fGenJetPhiBranch;
  fGenJetEtaBranch = in.fGenJetEtaBranch;
  fJetFilterBranch = in.fJetFilterBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fBeamScrapingBranch = in.fBeamScrapingBranch;
  fHfCoincidenceBranch = in.fHfCoincidenceBranch;
  fClusterCompatibilityBranch = in.fClusterCompatibilityBranch;
  fnTracksBranch = in.fnTracksBranch;
  fTrackPtBranch = in.fTrackPtBranch;
  fTrackPtErrorBranch = in.fTrackPtErrorBranch;
  fTrackPhiBranch = in.fTrackPhiBranch;
  fTrackEtaBranch = in.fTrackEtaBranch;
  fHighPurityTrackBranch = in.fHighPurityTrackBranch;
  fTrackVertexDistanceZBranch = in.fTrackVertexDistanceZBranch;
  fTrackVertexDistanceZErrorBranch = in.fTrackVertexDistanceZErrorBranch;
  fTrackVertexDistanceXYBranch = in.fTrackVertexDistanceXYBranch;
  fTrackVertexDistanceXYErrorBranch = in.fTrackVertexDistanceXYErrorBranch;
  fTrackChi2Branch = in.fTrackChi2Branch;
  fnTrackDegreesOfFreedomBranch = in.fnTrackDegreesOfFreedomBranch;
  fnHitsTrackerLayerBranch = in.fnHitsTrackerLayerBranch;
  fnHitsTrackBranch = in.fnHitsTrackBranch;
  fTrackEnergyEcalBranch = in.fTrackEnergyEcalBranch;
  fTrackEnergyHcalBranch = in.fTrackEnergyHcalBranch;
  fGenParticlePtBranch = in.fGenParticlePtBranch;
  fGenParticlePhiBranch = in.fGenParticlePhiBranch;
  fGenParticleEtaBranch = in.fGenParticleEtaBranch;
  fGenParticleChargeBranch = in.fGenParticleChargeBranch;
  fGenParticleSubeventBranch = in.fGenParticleSubeventBranch;
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fPtHat = in.fPtHat;
  fnJets = in.fnJets;
  fnGenJets = in.fnGenJets;
  fEventWeight = in.fEventWeight;
  fJetFilterBit = in.fJetFilterBit;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fHfCoincidenceFilterBit = in.fHfCoincidenceFilterBit;
  fClusterCompatibilityFilterBit = in.fClusterCompatibilityFilterBit;
  fnTracks = in.fnTracks;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
    fGenJetPtArray[i] = in.fGenJetPtArray[i];
    fGenJetPhiArray[i] = in.fGenJetPhiArray[i];
    fGenJetEtaArray[i] = in.fGenJetEtaArray[i];
  }
  
  // Copy the track arrays
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
  }
  
  // Copy the track vectors
  fTrackPtVector = in.fTrackPtVector;
  fTrackPtVector = in.fTrackPtVector;
  fTrackPhiVector = in.fTrackPhiVector;
  fTrackEtaVector = in.fTrackEtaVector;
  fHighPurityTrackVector = in.fHighPurityTrackVector;
  fTrackVertexDistanceZVector = in.fTrackVertexDistanceZVector;
  fTrackVertexDistanceZErrorVector = in.fTrackVertexDistanceZErrorVector;
  fTrackVertexDistanceXYVector = in.fTrackVertexDistanceXYVector;
  fTrackVertexDistanceXYErrorVector = in.fTrackVertexDistanceXYErrorVector;
  fTrackNormalizedChi2Vector = in.fTrackNormalizedChi2Vector;
  fnHitsTrackerLayerVector = in.fnHitsTrackerLayerVector;
  fnHitsTrackVector = in.fnHitsTrackVector;
  fTrackEnergyEcalVector = in.fTrackEnergyEcalVector;
  fTrackEnergyHcalVector = in.fTrackEnergyHcalVector;
  
  // Copy the generator level particle vectors
  fnGenParticles = in.fnGenParticles;
  fGenParticlePtArray = in.fGenParticlePtArray;
  fGenParticlePhiArray = in.fGenParticlePhiArray;
  fGenParticleEtaArray = in.fGenParticleEtaArray;
  fGenParticleChargeArray = in.fGenParticleChargeArray;
  fGenParticleSubeventArray = in.fGenParticleSubeventArray;
  
  return *this;
}

/*
 * Destructor
 */
ForestReader::~ForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void ForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fHeavyIonTree->SetBranchStatus("pthat",1);
    fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
    fHeavyIonTree->SetBranchStatus("weight",1);
    fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch); // event weight only for MC
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
    fEventWeight = 1;
  }
  
  // Connect the branches to the jet tree
  const char *jetAxis[2] = {"jt", "WTA"};
  const char *genJetAxis[2] = {"", "WTA"};
  char branchName[30];
  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("jtpt",1);
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  
  // If specified, select WTA axis for jet phi
  sprintf(branchName,"%sphi",jetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetPhiArray,&fJetPhiBranch);
  
  // If specified, select WTA axis for jet eta
  sprintf(branchName,"%seta",jetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetEtaArray,&fJetEtaBranch);
  
  fJetTree->SetBranchStatus("nref",1);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchStatus("rawpt",1);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchStatus("trackMax",1);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
  // If we are looking at Monte Carlo, connect the reference pT and parton arrays
  if(fDataType > kPbPb){
    fJetTree->SetBranchStatus("genpt",1);
    fJetTree->SetBranchAddress("genpt",&fGenJetPtArray,&fGenJetPtBranch);
    
    // If specified, select WTA axis for jet phi
    sprintf(branchName,"%sgenphi",genJetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fGenJetPhiArray,&fGenJetPhiBranch);
    
    // If specified, select WTA axis for jet eta
    sprintf(branchName,"%sgeneta",genJetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fGenJetEtaArray,&fGenJetEtaBranch);
    
    fJetTree->SetBranchStatus("ngen",1);
    fJetTree->SetBranchAddress("ngen",&fnGenJets,&fnGenJetsBranch);
  }
  
  // Event selection summary
  //
  //         tree                      branch                         What it is
  //  hltanalysis/HltTree   HLT_HIPuAK4CaloJet100_Eta5p1_v1      Event selection for PbPb
  //  hltanalysis/HltTree      HLT_AK4CaloJet80_Eta5p1_v1         Event selection for pp
  // skimanalysis/HltTree         pprimaryVertexFilter           Event selection for PbPb
  // skimanalysis/HltTree    HBHENoiseFilterResultRun2Loose   Event selection for pp and PbPb
  // skimanalysis/HltTree         pPAprimaryVertexFilter          Event selection for pp
  // skimanalysis/HltTree           pBeamScrapingFilter           Event selection for pp
  
  // Connect the branches to the HLT tree
  if(fUseTrigger){
    fHltTree->SetBranchStatus("*",0);
    
    if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
      
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1",1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1", &fJetFilterBit, &fJetFilterBranch);
      
    } else { // PbPb data or MC
      
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1",1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1", &fJetFilterBit, &fJetFilterBranch);
    }
  } else {
    fJetFilterBit = 1; // Jet trigger filter disabled
  }
  
  
  // Connect the branches to the skim tree (different for pp and PbPb data and Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  // pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter
  if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else { // PbPb data or MC
    
    // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    
    // Cut on noise on HCAL
    if(fIsMiniAOD){
      // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
      fSkimTree->SetBranchStatus("pphfCoincFilter2Th4",1);
      fSkimTree->SetBranchAddress("pphfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
    } else {
      
      // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
      fSkimTree->SetBranchStatus("phfCoincFilter2Th4",1);
      fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
    }
    
    // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  }
  
  // Connect the branches to the track tree
  
  fTrackTree->SetBranchStatus("*",0);
  
  // We need to read the forest to vectors for MiniAODs and to arrays for AODs
  if(fIsMiniAOD){
    
    fTrackTree->SetBranchStatus("trkPt",1);
    fTrackTree->SetBranchAddress("trkPt",&fTrackPtVector,&fTrackPtBranch);
    fTrackTree->SetBranchStatus("trkPtError",1);
    fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorVector,&fTrackPtErrorBranch);
    fTrackTree->SetBranchStatus("trkPhi",1);
    fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiVector,&fTrackPhiBranch);
    fTrackTree->SetBranchStatus("trkEta",1);
    fTrackTree->SetBranchAddress("trkEta",&fTrackEtaVector,&fTrackEtaBranch);
    fTrackTree->SetBranchStatus("nTrk",1);
    fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
    fTrackTree->SetBranchStatus("highPurity",1);
    fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackVector,&fHighPurityTrackBranch);
    fTrackTree->SetBranchStatus("trkDzFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDzFirstVtx",&fTrackVertexDistanceZVector,&fTrackVertexDistanceZBranch);
    fTrackTree->SetBranchStatus("trkDzErrFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDzErrFirstVtx",&fTrackVertexDistanceZErrorVector,&fTrackVertexDistanceZErrorBranch);
    fTrackTree->SetBranchStatus("trkDxyFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDxyFirstVtx",&fTrackVertexDistanceXYVector,&fTrackVertexDistanceXYBranch);
    fTrackTree->SetBranchStatus("trkDxyErrFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDxyErrFirstVtx",&fTrackVertexDistanceXYErrorVector,&fTrackVertexDistanceXYErrorBranch);
    fTrackTree->SetBranchStatus("trkNormChi2",1);
    fTrackTree->SetBranchAddress("trkNormChi2",&fTrackNormalizedChi2Vector,&fTrackChi2Branch);
    fTrackTree->SetBranchStatus("trkNLayers",1);
    fTrackTree->SetBranchAddress("trkNLayers",&fnHitsTrackerLayerVector,&fnHitsTrackerLayerBranch);
    fTrackTree->SetBranchStatus("trkNHits",1);
    fTrackTree->SetBranchAddress("trkNHits",&fnHitsTrackVector,&fnHitsTrackBranch);
    fTrackTree->SetBranchStatus("pfEcal",1);
    fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalVector,&fTrackEnergyEcalBranch);
    fTrackTree->SetBranchStatus("pfHcal",1);
    fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalVector,&fTrackEnergyHcalBranch);
    
  } else { // Read the tree from AOD files
    
    fTrackTree->SetBranchStatus("trkPt",1);
    fTrackTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
    fTrackTree->SetBranchStatus("trkPtError",1);
    fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
    fTrackTree->SetBranchStatus("trkPhi",1);
    fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
    fTrackTree->SetBranchStatus("trkEta",1);
    fTrackTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
    fTrackTree->SetBranchStatus("nTrk",1);
    fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
    fTrackTree->SetBranchStatus("highPurity",1);
    fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
    fTrackTree->SetBranchStatus("trkDz1",1);
    fTrackTree->SetBranchAddress("trkDz1",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
    fTrackTree->SetBranchStatus("trkDzError1",1);
    fTrackTree->SetBranchAddress("trkDzError1",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
    fTrackTree->SetBranchStatus("trkDxy1",1);
    fTrackTree->SetBranchAddress("trkDxy1",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
    fTrackTree->SetBranchStatus("trkDxyError1",1);
    fTrackTree->SetBranchAddress("trkDxyError1",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
    fTrackTree->SetBranchStatus("trkChi2",1);
    fTrackTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
    fTrackTree->SetBranchStatus("trkNdof",1);
    fTrackTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
    fTrackTree->SetBranchStatus("trkNlayer",1);
    fTrackTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
    fTrackTree->SetBranchStatus("trkNHit",1);
    fTrackTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
    fTrackTree->SetBranchStatus("pfEcal",1);
    fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
    fTrackTree->SetBranchStatus("pfHcal",1);
    fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
  }
  
  // Connect the branches to the generator level particle tree
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fGenParticleTree->SetBranchStatus("*",0);
    fGenParticleTree->SetBranchStatus("pt",1);
    fGenParticleTree->SetBranchAddress("pt",&fGenParticlePtArray,&fGenParticlePtBranch);
    fGenParticleTree->SetBranchStatus("phi",1);
    fGenParticleTree->SetBranchAddress("phi",&fGenParticlePhiArray,&fGenParticlePhiBranch);
    fGenParticleTree->SetBranchStatus("eta",1);
    fGenParticleTree->SetBranchAddress("eta",&fGenParticleEtaArray,&fGenParticleEtaBranch);
    fGenParticleTree->SetBranchStatus("chg",1);
    fGenParticleTree->SetBranchAddress("chg",&fGenParticleChargeArray,&fGenParticleChargeBranch);
    fGenParticleTree->SetBranchStatus("sube",1);
    fGenParticleTree->SetBranchAddress("sube",&fGenParticleSubeventArray,&fGenParticleSubeventBranch);
  } // Reading track trees
  
}


/*
 * Setter for fDataType
 */
void ForestReader::SetDataType(Int_t dataType){
  
  //Sanity check for given data type
  if(dataType < 0 || dataType > knDataTypes-1){
    cout << "ERROR: Data type input " << dataType << " is invalid in ForestReader.cxx!" << endl;
    cout << "Please give integer between 0 and " << knDataTypes-1 << "." << endl;
    cout << "Setting data type to 0 (pp)." << endl;
    fDataType = 0;
  } else {
    
    // If the sanity check passes, set the given data type
    fDataType = dataType;
  }
}

/*
 * Connect a new tree to the reader
 */
void ForestReader::ReadForestFromFile(TFile *inputFile){
  
  // When reading a forest, we need to check if it is AOD or MiniAOD forest as there are some differences
  // The HiForest tree is renamed to HiForestInfo in MiniAODs, so we can determine the forest type from this.
  TTree* miniAODcheck = (TTree*)inputFile->Get("HiForestInfo/HiForest");
  fIsMiniAOD = !(miniAODcheck == NULL);
  
  // Helper variable for finding the correct tree
  const char *treeName[4] = {"none","none","none","none"};
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  if(fUseTrigger) fHltTree = (TTree*)inputFile->Get("hltanalysis/HltTree");
  fSkimTree = (TTree*)inputFile->Get("skimanalysis/HltTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC){
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for puPF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  }
  
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  
  // The track tree has different name for pp and PbPb
  if(fIsMiniAOD && (fDataType == kPbPb || fDataType == kPbPbMC)){
    fTrackTree = (TTree*)inputFile->Get("PbPbTracks/trackTree");
  } else {
    fTrackTree = (TTree*)inputFile->Get("ppTrack/trackTree");
  }
  
  // Read the generator level particle tree if we are analyzing simulation
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fGenParticleTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
  }
  
  Initialize();
}

/*
 * Connect a new tree to the reader
 */
void ForestReader::ReadForestFromFileList(std::vector<TString> fileList){
  TFile *inputFile = TFile::Open(fileList.at(0));
  ReadForestFromFile(inputFile);
}

/*
 * Burn the current forest.
 */
void ForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  if(fUseTrigger) fHltTree->Delete();
  fSkimTree->Delete();
  fJetTree->Delete();
  fTrackTree->Delete();
  if(fDataType == kPpMC || fDataType == kPbPbMC) fGenParticleTree->Delete();
}

/*
 * Load an event to memory
 */
void ForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  if(fUseTrigger) fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
  if(fDataType == kPpMC || fDataType == kPbPbMC) {
    fGenParticleTree->GetEntry(nEvent);
   
    // Read the numbers of generator level particles for this event
    fnGenParticles = fGenParticlePtArray->size();
  }
}

// Getter for number of events in the tree
Int_t ForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for number of jets in an event
Int_t ForestReader::GetNJets() const{
  return fnJets;
}

// Getter for number of jets in an event
Int_t ForestReader::GetNGeneratorJets() const{
  return fnGenJets;
}

// Getter for jet pT
Float_t ForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t ForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t ForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for jet raw pT
Float_t ForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t ForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Getter for generator level jet pT
Float_t ForestReader::GetGeneratorJetPt(Int_t iJet) const{
  return fGenJetPtArray[iJet];
}

// Getter for generator level jet phi
Float_t ForestReader::GetGeneratorJetPhi(Int_t iJet) const{
  return fGenJetPhiArray[iJet];
}

// Getter for generator level jet eta
Float_t ForestReader::GetGeneratorJetEta(Int_t iJet) const{
  return fGenJetEtaArray[iJet];
}

// Getter for vertex z position
Float_t ForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. CMS has integer centrality bins from 0 to 200, thus division by 2.
Float_t ForestReader::GetCentrality() const{
  return fHiBin/2.0;
}

// Getter for hiBin. Return 1 for negative values (for easier handling of tracking efficiency correction)
Int_t ForestReader::GetHiBin() const{
  if(fHiBin < 0) return 1;
  return fHiBin;
}

// Getter for pT hat
Float_t ForestReader::GetPtHat() const{
  return fPtHat;
}

// Getter for pT hat
Float_t ForestReader::GetEventWeight() const{
  return fEventWeight;
}

// Getter for the selected jet filter bit
Int_t ForestReader::GetJetFilterBit() const{
  return fJetFilterBit;
}

// Getter for primary vertex filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetPrimaryVertexFilterBit() const{
  return fPrimaryVertexFilterBit;
}

// Getter for beam scraping filter bit. Always 1 for MC and PbPb (set in the initializer).
Int_t ForestReader::GetBeamScrapingFilterBit() const{
  return fBeamScrapingFilterBit;
}

// Getter for HF energy coincidence filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetHfCoincidenceFilterBit() const{
  return fHfCoincidenceFilterBit;
}

// Getter for cluster compatibility filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetClusterCompatibilityFilterBit() const{
  return fClusterCompatibilityFilterBit;
}

// Getter for number of tracks in an event
Int_t ForestReader::GetNTracks() const{
  return fnTracks;
}

// Getter for track pT
Float_t ForestReader::GetTrackPt(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPtVector->at(iTrack);
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t ForestReader::GetTrackPtError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPtErrorVector->at(iTrack);
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t ForestReader::GetTrackPhi(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPhiVector->at(iTrack);
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t ForestReader::GetTrackEta(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEtaVector->at(iTrack);
  return fTrackEtaArray[iTrack];
}

// Getter for high purity of the track
Bool_t ForestReader::GetTrackHighPurity(Int_t iTrack) const{
  if(fIsMiniAOD) return fHighPurityTrackVector->at(iTrack);
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t ForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceZVector->at(iTrack);
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t ForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceZErrorVector->at(iTrack);
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t ForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceXYVector->at(iTrack);
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t ForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceXYErrorVector->at(iTrack);
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for normalized track chi2 value from reconstruction fit
Float_t ForestReader::GetTrackNormalizedChi2(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackNormalizedChi2Vector->at(iTrack);
  return GetTrackChi2(iTrack) / (1.0*GetNTrackDegreesOfFreedom(iTrack));
}

// Getter for track chi2 value from reconstruction fit
Float_t ForestReader::GetTrackChi2(Int_t iTrack) const{
  if(fIsMiniAOD) return -1; // Does not exist in MiniAOD forest
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t ForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  if(fIsMiniAOD) return -1; // Does not exist in MiniAOD forest
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
Int_t ForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  if(fIsMiniAOD) return fnHitsTrackerLayerVector->at(iTrack);
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
Int_t ForestReader::GetNHitsTrack(Int_t iTrack) const{
  if(fIsMiniAOD) return fnHitsTrackVector->at(iTrack);
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t ForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEnergyEcalVector->at(iTrack);
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t ForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEnergyHcalVector->at(iTrack);
  return fTrackEnergyHcalArray[iTrack];
}

// Getter for number of generator level particles
Int_t ForestReader::GetNGenParticles() const{
  return fnGenParticles;
}

// Getter for generator level particle pT
Float_t ForestReader::GetGenParticlePt(Int_t iTrack) const{
  return fGenParticlePtArray->at(iTrack);
}

// Getter for generator level particle phi
Float_t ForestReader::GetGenParticlePhi(Int_t iTrack) const{
  return fGenParticlePhiArray->at(iTrack);
}

// Getter for generator level particle eta
Float_t ForestReader::GetGenParticleEta(Int_t iTrack) const{
  return fGenParticleEtaArray->at(iTrack);
}

// Getter for generator level particle charge
Int_t ForestReader::GetGenParticleCharge(Int_t iTrack) const{
  return fGenParticleChargeArray->at(iTrack);
}

// Getter for generator level particle subevent index
Int_t ForestReader::GetGenParticleSubevent(Int_t iTrack) const{
  return fGenParticleSubeventArray->at(iTrack);
}
