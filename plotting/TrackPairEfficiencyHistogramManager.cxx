/*
 * Implementation of TrackPairEfficiencyHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "TrackPairEfficiencyHistogramManager.h"

/*
 * Default constructor
 */
TrackPairEfficiencyHistogramManager::TrackPairEfficiencyHistogramManager() :
  fInputFile(NULL),
  fCard(NULL),
  fSystemAndEnergy(""),
  fCompactSystemAndEnergy(""),
  fLoadEventInformation(false),
  fLoadJets(false),
  fLoad2DHistograms(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(1),
  fFirstLoadedTrackPtBin(0),
  fLastLoadedTrackPtBin(1),
  fFirstLoadedTrackPairPtBin(0),
  fLastLoadedTrackPairPtBin(1),
  fFirstLoadedJetPtBin(0),
  fLastLoadedJetPtBin(1),
  fFirstLoadedAverageEtaBin(0),
  fLastLoadedAverageEtaBin(1),
  fnCentralityBins(kMaxCentralityBins),
  fnTrackPtBins(kMaxTrackPtBins),
  fnTrackPairPtBins(kMaxTrackPtBins),
  fnJetPtBins(kMaxJetPtBinsEEC),
  fnAverageEtaBins(kMaxAverageEtaBins)
{
  
  // Do not load tracks by default
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = false;
  }
  
  // Do not load track pairs by default
  for(int iDataType = 0; iDataType < TrackPairEfficiencyHistograms::knDataLevels; iDataType++){
    fLoadTrackPairs[iDataType] = false;
    fLoadTrackPairsCloseToJets[iDataType] = false;
  }
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPtBinBorders[iTrackPt] = 0;
    fTrackPairPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPairPtBinBorders[iTrackPt] = 0;
  }

  // Default binning for jet pT
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC + 1; iJetPt++){
    fJetPtBinIndices[iJetPt] = iJetPt + 1;
    fJetPtBinBorders[iJetPt] = 0;
  }
  
  // Default binning for average pair eta
  for(int iAverageEta = 0; iAverageEta < kMaxAverageEtaBins + 1; iAverageEta++){
    fAverageEtaBinIndices[iAverageEta] = iAverageEta+1;
    fAverageEtaBinBorders[iAverageEta] = 0;
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;               // Vertex z position
  fhVertexZWeighted = NULL;       // Weighted vertex z-position (only meaningfull for MC)
  fhEvents = NULL;                // Number of events surviving different event cuts
  fhCentrality = NULL;            // Centrality of all events
  fhCentralityWeighted = NULL;    // Weighted centrality distribution in all events (only meaningful for MC)
  fhPtHat = NULL;                 // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;         // Weighted pT hat distribution (only meaningful for MC)
  fhTrackCuts = NULL;             // Number of tracks surviving different track cuts
  fhGenParticleSelections = NULL; // Number of generator level particles surviving different selections
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
      
      // Jet histograms
      fhJetPt[iCentrality][iDataLevel] = NULL;      // Jet pT histograms
      fhJetPhi[iCentrality][iDataLevel] = NULL;     // Jet phi histograms
      fhJetEta[iCentrality][iDataLevel] = NULL;     // Jet eta histograms
      fhJetEtaPhi[iCentrality][iDataLevel] = NULL;  // 2D eta-phi histogram for jets
      
    } // Data level loop
    
    // Track histograms
    for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
      fhTrackPt[iTrackType][iCentrality] = NULL;   // Track pT histograms
      
      // Loop over track pT bins
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = NULL;    // Track phi histograms
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = NULL;    // Track eta histograms
        fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = NULL; // 2D eta-phi histogram for track
      } // Track pT loop
      
    } // Track category loop
    
    // Track pair histograms
    for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
      for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins + 1; iTriggerPt++){
        for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins + 1; iAssociatedPt++){
          for(int iAverageEta = 0; iAverageEta < kMaxAverageEtaBins+1; iAverageEta++){
            fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel] = NULL;
          }
        } // Associated pT loop
      } // Trigger pT loop
    } // Data level loop

    // Track pair histograms close to jets
    for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
      for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
        for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins + 1; iTriggerPt++){
          for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins + 1; iAssociatedPt++){
            for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC+1; iJetPt++){
              fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = NULL;
            } // Jet pT loop
          } // Associated pT loop
        } // Trigger pT loop
      } // Data level loop for tracks
    } // Data level loop for jets
    
  } // Centrality loop
}

/*
 * Constructor
 */
TrackPairEfficiencyHistogramManager::TrackPairEfficiencyHistogramManager(TFile *inputFile) :
  TrackPairEfficiencyHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile
  fCard = new TrackPairEfficiencyCard(inputFile);
  
  // Initialize values using the information in card
  InitializeFromCard();
  
}

/*
 * Constructor
 */
TrackPairEfficiencyHistogramManager::TrackPairEfficiencyHistogramManager(TFile *inputFile, TrackPairEfficiencyCard *card) :
  TrackPairEfficiencyHistogramManager()
{
  fInputFile = inputFile;
  
  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Initialize several member variables from TrackPairEfficiencyCard
 */
void TrackPairEfficiencyHistogramManager::InitializeFromCard(){
  
  // Read the collision system from the card
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Read bins for centrality and track pT from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPtBins();
  fnTrackPairPtBins = fCard->GetNTrackPairPtBins();
  fnJetPtBins = fCard->GetNJetPtBins();
  fnAverageEtaBins = fCard->GetNAverageEtaBins();
  
  // Centrality binning
  for(int iCentrality = 0; iCentrality <= fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fLastLoadedCentralityBin = fnCentralityBins-1;
  
  // Track pT binning
  for(int iTrackPt = 0; iTrackPt <= fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPt(iTrackPt);
  }
  fLastLoadedTrackPtBin = fnTrackPtBins-1;
  
  // Track pT binning in track pair histograms
  for(int iTrackPt = 0; iTrackPt <= fnTrackPairPtBins; iTrackPt++){
    fTrackPairPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPairPt(iTrackPt);
  }
  fLastLoadedTrackPairPtBin = fnTrackPairPtBins-1;
  
  // Jet pT binning
  for(int iJetPt = 0; iJetPt <= fnJetPtBins; iJetPt++){
    fJetPtBinBorders[iJetPt] = fCard->GetLowBinBorderJetPt(iJetPt);
  }
  fLastLoadedJetPtBin = fnJetPtBins-1;

  // Average pair eta binning
  for(int iAverageEta = 0; iAverageEta <= fnAverageEtaBins; iAverageEta++){
    fAverageEtaBinBorders[iAverageEta] = fCard->GetLowBinBorderAverageEta(iAverageEta);
  }
  fLastLoadedAverageEtaBin = fnAverageEtaBins;
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    fLastLoadedCentralityBin = 0;
    fCentralityBinBorders[0] = -0.5;
  }
  
}

/*
 * Copy constructor
 */
TrackPairEfficiencyHistogramManager::TrackPairEfficiencyHistogramManager(const TrackPairEfficiencyHistogramManager& in) :
  fInputFile(in.fInputFile),
  fCard(in.fCard),
  fSystemAndEnergy(in.fSystemAndEnergy),
  fCompactSystemAndEnergy(in.fCompactSystemAndEnergy),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadJets(in.fLoadJets),
  fLoad2DHistograms(in.fLoad2DHistograms),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fFirstLoadedTrackPtBin(in.fFirstLoadedTrackPtBin),
  fLastLoadedTrackPtBin(in.fLastLoadedTrackPtBin),
  fFirstLoadedTrackPairPtBin(in.fFirstLoadedTrackPairPtBin),
  fLastLoadedTrackPairPtBin(in.fLastLoadedTrackPairPtBin),
  fFirstLoadedJetPtBin(in.fFirstLoadedJetPtBin),
  fLastLoadedJetPtBin(in.fLastLoadedJetPtBin),
  fFirstLoadedAverageEtaBin(in.fFirstLoadedAverageEtaBin),
  fLastLoadedAverageEtaBin(in.fLastLoadedAverageEtaBin),
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted)
{
  // Copy constructor
  
  // Do not load tracks by default
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = in.fLoadTracks[iTrackType];
  }
  
  // Do not load track pairs by default
  for(int iDataType = 0; iDataType < TrackPairEfficiencyHistograms::knDataLevels; iDataType++){
    fLoadTrackPairs[iDataType] = in.fLoadTrackPairs[iDataType];
    fLoadTrackPairsCloseToJets[iDataType] = in.fLoadTrackPairsCloseToJets[iDataType];
  }
  
  // Copy binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }
  
  // Copy binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = in.fTrackPtBinIndices[iTrackPt];
    fTrackPtBinBorders[iTrackPt] = in.fTrackPtBinBorders[iTrackPt];
    fTrackPairPtBinIndices[iTrackPt] = in.fTrackPairPtBinIndices[iTrackPt];
    fTrackPairPtBinBorders[iTrackPt] = in.fTrackPairPtBinBorders[iTrackPt];
  }

  // Copy binning for jet pT
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC + 1; iJetPt++){
    fJetPtBinIndices[iJetPt] = in.fJetPtBinIndices[iJetPt];
    fJetPtBinBorders[iJetPt] = in.fJetPtBinBorders[iJetPt];
  }
  
  // Copy binning for average pair eta
  for(int iAverageEta = 0; iAverageEta < kMaxAverageEtaBins + 1; iAverageEta++){
    fAverageEtaBinIndices[iAverageEta] = in.fAverageEtaBinIndices[iAverageEta];
    fAverageEtaBinBorders[iAverageEta] = in.fAverageEtaBinBorders[iAverageEta];
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
      
      // Jet histograms
      fhJetPt[iCentrality][iDataLevel] = in.fhJetPt[iCentrality][iDataLevel];         // Jet pT histograms
      fhJetPhi[iCentrality][iDataLevel] = in.fhJetPhi[iCentrality][iDataLevel];       // Jet phi histograms
      fhJetEta[iCentrality][iDataLevel] = in.fhJetEta[iCentrality][iDataLevel];       // Jet eta histograms
      fhJetEtaPhi[iCentrality][iDataLevel] = in.fhJetEtaPhi[iCentrality][iDataLevel]; // 2D eta-phi histogram for jets
      
    } // Data level loop
    
    // Track histograms
    for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
      fhTrackPt[iTrackType][iCentrality] = in.fhTrackPt[iTrackType][iCentrality];   // Track pT histograms
      
      // Loop over track pT bins
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = in.fhTrackPhi[iTrackType][iCentrality][iTrackPt];       // Track phi histograms
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = in.fhTrackEta[iTrackType][iCentrality][iTrackPt];       // Track eta histograms
        fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = in.fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt]; // 2D eta-phi histogram for track
      } // Track pT loop
      
    } // Track category loop
    
    // Track pair histograms
    for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
      for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins + 1; iTriggerPt++){
        for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins + 1; iAssociatedPt++){
          for(int iAverageEta = 0; iAverageEta < kMaxAverageEtaBins + 1; iAverageEta++){
            fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel] = in.fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel];
          }
        } // Associated pT loop
      } // Trigger pT loop
    } // Data level loop

    // Track pair histograms close to jets
    for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
      for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
        for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins + 1; iTriggerPt++){
          for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins + 1; iAssociatedPt++){
            for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC+1; iJetPt++){
              fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = in.fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt];
            } // Jet pT loop
          } // Associated pT loop
        } // Trigger pT loop
      } // Data level loop for tracks
    } // Data level loop for jets
    
  } // Centrality loop
}

/*
 * Destructor
 */
TrackPairEfficiencyHistogramManager::~TrackPairEfficiencyHistogramManager(){
  delete fCard;
}

/*
 * Load all the selected histograms from the inputfile
 */
void TrackPairEfficiencyHistogramManager::LoadHistograms(){
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                                // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                             // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");             // MC weighted vertex z position
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                       // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");       // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                                 // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");                 // Weighted pT hat for MC events
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                         // Number of tracks surviving different track cuts
    fhGenParticleSelections = (TH1D*) fInputFile->Get("genParticleSelections"); // Number of generator level particles surviving different selections
  }
  
  // Load jet histograms
  LoadJetHistograms();
  
  // Load track histograms
  LoadTrackHistograms();
  
  // Load track pair histograms
  LoadTrackPairHistograms();

  // Load track pair histograms close to jets
  LoadTrackPairHistogramsCloseToJets();
  
}

/*
 * Loader for jet histograms
 *
 * THnSparse for jets:
 *
 *   Histogram name: inclusiveJet
 *
 *     Axis index               Content of axis
 * --------------------------------------------------------
 *       Axis 0                     Jet pT
 *       Axis 1                     Jet phi
 *       Axis 2                     Jet eta
 *       Axis 3                    Centrality
 *       Axis 4     Data level (Reconstructed / Generator level)
 *       Axis 5                Trigger selection
 */
void TrackPairEfficiencyHistogramManager::LoadJetHistograms(){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  int nAxes = 2;           // Number of constraining axes for this iteration
  
  // Find the histogram array from which the projections are made
  THnSparseD *histogramArray = (THnSparseD*) fInputFile->Get(fJetHistogramName);
  
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
    
    // Select the bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemoverCentrality;
    
    axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;  // Centrality
    
    for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
      
      // Select the correct data level
      axisIndices[1] = 4; lowLimits[1] = iDataLevel+1; highLimits[1] = iDataLevel+1;
      
      // Always load jet pT histograms
      fhJetPt[iCentrality][iDataLevel] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);
      
      if(!fLoadJets) continue;  // Only load the remaining jet histograms if selected
      
      fhJetPhi[iCentrality][iDataLevel] = FindHistogram(histogramArray,1,nAxes,axisIndices,lowLimits,highLimits);
      fhJetEta[iCentrality][iDataLevel] = FindHistogram(histogramArray,2,nAxes,axisIndices,lowLimits,highLimits);
      if(fLoad2DHistograms) fhJetEtaPhi[iCentrality][iDataLevel] = FindHistogram2D(histogramArray,1,2,nAxes,axisIndices,lowLimits,highLimits);
      
    } // Data level loop
    
  } // Loop over centrality bins
}

/*
 * Loader for track histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: track/trackUncorrected/genParticle
 *
 *     Axis index       Content of axis
 * ----------------------------------------
 *       Axis 0            Track pT
 *       Axis 1            Track phi
 *       Axis 2            Track eta
 *       Axis 3            Centrality
 */
void TrackPairEfficiencyHistogramManager::LoadTrackHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Loop over all track histograms
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    
    // Find the histogram array from which the projections are made
    THnSparseD *histogramArray = (THnSparseD*) fInputFile->Get(fTrackHistogramNames[iTrackType]);
    
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Reset the ranges for all the axes in the histogram array
      for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
        histogramArray->GetAxis(iAxis)->SetRange(0,0);
      }
      
      // Select the bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemoverCentrality;
      
      // Setup axes with restrictions, (3 = centrality)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      
      fhTrackPt[iTrackType][iCentrality] = FindHistogram(histogramArray,0,1,axisIndices,lowLimits,highLimits);
      fhTrackPhi[iTrackType][iCentrality][fnTrackPtBins] = FindHistogram(histogramArray,1,1,axisIndices,lowLimits,highLimits);
      fhTrackEta[iTrackType][iCentrality][fnTrackPtBins] = FindHistogram(histogramArray,2,1,axisIndices,lowLimits,highLimits);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentrality][fnTrackPtBins] = FindHistogram2D(histogramArray,1,2,1,axisIndices,lowLimits,highLimits);
      
      for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
        
        // Select the bin indices for track pT
        lowerTrackPtBin = fTrackPtBinIndices[iTrackPt];
        higherTrackPtBin = fTrackPtBinIndices[iTrackPt+1]+duplicateRemoverTrackPt;
        
        // Add restriction for pT axis (0)
        axisIndices[1] = 0; lowLimits[1] = lowerTrackPtBin; highLimits[1] = higherTrackPtBin;
        
        // Read the angle histograms in track pT bins
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = FindHistogram(histogramArray,1,2,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = FindHistogram(histogramArray,2,2,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = FindHistogram2D(histogramArray,1,2,2,axisIndices,lowLimits,highLimits);
        
      } // Track pT loop
    } // Centrality loop
  } // Track category loop
}

/*
 * Loader for track pair histograms
 *
 * THnSparse for track pairs:
 *
 *   Histogram name: trackPairs/genParticlePairs
 *
 *     Axis index              Content of axis
 * -----------------------------------------------------
 *       Axis 0            DeltaR between particles
 *       Axis 1               Trigger particle pT
 *       Axis 2             Associated particle pT
 *       Axis 3                Average pair phi
 *       Axis 4                Average pair eta
 *       Axis 5                   Centrality
 */
void TrackPairEfficiencyHistogramManager::LoadTrackPairHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  int nProjectionAxes = 3;
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerTriggerPtBin = 0;
  int higherTriggerPtBin = 0;
  int lowerAssociatedPtBin = 0;
  int higherAssociatedPtBin = 0;
  int lowerAverageEtaBin = 0;
  int higherAverageEtaBin = 0;
  
  // Loop over all track pair histograms
  for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
    if(!fLoadTrackPairs[iDataLevel]) continue;  // Only load the selected track pair types
    
    // Find the histogram array from which the projections are made
    THnSparseD *histogramArray = (THnSparseD*) fInputFile->Get(fTrackPairHistogramNames[iDataLevel]);
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Select the bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;
      
      // Setup axes with restrictions, (5 = centrality)
      axisIndices[0] = 5; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      
      // Trigger pT loop
      for(int iTriggerPt = fFirstLoadedTrackPairPtBin; iTriggerPt <= fLastLoadedTrackPairPtBin; iTriggerPt++){
        
        // Select the bin indices for track pT
        lowerTriggerPtBin = fTrackPairPtBinIndices[iTriggerPt];
        higherTriggerPtBin = fTrackPairPtBinIndices[iTriggerPt+1]+duplicateRemover;
        
        // Add restriction for trigger pT axis (1)
        axisIndices[1] = 1; lowLimits[1] = lowerTriggerPtBin; highLimits[1] = higherTriggerPtBin;
        
        // Associated pT loop
        for(int iAssociatedPt = fFirstLoadedTrackPairPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          
          // Select the bin indices for track pT
          lowerAssociatedPtBin = fTrackPairPtBinIndices[iAssociatedPt];
          higherAssociatedPtBin = fTrackPairPtBinIndices[iAssociatedPt+1]+duplicateRemover;
          
          // Add restriction for associated pT axis (2)
          axisIndices[2] = 2; lowLimits[2] = lowerAssociatedPtBin; highLimits[2] = higherAssociatedPtBin;
          
          for(int iAverageEta = fFirstLoadedAverageEtaBin; iAverageEta <= fLastLoadedAverageEtaBin; iAverageEta++){
            
            // Add a new projection axis if we are not in the average eta integral bin index
            if(iAverageEta < fnAverageEtaBins){
              nProjectionAxes = 4;
              lowerAverageEtaBin = fAverageEtaBinIndices[iAverageEta];
              higherAverageEtaBin = fAverageEtaBinIndices[iAverageEta+1]+duplicateRemover;
              axisIndices[3] = 4; lowLimits[3] = lowerAverageEtaBin; highLimits[3] = higherAverageEtaBin;
            } else {
              // If we are not projecting, reset the range in the average eta axis of the histogram
              nProjectionAxes = 3;
              histogramArray->GetAxis(4)->SetRange(0,0);
            }
            
            // Read the deltaR distribution histograms
            fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel] = FindHistogram(histogramArray, 0, nProjectionAxes, axisIndices, lowLimits, highLimits);
            
          }
          
        } // Associated pT loop
      } // Trigger pT loop
    } // Centrality loop
  } // Data level loop
}

/*
 * Loader for track pair histograms with jet pT binning
 *
 * THnSparse for track pairs:
 *
 *   Histogram name: trackPairsCloseToJet/genParticlePairsCloseToJet
 *
 *     Axis index              Content of axis
 * --------------------------------------------------------
 *       Axis 0            DeltaR between particles
 *       Axis 1               Trigger particle pT
 *       Axis 2             Associated particle pT
 *       Axis 3                     Jet pT
 *       Axis 4        Reconstructed / generator level jet
 *       Axis 5                   Centrality
 */
void TrackPairEfficiencyHistogramManager::LoadTrackPairHistogramsCloseToJets(){
  
  // Define arrays to help find the histograms
  int axisIndices[5] = {0};
  int lowLimits[5] = {0};
  int highLimits[5] = {0};
  int nProjectionAxes = 4;
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerTriggerPtBin = 0;
  int higherTriggerPtBin = 0;
  int lowerAssociatedPtBin = 0;
  int higherAssociatedPtBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;

  // Loop over all track pair histograms close to jets
  for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
    if(!fLoadTrackPairsCloseToJets[iDataLevelTracks]) continue;  // Only load the selected track pair types

    // Find the histogram array from which the projections are made
    THnSparseD *histogramArray = (THnSparseD *)fInputFile->Get(fTrackPairHistogramCloseToJetNames[iDataLevelTracks]);

    // Loop over data level for jets close to which the tracks are paired
    for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){

      // Setup the jet data level axis restrictions
      axisIndices[0] = 4; lowLimits[0] = iDataLevelJets+1; highLimits[0] = iDataLevelJets+1;

      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // Select the bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentrality];
        higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemover;

        // Setup axes with restrictions, (5 = centrality)
        axisIndices[1] = 5; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;

        // Trigger pT loop
        for(int iTriggerPt = fFirstLoadedTrackPairPtBin; iTriggerPt <= fLastLoadedTrackPairPtBin; iTriggerPt++){

          // Select the bin indices for track pT
          lowerTriggerPtBin = fTrackPairPtBinIndices[iTriggerPt];
          higherTriggerPtBin = fTrackPairPtBinIndices[iTriggerPt+1] + duplicateRemover;

          // Add restriction for trigger pT axis (1)
          axisIndices[2] = 1; lowLimits[2] = lowerTriggerPtBin; highLimits[2] = higherTriggerPtBin;

          // Associated pT loop
          for(int iAssociatedPt = fFirstLoadedTrackPairPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){

            // Select the bin indices for track pT
            lowerAssociatedPtBin = fTrackPairPtBinIndices[iAssociatedPt];
            higherAssociatedPtBin = fTrackPairPtBinIndices[iAssociatedPt+1] + duplicateRemover;

            // Add restriction for associated pT axis (2)
            axisIndices[3] = 2; lowLimits[3] = lowerAssociatedPtBin; highLimits[3] = higherAssociatedPtBin;

            for(int iJetPt = fFirstLoadedJetPtBin; iJetPt <= fLastLoadedJetPtBin; iJetPt++) {

              // Add a new projection axis if we are not averaging over all jet pT:s
              if(iJetPt < fnJetPtBins) {
                nProjectionAxes = 5;
                lowerJetPtBin = fJetPtBinIndices[iJetPt];
                higherJetPtBin = fJetPtBinIndices[iJetPt+1] + duplicateRemover;
                axisIndices[4] = 3; lowLimits[4] = lowerJetPtBin; highLimits[4] = higherJetPtBin;
              } else {
                // If we do not care about the jet pT, reset the range in the jet pT axis of the histogram
                nProjectionAxes = 4;
                histogramArray->GetAxis(3)->SetRange(0,0);
              }

              // Read the deltaR distribution histograms
              fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = FindHistogram(histogramArray, 0, nProjectionAxes, axisIndices, lowLimits, highLimits);
            }

          } // Associated pT loop
        } // Trigger pT loop
      } // Centrality loop
    } // Jet data level loop
  } // Track data level loop
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD *histogramArray = Inputfile containing the THnSparse to be read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* TrackPairEfficiencyHistogramManager::FindHistogram2D(THnSparseD *histogramArray, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth){
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName.Data());
  
  // Apply bin width normalization to the projected histogram
  if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* TrackPairEfficiencyHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName.Data());
  
  // Apply bin width normalization to the projected histogram
  if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* TrackPairEfficiencyHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram2D(inputFile,name,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   THnSparseD *histogramArray = Histogram array from which the desired histograms are projected
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* TrackPairEfficiencyHistogramManager::FindHistogram(THnSparseD *histogramArray, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth){
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName.Data());
  
    // Apply bin width normalization to the projected histogram
    if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* TrackPairEfficiencyHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex, const bool normalizeToBinWidth){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName.Data());
  
    // Apply bin width normalization to the projected histogram
    if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* TrackPairEfficiencyHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram(inputFile,name,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Write all the loaded histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void TrackPairEfficiencyHistogramManager::Write(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile *outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString histogramNamer;
  
  // Write the event information histograms to the output file
  if(fLoadEventInformation){
    fhEvents->Write("",TObject::kOverwrite);                // Number of events surviving different event cuts
    fhVertexZ->Write("",TObject::kOverwrite);               // Vertex z position
    fhVertexZWeighted->Write("",TObject::kOverwrite);       // MC weighted vertex z position
    fhCentrality->Write("",TObject::kOverwrite);            // Centrality in all events
    fhCentralityWeighted->Write("",TObject::kOverwrite);    // MC weighted centrality in all events
    fhPtHat->Write("",TObject::kOverwrite);                 // pT hat for MC events (only meaningful for MC)
    fhPtHatWeighted->Write("",TObject::kOverwrite);         // Weighted pT hat distribution (only meaningful for MC)
    fhTrackCuts->Write("",TObject::kOverwrite);             // Number of tracks surviving different track cuts
    fhGenParticleSelections->Write("",TObject::kOverwrite); // Number of generator level particles surviving different selections
  }
 
  // Write the jet histograms to the output file
  WriteJetHistograms();
  
  // Write the track histograms to the output file
  WriteTrackHistograms();
  
  // Write the track pair histograms to the output file
  WriteTrackPairHistograms();

  // Write the track pair histograms close to jets to the output file
  WriteTrackPairHistogramsCloseToJets();
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the jet histograms to the file that is currently open
 */
void TrackPairEfficiencyHistogramManager::WriteJetHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the jet histograms to the output file
  if(!fLoadJets) return;  // Only write the jet histograms if they are loaded
  
  for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
    
    // Create a directory for the histograms if it does not already exist
    histogramNamer = Form("%s%s", fJetHistogramName, fDataLevelName[iDataLevel]);
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Check that the histograms are actually there before trying to save them.
      if(fhJetPt[iCentrality][iDataLevel] == NULL) {
        cout << "Could not find histograms of type " << fJetHistogramName << " to write. Will skip writing these." << endl;
        continue;
      }
      
      // Jet pT
      histogramNamer = Form("%sPt_C%d",fJetHistogramName, iCentrality);
      if(fhJetPt[iCentrality][iDataLevel]) fhJetPt[iCentrality][iDataLevel]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // Jet phi
      histogramNamer = Form("%sPhi_C%d",fJetHistogramName, iCentrality);
      if(fhJetPhi[iCentrality][iDataLevel]) fhJetPhi[iCentrality][iDataLevel]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // Jet eta
      histogramNamer = Form("%sEta_C%d",fJetHistogramName, iCentrality);
      if(fhJetEta[iCentrality][iDataLevel]) fhJetEta[iCentrality][iDataLevel]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // Jet eta-phi
      histogramNamer = Form("%sEtaPhi_C%d",fJetHistogramName, iCentrality);
      if(fLoad2DHistograms && fhJetEtaPhi[iCentrality][iDataLevel]) fhJetEtaPhi[iCentrality][iDataLevel]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
    } // Loop over centrality bins
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Data level loop
  
}

/*
 * Write the track histograms to the file that is currently open
 */
void TrackPairEfficiencyHistogramManager::WriteTrackHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the track histograms to the output file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only write the loaded track types
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fTrackHistogramNames[iTrackType])) gDirectory->mkdir(fTrackHistogramNames[iTrackType]);
    gDirectory->cd(fTrackHistogramNames[iTrackType]);
    
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Track pT
      histogramNamer = Form("%sPt_C%d", fTrackHistogramNames[iTrackType], iCentrality);
      fhTrackPt[iTrackType][iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // pT integrated track phi
      histogramNamer = Form("%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentrality, fnTrackPtBins);
      fhTrackPhi[iTrackType][iCentrality][fnTrackPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // pT integrated track eta
      histogramNamer = Form("%sEta_C%dT%d", fTrackHistogramNames[iTrackType], iCentrality, fnTrackPtBins);
      fhTrackEta[iTrackType][iCentrality][fnTrackPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // pT integrated track eta-phi
      histogramNamer = Form("%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentrality, fnTrackPtBins);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentrality][fnTrackPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
        
        // Track phi in track pT bins
        histogramNamer = Form("%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentrality, iTrackPt);
        fhTrackPhi[iTrackType][iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
        // Track eta in track pT bins
        histogramNamer = Form("%sEta_C%dT%d", fTrackHistogramNames[iTrackType], iCentrality, iTrackPt);
        fhTrackEta[iTrackType][iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
        // Track eta-phi in track pT bins
        histogramNamer = Form("%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentrality, iTrackPt);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Track category loop

}

/*
 * Write the track pair histograms to the file that is currently open
 */
void TrackPairEfficiencyHistogramManager::WriteTrackPairHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the track histograms to the output file
  for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
    if(!fLoadTrackPairs[iDataLevel]) continue;  // Only write the loaded track pair types
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fTrackPairHistogramNames[iDataLevel])) gDirectory->mkdir(fTrackPairHistogramNames[iDataLevel]);
    gDirectory->cd(fTrackPairHistogramNames[iDataLevel]);
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Trigger pT loop
      for(int iTriggerPt = fFirstLoadedTrackPairPtBin; iTriggerPt <= fLastLoadedTrackPairPtBin; iTriggerPt++){
        
        // Associated pT loop
        for(int iAssociatedPt = fFirstLoadedTrackPairPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          
          // DeltaR histograms in average pair eta bins
          for(int iAverageEta = fFirstLoadedAverageEtaBin; iAverageEta <= fLastLoadedAverageEtaBin; iAverageEta++){
            
            if(iAverageEta == fnAverageEtaBins){
              // DeltaR histograms without eta or phi selections
              histogramNamer = Form("%sDeltaR_C%dT%dA%d", fTrackPairHistogramNames[iDataLevel], iCentrality, iTriggerPt, iAssociatedPt);
            } else {
              // DeltaR histograms in average eta bins
              histogramNamer = Form("%sDeltaR_C%dT%dA%dE%d", fTrackPairHistogramNames[iDataLevel], iCentrality, iTriggerPt, iAssociatedPt, iAverageEta);
            }
            fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel]->Write(histogramNamer.Data(), TObject::kOverwrite);
          }
          
        } // Associated pT loop
      } // Trigger pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Data level loop

}

/*
 * Write the track pair histograms close to jets to the file that is currently open
 */
void TrackPairEfficiencyHistogramManager::WriteTrackPairHistogramsCloseToJets(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the track histograms to the output file
  for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
    if(!fLoadTrackPairsCloseToJets[iDataLevelTracks]) continue;  // Only write the loaded track pair types
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fTrackPairHistogramCloseToJetNames[iDataLevelTracks])) gDirectory->mkdir(fTrackPairHistogramCloseToJetNames[iDataLevelTracks]);
    gDirectory->cd(fTrackPairHistogramCloseToJetNames[iDataLevelTracks]);

    // Loop over jet data level (reconstructed / generator level)
    for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
    
      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
        // Trigger pT loop
        for(int iTriggerPt = fFirstLoadedTrackPairPtBin; iTriggerPt <= fLastLoadedTrackPairPtBin; iTriggerPt++){
        
          // Associated pT loop
          for(int iAssociatedPt = fFirstLoadedTrackPairPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          
            // DeltaR histograms in average pair eta bins
            for(int iJetPt = fFirstLoadedJetPtBin; iJetPt <= fLastLoadedJetPtBin; iJetPt++){
            
              if(iJetPt == fnJetPtBins){
                // DeltaR histograms without jet pT selection
                histogramNamer = Form("%sDeltaRCloseTo%sJets_C%dT%dA%d", fTrackPairHistogramNames[iDataLevelTracks], fDataLevelName[iDataLevelJets], iCentrality, iTriggerPt, iAssociatedPt);
              } else {
                // DeltaR histograms in jet pT bins
                histogramNamer = Form("%sDeltaRCloseTo%sJets_C%dT%dA%dJ%d", fTrackPairHistogramNames[iDataLevelTracks], fDataLevelName[iDataLevelJets], iCentrality, iTriggerPt, iAssociatedPt, iJetPt);
              }
              fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
            }
          
          } // Associated pT loop
        } // Trigger pT loop
      } // Centrality loop
    } // Data level for jets loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Data level loop

}

/*
 * Load the selected histograms from a file containing readily processed histograms
 */
void TrackPairEfficiencyHistogramManager::LoadProcessedHistograms(){
  
  // Helper variable for finding names of loaded histograms
  TString histogramNamer;
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                                // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                             // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");             // MC weighted vertex z position
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                       // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");       // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                                 // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");                 // Weighted pT hat for MC events
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                         // Number of tracks surviving different track cuts
    fhGenParticleSelections = (TH1D*) fInputFile->Get("genParticleSelections"); // Number of generator level particles surviving different selections
  }
  
  // Load the jet histograms from the input file
  for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Always load jet pT histograms
      histogramNamer = Form("%s%s/%sPt_C%d", fJetHistogramName, fDataLevelName[iDataLevel], fJetHistogramName, iCentrality);
      fhJetPt[iCentrality][iDataLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      if(!fLoadJets) continue;  // Only load the loaded the selected histograms
      
      // Jet phi
      histogramNamer = Form("%s%s/%sPhi_C%d", fJetHistogramName, fDataLevelName[iDataLevel], fJetHistogramName, iCentrality);
      fhJetPhi[iCentrality][iDataLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // Jet eta
      histogramNamer = Form("%s%s/%sEta_C%d", fJetHistogramName, fDataLevelName[iDataLevel], fJetHistogramName, iCentrality);
      fhJetEta[iCentrality][iDataLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // Jet eta-phi
      histogramNamer = Form("%s%s/%sEtaPhi_C%d", fJetHistogramName, fDataLevelName[iDataLevel], fJetHistogramName, iCentrality);
      if(fLoad2DHistograms) fhJetEtaPhi[iCentrality][iDataLevel] = (TH2D*) fInputFile->Get(histogramNamer.Data());
      
    } // Loop over centrality bins
  } // Data level loop
  
  // Load the track histograms from the input file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Track pT
      histogramNamer = Form("%s/%sPt_C%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality);
      fhTrackPt[iTrackType][iCentrality] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // pT integrated track phi
      histogramNamer = Form("%s/%sPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality,fnTrackPtBins);
      fhTrackPhi[iTrackType][iCentrality][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // pT integrated track eta
      histogramNamer = Form("%s/%sEta_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality,fnTrackPtBins);
      fhTrackEta[iTrackType][iCentrality][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // pT integrated track eta-phi
      histogramNamer = Form("%s/%sEtaPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality,fnTrackPtBins);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentrality][fnTrackPtBins] = (TH2D*) fInputFile->Get(histogramNamer.Data());
      
      for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
        
        // Track phi in track pT bins
        histogramNamer = Form("%s/%sPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality,iTrackPt);
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
        
        // Track eta in track pT bins
        histogramNamer = Form("%s/%sEta_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality,iTrackPt);
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
        
        // Track eta-phi in track pT bins
        histogramNamer = Form("%s/%sEtaPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentrality,iTrackPt);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = (TH2D*) fInputFile->Get(histogramNamer.Data());
        
      } // Track pT loop
    } // Centrality loop
  } // Track category loop
  
  // Load the track pair histograms from the input file
  for(int iDataLevel = 0; iDataLevel < TrackPairEfficiencyHistograms::knDataLevels; iDataLevel++){
    if(!fLoadTrackPairs[iDataLevel]) continue;  // Only load the selected track pair types
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTriggerPt = fFirstLoadedTrackPairPtBin; iTriggerPt <= fLastLoadedTrackPairPtBin; iTriggerPt++){
        for(int iAssociatedPt = fFirstLoadedTrackPairPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
          for(int iAverageEta = fFirstLoadedAverageEtaBin; iAverageEta <= fLastLoadedAverageEtaBin; iAverageEta++){
            // DeltaR histograms without eta or phi selection
            if(iAverageEta == fnAverageEtaBins){
              histogramNamer = Form("%s/%sDeltaR_C%dT%dA%d", fTrackPairHistogramNames[iDataLevel], fTrackPairHistogramNames[iDataLevel], iCentrality, iTriggerPt, iAssociatedPt);
            } else {
              histogramNamer = Form("%s/%sDeltaR_C%dT%dA%dE%d", fTrackPairHistogramNames[iDataLevel], fTrackPairHistogramNames[iDataLevel], iCentrality, iTriggerPt, iAssociatedPt, iAverageEta);
            }
            fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          } // Average pair eta loop
        } // Associated pT loop
      } // Trigger pT loop
    } // Centrality loop
  } // Data level loop

  // Load the track pair histograms close to jets from the input file
  for(int iDataLevelTracks = 0; iDataLevelTracks < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelTracks++){
    if(!fLoadTrackPairsCloseToJets[iDataLevelTracks]) continue;  // Only load the selected track pair types
    for(int iDataLevelJets = 0; iDataLevelJets < TrackPairEfficiencyHistograms::knDataLevels; iDataLevelJets++){
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iTriggerPt = fFirstLoadedTrackPairPtBin; iTriggerPt <= fLastLoadedTrackPairPtBin; iTriggerPt++){
          for(int iAssociatedPt = fFirstLoadedTrackPairPtBin; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
            for(int iJetPt = fFirstLoadedJetPtBin; iJetPt <= fLastLoadedJetPtBin; iJetPt++){
              // DeltaR histograms without eta or phi selection
              if(iJetPt == fnJetPtBins){
                histogramNamer = Form("%s/%sDeltaRCloseTo%sJets_C%dT%dA%d", fTrackPairHistogramCloseToJetNames[iDataLevelTracks], fTrackPairHistogramNames[iDataLevelTracks], fDataLevelName[iDataLevelJets], iCentrality, iTriggerPt, iAssociatedPt);
              } else {
                histogramNamer = Form("%s/%sDeltaRCloseTo%sJets_C%dT%dA%dJ%d", fTrackPairHistogramCloseToJetNames[iDataLevelTracks], fTrackPairHistogramNames[iDataLevelTracks], fDataLevelName[iDataLevelJets], iCentrality, iTriggerPt, iAssociatedPt, iJetPt);
              }
              fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
            } // Average pair eta loop
          } // Associated pT loop
        } // Trigger pT loop
      } // Centrality loop
    } // Data level loop for jets
  } // Data level loop for tracks
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void TrackPairEfficiencyHistogramManager::SetBinIndices(const char* histogramName, const int nBins, int *binIndices, const double *binBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   double *copyBinBorders = Array to which a copy of bin borders is made
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void TrackPairEfficiencyHistogramManager::SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices){
  TH1D* hBinner;
  if(setIndices) hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    if(setIndices) binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Set up generic bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const char* histogramName = Name of the histogram from which the bin indices are searched
 *  const int iAxis = Axis from which the set indices can be found
 *  int nSetBins = Number of bins that is set
 *  double* setBinBorders = Bin borders that are set
 *  int* setBinIndices = Bin indices that are set
 *  const int nBins = New number of bins that is given
 *  const double *binBorders = New bin borders that are given
 *  const char* errorMessage = Type of the set bins to be printed in possible error message
 *  const int maxBins = Maximum number of allowed bins of this type
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TrackPairEfficiencyHistogramManager::SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double *binBorders, const char* errorMessage, const int maxBins, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    if(setIndices) SetBinIndices(histogramName, nSetBins, setBinIndices, setBinBorders, iAxis);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= maxBins){
      nSetBins = nBins;
      SetBinBordersAndIndices(histogramName, nSetBins, setBinBorders, setBinIndices, binBorders, iAxis, setIndices);
    } else {
      cout << "Error! Too many " << errorMessage << " bins given. Maximum number is " << maxBins << ". Will not set bins." << endl;
    }
  }
}

/*
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given centrality bins
 *  const double *binBorders = New bin borders for centrality
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TrackPairEfficiencyHistogramManager::SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, fJetHistogramName, 3, fnCentralityBins, fCentralityBinBorders, fCentralityBinIndices, nBins, binBorders, "centrality", kMaxCentralityBins, setIndices);
  
}

/*
 * Set up track pT bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TrackPairEfficiencyHistogramManager::SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "track", 0, fnTrackPtBins, fTrackPtBinBorders, fTrackPtBinIndices, nBins, binBorders, "track pT", kMaxTrackPtBins, setIndices);
  
}

/*
 * Set up track pT bin borders and indices according to provided bin borders for track pair histograms
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TrackPairEfficiencyHistogramManager::SetTrackPairPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "trackPairs", 1, fnTrackPairPtBins, fTrackPairPtBinBorders, fTrackPairPtBinIndices, nBins, binBorders, "track pair pT", kMaxTrackPtBins, setIndices);
  
}

/*
 * Set up jet pT bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TrackPairEfficiencyHistogramManager::SetJetPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "trackPairsCloseToJet", 3, fnJetPtBins, fJetPtBinBorders, fJetPtBinIndices, nBins, binBorders, "track pT", kMaxJetPtBinsEEC, setIndices);
  
}

/*
 * Set up track pT bin borders and indices according to provided bin borders for track pair histograms
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void TrackPairEfficiencyHistogramManager::SetAverageEtaBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "trackPairs", 4, fnAverageEtaBins, fAverageEtaBinBorders, fAverageEtaBinIndices, nBins, binBorders, "average eta", kMaxAverageEtaBins, setIndices);
  
}

// Setter for loading event information
void TrackPairEfficiencyHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for loading jet histograms
void TrackPairEfficiencyHistogramManager::SetLoadJetHistograms(const bool loadOrNot){
  fLoadJets = loadOrNot;
}

// Setter for loading tracks
void TrackPairEfficiencyHistogramManager::SetLoadTracks(const bool loadOrNot){
  fLoadTracks[kTrack] = loadOrNot;
}

// Setter for loading uncorrected tracks
void TrackPairEfficiencyHistogramManager::SetLoadTracksUncorrected(const bool loadOrNot){
  fLoadTracks[kUncorrectedTrack] = loadOrNot;
}

// Setter for loading generator level particles
void TrackPairEfficiencyHistogramManager::SetLoadGenParticles(const bool loadOrNot){
  fLoadTracks[kGenParticle] = loadOrNot;
}

// Setter for loading all track histograms
void TrackPairEfficiencyHistogramManager::SetLoadAllTracks(const bool loadTracks, const bool loadUncorrected, const bool loadGenParticles){
  SetLoadTracks(loadTracks);
  SetLoadTracksUncorrected(loadUncorrected);
  SetLoadGenParticles(loadGenParticles);
}

// Setter for loading track pairs
void TrackPairEfficiencyHistogramManager::SetLoadTrackPairs(const bool loadOrNot){
  fLoadTrackPairs[TrackPairEfficiencyHistograms::kReconstructed] = loadOrNot;
}

// Setter for loading generator level particle pairs
void TrackPairEfficiencyHistogramManager::SetLoadGenParticlePairs(const bool loadOrNot){
  fLoadTrackPairs[TrackPairEfficiencyHistograms::kGeneratorLevel] = loadOrNot;
}

// Setter for loading all track pair histograms
void TrackPairEfficiencyHistogramManager::SetLoadAllTrackPairs(const bool loadTracks, const bool loadGenParticles){
  SetLoadTrackPairs(loadTracks);
  SetLoadGenParticlePairs(loadGenParticles);
}

// Setter for loading track pairs close to jets
void TrackPairEfficiencyHistogramManager::SetLoadTrackPairsCloseToJets(const bool loadOrNot){
  fLoadTrackPairsCloseToJets[TrackPairEfficiencyHistograms::kReconstructed] = loadOrNot;
}

// Setter for loading generator level particle pairs close to jets
void TrackPairEfficiencyHistogramManager::SetLoadGenParticlePairsCloseToJets(const bool loadOrNot){
  fLoadTrackPairsCloseToJets[TrackPairEfficiencyHistograms::kGeneratorLevel] = loadOrNot;
}

// Setter for loading all track pair histograms close to jets
void TrackPairEfficiencyHistogramManager::SetLoadAllTrackPairsCloseToJets(const bool loadTracks, const bool loadGenParticles){
  SetLoadTrackPairsCloseToJets(loadTracks);
  SetLoadGenParticlePairsCloseToJets(loadGenParticles);
}

 // Setter for loading two-dimensional histograms
void TrackPairEfficiencyHistogramManager::SetLoad2DHistograms(const bool loadOrNot){
  fLoad2DHistograms = loadOrNot;
}

// Setter for loaded centrality bins
void TrackPairEfficiencyHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for centrality bins
  BinSanityCheck(fnCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Setter for loaded track pT bins
void TrackPairEfficiencyHistogramManager::SetTrackPtBinRange(const int first, const int last){
  fFirstLoadedTrackPtBin = first;
  fLastLoadedTrackPtBin = last;
  
  // Sanity check for track pT bins
  BinSanityCheck(fnTrackPtBins,fFirstLoadedTrackPtBin,fLastLoadedTrackPtBin);
}

// Setter for loaded track pair pT bins
void TrackPairEfficiencyHistogramManager::SetTrackPairPtBinRange(const int first, const int last){
  fFirstLoadedTrackPairPtBin = first;
  fLastLoadedTrackPairPtBin = last;
  
  // Sanity check for track pT bins
  BinSanityCheck(fnTrackPairPtBins,fFirstLoadedTrackPairPtBin,fLastLoadedTrackPairPtBin);
}

// Setter for loaded jet pT bins
void TrackPairEfficiencyHistogramManager::SetJetPtBinRange(const int first, const int last){
  fFirstLoadedJetPtBin = first;
  fLastLoadedJetPtBin = last;
  
  // Sanity check for track pT bins
  BinSanityCheck(fnJetPtBins+1,fFirstLoadedJetPtBin,fLastLoadedJetPtBin);
}

// Setter for loaded average pair eta bins
void TrackPairEfficiencyHistogramManager::SetAverageEtaBinRange(const int first, const int last){
  fFirstLoadedAverageEtaBin = first;
  fLastLoadedAverageEtaBin = last;
  
  // Sanity check for average pair eta bins
  BinSanityCheck(fnAverageEtaBins+1,fFirstLoadedAverageEtaBin,fLastLoadedAverageEtaBin);
}

// Sanity check for set bins
void TrackPairEfficiencyHistogramManager::BinSanityCheck(const int nBins, int& first, int& last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Sanity check for input bin index
int TrackPairEfficiencyHistogramManager::BinIndexCheck(const int nBins, const int binIndex) const{
  if(binIndex < 0) return 0;
  if(binIndex > nBins-1) return nBins-1;
  return binIndex;
}

// Getter for the number of centrality bins
int TrackPairEfficiencyHistogramManager::GetNCentralityBins() const{
  return fnCentralityBins;
}

// Getter for the number of track pT bins
int TrackPairEfficiencyHistogramManager::GetNTrackPtBins() const{
  return fnTrackPtBins;
}

// Getter for the number of track pT bins in track pair histograms
int TrackPairEfficiencyHistogramManager::GetNTrackPairPtBins() const{
  return fnTrackPairPtBins;
}

// Getter for the number of jet pT bins
int TrackPairEfficiencyHistogramManager::GetNJetPtBins() const{
  return fnJetPtBins;
}

// Getter for the number of average eta bins in track pair histograms
int TrackPairEfficiencyHistogramManager::GetNAverageEtaBins() const{
  return fnAverageEtaBins;
}

// Getter for the jet histogram name
const char* TrackPairEfficiencyHistogramManager::GetJetHistogramName() const{
  return fJetHistogramName;
}

// Getter for name suitable for x-axis in a given jet histogram
const char* TrackPairEfficiencyHistogramManager::GetJetAxisName() const{
  return fJetAxisName;
}

// Getter for the name of data level
const char* TrackPairEfficiencyHistogramManager::GetDataLevelName(const int iDataLevel) const{
  if(iDataLevel < 0 || iDataLevel >= TrackPairEfficiencyHistograms::knDataLevels) return "Unknown";
  return fDataLevelName[iDataLevel];
}

// Getter for collision system
TString TrackPairEfficiencyHistogramManager::GetSystem() const{
  return fCard->GetDataType();
}

// Getter for i:th centrality bin border
double TrackPairEfficiencyHistogramManager::GetCentralityBinBorder(const int iCentrality) const{
  return fCentralityBinBorders[iCentrality];
}

// Getter for i:th track pT bin border
double TrackPairEfficiencyHistogramManager::GetTrackPtBinBorder(const int iTrackPt) const{
  return fTrackPtBinBorders[iTrackPt];
}

// Getter for i:th track pT bin border in track pair histograms
double TrackPairEfficiencyHistogramManager::GetTrackPairPtBinBorder(const int iTrackPt) const{
  return fTrackPairPtBinBorders[iTrackPt];
}

// Getter for i:th jet pT bin border
double TrackPairEfficiencyHistogramManager::GetJetPtBinBorder(const int iJetPt) const{
  return fJetPtBinBorders[iJetPt];
}

// Getter for i:th average eta bin border in track pair histograms
double TrackPairEfficiencyHistogramManager::GetAverageEtaBinBorder(const int iAverageEta) const{
  return fAverageEtaBinBorders[iAverageEta];
}

// Getters for event information histograms

// Getter for z-vertex histogram
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramVertexZ() const{
  return fhVertexZ;
}

// Getter for z-vertex histogram
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramVertexZWeighted() const{
  return fhVertexZWeighted;
}

// Getter for histogram for number of events surviving different event cuts
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramEvents() const{
  return fhEvents;
}

// Getter for centrality histogram in all events
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramCentrality() const{
  return fhCentrality;
}

// Getter for weighted centrality histogram in all events
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramCentralityWeighted() const{
  return fhCentralityWeighted;
}

// Getters for jet histograms

// Getter for jet pT histograms
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramJetPt(int iCentrality, const int iDataLevel) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetPt[iCentrality][iDataLevel];
}

// Getter for jet phi histograms
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramJetPhi(int iCentrality, const int iDataLevel) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetPhi[iCentrality][iDataLevel];
}

// Getter for jet eta histograms
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramJetEta(int iCentrality, const int iDataLevel) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetEta[iCentrality][iDataLevel];
}

// Getter for 2D eta-phi histogram for jets
TH2D* TrackPairEfficiencyHistogramManager::GetHistogramJetEtaPhi(int iCentrality, const int iDataLevel) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetEtaPhi[iCentrality][iDataLevel];
}

// Getters for track histograms

// Getter for track pT histograms
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramTrackPt(const int iTrackType, const int iCentrality) const{
  return fhTrackPt[iTrackType][iCentrality];
}

// Getter for track phi histograms
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramTrackPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const{
  return fhTrackPhi[iTrackType][iCentrality][iTrackPt];
}

// Getter for track eta histograms
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramTrackEta(const int iTrackType, const int iCentrality, const int iTrackPt) const{
  return fhTrackEta[iTrackType][iCentrality][iTrackPt];
}

// Getter for 2D eta-phi histogram for track
TH2D* TrackPairEfficiencyHistogramManager::GetHistogramTrackEtaPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const{
  return fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt];
}

// Getters for track pair histograms

// Getter for DeltaR between tracks in pT and centrality bins
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramTrackPairDeltaR(const int iCentrality, const int iTriggerPt, const int iAssociatedPt, const int iAverageEta, const int iDataLevel) const{
  return fhTrackPairDeltaR[iCentrality][iTriggerPt][iAssociatedPt][iAverageEta][iDataLevel];
}

// Getter for DeltaR between tracks close to jets in pT and centrality bins
TH1D* TrackPairEfficiencyHistogramManager::GetHistogramTrackPairDeltaRCloseToJets(const int iDataLevelJets, const int iDataLevelTracks, const int iCentrality, const int iTriggerPt, const int iAssociatedPt, const int iJetPt) const{
  return fhTrackPairDeltaRCloseToJets[iDataLevelJets][iDataLevelTracks][iCentrality][iTriggerPt][iAssociatedPt][iJetPt];
}

// The rest

// Get the first loaded centrality bin
int TrackPairEfficiencyHistogramManager::GetFirstCentralityBin() const{
  return fFirstLoadedCentralityBin;
}

// Get the last loaded centrality bin
int TrackPairEfficiencyHistogramManager::GetLastCentralityBin() const{
  return fLastLoadedCentralityBin;
}

// Getter for the number of events passing the cuts
int TrackPairEfficiencyHistogramManager::GetNEvents() const{
  return fhEvents->GetBinContent(fhEvents->FindBin(TrackPairEfficiencyHistograms::kVzCut));
}

// Getter for the JCard
TrackPairEfficiencyCard* TrackPairEfficiencyHistogramManager::GetCard() const{
  return fCard;
}

// Getter for integral over inclusive jet pT. Include the overflow bin in the integral.
double TrackPairEfficiencyHistogramManager::GetJetPtIntegral(const int iCentrality, const int iDataLevel) const{
  return fhJetPt[iCentrality][iDataLevel]->Integral(1,fhJetPt[iCentrality][iDataLevel]->GetNbinsX()+1,"width");
}

/*
 * Getter for integral over inclusive jet pT over specified range
 *
 *  const int iCentrality = Centrality bin
 *  const int iDataLevel = Select reconstructed or generator level jets
 *  const double minPt = Lower pT range for integral calculation
 *  const double maxPt = Higher pT range for integral calculation
 */
double TrackPairEfficiencyHistogramManager::GetJetPtIntegral(const int iCentrality, const int iDataLevel, const double minPt, const double maxPt) const{
  return fhJetPt[iCentrality][iDataLevel]->Integral(fhJetPt[iCentrality][iDataLevel]->FindBin(minPt+0.001), fhJetPt[iCentrality][iDataLevel]->FindBin(maxPt-0.001), "width");
}
