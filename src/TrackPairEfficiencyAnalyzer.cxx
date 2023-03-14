// Class for the main analysis algorithms for the energy-energy correlator analysis

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "TrackPairEfficiencyAnalyzer.h"

using namespace std;

/*
 * Default constructor
 */
TrackPairEfficiencyAnalyzer::TrackPairEfficiencyAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fVzWeightFunction(0),
  fCentralityWeightFunctionCentral(0),
  fCentralityWeightFunctionPeripheral(0),
  fPtWeightFunction(0),
  fTrackEfficiencyCorrector2018(),
  fDataType(-1),
  fJetType(0),
  fUseTrigger(false),
  fDebugLevel(0),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fJetAxis(0),
  fVzCut(0),
  fMinimumPtHat(0),
  fMaximumPtHat(0),
  fJetEtaCut(0),
  fJetMinimumPtCut(0),
  fJetMaximumPtCut(0),
  fCutBadPhiRegion(false),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fTrackEtaCut(0),
  fTrackMinPtCut(0),
  fTrackMaxPtCut(0),
  fMaxTrackPtRelativeError(0),
  fMaxTrackDistanceToVertex(0),
  fCalorimeterSignalLimitPt(0),
  fHighPtEtFraction(0),
  fChi2QualityCut(0),
  fMinimumTrackHits(0),
  fSubeventCut(0)
{
  // Default constructor
  fHistograms = new TrackPairEfficiencyHistograms();
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fEventReader = NULL;
  
}

/*
 * Custom constructor
 */
TrackPairEfficiencyAnalyzer::TrackPairEfficiencyAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1)
{
  // Custom constructor
  fHistograms = new TrackPairEfficiencyHistograms(fCard);
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fEventReader = NULL;
  
  // Configurure the analyzer from input card
  ReadConfigurationFromCard();
  
  // pT weight function for Pythia to match 2017 MC and data pT spectra. Derived from all jets above 120 GeV
  fPtWeightFunction = new TF1("fPtWeightFunction","pol3",0,500);
  fPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09); // From JECv6
  
  // Find the correct folder for track correction tables based on data type
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC){
    
    // Weight function for 2017 MC
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);  // Weight function for 2017 MC
    fVzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
    
    // Track correction for 2017 pp data
    fTrackEfficiencyCorrector2018 = new TrkEff2017pp(false, "trackCorrectionTables/pp2017/");
    
  } else if (fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
    
    // The vz weight function is rederived from the miniAOD dataset.
    // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: f95771aa3242a7a9ed385c1ff364500481088eec
    // Input files: eecAnalysis_akFlowJets_wtaAxis_cutBadPhi_miniAODtesting_processed_2023-01-30.root
    //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_mAOD_4pC_wtaAxis_jetTrig_cutBadPhi_processed_2023-02-10.root
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);
    fVzWeightFunction->SetParameters(1.0082, -0.0190011, 0.000779051, -2.15118e-05, -6.70894e-06, 1.47181e-07, 6.65274e-09);
    
    // The centrality weight function is rederived for the miniAOD dataset.
    // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: f95771aa3242a7a9ed385c1ff364500481088eec
    // Input files: eecAnalysis_akFlowJets_wtaAxis_cutBadPhi_miniAODtesting_processed_2023-01-30.root
    //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_mAOD_4pC_wtaAxis_jetTrig_cutBadPhi_processed_2023-02-10.root
    fCentralityWeightFunctionCentral = new TF1("fCentralWeight","pol6",0,30);
    fCentralityWeightFunctionCentral->SetParameters(4.44918,-0.0544424, -0.0248668,0.00254486,-0.000117819,2.65985e-06,-2.35606e-08);
    fCentralityWeightFunctionPeripheral = new TF1("fPeripheralWeight","pol6",30,90);
    fCentralityWeightFunctionPeripheral->SetParameters(3.41938,-0.0643178, -0.00186948,7.67356e-05,-1.06981e-06,7.04102e-09,-1.84554e-11);
    
    // Track correction for 2018 PbPb data
    fTrackEfficiencyCorrector2018 = new TrkEff2018PbPb("general", false, "trackCorrectionTables/PbPb2018/");
    
  } else {
    fVzWeightFunction = NULL;
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
  }
  
}

/*
 * Copy constructor
 */
TrackPairEfficiencyAnalyzer::TrackPairEfficiencyAnalyzer(const TrackPairEfficiencyAnalyzer& in) :
  fEventReader(in.fEventReader),
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms),
  fVzWeightFunction(in.fVzWeightFunction),
  fCentralityWeightFunctionCentral(in.fCentralityWeightFunctionCentral),
  fCentralityWeightFunctionPeripheral(in.fCentralityWeightFunctionPeripheral),
  fPtWeightFunction(in.fPtWeightFunction),
  fDataType(in.fDataType),
  fJetType(in.fJetType),
  fUseTrigger(in.fUseTrigger),
  fDebugLevel(in.fDebugLevel),
  fVzWeight(in.fVzWeight),
  fCentralityWeight(in.fCentralityWeight),
  fPtHatWeight(in.fPtHatWeight),
  fTotalEventWeight(in.fTotalEventWeight),
  fJetAxis(in.fJetAxis),
  fVzCut(in.fVzCut),
  fMinimumPtHat(in.fMinimumPtHat),
  fMaximumPtHat(in.fMaximumPtHat),
  fJetEtaCut(in.fJetEtaCut),
  fJetMinimumPtCut(in.fJetMinimumPtCut),
  fJetMaximumPtCut(in.fJetMaximumPtCut),
  fCutBadPhiRegion(in.fCutBadPhiRegion),
  fMinimumMaxTrackPtFraction(in.fMinimumMaxTrackPtFraction),
  fMaximumMaxTrackPtFraction(in.fMaximumMaxTrackPtFraction),
  fTrackEtaCut(in.fTrackEtaCut),
  fTrackMinPtCut(in.fTrackMinPtCut),
  fTrackMaxPtCut(in.fTrackMaxPtCut),
  fMaxTrackPtRelativeError(in.fMaxTrackPtRelativeError),
  fMaxTrackDistanceToVertex(in.fMaxTrackDistanceToVertex),
  fCalorimeterSignalLimitPt(in.fCalorimeterSignalLimitPt),
  fHighPtEtFraction(in.fHighPtEtFraction),
  fChi2QualityCut(in.fChi2QualityCut),
  fMinimumTrackHits(in.fMinimumTrackHits),
  fSubeventCut(in.fSubeventCut)
{
  // Copy constructor
  
}

/*
 * Assingment operator
 */
TrackPairEfficiencyAnalyzer& TrackPairEfficiencyAnalyzer::operator=(const TrackPairEfficiencyAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fEventReader = in.fEventReader;
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  fVzWeightFunction = in.fVzWeightFunction;
  fCentralityWeightFunctionCentral = in.fCentralityWeightFunctionCentral;
  fCentralityWeightFunctionPeripheral = in.fCentralityWeightFunctionPeripheral;
  fPtWeightFunction = in.fPtWeightFunction;
  fDataType = in.fDataType;
  fJetType = in.fJetType;
  fUseTrigger = in.fUseTrigger;
  fDebugLevel = in.fDebugLevel;
  fVzWeight = in.fVzWeight;
  fCentralityWeight = in.fCentralityWeight;
  fPtHatWeight = in.fPtHatWeight;
  fTotalEventWeight = in.fTotalEventWeight;
  fJetAxis = in.fJetAxis;
  fVzCut = in.fVzCut;
  fMinimumPtHat = in.fMinimumPtHat;
  fMaximumPtHat = in.fMaximumPtHat;
  fJetEtaCut = in.fJetEtaCut;
  fJetMinimumPtCut = in.fJetMinimumPtCut;
  fJetMaximumPtCut = in.fJetMaximumPtCut;
  fCutBadPhiRegion = in.fCutBadPhiRegion;
  fMinimumMaxTrackPtFraction = in.fMinimumMaxTrackPtFraction;
  fMaximumMaxTrackPtFraction = in.fMaximumMaxTrackPtFraction;
  fTrackEtaCut = in.fTrackEtaCut;
  fTrackMinPtCut = in.fTrackMinPtCut;
  fTrackMaxPtCut = in.fTrackMaxPtCut;
  fMaxTrackPtRelativeError = in.fMaxTrackPtRelativeError;
  fMaxTrackDistanceToVertex = in.fMaxTrackDistanceToVertex;
  fCalorimeterSignalLimitPt = in.fCalorimeterSignalLimitPt;
  fHighPtEtFraction = in.fHighPtEtFraction;
  fChi2QualityCut = in.fChi2QualityCut;
  fMinimumTrackHits = in.fMinimumTrackHits;
  fSubeventCut = in.fSubeventCut;
  
  return *this;
}

/*
 * Destructor
 */
TrackPairEfficiencyAnalyzer::~TrackPairEfficiencyAnalyzer(){
  // destructor
  delete fHistograms;
  if(fVzWeightFunction) delete fVzWeightFunction;
  if(fCentralityWeightFunctionCentral) delete fCentralityWeightFunctionCentral;
  if(fCentralityWeightFunctionPeripheral) delete fCentralityWeightFunctionPeripheral;
  if(fPtWeightFunction) delete fPtWeightFunction;
  if(fTrackEfficiencyCorrector2018) delete fTrackEfficiencyCorrector2018;
  if(fEventReader) delete fEventReader;
}

/*
 * Read all the configuration from the input card
 */
void TrackPairEfficiencyAnalyzer::ReadConfigurationFromCard(){
  
  //****************************************
  //     Analyzed data type and trigger
  //****************************************
  fDataType = fCard->Get("DataType");
  
  //****************************************
  //         Event selection cuts
  //****************************************
  
  fVzCut = fCard->Get("ZVertexCut");             // Event cut vor the z-position of the primary vertex
  fMinimumPtHat = fCard->Get("LowPtHatCut");     // Minimum accepted pT hat value
  fMaximumPtHat = fCard->Get("HighPtHatCut");    // Maximum accepted pT hat value
  fUseTrigger = (fCard->Get("UseTrigger") == 1); // Flag telling if jet trigger is used in the analysis
  
  //****************************************
  //          Jet selection cuts
  //****************************************
  
  fJetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  fJetMinimumPtCut = fCard->Get("MinJetPtCut");   // Minimum pT cut for jets
  fJetMaximumPtCut = fCard->Get("MaxJetPtCut");   // Maximum pT accepted for jets (and tracks)
  fMinimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  fMaximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  fCutBadPhiRegion = (fCard->Get("CutBadPhi") == 1);   // Flag for cutting the phi region with bad tracking efficiency from the analysis
  

  //****************************************
  //            Jet selection
  //****************************************
  fJetType = fCard->Get("JetType");              // Select the type of analyzed jets (Calo, CSPF, PuPF, FlowPF)
  fJetAxis = fCard->Get("JetAxis");              // Select between escheme and WTA axes

  //****************************************
  //        Track selection cuts
  //****************************************
  
  fTrackEtaCut = fCard->Get("TrackEtaCut");     // Eta cut around midrapidity
  fTrackMinPtCut = fCard->Get("MinTrackPtCut"); // Minimum track pT cut
  fTrackMaxPtCut = fCard->Get("MaxTrackPtCut"); // Maximum track pT cut
  fMaxTrackPtRelativeError = fCard->Get("MaxTrackPtRelativeError");   // Maximum relative error for pT
  fMaxTrackDistanceToVertex = fCard->Get("VertexMaxDistance");        // Maximum distance to primary vetrex
  fCalorimeterSignalLimitPt = fCard->Get("CalorimeterSignalLimitPt"); // Require signal in calorimeters for track above this pT
  fHighPtEtFraction = fCard->Get("HighPtEtFraction"); // For high pT tracks, minimum required Et as a fraction of track pT
  fChi2QualityCut = fCard->Get("Chi2QualityCut");     // Quality cut for track reconstruction
  fMinimumTrackHits = fCard->Get("MinimumTrackHits"); // Quality cut for track hits
  fSubeventCut = fCard->Get("SubeventCut");           // Cut on subevent index in MC
  
  //************************************************
  //              Debug messages
  //************************************************
  fDebugLevel = fCard->Get("DebugLevel");
}

/*
 * Main analysis loop
 */
void TrackPairEfficiencyAnalyzer::RunAnalysis(){
  
  //************************************************
  //  Define variables needed in the analysis loop
  //************************************************
  
  // Input files and forest readers for analysis
  TFile *inputFile;
  
  // Event variables
  Int_t nEvents = 0;                // Number of events
  Double_t vz = 0;                  // Vertex z-position
  Double_t centrality = 0;          // Event centrality
  Int_t hiBin = 0;                  // CMS hiBin (centrality * 2)
  Double_t ptHat = 0;               // pT hat for MC events
  
  // Variables for tracks
  Double_t fillerTrack[4];                   // Track histogram filler
  Double_t fillerTrackPair[6];               // Track pair histogram filler
  Int_t nTracks;                             // Number of tracks in an event
  Double_t averagePairEta = 0;               // Average eta of the track pair
  Double_t averagePairPhi = 0;               // Average phi of the track pair
  Double_t pairDeltaR = 0;                   // DeltaR between the two tracks in a pair
  
  // Vectors in attempt to make the track pairings faster
  vector<double> selectedTrackPt;         // Track pT for the tracks passing the tracking cuts
  vector<double> selectedTrackEta;        // Track eta for the tracks passing the tracking cuts
  vector<double> selectedTrackPhi;        // Track phi for the tracks passing the tracking cuts
  vector<double> selectedTrackEfficiency; // Track efficiency for the tracks passing the tracking cuts
  
  // Variables for jets
  Int_t nJets = 0;                  // Number of jets in an event
  Double_t jetPt = 0;               // pT of the i:th jet in the event
  Double_t jetPhi = 0;              // phi of the i:th jet in the event
  Double_t jetEta = 0;              // eta of the i:th jet in the event
  Double_t jetPtWeight = 1;         // Weighting for jet pT
  
  // File name helper variables
  TString currentFile;
  
  // Fillers for THnSparses
  const Int_t nFillJet = 5;
  Double_t fillerJet[nFillJet];

  
  //************************************************
  //      Define forest reader for data files
  //************************************************
  
  fEventReader = new ForestReader(fDataType, fJetType, fJetAxis, fUseTrigger);
  
  
  //************************************************
  //       Main analysis loop over all files
  //************************************************
  
  // Loop over files
  Int_t nFiles = fFileNames.size();
  for(Int_t iFile = 0; iFile < nFiles; iFile++) {
    
    //************************************************
    //              Find and open files
    //************************************************
    
    // Find the filename and open the input file
    currentFile = fFileNames.at(iFile);
    inputFile = TFile::Open(currentFile);
    
    // Check that the file exists
    if(!inputFile){
      cout << "Error! Could not find the file: " << currentFile.Data() << endl;
      assert(0);
    }

    // Check that the file is open
    if(!inputFile->IsOpen()){
      cout << "Error! Could not open the file: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Error! The following file is a zombie: " << currentFile.Data() << endl;
      assert(0);
    }
    

    // Print the used files
    if(fDebugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;

    
    //************************************************
    //            Read forest from file
    //************************************************
    
    // If file is good, read the forest from the file
    fEventReader->ReadForestFromFile(inputFile);  // There might be a memory leak in handling the forest...
    nEvents = fEventReader->GetNEvents();

    //************************************************
    //         Main event loop for each file
    //************************************************
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){ // nEvents
      
      //************************************************
      //         Read basic event information
      //************************************************
      
      // Print to console how the analysis is progressing
      if(fDebugLevel > 1 && iEvent % 1000 == 0) cout << "Analyzing event " << iEvent << endl;
      
      // Read the event to memory
      fEventReader->GetEvent(iEvent);

      // Get vz, centrality and pT hat information
      vz = fEventReader->GetVz();
      centrality = fEventReader->GetCentrality();
      hiBin = fEventReader->GetHiBin();
      ptHat = fEventReader->GetPtHat();
      
      // We need to apply pT hat cuts before getting pT hat weight. There might be rare events above the upper
      // limit from which the weights are calculated, which could cause the code to crash.
      if(ptHat < fMinimumPtHat || ptHat >= fMaximumPtHat) continue;
      
      // Get the weighting for the event
      fVzWeight = GetVzWeight(vz);
      fCentralityWeight = GetCentralityWeight(hiBin);
      fPtHatWeight = fEventReader->GetEventWeight();
      fTotalEventWeight = fVzWeight*fCentralityWeight*fPtHatWeight;
      
      // Fill event counter histogram
      fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kAll);          // All the events looped over
      
      //  ============================================
      //  ===== Apply all the event quality cuts =====
      //  ============================================
      
      if(!PassEventCuts(fEventReader)) continue;
      
      // Fill the event information histograms for the events that pass the event cuts
      fHistograms->fhVertexZ->Fill(vz);                            // z vertex distribution from all events
      fHistograms->fhVertexZWeighted->Fill(vz,fVzWeight);          // z-vertex distribution weighted with the weight function
      fHistograms->fhCentrality->Fill(centrality);                 // Centrality filled from all events
      fHistograms->fhCentralityWeighted->Fill(centrality,fCentralityWeight); // Centrality weighted with the centrality weighting function
      fHistograms->fhPtHat->Fill(ptHat);                           // pT hat histogram
      fHistograms->fhPtHatWeighted->Fill(ptHat,fPtHatWeight);      // pT het histogram weighted with corresponding cross section and event number
      
      // ======================================
      // ===== Event quality cuts applied =====
      // ======================================
      
      //***********************************************************************
      //             Collect basic track distribution hisotgrams
      //***********************************************************************
      
      // Clear the track vectors
      selectedTrackPt.clear();
      selectedTrackEta.clear();
      selectedTrackPhi.clear();
      selectedTrackEfficiency.clear();
      
      // Loop over all track in the event
      nTracks = fEventReader->GetNTracks();
      for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
        
        // Check that all the track cuts are passed
        if(!PassTrackCuts(fEventReader,iTrack,fHistograms->fhTrackCuts,false)) continue;
        
        // Get the track information and add it to vectors
        selectedTrackPt.push_back(fEventReader->GetTrackPt(iTrack));
        selectedTrackEta.push_back(fEventReader->GetTrackEta(iTrack));
        selectedTrackPhi.push_back(fEventReader->GetTrackPhi(iTrack));
        selectedTrackEfficiency.push_back(GetTrackEfficiencyCorrection(iTrack));
        
        // Fill track histograms
        fillerTrack[0] = selectedTrackPt.back();      // Axis 0: Track pT
        fillerTrack[1] = selectedTrackPhi.back();     // Axis 1: Track phi
        fillerTrack[2] = selectedTrackEta.back();     // Axis 2: Track eta
        fillerTrack[3] = centrality;                  // Axis 3: Centrality
        fHistograms->fhTrack->Fill(fillerTrack,selectedTrackEfficiency.back()*fTotalEventWeight);  // Fill the track histogram
        fHistograms->fhTrackUncorrected->Fill(fillerTrack,fTotalEventWeight);                      // Fill the uncorrected track histogram
        
      } // Track loop
      
      // Once we have looped over all the tracks, only loop over tracks that pass the cuts to construct all possible track pairings
      for(Int_t iTrack = 0; iTrack < selectedTrackPt.size(); iTrack++){
        for(Int_t jTrack = iTrack+1; jTrack < selectedTrackPt.size(); jTrack++){
          
          // Calculate the distance of the two tracks from each other
          pairDeltaR = GetDeltaR(selectedTrackEta.at(iTrack), selectedTrackPhi.at(iTrack), selectedTrackEta.at(jTrack), selectedTrackPhi.at(jTrack));
          
          // Fill the track pair histograms for tracks relatively close to each other
          if(pairDeltaR < 0.8){
            
            // Calculate the average pair eta and phi positions
            averagePairEta = (selectedTrackEta.at(iTrack)+selectedTrackEta.at(jTrack))/2.0;
            averagePairPhi = (selectedTrackPhi.at(iTrack)+selectedTrackPhi.at(jTrack))/2.0;
            
            fillerTrackPair[0] = pairDeltaR;                   // Axis 0: DeltaR between the two tracks
            fillerTrackPair[1] = selectedTrackPt.at(iTrack);   // Axis 1: First track pT
            fillerTrackPair[2] = selectedTrackPt.at(jTrack);   // Axis 2: Second track pT
            fillerTrackPair[3] = averagePairPhi;               // Axis 3: Average pair phi
            fillerTrackPair[4] = averagePairEta;               // Axis 4: Average pair eta
            fillerTrackPair[5] = centrality;                   // Axis 5: Centrality
            fHistograms->fhTrackPairs->Fill(fillerTrackPair, selectedTrackEfficiency.at(iTrack) * selectedTrackEfficiency.at(jTrack) * fTotalEventWeight);  // Fill the track pair histogram
          }
          
        } // Inner track loop
      } // Outer track loop
      
      // Do the same for generator level tracks in case for running with Monte Carlo
      if(fDataType == ForestReader::kPpMC || fDataType == ForestReader::kPbPbMC){
       
        // Clear the track vectors
        selectedTrackPt.clear();
        selectedTrackEta.clear();
        selectedTrackPhi.clear();
        selectedTrackEfficiency.clear();
        
        nTracks = fEventReader->GetNGenParticles();
        for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
          
          // Check that all the particle selections are passed
          if(!PassGenParticleSelection(fEventReader,iTrack,fHistograms->fhGenParticleSelections,false)) continue;
          
          // Get the efficiency correction
          selectedTrackPt.push_back(fEventReader->GetGenParticlePt(iTrack));
          selectedTrackEta.push_back(fEventReader->GetGenParticleEta(iTrack));
          selectedTrackPhi.push_back(fEventReader->GetGenParticlePhi(iTrack));
          
          // Fill track histograms
          fillerTrack[0] = selectedTrackPt.back();      // Axis 0: Generator level particle pT
          fillerTrack[1] = selectedTrackPhi.back();     // Axis 1: Generator level particle phi
          fillerTrack[2] = selectedTrackEta.back();     // Axis 2: Generator level particle eta
          fillerTrack[3] = centrality;                  // Axis 3: Centrality
          fHistograms->fhGenParticle->Fill(fillerTrack,fTotalEventWeight);  // Fill the generator level particle histogram
        }
        
        // Once we have looped over all the tracks, only loop over tracks that pass the cuts to construct all possible track pairings
        for(Int_t iTrack = 0; iTrack < selectedTrackPt.size(); iTrack++){
          for(Int_t jTrack = iTrack+1; jTrack < selectedTrackPt.size(); jTrack++){

            // Calculate the distance of the two tracks from each other
            pairDeltaR = GetDeltaR(selectedTrackEta.at(iTrack), selectedTrackPhi.at(iTrack), selectedTrackEta.at(jTrack), selectedTrackPhi.at(jTrack));

            // Fill the track pair histograms for tracks relatively close to each other
            if(pairDeltaR < 0.8){

              // Calculate the average pair eta and phi positions
              averagePairEta = (selectedTrackEta.at(iTrack)+selectedTrackEta.at(jTrack))/2.0;
              averagePairPhi = (selectedTrackPhi.at(iTrack)+selectedTrackPhi.at(jTrack))/2.0;

              fillerTrackPair[0] = pairDeltaR;                   // Axis 0: DeltaR between the two tracks
              fillerTrackPair[1] = selectedTrackPt.at(iTrack);   // Axis 1: First track pT
              fillerTrackPair[2] = selectedTrackPt.at(jTrack);   // Axis 2: Second track pT
              fillerTrackPair[3] = averagePairPhi;               // Axis 3: Average pair phi
              fillerTrackPair[4] = averagePairEta;               // Axis 4: Average pair eta
              fillerTrackPair[5] = centrality;                   // Axis 5: Centrality
              fHistograms->fhGenParticlePairs->Fill(fillerTrackPair,fTotalEventWeight);  // Fill the track pair histogram
            }

          } // Inner track loop
        } // Outer track loop
        
      }
      
      
      //***********************************************************************
      //    Loop over all jets and fill histograms for different triggers
      //***********************************************************************
      
      // Jet loop
      nJets = fEventReader->GetNJets();
      for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++) {
        
        jetPt = fEventReader->GetJetPt(jetIndex);
        jetPhi = fEventReader->GetJetPhi(jetIndex);
        jetEta = fEventReader->GetJetEta(jetIndex);
        
        //  ========================================
        //  ======== Apply jet quality cuts ========
        //  ========================================
        
        if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
        if(fCutBadPhiRegion && (jetPhi > -0.1 && jetPhi < 1.2)) continue; // Cut the area of large inefficiency in tracker
        
        if(fMinimumMaxTrackPtFraction >= fEventReader->GetJetMaxTrackPt(jetIndex)/fEventReader->GetJetRawPt(jetIndex)) {
          continue; // Cut for jets with only very low pT particles
        }
        if(fMaximumMaxTrackPtFraction <= fEventReader->GetJetMaxTrackPt(jetIndex)/fEventReader->GetJetRawPt(jetIndex)) {
          continue; // Cut for jets where all the pT is taken by one track
        }
        
        
        //  ========================================
        //  ======= Jet quality cuts applied =======
        //  ========================================
        
        // After the jet pT can been corrected, apply analysis jet pT cuts
        if(jetPt < fJetMinimumPtCut) continue;
        if(jetPt > fJetMaximumPtCut) continue;
        
        //************************************************
        //         Fill histograms for all jets
        //************************************************
        
        // Find the pT weight for the jet
        jetPtWeight = GetJetPtWeight(jetPt);
        
        // Fill the axes in correct order
        fillerJet[0] = jetPt;          // Axis 0 = jet pT
        fillerJet[1] = jetPhi;         // Axis 1 = jet phi
        fillerJet[2] = jetEta;         // Axis 2 = jet eta
        fillerJet[3] = centrality;     // Axis 3 = centrality
        fillerJet[4] = TrackPairEfficiencyHistograms::kReconstructed;  // Axis 4 = Reconstruction flag
        
        fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
        
      } // End of jet loop
      
      // For MC, do another jet loop using generator level jets
      if(fDataType == ForestReader::kPpMC || fDataType == ForestReader::kPbPbMC){
        
        // Generator level jet loop
        nJets = fEventReader->GetNGeneratorJets();
        for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++) {
          
          jetPt = fEventReader->GetGeneratorJetPt(jetIndex);
          jetPhi = fEventReader->GetGeneratorJetPhi(jetIndex);
          jetEta = fEventReader->GetGeneratorJetEta(jetIndex);
          
          //  ==========================================
          //  ======== Apply jet kinematic cuts ========
          //  ==========================================
          
          if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
          if(jetPt < fJetMinimumPtCut) continue;
          if(jetPt > fJetMaximumPtCut) continue;
          
          //************************************************
          //     Fill histograms for generator level jets
          //************************************************

          // Find the pT weight for the jet
          jetPtWeight = GetJetPtWeight(jetPt);
          
          // Fill the axes in correct order
          fillerJet[0] = jetPt;          // Axis 0 = generator level jet pT
          fillerJet[1] = jetPhi;         // Axis 1 = generator level jet phi
          fillerJet[2] = jetEta;         // Axis 2 = generator level jet eta
          fillerJet[3] = centrality;     // Axis 3 = centrality
          fillerJet[4] = TrackPairEfficiencyHistograms::kGeneratorLevel;   // Axis 4 = Generator level flag
          
          fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
          
        } // End of jet loop
        
      } // MC if
      
      
    } // Event loop
    
    //************************************************
    //      Cleanup at the end of the file loop
    //************************************************
    
    // Close the input files after the event has been read
    inputFile->Close();
    
  } // File loop
  
}


/*
 * Get the proper vz weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t vz = Vertex z position for the event
 *
 *   return: Multiplicative correction factor for vz
 */
Double_t TrackPairEfficiencyAnalyzer::GetVzWeight(const Double_t vz) const{
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) return 1;  // No correction for real data
  if(fDataType == ForestReader::kPbPbMC || fDataType == ForestReader::kPpMC) return fVzWeightFunction->Eval(vz); // Weight for 2018 MC
  return -1; // Return crazy value for unknown data types, so user will not miss it
}

/*
 * Get the proper centrality weighting depending on analyzed system
 *
 *  Arguments:
 *   const Int_t hiBin = CMS hiBin
 *
 *   return: Multiplicative correction factor for the given CMS hiBin
 */
Double_t TrackPairEfficiencyAnalyzer::GetCentralityWeight(const Int_t hiBin) const{
  if(fDataType != ForestReader::kPbPbMC) return 1;
  
  // No weighting for the most peripheral centrality bins. Different weight function for central and peripheral.
  if(hiBin < 60) return fCentralityWeightFunctionCentral->Eval(hiBin/2.0);
  return (hiBin < 194) ? fCentralityWeightFunctionPeripheral->Eval(hiBin/2.0) : 1;
}

/*
 * Get the proper jet pT weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t jetPt = Jet pT for the weighted jet
 *
 *   return: Multiplicative correction factor for the jet pT
 */
Double_t TrackPairEfficiencyAnalyzer::GetJetPtWeight(const Double_t jetPt) const{
  if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPp) return 1.0;  // No weight for data
  
  return fPtWeightFunction->Eval(jetPt);
}

/*
 * Check if the event passes all the event cuts
 *
 *  Arguments:
 *   ForestReader *eventReader = ForestReader containing the event information checked for event cuts
 *
 *   return = True if all event cuts are passed, false otherwise
 */
Bool_t TrackPairEfficiencyAnalyzer::PassEventCuts(ForestReader *eventReader){

  // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction. Only applied for data.
  if(eventReader->GetPrimaryVertexFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kPrimaryVertex);
  
  // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV. Only applied for PbPb data.
  if(eventReader->GetHfCoincidenceFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kHfCoincidence);
  
  // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible. Only applied for PbPb data.
  if(eventReader->GetClusterCompatibilityFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kClusterCompatibility);
  
  // Cut for beam scraping. Only applied for pp data.
  if(eventReader->GetBeamScrapingFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kBeamScraping);
  
  // Jet trigger requirement.
  if(eventReader->GetJetFilterBit() == 0) return false;
  fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kCaloJet);
  
  // Cut for vertex z-position
  if(TMath::Abs(eventReader->GetVz()) > fVzCut) return false;
  fHistograms->fhEvents->Fill(TrackPairEfficiencyHistograms::kVzCut);
  
  return true;
  
}

/*
 * Check if a track passes the generator level particle selection criteria
 *
 *  Arguments:
 *   ForestReader *trackReader = ForestReader from which the tracks are read
 *   const Int_t iTrack = Index of the checked track in reader
 *   TH1F *trackCutHistogram = Histogram to which the track cut performance is filled
 *   const Bool_t bypassFill = Pass filling the track cut histograms
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t TrackPairEfficiencyAnalyzer::PassGenParticleSelection(ForestReader *trackReader, const Int_t iTrack, TH1F *trackCutHistogram, const Bool_t bypassFill){
 
  // Only fill the track cut histograms for same event data
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kAllGenParticles);
  
  // Cuts specific to generator level MC tracks
  if(trackReader->GetGenParticleCharge(iTrack) == 0) return false;  // Require that the track is charged
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kMcCharge);
  
  if(!PassSubeventCut(trackReader->GetGenParticleSubevent(iTrack))) return false;  // Require desired subevent
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kMcSube);
  
  Double_t trackPt = trackReader->GetGenParticlePt(iTrack);
  
  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;         // Minimum track pT cut
  if(trackPt >= fTrackMaxPtCut) return false;         // Maximum track pT cut
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kGenParticlePtCut);
  
  Double_t trackEta = trackReader->GetGenParticleEta(iTrack);
  
  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kGenParticleEtaCut);
  
  // If passed all checks, return true
  return true;
  
}


/*
 * Check if a track passes all the track cuts
 *
 *  Arguments:
 *   ForestReader *trackReader = ForestReader from which the tracks are read
 *   const Int_t iTrack = Index of the checked track in reader
 *   TH1F *trackCutHistogram = Histogram to which the track cut performance is filled
 *   const Bool_t bypassFill = Pass filling the track cut histograms
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t TrackPairEfficiencyAnalyzer::PassTrackCuts(ForestReader *trackReader, const Int_t iTrack, TH1F *trackCutHistogram, const Bool_t bypassFill){
  
  // Only fill the track cut histograms for same event data
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kAllTracks);
  
  Double_t trackPt = trackReader->GetTrackPt(iTrack);
  Double_t trackEta = trackReader->GetTrackEta(iTrack);
  Double_t trackEt = (trackReader->GetTrackEnergyEcal(iTrack)+trackReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
  
  //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
  
  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;                 // Minimum track pT cut
  if(trackPt >= fTrackMaxPtCut) return false;                 // Maximum track pT cut
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kPtCuts);
  
  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kEtaCut);
  
  // Cut for high purity
  if(!trackReader->GetTrackHighPurity(iTrack)) return false;     // High purity cut
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kHighPurity);
  
  // Cut for relative error for track pT
  if(trackReader->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) return false; // Cut for track pT relative error
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kPtError);
  
  // Cut for track distance from primary vertex
  if(TMath::Abs(trackReader->GetTrackVertexDistanceZ(iTrack)/trackReader->GetTrackVertexDistanceZError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in z-direction
  if(TMath::Abs(trackReader->GetTrackVertexDistanceXY(iTrack)/trackReader->GetTrackVertexDistanceXYError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in xy-direction
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kVertexDistance);
  
  // Cut for energy deposition in calorimeters for high pT tracks
  if(!(trackPt < fCalorimeterSignalLimitPt || (trackEt >= fHighPtEtFraction*trackPt))) return false;  // For high pT tracks, require signal also in calorimeters
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kCaloSignal);
  
  // Cuts for track reconstruction quality
  //if( trackReader->GetTrackChi2(iTrack) / (1.0*trackReader->GetNTrackDegreesOfFreedom(iTrack)) / (1.0*trackReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  if(trackReader->GetTrackNormalizedChi2(iTrack) / (1.0*trackReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  if(trackReader->GetNHitsTrack(iTrack) < fMinimumTrackHits) return false; // Cut for minimum number of hits per track
  if(!bypassFill) trackCutHistogram->Fill(TrackPairEfficiencyHistograms::kReconstructionQuality);
  
  // If passed all checks, return true
  return true;
}

/*
 * Get the track efficiency correction for a given track
 *
 *  Arguments:
 *   const Int_t iTrack = Index of the track for which the efficiency correction is obtained
 *
 *   return: Multiplicative track efficiency correction
 */
Double_t TrackPairEfficiencyAnalyzer::GetTrackEfficiencyCorrection(const Int_t iTrack){
  
  // Get track information
  Float_t trackPt = fEventReader->GetTrackPt(iTrack);    // Track pT
  Float_t trackEta = fEventReader->GetTrackEta(iTrack);  // Track eta
  Int_t hiBin = fEventReader->GetHiBin();                // hiBin for 2018 track correction
  
  // Get the correction using the track and event information
  return GetTrackEfficiencyCorrection(trackPt, trackEta, hiBin);
  
}

/*
 * Get the track efficiency correction according to given information
 *
 *  Arguments:
 *   const Float_t trackPt = pT of the track
 *   const Float_t trackEta = Eta of the track
 *   const Int_t hiBin = CMS hiBin index. Basically centrality/2
 *
 *   return: Multiplicative track efficiency correction
 */
Double_t TrackPairEfficiencyAnalyzer::GetTrackEfficiencyCorrection(const Float_t trackPt, const Float_t trackEta, const Int_t hiBin){
  
  // Weight factor only for 2017 pp MC as instructed be the tracking group
  double preWeight = 1.0;
  if(fDataType == ForestReader::kPpMC) preWeight = 0.979;
  
  // For PbPb2018 and pp2017, there is an efficiency table from which the correction comes
  return preWeight * fTrackEfficiencyCorrector2018->getCorrection(trackPt, trackEta, hiBin);
  
}

/*
 * Check is a track passes the required subevent cut
 *
 *  Arguments:
 *   const Int_t subeventIndex = Subevent index for the track in consideration
 *
 *  return: true if subevent cut is passes, false if not
 */
Bool_t TrackPairEfficiencyAnalyzer::PassSubeventCut(const Int_t subeventIndex) const{
  if(fSubeventCut == kSubeventAny) return true;
  if((fSubeventCut == kSubeventZero) && (subeventIndex == 0)) return true;
  if((fSubeventCut == kSubeventNonZero) && (subeventIndex > 0)) return true;
  return false;
}

/*
 * Getter for trigger histograms
 */
TrackPairEfficiencyHistograms* TrackPairEfficiencyAnalyzer::GetHistograms() const{
  return fHistograms;
}

/*
 * Get deltaR between two objects
 *
 *  Arguments:
 *   const Double_t eta1 = Eta of the first object
 *   const Double_t phi1 = Phi of the first object
 *   const Double_t eta2 = Eta of the second object
 *   const Double_t phi2 = Phi of the second object
 *
 *  return: DeltaR between the two objects
 */
Double_t TrackPairEfficiencyAnalyzer::GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const{

  Double_t deltaEta = eta1 - eta2;
  Double_t deltaPhi = phi1 - phi2;
  
  // Transform deltaPhi to interval [-pi,pi]
  while(deltaPhi > TMath::Pi()){deltaPhi += -2*TMath::Pi();}
  while(deltaPhi < -TMath::Pi()){deltaPhi += 2*TMath::Pi();}
  
  // Return the distance between the objects
  return TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  
}
