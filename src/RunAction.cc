//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//      ----------------------------------------------------------------
//                      K600 Spectrometer (iThemba Labs)
//      ----------------------------------------------------------------
//
//      Github repository: https://www.github.com/KevinCWLi/K600
//
//      Main Author:    K.C.W. Li
//
//      email: likevincw@gmail.com
//

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
: G4UserRunAction()
{
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);
    
    // Create analysis manager
    // The choice of analysis technology is done via selectin of a namespace
    // in Analysis.hh
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
    
    // Create directories
    //analysisManager->SetHistoDirectoryName("histograms");
    //analysisManager->SetNtupleDirectoryName("ntuple");
    analysisManager->SetVerboseLevel(1);
    //analysisManager->SetFirstHistoId(1);
    
    // Book histograms, ntuple
    //
    
    // Creating 1D histograms
    /*
     ////    CAKE DETECTORS, Histograms 1-6
     analysisManager->CreateH1("CvsE_CAKE1","CAKE 1 - Counts versus Energy", 10000, 0., 10000.);
     analysisManager->CreateH1("CvsE_CAKE2","CAKE 2 - Counts versus Energy", 10000, 0., 10000.);
     analysisManager->CreateH1("CvsE_CAKE3","CAKE 3 - Counts versus Energy", 10000, 0., 10000.);
     analysisManager->CreateH1("CvsE_CAKE4","CAKE 4 - Counts versus Energy", 10000, 0., 10000.);
     analysisManager->CreateH1("CvsE_CAKE5","CAKE 5 - Counts versus Energy", 10000, 0., 10000.);
     analysisManager->CreateH1("CvsE_CAKE_EA","Entire CAKE Array - Counts versus Energy", 10000, 0., 10000.);
     
     ////    PADDLE DETECTORS, Histograms 7-9
     analysisManager->CreateH1("CvsE_PADDLE1","PADDLE 1 - Counts versus Energy", 1000, 0., 100.);
     analysisManager->CreateH1("CvsE_PADDLE2","PADDLE 2 - Counts versus Energy", 1000, 0., 100.);
     analysisManager->CreateH1("CvsE_PADDLE3","PADDLE 3 - Counts versus Energy", 1000, 0., 100.);
     
     ////    CLOVER DETECTORS, Histograms 10-19
     analysisManager->CreateH1("CvsE_CLOVER1","CLOVER 1 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER2","CLOVER 2 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER3","CLOVER 3 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER4","CLOVER 4 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER5","CLOVER 5 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER6","CLOVER 6 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER7","CLOVER 7 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER8","CLOVER 8 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER9","CLOVER 9 - Counts versus Energy", 3000, 0., 3000.);
     analysisManager->CreateH1("CvsE_CLOVER_EA","Entire CLOVER Array - Counts versus Energy", 3000, 0., 3000.);
     
     
     // Creating 2D histograms
     
     ////    PADDLE DETECTORS, Histograms 1-6
     analysisManager->CreateH2("PvsE_PADDLE1","PADDLE 1 - Position vs. Energy", 1400, -700., 700., 3, -3*(102/2), 3*(102/2));
     analysisManager->CreateH2("PvsE_PADDLE2","PADDLE 2 - Position vs. Energy", 1400, -700., 700., 3, -3*(102/2), 3*(102/2));
     analysisManager->CreateH2("PvsE_PADDLE3","PADDLE 3 - Position vs. Energy", 1400, -700., 700., 3, -3*(102/2), 3*(102/2));
     analysisManager->CreateH2("EvsTOF_PADDLE1","PADDLE 1 - Energy vs. Time of Flight", 300, 0., 500.*0.3, 40, 0., 40.);
     analysisManager->CreateH2("EvsTOF_PADDLE2","PADDLE 2 - Energy vs. Time of Flight", 300, 0., 500.*0.3, 40, 0., 40.);
     analysisManager->CreateH2("EvsTOF_PADDLE3","PADDLE 3 - Energy vs. Time of Flight", 300, 0., 500.*0.3, 40, 0., 40.);
     */
    
    
    ////////////////////////////////////////////////////
    //                  DataTreeSim
    ////////////////////////////////////////////////////
    
    // Creating ntuple
    analysisManager->CreateNtuple("DataTreeSim", "K600 Spectrometer - Coincident Events");
    
    analysisManager->CreateNtupleDColumn(0, "InitialParticleKineticEnergy");
    
    //--------------------------------
    //      CLOVER Detectors
    analysisManager->CreateNtupleIColumn(0, "CLOVER_EventFold");
    analysisManager->CreateNtupleIColumn(0, "CLOVER_NCrystalsTriggered", CLOVER_nCrystalsTriggered);
    analysisManager->CreateNtupleIColumn(0, "CLOVER_Number", CLOVER_iD);
    analysisManager->CreateNtupleDColumn(0, "CLOVER_Energy", CLOVER_energy);
    analysisManager->CreateNtupleDColumn(0, "CLOVER_DetectorTheta", CLOVER_detectorTheta);
    analysisManager->CreateNtupleDColumn(0, "CLOVER_DetectorPhi", CLOVER_detectorPhi);

    //--------------------------------
    //      LaBr3Ce Detectors
    analysisManager->CreateNtupleIColumn(0, "LaBr3Ce_EventFold");
    analysisManager->CreateNtupleIColumn(0, "LaBr3Ce_Number", laBr3Ce_iD);
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_Energy", laBr3Ce_energy);
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_DetectorTheta", laBr3Ce_detectorTheta);
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_DetectorPhi", laBr3Ce_detectorPhi);
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_Theta", laBr3Ce_theta);
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_Phi", laBr3Ce_phi);
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_xPos", laBr3Ce_xPos); // cm (relative to the target/origin)
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_yPos", laBr3Ce_yPos); // cm (relative to the target/origin)
    analysisManager->CreateNtupleDColumn(0, "LaBr3Ce_zPos", laBr3Ce_zPos); // cm (relative to the target/origin)
    
    
    analysisManager->FinishNtuple(0);
    
    
    ////////////////////////////////////////////////////////////
    //                  GeometryAnalysisTree
    ////////////////////////////////////////////////////////////
    
    // Creating ntuple
    
    analysisManager->CreateNtuple("GeometryAnalysisTree", "K600 Spectrometer - GeometryAnalysis");
    
    ////    CAKE Detectors
    analysisManager->CreateNtupleIColumn(1, "CAKENo");
    analysisManager->CreateNtupleIColumn(1, "CAKE_RowNo");
    analysisManager->CreateNtupleIColumn(1, "CAKE_SectorNo");
    analysisManager->CreateNtupleDColumn(1, "Theta");
    analysisManager->CreateNtupleDColumn(1, "Phi");
    
    analysisManager->FinishNtuple(1);
    
    
    ////////////////////////////////////////////////////////////
    //                  Input Variable Tree
    ////////////////////////////////////////////////////////////
    
    // Creating ntuple
    analysisManager->CreateNtuple("InputVariableTree", "K600 Spectrometer - InputVariable");
    
    ////    Initial Particle Angular Distributions
    analysisManager->CreateNtupleDColumn(2, "ThetaDist");
    analysisManager->CreateNtupleDColumn(2, "PhiDist");
    analysisManager->CreateNtupleDColumn(2, "ThetaDist_projX");
    analysisManager->CreateNtupleDColumn(2, "ThetaDist_projY");

    analysisManager->FinishNtuple(2);
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    //inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // Open an output file
    //
    G4String fileName = "K600Output";
    analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    /*
     // print histogram statistics
     //
     if( analysisManager->GetH1(1) ){
     G4cout << "\n ----> print histograms statistic ";
     if(isMaster) {
     G4cout << "for the entire run \n" << G4endl;
     }
     else {
     G4cout << "for the local thread \n" << G4endl;
     }
     
     G4cout << " EAbs : mean = "
     << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
     << " rms = "
     << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
     
     G4cout << " EGap : mean = "
     << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy")
     << " rms = "
     << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
     
     G4cout << " LAbs : mean = "
     << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
     << " rms = "
     << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
     
     G4cout << " LGap : mean = "
     << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length")
     << " rms = "
     << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
     
     }
     */
    
    // save histograms & ntuple
    //
    analysisManager->Write();
    analysisManager->CloseFile();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
