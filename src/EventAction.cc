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

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include <fstream>
#include <string>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction, DetectorConstruction* detectorConstruction)
: G4UserEventAction(),
fRunAction(runAction),
fEnergyAbs(0.),
fEnergyGap(0.),
fTrackLAbs(0.),
fTrackLGap(0.),
/////
GainCAKE(1.0),
OffsetCAKE(0.0),
GainCLOVER(1.0),
OffsetCLOVER(0.0),
GainPADDLE(0),
OffsetPADDLE(0),
GainLEPS(1.0),
OffsetLEPS(0.0)
{    
    angles_CLOVER = detectorConstruction->GetAngles_CLOVER();
    angles_ALBA_LaBr3Ce = detectorConstruction->GetAngles_ALBA_LaBr3Ce();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
    // initialisation per event
    //fEnergyAbs = 0.;
    //fEnergyGap = 0.;
    //fTrackLAbs = 0.;
    //fTrackLGap = 0.;
    
    evtNb = evt->GetEventID();
    
    GA_LineOfSight = true;
    
    if(GA_MODE)
    {
        if(evtNb==0)
        {
            for(G4int i=0; i<5; i++)
            {
                for(G4int j=0; j<16; j++)
                {
                    for(G4int k=0; k<8; k++)
                    {
                        GA_MMM_AngDist_counter[i][j][k]= 0;
                    }
                }
            }
            
            
            for(G4int i=0; i<640; i++)
            {
                for(G4int j=0; j<4; j++)
                {
                    GA_CAKE_AA_stor[i][j] = 0;
                }
            }
        }
        
        for(G4int i=0; i<640; i++)
        {
            for(G4int j=0; j<3; j++)
            {
                GA_CAKE_AA[i][j] = 0;
            }
        }
    }
    
    
    //------------------------------------------------
    CLOVER_Number_vec.clear();
    CLOVER_Energy_vec.clear();
    CLOVER_DetectorTheta_vec.clear();
    CLOVER_DetectorPhi_vec.clear();

    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        for (G4int k=0; k<CLOVER_TotalTimeSamples; k++)
        {
            CLOVER_EDep[i][k] = 0;
            
            for(G4int j=0; j<4; j++)
            {
                CLOVER_HPGeCrystal_EDep[i][j][k] = 0;
                CLOVER_BGO_EDep[i][j][k] = 0;
                CLOVER_HPGeCrystal_EDepVETO[i][j][k] = false;
            }
        }
        
        for (G4int m=0; m<CLOVER_Shield_BGO_TotalTimeSamples; m++)
        {
            for(G4int l=0; l<16; l++)
            {
                CLOVER_BGO_EDep[i][l][m] = 0;
            }
        }
    }
    
    for(G4int i=0; i<6; i++)
    {
        for(G4int k=0; k<LEPS_TotalTimeSamples; k++)
        {
            LEPS_EDep[i][k] = 0.;

            for(G4int j=0; j<4; j++)
            {
                LEPS_HPGeCrystal_EDep[i][j][k] = 0;
            }
        }
    }
    
    //------------------------------------------------
    LaBr3Ce_Number_vec.clear();
    CLOVER_NCrystalTriggered_vec.clear();
    LaBr3Ce_Energy_vec.clear();
    LaBr3Ce_DetectorTheta_vec.clear();
    LaBr3Ce_DetectorPhi_vec.clear();
    LaBr3Ce_Theta_vec.clear();
    LaBr3Ce_Phi_vec.clear();
    LaBr3Ce_xPos_vec.clear();
    LaBr3Ce_yPos_vec.clear();
    LaBr3Ce_zPos_vec.clear();
    
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        for(G4int k=0; k<LaBr3Ce_TotalTimeSamples; k++)
        {
            LaBr3Ce_EDep[i][k] = 0.0;
        }
    }
    
    //------------------------------------------------
    for(G4int i=0; i<5; i++)
    {
        for(G4int k = 0; k<CAKE_TotalTimeSamples; k++)
        {
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<8; l++)
                {
                    for(G4int m=0; m<3; m++)
                    {
                        CAKE_AA[i][j][l][m][k] = 0;
                    }
                }
            }
        }
    }
    
    for(G4int i=0; i<3; i++)
    {
        for (G4int k=0; k<PADDLE_TotalTimeSamples; k++)
        {
            PADDLE_Trig[i] = 0;
            
            PADDLE_EDep[i][k] = 0;
            PADDLE_TOF[i][k] = 0;
            
            PADDLE_EWpositionX[i][k] = 0;
            PADDLE_EWpositionY[i][k] = 0;
            
            PADDLE_positionX[i][k] = 0;
            PADDLE_positionY[i][k] = 0;
        }
    }
    
    for(G4int j=0; j<4; j++)
    {
        for (G4int k=0; k<hit_buffersize; k++)
        {
            if(j==0) VDC_Observables[j][k] = -1;
            else{VDC_Observables[j][k] = 0;}
        }
    }
    
    ////    Input Variables
    InputDist[0] = 0;
    InputDist[1] = 0;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Accumulate statistics
    //
    
    // get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    
    ////////////////////////////////////////////////////////
    ////            RECOIL EXCITATION ENERGY
    ////////////////////////////////////////////////////////
    
    //recoilExcitationEnergy = G4RandGauss::shoot(recoilExcitationEnergy, 0.080*(1.0/2.35));
    
    
    
    ////////////////////////////////////////////////////////
    //
    //                CAKE DETECTOR ARRAY
    //
    ////////////////////////////////////////////////////////
    
    for(G4int i=0; i<5; i++)
    {
        for(G4int k=0; k<CAKE_TotalTimeSamples; k++)
        {
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<8; l++)
                {
                    ////    Processing the Energies of the CAKE Detectors
                    //CAKE_AA[i][j][l][0][k] = G4RandGauss::shoot(CAKE_AA[i][j][l][0][k], (0.010/2.3548));
                    //if(CAKE_AA[i][j][l][0][k] >= G4RandGauss::shoot(CAKE_AA_ThresholdEnergy, 0.1))
                    
                    //if(CAKE_AA[i][j][l][0][k] >= CAKE_AA_ThresholdEnergy)
                    if(CAKE_AA[i][j][l][0][k] > 0.0)
                    {
                        ////    Gaussian Smearing
                        //CAKE_AA[i][j][l][0][k] = G4RandGauss::shoot(CAKE_AA[i][j][l][0][k], 0.040*(1.0/2.35));
                        CAKE_AA[i][j][l][0][k] = abs(G4RandGauss::shoot(CAKE_AA[i][j][l][0][k], 0.040*(1.0/2.35)));

                        ////    Gaussian Smearing
                        //CAKE_AA[i][j][l][0][k] = G4RandGauss::shoot(CAKE_AA[i][j][l][0][k], (((double) j*0.004)) + 0.040)*(1.0/2.35);

                        ////      Counts versus Energy for each CAKE
                        //analysisManager->FillH1(1+i, GainCAKE*CAKE_AA[i][j][l][0][k] + OffsetCAKE);
                        
                        ////      Counts versus Energy for the Entire CAKE Array
                        //analysisManager->FillH1(6, GainCAKE*CAKE_AA[i][j][l][0][k] +  OffsetCAKE);
                        
                        ////////////////////////////////////////////////////////////
                        ////                Filling DataTreeSim
                        //analysisManager->FillNtupleIColumn(0, 0, i);
                        //analysisManager->FillNtupleIColumn(0, 1, j);
                        //analysisManager->FillNtupleIColumn(0, 2, k);
                        
                        //      CAKE_No
                        analysisManager->FillNtupleIColumn(0, 0, i);
                        //      CAKE_RowNo
                        analysisManager->FillNtupleIColumn(0, 1, j);
                        //      CAKE_SectorNo
                        analysisManager->FillNtupleIColumn(0, 2, l);
                        //      Energy
                        analysisManager->FillNtupleDColumn(0, 3, CAKE_AA[i][j][l][0][k]);
                        //      Theta
                        analysisManager->FillNtupleDColumn(0, 4, CAKE_AA[i][j][l][1][k]);
                        //      Phi
                        analysisManager->FillNtupleDColumn(0, 5, CAKE_AA[i][j][l][2][k]);
                        //      Ex
                        analysisManager->FillNtupleDColumn(0, 6, recoilExcitationEnergy);
                        //      DecayModeName
                        analysisManager->FillNtupleSColumn(0, 7, decayModeName);

                        analysisManager->AddNtupleRow(0);
                        
                    }
                }
            }
        }
    }
    
    
    ////////////////////////////////////////////////////////
    //
    //          PADDLE DETECTORS - Plastic Scintillators
    //
    ////////////////////////////////////////////////////////
    
    GainPADDLE = 1.0;
    OffsetPADDLE = 0.0;
    
    
    for(G4int i=0; i<3; i++)
    {
        for (G4int k=0; k<PADDLE_TotalTimeSamples; k++)
        {
            ////              Calculating energy weighted positions
            PADDLE_positionX[i][k] = G4RandGauss::shoot(PADDLE_EWpositionX[i][k]/PADDLE_EDep[i][k], 4.8);
            PADDLE_positionY[i][k] = PADDLE_EWpositionY[i][k]/PADDLE_EDep[i][k];
            
            ////              Calculating a Gaussian Smeared Energy Deposition
            PADDLE_EDep[i][k] = G4RandGauss::shoot(PADDLE_EDep[i][k], 0.10*PADDLE_EDep[i][k]);
            
            if( PADDLE_EDep[i][k] >= G4RandGauss::shoot(PADDLE_ThresholdEnergy, 0.01*PADDLE_ThresholdEnergy))
            {
                ////////////////////////////////////////////////////////
                //      PADDLE DETECTORS - 1D, Counts versus Energy
                ////////////////////////////////////////////////////////
                
                //analysisManager->FillH1(i+7, GainPADDLE*PADDLE_EDep[i][k] + OffsetPADDLE, 1);
                
                PADDLE_TOF[i][k] = G4RandGauss::shoot(PADDLE_TOF[i][k], 0.05*PADDLE_TOF[i][k]);
                
                
                ////////////////////////////////////////////////////////////////////
                //              PADDLE DETECTORS - 2D, Position versus Energy
                ////////////////////////////////////////////////////////////////////
                //analysisManager->FillH2(i+1, PADDLE_positionX[i][k], PADDLE_positionY[i][k], PADDLE_EDep[i][k]);
                
                ////////////////////////////////////////////////////////////////////
                //              PADDLE DETECTORS - 2D, Energy versus T.O.F.
                ////////////////////////////////////////////////////////////////////
                
                //analysisManager->FillH2(i+4, PADDLE_TOF[i][k], GainPADDLE*PADDLE_EDep[i][k] + OffsetPADDLE, 1);
                
            }
        }
    }
    
    ////////////////////////////////////////////////////////
    //
    //          Initial Particle Kinetic Energy
    //
    ////////////////////////////////////////////////////////

    analysisManager->FillNtupleDColumn(0, 0, initialParticleKineticEnergy);

    ////////////////////////////////////////////////////////
    //
    //                CLOVER DETECTOR ARRAY
    //
    ////////////////////////////////////////////////////////
    int eventN_CLOVER = 0;

    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        for(G4int k=0; k<CLOVER_TotalTimeSamples; k++)
        {
            int nCrystalsTriggered = 0;
            bool triggered = false;
            
            for(G4int j=0; j<4; j++)
            {
                //  0.849257 corresponds to a 2 keV FWHM
                if(G4RandGauss::shoot(CLOVER_HPGeCrystal_EDep[i][j][k], 0.849257) >= CLOVER_HPGeCrystal_ThresholdEnergy)
                {
                    nCrystalsTriggered++;
                    triggered = true;
                    
                    CLOVER_HPGeCrystal_EDep[i][j][k] = G4RandGauss::shoot(CLOVER_HPGeCrystal_EDep[i][j][k], 1.7);
                    
                    if(Activate_CLOVER_ComptonSupression)
                    {
                        for(G4int l=0; l<CLOVER_ComptonSupression_TimeWindow; l++)
                        {
                            for(G4int m=0; m<16; m++)
                            {
                                //      COMPTON SUPRESSION - VETO CLOVER Energy Depositions in anti-coincidence with BGO Shield Energy Deposition
                                if (CLOVER_BGO_EDep[i][m][k+l] >= CLOVER_BGO_ThresholdEnergy)
                                {
                                    CLOVER_HPGeCrystal_EDepVETO[i][j][k] = true;
                                }
                            }
                        }
                        if (CLOVER_HPGeCrystal_EDepVETO[i][j][k]) CLOVER_HPGeCrystal_EDep[i][j][k] = 0;
                    }
                    
                    if(Activate_CLOVER_ADDBACK && CLOVER_HPGeCrystal_EDep[i][j][k] != 0)
                    {
                        //      ADDBACK
                        CLOVER_EDep[i][k] += CLOVER_HPGeCrystal_EDep[i][j][k];
                        
                        //      For each Clover
                        //analysisManager->FillH1(i+10, GainCLOVER*CLOVER_EDep[i][k] + OffsetCLOVER);
                        
                        //      For the Entire Clover Array
                        //analysisManager->FillH1(19, GainCLOVER*CLOVER_EDep[i][k] +  OffsetCLOVER);
                    }
                    
                    /*
                    else if(CLOVER_HPGeCrystal_EDep[i][j][k] != 0)
                    {
                        //      For each Clover
                        //analysisManager->FillH1(i+10, GainCLOVER*CLOVER_HPGeCrystal_EDep[i][j][k] + OffsetCLOVER);
                        
                        //      For the Entire Clover Array
                        //analysisManager->FillH1(19, GainCLOVER*CLOVER_HPGeCrystal_EDep[i][j][k] +  OffsetCLOVER);
                    }
                    */
                    
                    /*
                    CLOVER_Number_vec.push_back(i);
                    CLOVER_Energy_vec.push_back(CLOVER_EDep[i][k]);
                    CLOVER_DetectorTheta_vec.push_back(std::get<1>(angles_CLOVER[i]));
                    CLOVER_DetectorPhi_vec.push_back(std::get<2>(angles_CLOVER[i]));
                                    
                    eventN_CLOVER++;
                    */
                }
            }
            
            if(triggered)
            {
                CLOVER_Number_vec.push_back(i);
                CLOVER_NCrystalTriggered_vec.push_back(nCrystalsTriggered);
                CLOVER_Energy_vec.push_back(CLOVER_EDep[i][k]);
                CLOVER_DetectorTheta_vec.push_back(std::get<1>(angles_CLOVER[i]));
                CLOVER_DetectorPhi_vec.push_back(std::get<2>(angles_CLOVER[i]));
                
                eventN_CLOVER++;
            }
        }
    }
    
    if(eventN_CLOVER>0)
    {
        analysisManager->FillNtupleIColumn(0, 1, eventN_CLOVER);

        fRunAction->SetCLOVER_IDs(CLOVER_Number_vec);
        fRunAction->SetCLOVER_NCrystalsTriggered(CLOVER_NCrystalTriggered_vec);
        fRunAction->SetCLOVER_Energies(CLOVER_Energy_vec);
        fRunAction->SetCLOVER_DetectorThetas(CLOVER_DetectorTheta_vec);
        fRunAction->SetCLOVER_DetectorPhis(CLOVER_DetectorPhi_vec);
    }
    
    
    ////////////////////////////////////////////////////
    //
    //              LEPS DETECTOR ARRAY
    //
    ////////////////////////////////////////////////////
    /*
    GainLEPS = 1.0;
    OffsetLEPS = 0.0;
    bool eventTriggered_LEPS = false;
    
    for(G4int i=0; i<6; i++)
    {
        for(G4int k=0; k<LEPS_TotalTimeSamples; k++)
        {
            for(G4int j=0; j<4; j++)
            {
                if(G4RandGauss::shoot(LEPS_HPGeCrystal_EDep[i][j][k], 0.7) >= LEPS_HPGeCrystal_ThresholdEnergy)
                {
                    LEPS_HPGeCrystal_EDep[i][j][k] = abs(G4RandGauss::shoot(LEPS_HPGeCrystal_EDep[i][j][k], 1.7));
                    
                    //      ADDBACK
                    if(Activate_LEPS_ADDBACK)
                    {
                        LEPS_EDep[i][k] += LEPS_HPGeCrystal_EDep[i][j][k];
                    }
                }
            }
            
            if(Activate_LEPS_ADDBACK && LEPS_EDep[i][k] >= LEPS_HPGeCrystal_ThresholdEnergy)
            {
                analysisManager->FillNtupleIColumn(0, i, 1);
                analysisManager->FillNtupleDColumn(0, i+2, GainLEPS*LEPS_EDep[i][k] + OffsetLEPS);
                eventTriggered_LEPS = true;
            }
        }
    }
     
    if(eventTriggered_LEPS) analysisManager->AddNtupleRow(0);
    */
    
    ////////////////////////////////////////////////////
    //
    //              LaBr3Ce DETECTOR ARRAY
    //
    ////////////////////////////////////////////////////
    
    //GainLaBr3Ce = 1.0;
    //OffsetLaBr3Ce = 0.0;
    int eventN_LaBr3Ce = 0;

    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        for(G4int k=0; k<LaBr3Ce_TotalTimeSamples; k++)
        {
            if(G4RandGauss::shoot(LaBr3Ce_EDep[i][k], 0.7) >= LaBr3Ce_LaBr3CeCrystal_ThresholdEnergy)
            {
                //------------------------------------------------
                LaBr3Ce_Number_vec.push_back(i);
                LaBr3Ce_DetectorTheta_vec.push_back(std::get<1>(angles_ALBA_LaBr3Ce[i]));
                LaBr3Ce_DetectorPhi_vec.push_back(std::get<2>(angles_ALBA_LaBr3Ce[i]));
                
                //------------------------------------------------
                LaBr3Ce_EWpositionX[i][k] *= (1.0/LaBr3Ce_EDep[i][k]);
                LaBr3Ce_EWpositionY[i][k] *= (1.0/LaBr3Ce_EDep[i][k]);
                LaBr3Ce_EWpositionZ[i][k] *= (1.0/LaBr3Ce_EDep[i][k]);
                
                LaBr3Ce_xPos_vec.push_back(LaBr3Ce_EWpositionX[i][k]);
                LaBr3Ce_yPos_vec.push_back(LaBr3Ce_EWpositionY[i][k]);
                LaBr3Ce_zPos_vec.push_back(LaBr3Ce_EWpositionZ[i][k]);
                
                //------------------------------------------------
                double normVector, theta, phi;
                
                normVector = pow(pow(LaBr3Ce_EWpositionX[i][k],2) + pow(LaBr3Ce_EWpositionY[i][k],2) + pow(LaBr3Ce_EWpositionZ[i][k],2) , 0.5);
                theta = acos(LaBr3Ce_EWpositionZ[i][k]/normVector)/deg;
                
                if(LaBr3Ce_EWpositionX[i][k]==0)
                {
                    if(LaBr3Ce_EWpositionY[i][k]==0) phi = 0;
                    if(LaBr3Ce_EWpositionY[i][k]>0) phi = 90;
                    if(LaBr3Ce_EWpositionY[i][k]<0) phi = 270;
                }
                else
                {
                    phi = atan(LaBr3Ce_EWpositionY[i][k]/LaBr3Ce_EWpositionX[i][k])/deg;
                    
                    if(LaBr3Ce_EWpositionX[i][k]>0 && LaBr3Ce_EWpositionY[i][k]>0) phi = phi; // deg
                    if(LaBr3Ce_EWpositionX[i][k]<0 && LaBr3Ce_EWpositionY[i][k]>0) phi = phi + 180.; // deg
                    if(LaBr3Ce_EWpositionX[i][k]<0 && LaBr3Ce_EWpositionY[i][k]<0) phi = phi + 180.; // deg
                    if(LaBr3Ce_EWpositionX[i][k]>0 && LaBr3Ce_EWpositionY[i][k]<0) phi = phi + 360.; // deg
                }

                LaBr3Ce_Theta_vec.push_back(theta);
                LaBr3Ce_Phi_vec.push_back(phi);

                //------------------------------------------------
                LaBr3Ce_EDep[i][k] = abs(G4RandGauss::shoot(LaBr3Ce_EDep[i][k], 1.7));
                LaBr3Ce_Energy_vec.push_back(LaBr3Ce_EDep[i][k]);

                eventN_LaBr3Ce++;
            }
        }
    }
    
    if(eventN_LaBr3Ce>0)
    {
        analysisManager->FillNtupleIColumn(0, 7, eventN_LaBr3Ce);

        fRunAction->SetLaBr3Ce_IDs(LaBr3Ce_Number_vec);
        fRunAction->SetLaBr3Ce_Energies(LaBr3Ce_Energy_vec);
        fRunAction->SetLaBr3Ce_DetectorThetas(LaBr3Ce_DetectorTheta_vec);
        fRunAction->SetLaBr3Ce_DetectorPhis(LaBr3Ce_DetectorPhi_vec);
        fRunAction->SetLaBr3Ce_Thetas(LaBr3Ce_Theta_vec);
        fRunAction->SetLaBr3Ce_Phis(LaBr3Ce_Phi_vec);
        fRunAction->SetLaBr3Ce_xPos(LaBr3Ce_xPos_vec);
        fRunAction->SetLaBr3Ce_yPos(LaBr3Ce_yPos_vec);
        fRunAction->SetLaBr3Ce_zPos(LaBr3Ce_zPos_vec);
    }
    
    //--------------------------------------------------------------------------------
    //      Combined data taking for both the LaBr3Ce and CLOVER detectors
    if(eventN_LaBr3Ce>0 || eventN_CLOVER>0)
    {
        analysisManager->AddNtupleRow(0);
    }
    

    
    ////////////////////////////////////////////////////////
    //
    //                VDC, VERTICAL DRIFT CHAMBER
    //
    ////////////////////////////////////////////////////////
    
    ////    VDC 1
    RayTrace(0, 0);     //RayTrace(VDCNo, XU_Wireplane)
    RayTrace(0, 1);
    CalcYFP(0);
    //G4cout << "Here is the Xpos[0] (VDC1)     -->     "<< Xpos[0] << G4endl;
    //G4cout << "Here is the ThetaFP[0] (VDC1)     -->     "<< ThetaFP[0] << G4endl;
    
    ////    VDC 2
    RayTrace(1, 0);     //RayTrace(VDCNo, XU_Wireplane)
    RayTrace(1, 1);
    CalcYFP(1);
    //G4cout << "Here is the Xpos[1] (VDC2)     -->     "<< Xpos[1] << G4endl;
    //G4cout << "Here is the ThetaFP[1] (VDC2)     -->     "<< ThetaFP[1] << G4endl;
    ////    Calculating the "True" positions through which the primary particle traverses the Wirechambers
    
    /*
    for(G4int i=0; i<4; i++)
    {
        ////    This condition checks if the PRE and POST points for each wireplane have been accounted for
        ////    PRE point, the last step point before traversing Wireplane
        ////    POST point, the first step point after traversing Wireplane
        
        if(WireplaneTraversePos[i][0][2]<0. && WireplaneTraversePOST[i])
        {
            //G4cout << "Here we are in the Event Action" << G4endl;
            
            ////    j==0 and j==1 searches for the X and Y components of the local position of traversal accross each wireplane
            for(G4int j=0; j<2; j++)
            {
                a = (WireplaneTraversePos[i][1][j] - WireplaneTraversePos[i][0][j])/(WireplaneTraversePos[i][1][2] - WireplaneTraversePos[i][0][2]);
                b = WireplaneTraversePos[i][0][j] - a*WireplaneTraversePos[i][0][2];
                
                //WireplaneTraversePos[i][2][j] = b;
                WireplaneTraversePos[i][2][j] = b + (936./2); // Offset to be aligned with the Raytraced position
                
            }
        }
    }
    */
    
    
    ////////////////////////////////////////////////////
    ////                                            ////
    ////                DataTreeSim                 ////
    ////                                            ////
    ////////////////////////////////////////////////////
    /*
    analysisManager->FillNtupleDColumn(0, 0, Xpos[0]);
    analysisManager->FillNtupleDColumn(0, 1, Y[0]);
    analysisManager->FillNtupleDColumn(0, 2, ThetaFP[0]);
    analysisManager->FillNtupleDColumn(0, 3, ThetaSCAT[0]);
    
    analysisManager->FillNtupleDColumn(0, 4, Xpos[1]);
    analysisManager->FillNtupleDColumn(0, 5, Y[1]);
    analysisManager->FillNtupleDColumn(0, 6, ThetaFP[1]);
    analysisManager->FillNtupleDColumn(0, 7, ThetaSCAT[1]);
    
    ////    Points of traversal
    ////    VDC1
    analysisManager->FillNtupleDColumn(0, 8, WireplaneTraversePos[0][2][0]);
    analysisManager->FillNtupleDColumn(0, 9, WireplaneTraversePos[0][2][1]);
    analysisManager->FillNtupleDColumn(0, 10, WireplaneTraversePos[1][2][0]);
    analysisManager->FillNtupleDColumn(0, 11, WireplaneTraversePos[1][2][1]);
    
    ////    VDC2
    analysisManager->FillNtupleDColumn(0, 12, WireplaneTraversePos[2][2][0]);
    analysisManager->FillNtupleDColumn(0, 13, WireplaneTraversePos[2][2][1]);
    analysisManager->FillNtupleDColumn(0, 14, WireplaneTraversePos[3][2][0]);
    analysisManager->FillNtupleDColumn(0, 15, WireplaneTraversePos[3][2][1]);
    
    
    analysisManager->AddNtupleRow(0);
    */
    
    
    
    
    
    ////////////////////////////////////////////////////////////////
    ////                                                        ////
    ////                    GEOMETRY ANALYSIS                   ////
    ////                                                        ////
    ////////////////////////////////////////////////////////////////
    
    if(GA_MODE)
    {
        if(GA_GenInputVar)
        {
            
            double theta_projX = atan((tan(InputDist[0]*deg)/rad)*rad*(cos(InputDist[1]*deg)/rad)*rad)/deg; // deg
            double theta_projY = atan((tan(InputDist[0]*deg)/rad)*rad*(sin(InputDist[1]*deg)/rad)*rad)/deg; // deg

            //if(InputDist[1]>180.0) theta_projY = -theta_projY;
                
            ////////////////////////////////
            ////    Input Variables
            analysisManager->FillNtupleDColumn(2, 0, InputDist[0]);
            analysisManager->FillNtupleDColumn(2, 1, InputDist[1]);
            analysisManager->FillNtupleDColumn(2, 2, theta_projX);
            analysisManager->FillNtupleDColumn(2, 3, theta_projY);
            

            
            analysisManager->AddNtupleRow(2);
            
            //G4cout << "Here is the value of InputDist[0]:    -->     " << InputDist[0] << G4endl;
            //G4cout << "Here is the value of InputDist[1]:    -->     " << InputDist[1] << G4endl;
        }
        
        
        ////////////////////////////////////////////////////////////
        ////    Creating Distribution txt file for Data Sorting
        
        if(GA_GenAngDist)
        {
            if(GA_CAKE)
            {
                if(evtNb==0)
                {
                    sprintf(filenameV,"K600Veridical_CAKE_AngDist.txt");
                    fileNameHolder = string(filenameV);
                    fileV_MMM.open(fileNameHolder, std::ios_base::app);
                }
                
                for(G4int i=0; i<512; i++)  //  i<640 for all 5 silicons
                {
                    if(GA_CAKE_AA[i][0]!=0)
                    {
                        CAKE_No = i/128;
                        CAKE_RowNo = (i - (CAKE_No*128))/8;
                        CAKE_SectorNo = (i - (CAKE_No*128))%8;
                        
                        ////////////////////////////////////////////////////////////
                        ////            Filling GeometryAnalysisTree
                        
                        analysisManager->FillNtupleIColumn(1, 0, CAKE_No);
                        analysisManager->FillNtupleIColumn(1, 1, CAKE_RowNo);
                        analysisManager->FillNtupleIColumn(1, 2, CAKE_SectorNo);
                        
                        //      Theta
                        analysisManager->FillNtupleDColumn(1, 3, GA_CAKE_AA[i][1]);
                        //      Phi
                        analysisManager->FillNtupleDColumn(1, 4, GA_CAKE_AA[i][2]);
                        
                        analysisManager->AddNtupleRow(1);
                        
                        fileV_MMM << CAKE_No << "    " << CAKE_RowNo << "    " << CAKE_SectorNo << "    " << GA_CAKE_AA[i][1] << "    " << GA_CAKE_AA[i][2] << endl;
                        
                    }
                }
                
                ////    Closing the Distribution File
                if((evtNb+1)%GA_numberOfEvents==0 )
                {
                    fileV_MMM.close();
                }

            }
            
            if(GA_W1)
            {
                if(evtNb==0)
                {
                    sprintf(filenameV,"K600Veridical_W1_AngDist.txt");
                    fileNameHolder = string(filenameV);
                    fileV_MMM.open(fileNameHolder, std::ios_base::app);
                }
                
                for(G4int i=0; i<1024; i++)  //
                {
                    if(GA_W1_AA[i][0]!=0)
                    {
                        W1_No = i/256;
                        W1_RowNo = (i - (W1_No*256))/16;
                        W1_ColumnNo = (i - (W1_No*256))%16;
                        
                        ////////////////////////////////////////////////////////////
                        ////            Filling GeometryAnalysisTree
                        
                        analysisManager->FillNtupleIColumn(1, 0, W1_No);
                        analysisManager->FillNtupleIColumn(1, 1, W1_RowNo);
                        analysisManager->FillNtupleIColumn(1, 2, W1_ColumnNo);
                        
                        //      Theta
                        analysisManager->FillNtupleDColumn(1, 3, GA_W1_AA[i][1]);
                        //      Phi
                        analysisManager->FillNtupleDColumn(1, 4, GA_W1_AA[i][2]);
                        
                        analysisManager->AddNtupleRow(1);
                        
                        fileV_MMM << W1_No << "    " << W1_RowNo << "    " << W1_ColumnNo << "    " << GA_W1_AA[i][1] << "    " << GA_W1_AA[i][2] << endl;
                        
                    }
                }
                
                ////    Closing the Distribution File
                if((evtNb+1)%GA_numberOfEvents==0 )
                {
                    fileV_MMM.close();
                }
                
            }

            
            
        }
        
        
        
        
        /////////////////////////////////////////////////////////////////////
        ////    Writing the Output files for the GEOMETRY ANALYSIS mode
        if(evtNb == (GA_numberOfEvents-1))
        {
            G4double av_xPos, av_yPos, av_zPos;
            G4double normVector, theta, phi, solidAngle;
            
            if(GA_CAKE)
            {
                // append text file
                G4String fileName1 = "K600SimOutput_CAKE.txt";
                G4String fileName2 = "K600SimOutput_CAKE.h";
                
                std::ofstream file1;
                file1.open(fileName1, std::ios_base::app);
                
                std::ofstream file2;
                file2.open(fileName2, std::ios_base::app);
                
                file1 << "(CAKE NUMBER)  (ROW NUMBER)  (SECTOR NUMBER)  (THETA)      (PHI)      (SOLID ANGLE)"<< endl;
                file2 << "                  "<< endl;
                file2 << "Double_t GA_CAKE[5][16][8][3];"<< endl;
                file2 << "                  "<< endl;
                file2 << "void initialize_GA()"<< endl;
                file2 << "{"<< endl;
                
                G4double GA_numberOfEvents_double = GA_numberOfEvents;
                
                ////    For the Silicon Array
                //for(G4int i=0; i<640; i++)  //  i<640 for all 5 silicons
                for(G4int i=0; i<512; i++)  //  i<640 for all 5 silicons
                {
                    if(i%128 == 0) file1 << "   " << endl;
                    
                    CAKE_No = i/128;
                    CAKE_RowNo = (i - (CAKE_No*128))/8;
                    CAKE_SectorNo = (i - (CAKE_No*128))%8;
                    
                    av_xPos = GA_CAKE_AA_stor[i][0]/GA_CAKE_AA_stor[i][3];
                    av_yPos = GA_CAKE_AA_stor[i][1]/GA_CAKE_AA_stor[i][3];
                    av_zPos = GA_CAKE_AA_stor[i][2]/GA_CAKE_AA_stor[i][3];
                    
                    normVector = pow(pow(av_xPos,2) + pow(av_yPos,2) + pow(av_zPos,2) , 0.5);
                    theta = acos(av_zPos/normVector)/deg;
                    //solidAngle = (GA_CAKE_AA_stor[i][3]/GA_numberOfEvents_double);
                    
                    ////    The 0.5 factor is correct for a biased calculation where the primary particle vector only spans 1 hemisphere
                    solidAngle = (0.5)*(GA_CAKE_AA_stor[i][3]/GA_numberOfEvents_double);
                    
                    /*
                     G4cout << "" << G4endl;
                     G4cout << "Here is the GA_CAKE_AA_stor[3] value     -->     "<< GA_CAKE_AA_stor[i][3] << G4endl;
                     G4cout << "Here is the GA_numberOfEvents_double value     -->     "<< GA_numberOfEvents_double << G4endl;
                     G4cout << "Here is the solidAngle value     -->     "<< solidAngle << G4endl;
                     G4cout << "" << G4endl;
                     */
                    
                    if(av_xPos==0)
                    {
                        if(av_yPos==0) phi = 0;
                        if(av_yPos>0) phi = 90;
                        if(av_yPos<0) phi = 270;
                    }
                    else
                    {
                        phi = atan(av_yPos/av_xPos)/deg;
                        
                        if(av_xPos>0 && av_yPos>0) phi = phi; // deg
                        if(av_xPos<0 && av_yPos>0) phi = phi + 180.; // deg
                        if(av_xPos<0 && av_yPos<0) phi = phi + 180.; // deg
                        if(av_xPos>0 && av_yPos<0) phi = phi + 360.; // deg
                    }
                    
                    ////////////////////////////////////////////////////////////////////////////////////////
                    file1 << CAKE_No << ",              " << CAKE_RowNo << ",            " << CAKE_SectorNo << ",               " << theta << ",     " << phi <<  ",   " << solidAngle << endl;
                    
                    file2 << "    GA_CAKE[" << CAKE_No << "][" << CAKE_RowNo << "][" << CAKE_SectorNo << "][0]=" << theta << ";   GA_CAKE[" << CAKE_No << "][" << CAKE_RowNo << "][" << CAKE_SectorNo << "][1]=" << phi << ";   GA_CAKE[" << CAKE_No << "][" << CAKE_RowNo << "][" << CAKE_SectorNo << "][2]=" << solidAngle << ";" << endl;
                    ////////////////////////////////////////////////////////////////////////////////////////
                    
                }
                
                file2 << "}"<< endl;
                
                file1.close();
                file2.close();

            }
            
            if(GA_W1)
            {
                // append text file
                G4String fileName1 = "K600SimOutput_W1.txt";
                G4String fileName2 = "K600SimOutput_W1.h";
                
                std::ofstream file1;
                file1.open(fileName1, std::ios_base::app);
                
                std::ofstream file2;
                file2.open(fileName2, std::ios_base::app);
                
                file1 << "(W1 NUMBER)  (ROW NUMBER)  (SECTOR NUMBER)  (THETA)      (PHI)      (SOLID ANGLE)"<< endl;
                file2 << "                  "<< endl;
                file2 << "Double_t GA_W1[4][16][16][3];"<< endl;
                file2 << "                  "<< endl;
                file2 << "void initialize_GA()"<< endl;
                file2 << "{"<< endl;
                
                G4double GA_numberOfEvents_double = GA_numberOfEvents;
                
                ////    For the Silicon Array
                for(G4int i=0; i<1024; i++)  //  i<640 for all 5 silicons
                {
                    if(i%256 == 0) file1 << "   " << endl;
                    
                    W1_No = i/256;
                    W1_RowNo = (i - (W1_No*256))/16;
                    W1_ColumnNo = (i - (W1_No*256))%16;
                    
                    av_xPos = GA_W1_AA_stor[i][0]/GA_W1_AA_stor[i][3];
                    av_yPos = GA_W1_AA_stor[i][1]/GA_W1_AA_stor[i][3];
                    av_zPos = GA_W1_AA_stor[i][2]/GA_W1_AA_stor[i][3];
                    
                    normVector = pow(pow(av_xPos,2) + pow(av_yPos,2) + pow(av_zPos,2) , 0.5);
                    theta = acos(av_zPos/normVector)/deg;
                    //solidAngle = (GA_W1_AA_stor[i][3]/GA_numberOfEvents_double);
                    ////    The 0.5 factor is correct for a biased calculation where the primary particle vector only spans 1 hemisphere
                    solidAngle = (0.5)*(GA_W1_AA_stor[i][3]/GA_numberOfEvents_double);
                    
                    
                    if(av_xPos==0)
                    {
                        if(av_yPos==0) phi = 0;
                        if(av_yPos>0) phi = 90;
                        if(av_yPos<0) phi = 270;
                    }
                    else
                    {
                        phi = atan(av_yPos/av_xPos)/deg;
                        
                        if(av_xPos>0 && av_yPos>0) phi = phi; // deg
                        if(av_xPos<0 && av_yPos>0) phi = phi + 180.; // deg
                        if(av_xPos<0 && av_yPos<0) phi = phi + 180.; // deg
                        if(av_xPos>0 && av_yPos<0) phi = phi + 360.; // deg
                    }
                    
                    ////////////////////////////////////////////////////////////////////////////////////////
                    file1 << W1_No << ",              " << W1_RowNo << ",            " << W1_ColumnNo << ",               " << theta << ",     " << phi <<  ",   " << solidAngle << endl;
                    
                    file2 << "    GA_W1[" << W1_No << "][" << W1_RowNo << "][" << W1_ColumnNo << "][0]=" << theta << ";   GA_W1[" << W1_No << "][" << W1_RowNo << "][" << W1_ColumnNo << "][1]=" << phi << ";   GA_W1[" << W1_No << "][" << W1_RowNo << "][" << W1_ColumnNo << "][2]=" << solidAngle << ";" << endl;
                    ////////////////////////////////////////////////////////////////////////////////////////
                    
                }
                
                file2 << "}"<< endl;
                
                file1.close();
                file2.close();
                
            }

            
        }
    }
    
    
    
    
    
    ///////////////////////////////////////////////////////
    //          PRINTING the VDC Values
    /*
     for (G4int k=0; k<hit_buffersize; k++)
     {
     if(VDC_Observables[0][k]>=0 && VDC_Observables[1][k]>20)
     {
     G4cout << "                         " << G4endl;
     G4cout << "Here is the k (index) value     -->     "<< k << G4endl;
     G4cout << "Here is the VDC_Observables[0][k], Channel ID     -->     "<< VDC_Observables[0][k] << G4endl;
     G4cout << "Here is the VDC_Observables[1][k], Edep     -->     "<< VDC_Observables[1][k] << G4endl;
     G4cout << "Here is the VDC_Observables[2][k], EW_zpos     -->     "<< VDC_Observables[2][k] << G4endl;
     G4cout << "Here is the VDC_Observables[3][k], EW_t     -->     "<< VDC_Observables[3][k] << G4endl;
     G4cout << "Here is the VDC_Observables[2][k], zpos     -->     "<< VDC_Observables[2][k]/VDC_Observables[1][k] << G4endl;
     G4cout << "Here is the VDC_Observables[3][k], t     -->     "<< VDC_Observables[3][k]/VDC_Observables[1][k] << G4endl;
     }
     }
     */
    
    /*
     for (G4int k=0; k<hit_buffersize; k++)
     {
     if(VDC_Observables[1][k]>20 && (VDC_Observables[0][k]>=0) && (VDC_Observables[0][k]<=142))
     {
     G4cout << "Here is the VDC_Observables[0][k] value     -->     "<< VDC_Observables[0][k] << G4endl;
     
     }
     }
     G4cout << "Here is the Upos[0] value     -->     "<< Upos[0] << G4endl;
     */
    
    
    
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
