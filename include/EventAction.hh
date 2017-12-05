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

#ifndef EventAction_h
#define EventAction_h 1

#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "RunAction.hh"

#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////
//                          OPERATION MODES
//////////////////////////////////////////////////////////////////////////

///////////////     GEOMETRY ANALYSIS
const G4bool        GA_MODE = true;
const G4bool        GA_CAKE = true;
const G4bool        GA_W1 = false;
const G4bool        GA_LineOfSightMODE = true;
const G4int         GA_numberOfEvents = 100000000;

const G4bool        GA_GenInputVar = true;
const G4bool        GA_GenAngDist = true;
//const G4int         GA_GenAngDist_buffer = 5000;


//////////////////////////////////////////////////////////////////////////

///////////////     CAKE Detectors - PIXIE16 Sampling     ///////////////////
const G4double      CAKE_SamplingTime = 200000; // ns
const G4int         CAKE_TotalTimeSamples = 1; //
const G4double      CAKE_TotalSampledTime = CAKE_SamplingTime * CAKE_TotalTimeSamples; // ns
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///////////////     VDC Detectors       ///////////////////
const G4int         hit_buffersize = 100;
const G4double      VDC_SamplingTime = 10; // ns
const G4int         VDC_TotalTimeSamples = 15; //
const G4double      VDC_TotalSampledTime = VDC_SamplingTime * VDC_TotalTimeSamples; // ns
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///////////////     PADDLE Detectors - Analogue Sampling    ///////////////////
const G4double      PADDLE_SamplingTime = 10; // ns
const G4int         PADDLE_TotalTimeSamples = 15; //
const G4double      PADDLE_TotalSampledTime = PADDLE_SamplingTime * PADDLE_TotalTimeSamples; // ns
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///////////////     CLOVER Detectors - PIXIE16 Sampling     ///////////////////
const G4bool        Activate_CLOVER_ADDBACK = true;
const G4bool        Activate_CLOVER_ComptonSupression = false;

const G4double      CLOVER_SamplingTime = 10; // ns
const G4int         CLOVER_TotalTimeSamples = 10; //
const G4double      CLOVER_TotalSampledTime = CLOVER_SamplingTime * CLOVER_TotalTimeSamples; // ns
const G4int         CLOVER_ComptonSupression_TimeWindow = 3; // Amount of CLOVER Time Samples

///////////////     CLOVER BGO Anti-Compton Shield - PIXIE16 Sampling    ///////////////////
const G4double      CLOVER_Shield_BGO_SamplingTime = CLOVER_SamplingTime; // ns
const G4int         CLOVER_Shield_BGO_TotalTimeSamples = CLOVER_TotalTimeSamples + CLOVER_ComptonSupression_TimeWindow; //
const G4double      CLOVER_Shield_BGO_TotalSampledTime = CLOVER_Shield_BGO_SamplingTime * CLOVER_Shield_BGO_TotalTimeSamples; // ns
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

///////////////     LEPS Detectors - PIXIE16 Sampling     ///////////////////
const G4bool        Activate_LEPS_ADDBACK = true;

const G4double      LEPS_SamplingTime = 10; // ns
const G4int         LEPS_TotalTimeSamples = 10; //
const G4double      LEPS_TotalSampledTime = LEPS_SamplingTime * LEPS_TotalTimeSamples; // ns

///////////////     LaBr3Ce Detectors - PIXIE16 Sampling     ///////////////////
const G4double      LaBr3Ce_SamplingTime = 10; // ns
const G4int         LaBr3Ce_TotalTimeSamples = 100; //
const G4double      LaBr3Ce_TotalSampledTime = LaBr3Ce_SamplingTime * LaBr3Ce_TotalTimeSamples; // ns

///////////////     VDC Signal Wires - Energy Threshold     ///////////////////
const G4double      VDC1_U_WIRE_ThresholdEnergy = 10.;   // keV
const G4double      VDC1_X_WIRE_ThresholdEnergy = 10.;   // keV

const G4double      VDC2_U_WIRE_ThresholdEnergy = 10.;   // keV
const G4double      VDC2_X_WIRE_ThresholdEnergy = 10.;   // keV

///////////////     CAKE - Energy Threshold     ///////////////////
const G4double      CAKE_AA_ThresholdEnergy = .5;   // MeV

///////////////     CLOVER - Energy Threshold     ///////////////////
const G4double      CLOVER_HPGeCrystal_ThresholdEnergy = 6.;   // keV

///////////////     CLOVER BGO Anti-Compton Shield - Energy Threshold     ///////////////////
const G4double      CLOVER_BGO_ThresholdEnergy = 5.;  //keV

///////////////     LEPS - Energy Threshold     ///////////////////
const G4double      LEPS_HPGeCrystal_ThresholdEnergy = 6.;   // keV

///////////////     LaBr3Ce - Energy Threshold     ///////////////////
const G4double      LaBr3Ce_LaBr3CeCrystal_ThresholdEnergy = 6.;   // keV

///////////////     PADDLE, Plastic Scintillators - Energy Threshold     ///////////////////
const G4double      PADDLE_ThresholdEnergy = 0.5;  //  MeV

///////////////     Average particles per packet, (from beam intensity and frequency)     ///////
const G4bool        Activate_CyclotronBeam_Timing = false;
const G4int         Particles_per_Bunch = 100;  // Particles per Bunch



////////////////////////////////////////////////////
////    Variables for calculating VDC observables

const G4double a0 = -1.01703, a1 = -6.25653e-05, a2 = 0.;
const G4double b0 = 33.6679, b1 = -0.0025703, b2 = 0.;

//  Variables for CalcYFP
const G4double sinThetaU = 0.766044443;
const G4double tanThetaU = 1.191753593;


class EventAction : public G4UserEventAction
{
public:
    EventAction(RunAction *runAction, DetectorConstruction *detectorConstruction);
    virtual ~EventAction();
    
    RunAction*  fRunAction;

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddAbs(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);
    
    G4int       evtNb;
    
    G4double initialParticleKineticEnergy;
    void SetInitialParticleKineticEnergy(G4double a) {initialParticleKineticEnergy = a;};
    
    /////////////////////
    //      CAKE
    G4double GainCAKE;
    G4double OffsetCAKE;
    
    G4double    CAKE_AA[5][16][8][3][CAKE_TotalTimeSamples];
    //  First index designates the CAKENo
    //  Second index designates the CAKE_RowNo
    //  Third index designates the CAKE_SectorNo
    //  Fourth index designates:
    //  0 -> Energy, 1 -> Theta (of first interaction), 2 -> Phi (of first interaction)
    
    
    void FillVar_CAKE_AA(G4int i, G4int j, G4int l, G4int m, G4int k, G4double a)
    {CAKE_AA[i][j][l][m][k] += a;};
    
    void SetVar_CAKE_AA(G4int i, G4int j, G4int l, G4int m, G4int k, G4double a)
    {CAKE_AA[i][j][l][m][k] = a;};
    
    G4double GetVar_CAKE_AA(G4int i, G4int j, G4int l, G4int m, G4int k)
    {return CAKE_AA[i][j][l][m][k];};
    
    
    ////////////////////////
    //      CLOVERS
    //G4int nCLOVER_detectors;
    G4double GainCLOVER;
    G4double OffsetCLOVER;
    
    std::vector<int> CLOVER_Number_vec;
    std::vector<int> CLOVER_NCrystalTriggered_vec;
    std::vector<double> CLOVER_Energy_vec;
    std::vector<double> CLOVER_DetectorTheta_vec;
    std::vector<double> CLOVER_DetectorPhi_vec;
    
    G4double    CLOVER_HPGeCrystal_EDep[numberOf_CLOVER][4][CLOVER_TotalTimeSamples];
    G4bool      CLOVER_HPGeCrystal_EDepVETO[numberOf_CLOVER][4][CLOVER_TotalTimeSamples];
    G4double    CLOVER_EDep[numberOf_CLOVER][CLOVER_TotalTimeSamples];
    
    void AddEnergyCLOVER_HPGeCrystal(G4int i, G4int j, G4int k, G4double a)	{CLOVER_HPGeCrystal_EDep[i][j][k] += a;};
    
    std::vector<std::tuple<int, double, double>> angles_CLOVER;
    
    /////////////////////////////////////////
    //      CLOVER Shield BGO Crystals
    G4double    CLOVER_BGO_EDep[numberOf_CLOVER_Shields][16][CLOVER_Shield_BGO_TotalTimeSamples+CLOVER_ComptonSupression_TimeWindow];
    
    void AddEnergyBGODetectors(G4int i, G4int j, G4int k, G4double a)	{CLOVER_BGO_EDep[i][j][k] += a;};
    
    
    
    ////////////////////////
    //      LEPS
    G4double GainLEPS;
    G4double OffsetLEPS;
    // Previous versions, moved declaration to EventAction.cc constructor
    //G4double GainLEPS = 1.0;
    //G4double OffsetLEPS = 0.0;
    
    G4double    LEPS_HPGeCrystal_EDep[6][4][LEPS_TotalTimeSamples];
    G4double    LEPS_EDep[6][LEPS_TotalTimeSamples];
    
    void AddEnergyLEPS_HPGeCrystals(G4int i, G4int j, G4int k, G4double a)	{LEPS_HPGeCrystal_EDep[i][j][k] += a;};
    
    ////////////////////////
    //      LaBr3Ce
    //G4int nLaBr3Ce_detectors;
    G4double GainLaBr3Ce;
    G4double OffsetLaBr3Ce;

    std::vector<int> LaBr3Ce_Number_vec;
    std::vector<double> LaBr3Ce_Energy_vec;
    std::vector<double> LaBr3Ce_DetectorTheta_vec;
    std::vector<double> LaBr3Ce_DetectorPhi_vec;
    std::vector<double> LaBr3Ce_Theta_vec;
    std::vector<double> LaBr3Ce_Phi_vec;
    std::vector<double> LaBr3Ce_xPos_vec;
    std::vector<double> LaBr3Ce_yPos_vec;
    std::vector<double> LaBr3Ce_zPos_vec;

    G4double    LaBr3Ce_EDep[20][LaBr3Ce_TotalTimeSamples];
    G4double    LaBr3Ce_EWpositionX[20][LaBr3Ce_TotalTimeSamples];
    G4double    LaBr3Ce_EWpositionY[20][LaBr3Ce_TotalTimeSamples];
    G4double    LaBr3Ce_EWpositionZ[20][LaBr3Ce_TotalTimeSamples];
    
    void AddEnergyLaBr3Ce_LaBr3CeCrystal(G4int i, G4int k, G4double a)	{LaBr3Ce_EDep[i][k] += a;};
    void AddEWpositionX_LaBr3Ce_LaBr3CeCrystal(G4int i, G4int k, G4double a) {LaBr3Ce_EWpositionX[i][k] += a;};
    void AddEWpositionY_LaBr3Ce_LaBr3CeCrystal(G4int i, G4int k, G4double a) {LaBr3Ce_EWpositionY[i][k] += a;};
    void AddEWpositionZ_LaBr3Ce_LaBr3CeCrystal(G4int i, G4int k, G4double a) {LaBr3Ce_EWpositionZ[i][k] += a;};

    std::vector<std::tuple<int, double, double>> angles_ALBA_LaBr3Ce;
    
    /////////////////////////////////////////
    //          PADDLE DETECTORS
    G4double GainPADDLE;
    G4double OffsetPADDLE;
    
    G4double    PADDLE_EDep[3][PADDLE_TotalTimeSamples];
    G4double    PADDLE_TOF[3][PADDLE_TotalTimeSamples];
    G4bool      PADDLE_Trig[3];
    G4int       PADDLE_numberDetTrig;
    
    //      Energy Weighted Positioning
    G4double       PADDLE_EWpositionX[3][PADDLE_TotalTimeSamples];
    G4double       PADDLE_EWpositionY[3][PADDLE_TotalTimeSamples];
    G4double       PADDLE_positionX[3][PADDLE_TotalTimeSamples];
    G4double       PADDLE_positionY[3][PADDLE_TotalTimeSamples];
    
    void AddEnergy_PADDLE(G4int i, G4int j, G4double a)	{PADDLE_EDep[i][j] += a;};
    void TagTOF_PADDLE(G4int i, G4int j, G4double a)	{PADDLE_TOF[i][j] = a;};
    void AddEWpositionX_PADDLE(G4int i, G4int j, G4double a)  {PADDLE_EWpositionX[i][j] += a;};
    void AddEWpositionY_PADDLE(G4int i, G4int j, G4double a)  {PADDLE_EWpositionY[i][j] += a;};
    void Set_PADDLE_Trig(G4int i, G4bool b) {PADDLE_Trig[i] = b;};
    G4bool Get_PADDLE_Trig(G4int i) {return PADDLE_Trig[i];};
    
    
    /////////////////////////////////////////
    //          VDC DETECTORS
    G4double GainVDC;
    G4double OffsetVDC;
    
    // VDC_Observables[i][k]
    G4double    VDC_Observables[4][hit_buffersize]; // buffer of approx. 50 possible hits
    // i==0 => CELL NUMBER, where U1:(0->142), X1:(143->340), U2:(341->483), X2:(484->681)
    // i==1 => Edep
    // i==2 => E-weighted z-position
    // i==3 => E-weighted time
    // k => cell hits (not yet verified to be valid hits)
    
    //  Variables for RayTrace
    G4double Xpos[2], Upos[2], Y[2];
    G4double ThetaFP[2];
    G4double ThetaSCAT[2];
    
    G4double EnergyThreshold;
    
    //  Variables for CalcYFP
    G4double tmp1,tmp2;
    G4double tanThetaFP;
    
    
    
    void FillVDC_Observables(G4int k, G4int channelID, G4double Edep, G4double EW_zpos, G4double EW_t)
    {
        if(VDC_Observables[0][k] == -1)
        {
            VDC_Observables[0][k] = channelID;
        }
        
        VDC_Observables[1][k] += Edep;
        VDC_Observables[2][k] += EW_zpos;
        VDC_Observables[3][k] += EW_t;
    };
    
    G4double GetVDC_ObservablesChannelID(G4int k)    {return VDC_Observables[0][k];};
    
    void RayTrace(G4int VDCNo, G4int XU_Wireplane);
    void CalcYFP(G4int VDCNo);
    
    ////    WireplaneTraversePos[A][B][C]
    ////    A -> Wireplane Number. 0,1->VDC1 and 2,3->VDC2
    ////    B -> 0: PRE point, the last step point before traversing Wireplane
    ////    B -> 1: POST point, the first step point after traversing Wireplane
    ////    B -> 2: TRAVERSAL point accross the wireplane
    ////    C -> 0, 1, 2: x, y and z positions respectively
    G4double    WireplaneTraversePos[4][3][3];
    
    void SetVDC_WireplaneTraversePos(G4int WireplaneNumber, G4int i, G4int component, G4double componentPosition)
    {
        WireplaneTraversePos[WireplaneNumber][i][component] = componentPosition;
    }
    
    ////    WireplaneTraversePOST[A]
    ////    A -> Wireplane Number
    ////    True implies that the POST point has been accounted for, False it is unnacounted for
    G4bool      WireplaneTraversePOST[4];
    
    void SetVDC_WireplaneTraversePOST(G4int WireplaneNumber, G4bool decision)
    {
        WireplaneTraversePOST[WireplaneNumber] = decision;
    }
    
    G4bool GetVDC_WireplaneTraversePOST(G4int WireplaneNumber)
    {
        return WireplaneTraversePOST[WireplaneNumber];
    }
    
    /////////////////////////////////
    //      GEOMETRY ANALYSIS
    /////////////////////////////////
    
    G4bool      GA_LineOfSight;  //  LOF -> Line of Sight
    
    G4bool  GA_GetLineOfSight()	{return GA_LineOfSight;};
    void GA_SetLineOfSight(G4bool b)	{GA_LineOfSight = b;};
    
    ////    Input Variables
    G4double    InputDist[2];
    void SetInputDist(G4int i, G4double a)	{InputDist[i] = a;};
    
    ////    CAKE
    G4int CAKE_No, CAKE_RowNo, CAKE_SectorNo;
    G4double    GA_CAKE_AA_stor[640][4];
    
    ////    W1
    G4int W1_No, W1_RowNo, W1_ColumnNo;
    G4double    GA_W1_AA_stor[1024][4];

    
    ////    TEST
    std::ofstream fileV_MMM;
    char filenameV[512];
    G4String fileNameHolder;
    
    
    ////    Angular Distribution for Data Sorting
    G4int       GA_MMM_AngDist_counter[5][16][8];
    G4double    GA_MMM_AngDist[4][16][8][2][100];
    
    ////////////////////////////////////
    ////    GA mode for CAKE
    //  First index designates channel
    //  Second index:
    //  Indices 0, 1 and 2 designates summed x, y and z positions respectively whilst an index of 3 designates the number of valid hits
    void FillGA_CAKEstor(G4int i, G4int j, G4double a)	{GA_CAKE_AA_stor[i][j] += a;};
    
    
    G4double    GA_CAKE_AA[640][3];
    //  First index designates channel
    //  Second index:
    //  An index of 0, designates the theta value of the first relevant interaction
    //  An index of 1, designates the phi value of the first relevant interaction
    //  An index of 2, designates whether the volume of interest has been hit or not: 0=>!hit, 1=>hit
    void SetGA_CAKE(G4int i, G4int j, G4double a)	{GA_CAKE_AA[i][j] = a;};
    double GetGA_CAKE(G4int i, G4int j)	{return GA_CAKE_AA[i][j];};
    

    ////////////////////////////////////
    ////    GA mode for W1
    //  First index designates channel
    //  Second index:
    //  Indices 0, 1 and 2 designates summed x, y and z positions respectively whilst an index of 3 designates the number of valid hits
    void FillGA_W1stor(G4int i, G4int j, G4double a)	{GA_W1_AA_stor[i][j] += a;};

    G4double    GA_W1_AA[1024][3];
    //  First index designates channel
    //  Second index:
    //  An index of 0, designates the theta value of the first relevant interaction
    //  An index of 1, designates the phi value of the first relevant interaction
    //  An index of 2, designates whether the volume of interest has been hit or not: 0=>!hit, 1=>hit
    void SetGA_W1(G4int i, G4int j, G4double a)	{GA_W1_AA[i][j] = a;};
    double GetGA_W1(G4int i, G4int j)	{return GA_W1_AA[i][j];};

    
    ////    Primary Generator Action Variables
    G4double recoilExcitationEnergy;
    void SetRecoilExcitationEnergy(G4double initialExcitationEnergy) {recoilExcitationEnergy = initialExcitationEnergy;};
    
    ////    Primary Generator Action Variables
    G4String decayModeName;
    void SetDecayMode(G4String string_decayMode) {decayModeName = string_decayMode;};

    
private:
    G4double  fEnergyAbs;
    G4double  fEnergyGap;
    G4double  fTrackLAbs;
    G4double  fTrackLGap;
    
    
    
    
    
};


inline void EventAction::RayTrace(G4int VDCNo, G4int XU_Wireplane)
{
    G4double signalWirePos, z_dd, sum_n=0.0, sum_x=0.0, sum_z=0.0, sum_xz=0.0, sum_x2=0.0;
    
    G4int wireChannelMin, wireChannelMax, wireOffset;
    
    ////////////////    Wire channel mapping for the case when the X wireframe is upstream of the U wireframe
    ////    VDC 1
    if(VDCNo==0 && XU_Wireplane==0) wireChannelMin = 0, wireChannelMax = 197, wireOffset = 0, EnergyThreshold = VDC1_X_WIRE_ThresholdEnergy;
    if(VDCNo==0 && XU_Wireplane==1) wireChannelMin = 198, wireChannelMax = 340, wireOffset = 143, EnergyThreshold = VDC1_U_WIRE_ThresholdEnergy;
    
    ////    VDC 2
    if(VDCNo==1 && XU_Wireplane==0) wireChannelMin = 341, wireChannelMax = 538, wireOffset = 341, EnergyThreshold = VDC2_X_WIRE_ThresholdEnergy;
    if(VDCNo==1 && XU_Wireplane==1) wireChannelMin = 539, wireChannelMax = 681, wireOffset = 484, EnergyThreshold = VDC2_U_WIRE_ThresholdEnergy;
    
    /*
     ////////////////    Wire channel mapping for the case when the U wireframe is upstream of the X wireframe
     ////    VDC 1
     if(VDCNo==0 && XU_Wireplane==0) wireChannelMin = 0, wireChannelMax = 142, wireOffset = 0, EnergyThreshold = VDC1_U_WIRE_ThresholdEnergy;
     if(VDCNo==0 && XU_Wireplane==1) wireChannelMin = 143, wireChannelMax = 340, wireOffset = 143, EnergyThreshold = VDC1_X_WIRE_ThresholdEnergy;
     
     ////    VDC 2
     if(VDCNo==1 && XU_Wireplane==0) wireChannelMin = 341, wireChannelMax = 483, wireOffset = 341, EnergyThreshold = VDC2_U_WIRE_ThresholdEnergy;
     if(VDCNo==1 && XU_Wireplane==1) wireChannelMin = 484, wireChannelMax = 681, wireOffset = 484, EnergyThreshold = VDC2_X_WIRE_ThresholdEnergy;
     */
    
    
    
    for(G4int k=0; k<hit_buffersize; k++)
    {
        if( (VDC_Observables[0][k]>=wireChannelMin) && (VDC_Observables[0][k]<=wireChannelMax) && (VDC_Observables[1][k]>EnergyThreshold) )
        {
            signalWirePos = 4.0*(VDC_Observables[0][k] - wireOffset);  // mm
            z_dd = VDC_Observables[2][k]/VDC_Observables[1][k];
            
            sum_n  += 1.0;
            sum_x  += signalWirePos;
            sum_z  += z_dd;
            sum_xz += signalWirePos*z_dd;
            sum_x2 += pow(signalWirePos,2);
        }
    }
    
    G4double a;
    G4double b;

    // Equation is of the form: z = ax + b
    a = (sum_x*sum_z-sum_n*sum_xz)/(pow(sum_x,2)-sum_n*sum_x2);
    b = (sum_x*sum_xz-sum_x2*sum_z)/(pow(sum_x,2)-sum_n*sum_x2);
    
    
    if(XU_Wireplane==0)
    {
        Upos[VDCNo]  = (-1.)*b/a; // X position at the X Wireframe, mm
    }
    
    
    if(XU_Wireplane==1)
    {
        Xpos[VDCNo]  = (-1.)*b/a; // X position at the X Wireframe, mm
        ThetaFP[VDCNo] = (-1.)*atan(a)/deg;
        //G4cout << "Here is the ThetaFP[VDCNo]     -->     "<< ThetaFP[VDCNo] << G4endl;
        ThetaSCAT[VDCNo] = (a0 + a1*Xpos[VDCNo])*ThetaFP[VDCNo] + (b0 + b1*Xpos[VDCNo]);
    }
    
}



inline void EventAction::CalcYFP(G4int VDCNo)
{
    tanThetaFP = tan(ThetaFP[VDCNo]*deg);
    
    //G4cout << "Here is the ThetaFP[VDCNo]     -->     "<< ThetaFP[VDCNo] << G4endl;
    //G4cout << "Here is the tanThetaFP     -->     "<< tanThetaFP << G4endl;
    
    /*
     G4double sinThetaU = 0.766044443; //sin(U_WIRE_ANGLE/57.2957);
     G4double tanThetaU = 1.191753593; //tan(U_WIRE_ANGLE/57.2957);
     G4double tmp1,tmp2;
     G4double tanThetaFP;
     */
    
    // for UX configuration. See RN K600 book6 p20-23
    //  tmp1=(u*tanfp+sinu*16);
    //  tmp2=sinu*tanfp;
    //  *y=(tmp1/tmp2-x)*tanu+76.27 -50;  // the 76.27 is the offset due to first u and x wires not sharing the same origin
    // the -50 is to put it around zero
    // for XU configuration
    tmp1 = (Upos[VDCNo]*tanThetaFP - sinThetaU*16);
    tmp2 = sinThetaU*tanThetaFP;
    //Y[VDCNo] = -1*((tmp1/tmp2-Xpos[VDCNo])*tanThetaU + 26.21);
    Y[VDCNo] = (Xpos[VDCNo] - (tmp1/tmp2))*tanThetaU + 35.;
    
    /*
     tmp1=(u*tanfp-sinu*16);
     tmp2=sinu*tanfp;
     *y=-1*((tmp1/tmp2-x)*tanu+26.21);
     */
    
}




inline void EventAction::AddAbs(G4double de, G4double dl) {
    fEnergyAbs += de;
    fTrackLAbs += dl;
}

inline void EventAction::AddGap(G4double de, G4double dl) {
    fEnergyGap += de;
    fTrackLGap += dl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


