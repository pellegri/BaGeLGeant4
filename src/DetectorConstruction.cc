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

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Trd.hh"
#include "Randomize.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4Isotope.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalConstants.hh"

#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"

#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "CADMesh.hh"
#include "MagneticFieldMapping.hh"
//#include "G4BlineTracer.hh"

#include "GeometryConstructionDANDELION3.hh"
#include "GeoConstruct_22_03_18.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

G4ThreadLocal G4QuadrupoleMagField* DetectorConstruction::MagneticField_K600_Q = 0;
G4ThreadLocal G4UniformMagField* DetectorConstruction::MagneticField_K600_D1 = 0;
G4ThreadLocal G4UniformMagField* DetectorConstruction::MagneticField_K600_D2 = 0;

G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_K600_Q = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_K600_D1 = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fieldManagerMagneticField_K600_D2 = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
fAbsorberPV(0), fGapPV(0), fCheckOverlaps(false), PhysiCLOVER_HPGeCrystal(0), PhysiCLOVER_Shield_BGOCrystal(0), PhysiCLOVER_Shield_PMT(0), PhysiCAKE_AA_RS(0), PhysiPADDLE(0), PhysiK600_Quadrupole(0), PhysiK600_Dipole1(0), PhysiK600_Dipole2(0), PhysiHAGAR_NaICrystal(0), PhysiHAGAR_Annulus(0), PhysiHAGAR_FrontDisc(0), Physical_LEPS_HPGeCrystal(0)
{
    WorldSize = 15.*m;
    setPreconfiguredVersion = false;
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        LaBr3Ce_Presence[i] = false;
    }
    LaBR3Ce_SetGlobalDistance = false;
    LaBR3Ce_automaticOrientation = false;
    configuration_truncatedIcosahedron_hexagons = false;
    SetupTruncatedIcosahedron();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    
    //////////////////////////////////////////////////////////
    ////                                                  ////
    ////            DETECTOR ARRAY POSITIONING            ////
    ////                                                  ////
    //////////////////////////////////////////////////////////
    
    
    ////////////////////////////
    ////    VDC SETUP
    VDC_AllPresent_Override = false;
    VDC_AllAbsent_Override = true;
    
    //  VDC 1
    VDC_Presence[0] = true;
    VDC_CentrePositionX[0] = 481.93; // cm
    VDC_CentrePositionZ[0] = 352.050; // cm
    VDC_CentrePosition[0] = G4ThreeVector(VDC_CentrePositionX[0]*cm - 200.*cm, 0., VDC_CentrePositionZ[0]*cm - 100.*cm);
    VDC_RotationY[0] = -14.03; // deg;
    VDC_rotm[0].rotateY(VDC_RotationY[0]*deg);
    
    //  VDC 2
    VDC_Presence[1] = false;
    VDC_CentrePositionX[1] = 531.23; // cm
    VDC_CentrePositionZ[1] = 337.40; // cm
    VDC_CentrePosition[1] = G4ThreeVector(VDC_CentrePositionX[1]*cm - 200.*cm, 0., VDC_CentrePositionZ[1]*cm - 100.*cm);
    VDC_RotationY[1] = -14.03; // deg;
    VDC_rotm[1].rotateY(VDC_RotationY[1]*deg);
    
    
    for(G4int i=0; i<numberOf_VDC; i++)
    {
        if(VDC_AllPresent_Override) VDC_Presence[i] = true;
        if(VDC_AllAbsent_Override) VDC_Presence[i] = false;
        if(VDC_AllPresent_Override && VDC_AllAbsent_Override) VDC_Presence[i] = false;
    }
    
    
    /////////////////////////////
    ////    PADDLE SETUP
    
    PADDLE_AllPresent_Override = false;
    PADDLE_AllAbsent_Override = true;
    
    //  PADDLE 1
    PADDLE_Presence[0] = true;
    PADDLE_CentrePositionX[0] = 513.726; // cm
    PADDLE_CentrePositionZ[0] = 324.573; // cm
    PADDLE_CentrePosition[0] = G4ThreeVector(PADDLE_CentrePositionX[0]*cm - 200.*cm, 0., PADDLE_CentrePositionZ[0]*cm - 100.*cm);
    PADDLE_RotationY[0] = -14.03; // deg
    PADDLE_rotm[0].rotateY(PADDLE_RotationY[0]*deg);
    
    //  PADDLE 2
    PADDLE_Presence[1] = false;
    PADDLE_CentrePositionX[1] = 529.38; // cm
    PADDLE_CentrePositionZ[1] = 320.34; // cm
    PADDLE_CentrePosition[1] = G4ThreeVector(PADDLE_CentrePositionX[1]*cm - 200.*cm, 0., PADDLE_CentrePositionZ[1]*cm - 100.*cm);
    PADDLE_RotationY[1] = -14.03; // deg
    PADDLE_rotm[1].rotateY(PADDLE_RotationY[1]*deg);
    
    //  PADDLE 3
    PADDLE_Presence[2] = true;
    PADDLE_CentrePositionX[2] = 529.38; // cm
    PADDLE_CentrePositionZ[2] = 320.34; // cm
    PADDLE_CentrePosition[2] = G4ThreeVector(PADDLE_CentrePositionX[2]*cm - 200.*cm, 0., PADDLE_CentrePositionZ[2]*cm - 100.*cm);
    PADDLE_RotationY[2] = -14.03; // deg
    PADDLE_rotm[2].rotateY(PADDLE_RotationY[2]*deg);
    
    for(G4int i=0; i<3; i++)
    {
        if(PADDLE_AllPresent_Override) PADDLE_Presence[i] = true;
        if(PADDLE_AllAbsent_Override) PADDLE_Presence[i] = false;
        if(PADDLE_AllPresent_Override && PADDLE_AllAbsent_Override) PADDLE_Presence[i] = false;
    }
    
    
    ////////////////////////////
    ////    CAKE SETUP
    
    CAKE_AllPresent_Override = false;
    CAKE_AllAbsent_Override = true;
    
    ////    STANDARD PR226 - Assuming target at origin
    //offset_CAKE_BeamAxis = -131.3217600; // mm
    
    ////    PR226 - Error estimation
    //offset_CAKE_BeamAxis = -131.3217600 + 3.0; // mm, Standard with the +3mm target frame offset

    ////    Phil's experiment
    //offset_CAKE_BeamAxis = -130.5217600 + 2.0; // mm

    ////    Fhumi's experiment
    offset_CAKE_BeamAxis = -145.07 + 2.0; // mm, Standard with the +3mm target frame offset
    //  14/12/16
    //  MT:     +0.0
    //  MT2:    -2.0
    //  MT3:    +2.0


    //  CAKE 1
    CAKE_Presence[0] = true;
    CAKE_AA_CentrePosition[0] = G4ThreeVector(1.9682273*mm, 0.1740583*mm, offset_CAKE_BeamAxis*mm);
    CAKE_rotm[0].rotateX(22.7683986*deg);
    CAKE_rotm[0].rotateY(-35.0661884*deg);
    CAKE_rotm[0].rotateZ((-54.1496802 + (0*72.))*deg);

    //  For Fhumi
    CAKE_rotm[0].rotateZ((192.0)*deg);
    
    ////    This positioning if for viewing the CAKE
    //CAKE_rotm[0].rotateZ(59.63*deg);
    
    //  CAKE 2
    CAKE_Presence[1] = false;
    CAKE_AA_CentrePosition[1] = G4ThreeVector(0.4535687*mm, 1.9147014*mm, offset_CAKE_BeamAxis*mm);
    CAKE_rotm[1].rotateX(22.7683986*deg);
    CAKE_rotm[1].rotateY(-35.0661884*deg);
    CAKE_rotm[1].rotateZ((-54.1496802 + (1*72.))*deg);
    
    
    //  CAKE 3
    CAKE_Presence[2] = false;
    CAKE_AA_CentrePosition[2] = G4ThreeVector(-1.6699365*mm, 1.0120638*mm, offset_CAKE_BeamAxis*mm);
    CAKE_rotm[2].rotateX(22.7683986*deg);
    CAKE_rotm[2].rotateY(-35.0661884*deg);
    CAKE_rotm[2].rotateZ((-54.1496802 + (2*72.))*deg);
    
    
    //  CAKE 4
    CAKE_Presence[3] = false;
    CAKE_AA_CentrePosition[3] = G4ThreeVector(-1.4676763*mm, -1.2864401*mm, offset_CAKE_BeamAxis*mm);
    CAKE_rotm[3].rotateX(22.7683986*deg);
    CAKE_rotm[3].rotateY(-35.0661884*deg);
    CAKE_rotm[3].rotateZ((-54.1496802 + (3*72.))*deg);
    
    
    //  CAKE 5
    //  This CAKE was not present for PR226 and therefore the CAKE_AA_CentrePosition[4] was not yet measured from SolidEdge. Still to do for COMPLETENESS.
    //  In principle, one only has to rotate the CAKE_AA_CentrePosition[] for a particular MMM about the Z-axis for all 5 MMM's in the CAKE array configuration
    CAKE_Presence[4] = false;
    CAKE_AA_CentrePosition[4] = G4ThreeVector(0.773759*mm, -1.81811*mm, offset_CAKE_BeamAxis*mm);
    CAKE_rotm[4].rotateX(22.7683986*deg);
    CAKE_rotm[4].rotateY(-35.0661884*deg);
    CAKE_rotm[4].rotateZ((-54.1496802 + (4*72.))*deg);

    
    /*
     ////   Old
     CAKE_AA_CentrePosition[4] = G4ThreeVector(0.90*mm, -1.50*mm, offset_CAKE_BeamAxis*mm);
     CAKE_rotm[4].rotateX(-37.891*deg);
     CAKE_rotm[4].rotateY(17.002*deg);
     CAKE_rotm[4].rotateZ((-38.571 + 4*72)*deg);
     */
    
    for(G4int i=0; i<numberOf_CAKE; i++)
    {
        if(CAKE_AllPresent_Override) CAKE_Presence[i] = true;
        if(CAKE_AllAbsent_Override) CAKE_Presence[i] = false;
        if(CAKE_AllPresent_Override && CAKE_AllAbsent_Override) CAKE_Presence[i] = false;
    }
    
    
    ////////////////////////////
    ////    W1 SETUP
    
    W1_AllPresent_Override = false;
    W1_AllAbsent_Override = true;
    
    offset_W1_BeamAxis = -65.9; // mm
    
    //  W1 1
    W1_Presence[0] = false;
    W1_AA_CentrePosition[0] = G4ThreeVector(0.0*mm, 34.0*mm, offset_W1_BeamAxis*mm);
    //W1_AA_CentrePosition[0] = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    W1_rotm[0].rotateX(90.0*deg);
    
    //  W1 2
    W1_Presence[1] = false;
    W1_AA_CentrePosition[1] = G4ThreeVector(-34.0*mm, 0.0*mm, offset_W1_BeamAxis*mm);
    W1_rotm[1].rotateX(90.0*deg);
    W1_rotm[1].rotateZ(90.0*deg);

    //  W1 3
    W1_Presence[2] = false;
    W1_AA_CentrePosition[2] = G4ThreeVector(0.0*mm, -34.0*mm, offset_W1_BeamAxis*mm);
    W1_rotm[2].rotateX(90.0*deg);
    W1_rotm[2].rotateZ(180.0*deg);

    //  W1 4
    W1_Presence[3] = false;
    W1_AA_CentrePosition[3] = G4ThreeVector(34.0*mm, 0.0*mm, offset_W1_BeamAxis*mm);
    W1_rotm[3].rotateX(90.0*deg);
    W1_rotm[3].rotateZ(270.0*deg);
    
    for(G4int i=0; i<numberOf_W1; i++)
    {
        if(W1_AllPresent_Override) W1_Presence[i] = true;
        if(W1_AllAbsent_Override) W1_Presence[i] = false;
        if(W1_AllPresent_Override && W1_AllAbsent_Override) W1_Presence[i] = false;
    }

    
    ////////////////////////////
    ////    HAGAR SETUP
    
    //  HAGAR NaI Crystal
    HAGAR_NaICrystal_Presence = false;
    HAGAR_NaICrystal_CentrePosition = G4ThreeVector(-100.*cm, 0.*cm, 0.*cm);
    HAGAR_rotm.rotateY(-90*deg);
    
    //  HAGAR Annulus
    HAGAR_Annulus_Presence = false;
    HAGAR_Annulus_CentrePosition = G4ThreeVector(-105.*cm, 0.*cm, 0.*cm);
    
    //  HAGAR Front
    HAGAR_FrontDisc_Presence = false;
    HAGAR_FrontDisc_CentrePosition = HAGAR_Annulus_CentrePosition + G4ThreeVector((61./2 +8/2)*cm, 0.*cm, 0.*cm);
    
    
    ////////////////////////////
    ////    CLOVER SETUP
    
    CLOVER_AllPresent_Override = false;
    CLOVER_AllAbsent_Override = false;
    
    CLOVER_Shield_AllPresent_Override = false;
    CLOVER_Shield_AllAbsent_Override = true;
    
    //--------------------------------
    useCLOVER_Walid = false;
    
    if(useCLOVER_Walid)
    {
        DefineHPGeCrystal_Walid();
        DefineHPGeCrystal_Walid_2();
    }
    
    /*
    //  CLOVER 1
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = true;
    CLOVER_Distance[0] = 10*cm;
    CLOVER_phi[0] = 90*deg;
    CLOVER_theta[0] = 135*deg;
    
    //  CLOVER 2
    CLOVER_Presence[1] = true;
    CLOVER_Shield_Presence[1] = true;
    CLOVER_Distance[1] = 10*cm;
    CLOVER_phi[1] = 0*deg;
    CLOVER_theta[1] = 135*deg;
    
    //  CLOVER 3
    CLOVER_Presence[2] = true;
    CLOVER_Shield_Presence[2] = true;
    CLOVER_Distance[2] = 10*cm;
    CLOVER_phi[2] = 270*deg;
    CLOVER_theta[2] = 135*deg;
    
    //  CLOVER 4
    CLOVER_Presence[3] = true;
    CLOVER_Shield_Presence[3] = true;
    CLOVER_Distance[3] = 10*cm;
    CLOVER_phi[3] = 180*deg;
    CLOVER_theta[3] = 135*deg;
    
    //  CLOVER 5
    CLOVER_Presence[4] = true;
    CLOVER_Shield_Presence[4] = true;
    CLOVER_Distance[4] = 10*cm;
    CLOVER_phi[4] = 45*deg;
    CLOVER_theta[4] = 90*deg;
    
    //  CLOVER 6
    CLOVER_Presence[5] = true;
    CLOVER_Shield_Presence[5] = true;
    CLOVER_Distance[5] = 10*cm;
    CLOVER_phi[5] = 0*deg;
    CLOVER_theta[5] = 90*deg;
    
    //  CLOVER 7
    CLOVER_Presence[6] = true;
    CLOVER_Shield_Presence[6] = true;
    CLOVER_Distance[6] = 10*cm;
    CLOVER_phi[6] = 180*deg;
    CLOVER_theta[6] = 90*deg;
    
    //  CLOVER 8
    CLOVER_Presence[7] = true;
    CLOVER_Shield_Presence[7] = true;
    CLOVER_Distance[7] = 10*cm;
    CLOVER_phi[7] = 135*deg;
    CLOVER_theta[7] = 90*deg;
    
    //  CLOVER 9
    CLOVER_Presence[8] = true;
    CLOVER_Shield_Presence[8] = true;
    CLOVER_Distance[8] = 10*cm;
    CLOVER_phi[8] = -45*deg;
    CLOVER_theta[8] = 90*deg;
    */
    
    //  CLOVER 1 - R1
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = false;
    CLOVER_Distance[0] = 16.6*cm;
    CLOVER_phi[0] = 180*deg;
    CLOVER_theta[0] = 90*deg;
    
    //  CLOVER 2 - R2
    CLOVER_Presence[1] = true;
    CLOVER_Shield_Presence[1] = false;
    CLOVER_Distance[1] = 15.9*cm;
    CLOVER_phi[1] = 180*deg;
    CLOVER_theta[1] = 135*deg;
    
    //  CLOVER 3 - R3
    CLOVER_Presence[2] = true;
    CLOVER_Shield_Presence[2] = false;
    CLOVER_Distance[2] = 16.6*cm;
    CLOVER_phi[2] = 135*deg;
    CLOVER_theta[2] = 125.2643897*deg;
    
    //  CLOVER 4 - R4
    CLOVER_Presence[3] = true;
    CLOVER_Shield_Presence[3] = false;
    CLOVER_Distance[3] = 16.5*cm;
    CLOVER_phi[3] = 145*deg;
    CLOVER_theta[3] = 90*deg;
    
    //  CLOVER 5 - R5
    CLOVER_Presence[4] = false;
    CLOVER_Shield_Presence[4] = false;
    CLOVER_Distance[4] = 16.3*cm;
    CLOVER_phi[4] = 216*deg;
    CLOVER_theta[4] = 72.19*deg;
    
    //  CLOVER 6 - R6
    CLOVER_Presence[5] = false;
    CLOVER_Shield_Presence[5] = false;
    CLOVER_Distance[5] = 16.6*cm;
    CLOVER_phi[5] = 288*deg;
    CLOVER_theta[5] = 72.19*deg;
    
    //  CLOVER 7 - L1
    CLOVER_Presence[6] = true;
    CLOVER_Shield_Presence[6] = false;
    CLOVER_Distance[6] = 16.3*cm;
    CLOVER_phi[6] = 45*deg;
    CLOVER_theta[6] = 125.2643897*deg;
    
    //  CLOVER 8 - L2
    CLOVER_Presence[7] = true;
    CLOVER_Shield_Presence[7] = false;
    CLOVER_Distance[7] = 16.6*cm;
    CLOVER_phi[7] = 0*deg;
    CLOVER_theta[7] = 135*deg;
    
    //  CLOVER 9 - L3
    CLOVER_Presence[8] = true;
    CLOVER_Shield_Presence[8] = false;
    CLOVER_Distance[8] = 16.7*cm;
    CLOVER_phi[8] = 0*deg;
    CLOVER_theta[8] = 90*deg;
    
    //  CLOVER 10 - L4
    CLOVER_Presence[9] = true;
    CLOVER_Shield_Presence[9] = false;
    CLOVER_Distance[9] = 16.9*cm;
    CLOVER_phi[9] = 35*deg;
    CLOVER_theta[9] = 90*deg;
    
    //  CLOVER 11 - L5
    CLOVER_Presence[10] = false;
    CLOVER_Shield_Presence[10] = false;
    CLOVER_Distance[10] = 10.6*cm;
    CLOVER_phi[10] = 324*deg;
    CLOVER_theta[10] = 107.81*deg;
    
    //  CLOVER 12 - L6
    CLOVER_Presence[11] = false;
    CLOVER_Shield_Presence[11] = false;
    CLOVER_Distance[11] = 10.6*cm;
    CLOVER_phi[11] = 36*deg;
    CLOVER_theta[11] = 135.62*deg;
    
    //  CLOVER 13
    CLOVER_Presence[12] = false;
    CLOVER_Shield_Presence[12] = false;
    CLOVER_Distance[12] = 10.6*cm;
    CLOVER_phi[12] = 108*deg;
    CLOVER_theta[12] = 135.62*deg;
    
    //  CLOVER 14
    CLOVER_Presence[13] = false;
    CLOVER_Shield_Presence[13] = false;
    CLOVER_Distance[13] = 10.6*cm;
    CLOVER_phi[13] = 180.*deg;
    CLOVER_theta[13] = 149.37*deg;
    
    //  CLOVER 15
    CLOVER_Presence[14] = false;
    CLOVER_Shield_Presence[14] = false;
    CLOVER_Distance[14] = 10.6*cm;
    CLOVER_phi[14] = 252*deg;
    CLOVER_theta[14] = 135.62*deg;
    
    //  CLOVER 16
    CLOVER_Presence[15] = false;
    CLOVER_Shield_Presence[15] = false;
    CLOVER_Distance[15] = 10.6*cm;
    CLOVER_phi[15] = 324*deg;
    CLOVER_theta[15] = 149.37*deg;
    
    //  CLOVER 17
    CLOVER_Presence[16] = false;
    CLOVER_Shield_Presence[16] = false;
    CLOVER_Distance[16] = 10.6*cm;
    CLOVER_phi[16] = 216*deg;
    CLOVER_theta[16] = 30.63*deg;

    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        double smallestPolarAngleDifference;
        std::vector<int> detNs;
        
        double x, y, z;
        
        x = sin(CLOVER_theta[i])*cos(CLOVER_phi[i]);
        y = sin(CLOVER_theta[i])*sin(CLOVER_phi[i]);
        z = cos(CLOVER_theta[i]);
        
        G4ThreeVector v1(x, y, z);
        
        bool initialComparison = false;
        
        for(G4int j=0; j<numberOf_CLOVER; j++)
        {
            double angle;
            G4ThreeVector v2;
            
            if(i!=j)
            {
                x = sin(CLOVER_theta[j])*cos(CLOVER_phi[j]);
                y = sin(CLOVER_theta[j])*sin(CLOVER_phi[j]);
                z = cos(CLOVER_theta[j]);
                
                v2 = G4ThreeVector(x, y, z);
            }
            
            angle = v1.angle(v2)/deg;

            if(!initialComparison || angle<smallestPolarAngleDifference)
            {
                detNs.clear();
                detNs.push_back(j);
                smallestPolarAngleDifference = angle;
                
                initialComparison = true;
            }
            else if(abs(angle-smallestPolarAngleDifference)<=0.1)
            {
                detNs.push_back(j);
                
                if(angle<smallestPolarAngleDifference)
                {
                    smallestPolarAngleDifference = angle;
                }
            }
        }
        
        //std::cout << "CLOVER " << i << ", smallestAngle (with detector " << detN_1 << "): " << smallestPolarAngleDifference << std::endl;
        
        std::cout << "CLOVER " << i << ", smallestAngle (with detectors";
        
        for(int j=0; j<(int) detNs.size(); j++)
        {
            std::cout << ", " << detNs[j] ;
        }
        std::cout << "): " << smallestPolarAngleDifference << std::endl;
    }
    
    /*
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        CLOVER_Presence[i] = false;
        CLOVER_Shield_Presence[i] = false;
    }
    */
    
    //----------------------------------------------------------------------------------------------------------------
    //      Used for 25.0 cm displacement (from target to face of CLOVER detector) for benchmarking against Walid
    /*
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = false;
    CLOVER_Distance[0] = (25.0-7.3)*cm;
    CLOVER_phi[0] = 0*deg;
    CLOVER_theta[0] = 0*deg;
    */
    
    
    //----------------------------------------------------------------------------------------------------------------
    //      Used for 10.6 cm displacement (from the front face of the HEAVIMET collimator to the source)
    /*
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = true;
    CLOVER_Distance[0] = 10.6*cm;
    CLOVER_phi[0] = 0*deg;
    CLOVER_theta[0] = 0*deg;
    */
    
    /*
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = true;
    CLOVER_Distance[0] = 0.0*cm;
    CLOVER_phi[0] = 0*deg;
    CLOVER_theta[0] = 0*deg;
    */
    
    
    /*
    //  Dandelion equivalent crystal distance
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = true;
    CLOVER_Distance[0] = 8.8*cm;
    CLOVER_phi[0] = 0*deg;
    CLOVER_theta[0] = 0*deg;
    */
    
    /*
    //  PK
    CLOVER_Presence[0] = true;
    CLOVER_Shield_Presence[0] = true;
    CLOVER_Distance[0] = 10.6*cm;
    CLOVER_phi[0] = 0*deg;
    CLOVER_theta[0] = 0*deg;
    */
    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        if(CLOVER_AllPresent_Override) CLOVER_Presence[i] = true;
        if(CLOVER_AllAbsent_Override) CLOVER_Presence[i] = false;
        if(CLOVER_AllPresent_Override && CLOVER_AllAbsent_Override) CLOVER_Presence[i] = false;
        
        if(CLOVER_Shield_AllPresent_Override) CLOVER_Shield_Presence[i] = true;
        if(CLOVER_Shield_AllAbsent_Override) CLOVER_Shield_Presence[i] = false;
        if(CLOVER_Shield_AllPresent_Override && CLOVER_Shield_AllAbsent_Override) CLOVER_Shield_Presence[i] = false;
    }
    
    
    ////////////////////////////////
    ////        LEPS SETUP
    
    LEPS_AllPresent_Override = false;
    LEPS_AllAbsent_Override = true;
    
    
    //  LEPS 1
    LEPS_Presence[0] = true;
    LEPS_Distance[0] = 4.5*cm;
    LEPS_phi[0] = 0.*deg;
    LEPS_theta[0] = 0.*deg;
    LEPS_rotm[0].rotateX(180.*deg);
    
    //  LEPS 2
    LEPS_Presence[1] = true;
    LEPS_Distance[1] = 4.5*cm;
    LEPS_phi[1] = 0.*deg;
    LEPS_theta[1] = 90.*deg;
    LEPS_rotm[1].rotateY(-90.*deg);
    
    //  LEPS 3
    LEPS_Presence[2] = false;
    LEPS_Distance[2] = 4.5*cm;
    LEPS_phi[2] = 90.*deg;
    LEPS_theta[2] = 90.*deg;
    LEPS_rotm[2].rotateX(90.*deg);
    
    //  LEPS 4
    LEPS_Presence[3] = false;
    LEPS_Distance[3] = 4.5*cm;
    LEPS_phi[3] = 180.*deg;
    LEPS_theta[3] = 90*deg;
    LEPS_rotm[3].rotateY(90.*deg);
    
    //  LEPS 5
    LEPS_Presence[4] = false;
    LEPS_Distance[4] = 4.5*cm;
    LEPS_phi[4] = 270.*deg;
    LEPS_theta[4] = 90*deg;
    LEPS_rotm[4].rotateX(-90.*deg);
    
    
    //  LEPS 6
    LEPS_Presence[5] = false;
    LEPS_Distance[5] = 4.5*cm;
    LEPS_phi[5] = 0*deg;
    LEPS_theta[5] = 180*deg;
    LEPS_rotm[5].rotateY(0.*deg);
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        if( LEPS_AllPresent_Override == true ) LEPS_Presence[i] = true;
        if( LEPS_AllAbsent_Override == true ) LEPS_Presence[i] = false;
        if( LEPS_AllPresent_Override == true && LEPS_AllAbsent_Override == true ) LEPS_Presence[i] = false;
    }

    
    ////////////////////////////////
    ////        LaBr3Ce SETUP
    
    //------------------------------------------------
    //      MANUAL SETUP: Unneccesary if SetupPreconfiguredVersion() is used
    
    
    LaBr3Ce_AllPresent_Override = false;
    LaBr3Ce_AllAbsent_Override = false;
    
    /*
    LaBr3CeSetupVersion = 1;
    LaBR3Ce_SetGlobalDistance = false;
    LaBR3Ce_automaticOrientation = false;
    configuration_truncatedIcosahedron_hexagons = false;
    
    G4double LaBR3Ce_GlobalDistance = 0.0;
    
    //      Full cylindrical crystal
    LaBR3Ce_GlobalDistance = 13.1*cm;
    
    //  LaBr3Ce 1
    LaBr3Ce_Presence[0] = true;
    LaBr3Ce_Distance[0] = 0.0*cm;
    */
    
    
    //  LaBr3Ce 1
    LaBr3Ce_Presence[0] = true;
    LaBr3Ce_Distance[0] = 281.0*mm;
    LaBr3Ce_theta[0] = 315*deg;
    LaBr3Ce_phi[0] = 0.0*deg;
    
    //  LaBr3Ce 2
    LaBr3Ce_Presence[1] = true;
    LaBr3Ce_Distance[1] = 243.0*mm;
    LaBr3Ce_theta[1] = 270.0*deg;
    LaBr3Ce_phi[1] = 35.0*deg;
    
    //  LaBr3Ce 3
    LaBr3Ce_Presence[2] = false;
    LaBr3Ce_Distance[2] = 217.8*mm;
    LaBr3Ce_theta[2] = 135.0*deg;
    LaBr3Ce_phi[2] = 0.0*deg;
    
    //  LaBr3Ce 4
    LaBr3Ce_Presence[3] = false;
    LaBr3Ce_Distance[3] = 236.8*mm;
    LaBr3Ce_theta[3] = 90.0*deg;
    LaBr3Ce_phi[3] = 62.0*deg;
    
    
    //------------------------------------------------
    //      PRECONFIGURED SETUPS
    
    //      Work for Christiaan/Mathis/Katarzyna
    //SetupPreconfiguredVersion(-1);

    //  No tapering
    //  LaBR3Ce_GlobalDistance = 13.2*cm;
    //SetupPreconfiguredVersion(0);

    //  LaBR3Ce_GlobalDistance = 12.0*cm;
    //SetupPreconfiguredVersion(1);
    
    //  LaBR3Ce_GlobalDistance = 11.0*cm;
    //SetupPreconfiguredVersion(2);

    //  LaBR3Ce_GlobalDistance = 10.0*cm;
    //SetupPreconfiguredVersion(3);

    //  LaBR3Ce_GlobalDistance = 9.0*cm;
    //SetupPreconfiguredVersion(4);

    //  LaBR3Ce_GlobalDistance = 8.0*cm;
    //SetupPreconfiguredVersion(5);

    //  LaBR3Ce_GlobalDistance = 7.0*cm;
    //SetupPreconfiguredVersion(6);

    //  LaBR3Ce_GlobalDistance = 6.0*cm;
    //SetupPreconfiguredVersion(7);

    //  LaBR3Ce_GlobalDistance = 8.9632*cm;
    //SetupPreconfiguredVersion(8);
    
    //  LaBR3Ce_GlobalDistance = 10.0*cm;
    //SetupPreconfiguredVersion(9);

    //  LaBR3Ce_GlobalDistance = 12.0*cm;
    //SetupPreconfiguredVersion(10);

    //------------------------------------------------
    
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        if( LaBr3Ce_AllPresent_Override == true ) LaBr3Ce_Presence[i] = true;
        if( LaBr3Ce_AllAbsent_Override == true ) LaBr3Ce_Presence[i] = false;
        if( LaBr3Ce_AllPresent_Override == true && LaBr3Ce_AllAbsent_Override == true ) LaBr3Ce_Presence[i] = false;
        
        if(LaBR3Ce_SetGlobalDistance)
        {
            LaBr3Ce_Distance[i] = LaBR3Ce_GlobalDistance;
        }
        
        double x, y, z;
        double theta, phi, r;
        
        if(i<20)
        {
            x = vertex_hexagonFaces_truncatedIcosahedron[i].x();
            y = vertex_hexagonFaces_truncatedIcosahedron[i].y();
            z = vertex_hexagonFaces_truncatedIcosahedron[i].z();
            
            r = vertex_hexagonFaces_truncatedIcosahedron[i].mag();
        }
        else if(i<32)
        {
            int pentIndex = i-20;
            x = vertex_pentagonFaces_truncatedIcosahedron[pentIndex].x();
            y = vertex_pentagonFaces_truncatedIcosahedron[pentIndex].y();
            z = vertex_pentagonFaces_truncatedIcosahedron[pentIndex].z();
            
            r = vertex_pentagonFaces_truncatedIcosahedron[pentIndex].mag();
        }
        
        theta = acos(z/r)/deg;
        
        if(x>0.0)
        {
            phi = atan(y/x)/deg;
        }
        else if(x<0.0 && y>=0.0)
        {
            phi = atan(y/x)/deg + 180.0;
        }
        else if(x<0.0 && y<0.0)
        {
            phi = atan(y/x)/deg - 180.0;
        }
        else if(x==0.0 && y>0.0)
        {
            phi = 90.0;
        }
        else if(x==0.0 && y<0.0)
        {
            phi = -90.0;
        }
        else if(x==0.0 && y==0.0)
        {
            phi = G4RandFlat::shoot(0.0, 360.0);
        }
        
        if(configuration_truncatedIcosahedron_hexagons)
        {
            LaBr3Ce_theta[i] = theta*deg;
            LaBr3Ce_phi[i] = phi*deg;
        }
    }

    
    ////////////////////////////////////////////
    ////                                    ////
    ////        K600 SPECTROMETER SETUP     ////
    ////                                    ////
    ////////////////////////////////////////////
    
    //  K600 Quadrapole
    K600_Quadrupole = false;
    Ideal_Quadrupole = false;
    Mapped_Quadrupole = true;
    //K600_Q_gradient = 0.030*tesla/cm;  // gradient = dB/dr for Ideal Quadrupole
    K600_Q_gradient = 0.001*tesla/cm;  // gradient = dB/dr for Ideal Quadrupole
    K600_Quadrupole_CentrePosition = G4ThreeVector(0.*cm, 0.*cm,100.*cm);
    //K600_Quadrupole_rotm.rotateZ(90.*deg);
    //K600_Quadrupole_rotm.rotateY(90*deg);
    //K600_Quadrupole_rotm2.rotateX(90*deg);
    
    //  K600 Dipole 1
    K600_Dipole1 = false;
    K600_Dipole1_BZ = -2.30*tesla;
    //K600_Dipole1_BZ = -3.30*tesla;
    K600_Dipole1_CentrePosition = G4ThreeVector(75*cm, 0.*cm,250*cm);
    K600_Dipole1_rotm.rotateX(-90.*deg);
    K600_Dipole1_rotm.rotateY(180.*deg);
    
    //  K600 Dipole 2
    K600_Dipole2 = false;
    K600_Dipole2_BZ = -2.40*tesla;
    //K600_Dipole2_BZ = -3.40*tesla;
    K600_Dipole2_CentrePosition = G4ThreeVector(75*cm, 0.*cm,250*cm);
    K600_Dipole2_rotm.rotateX(-90*deg);
    K600_Dipole2_rotm.rotateY(180.*deg);
    
    if(Ideal_Quadrupole && Mapped_Quadrupole || (!Ideal_Quadrupole && !Mapped_Quadrupole))
    {
        Ideal_Quadrupole = false;
        Mapped_Quadrupole = false;
    }
    
    ////////////////////////////////////////////////////
    ////                                            ////
    ////        K600 VAULT - STRUCTURE SETUP        ////
    ////                                            ////
    ////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, both sides on
    K600_BACTAR_sidesOn_Presence = true;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, both sides off
    K600_BACTAR_sidesOff_Presence = false;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, beam left side side off
    K600_BACTAR_beamRightSideOff_Presence = false;
    
    ////////////////////////////////////////////////
    ////    New K600 Target Chamber - New scattering chamber, beam right side side off
    K600_BACTAR_beamLeftSideOff_Presence = false;

    ////////////////////////////////////////////////
    ////    ALBA shielding - Truncated Icosahedron
    K600_ALBA_TruncIcos_Shielding_Presence = false;
    
    /////////////////////////////////////
    //  K600 Target
    K600_Target_Presence = false;
    
    /////////////////////////////////////
    //  K600 Target Backing
    K600_TargetBacking_Presence = false;
    
    ////////////////////////////////////////////////
    ////    New AFRODITE Target Chamber by Mathis
    AFRODITE_MathisTC_Presence = false;

    
    // Define materials
    DefineMaterials();
    
    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    //  NIST Material Database - Materials
    nistManager->FindOrBuildMaterial("G4_Galactic");
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_BGO");
    nistManager->FindOrBuildMaterial("G4_Ge");
    nistManager->FindOrBuildMaterial("G4_Al");
    nistManager->FindOrBuildMaterial("G4_Si");
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    nistManager->FindOrBuildMaterial("G4_MYLAR");
    nistManager->FindOrBuildMaterial("G4_W");
    nistManager->FindOrBuildMaterial("G4_Ar");
    nistManager->FindOrBuildMaterial("G4_Pb");
    nistManager->FindOrBuildMaterial("G4_C");
    nistManager->FindOrBuildMaterial("G4_Li");
    nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    nistManager->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    nistManager->FindOrBuildMaterial("G4_LITHIUM_CARBONATE");
    nistManager->FindOrBuildMaterial("G4_Be");

    //  NIST Elementary Material Database - ELEMENTS
    nistManager->FindOrBuildElement("H");
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("N");
    nistManager->FindOrBuildElement("O");
    nistManager->FindOrBuildElement("Fe");
    nistManager->FindOrBuildElement("Co");
    nistManager->FindOrBuildElement("Ni");
    nistManager->FindOrBuildElement("Cu");
    nistManager->FindOrBuildElement("Pb");
    nistManager->FindOrBuildElement("W");
    nistManager->FindOrBuildElement("Li");
    nistManager->FindOrBuildElement("Ar");
    nistManager->FindOrBuildElement("La");
    nistManager->FindOrBuildElement("Br");
    nistManager->FindOrBuildElement("Ce");
    
    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    //////////////////////////////////////
    //          Get Elements            //
    //////////////////////////////////////
    
    G4Element* Hydrogen = G4Element::GetElement("H");
    G4Element* Carbon = G4Element::GetElement("C");
    G4Element* Nitrogen = G4Element::GetElement("N");
    G4Element* Oxygen = G4Element::GetElement("O");
    G4Element* Iron = G4Element::GetElement("Fe");
    G4Element* Cobalt = G4Element::GetElement("Co");
    G4Element* Nickel = G4Element::GetElement("Ni");
    G4Element* Copper = G4Element::GetElement("Cu");
    G4Element* Lead = G4Element::GetElement("Pb");
    G4Element* Tungsten = G4Element::GetElement("W");
    G4Element* Lithium = G4Element::GetElement("Li");
    G4Element* Argon = G4Element::GetElement("Ar");
    G4Element* Lanthanum = G4Element::GetElement("La");
    G4Element* Bromine = G4Element::GetElement("Br");
    G4Element* Cerium = G4Element::GetElement("Ce");
    
    
    //////////////////////////////////////
    //          Get Materials           //
    //////////////////////////////////////
    
    ////    NIST Defined Elemental Material
    G4Material* G4_Ge_Material = G4Material::GetMaterial("G4_Ge");
    G4Material* G4_Al_Material  = G4Material::GetMaterial("G4_Al");
    G4Material* G4_Si_Material = G4Material::GetMaterial("G4_Si");
    G4Material* G4_W_Material = G4Material::GetMaterial("G4_W");
    G4Material* G4_Ar_Material = G4Material::GetMaterial("G4_Ar");
    G4Material* G4_C_Material = G4Material::GetMaterial("G4_C");
    G4Material* G4_Li_Material = G4Material::GetMaterial("G4_Li");
    G4Material* G4_Pb_Material = G4Material::GetMaterial("G4_Pb");
    
    ////    NIST Defined Materials and Compounds
    G4Material* G4_Galactic_Material = G4Material::GetMaterial("G4_Galactic");
    G4Material* G4_AIR_Material = G4Material::GetMaterial("G4_AIR");
    G4Material* G4_BGO_Material = G4Material::GetMaterial("G4_BGO");
    G4Material* G4_PLASTIC_SC_VINYLTOLUENE_Material = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* G4_Mylar_Material = G4Material::GetMaterial("G4_MYLAR");
    G4Material* G4_CARBON_DIOXIDE_Material = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
    G4Material* G4_SODIUM_IODIDE_Material = G4Material::GetMaterial("G4_SODIUM_IODIDE");
    G4Material* G4_LITHIUM_CARBONATE_Material = G4Material::GetMaterial("G4_LITHIUM_CARBONATE");
    G4Material* G4_Be_Material = G4Material::GetMaterial("G4_Be");

    ////    CLOVER Detector Shield, HEAVIMET Material
    G4Material* Heavimet_Material = new G4Material("Heavimet_Material",19.25*g/cm3, 5);
    Heavimet_Material->AddElement(Tungsten, 94.20*perCent);
    Heavimet_Material->AddElement(Nickel, 4.35*perCent);
    Heavimet_Material->AddElement(Iron, 0.85*perCent);
    Heavimet_Material->AddElement(Cobalt, 0.50*perCent);
    Heavimet_Material->AddElement(Copper, 0.10*perCent);
    
    ////    LaBr3Ce3Ce Detector, LaBr3Ce3C Material
    G4Material* LaBr3Ce3Ce_Material = new G4Material("LaBr3Ce3Ce_Material",5.08*g/cm3, 2);
    //      Version 1: LaBr3Ce3Ce
    //LaBr3Ce3Ce_Material->AddElement(Lanthanum, 26.77781134886090797*perCent);
    //LaBr3Ce3Ce_Material->AddElement(Bromine, 46.211014685441440975*perCent);
    //LaBr3Ce3Ce_Material->AddElement(Cerium, 27.011173965697647503*perCent);
    //      Version 2: LaBr3Ce3 (no Cerium doping)
    LaBr3Ce3Ce_Material->AddElement(Lanthanum, 36.687549045214417731*perCent);
    LaBr3Ce3Ce_Material->AddElement(Bromine, 63.312450954785568058*perCent);
    
    ////    CAKE Detector Materials
    //G4Material*         CAKE_PCBMaterial;
    
    ////    PADDLE Detector Materials
    G4Material* BC408_Material = new G4Material("BC408_Material",1.032*g/cm3, 2);
    BC408_Material->AddElement(Hydrogen, 8.4748*perCent);
    BC408_Material->AddElement(Carbon, 91.5252*perCent);
    
    ////    VDC Detector Materials
    G4Material* VDC_SR_Gas_Material = new G4Material("VDC_SR_Gas_Material", 0.00184212*g/cm3, 2);
    VDC_SR_Gas_Material->AddMaterial(G4_Ar_Material, 90.*perCent);
    VDC_SR_Gas_Material->AddMaterial(G4_CARBON_DIOXIDE_Material, 10.*perCent);
    
    
    ////    Vacuum - Target Chamber
    G4Material* Vacuum_TC_Material = new G4Material("Vacuum_TC", 1.2e-14*g/cm3, 4);
    Vacuum_TC_Material->AddElement(Carbon, 0.0124*perCent);
    Vacuum_TC_Material->AddElement(Nitrogen, 75.5268*perCent);
    Vacuum_TC_Material->AddElement(Oxygen, 23.1781*perCent);
    Vacuum_TC_Material->AddElement(Argon, 1.2827*perCent);
    
    
    ////    Target
    //G4Isotope* Be9 = new G4Isotope(name="Be9", iz=5, n=9, a=235.01*g/mole);

    
    /*
     4     G4_AIR                     0.00120479        85.7
     6     0.000124
     7     0.755268
     8     0.231781
     18    0.012827
     */
    
    
    ////    Variables for CADMesh
    char    meshPath[512];
    char    meshType[] = "PLY";

    
    //////////////////////////////////////////////////////////
    //                      WORLD                           //
    //////////////////////////////////////////////////////////
    
    G4Box* SolidWorld = new G4Box("World", WorldSize/2, WorldSize/2, WorldSize/2);
    
    G4LogicalVolume*
    LogicWorld = new G4LogicalVolume(SolidWorld,                //its solid
                                     G4_Galactic_Material,      //its material
                                     "World");                  //its name
    G4VPhysicalVolume*
    PhysiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),          //at (0,0,0)
                                   LogicWorld,               //its logical volume
                                   "World",                  //its name
                                   0,                        //its mother  volume
                                   false,                    //no boolean operation
                                   0);                       //copy number
    
    
    
    //////////////////////////////////////////////////////////
    //                      VACUUM CHAMBER
    //////////////////////////////////////////////////////////
    
    G4ThreeVector positionVacuumChamber = G4ThreeVector(0,0,0);
    
    G4Box* SolidVacuumChamber = new G4Box("VacuumChamber", (200./2)*cm, (200./2)*cm, (200./2)*cm);
    
    G4LogicalVolume* LogicVacuumChamber = new G4LogicalVolume(SolidVacuumChamber, G4_Galactic_Material,"VacuumChamber",0,0,0);
    
    
    new G4PVPlacement(0,               // no rotation
                      positionVacuumChamber, // at (x,y,z)
                      LogicVacuumChamber,       // its logical volume
                      "VacuumChamber",       // its name
                      LogicWorld,         // its mother  volume
                      false,           // no boolean operations
                      0,               // copy number
                      fCheckOverlaps); // checking overlaps
    
    
    
    
    
    //////////////////////////////////////////////////////////
    //              Scattering Chamber - CADMesh
    //////////////////////////////////////////////////////////
    
    if(K600_BACTAR_sidesOn_Presence)
    {
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/STRUCTURES/BACTAR/BACTAR_sidesOn.ply");
    }
    if(K600_BACTAR_sidesOff_Presence)
    {
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/STRUCTURES/BACTAR/BACTAR_sidesOff.ply");
    }
    if(K600_BACTAR_beamRightSideOff_Presence)
    {
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/STRUCTURES/BACTAR/BACTAR_beamRightSideOff.ply");
    }
    if(K600_BACTAR_beamLeftSideOff_Presence)
    {
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/STRUCTURES/BACTAR/BACTAR_beamLightSideOff.ply");
    }

    
    if(K600_BACTAR_sidesOn_Presence || K600_BACTAR_sidesOff_Presence || K600_BACTAR_beamRightSideOff_Presence || K600_BACTAR_beamLeftSideOff_Presence)
    {
        G4ThreeVector offset_BACTAR = G4ThreeVector(0*cm, 0*cm, 0*cm);
        
        CADMesh * mesh_BACTAR = new CADMesh(meshPath, meshType, mm, offset_BACTAR, false);
        
        G4VSolid * SolidBACTAR = mesh_BACTAR->TessellatedMesh();
        
        G4LogicalVolume* LogicBACTAR = new G4LogicalVolume(SolidBACTAR, G4_Al_Material, "BACTAR", 0, 0, 0);
        
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (x,y,z)
                          LogicBACTAR,       // its logical volume
                          "BACTAR",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
        //  Visualisation
        G4VisAttributes* K600_BACTAR_TC_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
        K600_BACTAR_TC_VisAtt->SetForceSolid(true);
        LogicBACTAR ->SetVisAttributes(K600_BACTAR_TC_VisAtt);
    }
    
    //////////////////////////////////////////////////////////
    //              Scattering Chamber - CADMesh
    //////////////////////////////////////////////////////////
    
    if(K600_ALBA_TruncIcos_Shielding_Presence)
    {
        //  Without tolerances
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/STRUCTURES/ALBA/ALBA_shield.ply");
        //  With tolerances
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/STRUCTURES/ALBA/ALBA_shield_mod.ply");
    }
    
    if(K600_ALBA_TruncIcos_Shielding_Presence)
    {
        G4ThreeVector offset_ALBA_shield = G4ThreeVector(0*cm, 0*cm, 0*cm);
        
        CADMesh * mesh_ALBA_shield = new CADMesh(meshPath, meshType, mm, offset_ALBA_shield, false);
        
        G4VSolid * SolidALBA_shield = mesh_ALBA_shield->TessellatedMesh();
        
        G4LogicalVolume* LogicALBA_shield = new G4LogicalVolume(SolidALBA_shield, G4_Pb_Material, "ALBA_shield", 0, 0, 0);
        
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (x,y,z)
                          LogicALBA_shield,       // its logical volume
                          "ALBA_shield",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
        //  Visualisation
        G4VisAttributes* K600_ALBA_shield_TruncIcos_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
        K600_ALBA_shield_TruncIcos_VisAtt->SetForceSolid(true);
        LogicALBA_shield ->SetVisAttributes(K600_ALBA_shield_TruncIcos_VisAtt);
    }
    
    
    //////////////////////////////////////////////////////////
    //              Scattering Chamber - CADMesh
    //////////////////////////////////////////////////////////
    
    if(AFRODITE_MathisTC_Presence)
    {
        G4ThreeVector offset_MathisTC = G4ThreeVector(0*cm, 0*cm, 0*cm);
        
        CADMesh * mesh_MathisTC = new CADMesh("../K600-ALBA/Mesh-Models/STRUCTURES/MathisTC/MathisTC.ply", "PLY", mm, offset_MathisTC, false);
        
        G4VSolid * SolidMathisTC = mesh_MathisTC->TessellatedMesh();
        
        G4LogicalVolume* LogicMathisTC = new G4LogicalVolume(SolidMathisTC, G4_Al_Material, "BACTAR", 0, 0, 0);
        
        /*
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(), // at (x,y,z)
                          LogicMathisTC,       // its logical volume
                          "BACTAR",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        */
        
        //------------------------------------------------
        //      Work for Christiaan/Mathis/Katarzyna

        G4RotationMatrix    AFRODITE_MathisTC_rotm;
        G4ThreeVector       AFRODITE_MathisTC_position = G4ThreeVector();

        AFRODITE_MathisTC_rotm.rotateZ(-90.0*deg);
        
        G4Transform3D AFRODITE_MathisTC_transform = G4Transform3D(AFRODITE_MathisTC_rotm, AFRODITE_MathisTC_position);

        new G4PVPlacement(AFRODITE_MathisTC_transform,
                          LogicMathisTC,       // its logical volume
                          "BACTAR",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
        //  Visualisation
        G4VisAttributes* AFRODITE_MathisTC_VisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
        AFRODITE_MathisTC_VisAtt->SetForceSolid(true);
        LogicMathisTC->SetVisAttributes(AFRODITE_MathisTC_VisAtt);
    }

    
    //////////////////////////////////////////////////////////
    //                  TARGET DEFINITION
    //////////////////////////////////////////////////////////
    
    G4double targetThickness = 2.42; // um
    //G4double targetThickness = 2.18; // um
    //G4double targetThickness = 3.5; // um
    //G4double targetThickness = 4.5; // um
    
    
    ////    For Daniel Lombarri's Experiment
    //G4double targetThickness = 1.5; // um
    //G4double targetThickness = 100.; // um
    
    G4ThreeVector positionTarget = G4ThreeVector(0,0,0);
    
    G4Box* SolidTarget = new G4Box("Target", (20./2)*mm, (20./2)*mm, (targetThickness/2)*um);
    
    
    //G4LogicalVolume* LogicTarget = new G4LogicalVolume(SolidTarget, G4_LITHIUM_CARBONATE_Material,"Target",0,0,0);
    G4LogicalVolume* LogicTarget = new G4LogicalVolume(SolidTarget, G4_Be_Material,"Target",0,0,0);

    if(K600_Target_Presence)
    {
        new G4PVPlacement(0,               // no rotation
                          positionTarget, // at (x,y,z)
                          LogicTarget,       // its logical volume
                          "Target",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
    }
    
    
    //////////////////////////////////////////////////////
    //                  TARGET BACKING
    //////////////////////////////////////////////////////
    
    G4double targetBackingThickness = 0.25; // um
    
    G4ThreeVector positionTargetBacking = G4ThreeVector(0,0,-((targetThickness/2) + (targetBackingThickness/2))*um);
    
    G4Box* SolidTargetBacking = new G4Box("TargetBacking", (20./2)*mm, (20./2)*mm, (targetBackingThickness/2)*um);
    
    G4LogicalVolume* LogicTargetBacking = new G4LogicalVolume(SolidTargetBacking, G4_C_Material,"TargetBacking",0,0,0);
    
    if(K600_TargetBacking_Presence)
    {
        new G4PVPlacement(0,               // no rotation
                          positionTargetBacking, // at (x,y,z)
                          LogicTargetBacking,       // its logical volume
                          "TargetBacking",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
    }
    
    
    
    
    /////////////////////////////////////////////
    //             CAKE DEFINITION
    /////////////////////////////////////////////
    
    
    /////////////////////////////////////////
    //          CAKE - PCB
    
    G4double thickness_CAKE_PCB = 1.6e+3; // um
    
    G4Tubs* Solid_CAKE_PCB_filled = new G4Tubs("CAKE_PCB",0.0*mm, 170.0*mm, (thickness_CAKE_PCB/2)*um, 0*deg, (57.74)*deg);
    
    G4Box* PCB_razor1 = new G4Box("PCB_razor1", (200/2)*mm, (20.)*mm, (100/2)*mm);
    
    G4ThreeVector   position_CAKE_PCB_razor1[4];
    position_CAKE_PCB_razor1[0] = G4ThreeVector(7.0236770*mm, 2.0832170*mm, 0.0*mm);
    position_CAKE_PCB_razor1[1] = G4ThreeVector(171.1783094*mm, 25.6019696*mm, 0.0*mm);
    position_CAKE_PCB_razor1[2] = G4ThreeVector(154.2410751*mm, 99.9806850*mm,0.0*mm);
    position_CAKE_PCB_razor1[3] = G4ThreeVector(107.0949198*mm, 135.0207195*mm,0.0*mm);
    
    G4RotationMatrix* rm1_CAKE_PCB = new G4RotationMatrix();
    G4RotationMatrix* rm2_CAKE_PCB = new G4RotationMatrix();
    G4RotationMatrix* rm3_CAKE_PCB = new G4RotationMatrix();
    G4RotationMatrix* rm4_CAKE_PCB = new G4RotationMatrix();
    
    rm1_CAKE_PCB->rotateZ(59.63*deg);
    rm2_CAKE_PCB->rotateZ(101.2030637*deg);
    rm3_CAKE_PCB->rotateZ(59.63*deg);
    rm4_CAKE_PCB->rotateZ(15.8652730*deg);
    
    G4VSolid* Solid_CAKE_PCB_shear1 = new G4SubtractionSolid("CAKE_PCB_shear1", Solid_CAKE_PCB_filled, PCB_razor1, rm1_CAKE_PCB, position_CAKE_PCB_razor1[0]);
    G4VSolid* Solid_CAKE_PCB_shear2 = new G4SubtractionSolid("CAKE_PCB_shear2", Solid_CAKE_PCB_shear1, PCB_razor1, rm2_CAKE_PCB, position_CAKE_PCB_razor1[1]);
    G4VSolid* Solid_CAKE_PCB_shear3 = new G4SubtractionSolid("CAKE_PCB_shear3", Solid_CAKE_PCB_shear2, PCB_razor1, rm3_CAKE_PCB, position_CAKE_PCB_razor1[2]);
    G4VSolid* Solid_CAKE_PCB_shear4 = new G4SubtractionSolid("CAKE_PCB_shear4", Solid_CAKE_PCB_shear3, PCB_razor1, rm4_CAKE_PCB, position_CAKE_PCB_razor1[3]);
    
    /////////////////////////////////////////
    //          CAKE - PCB PUNCH
    
    G4double thickness_CAKE_PCBpunch = 2.0e+3; // um
    
    G4Tubs* Solid_CAKE_PCBpunch_filled = new G4Tubs("CAKE_PCBpunch_filled",0.0*mm, 160.0*mm, (thickness_CAKE_PCBpunch/2)*um, 0.*deg, (60.)*deg);
    
    G4Box* PCBpunch_shear1 = new G4Box("PCBpunch_shear1", (300/2)*mm, (20.)*mm, (100/2)*mm);
    
    G4ThreeVector   position_CAKE_PCBpunch_shear[6];
    position_CAKE_PCBpunch_shear[0] = G4ThreeVector(7.6909778*mm, 5.8892601*mm, 0.0*mm);
    position_CAKE_PCBpunch_shear[1] = G4ThreeVector(82.8083631*mm, -15.9622884*mm, 0.0*mm);
    position_CAKE_PCBpunch_shear[2] = G4ThreeVector(152.9061303*mm, 34.2318836*mm, 0.0*mm);
    position_CAKE_PCBpunch_shear[3] = G4ThreeVector(133.6409271*mm, 81.5333641*mm, 0.0*mm);
    position_CAKE_PCBpunch_shear[4] = G4ThreeVector(105.9034003*mm, 114.8979515*mm, 0.0*mm);
    position_CAKE_PCBpunch_shear[5] = G4ThreeVector(26.8609291*mm, 78.1519180*mm, 0.0*mm);
    
    G4RotationMatrix* rm_PCBpunch_shear1 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear2 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear3 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear4 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear5 = new G4RotationMatrix();
    G4RotationMatrix* rm_PCBpunch_shear6 = new G4RotationMatrix();
    
    rm_PCBpunch_shear1->rotateZ(59.6297466*deg);
    rm_PCBpunch_shear2->rotateZ(-3.0000000*deg);
    rm_PCBpunch_shear3->rotateZ(76.7362500*deg);
    rm_PCBpunch_shear4->rotateZ(59.6300000*deg);
    rm_PCBpunch_shear5->rotateZ(42.5237500*deg);
    rm_PCBpunch_shear6->rotateZ(-57.7400000*deg);
    
    G4VSolid* Solid_CAKE_PCB_punch1 = new G4SubtractionSolid("CAKE_PCB_punch1", Solid_CAKE_PCBpunch_filled, PCBpunch_shear1, rm_PCBpunch_shear1, position_CAKE_PCBpunch_shear[0]);
    G4VSolid* Solid_CAKE_PCB_punch2 = new G4SubtractionSolid("CAKE_PCB_punch2", Solid_CAKE_PCB_punch1, PCBpunch_shear1, rm_PCBpunch_shear2, position_CAKE_PCBpunch_shear[1]);
    G4VSolid* Solid_CAKE_PCB_punch3 = new G4SubtractionSolid("CAKE_PCB_punch3", Solid_CAKE_PCB_punch2, PCBpunch_shear1, rm_PCBpunch_shear3, position_CAKE_PCBpunch_shear[2]);
    G4VSolid* Solid_CAKE_PCB_punch4 = new G4SubtractionSolid("CAKE_PCB_punch4", Solid_CAKE_PCB_punch3, PCBpunch_shear1, rm_PCBpunch_shear4, position_CAKE_PCBpunch_shear[3]);
    G4VSolid* Solid_CAKE_PCB_punch5 = new G4SubtractionSolid("CAKE_PCB_punch5", Solid_CAKE_PCB_punch4, PCBpunch_shear1, rm_PCBpunch_shear5, position_CAKE_PCBpunch_shear[4]);
    G4VSolid* Solid_CAKE_PCB_punch6 = new G4SubtractionSolid("CAKE_PCB_punch6", Solid_CAKE_PCB_punch5, PCBpunch_shear1, rm_PCBpunch_shear6, position_CAKE_PCBpunch_shear[5]);
    
    ////    This offset is the offset between the origins of the created PCB and the PCBpunch
    //G4ThreeVector position_CAKE_PCB_perf = G4ThreeVector(1.0437151*mm, -0.2199451*mm, 0.0*mm);
    G4ThreeVector position_CAKE_PCB_perf = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    
    G4RotationMatrix* rm_perf = new G4RotationMatrix();
    //rm_perf->rotateZ(-3.*deg);
    
    G4VSolid* Solid_CAKE_PCB_punched = new G4SubtractionSolid("CAKE_PCB_punched", Solid_CAKE_PCB_shear4, Solid_CAKE_PCB_punch6, rm_perf, position_CAKE_PCB_perf);
    
    
    /////////////////////////////////////////
    //          CAKE - Silicon Wafer
    
    G4double thickness_SiliconWafer_MMM = 400.0; // um
    
    G4double thickness_DL_MMM_Front = 2.0; // um
    G4double thickness_DL_MMM_Back = 2.0; // um
    
    G4double thickness_2M_MMM_Front = 0.3; // um
    G4double thickness_2M_MMM_Back = 0.3; // um
    
    ////    Active Area (Silicon) thickness
    G4double thickness_MMM_AA = thickness_SiliconWafer_MMM - thickness_DL_MMM_Front - thickness_DL_MMM_Back; // um
    
    G4Tubs* Solid_CAKE_SiliconWafer_filled = new G4Tubs("CAKE_SiliconWafer_filled", 0.0*mm, 160.0*mm, (thickness_MMM_AA/2)*um, 0.*deg, (60.)*deg);
    
    G4Box* SiliconWafer_shear1 = new G4Box("SiliconWafer_shear1", (300/2)*mm, (20.)*mm, (100/2)*mm);
    
    G4ThreeVector   position_CAKE_PCB_shear[6];
    position_CAKE_PCB_shear[0] = G4ThreeVector(7.5760469*mm, 4.1074720*mm, 0.0*mm);
    position_CAKE_PCB_shear[1] = G4ThreeVector(86.1695456*mm, -16.7875086*mm, 0.0*mm);
    position_CAKE_PCB_shear[2] = G4ThreeVector(152.5181096*mm, 40.2365300*mm, 0.0*mm);
    position_CAKE_PCB_shear[3] = G4ThreeVector(135.4631769*mm, 80.4016034*mm, 0.0*mm);
    position_CAKE_PCB_shear[4] = G4ThreeVector(105.7796837*mm, 116.3682685*mm, 0.0*mm);
    position_CAKE_PCB_shear[5] = G4ThreeVector(28.3491241*mm, 82.3831462*mm, 0.0*mm);
    
    G4RotationMatrix* rm_SiliconWafer_shear1 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear2 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear3 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear4 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear5 = new G4RotationMatrix();
    G4RotationMatrix* rm_SiliconWafer_shear6 = new G4RotationMatrix();
    
    rm_SiliconWafer_shear1->rotateZ(59.6300000*deg);
    rm_SiliconWafer_shear2->rotateZ(-3.0000000*deg);
    rm_SiliconWafer_shear3->rotateZ(76.7362500*deg);
    rm_SiliconWafer_shear4->rotateZ(59.6300000*deg);
    rm_SiliconWafer_shear5->rotateZ(42.5237500*deg);
    rm_SiliconWafer_shear6->rotateZ(-57.7400000*deg);
    
    G4VSolid* Solid_CAKE_SiliconWafer_shear1 = new G4SubtractionSolid("CAKE_SiliconWafer_shear1", Solid_CAKE_SiliconWafer_filled, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_CAKE_PCB_shear[0]);
    G4VSolid* Solid_CAKE_SiliconWafer_shear2 = new G4SubtractionSolid("CAKE_SiliconWafer_shear2", Solid_CAKE_SiliconWafer_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_CAKE_PCB_shear[1]);
    G4VSolid* Solid_CAKE_SiliconWafer_shear3 = new G4SubtractionSolid("CAKE_SiliconWafer_shear3", Solid_CAKE_SiliconWafer_shear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_CAKE_PCB_shear[2]);
    G4VSolid* Solid_CAKE_SiliconWafer_shear4 = new G4SubtractionSolid("CAKE_SiliconWafer_shear3", Solid_CAKE_SiliconWafer_shear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_CAKE_PCB_shear[3]);
    G4VSolid* Solid_CAKE_SiliconWafer_shear5 = new G4SubtractionSolid("CAKE_SiliconWafer_shear3", Solid_CAKE_SiliconWafer_shear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_CAKE_PCB_shear[4]);
    G4VSolid* Solid_CAKE_SiliconWafer_Final = new G4SubtractionSolid("CAKE_SiliconWafer_Final", Solid_CAKE_SiliconWafer_shear5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_CAKE_PCB_shear[5]);
    
    ////    This offset is the offset between the origins of the created PCB and the PCBpunch
    G4ThreeVector position_CAKE_SiliconWafer_shear = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    
    G4RotationMatrix* rm_SiliconWafer_shear = new G4RotationMatrix();
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////    Final PCB
    ////    This includes the final PCB shear, which makes space for the Silicon Wafer
    
    // The "+10" for the thickness is to prevent logical issues when using boolean subtraction
    //G4Tubs* Solid_CAKE_SiliconWafer_finalshear_filled = new G4Tubs("CAKE_SiliconWafer_finalshear_filled",0.0*mm, 160.0*mm, ((thickness_MMM_AA/2) + 10)*um, 0.*deg, 60.*deg);
    G4Tubs* Solid_CAKE_SiliconWafer_finalshear_filled = new G4Tubs("CAKE_SiliconWafer_finalshear_filled",0.0*mm, 160.0*mm, ((400./2) + 10.)*um, 0.*deg, 60.*deg);
    
    G4ThreeVector position_CAKE_SiliconWafer_finalshear = G4ThreeVector(0.0*mm, 0.0*mm, ((thickness_CAKE_PCB/2.) - (thickness_SiliconWafer_MMM/2.) + 10.)*um);
    
    G4VSolid* Solid_CAKE_finalshear1 = new G4SubtractionSolid("CAKE_finalshear1", Solid_CAKE_SiliconWafer_finalshear_filled, SiliconWafer_shear1, rm_SiliconWafer_shear1, position_CAKE_PCB_shear[0]);
    G4VSolid* Solid_CAKE_finalshear2 = new G4SubtractionSolid("CAKE_finalshear2", Solid_CAKE_finalshear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_CAKE_PCB_shear[1]);
    G4VSolid* Solid_CAKE_finalshear3 = new G4SubtractionSolid("CAKE_finalshear3", Solid_CAKE_finalshear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_CAKE_PCB_shear[2]);
    G4VSolid* Solid_CAKE_finalshear4 = new G4SubtractionSolid("CAKE_finalshear4", Solid_CAKE_finalshear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_CAKE_PCB_shear[3]);
    G4VSolid* Solid_CAKE_finalshear5 = new G4SubtractionSolid("CAKE_finalshear5", Solid_CAKE_finalshear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_CAKE_PCB_shear[4]);
    G4VSolid* Solid_CAKE_finalshear6 = new G4SubtractionSolid("CAKE_finalshear6", Solid_CAKE_finalshear5, SiliconWafer_shear1, rm_SiliconWafer_shear6, position_CAKE_PCB_shear[5]);
    
    //  purposely use Solid_CAKE_finalshear5 and not Solid_CAKE_finalshear6 to prevent logical error through boolean operation
    G4VSolid* Solid_CAKE_PCB_Final = new G4SubtractionSolid("CAKE_PCB_Final", Solid_CAKE_PCB_punched, Solid_CAKE_finalshear5, rm_SiliconWafer_shear, position_CAKE_SiliconWafer_finalshear);
    
    G4LogicalVolume* Logic_CAKE_PCB = new G4LogicalVolume(Solid_CAKE_PCB_Final, G4_PLASTIC_SC_VINYLTOLUENE_Material, "CAKE_PCB_Final", 0, 0, 0);
    
    
    
    /////////////////////////////////////////
    //          CAKE - Active Area
    
    G4double rotationphi;
    
    G4Tubs*             Solid_CAKE_RS_filled[16][8];
    G4VSolid*           Solid_CAKE_RS_shear[16][8];
    G4VSolid*           Solid_CAKE_RS[16][8];
    G4LogicalVolume*    Logic_CAKE_AA_RS[16][8];
    
    G4ThreeVector position_CAKE_Sector_shear[8][2];
    G4RotationMatrix* rm_CAKE_Sector_shear[8][2];
    
    for(G4int j=0; j<8; j++)
    {
        for(G4int l=0; l<2; l++)
        {
            rm_CAKE_Sector_shear[j][l] = new G4RotationMatrix();
        }
    }
    
    G4ThreeVector position_CAKE_Sector = G4ThreeVector(1.4448593*mm,  -0.0493525*mm, 0.0*mm);
    
    G4Box* CAKE_Sector_razor = new G4Box("CAKE_Sector_razor", (150.)*mm, (150.)*mm, (100/2)*mm);
    
    position_CAKE_Sector_shear[0][0] = G4ThreeVector(93.0470539*mm, -145.2243182*mm, 0.0*mm);
    position_CAKE_Sector_shear[0][1] = G4ThreeVector(58.2842550*mm,  162.2461746*mm, 0.0*mm);
    rm_CAKE_Sector_shear[0][0]->rotateZ(-3.0000732*deg);
    rm_CAKE_Sector_shear[0][1]->rotateZ(-9.8424272*deg);
    
    position_CAKE_Sector_shear[1][0] = G4ThreeVector(106.8796942*mm, -133.5912908*mm, 0.0*mm);
    position_CAKE_Sector_shear[1][1] = G4ThreeVector(35.1278912*mm,  167.0122283*mm, 0.0*mm);
    rm_CAKE_Sector_shear[1][0]->rotateZ(-9.8425735*deg);
    rm_CAKE_Sector_shear[1][1]->rotateZ(-16.6849271*deg);
    
    position_CAKE_Sector_shear[2][0] = G4ThreeVector(126.6444045*mm, -118.5243029*mm, 0.0*mm);
    position_CAKE_Sector_shear[2][1] = G4ThreeVector(19.1283325*mm,  171.8139616*mm, 0.0*mm);
    rm_CAKE_Sector_shear[2][0]->rotateZ(-16.6850728*deg);
    rm_CAKE_Sector_shear[2][1]->rotateZ(-23.5274272*deg);
    
    position_CAKE_Sector_shear[3][0] = G4ThreeVector(138.3487001*mm, -103.2511752*mm, 0.0*mm);
    position_CAKE_Sector_shear[3][1] = G4ThreeVector(-1.5603892*mm,  172.8208992*mm, 0.0*mm);
    rm_CAKE_Sector_shear[3][0]->rotateZ(-23.5275730*deg);
    rm_CAKE_Sector_shear[3][1]->rotateZ(-30.3700000*deg);
    
    position_CAKE_Sector_shear[4][0] = G4ThreeVector(147.5802745*mm, -87.2542869*mm, 0.0*mm);
    position_CAKE_Sector_shear[4][1] = G4ThreeVector(-24.6202415*mm,  169.5199469*mm, 0.0*mm);
    rm_CAKE_Sector_shear[4][0]->rotateZ(-30.3700732*deg);
    rm_CAKE_Sector_shear[4][1]->rotateZ(-37.2124996*deg);
    
    position_CAKE_Sector_shear[5][0] = G4ThreeVector(164.4752025*mm, -63.3162228*mm, 0.0*mm);
    position_CAKE_Sector_shear[5][1] = G4ThreeVector(-38.4388028*mm,  171.3807948*mm, 0.0*mm);
    rm_CAKE_Sector_shear[5][0]->rotateZ(-37.2124996*deg);
    rm_CAKE_Sector_shear[5][1]->rotateZ(-44.0549272*deg);
    
    position_CAKE_Sector_shear[6][0] = G4ThreeVector(169.2760955*mm, -44.7896689*mm, 0.0*mm);
    position_CAKE_Sector_shear[6][1] = G4ThreeVector(-59.5849311*mm,  164.3482428*mm, 0.0*mm);
    rm_CAKE_Sector_shear[6][0]->rotateZ(-44.0550729*deg);
    rm_CAKE_Sector_shear[6][1]->rotateZ(-50.8974265*deg);
    
    position_CAKE_Sector_shear[7][0] = G4ThreeVector(176.8824451*mm, -20.0264597*mm, 0.0*mm);
    position_CAKE_Sector_shear[7][1] = G4ThreeVector(-75.2407310*mm,  161.6241136*mm, 0.0*mm);
    rm_CAKE_Sector_shear[7][0]->rotateZ(-50.8974999*deg);
    rm_CAKE_Sector_shear[7][1]->rotateZ(-57.7399268*deg);
    
    G4double innerRadius, outerRadius;
    
    for(G4int j=0; j<16; j++)
    {
        innerRadius = 32.60 + ((6.3125 + 0.1)*j); // mm
        outerRadius = innerRadius + 6.3125; // mm
        
        for(G4int l=0; l<8; l++)
        {
            Solid_CAKE_RS_filled[j][l] = new G4Tubs("CAKE_RS_filled", innerRadius*mm, outerRadius*mm, (thickness_MMM_AA/2)*um, 0.*deg, 60.*deg);
            
            Solid_CAKE_RS_shear[j][l] = new G4SubtractionSolid("CAKE_RS_shear", Solid_CAKE_RS_filled[j][l], CAKE_Sector_razor, rm_CAKE_Sector_shear[l][0], position_CAKE_Sector_shear[l][0]);
            
            Solid_CAKE_RS[j][l] = new G4SubtractionSolid("CAKE_RS", Solid_CAKE_RS_shear[j][l], CAKE_Sector_razor, rm_CAKE_Sector_shear[l][1], position_CAKE_Sector_shear[l][1]);
            
            Logic_CAKE_AA_RS[j][l] = new G4LogicalVolume(Solid_CAKE_RS[j][l], G4_Si_Material, "Logic_CAKE_AA_RS", 0, 0, 0);
            
            
        }
    }
    
    
    /////////////////////////////////////////////
    //          CAKE - Active Area Punch
    
    G4Tubs*             Solid_CAKE_RS_punch_filled[16][8];
    G4VSolid*           Solid_CAKE_RS_punch_shear[16][8];
    G4VSolid*           Solid_CAKE_RS_punch[16][8];
    G4LogicalVolume*    Logic_CAKE_RS_punch[16][8];
    
    
    for(G4int j=0; j<16; j++)
    {
        innerRadius = 32.60 + ((6.3125 + 0.1)*j); // mm
        outerRadius = innerRadius + 6.3125; // mm
        
        for(G4int l=0; l<8; l++)
        {
            //  The "+ 300" is to avoid a logical error occuring during the boolean subtraction
            //Solid_CAKE_RS_punch_filled[j][l] = new G4Tubs("CAKE_RS_punch_filled", innerRadius*mm, outerRadius*mm, ((thickness_MMM_AA/2) + 1000)*um, 0.*deg, 60.*deg);
            Solid_CAKE_RS_punch_filled[j][l] = new G4Tubs("CAKE_RS_punch_filled", innerRadius*mm, outerRadius*mm, 2.*mm, 0.*deg, 60.*deg);
            
            Solid_CAKE_RS_punch_shear[j][l] = new G4SubtractionSolid("CAKE_RS_punch_shear", Solid_CAKE_RS_punch_filled[j][l], CAKE_Sector_razor, rm_CAKE_Sector_shear[l][0], position_CAKE_Sector_shear[l][0]);
            
            Solid_CAKE_RS_punch[j][l] = new G4SubtractionSolid("CAKE_RS_punch", Solid_CAKE_RS_punch_shear[j][l], CAKE_Sector_razor, rm_CAKE_Sector_shear[l][1], position_CAKE_Sector_shear[l][1]);
            
            Logic_CAKE_RS_punch[j][l] = new G4LogicalVolume(Solid_CAKE_RS_punch[j][l], G4_Si_Material, "Logic_CAKE_RS_punch", 0, 0, 0);
            
        }
    }
    
    
    
    /////////////////////////////////////////////////
    //          CAKE - Inter-strip Region
    
    //CAKE_Interstrip_rotm.rotateZ(3.*deg);
    
    G4ThreeVector position_CAKE_Interstrip = position_CAKE_Sector;
    
    G4RotationMatrix* rm_CAKE_Interstrip = new G4RotationMatrix();
    //rm_CAKE_Interstrip->rotateZ(3.*deg);
    
    //G4Box* Solid_CAKE_Interstrip_filled = new G4Box("CAKE_Interstrip_filled", (40/2)*cm, (40/2)*cm, (thickness_MMM_AA/2)*um);
    G4Tubs* Solid_CAKE_Interstrip_filled = new G4Tubs("CAKE_Interstrip_filled", 0.0*mm, 160.0*mm, (thickness_MMM_AA/2)*um, 0.*deg, (60.)*deg);
    
    G4VSolid* Solid_CAKE_Interstrip_Final[16][8];
    
    for(G4int j=0; j<16; j++)
    {
        for(G4int l=0; l<8; l++)
        {
            if(j==0 && l==0)
            {
                Solid_CAKE_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_CAKE_Interstrip_Final", Solid_CAKE_Interstrip_filled, Solid_CAKE_RS_punch[j][l], rm_CAKE_Interstrip, position_CAKE_Interstrip);
                
            }
            if(j!=0 && l==0)
            {
                Solid_CAKE_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_CAKE_Interstrip_Final", Solid_CAKE_Interstrip_Final[j-1][7], Solid_CAKE_RS_punch[j][l], rm_CAKE_Interstrip, position_CAKE_Interstrip);
                
            }
            if(l!=0)
            {
                Solid_CAKE_Interstrip_Final[j][l] = new G4SubtractionSolid("Solid_CAKE_Interstrip_Final", Solid_CAKE_Interstrip_Final[j][l-1], Solid_CAKE_RS_punch[j][l], rm_CAKE_Interstrip, position_CAKE_Interstrip);
                
            }
        }
    }
    

    
    ////    STD Shearing attempt
    G4VSolid* Solid_CAKE_Interstrip_shear1 = new G4SubtractionSolid("Solid_CAKE_Interstrip_shear1", Solid_CAKE_Interstrip_Final[15][7], SiliconWafer_shear1, rm_SiliconWafer_shear1, position_CAKE_PCB_shear[0]);
    G4VSolid* Solid_CAKE_Interstrip_shear2 = new G4SubtractionSolid("Solid_CAKE_Interstrip_shear2", Solid_CAKE_Interstrip_shear1, SiliconWafer_shear1, rm_SiliconWafer_shear2, position_CAKE_PCB_shear[1]);
    G4VSolid* Solid_CAKE_Interstrip_shear3 = new G4SubtractionSolid("Solid_CAKE_Interstrip_shear3", Solid_CAKE_Interstrip_shear2, SiliconWafer_shear1, rm_SiliconWafer_shear3, position_CAKE_PCB_shear[2]);
    G4VSolid* Solid_CAKE_Interstrip_shear4 = new G4SubtractionSolid("Solid_CAKE_Interstrip_shear4", Solid_CAKE_Interstrip_shear3, SiliconWafer_shear1, rm_SiliconWafer_shear4, position_CAKE_PCB_shear[3]);
    G4VSolid* Solid_CAKE_Interstrip_shear5 = new G4SubtractionSolid("Solid_CAKE_Interstrip_shear5", Solid_CAKE_Interstrip_shear4, SiliconWafer_shear1, rm_SiliconWafer_shear5, position_CAKE_PCB_shear[4]);
    
    G4VSolid* Solid_CAKE_Interstrip1 = new G4SubtractionSolid("Solid_CAKE_Interstrip_shear3", Solid_CAKE_Interstrip_Final[15][7], SiliconWafer_shear1, rm_SiliconWafer_shear6, position_CAKE_PCB_shear[5]);

    
    G4LogicalVolume * Logic_CAKE_SiliconWafer[numberOf_CAKE];
    
    for(G4int i=0; i<numberOf_CAKE; i++)
    {
        Logic_CAKE_SiliconWafer[i] = new G4LogicalVolume(Solid_CAKE_SiliconWafer_Final, G4_Si_Material, "Logic_CAKE_Interstrip", 0, 0, 0);
    }
    
    
    
    
    
    /////////////////////////////////////////////////////
    //          CAKE 2M Windows, Front and Back
    
    G4Tubs* Solid_CAKE_2M_Front = new G4Tubs("CAKE_2M_Front",0.0*mm, 140.0*mm, (thickness_2M_MMM_Front/2)*um, 0.*deg, (54.74)*deg);
    G4Tubs* Solid_CAKE_2M_Back = new G4Tubs("CAKE_2M_Back",0.0*mm, 140.0*mm, (thickness_2M_MMM_Back/2)*um, 0.*deg, (54.74)*deg);
    
    G4LogicalVolume* Logic_CAKE_2M_Front = new G4LogicalVolume(Solid_CAKE_2M_Front, G4_Al_Material, "Logic_CAKE_DL_Front", 0, 0, 0);
    G4LogicalVolume* Logic_CAKE_2M_Back = new G4LogicalVolume(Solid_CAKE_2M_Back, G4_Al_Material, "Logic_CAKE_DL_Back", 0, 0, 0);
    

    
    
    /////////////////////////////////////////
    //          CAKE Asm
    
    G4Tubs* Solid_CAKE_Asm_filled = new G4Tubs("Solid_CAKE_Asm_filled",0.0*mm, 170.0*mm, (thickness_CAKE_PCB/2)*um, 0*deg, (57.74)*deg);
    
    G4VSolid* Solid_CAKE_Asm_shear1 = new G4SubtractionSolid("CAKE_Asm_shear1", Solid_CAKE_Asm_filled, PCB_razor1, rm1_CAKE_PCB, position_CAKE_PCB_razor1[0]);
    G4VSolid* Solid_CAKE_Asm_shear2 = new G4SubtractionSolid("CAKE_Asm_shear2", Solid_CAKE_Asm_shear1, PCB_razor1, rm2_CAKE_PCB, position_CAKE_PCB_razor1[1]);
    G4VSolid* Solid_CAKE_Asm_shear3 = new G4SubtractionSolid("CAKE_Asm_shear3", Solid_CAKE_Asm_shear2, PCB_razor1, rm3_CAKE_PCB, position_CAKE_PCB_razor1[2]);
    G4VSolid* Solid_CAKE_Asm = new G4SubtractionSolid("Solid_CAKE_Asm", Solid_CAKE_Asm_shear3, PCB_razor1, rm4_CAKE_PCB, position_CAKE_PCB_razor1[3]);
    
    G4LogicalVolume*    Logic_CAKE_Asm[numberOf_CAKE];
    
    for(G4int i=0; i<numberOf_CAKE; i++)
    {
        Logic_CAKE_Asm[i] = new G4LogicalVolume(Solid_CAKE_Asm, G4_Galactic_Material, "Logic_CAKE_Asm", 0, 0, 0);
    }
    
    
    
    
    
    ////////////////////////////////////////////////////
    //               CAKE INITIALIZATION             //
    ////////////////////////////////////////////////////
    
    //G4ThreeVector offset_CAKE_AA = G4ThreeVector(0*cm, 0*cm, (thickness_CAKE_PCB/2 - thickness_DL_MMM_Front - thickness_MMM_AA/2)*um);
    //G4ThreeVector offset_CAKE_AA = position_CAKE_Sector + G4ThreeVector(0*cm, 0*cm, (thickness_CAKE_PCB/2 - thickness_DL_MMM_Front - (thickness_MMM_AA/2))*um);
    
    //G4ThreeVector offset_CAKE_SiliconWafer = G4ThreeVector(0*cm, 0*cm, (thickness_CAKE_PCB/2 - thickness_DL_MMM_Front - (thickness_MMM_AA/2))*um);
    G4ThreeVector offset_CAKE_SiliconWafer = G4ThreeVector(0*cm, 0*cm, (thickness_CAKE_PCB/2. - (400./2.))*um);
    
    G4ThreeVector offset_CAKE_AA = position_CAKE_Sector + G4ThreeVector(0*cm, 0*cm, ((400./2) - (thickness_MMM_AA/2) - thickness_DL_MMM_Front)*um);
    
    G4ThreeVector offset_CAKE_DL_Front = offset_CAKE_AA + G4ThreeVector(0.0*mm,0.0*mm, ((thickness_MMM_AA/2 + thickness_DL_MMM_Front/2))*um);
    G4ThreeVector offset_CAKE_DL_Back = offset_CAKE_AA + G4ThreeVector(0.0*mm,0.0*mm, -((thickness_MMM_AA/2 + thickness_DL_MMM_Front/2))*um);
    
    G4ThreeVector offset_CAKE_2M_Front = offset_CAKE_DL_Front + G4ThreeVector(0.0*mm,0.0*mm, ((thickness_DL_MMM_Front/2 + thickness_2M_MMM_Front/2))*um);
    G4ThreeVector offset_CAKE_2M_Back = offset_CAKE_DL_Back + G4ThreeVector(0.0*mm,0.0*mm, -((thickness_DL_MMM_Back/2 + thickness_2M_MMM_Back/2))*um);
    
    //G4ThreeVector offset_CAKE_Interstrip = position_CAKE_Sector + G4ThreeVector(0.0*mm,0.0*mm, ((thickness_MMM_AA/2 + thickness_DL_MMM_Front/2))*um);
    //G4ThreeVector offset_CAKE_Interstrip = position_CAKE_Sector + G4ThreeVector(0*cm, 0*cm, (thickness_CAKE_PCB/2 - thickness_DL_MMM_Front - (thickness_MMM_AA/2))*um);
    
    //    G4ThreeVector offset_CAKE_Interstrip = G4ThreeVector(2.4748830*mm, 1.3857888*mm, (thickness_CAKE_PCB/2 - thickness_DL_MMM_Front - (thickness_MMM_AA/2))*um);
    
    CAKE_SiliconWafer_transform = G4Transform3D(CAKE_SiliconWafer_rotm, offset_CAKE_SiliconWafer);
    
    
    for(G4int i=0; i<numberOf_CAKE; i++)
    {
        CAKE_transform[i] = G4Transform3D(CAKE_rotm[i],CAKE_AA_CentrePosition[i]);
        
        
        if(CAKE_Presence[i])
        {
            ///////////////////////////////////////
            //      CAKE AA - Rings and Sectors
            
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<8; l++)
                {
                    CAKE_AA_RS_transform[l] = G4Transform3D(CAKE_AA_RS_rotm[l], offset_CAKE_AA);
                    
                    PhysiCAKE_AA_RS = new G4PVPlacement(CAKE_AA_RS_transform[l],
                                                         Logic_CAKE_AA_RS[j][l],
                                                         "CAKE_AA_RS", // its name
                                                         Logic_CAKE_SiliconWafer[i],
                                                         false, // no boolean operations
                                                         i*128 + l + j*8,  // copy number
                                                         fCheckOverlaps); // checking overlaps
                    
                }
            }
            
            
            ////////////////////////////////
            //      CAKE INTERSTRIP
            
            PhysiCAKE_SiliconWafer = new G4PVPlacement(CAKE_SiliconWafer_transform,
                                                        Logic_CAKE_SiliconWafer[i],
                                                        "CAKE_SiliconWafer", // its name
                                                        Logic_CAKE_Asm[i],
                                                        false, // no boolean operations
                                                        i,  // copy number
                                                        fCheckOverlaps); // checking overlaps
            
            
            ////////////////////////////////
            //      CAKE PCB
            
            new G4PVPlacement(0,               // no rotation
                              G4ThreeVector(), // at (x,y,z)
                              Logic_CAKE_PCB,       // its logical volume
                              "CAKE_PCB",       // its name
                              Logic_CAKE_Asm[i],         // its mother  volume
                              false,           // no boolean operations
                              0,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            /*
             ////////////////////////////////////////////
             //      CAKE 2M Windows, Front and Back
             
             CAKE_2M_transform[0] = G4Transform3D(CAKE_DL_2M_rotm, offset_CAKE_2M_Front);
             CAKE_2M_transform[1] = G4Transform3D(CAKE_DL_2M_rotm, offset_CAKE_2M_Back);
             
             new G4PVPlacement(CAKE_2M_transform[0],
             Logic_CAKE_2M_Front,       // its logical volume
             "Physical_CAKE_2M_Front",       // its name
             Logic_CAKE_Asm[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             
             new G4PVPlacement(CAKE_2M_transform[1],
             Logic_CAKE_2M_Back,       // its logical volume
             "Physical_CAKE_2M_Back",       // its name
             Logic_CAKE_Asm[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             
             
             ////////////////////////////////////////////
             //      CAKE - Dead Layers, Front and Back
             
             CAKE_DL_transform[0] = G4Transform3D(CAKE_DL_2M_rotm, offset_CAKE_DL_Front);
             CAKE_DL_transform[1] = G4Transform3D(CAKE_DL_2M_rotm, offset_CAKE_DL_Back);
             
             new G4PVPlacement(CAKE_DL_transform[0],
             Logic_CAKE_DL_Front,       // its logical volume
             "Physical_CAKE_DL_Front",       // its name
             Logic_CAKE_Asm[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             
             new G4PVPlacement(CAKE_DL_transform[1],
             Logic_CAKE_DL_Front,       // its logical volume
             "Physical_CAKE_DL_Back",       // its name
             Logic_CAKE_Asm[i],         // its mother  volume
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             */
            ////////////////////////////////
            //      CAKE Asm
            
            new G4PVPlacement(CAKE_transform[i], // transformation
                              Logic_CAKE_Asm[i],
                              "Physical_CAKE_Asm",
                              LogicVacuumChamber,
                              false,           // no boolean
                              0,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
    }
    
    
    
    /////////////////////////////////////////////
    //             W1 DEFINITION               //
    /////////////////////////////////////////////
    
    
    /////////////////////////////////////////
    //          W1 - PCB
    
    G4double thickness_W1_PCB = 1.6e+3; // um
    G4double W1_AA_length = 50.0; // mm: the active area is 50.0mm x 50.0mm
    
    G4Box* Solid_W1_PCB_filled = new G4Box("Solid_W1_PCB_filled", (60.0/2.0)*mm, (60.0/2.0)*mm, (thickness_W1_PCB/2.0)*um);

    G4Box* Solid_W1_PCB_punch = new G4Box("Solid_W1_PCB_punch", (W1_AA_length/2.0)*mm, (W1_AA_length/2.0)*mm, thickness_W1_PCB*um);

    G4ThreeVector position_W1_PCB_punch = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    G4RotationMatrix* rm1_W1_PCB = new G4RotationMatrix();
    
    G4VSolid* Solid_W1_PCB = new G4SubtractionSolid("Solid_W1_PCB", Solid_W1_PCB_filled, Solid_W1_PCB_punch, rm1_W1_PCB, position_W1_PCB_punch);

    G4LogicalVolume* Logic_W1_PCB = new G4LogicalVolume(Solid_W1_PCB, G4_PLASTIC_SC_VINYLTOLUENE_Material, "Logic_W1_PCB", 0, 0, 0);
    
    ////    W1 - Asm
    G4LogicalVolume*    Logic_W1_Asm[numberOf_CAKE];
    
    for(G4int i=0; i<numberOf_W1; i++)
    {
        Logic_W1_Asm[i] = new G4LogicalVolume(Solid_W1_PCB_filled, G4_Galactic_Material, "Logic_W1_Asm", 0, 0, 0);
    }

    
    /////////////////////////////////////////
    //          W1 - Silicon Wafer
    
    G4double thickness_SiliconWafer_W1 = 140.0; // um
    
    G4double thickness_DL_W1_Front = 0.3; // um
    G4double thickness_DL_W1_Back = 0.3; // um
    
    G4double thickness_2M_W1_Front = 0.3; // um
    G4double thickness_2M_W1_Back = 0.3; // um
    
    G4double W1_element_pitch = 3100.0; // um
    G4double W1_element_width = 3000.0; // um
    G4double interstrip_width = 100.0; // um
    
    G4Box* Solid_W1_SiliconWafer_filled = new G4Box("Solid_W1_SiliconWafer_filled", (W1_AA_length/2.0)*mm, (W1_AA_length/2.0)*mm, (thickness_SiliconWafer_W1/2.0)*um);

    ////    Active Area (Silicon) thickness
    G4double thickness_W1_AA = thickness_SiliconWafer_W1 - thickness_DL_W1_Front - thickness_DL_W1_Back; // um
    G4Box* Solid_W1_AA = new G4Box("Solid_W1_AA", (W1_element_width/2.0)*um, (W1_element_width/2.0)*um, (thickness_W1_AA/2.0)*um);
    G4LogicalVolume* Logic_W1_AA = new G4LogicalVolume(Solid_W1_AA, G4_Si_Material, "Logic_W1_AA", 0, 0, 0);

    ////    Dead Layers
    G4Box* Solid_W1_DL_Front = new G4Box("Solid_W1_DL_Front", (W1_AA_length/2.0)*mm, (W1_AA_length/2.0)*mm, (thickness_DL_W1_Front/2.0)*um);
    G4LogicalVolume* Logic_W1_DL_Front = new G4LogicalVolume(Solid_W1_DL_Front, G4_Si_Material, "Logic_W1_DL_Front", 0, 0, 0);

    G4Box* Solid_W1_DL_Back = new G4Box("Solid_W1_DL_Back", (W1_AA_length/2.0)*mm, (W1_AA_length/2.0)*mm, (thickness_DL_W1_Back/2.0)*um);
    G4LogicalVolume* Logic_W1_DL_Back = new G4LogicalVolume(Solid_W1_DL_Back, G4_Si_Material, "Logic_W1_DL_Back", 0, 0, 0);

    
    ////    W1 - Silicon Wafer Asm (includes dead layers and windows)
    G4LogicalVolume*    Logic_W1_SiliconWafer_Asm[numberOf_W1];
    
    for(G4int i=0; i<numberOf_W1; i++)
    {
        Logic_W1_SiliconWafer_Asm[i] = new G4LogicalVolume(Solid_W1_SiliconWafer_filled, G4_Si_Material, "Logic_W1_SiliconWafer_Asm", 0, 0, 0);
    }



    //////////////////////////////////////////////////
    //              W1 INITIALIZATION               //
    //////////////////////////////////////////////////

    
    G4ThreeVector offset_W1_SiliconWafer = G4ThreeVector(0.0*mm, 0.0*mm, ((thickness_W1_PCB/2.0) - (thickness_SiliconWafer_W1/2.0))*um);
    //G4ThreeVector offset_W1_SiliconWafer = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*um);

    G4ThreeVector offset_W1_AA = G4ThreeVector(0.0*mm, 0.0*mm, ((-(thickness_DL_W1_Front) + (thickness_DL_W1_Back))/2.0)*um);
    
    G4ThreeVector offset_W1_DL_Front = G4ThreeVector(0.0*mm, 0.0*mm, ((thickness_SiliconWafer_W1/2.0) - (thickness_DL_W1_Front/2.0))*um);
    G4ThreeVector offset_W1_DL_Back = G4ThreeVector(0.0*mm, 0.0*mm, (-(thickness_SiliconWafer_W1/2.0) + (thickness_DL_W1_Back/2.0))*um);
    
    
    W1_SiliconWafer_transform = G4Transform3D(W1_SiliconWafer_rotm, offset_W1_SiliconWafer);

    G4ThreeVector position_W1_AA[16][16];
    
    for(G4int i=0; i<numberOf_W1; i++)
    {
        W1_transform[i] = G4Transform3D(W1_rotm[i],W1_AA_CentrePosition[i]);
        
        if(W1_Presence[i])
        {
            ////    W1 - Active Area
            for(G4int j=0; j<16; j++)
            {
                for(G4int l=0; l<16; l++)
                {
                    
                    position_W1_AA[j][l] = G4ThreeVector((-7.5 + l)*(W1_element_width + interstrip_width)*um, (-7.5 + j)*(W1_element_width + interstrip_width)*um, 0.0) + offset_W1_AA;
                    
                    PhysiW1_AA = new G4PVPlacement(0,               // no rotation
                                                   position_W1_AA[j][l], // at (x,y,z)
                                                   Logic_W1_AA,       // its logical volume
                                                   "W1_AA",       // its name
                                                   Logic_W1_SiliconWafer_Asm[i], // its mother  volume
                                                   false,           // no boolean operations
                                                   (i*16*16) + (j*16) + l, // copy number
                                                   fCheckOverlaps); // checking overlaps
                    
                    
                    /*
                    if(j!=2 && l>3)
                    {
                        position_W1_AA[j][l] = G4ThreeVector((-7.5 + l)*(W1_element_width + interstrip_width)*um, (-7.5 + j)*(W1_element_width + interstrip_width)*um, 0.0) + offset_W1_AA;
                        
                        PhysiW1_AA = new G4PVPlacement(0,               // no rotation
                                                       position_W1_AA[j][l], // at (x,y,z)
                                                       Logic_W1_AA,       // its logical volume
                                                       "W1_AA",       // its name
                                                       Logic_W1_SiliconWafer_Asm[i], // its mother  volume
                                                       false,           // no boolean operations
                                                       (i*16*16) + (j*16) + l, // copy number
                                                       fCheckOverlaps); // checking overlaps
                        
                    }
                    */
                    
                }
            }
            
            ////    W1 - Dead Layer, Front
            new G4PVPlacement(0,               // no rotation
                              offset_W1_DL_Front, // at (x,y,z)
                              Logic_W1_DL_Front, // its logical volume
                              "W1_DL_Front",       // its name
                              Logic_W1_SiliconWafer_Asm[i], // its mother  volume
                              false,           // no boolean operations
                              i, // copy number
                              fCheckOverlaps); // checking overlaps

            ////    W1 - Dead Layer, Back
            new G4PVPlacement(0,               // no rotation
                              offset_W1_DL_Back, // at (x,y,z)
                              Logic_W1_DL_Back, // its logical volume
                              "W1_DL_Back",       // its name
                              Logic_W1_SiliconWafer_Asm[i], // its mother  volume
                              false,           // no boolean operations
                              i, // copy number
                              fCheckOverlaps); // checking overlaps

            
            ////    W1 - Combined Silicon Wafer Assembly
            PhysiW1_SiliconWafer = new G4PVPlacement(0,               // no rotation
                                                     offset_W1_SiliconWafer, // at (x,y,z)
                                                     Logic_W1_SiliconWafer_Asm[i], // its logical volume
                                                     "PhysiW1_SiliconWafer",       // its name
                                                     Logic_W1_Asm[i], // its mother  volume
                                                     false,           // no boolean operations
                                                     i, // copy number
                                                     fCheckOverlaps); // checking overlaps

            
            
            
            new G4PVPlacement(W1_transform[i], // transformation
                              Logic_W1_Asm[i],
                              "Physical_W1_Asm",
                              LogicVacuumChamber,
                              false,           // no boolean
                              0,               // copy number
                              fCheckOverlaps); // checking overlaps

            
            
        }
    }
    
    
    
    
    
    
    
    
    //////////////////////////////////////
    //          VDC DEFINITION          //
    //////////////////////////////////////
    
    ///////////////////////////////
    //      VDC - Asm
    G4Box* Solid_VDC_Asm_USDS = new G4Box("Solid_Stesalit_Frame", (936./2)*mm, (240./2)*mm, (24./2)*mm);
    G4Box* Solid_VDC_Asm = new G4Box("Solid_Stesalit_Frame", (936./2)*mm, (240./2)*mm, 78./2*mm);
    
    G4LogicalVolume* Logic_VDC_SenseRegion_USDS[numberOf_VDC][2];
    G4LogicalVolume* Logic_VDC_Asm[numberOf_VDC];
    
    for(G4int i=0; i<2; i++)
    {
        Logic_VDC_Asm[i] = new G4LogicalVolume(Solid_VDC_Asm, G4_AIR_Material, "Logic_VDC_Asm",0,0,0);
        
        for(G4int j=0; j<2; j++)
        {
            Logic_VDC_SenseRegion_USDS[i][j] = new G4LogicalVolume(Solid_VDC_Asm_USDS, VDC_SR_Gas_Material, "Logic_VDC_SenseRegion_USDS",0,0,0);
        }
    }
    
    ///////////////////////////////////////////////
    //      VDC - Stesalit Frame Templates
    G4double thickness_Stesalit_Frame = 8.0; // mm
    G4double thickness_StesalitPCB_Frame = 5.5; // mm
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT STANDARD FRAME
    G4Box* Solid_Stesalit_StdFrame = new G4Box("Solid_Stesalit_StdFrame", (936./2)*mm, (240./2)*mm, (8./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT STANDARD PCB FRAME
    G4Box* Solid_StesalitPCB_StdFrame = new G4Box("Solid_StesalitPCB_StdFrame", (936./2)*mm, (240./2)*mm, (5.5/2)*mm);
    G4Box* Solid_Wire_subtraction = new G4Box("Solid_Wire_subtraction", (2000./2)*mm, (1000./2)*mm, (10./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - PCB FRAME
    G4Box* Solid_PCB_StdFrame = new G4Box("Solid_PCB_StdFrame", (936./2)*mm, (240./2)*mm, (2.5/2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - ALUMINIUM OUTER FRAME
    G4Box* Solid_Al_StdFrame = new G4Box("Solid_Al_StdFrame", (936./2)*mm, (240./2)*mm, (15./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - FRAME FINISHING
    G4ThreeVector   position_punch[4];
    
    G4Box* punch1 = new G4Box("punch1", (825./2)*mm, (100./2)*mm, (20./2)*mm);
    G4Box* punch2 = new G4Box("punch2", (809./2)*mm, (131./2)*mm, (20./2)*mm);
    G4Box* punch3 = new G4Box("punch3", (800./2)*mm, (100./2)*mm, (20./2)*mm);
    G4Box* punch4 = new G4Box("punch4", (852./2)*mm, (100./2)*mm, (20./2)*mm);
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT HV Frames, Upstream and Downstream
    position_punch[0] = G4ThreeVector((-9.5)*mm, 0.0*mm, 0.0*mm);
    
    G4VSolid* Solid_Stesalit_HV_USDS_Frame = new G4SubtractionSolid("Solid_Stesalit_HV_USDS_Frame", Solid_Stesalit_StdFrame, punch1, 0, position_punch[0]);
    
    ///////////////////////////////////////////////
    position_punch[0] = G4ThreeVector(-(825./2 - 8.7)*mm, 0.0*mm, (0.)*mm);
    
    G4Box* extrusion1 = new G4Box("extrusion1", (830./2)*mm, (113./2)*mm, (20./2)*mm);
    
    G4VSolid* Solid_VDC_HVFrame_USDS2 = new G4SubtractionSolid("Solid_VDC_HVFrame_USDS2", Solid_Stesalit_HV_USDS_Frame, extrusion1, 0, G4ThreeVector((-9.5-2.5)*mm, 0.0*mm, (4+10-2)*mm));
    
    ///////////////////////////////////////////////
    
    G4RotationMatrix* rm_punch = new G4RotationMatrix();
    rm_punch->rotateY(-35.*deg);
    
    G4VSolid* Solid_VDC_HVFrame_USDS = new G4SubtractionSolid("Solid_VDC_HVFrame_USDS", Solid_VDC_HVFrame_USDS2, punch1, rm_punch, position_punch[0]);
    
    G4LogicalVolume*    Logic_VDC_HV_Frame[3];
    
    Logic_VDC_HV_Frame[0] =  new G4LogicalVolume(Solid_VDC_HVFrame_USDS, G4_Al_Material, "Logic_VDC_HV_Frame",0,0,0);
    Logic_VDC_HV_Frame[2] = new G4LogicalVolume(Solid_VDC_HVFrame_USDS, G4_Al_Material, "Logic_VDC_HV_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT GAS CELL
    position_punch[1] = G4ThreeVector(-1.5*mm, -1.5*mm, 0.0*mm);
    
    G4VSolid* Solid_Stesalit_GasFrame = new G4SubtractionSolid("Solid_Stesalit_GasFrame", Solid_Stesalit_StdFrame, punch2, 0, position_punch[1]);
    
    position_punch[1] = G4ThreeVector(-(809./2 + 0.1 - 8.192)*mm, -1.5*mm, 1.3*mm);
    
    G4VSolid* Solid_VDC_GasFrame = new G4SubtractionSolid("Solid_VDC_GasFrame", Solid_Stesalit_GasFrame, punch2, rm_punch, position_punch[1]);
    
    G4LogicalVolume* Logic_VDC_GasFrame = new G4LogicalVolume(Solid_VDC_GasFrame, G4_Al_Material, "Logic_VDC_GasFrame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT X-WIRE/U-WIRE FRAME
    position_punch[2] = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    
    G4VSolid* Solid_VDC_Stesalit_XU_Frame = new G4SubtractionSolid("Solid_VDC_Stesalit_XU_Frame", Solid_StesalitPCB_StdFrame, punch3, 0, position_punch[2]);
    
    G4VSolid* Solid_Wire_subtraction2 = new G4SubtractionSolid("Solid_VDC_Stesalit_XU_Frame", Solid_Wire_subtraction, punch3, 0, position_punch[2]);
    
    G4LogicalVolume* Logic_VDC_XU_Frame = new G4LogicalVolume(Solid_VDC_Stesalit_XU_Frame, G4_Al_Material, "Logic_VDC_XU_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - PCB X-WIRE/U-WIRE FRAME
    position_punch[2] = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    
    G4VSolid* Solid_VDC_PCB_XU_Frame = new G4SubtractionSolid("Solid_VDC_PCB_XU_Frame", Solid_PCB_StdFrame, punch3, 0, position_punch[2]);
    
    G4LogicalVolume* Logic_VDC_XU_PCBFrame = new G4LogicalVolume(Solid_VDC_PCB_XU_Frame, G4_Al_Material, "Logic_VDC_XU_PCBFrame",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - STESALIT HV Frame, MIDDLE
    position_punch[3] = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
    
    G4VSolid* Solid_Stesalit_HV_MID_Frame = new G4SubtractionSolid("Solid_Stesalit_HV_MID_Frame", Solid_Stesalit_StdFrame, punch4, 0, position_punch[3]);
    
    Logic_VDC_HV_Frame[1] = new G4LogicalVolume(Solid_Stesalit_HV_MID_Frame, G4_Al_Material, "Logic_VDC_HV_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - ALUMINIUM OUTER FRAME
    G4Box* AlFrame_punch = new G4Box("AlFrame_punch", (834./2)*mm, (100./2)*mm, (20./2)*mm);
    
    G4Box* Solid_Aluminium_Frame = new G4Box("Solid_Stesalit_Frame", (936./2)*mm, (240./2)*mm, (15./2)*mm);
    
    G4VSolid* Solid_VDC_Al_Frame_shear2perf = new G4SubtractionSolid("Solid_VDC_Al_Frame_shear2perf", Solid_Aluminium_Frame, AlFrame_punch, 0, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm));
    
    G4VSolid* Solid_VDC_Al_Frame = new G4SubtractionSolid("Solid_VDC_Al_Frame", Solid_VDC_Al_Frame_shear2perf, AlFrame_punch, rm_punch, G4ThreeVector(-(834./2 + 14.0 - 8.192)*mm, 0.0*mm, 13.0*mm));
    
    G4LogicalVolume* Logic_VDC_Al_Frame = new G4LogicalVolume(Solid_VDC_Al_Frame, G4_Al_Material, "Logic_VDC_Al_Frame",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - HIGH VOLTAGE PLANE
    G4Box* Solid_VDC_HV_Plane_US = new G4Box("Solid_VDC_HV_Plane_US", (809./2)*mm, (100./2)*mm, (20./2)*um);
    G4Box* Solid_VDC_HV_Plane_MID = new G4Box("Solid_VDC_HV_Plane_MID", (800./2)*mm, (100./2)*mm, (20./2)*um);
    G4Box* Solid_VDC_HV_Plane_DS = new G4Box("Solid_VDC_HV_Plane_DS", (800./2)*mm, (100./2)*mm, (20./2)*um);
    
    G4LogicalVolume*    Logic_VDC_HV_Plane[3];
    
    Logic_VDC_HV_Plane[0] = new G4LogicalVolume(Solid_VDC_HV_Plane_US, G4_Al_Material, "Logic_VDC_HV_Plane_US",0,0,0);
    Logic_VDC_HV_Plane[1] = new G4LogicalVolume(Solid_VDC_HV_Plane_MID, G4_Al_Material, "Logic_VDC_HV_Plane_MID",0,0,0);
    Logic_VDC_HV_Plane[2] = new G4LogicalVolume(Solid_VDC_HV_Plane_DS, G4_Al_Material, "Logic_VDC_HV_Plane_DS",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - MYLAR WINDOW
    G4Box* Solid_VDC_MYLAR_Plane = new G4Box("Solid_VDC_MYLAR_Plane", (834./2)*mm, (100./2)*mm, (25./2)*um);
    
    G4LogicalVolume* Logic_VDC_MYLAR_Plane = new G4LogicalVolume(Solid_VDC_MYLAR_Plane, G4_Mylar_Material, "Logic_VDC_MYLAR_Plane",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - X WIRES
    G4Tubs* Solid_VDC_X_WIRE = new G4Tubs("Solid_VDC_X_WIRE", 0.*um, 20.*um, 100./2*mm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_VDC_X_WIRE = new G4LogicalVolume(Solid_VDC_X_WIRE, G4_W_Material, "Logic_VDC_X_WIRE",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - X GUARD WIRES
    G4Tubs* Solid_VDC_X_GUARDWIRE = new G4Tubs("Solid_VDC_X_GUARDWIRE", 0.*um, 50.*um, 100./2*mm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_VDC_X_GUARDWIRE = new G4LogicalVolume(Solid_VDC_X_GUARDWIRE, G4_W_Material, "Logic_VDC_X_GUARDWIRE",0,0,0);
    
    
    ///////////////////////////////////////////////
    //      VDC - X GUARD WIRES, THICK
    G4Tubs* Solid_VDC_X_GUARDWIRE_THICK = new G4Tubs("Solid_VDC_X_GUARDWIRE_THICK", 0.*um, 100.*um, 100./2*mm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_VDC_X_GUARDWIRE_Thick = new G4LogicalVolume(Solid_VDC_X_GUARDWIRE_THICK, G4_W_Material, "Logic_VDC_X_GUARDWIRE_Thick",0,0,0);
    
    ///////////////////////////////////////////////
    //      VDC - U WIRES
    G4Tubs* Solid_VDC_U_WIRE_Full = new G4Tubs("Solid_VDC_U_WIRE", 0.*um, 20.*um, 150./2*mm, 0.*deg, 360.*deg);
    
    G4RotationMatrix* rm_U_WIRE = new G4RotationMatrix();
    rm_U_WIRE->rotateY(90.*deg);
    rm_U_WIRE->rotateZ((90.-40.)*deg);
    
    G4VSolid* Solid_VDC_U_WIRE;
    G4ThreeVector   offset_VDC_U_WIRE_subtraction;
    
    G4LogicalVolume* Logic_VDC_U_WIRE[143];
    
    for(G4int k=0; k<143; k++)
    {
        offset_VDC_U_WIRE_subtraction = G4ThreeVector(0.0*mm, (-k + 71)*4*mm, (k - 71)*(4/tan(50*deg))*mm);
        
        Solid_VDC_U_WIRE = new G4SubtractionSolid("Solid_VDC_U_WIRE", Solid_VDC_U_WIRE_Full, Solid_Wire_subtraction2, rm_U_WIRE, offset_VDC_U_WIRE_subtraction);
        
        Logic_VDC_U_WIRE[k] = new G4LogicalVolume(Solid_VDC_U_WIRE, G4_W_Material, "Logic_VDC_U_WIRE",0,0,0);
    }
    
    ///////////////////////////////////////////////
    //      VDC - U GUARD WIRES
    G4Tubs* Solid_VDC_U_GUARDWIRE_Full = new G4Tubs("Solid_VDC_U_GUARDWIRE", 0.*um, 50.*um, 150./2*mm, 0.*deg, 360.*deg);
    
    G4VSolid* Solid_VDC_U_GUARDWIRE;
    G4LogicalVolume* Logic_VDC_U_GUARDWIRE[144];
    
    for(G4int k=0; k<144; k++)
    {
        offset_VDC_U_WIRE_subtraction = G4ThreeVector(0.0*mm, (-k + 71.5)*4*mm, (k - 71.5)*(4/tan(50*deg))*mm);
        
        Solid_VDC_U_GUARDWIRE = new G4SubtractionSolid("Solid_VDC_U_WIRE", Solid_VDC_U_GUARDWIRE_Full, Solid_Wire_subtraction2, rm_U_WIRE, offset_VDC_U_WIRE_subtraction);
        
        Logic_VDC_U_GUARDWIRE[k] = new G4LogicalVolume(Solid_VDC_U_GUARDWIRE, G4_W_Material, "Logic_VDC_U_GUARDWIRE",0,0,0);
        
    }
    
    ////////////////////////////////////////////////
    //      VDC - U GUARD WIRES, Thick
    G4Tubs* Solid_VDC_U_GUARDWIRE_Thick_Full = new G4Tubs("Solid_VDC_U_GUARDWIRE_Thick", 0.*um, 100.*um, 150./2*mm, 0.*deg, 360.*deg);
    
    G4VSolid* Solid_VDC_U_GUARDWIRE_Thick;
    
    G4LogicalVolume*    Logic_VDC_U_GUARDWIRE_Thick[2];
    
    for(G4int k=0; k<2; k++)
    {
        offset_VDC_U_WIRE_subtraction = G4ThreeVector(0.0*mm, (-k*144. + 72.)*4*mm, (k*144. - 72.)*(4/tan(50*deg))*mm);
        
        Solid_VDC_U_GUARDWIRE_Thick = new G4SubtractionSolid("Solid_VDC_U_WIRE_Thick", Solid_VDC_U_GUARDWIRE_Thick_Full, Solid_Wire_subtraction2, rm_U_WIRE, offset_VDC_U_WIRE_subtraction);
        
        Logic_VDC_U_GUARDWIRE_Thick[k] = new G4LogicalVolume(Solid_VDC_U_GUARDWIRE_Thick, G4_W_Material, "Logic_VDC_U_GUARDWIRE_Thick",0,0,0);
    }
    
    
    ////////////////////////////////////////////////
    //             VDC INITIALISATION             //
    ////////////////////////////////////////////////
    
    G4ThreeVector   offset_VDC_SenseRegion_USDS[3];
    offset_VDC_SenseRegion_USDS[0] = G4ThreeVector(0.0*mm,0.0*mm, 1.5*thickness_Stesalit_Frame*mm);
    offset_VDC_SenseRegion_USDS[1] = G4ThreeVector(0.0*mm,0.0*mm, -1.5*thickness_Stesalit_Frame*mm);
    
    G4ThreeVector   offset_VDC_HVFrame[3];
    offset_VDC_HVFrame[0] = G4ThreeVector(0.0*mm,0.0*mm, thickness_Stesalit_Frame*mm);
    offset_VDC_HVFrame[1] = G4ThreeVector(0.0*mm,0.0*mm, thickness_Stesalit_Frame*mm);
    offset_VDC_HVFrame[2] = G4ThreeVector(0.0*mm,0.0*mm, -thickness_Stesalit_Frame*mm);
    VDC_HVFrame_rotm[2].rotateY(180.*deg);
    
    G4ThreeVector   offset_VDC_HV_Plane[3];
    offset_VDC_HV_Plane[0] = G4ThreeVector(-1.5*mm, 0.0*mm, (4000.)*um);
    offset_VDC_HV_Plane[1] = G4ThreeVector(0.0*mm, 0.0*mm, (12000. - 10.)*um);
    offset_VDC_HV_Plane[2] = G4ThreeVector(0.0*mm, 0.0*mm, -(3990.)*um);
    
    G4ThreeVector   offset_VDC_GasFrame;
    offset_VDC_GasFrame = G4ThreeVector(0.0*mm,0.0*mm, 0.0*mm);
    
    G4ThreeVector   offset_VDC_XU_Frame[2];
    offset_VDC_XU_Frame[0] = G4ThreeVector(0.0*mm, 0.0*mm, (-1.5*thickness_Stesalit_Frame + 5.5/2)*mm);
    offset_VDC_XU_Frame[1] = G4ThreeVector(0.0*mm, 0.0*mm, -(2.5/2)*mm);
    
    G4ThreeVector   offset_VDC_XU_PCBFrame[2];
    offset_VDC_XU_PCBFrame[0] = G4ThreeVector(0.0*mm, 0.0*mm, (-thickness_Stesalit_Frame/2 - 2.5/2)*mm);
    offset_VDC_XU_PCBFrame[1] = G4ThreeVector(0.0*mm, 0.0*mm, (thickness_Stesalit_Frame/2 - 2.5/2)*mm);
    
    G4ThreeVector   offset_VDC_Al_Frame[2];
    offset_VDC_Al_Frame[0] = G4ThreeVector(0.0*mm, 0.0*mm, (3.*thickness_Stesalit_Frame  + 7.5)*mm);
    offset_VDC_Al_Frame[1] = G4ThreeVector(0.0*mm, 0.0*mm, -(3.*thickness_Stesalit_Frame + 7.5)*mm);
    VDC_Al_Frame_rotm[1].rotateY(180.*deg);
    
    G4ThreeVector   offset_VDC_MYLAR_Plane[2];
    offset_VDC_MYLAR_Plane[0] = G4ThreeVector(0.0*mm, 0.0*mm, (24000. + 25./2)*um);
    offset_VDC_MYLAR_Plane[1] = G4ThreeVector(0.0*mm, 0.0*mm, -(24000. + 25./2)*um);
    
    G4RotationMatrix* rm_GAS_VDC_HV = new G4RotationMatrix();
    rm_GAS_VDC_HV->rotateY(180.*deg);
    
    G4int usds;
    
    G4ThreeVector   offset_VDC_X_WIRE;
    VDC_X_WIRE_rotm.rotateX(90.*deg);
    
    G4ThreeVector   offset_VDC_U_WIRE;
    VDC_U_WIRE_rotm.rotateX(-90.*deg);
    VDC_U_WIRE_rotm.rotateY(-90.*deg);
    VDC_U_WIRE_rotm.rotateZ(40.*deg);
    
    
    for(G4int i=0; i<2; i++)
    {
        if(VDC_Presence[i])
        {
            ///////////////////////////////////////////
            //      VDC - STESALIT STANDARD FRAME
            for(G4int j=0; j<3; j++)
            {
                VDC_HVFrame_transform = G4Transform3D(VDC_HVFrame_rotm[j], offset_VDC_HVFrame[j]);
                
                if(j==0) usds = 0;
                if(j==1 || j ==2) usds = 1;
                
                
                new G4PVPlacement(VDC_HVFrame_transform,
                                  Logic_VDC_HV_Frame[j],
                                  "VDC_HV_Frame",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                new G4PVPlacement(0,
                                  offset_VDC_HV_Plane[j],
                                  Logic_VDC_HV_Plane[j],
                                  "VDC_HV_Plane",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
            }
            
            ///////////////////////////////////////
            //      VDC - GAS CELL
            VDC_GasFrame_transform = G4Transform3D(VDC_GasFrame_rotm, offset_VDC_GasFrame);
            
            new G4PVPlacement(VDC_GasFrame_transform,
                              Logic_VDC_GasFrame,
                              "VDC_GasFrame",
                              Logic_VDC_SenseRegion_USDS[i][0],
                              false,    // no boolean operations
                              0,    // copy number
                              fCheckOverlaps); // checking overlaps
            
            //////////////////////////////////////////////////////
            //     PLACING OBJECTS within each SENSE REGION
            for(G4int j=0; j<2; j++)
            {
                if(j==0) usds = 0;
                if(j==1) usds = 1;
                
                if(j==1)
                {
                    new G4PVPlacement(0,    // no rotation
                                      offset_VDC_XU_PCBFrame[j],
                                      Logic_VDC_XU_PCBFrame,
                                      "VDC_XU_PCBFrame",
                                      Logic_VDC_SenseRegion_USDS[i][usds],
                                      false,    // no boolean operations
                                      i*2 + j,    // copy number
                                      fCheckOverlaps); // checking overlaps
                    
                }
                
                /////////////////////////////////////////////////
                //      VDC - STESALIT X-WIRE/U-WIRE FRAME
                new G4PVPlacement(0,    // no rotation
                                  offset_VDC_XU_Frame[j],
                                  Logic_VDC_XU_Frame,
                                  "VDC_XU_Frame",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                ///////////////////////////////////////////
                //      VDC - PCB X-WIRE/U-WIRE FRAME
                new G4PVPlacement(0,    // no rotation
                                  offset_VDC_XU_PCBFrame[j],
                                  Logic_VDC_XU_PCBFrame,
                                  "VDC_XU_PCBFrame",
                                  Logic_VDC_SenseRegion_USDS[i][usds],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                if(j==0)
                {
                    //////////////////////////////////////////////
                    //      VDC - U WIRES
                    for(G4int k=0; k<143; k++)
                    {
                        //offset_VDC_U_WIRE = G4ThreeVector((k - 71)*4*tan(40*deg)*mm, 0.0*mm, (4000. - 20./2)*um);
                        offset_VDC_U_WIRE = G4ThreeVector((k - 71)*(4/sin(50.*deg))*mm, 0.0*mm, (-4000. - 20./2)*um);
                        
                        VDC_U_WIRE_transform = G4Transform3D(VDC_U_WIRE_rotm, offset_VDC_U_WIRE);
                        
                        PhysiVDC_U_WIRE = new G4PVPlacement(VDC_U_WIRE_transform,
                                                            Logic_VDC_U_WIRE[k],
                                                            "VDC_U_WIRE",
                                                            Logic_VDC_SenseRegion_USDS[i][j],
                                                            false,    // no boolean operations
                                                            i*143 + k,    // copy number
                                                            fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - U GUARD WIRES
                    for(G4int k=0; k<144; k++)
                    {
                        offset_VDC_U_WIRE = G4ThreeVector((k - 71.5)*(4/sin(50.*deg))*mm, 0.0*mm, (-4000. - 20./2)*um);
                        VDC_U_WIRE_transform = G4Transform3D(VDC_U_WIRE_rotm, offset_VDC_U_WIRE);
                        
                        new G4PVPlacement(VDC_U_WIRE_transform,
                                          Logic_VDC_U_GUARDWIRE[k],
                                          "VDC_U_GUARDWIRE",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*144 + k,    // copy number
                                          fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - U GUARD WIRES, Thick
                    for(G4int k=0; k<2; k++)
                    {
                        offset_VDC_U_WIRE = G4ThreeVector((k*144. - 72.)*(4/sin(50.*deg))*mm, 0.0*mm, (-4000. - 20./2)*um);
                        VDC_U_WIRE_transform = G4Transform3D(VDC_U_WIRE_rotm, offset_VDC_U_WIRE);
                        
                        new G4PVPlacement(VDC_U_WIRE_transform,
                                          Logic_VDC_U_GUARDWIRE_Thick[k],
                                          "VDC_U_GUARDWIRE_Thick",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*2 + k,    // copy number
                                          fCheckOverlaps);
                        
                    }
                }
                
                if(j==1)
                {
                    //////////////////////////////////////////////
                    //      VDC - X WIRES
                    for(G4int k=0; k<198; k++)
                    {
                        offset_VDC_X_WIRE = G4ThreeVector(( (k*4.) - 394.)*mm, 0.0*mm, (4000. - 20./2)*um);
                        VDC_X_WIRE_transform = G4Transform3D(VDC_X_WIRE_rotm, offset_VDC_X_WIRE);
                        
                        PhysiVDC_X_WIRE = new G4PVPlacement(VDC_X_WIRE_transform,
                                                            Logic_VDC_X_WIRE,
                                                            "VDC_X_WIRE",
                                                            Logic_VDC_SenseRegion_USDS[i][j],
                                                            false,    // no boolean operations
                                                            i*198 + k,    // copy number
                                                            fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - X GUARD WIRES
                    for(G4int k=0; k<199; k++)
                    {
                        offset_VDC_X_WIRE = G4ThreeVector(( (k*4.) - 396.)*mm, 0.0*mm, (4000. - 20./2)*um);
                        VDC_X_WIRE_transform = G4Transform3D(VDC_X_WIRE_rotm, offset_VDC_X_WIRE);
                        
                        new G4PVPlacement(VDC_X_WIRE_transform,
                                          Logic_VDC_X_GUARDWIRE,
                                          "VDC_X_GUARDWIRE",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*199 + k,    // copy number
                                          fCheckOverlaps); // checking overlaps
                        
                    }
                    
                    //////////////////////////////////////////////
                    //      VDC - X GUARD WIRES, Thick
                    for(G4int k=0; k<2; k++)
                    {
                        offset_VDC_X_WIRE = G4ThreeVector(( (k*2.*398.) - 398.)*mm, 0.0*mm, (4000. - 20./2)*um);
                        VDC_X_WIRE_transform = G4Transform3D(VDC_X_WIRE_rotm, offset_VDC_X_WIRE);
                        
                        new G4PVPlacement(VDC_X_WIRE_transform,
                                          Logic_VDC_X_GUARDWIRE_Thick,
                                          "VDC_X_GUARDWIRE_Thick",
                                          Logic_VDC_SenseRegion_USDS[i][j],
                                          false,    // no boolean operations
                                          i*2 + k,    // copy number
                                          fCheckOverlaps);
                        
                    }
                }
            }
        }
    }
    
    
    for(G4int i=0; i<2; i++)
    {
        if(VDC_Presence[i])
        {
            VDC_transform[i] = G4Transform3D(VDC_rotm[i],VDC_CentrePosition[i]);
            
            for(G4int j=0; j<2; j++)
            {
                
                ////////////////////////////////////////////////////
                //      VDC - US and DS components of Asm
                
                Physical_VDC_SenseRegion_USDS = new G4PVPlacement(0,
                                                                  offset_VDC_SenseRegion_USDS[j],
                                                                  Logic_VDC_SenseRegion_USDS[i][j],
                                                                  "VDC_SenseRegion_USDS",
                                                                  Logic_VDC_Asm[i],
                                                                  false,           // no boolean
                                                                  i*2 + j,               // copy number
                                                                  fCheckOverlaps); // checking overlaps
                
                ////////////////////////////////////////
                //      VDC - ALUMINIUM FRAME
                VDC_Al_Frame_transform = G4Transform3D(VDC_Al_Frame_rotm[j], offset_VDC_Al_Frame[j]);
                
                new G4PVPlacement(VDC_Al_Frame_transform,
                                  Logic_VDC_Al_Frame,
                                  "VDC_Al_Frame",
                                  Logic_VDC_Asm[i],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
                ////////////////////////////////////////
                //      VDC - MYLAR WINDOW
                
                new G4PVPlacement(0,
                                  offset_VDC_MYLAR_Plane[j],
                                  Logic_VDC_MYLAR_Plane,
                                  "VDC_MYLAR_Plane",
                                  Logic_VDC_Asm[i],
                                  false,    // no boolean operations
                                  i*2 + j,    // copy number
                                  fCheckOverlaps); // checking overlaps
                
            }
            
            ////////////////////////////////
            //      VDC Asm
            
            Physical_VDC_Asm = new G4PVPlacement(VDC_transform[i], // transformation
                                                      Logic_VDC_Asm[i],
                                                      "VDC_Asm",
                                                      LogicWorld,
                                                      false,           // no boolean
                                                      i,               // copy number
                                                      fCheckOverlaps); // checking overlaps
            
        }
    }
    
    
    
    
    
    
    
    //////////////////////////////////
    //      PADDLE DEFINITION       //
    //////////////////////////////////
    
    G4Box* Solid_PADDLE_1 = new G4Box("Solid_PADDLE_1", (122/2)*cm, (10.2/2)*cm, (0.3175/2)*cm);
    G4Box* Solid_PADDLE_2 = new G4Box("Solid_PADDLE_2", (122/2)*cm, (10.2/2)*cm, (0.6350/2)*cm);
    G4Box* Solid_PADDLE_3 = new G4Box("Solid_PADDLE_3", (122/2)*cm, (10.2/2)*cm, (1.2700/2)*cm);
    
    G4LogicalVolume * Logic_PADDLE[numberOf_PADDLE];
    
    Logic_PADDLE[0] = new G4LogicalVolume(Solid_PADDLE_1, BC408_Material,"PADDLE_1",0,0,0);
    Logic_PADDLE[1] = new G4LogicalVolume(Solid_PADDLE_2, BC408_Material,"PADDLE_2",0,0,0);
    Logic_PADDLE[2] = new G4LogicalVolume(Solid_PADDLE_3, BC408_Material,"PADDLE_3",0,0,0);
    
    
    //////////////////////////////////////
    //      PADDLE INITIALISATION       //
    //////////////////////////////////////
    
    for(G4int i=0; i<3; i++)
    {
        if(PADDLE_Presence[i])
        {
            PADDLE_transform[i] = G4Transform3D(PADDLE_rotm[i],PADDLE_CentrePosition[i]);
            
            PhysiPADDLE = new G4PVPlacement(PADDLE_transform[i],
                                            Logic_PADDLE[i],    // its logical volume
                                            "PADDLE",           // its name
                                            LogicWorld,         // its mother  volume
                                            false,              // no boolean operations
                                            i,                  // copy number
                                            fCheckOverlaps);    // checking overlaps
            
        }
    }
    
    
    
    //////////////////////////////////
    //      HAGAR DEFINITION        //
    //////////////////////////////////
    
    ////////////////////////////////
    //      HAGAR - NaI Crystal
    G4Tubs* Solid_HAGAR_NaICrystal = new G4Tubs("HAGAR_NaICrystal", 0.*cm, (23.8/2)*cm, (35.6/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_NaICrystal = new G4LogicalVolume(Solid_HAGAR_NaICrystal, G4_SODIUM_IODIDE_Material, "Logic_HAGAR_NaICrystal", 0, 0, 0);
    
    ////////////////////////////////
    //      HAGAR - Annulus
    G4Tubs* Solid_HAGAR_Anulus = new G4Tubs("HAGAR_Anulus", (28.86/2)*cm, ((28.86/2) + 9.84)*cm, (61/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_Annulus = new G4LogicalVolume(Solid_HAGAR_Anulus, BC408_Material, "Logic_HAGAR_Annulus", 0, 0, 0);
    
    ////////////////////////////////
    //      HAGAR - Front Disc
    G4Tubs* Solid_HAGAR_FrontDisc = new G4Tubs("HAGAR_FrontDisc", 0.*cm, (48.58/2)*cm, (8/2)*cm, 0.*deg, 360.*deg);
    
    G4LogicalVolume* Logic_HAGAR_FrontDisc = new G4LogicalVolume(Solid_HAGAR_FrontDisc, BC408_Material, "Logic_HAGAR_FrontDisc", 0, 0, 0);
    
    
    
    ////////////////////////////////////////////////////
    //               HAGAR INITIALIZATION             //
    ////////////////////////////////////////////////////
    
    
    ///////////////////////////////
    //      HAGAR - NaI Crystal
    if(HAGAR_NaICrystal_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_NaICrystal_CentrePosition);
        
        PhysiHAGAR_NaICrystal = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                                  Logic_HAGAR_NaICrystal,       // its logical volume
                                                  "HAGAR_NaICrystal",       // its name
                                                  LogicWorld,         // its mother  volume
                                                  false,           // no boolean operations
                                                  0,               // copy number
                                                  fCheckOverlaps); // checking overlaps
    }
    
    /////////////////////////////
    //      HAGAR - Annulus
    if(HAGAR_Annulus_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_Annulus_CentrePosition);
        
        PhysiHAGAR_Annulus = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                               Logic_HAGAR_Annulus,       // its logical volume
                                               "HAGAR_Annulus",       // its name
                                               LogicWorld,         // its mother  volume
                                               false,           // no boolean operations
                                               0,               // copy number
                                               fCheckOverlaps); // checking overlaps
    }
    
    ///////////////////////////////
    //      HAGAR - Front Disc
    if(HAGAR_FrontDisc_Presence)
    {
        HAGAR_transform = G4Transform3D(HAGAR_rotm, HAGAR_FrontDisc_CentrePosition);
        
        PhysiHAGAR_FrontDisc = new G4PVPlacement(HAGAR_transform,   // transformation matrix
                                                 Logic_HAGAR_FrontDisc,       // its logical volume
                                                 "HAGAR_FrontDisc",       // its name
                                                 LogicWorld,         // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps
    }
    
    
    
    
    
    
    
    /////////////////////////////////////////////
    //             CLOVER DEFINITION           //
    /////////////////////////////////////////////
    
    //G4double CLOVERtoShield_displacement = 10.0;  // cm

    //  Minimum closest distance the CLOVER can placed away from the HEAVIMET to not clash (with CloverEncasement.ply)
    //G4double CLOVERtoShield_displacement = 7.8;  // cm
    
    //  Measured closest distance the CLOVER can placed away from the HEAVIMET (to be used with CloverEncasement_approx.ply
   // G4double CLOVERtoShield_displacement = 7.3;  // cm
    G4double CLOVERtoShield_displacement = 0;  // cm

    G4ThreeVector offset_CLOVERInternalVacuum = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVEREncasement = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal1 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal2 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal3 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    G4ThreeVector offset_CLOVERHPGeCrystal4 = G4ThreeVector(0*cm, 0*cm, -CLOVERtoShield_displacement*cm);
    
    G4LogicalVolume * Logic_CLOVER_InternalVacuum[numberOf_CLOVER];
    G4LogicalVolume * Logic_CLOVER_Encasement;
    G4LogicalVolume * Logic_CLOVER_HPGeCrystal[4];
    G4LogicalVolume * Logic_CLOVER_HPGeCrystal_LithiumContact[4];
    
    bool useCLOVER = false;
    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        if(CLOVER_Presence[i])
        {
            useCLOVER = true;
        }
    }
    
    if(useCLOVER)
    {
        //////////////////////////////////////////////////////////
        //              CLOVER Internal Vacuum - CADMesh
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/CLOVER-InternalVacuum/CloverInternalVacuum.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/CLOVER-InternalVacuum/CloverInternalVacuum_approx.ply");

        CADMesh * mesh_CLOVERInternalVacuum = new CADMesh(meshPath, meshType, mm, offset_CLOVERInternalVacuum, false);
        
        G4VSolid * Solid_CLOVERInternalVacuum = mesh_CLOVERInternalVacuum->TessellatedMesh();
        
        G4VisAttributes* CLOVER_InternalVacuum_VisAtt = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7));
        CLOVER_InternalVacuum_VisAtt->SetVisibility(false);

        for(G4int i=0; i<numberOf_CLOVER; i++)
        {
            Logic_CLOVER_InternalVacuum[i] = new G4LogicalVolume(Solid_CLOVERInternalVacuum, G4_Galactic_Material, "LogicCLOVERInternalVacuum", 0, 0, 0);
            Logic_CLOVER_InternalVacuum[i]->SetVisAttributes(CLOVER_InternalVacuum_VisAtt);
        }
        
        
        ///////////////////////////////////////////////////////
        //              CLOVER Encasement - CADMesh
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Clover-Encasement/CloverEncasement.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Clover-Encasement/CloverEncasement_new_approx.ply");

        CADMesh * mesh_CLOVEREncasement = new CADMesh(meshPath, meshType, mm, offset_CLOVEREncasement, false);
        
        G4VSolid * Solid_CLOVEREncasement = mesh_CLOVEREncasement->TessellatedMesh();
        
        Logic_CLOVER_Encasement = new G4LogicalVolume(Solid_CLOVEREncasement, G4_Al_Material, "LogicCLOVERCloverEncasement", 0, 0, 0);
        
        G4VisAttributes* CLOVER_Encasement_VisAtt = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7));
        //CLOVER_Encasement_VisAtt->SetForceSolid(true);
        Logic_CLOVER_Encasement->SetVisAttributes(CLOVER_Encasement_VisAtt);

        
        //////////////////////////////////////////////////////////
        //              CLOVER HPGeCrystals - CADMesh
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal1.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-RoundedCrystal1_10um.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-RoundedCrystal1_10um_invertedNormals.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_RoundedCrystal1_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal1 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal1, false);
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal2.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-RoundedCrystal2_10um.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_RoundedCrystal2_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal2 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal2, false);
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal3.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-RoundedCrystal3_10um.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_RoundedCrystal3_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal3 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal3, false);
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-Crystal4.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe-RoundedCrystal4_10um.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_RoundedCrystal4_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal4 = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal4, false);
        
        G4VSolid * Solid_HPGeCrystal1 = mesh_CLOVERHPGeCrystal1->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal2 = mesh_CLOVERHPGeCrystal2->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal3 = mesh_CLOVERHPGeCrystal3->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal4 = mesh_CLOVERHPGeCrystal4->TessellatedMesh();
        
        Logic_CLOVER_HPGeCrystal[0] = new G4LogicalVolume(Solid_HPGeCrystal1, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[1] = new G4LogicalVolume(Solid_HPGeCrystal2, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[2] = new G4LogicalVolume(Solid_HPGeCrystal3, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        Logic_CLOVER_HPGeCrystal[3] = new G4LogicalVolume(Solid_HPGeCrystal4, G4_Ge_Material,"LogicCLOVERHPGeCrystal",0,0,0);
        
        G4VisAttributes* CLOVER_HPGeCrystals_VisAtt = new G4VisAttributes(G4Colour(0.9, 0.9, 0.0));
        //CLOVER_HPGeCrystals_VisAtt->SetForceSolid(true);

        Logic_CLOVER_HPGeCrystal[0]->SetVisAttributes(CLOVER_HPGeCrystals_VisAtt);
        Logic_CLOVER_HPGeCrystal[1]->SetVisAttributes(CLOVER_HPGeCrystals_VisAtt);
        Logic_CLOVER_HPGeCrystal[2]->SetVisAttributes(CLOVER_HPGeCrystals_VisAtt);
        Logic_CLOVER_HPGeCrystal[3]->SetVisAttributes(CLOVER_HPGeCrystals_VisAtt);
        
        //////////////////////////////////////////////////////////////////////
        //              CLOVER HPGeCrystals - Lithium contacts - CADMesh

        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_LithiumContact1_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal1_LithiumContact = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal1, false);

        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_LithiumContact2_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal2_LithiumContact = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal2, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_LithiumContact3_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal3_LithiumContact = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal3, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/HPGe-Crystals/HPGe_pureCylindricalBorehole_LithiumContact4_10um.ply");
        CADMesh * mesh_CLOVERHPGeCrystal4_LithiumContact = new CADMesh(meshPath, meshType, mm, offset_CLOVERHPGeCrystal4, false);

        G4VSolid * Solid_HPGeCrystal1_LithiumContact = mesh_CLOVERHPGeCrystal1_LithiumContact->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal2_LithiumContact = mesh_CLOVERHPGeCrystal2_LithiumContact->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal3_LithiumContact = mesh_CLOVERHPGeCrystal3_LithiumContact->TessellatedMesh();
        G4VSolid * Solid_HPGeCrystal4_LithiumContact = mesh_CLOVERHPGeCrystal4_LithiumContact->TessellatedMesh();

        Logic_CLOVER_HPGeCrystal_LithiumContact[0] = new G4LogicalVolume(Solid_HPGeCrystal1_LithiumContact, G4_Li_Material,"LogicCLOVERHPGeCrystal_LithiumContact",0,0,0);
        Logic_CLOVER_HPGeCrystal_LithiumContact[1] = new G4LogicalVolume(Solid_HPGeCrystal2_LithiumContact, G4_Li_Material,"LogicCLOVERHPGeCrystal_LithiumContact",0,0,0);
        Logic_CLOVER_HPGeCrystal_LithiumContact[2] = new G4LogicalVolume(Solid_HPGeCrystal3_LithiumContact, G4_Li_Material,"LogicCLOVERHPGeCrystal_LithiumContact",0,0,0);
        Logic_CLOVER_HPGeCrystal_LithiumContact[3] = new G4LogicalVolume(Solid_HPGeCrystal4_LithiumContact, G4_Li_Material,"LogicCLOVERHPGeCrystal_LithiumContact",0,0,0);

        
        
        //------------------------------------
        //      Lithium doped layer
        for(G4int j=0; j<4; j++)
        {
            //------------------------------------------------
            G4Tubs* outerCylinder = new G4Tubs("outerCylinder", 4.0*mm, 5.5*mm, 0.5*55.5*mm, 0.*deg, 360.*deg);
            G4LogicalVolume *Logic_CLOVER_LithiumDeadLayer_cylinder = new G4LogicalVolume(outerCylinder, G4_Li_Material, "Logic_CLOVER_LithiumDeadLayer_cylinder");
            
            G4Tubs* flatCylinderCap = new G4Tubs("flatCylinderCap", 0.0*mm, 5.1*mm, 0.5*0.5*mm, 0.*deg, 360.*deg);
            G4LogicalVolume *Logic_CLOVER_LithiumDeadLayer_flatCap = new G4LogicalVolume(flatCylinderCap, G4_Li_Material, "Logic_CLOVER_LithiumDeadLayer_flatCap");
            
            G4VisAttributes* CLOVER_DeadLayer_Lithium_VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
            CLOVER_DeadLayer_Lithium_VisAtt->SetForceSolid(true);
            Logic_CLOVER_LithiumDeadLayer_cylinder->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);
            Logic_CLOVER_LithiumDeadLayer_flatCap->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);
            Logic_CLOVER_HPGeCrystal_LithiumContact[j]->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);
            
            //------------------------------------------------
            G4ThreeVector position_LithiumDeadLayer_cylinder;
            
            if(j==0)
            {
                position_LithiumDeadLayer_cylinder = G4ThreeVector(-20.5*mm, 20.5*mm, -(20.0+70.0)*mm - (CLOVERtoShield_displacement*10.0)*mm + 0.5*55.5*mm);
            }
            else if(j==1)
            {
                position_LithiumDeadLayer_cylinder = G4ThreeVector(-20.5*mm, -20.5*mm, -(20.0+70.0)*mm - (CLOVERtoShield_displacement*10.0)*mm + 0.5*55.5*mm);
            }
            else if(j==2)
            {
                position_LithiumDeadLayer_cylinder = G4ThreeVector(20.5*mm, -20.5*mm, -(20.0+70.0)*mm - (CLOVERtoShield_displacement*10.0)*mm + 0.5*55.5*mm);
            }
            else if(j==3)
            {
                position_LithiumDeadLayer_cylinder = G4ThreeVector(20.5*mm, 20.5*mm, -(20.0+70.0)*mm - (CLOVERtoShield_displacement*10.0)*mm + 0.5*55.5*mm);
            }
            
            //------------------------------------------------
            G4RotationMatrix* rm = new G4RotationMatrix();
            G4VSolid* CLOVER_LithiumDeadLayer_cylinderAndFlatCap = new G4UnionSolid("CLOVER_LithiumDeadLayer_cylinderAndFlatCap", outerCylinder, flatCylinderCap, rm, G4ThreeVector(0.,0.,0.5*55.5*mm - 0.5*0.5*mm));
            G4LogicalVolume *Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap = new G4LogicalVolume(CLOVER_LithiumDeadLayer_cylinderAndFlatCap, G4_Li_Material, "Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap");
            Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);
            
            //------------------------------------------------
            G4VSolid* CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge;
            
            if(j==0)
            {
                CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge = new G4IntersectionSolid("CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge", CLOVER_LithiumDeadLayer_cylinderAndFlatCap, Solid_HPGeCrystal1, rm, -position_LithiumDeadLayer_cylinder);
            }
            else if(j==1)
            {
                CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge = new G4IntersectionSolid("CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge", CLOVER_LithiumDeadLayer_cylinderAndFlatCap, Solid_HPGeCrystal2, rm, -position_LithiumDeadLayer_cylinder);
            }
            else if(j==2)
            {
                CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge = new G4IntersectionSolid("CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge", CLOVER_LithiumDeadLayer_cylinderAndFlatCap, Solid_HPGeCrystal3, rm, -position_LithiumDeadLayer_cylinder);
            }
            else if(j==3)
            {
                CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge = new G4IntersectionSolid("CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge", CLOVER_LithiumDeadLayer_cylinderAndFlatCap, Solid_HPGeCrystal4, rm, -position_LithiumDeadLayer_cylinder);
            }
            
            G4LogicalVolume *Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge = new G4LogicalVolume(CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge, G4_Li_Material, "Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge");
            Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);

            //------------------------------------------------
            /*
            static int test = 0;
            
            if(test==0)
            {
                new G4PVPlacement(0,               // no rotation
                                  G4ThreeVector(), // at (x,y,z)
                                  Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge,
                                  "CLOVER_HPGeLithiumDopedDeadlayer", // its name
                                  LogicWorld,
                                  false,           // no boolean operations
                                  0,               // copy number
                                  fCheckOverlaps); // checking overlaps
            }
            
            test++;
            */
            
            //----------------------------------------------------------------------------
            //      Placement for lithium contact created through boolean operations
            /*
             new G4PVPlacement(0,               // no rotation
             position_LithiumDeadLayer_cylinder, // at (x,y,z)
             Logic_CLOVER_LithiumDeadLayer_cylinderAndFlatCap_roundedEdge,
             "CLOVER_HPGeLithiumDopedDeadlayer", // its name
             Logic_CLOVER_HPGeCrystal[j],
             false,           // no boolean operations
             0,               // copy number
             fCheckOverlaps); // checking overlaps
             */
            
            //----------------------------------------------------------------
            //      Placement for lithium contact created through CADMesh
            new G4PVPlacement(0,               // no rotation
                              G4ThreeVector(), // at (x,y,z)
                              Logic_CLOVER_HPGeCrystal_LithiumContact[j],
                              "CLOVER_HPGeLithiumDopedDeadlayer", // its name
                              Logic_CLOVER_HPGeCrystal[j],
                              false,           // no boolean operations
                              0,               // copy number
                              fCheckOverlaps); // checking overlaps
        }
    }
    
    
    ////////////////////////////////////////////
    //              CLOVER SHIELD
    ////////////////////////////////////////////
    
    G4ThreeVector offset_CLOVER_Shield_Body = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_Heavimet = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_BGOCrystals = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_PMT = G4ThreeVector(0*cm, 0*cm, 0*cm);
    G4ThreeVector offset_CLOVER_Shield_PMTConArray = G4ThreeVector(0*cm, 0*cm, 0*cm);
    
    G4LogicalVolume* Logic_CLOVER_Shield_Body;
    G4LogicalVolume* Logic_CLOVER_Shield_Heavimet;
    G4LogicalVolume* Logic_CLOVER_Shield_PMTConArray;
    G4LogicalVolume* Logic_CLOVER_Shield_BGOCrystal[16];
    G4LogicalVolume* Logic_CLOVER_Shield_PMT[16];
    
    bool useCLOVER_Shield = false;
    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        if(CLOVER_Shield_Presence[i])
        {
            useCLOVER_Shield = true;
        }
    }

    if(useCLOVER_Shield)
    {
        ///////////////////////////////////////////////////////
        //              CLOVER Shield Body - CADMesh
        ///////////////////////////////////////////////////////
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body_Modified.ply");
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body_Modified2_tol_10um.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body_Modified3_tol_10um.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Body/Body_Modified4_tol_10um.ply");
        
        CADMesh * mesh_CLOVER_Shield_Body = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_Body, false);
        
        G4VSolid * Solid_CLOVER_Shield_Body = mesh_CLOVER_Shield_Body->TessellatedMesh();
        
        Logic_CLOVER_Shield_Body = new G4LogicalVolume(Solid_CLOVER_Shield_Body, G4_Al_Material, "LogicCLOVERShieldBody", 0, 0, 0);
        
        G4VisAttributes* CLOVER_ShieldBody_VisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6));
        CLOVER_ShieldBody_VisAtt->SetForceSolid(true);
        Logic_CLOVER_Shield_Body->SetVisAttributes(CLOVER_ShieldBody_VisAtt);
        
        ///////////////////////////////////////////////////////
        //              CLOVER Shield Heavimet - CADMesh
        ///////////////////////////////////////////////////////
        
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HeavimetShield.ply");
        
        //--------------------------------------------------------------------------------------------------------------------
        //      The "original" HEAVIMET shield for the AFRODITE array.
        //      The "_Modified" suffix of the .ply filename is because it is has been negligibly altered to avoid geometry overlaps (with a tolerance of 10um)
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HeavimetShield_Modified.ply");
        
        //--------------------------------------------------------------------------------------------------------------------
        //      The "40mm" HEAVIMET shield Paul drew before 05/05/18
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_40mm_mod_10um.ply");
        
        //--------------------------------------------------------------------------------------------------------------------
        //      The new HEAVIMET drawings to optomize the HEAVIMET acceptance with respec to the peak-to-total ratios
        //      Models implemented on 05/05/18 and range from 30mm to 50mm in 2.5mm steps
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_30mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_32_5mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_35mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_37_5mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_40mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_42_5mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_45mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_47_5mm.ply");
        //sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/Heavimet-Shield/HEAVIMET_50mm.ply");
        
        CADMesh * mesh_CLOVER_Shield_Heavimet = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_Heavimet, false);
        
        G4VSolid * Solid_CLOVER_Shield_Heavimet = mesh_CLOVER_Shield_Heavimet->TessellatedMesh();
        
        Logic_CLOVER_Shield_Heavimet = new G4LogicalVolume(Solid_CLOVER_Shield_Heavimet, Heavimet_Material, "LogicCLOVERShieldHeavimet", 0, 0, 0);
        
        G4VisAttributes* CLOVER_Shield_HEAVIMET_VisAtt = new G4VisAttributes(G4Colour(0.4, 0.2, 0.2));
        CLOVER_Shield_HEAVIMET_VisAtt->SetForceSolid(true);
        Logic_CLOVER_Shield_Heavimet->SetVisAttributes(CLOVER_Shield_HEAVIMET_VisAtt);

        ///////////////////////////////////////////////////////
        //      CLOVER Shield PMT Connecter Array - CADMesh
        ///////////////////////////////////////////////////////
        
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMT-Connectors/PMT-ConnecterArray.ply");
        CADMesh * mesh_CLOVER_Shield_PMTConArray = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMTConArray, false);
        
        G4VSolid * Solid_CLOVER_Shield_PMTConArray = mesh_CLOVER_Shield_PMTConArray->TessellatedMesh();
        
        Logic_CLOVER_Shield_PMTConArray = new G4LogicalVolume(Solid_CLOVER_Shield_PMTConArray, G4_Al_Material, "LogicCLOVERShieldHeavimet", 0, 0, 0);
        
        
        ///////////////////////////////////////////////////////
        //      CLOVER Shield BGO Crystals - CADMesh
        ///////////////////////////////////////////////////////
        
        /*
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal1.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal1 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal2.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal2 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal3.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal3 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal4.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal4 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal5.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal5 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal6.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal6 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal7.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal7 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal8.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal8 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal9.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal9 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal10.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal10 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal11.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal11 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal12.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal12 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal13.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal13 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal14.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal14 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal15.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal15 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal16.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal16 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        */
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_1.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal1 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_2.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal2 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_3.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal3 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_4.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal4 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_5.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal5 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_6.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal6 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_7.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal7 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_8.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal8 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_9.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal9 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_10.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal10 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_11.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal11 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_12.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal12 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_13.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal13 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_14.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal14 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_15.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal15 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/BGO-Crystals/BGO-Crystal_Modified_16.ply");
        CADMesh * mesh_CLOVER_Shield_BGOCrystal16 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_BGOCrystals, false);

        
        G4VSolid * Solid_CLOVER_Shield_BGOCrystal[16];
        
        Solid_CLOVER_Shield_BGOCrystal[0] = mesh_CLOVER_Shield_BGOCrystal1->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[1] = mesh_CLOVER_Shield_BGOCrystal2->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[2] = mesh_CLOVER_Shield_BGOCrystal3->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[3] = mesh_CLOVER_Shield_BGOCrystal4->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[4] = mesh_CLOVER_Shield_BGOCrystal5->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[5] = mesh_CLOVER_Shield_BGOCrystal6->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[6] = mesh_CLOVER_Shield_BGOCrystal7->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[7] = mesh_CLOVER_Shield_BGOCrystal8->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[8] = mesh_CLOVER_Shield_BGOCrystal9->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[9] = mesh_CLOVER_Shield_BGOCrystal10->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[10] = mesh_CLOVER_Shield_BGOCrystal11->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[11] = mesh_CLOVER_Shield_BGOCrystal12->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[12] = mesh_CLOVER_Shield_BGOCrystal13->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[13] = mesh_CLOVER_Shield_BGOCrystal14->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[14] = mesh_CLOVER_Shield_BGOCrystal15->TessellatedMesh();
        Solid_CLOVER_Shield_BGOCrystal[15] = mesh_CLOVER_Shield_BGOCrystal16->TessellatedMesh();
        
        
        for(G4int k=0; k<16; k++)
        {
            Logic_CLOVER_Shield_BGOCrystal[k] = new G4LogicalVolume(Solid_CLOVER_Shield_BGOCrystal[k], G4_BGO_Material,"LogicCLOVERShieldBGOCrystal",0,0,0);
        }
        
        ////////////////////////////////////
        //      CLOVER Shield PMT's
        ////////////////////////////////////
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT1.ply");
        CADMesh * mesh_CLOVER_Shield_PMT1 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT2.ply");
        CADMesh * mesh_CLOVER_Shield_PMT2 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT3.ply");
        CADMesh * mesh_CLOVER_Shield_PMT3 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT4.ply");
        CADMesh * mesh_CLOVER_Shield_PMT4 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT5.ply");
        CADMesh * mesh_CLOVER_Shield_PMT5 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT6.ply");
        CADMesh * mesh_CLOVER_Shield_PMT6 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT7.ply");
        CADMesh * mesh_CLOVER_Shield_PMT7 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT8.ply");
        CADMesh * mesh_CLOVER_Shield_PMT8 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT9.ply");
        CADMesh * mesh_CLOVER_Shield_PMT9 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT10.ply");
        CADMesh * mesh_CLOVER_Shield_PMT10 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT11.ply");
        CADMesh * mesh_CLOVER_Shield_PMT11 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT12.ply");
        CADMesh * mesh_CLOVER_Shield_PMT12 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT13.ply");
        CADMesh * mesh_CLOVER_Shield_PMT13 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT14.ply");
        CADMesh * mesh_CLOVER_Shield_PMT14 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT15.ply");
        CADMesh * mesh_CLOVER_Shield_PMT15 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        sprintf(meshPath, "../K600-ALBA/Mesh-Models/DETECTORS/CLOVER/Shield/PMTs/PMT16.ply");
        CADMesh * mesh_CLOVER_Shield_PMT16 = new CADMesh(meshPath, meshType, mm, offset_CLOVER_Shield_PMT, false);
        
        
        G4VSolid * Solid_CLOVER_Shield_PMT[16];
        
        Solid_CLOVER_Shield_PMT[0] = mesh_CLOVER_Shield_PMT1->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[1] = mesh_CLOVER_Shield_PMT2->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[2] = mesh_CLOVER_Shield_PMT3->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[3] = mesh_CLOVER_Shield_PMT4->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[4] = mesh_CLOVER_Shield_PMT5->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[5] = mesh_CLOVER_Shield_PMT6->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[6] = mesh_CLOVER_Shield_PMT7->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[7] = mesh_CLOVER_Shield_PMT8->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[8] = mesh_CLOVER_Shield_PMT9->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[9] = mesh_CLOVER_Shield_PMT10->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[10] = mesh_CLOVER_Shield_PMT11->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[11] = mesh_CLOVER_Shield_PMT12->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[12] = mesh_CLOVER_Shield_PMT13->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[13] = mesh_CLOVER_Shield_PMT14->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[14] = mesh_CLOVER_Shield_PMT15->TessellatedMesh();
        Solid_CLOVER_Shield_PMT[15] = mesh_CLOVER_Shield_PMT16->TessellatedMesh();
        
        for(G4int k=0; k<16; k++)
        {
            Logic_CLOVER_Shield_PMT[k] = new G4LogicalVolume(Solid_CLOVER_Shield_PMT[k], G4_Al_Material,"LogicCLOVERShieldPMT",0,0,0);
        }
    }
    
    
    ////////////////////////////////////////////////////
    //               CLOVER INITIALIZATION            //
    ////////////////////////////////////////////////////
    
    //------------------------------------------------------------------------------------
    //      Offsets for Walid's crystal geometry
    G4double length_head = (5./atan(7.1*deg)); // 5 = total enlevé vers lavant de la face tapered
    //double offset = (-35.0-20.0); // mm
    //double offset = (-35.0-20.0)-(CLOVERtoShield_displacement*10.0); // mm
    
    //  Properly centered Logic_HPGeCrystal_Walid
    //double offset = (-35.0-20.0)-(CLOVERtoShield_displacement*10.0); // mm
    
    //  Walid's mistake
    //double offset = (28.5-length_head-20.0)-20.0-(CLOVERtoShield_displacement*10.0); // mm
    //double offset = (28.5-length_head-20); // mm
    
    //double offset = (28.5-length_head-20.0)-(CLOVERtoShield_displacement*10.0); // mm
    //double offset = (-35.0+length_head-20.0)-(CLOVERtoShield_displacement*10.0); // mm
    //double offset = (28.5-length_head-20+1.5)-(CLOVERtoShield_displacement*10.0); // mm
    
    //------------------------------------------------
    //      New working for Logic_HPGeCrystal_Walid_2
    //      This should be the working, but it is slightly off - likely due to the origin not being defined well (clearly) for the volume and/or approximations
    //double offset = (-length_head*mm-20.0*mm)-(CLOVERtoShield_displacement*10.0)*mm; // mm
    
    //      Modified to match by eye - matched with our crystals
    double walidCrystaloffset = (-length_head*mm-20.0*mm+0.41*mm)-(CLOVERtoShield_displacement*10.0)*mm; // mm
    
    //      Modified to match by eye - to match the front face with the X-Y plane
    //double offset = (-length_head-(1.5/2.0))*mm; // mm
    
    //------------------------------------
    //      Lithium doped layer
    if(useCLOVER_Walid)
    {
        G4Tubs* outerCylinder = new G4Tubs("outerCylinder", 5.0*mm, 5.5*mm, 0.5*55.5*mm, 0.*deg, 360.*deg);
        G4LogicalVolume *Logic_CLOVER_LithiumDeadLayer_cylinder = new G4LogicalVolume(outerCylinder, G4_Li_Material, "Logic_CLOVER_LithiumDeadLayer_cylinder");
        
        G4Tubs* flatCylinderCap = new G4Tubs("flatCylinderCap", 0.0*mm, 5.0*mm, 0.5*0.5*mm, 0.*deg, 360.*deg);
        G4LogicalVolume *Logic_CLOVER_LithiumDeadLayer_flatCap = new G4LogicalVolume(flatCylinderCap, G4_Li_Material, "Logic_CLOVER_LithiumDeadLayer_flatCap");
        
        G4VisAttributes* CLOVER_DeadLayer_Lithium_VisAtt = new G4VisAttributes(G4Colour(0.7, 0.7, 0.0));
        CLOVER_DeadLayer_Lithium_VisAtt->SetForceSolid(true);
        Logic_CLOVER_LithiumDeadLayer_cylinder->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);
        Logic_CLOVER_LithiumDeadLayer_flatCap->SetVisAttributes(CLOVER_DeadLayer_Lithium_VisAtt);

        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(0,0,-(70.0-length_head)*mm + (0.5*55.5)*mm - 0.41*mm), // at (x,y,z)
                          Logic_CLOVER_LithiumDeadLayer_cylinder,
                          "CLOVER_HPGeLithiumDopedDeadlayer", // its name
                          Logic_HPGeCrystal_Walid_2,
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
        
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(0,0,-(70.0-length_head)*mm + (0.5*55.5)*mm - 0.41*mm + (0.5*55.5) - (0.5*0.5)*mm), // at (x,y,z)
                          Logic_CLOVER_LithiumDeadLayer_flatCap,
                          "CLOVER_HPGeLithiumDopedDeadlayer", // its name
                          Logic_HPGeCrystal_Walid_2,
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        
    }

    //--------------------------------------------------------------------------------
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        //--------------------------------------------
        CLOVER_position[i] = (CLOVER_Distance[i])*G4ThreeVector(sin(CLOVER_theta[i]) * cos(CLOVER_phi[i]), sin(CLOVER_theta[i]) * sin(CLOVER_phi[i]), cos(CLOVER_theta[i]));
        CLOVER_Shield_position[i] = (CLOVER_Distance[i])*G4ThreeVector(sin(CLOVER_theta[i]) * cos(CLOVER_phi[i]), sin(CLOVER_theta[i]) * sin(CLOVER_phi[i]), cos(CLOVER_theta[i]));
        
        //--------------------------------------------
        G4ThreeVector positionVector = CLOVER_position[i].unit();
        
        
        G4ThreeVector positionVector_z = positionVector.unit();
        G4ThreeVector positionVector_y = (positionVector.orthogonal()).unit();
        G4ThreeVector positionVector_x = (positionVector_y.cross(positionVector_z)).unit();
        positionVector_z = -positionVector_z;
        positionVector_y = -positionVector_y;
        
        G4RotationMatrix rotmPrime(positionVector_x, positionVector_y, positionVector_z);
        //G4RotationMatrix rotmPrime(-positionVector);
        
        
        //--------------------------------------------
        
        G4ThreeVector initialXAxis(1.0, 0.0, 0.0); //
        G4ThreeVector initialYAxis(0.0, 1.0, 0.0); //
        G4ThreeVector rotatedXAxis = rotmPrime*initialXAxis;
        
        G4ThreeVector requiredFinalXaxis;
        
        //      Sides of the detector to be parallel to the vertical axis
        requiredFinalXaxis = positionVector.cross(initialYAxis).unit();
        //      Detector orientated to be along constant-phi lines
        //requiredFinalXaxis = G4ThreeVector(-sin(CLOVER_phi[i]*deg), cos(CLOVER_phi[i]*deg), 0.0); // Unit vector of phi
       // requiredFinalXaxis = G4ThreeVector(cos(CLOVER_theta[i])*cos(CLOVER_phi[i]), cos(CLOVER_theta[i])*sin(CLOVER_phi[i]), -sin(CLOVER_theta[i])); // unit vector of theta

        G4double requiredRotationAngle = rotatedXAxis.angle(requiredFinalXaxis)/deg;
        
        
        
        if(CLOVER_Presence[i])
        {
            //G4cout << "requiredFinalXaxis.x(): " << requiredFinalXaxis.x() << G4endl;
           // G4cout << "requiredFinalXaxis.y(): " << requiredFinalXaxis.y() << G4endl;
            //G4cout << "requiredFinalXaxis.z(): " << requiredFinalXaxis.z() << G4endl;
            G4cout << "requiredRotationAngle: " << requiredRotationAngle << G4endl;
            G4cout << G4endl;
        }
        
        
        
        G4RotationMatrix rotm_AlongDetectorAxis(positionVector, (requiredRotationAngle)*deg);
        
        //CLOVER_rotm[i] = rotm_AlongDetectorAxis*rotmPrime;
        CLOVER_rotm[i] = rotmPrime.transform(rotm_AlongDetectorAxis);
        
        
        //--------------------------------------------
       // G4ThreeVector requiredFinalXaxis = G4ThreeVector(cos(CLOVER_theta[i])*cos(CLOVER_phi[i]), cos(CLOVER_theta[i])*sin(CLOVER_phi[i]), -sin(CLOVER_theta[i])); // unit vector of theta
        G4ThreeVector requiredFinalZaxis = -positionVector;
        G4ThreeVector requiredFinalYaxis = requiredFinalZaxis.cross(requiredFinalXaxis).unit();

        G4RotationMatrix rotmPrime2;
        rotmPrime2.rotateAxes(requiredFinalXaxis, requiredFinalYaxis, requiredFinalZaxis);
        //rotmPrime2.rotateUz(requiredFinalZaxis);
        
        //G4double requiredRotationAngle = rotatedXAxis.angle(requiredFinalXaxis)/deg;
        
        /*
         G4RotationMatrix rotm_AlongDetectorAxis(-positionVector, (requiredRotationAngle)*deg);
         
         CLOVER_rotm[i] = rotm_AlongDetectorAxis*rotmPrime;
         */
        
        
        //G4RotationMatrix *rotmPrime2 = new G4RotationMatrix;
        //(initialXAxis, initialYAxis, initialZAxis);
        //rotmPrime2->rotateAxes(requiredFinalXaxis, requiredFinalYaxis, requiredFinalZaxis);
        
        //--------------------------------------------
        //CLOVER_rotm[i] = rotmPrime;
        //CLOVER_rotm[i] = rotm_AlongDetectorAxis*rotmPrime;
        //CLOVER_rotm[i] = rotmPrime;
        //CLOVER_rotm[i] = rotmPrime2;
        CLOVER_transform[i] = G4Transform3D(CLOVER_rotm[i],CLOVER_position[i]);
        
        
      //  CLOVER_transform[i] = G4Transform3D(rotmPrime2, CLOVER_position[i]);
        
        CLOVER_Shield_transform[i] = CLOVER_transform[i];

        /////////////////////////////
        //          CLOVER
        if(CLOVER_Presence[i])
        {
            new G4PVPlacement(CLOVER_transform[i],   // transformation matrix
                              Logic_CLOVER_Encasement,       // its logical volume
                              "CLOVER_Encasement",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
        
            if(useCLOVER_Walid)
            {
                G4VisAttributes* CLOVER_HPGeCrystals_Walid_VisAtt = new G4VisAttributes(G4Colour(0.9, 0.9, 0.0));
                CLOVER_HPGeCrystals_Walid_VisAtt->SetForceSolid(true);
                Logic_HPGeCrystal_Walid->SetVisAttributes(CLOVER_HPGeCrystals_Walid_VisAtt);
                //Logic_HPGeCrystal_Walid_2->SetVisAttributes(CLOVER_HPGeCrystals_Walid_VisAtt);
                
                //G4double walid_crystalBackShift = abs(28.5-length_head-20+1.5);
                
                G4cout << "length_head: " << length_head << G4endl;
                //--------------------
                //      STANDARD
                
                bool walidComparison = false;
                G4int nWalidCrystals = 0;
                
                if(walidComparison)
                {
                    nWalidCrystals = 3;
                }
                else
                {
                    nWalidCrystals = 4;
                }
                
                for(int j=0; j<nWalidCrystals; j++)
                {
                    //------------------------------------------------
                    G4ThreeVector position_HPGeCrystal_Walid;
                    
                    if(j==0)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(-20.501*mm, 20.501*mm, walidCrystaloffset*mm);
                    }
                    else if(j==1)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(-20.501*mm, -20.501*mm, walidCrystaloffset*mm);
                    }
                    else if(j==2)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(20.501*mm, -20.501*mm, walidCrystaloffset*mm);
                    }
                    else if(j==3)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(20.501*mm, 20.501*mm, walidCrystaloffset*mm);
                    }
                    
                    //------------------------------------------------
                    G4RotationMatrix* rm_HPGeCrystal_Walid = new G4RotationMatrix();
                    rm_HPGeCrystal_Walid->rotateZ(((j+1)*90.0)*deg);
                    
                    //------------------------------------------------
                    G4Transform3D CLOVER_HPGeCrystal_Walid_transform = G4Transform3D(*rm_HPGeCrystal_Walid, position_HPGeCrystal_Walid);
                    
                    //------------------------------------------------
                    
                    
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(CLOVER_HPGeCrystal_Walid_transform,
                                                                Logic_HPGeCrystal_Walid_2,
                                                                "CLOVER_HPGeCrystal", // its name
                                                                Logic_CLOVER_InternalVacuum[i],
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                    
                    
                    /*
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(CLOVER_HPGeCrystal_Walid_transform,
                                                                Logic_HPGeCrystal_Walid_2,
                                                                "CLOVER_HPGeCrystal", // its name
                                                                LogicWorld,
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                    */
                    
                    
                }

                for(int j=nWalidCrystals; j<4; j++)
                {
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(0,               // no rotation
                                                                G4ThreeVector(0,0,0), // at (x,y,z)
                                                                Logic_CLOVER_HPGeCrystal[j],
                                                                "CLOVER_HPGeCrystal", // its name
                                                                Logic_CLOVER_InternalVacuum[i],
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                    
                }
                
                
                new G4PVPlacement(0,               // no rotation
                                  G4ThreeVector(0,0,-20.0*mm-70.0*mm-10.0*mm-(CLOVERtoShield_displacement*10.0)*mm), // at (x,y,z)
                                  Logic_CLOVER_HPGeAluminiumBackingPlate,
                                  "CLOVER_HPGeAluminiumBackingPlate", // its name
                                  Logic_CLOVER_InternalVacuum[i],
                                  false,           // no boolean operations
                                  i,               // copy number
                                  fCheckOverlaps); // checking overlaps
                
                
                
                /*
                //--------------------------------------------
                //      Test at origin (of VacuumChamber)
                for(int j=0; j<1; j++)
                {
                    //------------------------------------------------
                    G4ThreeVector position_HPGeCrystal_Walid;
                    //double offset = -(28.5-length_head); // mm
                    double offset = (28.5-length_head-20); // mm
                    //double offset = (-35.0-20.0)-(CLOVERtoShield_displacement*10.0); // mm
                    //double offset = (-35.0-20.0)-(CLOVERtoShield_displacement*10.0); // mm
                    
                    if(j==0)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(-20.501*mm, 20.501*mm, offset*mm);
                    }
                    else if(j==1)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(-20.501*mm, -20.501*mm, offset*mm);
                    }
                    else if(j==2)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(20.501*mm, -20.501*mm, offset*mm);
                    }
                    else if(j==3)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(20.501*mm, 20.501*mm, offset*mm);
                    }
                    
                    //------------------------------------------------
                    G4RotationMatrix* rm_HPGeCrystal_Walid = new G4RotationMatrix();
                    rm_HPGeCrystal_Walid->rotateZ(((j+1)*90.0)*deg);
                    
                    //------------------------------------------------
                    G4Transform3D CLOVER_HPGeCrystal_Walid_transform = G4Transform3D(*rm_HPGeCrystal_Walid, position_HPGeCrystal_Walid);
                    
                    //------------------------------------------------
                    
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(CLOVER_HPGeCrystal_Walid_transform,
                                                                Logic_HPGeCrystal_Walid,
                                                                "CLOVER_HPGeCrystal", // its name
                                                                LogicVacuumChamber,
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                    
                    
                    
                    offset = (CLOVERtoShield_displacement*cm+20*mm); // mm
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(0,               // no rotation
                                                                G4ThreeVector(0,0,offset), // at (x,y,z)
                                                                Logic_CLOVER_HPGeCrystal[j],
                                                                "CLOVER_HPGeCrystal", // its name
                                                                Logic_CLOVER_InternalVacuum[i],
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                }
                */
                
                //--------------------
                //      TEMP
                /*
                for(int j=0; j<4; j++)
                {
                    //------------------------------------------------
                    G4ThreeVector position_HPGeCrystal_Walid;
                    
                    if(j==0)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(0.0, 0.0, 0.0);
                    }
                    else if(j==1)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(0.0, 0.0, 0.0);
                    }
                    else if(j==2)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(0.0, 0.0, 0.0);
                    }
                    else if(j==3)
                    {
                        position_HPGeCrystal_Walid = G4ThreeVector(0.0, 0.0, 0.0);
                    }
                    
                    //------------------------------------------------
                    G4RotationMatrix* rm_HPGeCrystal_Walid = new G4RotationMatrix();
                    rm_HPGeCrystal_Walid->rotateZ((i*90.0)*deg);
                    
                    //------------------------------------------------
                    G4Transform3D CLOVER_HPGeCrystal_Walid_transform = G4Transform3D(*rm_HPGeCrystal_Walid, position_HPGeCrystal_Walid);

                    //------------------------------------------------
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(CLOVER_HPGeCrystal_Walid_transform,
                                                                Logic_HPGeCrystal_Walid,
                                                                "CLOVER_HPGeCrystal", // its name
                                                                Logic_CLOVER_InternalVacuum[i],
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                    
                }
                */
            }
            else
            {
                for(int j=0; j<4; j++)
                {
                    PhysiCLOVER_HPGeCrystal = new G4PVPlacement(0,               // no rotation
                                                                G4ThreeVector(0,0,0), // at (x,y,z)
                                                                Logic_CLOVER_HPGeCrystal[j],
                                                                "CLOVER_HPGeCrystal", // its name
                                                                Logic_CLOVER_InternalVacuum[i],
                                                                false,           // no boolean operations
                                                                i*4 + j,               // copy number
                                                                fCheckOverlaps); // checking overlaps
                    
                }
            }
            
            new G4PVPlacement(CLOVER_transform[i],
                              Logic_CLOVER_InternalVacuum[i],
                              "CLOVER_InternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
        }
        
        /*
        //------------------------------------------------
        //      TEST
        new G4PVPlacement(0,               // no rotation
                          G4ThreeVector(0,0,-(5.0)*mm), // at (x,y,z)
                          Logic_CLOVER_HPGeAluminiumBackingPlate,
                          "CLOVER_HPGeCrystal", // its name
                          LogicWorld,
                          false,           // no boolean operations
                          i,               // copy number
                          fCheckOverlaps); // checking overlaps
        */
        
        
        /////////////////////////////
        //      CLOVER SHIELD
        if(CLOVER_Shield_Presence[i])
        {
            char volumeName[512];
            sprintf(volumeName, "CLOVER_Shield_Body_%d", i);
            
            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_Body,
                              volumeName,       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            for(int j=0; j<16; j++)
            {
                //------------------------------------------------
                sprintf(volumeName, "CLOVER_%d_Shield_BGOCrystal%d", i, j);
                
                /*
                PhysiCLOVER_Shield_BGOCrystal = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                                  Logic_CLOVER_Shield_BGOCrystal[j],
                                                                  volumeName,
                                                                  LogicVacuumChamber,
                                                                  false,
                                                                  i*16 + j,
                                                                  fCheckOverlaps);
                */
                
                PhysiCLOVER_Shield_BGOCrystal = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                                  Logic_CLOVER_Shield_BGOCrystal[j],
                                                                  "CLOVER_Shield_BGOCrystal",
                                                                  LogicVacuumChamber,
                                                                  false,
                                                                  i*16 + j,
                                                                  fCheckOverlaps);

                //------------------------------------------------
                sprintf(volumeName, "CLOVER_%d_Shield_PMT%d", i, j);
                
                PhysiCLOVER_Shield_PMT = new G4PVPlacement(CLOVER_Shield_transform[i],
                                                           Logic_CLOVER_Shield_PMT[j],
                                                           "CLOVER_Shield_PMT", // its name
                                                           LogicVacuumChamber,
                                                           false, // no boolean operations
                                                           i*16 + j,  // copy number
                                                           fCheckOverlaps); // checking overlaps
                
            }
            
            /*
            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_PMTConArray,
                              "CLOVER_Shield_PMTConArray",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            */
            
            
            sprintf(volumeName, "CLOVER_Shield_Heavimet_%d", i);

            new G4PVPlacement(CLOVER_Shield_transform[i],
                              Logic_CLOVER_Shield_Heavimet,
                              volumeName,       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
    }
    
    
    ////////////////////////////////////////
    ////        LEPS DEFINITION         ////
    ////////////////////////////////////////
    
    //////////////////////////////////////////////////////////
    //              LEPS Internal Vacuum - CADMesh
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_InternalVacuum = new G4Tubs("Solid_LEPSInternalVacuum", 0.0*mm, 38.4*mm, 45.0*mm, 0.*deg, 360*deg);
    G4LogicalVolume* Logic_LEPS_InternalVacuum[numberOf_LEPS];
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        Logic_LEPS_InternalVacuum[i] = new G4LogicalVolume(Solid_LEPS_InternalVacuum, G4_Galactic_Material, "LogicLEPSInternalVacuum", 0, 0, 0);
    }
    
    ///////////////////////////////////////////////////////
    //              LEPS Encasement - CADMesh
    ///////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_Encasement = new G4Tubs("Solid_LEPSEncasement", 38.5*mm, 40.0*mm, (90./2)*mm, 0.*deg, 360*deg);
    
    G4LogicalVolume* Logic_LEPS_Encasement = new G4LogicalVolume(Solid_LEPS_Encasement, G4_Al_Material, "LogicLEPSLEPSEncasement", 0, 0, 0);
    
    
    ///////////////////////////////////////////////////////
    //              LEPS Beryllium Window - CADMesh
    ///////////////////////////////////////////////////////
    
    G4Tubs* Solid_LEPS_Window = new G4Tubs("Solid_LEPSWindow", 0.0*mm, 38.5*mm, (0.3/2)*mm, 0.*deg, 360*deg);
    
    G4LogicalVolume* Logic_LEPS_Window = new G4LogicalVolume(Solid_LEPS_Window, G4_Be_Material, "Logic_LEPS_Window", 0, 0, 0);
    
    
    //////////////////////////////////////////////////////////
    //              LEPS HPGeCrystals - CADMesh
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_HPGeCrystal = new G4Tubs("Solid_HPGeCrystal1", 0.0*mm, 33.0*mm, 5.5*mm, 0.*deg, 90.*deg);
    
    G4LogicalVolume* Logic_LEPS_HPGeCrystal = new G4LogicalVolume(Solid_HPGeCrystal, G4_Ge_Material,"LogicLEPSHPGeCrystal",0,0,0);
    
    LEPS_HPGeCrystal_rotm[0].rotateZ(0.*deg);
    LEPS_HPGeCrystal_rotm[1].rotateZ(90.*deg);
    LEPS_HPGeCrystal_rotm[2].rotateZ(180.*deg);
    LEPS_HPGeCrystal_rotm[3].rotateZ(270.*deg);
    
    for(G4int i=0; i<4; i++)
    {
        LEPS_HPGeCrystal_transform[i] = G4Transform3D(LEPS_HPGeCrystal_rotm[i], G4ThreeVector(0,0,(29.0-0.5)*mm));
    }
    
    ////////////////////////////////////////////////////
    //               LEPS INITIALIZATION
    ////////////////////////////////////////////////////
    
    
    for(G4int i=0; i<numberOf_LEPS; i++)
    {
        LEPS_position[i] = (LEPS_Distance[i] + 4.5*cm)*G4ThreeVector(std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        
        LEPS_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_position[i]);
        
        LEPS_InternalVacuum_position[i] = (LEPS_Distance[i]+ 4.5*cm + 0.5*mm)*G4ThreeVector(std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        LEPS_InternalVacuum_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_InternalVacuum_position[i]);
        
        LEPS_Window_position[i] = (LEPS_Distance[i] + 4.5*cm +(-45.0+0.15)*mm)*G4ThreeVector(std::sin(LEPS_theta[i]) * std::cos(LEPS_phi[i]), std::sin(LEPS_theta[i]) * std::sin(LEPS_phi[i]), std::cos(LEPS_theta[i]));
        LEPS_Window_transform[i] = G4Transform3D(LEPS_rotm[i],LEPS_Window_position[i]);
        
        /////////////////////////////
        //          LEPS
        if(LEPS_Presence[i] == true)
        {
            
            new G4PVPlacement(LEPS_transform[i],   // transformation matrix
                              Logic_LEPS_Encasement,       // its logical volume
                              "LEPSEncasement",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            new G4PVPlacement(LEPS_Window_transform[i],   // transformation matrix
                              Logic_LEPS_Window,       // its logical volume
                              "LEPSWindow",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
            
            for(int j=0; j<4; j++)
            {
                Physical_LEPS_HPGeCrystal = new G4PVPlacement(LEPS_HPGeCrystal_transform[j],
                                                              Logic_LEPS_HPGeCrystal,       // its logical volume
                                                              "LEPSHPGeCrystal",       // its name
                                                              Logic_LEPS_InternalVacuum[i],    // its mother  volume
                                                              false,           // no boolean operations
                                                              j + (i*4),               // copy number
                                                              fCheckOverlaps); // checking overlaps
                
            }
            
            new G4PVPlacement(LEPS_InternalVacuum_transform[i],
                              Logic_LEPS_InternalVacuum[i],
                              "LEPSInternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            
        }
        
    }
    
    
    ////////////////////////////////////////////
    ////        LaBr3Ce DEFINITION          ////
    ////////////////////////////////////////////

    /*
    double LaBr3Ce_encasement_innerRadius = 45.0; // mm
    double LaBr3Ce_encasement_outerRadius = 50.0; // mm

    double LaBe3CeCrystal_cone_innerRadius;
    double LaBe3CeCrystal_cone_outerRadius;
    
    double LaBr3Ce_encasement_cone_innerRadius;
    double LaBr3Ce_encasement_cone_outerRadius = LaBr3Ce_encasement_outerRadius;
    
    double crystalBase_outerRadius = 45.0; // mm
    double LaBr3Ce_window_axialLength = 5.0; // mm
    */
    
    //----------------------------------------------------
    //      Work for Christiaan/Mathis/Katarzyna
    
    double LaBr3Ce_encasement_innerRadius = 45.1; // mm
    double LaBr3Ce_encasement_outerRadius = 50.0; // mm
    
    double LaBe3CeCrystal_cone_innerRadius;
    double LaBe3CeCrystal_cone_outerRadius;
    
    double LaBr3Ce_encasement_cone_innerRadius;
    double LaBr3Ce_encasement_cone_outerRadius = LaBr3Ce_encasement_outerRadius;
    
    double crystalBase_outerRadius = 45.0; // mm
    double LaBr3Ce_window_axialLength = 5.0; // mm
    
    
    /*
    LaBr3Ce_encasement_innerRadius = ;
    LaBr3Ce_encasement_outerRadius = (89.0/2.0) + 5.0;; // mm
    crystalBase_outerRadius = (89.0/2.0); // mm
    */
    
    //----------------------------------------------------
    //  This sets the tapering
    //----------------------------------------------------
    
    if(LaBr3CeSetupVersion==-1)
    {
        //      No tapering
        LaBe3CeCrystal_cone_outerRadius = crystalBase_outerRadius; // mm
    }
    else if(LaBr3CeSetupVersion==0)
    {
        //      No tapering
        LaBe3CeCrystal_cone_outerRadius = crystalBase_outerRadius; // mm
    }
    else if(LaBr3CeSetupVersion>0)
    {
        double theta = 41.8103/2.0; // deg
        double a, b;
        
        b = LaBR3Ce_GlobalDistance/mm; // mm
        a = b*tan(theta*deg); // mm
        
        LaBe3CeCrystal_cone_outerRadius = a - 3.4; // mm
    }
    
    //----------------------------------------------------
    /*
    if(LaBe3CeCrystal_cone_outerRadius > LaBr3Ce_encasement_innerRadius)
    {
        LaBe3CeCrystal_cone_outerRadius = LaBr3Ce_encasement_innerRadius;
    }
    */


    //----------------------------------------------------

    double LaBr3Ce_window_smallerRadius;
    //double LaBr3Ce_window_largerRadius = (LaBe3CeCrystal_cone_outerRadius+LaBr3Ce_window_axialLength); // mm
    LaBr3Ce_encasement_cone_innerRadius = LaBe3CeCrystal_cone_outerRadius+LaBr3Ce_window_axialLength; // mm

    double LaBr3Ce_window_largerRadius = (LaBr3Ce_encasement_cone_innerRadius); // mm



    //////////////////////////////////////////////////////
    //              LaBr3Ce Crystals
    //////////////////////////////////////////////////////
    double LaBr3Ce_crystal_axialLength = 203.0; // mm
    double LaBr3Ce_crystalCone_axialLength;

    if(LaBr3CeSetupVersion==-1)
    {
        LaBr3Ce_crystalCone_axialLength = 20.0; // mm
    }
    else if(LaBr3CeSetupVersion==0)
    {
        LaBr3Ce_crystalCone_axialLength = 20.0; // mm
    }
    else
    {
        double a, b;
        
        b = crystalBase_outerRadius - LaBe3CeCrystal_cone_outerRadius;
        a = b/(tan((41.8103/2.0)*deg));
        
        LaBr3Ce_crystalCone_axialLength = a;
    }

    bool automaticWindowConeRadii = true;
    if(automaticWindowConeRadii)
    {
        double a, b;
        
        b = LaBr3Ce_window_axialLength;
        a = b*tan((41.8103/2.0)*deg);
        
        LaBr3Ce_window_smallerRadius = LaBr3Ce_window_largerRadius - a;
    }

    //----------------------------------------------------------------
    double LaBr3Ce_crystalBase_axialLength = LaBr3Ce_crystal_axialLength-LaBr3Ce_crystalCone_axialLength; // mm

    G4Tubs* Solid_LaBe3CeCrystal_base = new G4Tubs("Solid_LaBe3CeCrystal_base", 0.0*mm, crystalBase_outerRadius*mm, (LaBr3Ce_crystalBase_axialLength/2.0)*mm, 0.*deg, 360.0*deg);
    G4Cons* Solid_LaBe3CeCrystal_cone = new G4Cons("Solid_LaBe3CeCrystal_cone", 0.0*mm, crystalBase_outerRadius*mm, 0.0*mm, LaBe3CeCrystal_cone_outerRadius*mm, (LaBr3Ce_crystalCone_axialLength/2.0)*mm, 0.*deg, 360.0*deg);
    
    G4RotationMatrix* rm_LaBr3CeCrystal = new G4RotationMatrix();
    G4ThreeVector position_LaBr3CeCrystal(0.0, 0.0, ((LaBr3Ce_crystalBase_axialLength+LaBr3Ce_crystalCone_axialLength)/2.0 - 1.0)*mm);
    
    G4LogicalVolume* Logic_LaBr3Ce_LaBr3CeCrystal_base = new G4LogicalVolume(Solid_LaBe3CeCrystal_base, LaBr3Ce3Ce_Material,"Logic_LaBr3Ce_LaBr3CeCrystal_base",0,0,0);
    G4LogicalVolume* Logic_LaBr3Ce_LaBr3CeCrystal_cone = new G4LogicalVolume(Solid_LaBe3CeCrystal_cone, LaBr3Ce3Ce_Material,"Logic_LaBr3Ce_LaBr3CeCrystal_cone",0,0,0);
    
    G4LogicalVolume* Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalBase = new G4LogicalVolume(Solid_LaBe3CeCrystal_base, LaBr3Ce3Ce_Material,"Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalBase",0,0,0);
    G4LogicalVolume* Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalCone = new G4LogicalVolume(Solid_LaBe3CeCrystal_cone, LaBr3Ce3Ce_Material,"Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalCone",0,0,0);

    double vol1 = Solid_LaBe3CeCrystal_base->GetCubicVolume()/cm3;
    double vol2 = Solid_LaBe3CeCrystal_cone->GetCubicVolume()/cm3;
    double vol_LaBr3Ce = vol1 + vol2;
    
    ///////////////////////////////////////////////////////
    //              LaBr3Ce Aluminium Window
    ///////////////////////////////////////////////////////
    
    if(LaBr3Ce_window_largerRadius > LaBr3Ce_encasement_outerRadius)
    {
        LaBr3Ce_window_largerRadius = LaBr3Ce_encasement_outerRadius;
    }
    
    G4Tubs* Solid_LaBr3Ce_Window = new G4Tubs("Solid_LaBr3CeWindow", 0.0*mm, LaBr3Ce_window_largerRadius*mm, (LaBr3Ce_window_axialLength/2.0)*mm, 0.*deg, 360*deg);
    G4Cons* Solid_LaBr3Ce_Window_cone;
    
    if(LaBr3CeSetupVersion>0)
    {
        Solid_LaBe3CeCrystal_cone = new G4Cons("Solid_LaBe3CeCrystal_cone", 0.0*mm, LaBr3Ce_window_largerRadius*mm, 0.0*mm, LaBr3Ce_window_smallerRadius*mm, (LaBr3Ce_window_axialLength/2.0)*mm, 0.*deg, 360.0*deg);
    }
    
    G4LogicalVolume* Logic_LaBr3Ce_Window;
    if(LaBr3CeSetupVersion==-1)
    {
        Logic_LaBr3Ce_Window = new G4LogicalVolume(Solid_LaBr3Ce_Window, G4_Al_Material, "Logic_LaBr3Ce_Window", 0, 0, 0);
    }
    else if(LaBr3CeSetupVersion==0)
    {
        Logic_LaBr3Ce_Window = new G4LogicalVolume(Solid_LaBr3Ce_Window, G4_Al_Material, "Logic_LaBr3Ce_Window", 0, 0, 0);
    }
    else
    {
        Logic_LaBr3Ce_Window = new G4LogicalVolume(Solid_LaBe3CeCrystal_cone, G4_Al_Material, "Logic_LaBr3Ce_Window", 0, 0, 0);
    }
    
    ///////////////////////////////////////////////////////
    //              LaBr3Ce Encasement
    ///////////////////////////////////////////////////////
    
    //G4Tubs* Solid_LaBr3Ce_Encasement = new G4Tubs("Solid_LaBr3CeEncasement", LaBr3Ce_encasement_innerRadius*mm, LaBr3Ce_encasement_outerRadius*mm, ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm, 0.*deg, 360*deg);
    G4Tubs* Solid_LaBr3Ce_Encasement_base = new G4Tubs("Solid_LaBr3Ce_Encasement_base", LaBr3Ce_encasement_innerRadius*mm, LaBr3Ce_encasement_outerRadius*mm, (LaBr3Ce_crystalBase_axialLength/2.0)*mm, 0.*deg, 360*deg);
    G4Cons* Solid_LaBr3Ce_Encasement_cone = new G4Cons("Solid_LaBr3Ce_Encasement_cone", LaBr3Ce_encasement_innerRadius*mm, LaBr3Ce_encasement_outerRadius*mm, LaBe3CeCrystal_cone_outerRadius*mm, (LaBr3Ce_encasement_cone_innerRadius)*mm, (LaBr3Ce_crystalCone_axialLength/2.0)*mm, 0.*deg, 360.0*deg);

    G4LogicalVolume* Logic_LaBr3Ce_Encasement_base = new G4LogicalVolume(Solid_LaBr3Ce_Encasement_base, G4_Al_Material, "LogicLaBr3CeLaBr3CeEncasement", 0, 0, 0);
    G4LogicalVolume* Logic_LaBr3Ce_Encasement_cone = new G4LogicalVolume(Solid_LaBr3Ce_Encasement_cone, G4_Al_Material, "LogicLaBr3CeLaBr3CeEncasement", 0, 0, 0);

    //////////////////////////////////////////////////////////
    //              LaBr3Ce Internal Vacuum
    //////////////////////////////////////////////////////////
    
    G4Tubs* Solid_LaBr3Ce_InternalVacuum = new G4Tubs("Solid_LaBr3CeInternalVacuum", 0.0*mm, (LaBr3Ce_encasement_outerRadius)*mm, ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm, 0.*deg, 360*deg);

    G4LogicalVolume* Logic_LaBr3Ce_InternalVacuum[numberOf_LaBr3Ce];
    
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        if(setPreconfiguredVersion)
        {
            if(LaBr3CeSetupVersion==-1)
            {
                Logic_LaBr3Ce_InternalVacuum[i] = new G4LogicalVolume(Solid_LaBr3Ce_InternalVacuum, G4_Galactic_Material, "LogicLaBr3CeInternalVacuum", 0, 0, 0);
            }
            if(LaBr3CeSetupVersion==0)
            {
                Logic_LaBr3Ce_InternalVacuum[i] = new G4LogicalVolume(Solid_LaBr3Ce_InternalVacuum, G4_Galactic_Material, "LogicLaBr3CeInternalVacuum", 0, 0, 0);
            }
            else if(LaBr3CeSetupVersion>0)
            {
                double a = LaBr3Ce_window_axialLength+LaBr3Ce_crystalCone_axialLength;
                double b = LaBr3Ce_encasement_outerRadius-LaBr3Ce_window_smallerRadius;
                double c = (LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength);
                double d = c*(b/a);
                double largerRadius = LaBr3Ce_window_smallerRadius + d; // mm
                
                G4Cons* Solid_LaBr3CeInternalVacuum_cone = new G4Cons("Solid_LaBr3CeInternalVacuum_cone", 0.0*mm, largerRadius*mm, 0.0*mm, LaBr3Ce_window_smallerRadius*mm, ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm, 0.*deg, 360.0*deg);
                
                G4RotationMatrix* rm1_LaBr3CeInternalVacuum = new G4RotationMatrix();
                G4VSolid* Solid_LaBr3CeInternalVacuum_mod = new G4IntersectionSolid("Solid_LaBr3CeInternalVacuum_mod", Solid_LaBr3CeInternalVacuum_cone, Solid_LaBr3Ce_InternalVacuum, rm1_LaBr3CeInternalVacuum, G4ThreeVector());
                
                Logic_LaBr3Ce_InternalVacuum[i] = new G4LogicalVolume(Solid_LaBr3CeInternalVacuum_mod, G4_Galactic_Material, "LogicLaBr3CeInternalVacuum", 0, 0, 0);
            }
        }
        else
        {
            //--------------------
            //      OPTION 1
            Logic_LaBr3Ce_InternalVacuum[i] = new G4LogicalVolume(Solid_LaBr3Ce_InternalVacuum, G4_Galactic_Material, "LogicLaBr3CeInternalVacuum", 0, 0, 0);
            
            //--------------------
            //      OPTION2
            
         /*   double a = LaBr3Ce_window_axialLength+LaBr3Ce_crystalCone_axialLength;
            double b = LaBr3Ce_encasement_outerRadius-LaBr3Ce_window_smallerRadius;
            double c = (LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength);
            double d = c*(b/a);
            double largerRadius = LaBr3Ce_window_smallerRadius + d; // mm
            
            G4Cons* Solid_LaBr3CeInternalVacuum_cone = new G4Cons("Solid_LaBr3CeInternalVacuum_cone", 0.0*mm, largerRadius*mm, 0.0*mm, LaBr3Ce_window_smallerRadius*mm, ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm, 0.*deg, 360.0*deg);
            
            G4RotationMatrix* rm1_LaBr3CeInternalVacuum = new G4RotationMatrix();
            G4VSolid* Solid_LaBr3CeInternalVacuum_mod = new G4IntersectionSolid("Solid_LaBr3CeInternalVacuum_mod", Solid_LaBr3CeInternalVacuum_cone, Solid_LaBr3Ce_InternalVacuum, rm1_LaBr3CeInternalVacuum, G4ThreeVector());
            
            Logic_LaBr3Ce_InternalVacuum[i] = new G4LogicalVolume(Solid_LaBr3CeInternalVacuum_mod, G4_Galactic_Material, "LogicLaBr3CeInternalVacuum", 0, 0, 0);
          */
        }
    }
    
    
    ////////////////////////////////////////////////////
    //               LaBr3Ce INITIALIZATION
    ////////////////////////////////////////////////////
    
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        G4RotationMatrix    LaBr3Ce_relativeInternalVacuum_rotm;
        LaBr3Ce_Encasement_transform[i] = G4Transform3D(LaBr3Ce_relativeInternalVacuum_rotm, G4ThreeVector(0.0, 0.0, -((LaBr3Ce_window_axialLength+LaBr3Ce_crystalCone_axialLength)/2.0)*mm));

        LaBr3Ce_LaBr3CeCrystal_base_transform[i] = G4Transform3D(LaBr3Ce_relativeInternalVacuum_rotm, G4ThreeVector(0.0, 0.0, -((LaBr3Ce_window_axialLength+LaBr3Ce_crystalCone_axialLength)/2.0)*mm));
        LaBr3Ce_LaBr3CeCrystal_cone_transform[i] = G4Transform3D(LaBr3Ce_relativeInternalVacuum_rotm, G4ThreeVector(0.0, 0.0, +((LaBr3Ce_crystalBase_axialLength-LaBr3Ce_window_axialLength)/2.0)*mm));

        LaBr3Ce_Window_position[i] = G4ThreeVector(0.0, 0.0, ((LaBr3Ce_crystal_axialLength)/2.0)*mm);
        LaBr3Ce_Window_transform[i] = G4Transform3D(LaBr3Ce_relativeInternalVacuum_rotm, LaBr3Ce_Window_position[i]);
        
        ////    STANDARD
        //LaBr3Ce_InternalVacuum_position[i] = (LaBr3Ce_Distance[i] + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        LaBr3Ce_InternalVacuum_position[i] = (LaBr3Ce_Distance[i] + 0.4*mm + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        ////    TEMP
        //LaBr3Ce_InternalVacuum_position[i] = (17.0*cm + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        //  For comparison for rates
        //LaBr3Ce_InternalVacuum_position[i] = (13.3*cm + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        //  For comparison for peak-to-total
        //LaBr3Ce_InternalVacuum_position[i] = (13.2*cm + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        //LaBr3Ce_InternalVacuum_position[i] = (8.9632*cm + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        
        if(setPreconfiguredVersion && LaBr3CeSetupVersion==-1)
        {
            //  LaBr3Ce 1
            LaBr3Ce_Presence[0] = true;
            LaBr3Ce_Distance[0] = 246.4*mm;
            LaBr3Ce_theta[0] = 90.0*deg;
            LaBr3Ce_phi[0] = 151.0*deg;
            
            //  LaBr3Ce 2
            LaBr3Ce_Presence[1] = true;
            LaBr3Ce_Distance[1] = 247.65*mm;
            LaBr3Ce_theta[1] = 135.0*deg;
            LaBr3Ce_phi[1] = 180.0*deg;
            
            //  LaBr3Ce 3
            LaBr3Ce_Presence[2] = true;
            LaBr3Ce_Distance[2] = 217.8*mm;
            LaBr3Ce_theta[2] = 135.0*deg;
            LaBr3Ce_phi[2] = 0.0*deg;
            
            //  LaBr3Ce 4
            LaBr3Ce_Presence[3] = true;
            LaBr3Ce_Distance[3] = 236.8*mm;
            LaBr3Ce_theta[3] = 90.0*deg;
            LaBr3Ce_phi[3] = 62.0*deg;
            
            double x = std::sin(LaBr3Ce_theta[i])*std::cos(LaBr3Ce_phi[i]);
            double y = std::sin(LaBr3Ce_theta[i])*std::sin(LaBr3Ce_phi[i]);
            double z = std::cos(LaBr3Ce_theta[i]);
            G4ThreeVector positionDetector(x, y, z);
            
            LaBr3Ce_InternalVacuum_position[i] = (LaBr3Ce_Distance[i] + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*positionDetector.unit();
            
            //LaBr3Ce_InternalVacuum_position[i] = (LaBr3Ce_Distance[i] + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*G4ThreeVector(std::sin(LaBr3Ce_theta[i]) * std::cos(LaBr3Ce_phi[i]), std::sin(LaBr3Ce_theta[i]) * std::sin(LaBr3Ce_phi[i]), std::cos(LaBr3Ce_theta[i]));
        }
        else if(setPreconfiguredVersion && LaBr3CeSetupVersion>=0 && i>=20)
        {
            G4double pentagonalDistance;
            
            if(LaBr3CeSetupVersion==0)
            {
                pentagonalDistance = 14.8*cm; // mm
            }
            else if(LaBr3CeSetupVersion==8)
            {
                pentagonalDistance = 12.1*cm; // mm
            }
            else if(LaBr3CeSetupVersion==9)
            {
                pentagonalDistance = 12.2*cm; // mm
            }

            int pentIndex = i-20;
            LaBr3Ce_InternalVacuum_position[i] = (pentagonalDistance + ((LaBr3Ce_crystal_axialLength+LaBr3Ce_window_axialLength)/2.0)*mm)*vertex_pentagonFaces_truncatedIcosahedron[pentIndex].unit();
        }
        
        G4ThreeVector positionVector = LaBr3Ce_InternalVacuum_position[i].unit();
        
        G4ThreeVector positionVector_z = positionVector.unit();
        G4ThreeVector positionVector_y = (positionVector.orthogonal()).unit();
        G4ThreeVector positionVector_x = (positionVector_y.cross(positionVector_z)).unit();
        positionVector_z = -positionVector_z;
        positionVector_y = -positionVector_y;
        
        G4RotationMatrix rotmPrime(positionVector_x, positionVector_y, positionVector_z);
        
        LaBr3Ce_rotm[i] = rotmPrime;
        
        LaBr3Ce_InternalVacuum_transform[i] = G4Transform3D(LaBr3Ce_rotm[i],LaBr3Ce_InternalVacuum_position[i]);

        /////////////////////////////
        //          LaBr3Ce
        if(LaBr3Ce_Presence[i] == true)
        {
            
            new G4PVPlacement(LaBr3Ce_Encasement_transform[i],   // transformation matrix
                              Logic_LaBr3Ce_Encasement_base,       // its logical volume
                              "LaBr3CeEncasement",       // its name
                              Logic_LaBr3Ce_InternalVacuum[i],         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            new G4PVPlacement(LaBr3Ce_LaBr3CeCrystal_cone_transform[i],   // transformation matrix
                              Logic_LaBr3Ce_Encasement_cone,       // its logical volume
                              "LaBr3CeEncasement",       // its name
                              Logic_LaBr3Ce_InternalVacuum[i],         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps

            new G4PVPlacement(LaBr3Ce_Window_transform[i],   // transformation matrix
                              Logic_LaBr3Ce_Window,       // its logical volume
                              "LaBr3CeWindow",       // its name
                              Logic_LaBr3Ce_InternalVacuum[i],         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
            if(i<20)
            {
                new G4PVPlacement(LaBr3Ce_LaBr3CeCrystal_base_transform[i],
                                  Logic_LaBr3Ce_LaBr3CeCrystal_base,       // its logical volume
                                  "LaBr3CeHPGeCrystal",       // its name
                                  Logic_LaBr3Ce_InternalVacuum[i],    // its mother  volume
                                  false,           // no boolean operations
                                  i,               // copy number
                                  fCheckOverlaps); // checking overlaps
                
                new G4PVPlacement(LaBr3Ce_LaBr3CeCrystal_cone_transform[i],
                                  Logic_LaBr3Ce_LaBr3CeCrystal_cone,       // its logical volume
                                  "LaBr3CeHPGeCrystal",       // its name
                                  Logic_LaBr3Ce_InternalVacuum[i],    // its mother  volume
                                  false,           // no boolean operations
                                  i,               // copy number
                                  fCheckOverlaps); // checking overlaps
            }
            else
            {
                new G4PVPlacement(LaBr3Ce_LaBr3CeCrystal_base_transform[i],
                                  Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalBase,       // its logical volume
                                  "LaBr3CeHPGeCrystal",       // its name
                                  Logic_LaBr3Ce_InternalVacuum[i],    // its mother  volume
                                  false,           // no boolean operations
                                  i,               // copy number
                                  fCheckOverlaps); // checking overlaps
                
                new G4PVPlacement(LaBr3Ce_LaBr3CeCrystal_cone_transform[i],
                                  Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalCone,       // its logical volume
                                  "LaBr3CeHPGeCrystal",       // its name
                                  Logic_LaBr3Ce_InternalVacuum[i],    // its mother  volume
                                  false,           // no boolean operations
                                  i,               // copy number
                                  fCheckOverlaps); // checking overlaps
            }
            
            new G4PVPlacement(LaBr3Ce_InternalVacuum_transform[i],
                              Logic_LaBr3Ce_InternalVacuum[i],
                              "LaBr3CeInternalVacuum",       // its name
                              LogicVacuumChamber,         // its mother  volume
                              false,           // no boolean operations
                              i,               // copy number
                              fCheckOverlaps); // checking overlaps
            
        }
    }

    
    //////////////////////////////////////////////////////////
    //              LaBr3Ce - Pentagonal
    //////////////////////////////////////////////////////////
    
    //  Minimum angle between the normal vectors of hexagonal and pentagonal faces: (37.3774/2.0)
    
    int setupN_pentagonal = 2;
    
    G4Tubs *Solid_LaBr3Ce_pentagonal;
    if(setupN_pentagonal==0)
    {
        //      Setup 3
        Solid_LaBr3Ce_pentagonal = new G4Tubs("Solid_LaBr3Ce_pentagonal", 0.0*mm, (30.0)*mm, (200.0/2.0)*mm, 0.*deg, 360*deg);
    }
    else if(setupN_pentagonal==1)
    {
        //      Setup 8
        Solid_LaBr3Ce_pentagonal = new G4Tubs("Solid_LaBr3Ce_pentagonal", 0.0*mm, (27.0)*mm, (200.0/2.0)*mm, 0.*deg, 360*deg);
    }
    else if(setupN_pentagonal==2)
    {
        //      Setup *
        Solid_LaBr3Ce_pentagonal = new G4Tubs("Solid_LaBr3Ce_pentagonal", 0.0*mm, (39.0)*mm, (200.0/2.0)*mm, 0.*deg, 360*deg);
    }

    G4LogicalVolume *Logic_LaBr3Ce_pentagonal = new G4LogicalVolume(Solid_LaBr3Ce_pentagonal, G4_Galactic_Material, "Logic_LaBr3Ce_pentagonal", 0, 0, 0);
    
    for(int i=0; i<(int) vertex_pentagonFaces_truncatedIcosahedron.size(); i++)
    {
        double pentagonalDistance = ((sqrt(10.0*(125.0+41.0*sqrt(5.0)))/20.0)/((3.0*sqrt(3.0)+sqrt(15.0))/4.0))*LaBr3Ce_Distance[0]; // mm
        
        if(setupN_pentagonal==2)
        {
            pentagonalDistance = 13.2*10.0; // mm
        }
        
        pentagonalDistance = 13.1*10.0; // mm

        G4ThreeVector positionVector_LaBr3Ce_pentagonal = (pentagonalDistance*mm + (200.0/2.0)*mm)*vertex_pentagonFaces_truncatedIcosahedron[i].unit();

        G4ThreeVector positionVector_LaBr3Ce_pentagonal_z = positionVector_LaBr3Ce_pentagonal.unit();
        G4ThreeVector positionVector_LaBr3Ce_pentagonal_y = (positionVector_LaBr3Ce_pentagonal.orthogonal()).unit();
        G4ThreeVector positionVector_LaBr3Ce_pentagonal_x = (positionVector_LaBr3Ce_pentagonal_y.cross(positionVector_LaBr3Ce_pentagonal_z)).unit();
        positionVector_LaBr3Ce_pentagonal_z = -positionVector_LaBr3Ce_pentagonal_z;
        positionVector_LaBr3Ce_pentagonal_y = -positionVector_LaBr3Ce_pentagonal_y;
        
        G4RotationMatrix rotm_LaBr3Ce_pentagonal(positionVector_LaBr3Ce_pentagonal_x, positionVector_LaBr3Ce_pentagonal_y, positionVector_LaBr3Ce_pentagonal_z);
        
        G4Transform3D LaBr3Ce_pentagonal_transform = G4Transform3D(rotm_LaBr3Ce_pentagonal, positionVector_LaBr3Ce_pentagonal);
        
        /*
        new G4PVPlacement(LaBr3Ce_pentagonal_transform,
                          Logic_LaBr3Ce_InternalVacuum[0],
                          "LaBr3Ce_pentagonal",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        */
        
        /*
        new G4PVPlacement(LaBr3Ce_pentagonal_transform,
                          Logic_LaBr3Ce_pentagonal,
                          "LaBr3Ce_pentagonal",       // its name
                          LogicVacuumChamber,         // its mother  volume
                          false,           // no boolean operations
                          0,               // copy number
                          fCheckOverlaps); // checking overlaps
        */
        
    }
    
    
    
    //////////////////////////////////////////////////
    //      K600 SPECTROMETER INITIIALIZATION       //
    //////////////////////////////////////////////////
    
    
    // A field object is held by a field manager
    // Find the global Field Manager
    G4TransportationManager* tmanagerMagneticField = G4TransportationManager::GetTransportationManager();
    tmanagerMagneticField->GetPropagatorInField()->SetLargestAcceptableStep(1*mm);
    G4double minStepMagneticField = 0.0025*mm ;
    
    
    //////////////////////////////////////////////////////
    //              K600 - QUADRUPOLE
    //////////////////////////////////////////////////////
    
    if(K600_Quadrupole)
    {
        
        G4RotationMatrix* K600_Q_MagField_rotm = new G4RotationMatrix;
        //K600_Q_MagField_rotm->rotateX(90.*deg);
        
        K600_Quadrupole_transform = G4Transform3D(K600_Quadrupole_rotm, K600_Quadrupole_CentrePosition);
        
        G4Box* Solid_K600_Quadrupole = new G4Box("Solid_K600_Quadrupole", (50./2)*cm, (50./2)*cm, (30./2)*cm);
        
        G4LogicalVolume* Logic_K600_Quadrupole = new G4LogicalVolume(Solid_K600_Quadrupole, G4_Galactic_Material,"Logic_K600_Quadrupole",0,0,0);
        
        ////    IDEAL MAGNETIC FIELD for QUADRUPOLE
        if(Ideal_Quadrupole)
        {
            MagneticField_K600_Q = new G4QuadrupoleMagField(K600_Q_gradient, K600_Quadrupole_CentrePosition, K600_Q_MagField_rotm);
            fEquationMagneticField_K600_Q = new G4Mag_UsualEqRhs(MagneticField_K600_Q);
            
            fieldManagerMagneticField_K600_Q = new G4FieldManager(MagneticField_K600_Q);
            
            stepperMagneticField_K600_Q = new G4ClassicalRK4( fEquationMagneticField_K600_Q );
            fieldManagerMagneticField_K600_Q -> SetDetectorField(MagneticField_K600_Q);
            
            fChordFinder_K600_Q = new G4ChordFinder( MagneticField_K600_Q, minStepMagneticField, stepperMagneticField_K600_Q);
            
            Logic_K600_Quadrupole -> SetFieldManager(fieldManagerMagneticField_K600_Q, true) ;
            
            
            //G4BlineTracer* theBlineTool = new G4BlineTracer();
            //theBlineTool->ComputeBlines();
        }
        
        ////    MAPPED MAGNETIC FIELD for QUADRUPOLE
        if(Mapped_Quadrupole)
        {
            
            G4double z_Q_Offset = 4.4*mm+ 100*cm;
            
            G4MagneticField* PurgMagField = new MagneticFieldMapping("../K600-ALBA/MagneticFieldMaps/Quadrupole_MagneticFieldMap.TABLE", z_Q_Offset);
            fEquationMagneticField_K600_Q = new G4Mag_UsualEqRhs(PurgMagField);
            
            fieldManagerMagneticField_K600_Q = new G4FieldManager(PurgMagField);
            
            stepperMagneticField_K600_Q = new G4ClassicalRK4( fEquationMagneticField_K600_Q );
            fieldManagerMagneticField_K600_Q -> SetDetectorField(PurgMagField);
            
            fChordFinder_K600_Q = new G4ChordFinder( PurgMagField, minStepMagneticField, stepperMagneticField_K600_Q);
            
            Logic_K600_Quadrupole -> SetFieldManager(fieldManagerMagneticField_K600_Q, true) ;
            
            
            //G4BlineTracer* theBlineTool = new G4BlineTracer();
            
        }
        
        PhysiK600_Quadrupole = new G4PVPlacement(K600_Quadrupole_transform,
                                                 Logic_K600_Quadrupole,       // its logical volume
                                                 "K600_Quadrupole",       // its name
                                                 LogicWorld,         // its mother  volume
                                                 false,           // no boolean operations
                                                 0,               // copy number
                                                 fCheckOverlaps); // checking overlaps
        
    }
    
    //////////////////////////////////////////////////////
    //              K600 - DIPOLE 1
    //////////////////////////////////////////////////////
    
    /*
     // magnetic field ----------------------------------------------------------
     MagneticField_K600_D1 = new G4UniformMagField(G4ThreeVector(0., K600_Dipole1_BZ, 0.));
     fieldManagerMagneticField_K600_D1 = new G4FieldManager();
     fieldManagerMagneticField_K600_D1->SetDetectorField(MagneticField_K600_D1);
     fieldManagerMagneticField_K600_D1->CreateChordFinder(MagneticField_K600_D1);
     G4bool forceToAllDaughters = true;
     */
    
    if(K600_Dipole1)
    {
        ////    MAGNETIC FIELD for DIPOLE 1
        MagneticField_K600_D1 = new G4UniformMagField(G4ThreeVector(0., K600_Dipole1_BZ, 0.));
        fEquationMagneticField_K600_D1 = new G4Mag_UsualEqRhs(MagneticField_K600_D1);
        
        fieldManagerMagneticField_K600_D1 = new G4FieldManager(MagneticField_K600_D1);
        
        stepperMagneticField_K600_D1 = new G4ClassicalRK4( fEquationMagneticField_K600_D1 );
        fieldManagerMagneticField_K600_D1 -> SetDetectorField(MagneticField_K600_D1);
        
        fChordFinder_K600_D1 = new G4ChordFinder( MagneticField_K600_D1, minStepMagneticField, stepperMagneticField_K600_D1);
        
        
        /////////////////////////////////////////////
        
        K600_Dipole1_transform = G4Transform3D(K600_Dipole1_rotm, K600_Dipole1_CentrePosition);
        
        //G4Box* Solid_K600_Dipole1 = new G4Box("Solid_K600_Dipole1", (50./2)*cm, (50./2)*cm, (30./2)*cm);
        //G4Tubs* Solid_K600_Dipole1 = new G4Tubs("Solid_K600_Dipole1", 50.*cm, 100.0*cm, 30.*cm, 0.*deg, 40.*deg);
        G4Tubs* Solid_K600_Dipole1 = new G4Tubs("Solid_K600_Dipole1", 30.*cm, 150.0*cm, 30.*cm, 0.*deg, 40.*deg);
        
        G4LogicalVolume* Logic_K600_Dipole1 = new G4LogicalVolume(Solid_K600_Dipole1, G4_Galactic_Material,"Logic_K600_Dipole1",0,0,0);
        Logic_K600_Dipole1 -> SetFieldManager(fieldManagerMagneticField_K600_D1, true) ;
        
        // Register the field and its manager for deleting
        //G4AutoDelete::Register(MagneticField_K600_D1);
        //G4AutoDelete::Register(fieldManagerMagneticField_K600_D1);
        
        PhysiK600_Dipole1 = new G4PVPlacement(K600_Dipole1_transform,
                                              Logic_K600_Dipole1,       // its logical volume
                                              "K600_Dipole1",       // its name
                                              LogicWorld,         // its mother  volume
                                              false,           // no boolean operations
                                              0,               // copy number
                                              fCheckOverlaps); // checking overlaps
        
    }
    
    
    
    //////////////////////////////////////////////////////
    //              K600 - DIPOLE 2
    //////////////////////////////////////////////////////
    
    if(K600_Dipole2)
    {
        
        ////    MAGNETIC FIELD for DIPOLE 2
        MagneticField_K600_D2 = new G4UniformMagField(G4ThreeVector(0., K600_Dipole2_BZ, 0.));
        fEquationMagneticField_K600_D2 = new G4Mag_UsualEqRhs(MagneticField_K600_D2);
        
        fieldManagerMagneticField_K600_D2 = new G4FieldManager(MagneticField_K600_D2);
        
        stepperMagneticField_K600_D2 = new G4ClassicalRK4( fEquationMagneticField_K600_D2 );
        fieldManagerMagneticField_K600_D2 -> SetDetectorField(MagneticField_K600_D2);
        
        fChordFinder_K600_D2 = new G4ChordFinder( MagneticField_K600_D2, minStepMagneticField, stepperMagneticField_K600_D2);
        
        
        /////////////////////////////////////////////
        
        K600_Dipole2_transform = G4Transform3D(K600_Dipole2_rotm, K600_Dipole2_CentrePosition);
        
        //G4Tubs* Solid_K600_Dipole2 = new G4Tubs("Solid_K600_Dipole2", 50.*cm, 100.0*cm, 30.*cm, 50.*deg, 70.*deg);
        G4Tubs* Solid_K600_Dipole2 = new G4Tubs("Solid_K600_Dipole2", 30.*cm, 150.0*cm, 30.*cm, 50.*deg, 70.*deg);
        
        
        G4LogicalVolume* Logic_K600_Dipole2 = new G4LogicalVolume(Solid_K600_Dipole2, G4_Galactic_Material,"Logic_K600_Dipole2",0,0,0);
        Logic_K600_Dipole2 -> SetFieldManager(fieldManagerMagneticField_K600_D2, true) ;
        
        PhysiK600_Dipole2 = new G4PVPlacement(K600_Dipole2_transform,
                                              Logic_K600_Dipole2,       // its logical volume
                                              "K600_Dipole2",       // its name
                                              LogicWorld,         // its mother  volume
                                              false,           // no boolean operations
                                              0,               // copy number
                                              fCheckOverlaps); // checking overlaps
        
        
    }
    
    
    
    
    //////////////////////////////////////////////////////////
    //                      VISUALISATION
    //////////////////////////////////////////////////////////
    
    G4VisAttributes* World_VisAtt= new G4VisAttributes(G4Colour(0., 0., 0.));
    World_VisAtt->SetVisibility(false);
    LogicWorld->SetVisAttributes(World_VisAtt);
    
    //G4VisAttributes* Target_VisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //LogicTarget->SetVisAttributes(Target_VisAtt);
    
    
    ////////////////////////////
    //  Vacuum Chamber
    
    G4VisAttributes* VacuumChamber_VisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VacuumChamber_VisAtt->SetVisibility(false);
    LogicVacuumChamber->SetVisAttributes(VacuumChamber_VisAtt);
    
    
    ////////////////////////////////
    //      VDC VISUALIZATION
    ////////////////////////////////
    
    //  VDC - Asm
    G4VisAttributes* VDC_Asm_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    VDC_Asm_VisAtt->SetVisibility(false);
    
    //  VDC Asm USDS
    G4VisAttributes* VDC_Asm_USDS_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    VDC_Asm_USDS_VisAtt->SetVisibility(false);
    
    //  VDC - STESALIT HV Frames, Upstream, Middle and Downstream
    G4VisAttributes* VDC_HV_Frame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_HV_Frame_VisAtt->SetForceSolid(true);
    //VDC_HV_Frame_VisAtt->SetVisibility(false);
    
    //  VDC - HIGH VOLTAGE PLANE
    G4VisAttributes* VDC_HV_Plane_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    //VDC_HV_Plane_VisAtt->SetForceSolid(true);
    //VDC_HV_Plane_VisAtt->SetVisibility(false);
    
    //  VDC - GAS FRAME
    G4VisAttributes* VDC_GasFrame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_GasFrame_VisAtt->SetForceSolid(true);
    //VDC_GasFrame_VisAtt->SetVisibility(false);
    
    //  VDC - STESALIT X-WIRE/U-WIRE FRAME
    G4VisAttributes* VDC_XU_Frame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_XU_Frame_VisAtt->SetForceSolid(true);
    //VDC_XU_Frame_VisAtt->SetVisibility(false);
    
    //  VDC - PCB X-WIRE/U-WIRE FRAME
    G4VisAttributes* VDC_XU_PCBFrame_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.5, 0.35));
    VDC_XU_PCBFrame_VisAtt->SetForceSolid(true);
    //VDC_XU_PCBFrame_VisAtt->SetVisibility(false);
    
    //  VDC - ALUMINIUM OUTER FRAME
    G4VisAttributes* VDC_Al_Frame_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VDC_Al_Frame_VisAtt->SetForceSolid(true);
    
    //  VDC - MYLAR PLANE
    G4VisAttributes* VDC_MYLAR_Plane_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //VDC_MYLAR_Plane_VisAtt->SetForceSolid(true);
    
    //  VDC - X WIRES
    G4VisAttributes* VDC_X_WIRE_VisAtt = new G4VisAttributes(G4Colour(0., 0.7, 0.7));
    VDC_X_WIRE_VisAtt->SetForceSolid(true);
    
    //  VDC - X GUARD WIRES
    G4VisAttributes* VDC_X_GUARDWIRE_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VDC_X_GUARDWIRE_VisAtt->SetVisibility(false);
    VDC_X_GUARDWIRE_VisAtt->SetForceSolid(true);
    
    //  VDC - U WIRES
    G4VisAttributes* VDC_U_WIRE_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.));
    VDC_U_WIRE_VisAtt->SetForceSolid(true);
    
    //  VDC - U GUARD WIRES
    G4VisAttributes* VDC_U_GUARDWIRE_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    VDC_U_GUARDWIRE_VisAtt->SetVisibility(false);
    VDC_U_GUARDWIRE_VisAtt->SetForceSolid(true);
    
    Logic_VDC_GasFrame->SetVisAttributes(VDC_GasFrame_VisAtt);
    Logic_VDC_XU_Frame->SetVisAttributes(VDC_XU_Frame_VisAtt);
    Logic_VDC_XU_PCBFrame->SetVisAttributes(VDC_XU_PCBFrame_VisAtt);
    Logic_VDC_Al_Frame->SetVisAttributes(VDC_Al_Frame_VisAtt);
    Logic_VDC_MYLAR_Plane->SetVisAttributes(VDC_MYLAR_Plane_VisAtt);
    Logic_VDC_X_WIRE->SetVisAttributes(VDC_X_WIRE_VisAtt);
    Logic_VDC_X_GUARDWIRE->SetVisAttributes(VDC_X_GUARDWIRE_VisAtt);
    Logic_VDC_X_GUARDWIRE_Thick->SetVisAttributes(VDC_X_GUARDWIRE_VisAtt);
    
    for(G4int k=0; k<143; k++)
    {
        Logic_VDC_U_WIRE[k]->SetVisAttributes(VDC_U_WIRE_VisAtt);
    }
    
    for(G4int k=0; k<144; k++)
    {
        Logic_VDC_U_GUARDWIRE[k]->SetVisAttributes(VDC_U_GUARDWIRE_VisAtt);
    }
    
    for(G4int k=0; k<2; k++)
    {
        Logic_VDC_U_GUARDWIRE_Thick[k]->SetVisAttributes(VDC_U_GUARDWIRE_VisAtt);
    }
    
    for(G4int k=0; k<3; k++)
    {
        Logic_VDC_HV_Frame[k]->SetVisAttributes(VDC_HV_Frame_VisAtt);
        Logic_VDC_HV_Plane[k]->SetVisAttributes(VDC_HV_Plane_VisAtt);
    }
    
    for(G4int i=0; i<2; i++)
    {
        for(G4int j=0; j<2; j++)
        {
            Logic_VDC_SenseRegion_USDS[i][j]->SetVisAttributes(VDC_Asm_USDS_VisAtt);
        }
    }
    
    for(G4int i=0; i<numberOf_VDC; i++)
    {
        Logic_VDC_Asm[i]->SetVisAttributes(VDC_Asm_VisAtt);
    }
    
    
    
    
    //////////////////////////////////
    //      PADDLE VISUALIZATION
    //////////////////////////////////
    
    G4VisAttributes* scintillatorVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    for(G4int i=0; i<3; i++)
    {Logic_PADDLE[i] -> SetVisAttributes(scintillatorVisAtt);}
    
    
    
    //////////////////////////////////
    //      HAGAR VISUALIZATION
    //////////////////////////////////
    
    //  HAGAR - NaI Crystal
    G4VisAttributes* HAGAR_NaICrystal_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    HAGAR_NaICrystal_VisAtt->SetForceSolid(true);
    Logic_HAGAR_NaICrystal->SetVisAttributes(HAGAR_NaICrystal_VisAtt);
    
    //  HAGAR - Annulus
    G4VisAttributes* HAGAR_Annulus_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    Logic_HAGAR_Annulus->SetVisAttributes(HAGAR_Annulus_VisAtt);
    
    //  HAGAR - Front Disc
    G4VisAttributes* HAGAR_FrontDisc_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    HAGAR_FrontDisc_VisAtt->SetForceSolid(true);
    Logic_HAGAR_FrontDisc->SetVisAttributes(HAGAR_FrontDisc_VisAtt);
    
    
    /////////////////////////////////
    //      CAKE VISUALIZATION
    /////////////////////////////////
    
    //  Asm
    G4VisAttributes* CAKE_Asm_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    CAKE_Asm_VisAtt->SetVisibility(false);
    //CAKE_Asm_VisAtt->SetForceSolid(false);
    //CAKE_Asm_VisAtt->SetForceWireframe(false);
    
    //  PCB
    G4VisAttributes* CAKE_PCB_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    CAKE_PCB_VisAtt->SetForceSolid(true);
    
    //  AA - RS
    //G4VisAttributes* CAKE_AA_RS_VisAtt = new G4VisAttributes(G4Colour(0., 0.7, 0.7));
    G4VisAttributes* CAKE_AA_RS_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    //CAKE_AA_RS_VisAtt->SetForceSolid(true);
    CAKE_AA_RS_VisAtt->SetForceLineSegmentsPerCircle(100);
    
    //  Silicon Wafer
    G4VisAttributes* CAKE_SiliconWafer_VisAtt = new G4VisAttributes(G4Colour(0., 0.2, 0.2));
    //CAKE_SiliconWafer_VisAtt->SetForceSolid(true);
    CAKE_SiliconWafer_VisAtt->SetForceLineSegmentsPerCircle(100);
    
    
    //  DL
    G4VisAttributes* CAKE_DL_VisAtt = new G4VisAttributes(G4Colour(1., 0., 0.));
    //CAKE_DL_VisAtt->SetForceSolid(true);
    
    //  2M
    G4VisAttributes* CAKE_2M_VisAtt = new G4VisAttributes(G4Colour(1., 1., 1.));
    CAKE_2M_VisAtt->SetVisibility(false);
    //CAKE_2M_VisAtt->SetForceSolid(false);
    //CAKE_2M_VisAtt->SetForceWireframe(false);
    
    
    Logic_CAKE_PCB->SetVisAttributes(CAKE_PCB_VisAtt);
    //Logic_CAKE_DL_Front->SetVisAttributes(CAKE_DL_VisAtt);
    //Logic_CAKE_DL_Back->SetVisAttributes(CAKE_DL_VisAtt);
    //Logic_CAKE_2M_Front->SetVisAttributes(CAKE_2M_VisAtt);
    //Logic_CAKE_2M_Back->SetVisAttributes(CAKE_2M_VisAtt);
    
    for(G4int i=0; i<numberOf_CAKE; i++)
    {
        Logic_CAKE_Asm[i]->SetVisAttributes(CAKE_Asm_VisAtt);
        Logic_CAKE_SiliconWafer[i]->SetVisAttributes(CAKE_SiliconWafer_VisAtt);
        
    }
    
    for(G4int j=0; j<16; j++)
    {
        for(G4int l=0; l<8; l++)
        {
            Logic_CAKE_AA_RS[j][l]->SetVisAttributes(CAKE_AA_RS_VisAtt);
            Logic_CAKE_RS_punch[j][l]->SetVisAttributes(CAKE_AA_RS_VisAtt);
            
        }
    }
    
    
    //////////////////////////////////
    //      LaBr3Ce VISUALIZATION
    //////////////////////////////////
    
    //  LaBr3Ce - LaBr3Ce Crystal
    G4VisAttributes* LaBr3Ce_LaBr3CeCrystal_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
    LaBr3Ce_LaBr3CeCrystal_VisAtt->SetForceSolid(true);
    Logic_LaBr3Ce_LaBr3CeCrystal_base->SetVisAttributes(LaBr3Ce_LaBr3CeCrystal_VisAtt);
    Logic_LaBr3Ce_LaBr3CeCrystal_cone->SetVisAttributes(LaBr3Ce_LaBr3CeCrystal_VisAtt);
    
    G4VisAttributes* LaBr3Ce_LaBr3CeCrystal_pentagonal_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0));
    LaBr3Ce_LaBr3CeCrystal_pentagonal_VisAtt->SetForceSolid(true);
    Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalBase->SetVisAttributes(LaBr3Ce_LaBr3CeCrystal_pentagonal_VisAtt);
    Logic_LaBr3Ce_LaBr3CeCrystal_pentagonalCone->SetVisAttributes(LaBr3Ce_LaBr3CeCrystal_pentagonal_VisAtt);
    
    //  LaBr3Ce - Aluminium Window
    G4VisAttributes* LaBr3Ce_Window_VisAtt = new G4VisAttributes(G4Colour(0.0, 0.7, 0.0));
    LaBr3Ce_Window_VisAtt->SetForceSolid(true);
    Logic_LaBr3Ce_Window->SetVisAttributes(LaBr3Ce_Window_VisAtt);
    
    //  LaBr3Ce - Encasement
    G4VisAttributes* LaBr3Ce_Encasement_VisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    //LaBr3Ce_Encasement_VisAtt->SetForceSolid(true);
    Logic_LaBr3Ce_Encasement_base->SetVisAttributes(LaBr3Ce_Encasement_VisAtt);
    Logic_LaBr3Ce_Encasement_cone->SetVisAttributes(LaBr3Ce_Encasement_VisAtt);
    
    //  LaBr3Ce - InternalVacuumChamber
    G4VisAttributes* LaBr3Ce_InternalVacuumChamber_VisAtt = new G4VisAttributes(G4Colour(1., 1., 1.));
    LaBr3Ce_InternalVacuumChamber_VisAtt->SetVisibility(false);
    
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        Logic_LaBr3Ce_InternalVacuum[i]->SetVisAttributes(LaBr3Ce_InternalVacuumChamber_VisAtt);
    }

    
    //--------------------------------
    //  TEMP
    //  LaBr3Ce - LaBr3Ce Pentagonal
    G4VisAttributes* LaBr3Ce_pentagonal_VisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    LaBr3Ce_pentagonal_VisAtt->SetForceSolid(true);
    Logic_LaBr3Ce_pentagonal->SetVisAttributes(LaBr3Ce_pentagonal_VisAtt);
    
    
    G4cout << "//----------------------------------------" << G4endl;
    G4cout << "     LaBr3CeSetupVersion: " << LaBr3CeSetupVersion << G4endl;
    G4cout << "     LaBr3Ce detector volume [cm3]: " << vol_LaBr3Ce << G4endl;
    G4cout << "     Distance from target (origin) [cm]: " << LaBR3Ce_GlobalDistance/cm << G4endl;
    G4cout << "     Tapering pitch (LaBr3Ce_crystalCone_axialLength) [mm]: " << LaBr3Ce_crystalCone_axialLength << G4endl;
    G4cout << "     Tapering radius change [mm]: " << crystalBase_outerRadius << " -> " << LaBe3CeCrystal_cone_outerRadius << G4endl;
    G4cout << "     Tapering angle [deg]: " << (41.8103/2.0) << G4endl;
    G4cout << "//----------------------------------------" << G4endl;

    //
    //always return the physical World
    //
    return PhysiWorld;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<std::tuple<int, double, double>> DetectorConstruction::GetAngles_CLOVER()
{
    std::vector<std::tuple<int, double, double>> angles;
    
    for(G4int i=0; i<numberOf_CLOVER; i++)
    {
        angles.push_back(std::make_tuple(i, CLOVER_theta[i]/deg, CLOVER_phi[i]/deg));
    }
    
    return angles;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<std::tuple<int, double, double>> DetectorConstruction::GetAngles_ALBA_LaBr3Ce()
{
    std::vector<std::tuple<int, double, double>> angles;
    
    for(G4int i=0; i<numberOf_LaBr3Ce; i++)
    {
        angles.push_back(std::make_tuple(i, LaBr3Ce_theta[i]*57.2958, LaBr3Ce_phi[i]*57.2958));
    }
    
    return angles;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetupTruncatedIcosahedron()
{
    //------------------------------------
    //  Truncated Icosahedron Arrangement
    
    G4double C0 = (1.0 + sqrt(5.0))/4.0; // 0.809016994374947424102293417183
    G4double C1 = (1.0 + sqrt(5.0))/2.0; // 1.61803398874989484820458683437
    G4double C2 = (5.0 + sqrt(5.0))/4.0; // 1.80901699437494742410229341718
    G4double C3 = (2.0 + sqrt(5.0))/2.0; // 2.11803398874989484820458683437
    G4double C4 = 3.0 * (1.0 + sqrt(5.0))/4.0; // 2.427050983124842272306880251548
    
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.5,  0.0,   C4));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.5,  0.0,  -C4));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-0.5,  0.0,   C4));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-0.5,  0.0,  -C4));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C4,  0.5,  0.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C4, -0.5,  0.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C4,  0.5,  0.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C4, -0.5,  0.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.0,   C4,  0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.0,   C4, -0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.0,  -C4,  0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.0,  -C4, -0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(1.0,   C0,   C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(1.0,   C0,  -C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(1.0,  -C0,   C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(1.0,  -C0,  -C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-1.0,   C0,   C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-1.0,   C0,  -C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-1.0,  -C0,   C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-1.0,  -C0,  -C3));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C3,  1.0,   C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C3,  1.0,  -C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C3, -1.0,   C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C3, -1.0,  -C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C3,  1.0,   C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C3,  1.0,  -C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C3, -1.0,   C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C3, -1.0,  -C0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C0,   C3,  1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C0,   C3, -1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C0,  -C3,  1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C0,  -C3, -1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C0,   C3,  1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C0,   C3, -1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C0,  -C3,  1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C0,  -C3, -1.0));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.5,   C1,   C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.5,   C1,  -C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.5,  -C1,   C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(0.5,  -C1,  -C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-0.5,   C1,   C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-0.5,   C1,  -C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-0.5,  -C1,   C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-0.5,  -C1,  -C2));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C2,  0.5,   C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C2,  0.5,  -C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C2, -0.5,   C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C2, -0.5,  -C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C2,  0.5,   C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C2,  0.5,  -C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C2, -0.5,   C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C2, -0.5,  -C1));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C1,   C2,  0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C1,   C2, -0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C1,  -C2,  0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(C1,  -C2, -0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C1,   C2,  0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C1,   C2, -0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C1,  -C2,  0.5));
    vertices_truncatedIcosahedron.push_back(G4ThreeVector(-C1,  -C2, -0.5));
    
    //------------------------------------------------
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{0,  2, 18, 42, 38, 14});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{1,  3, 17, 41, 37, 13});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{2,  0, 12, 36, 40, 16});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{3,  1, 15, 39, 43, 19});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{4,  5, 23, 47, 45, 21});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{5,  4, 20, 44, 46, 22});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{6,  7, 26, 50, 48, 24});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{7,  6, 25, 49, 51, 27});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{8,  9, 33, 57, 56, 32});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{9,  8, 28, 52, 53, 29});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{10, 11, 31, 55, 54, 30});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{11, 10, 34, 58, 59, 35});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{12, 44, 20, 52, 28, 36});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{13, 37, 29, 53, 21, 45});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{14, 38, 30, 54, 22, 46});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{15, 47, 23, 55, 31, 39});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{16, 40, 32, 56, 24, 48});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{17, 49, 25, 57, 33, 41});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{18, 50, 26, 58, 34, 42});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{19, 43, 35, 59, 27, 51});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{0, 14, 46, 44, 12});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{1, 13, 45, 47, 15});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{2, 16, 48, 50, 18});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{3, 19, 51, 49, 17});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{4, 21, 53, 52, 20});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{5, 22, 54, 55, 23});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{6, 24, 56, 57, 25});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{7, 27, 59, 58, 26});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{8, 32, 40, 36, 28});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{9, 29, 37, 41, 33});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{10, 30, 38, 42, 34});
    verticesPerFace_truncatedIcosahedron.push_back(std::vector<double>{11, 35, 43, 39, 31});
    
    for(int i=0; i<(int) verticesPerFace_truncatedIcosahedron.size(); i++)
    {
        int nVerticesPerFace = verticesPerFace_truncatedIcosahedron[i].size();
        G4ThreeVector centralPoint(0.0, 0.0, 0.0);
        
        for(int j=0; j<(int) verticesPerFace_truncatedIcosahedron[i].size(); j++)
        {
            int vertexNumber = verticesPerFace_truncatedIcosahedron[i][j];
            centralPoint += vertices_truncatedIcosahedron[vertexNumber];
        }
        
        centralPoint = centralPoint.unit();
        vertex_allFaces_truncatedIcosahedron.push_back(std::make_tuple(nVerticesPerFace, centralPoint));
        
        if(nVerticesPerFace==5)
        {
            vertex_pentagonFaces_truncatedIcosahedron.push_back(centralPoint);
        }
        else if(nVerticesPerFace==6)
        {
            vertex_hexagonFaces_truncatedIcosahedron.push_back(centralPoint);
        }
    }
    
    //----------------------------------------------------------------------------
    double smallestAngle_hexgonal = 0.0; // deg
    for(int i=0; i<(int) vertex_hexagonFaces_truncatedIcosahedron.size(); i++)
    {
        for(int j=i+1; j<(int) vertex_hexagonFaces_truncatedIcosahedron.size(); j++)
        {
            double angle = vertex_hexagonFaces_truncatedIcosahedron[i].angle(vertex_hexagonFaces_truncatedIcosahedron[j])/deg;
            
            if(i==0 && j==1)
            {
                smallestAngle_hexgonal = angle;
            }
            else if(angle<smallestAngle_hexgonal)
            {
                smallestAngle_hexgonal = angle;
            }
        }
    }
    
    G4cout.precision(15);
    
    G4cout << "smallestAngle_hexgonal [deg]: " << smallestAngle_hexgonal << G4endl;
    
    //----------------------------------------------------------------------------
    double smallestAngle_pentagonal = 0.0;
    for(int i=0; i<(int) vertex_pentagonFaces_truncatedIcosahedron.size(); i++)
    {
        for(int j=i+1; j<(int) vertex_pentagonFaces_truncatedIcosahedron.size(); j++)
        {
            double angle = vertex_pentagonFaces_truncatedIcosahedron[i].angle(vertex_pentagonFaces_truncatedIcosahedron[j])/deg;
            
            if(i==0 && j==1)
            {
                smallestAngle_pentagonal = angle;
                G4cout << "i: " << i << G4endl;
                G4cout << "j: " << j << G4endl;
                G4cout << G4endl;
            }
            else if(angle<smallestAngle_pentagonal)
            {
                smallestAngle_pentagonal = angle;
                G4cout << "i: " << i << G4endl;
                G4cout << "j: " << j << G4endl;
                G4cout << G4endl;
            }
        }
    }
    
    G4cout << "smallestAngle_pentagonal [deg]: " << smallestAngle_pentagonal << G4endl;
    
    //----------------------------------------------------------------------------
    double smallestAngle_hexagonal_pentagonal = 0.0;
    for(int i=0; i<(int) vertex_hexagonFaces_truncatedIcosahedron.size(); i++)
    {
        for(int j=0; j<(int) vertex_pentagonFaces_truncatedIcosahedron.size(); j++)
        {
            double angle = vertex_hexagonFaces_truncatedIcosahedron[i].angle(vertex_pentagonFaces_truncatedIcosahedron[j])/deg;
            
            if(i==0 && j==0)
            {
                smallestAngle_hexagonal_pentagonal = angle;
            }
            else if(angle<smallestAngle_hexagonal_pentagonal)
            {
                smallestAngle_hexagonal_pentagonal = angle;
            }
        }
    }
    
    G4cout << "smallestAngle_hexagonal_pentagonal [deg]: " << smallestAngle_hexagonal_pentagonal << G4endl;

    //----------------------------------------------------------------------------
    //      Modifying the orientation of the truncated icosahedron
    for(int i=0; i<(int) vertex_allFaces_truncatedIcosahedron.size(); i++)
    {
        std::get<1>(vertex_allFaces_truncatedIcosahedron[i]).rotateY(90.0*deg);
        std::get<1>(vertex_allFaces_truncatedIcosahedron[i]).rotateX(smallestAngle_pentagonal/2.0*deg);
        std::get<1>(vertex_allFaces_truncatedIcosahedron[i]).rotateZ(90.0*deg);
    }
    
    for(int i=0; i<(int) vertex_pentagonFaces_truncatedIcosahedron.size(); i++)
    {
        vertex_pentagonFaces_truncatedIcosahedron[i].rotateY(90.0*deg);
        vertex_pentagonFaces_truncatedIcosahedron[i].rotateX(smallestAngle_pentagonal/2.0*deg);
        vertex_pentagonFaces_truncatedIcosahedron[i].rotateZ(90.0*deg);
        
    }

    for(int i=0; i<(int) vertex_hexagonFaces_truncatedIcosahedron.size(); i++)
    {
        vertex_hexagonFaces_truncatedIcosahedron[i].rotateY(90.0*deg);
        vertex_hexagonFaces_truncatedIcosahedron[i].rotateX(smallestAngle_pentagonal/2.0*deg);
        vertex_hexagonFaces_truncatedIcosahedron[i].rotateZ(90.0*deg);
    }
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    //G4ThreeVector fieldValue = G4ThreeVector(0., 5*tesla, 0.);
    
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);
    
    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetupPreconfiguredVersion(int a)
{
    setPreconfiguredVersion = true;
    LaBr3CeSetupVersion = a;
    
    LaBr3Ce_AllPresent_Override = true;
    LaBr3Ce_AllAbsent_Override = false;
    
    LaBR3Ce_SetGlobalDistance = true;
    LaBR3Ce_automaticOrientation = true;
    configuration_truncatedIcosahedron_hexagons = true;

    if(LaBr3CeSetupVersion==-1)
    {
        //--------------------------------
        //      Work for Christiaan/Mathis/Katarzyna
    }
    else if(LaBr3CeSetupVersion==0)
    {
        //  Full cylindrical crystal (untapered)
        LaBR3Ce_GlobalDistance = 13.2*cm;
    }
    else if(LaBr3CeSetupVersion==1)
    {
        LaBR3Ce_GlobalDistance = 12.0*cm;
    }
    else if(LaBr3CeSetupVersion==2)
    {
        LaBR3Ce_GlobalDistance = 11.0*cm;
    }
    else if(LaBr3CeSetupVersion==3)
    {
        LaBR3Ce_GlobalDistance = 10.0*cm;
    }
    else if(LaBr3CeSetupVersion==4)
    {
        LaBR3Ce_GlobalDistance = 9.0*cm;
    }
    else if(LaBr3CeSetupVersion==5)
    {
        LaBR3Ce_GlobalDistance = 8.0*cm;
    }
    else if(LaBr3CeSetupVersion==6)
    {
        LaBR3Ce_GlobalDistance = 7.0*cm;
    }
    else if(LaBr3CeSetupVersion==7)
    {
        LaBR3Ce_GlobalDistance = 6.0*cm;
    }
    else if(LaBr3CeSetupVersion==8)
    {
        LaBR3Ce_GlobalDistance = 8.9632*cm;
    }
    else if(LaBr3CeSetupVersion==9)
    {
        LaBR3Ce_GlobalDistance = 10.0*cm;
    }
    else if(LaBr3CeSetupVersion==10)
    {
        LaBR3Ce_GlobalDistance = 12.0*cm;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
