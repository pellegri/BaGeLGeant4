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

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

#include "G4IonTable.hh"

#include "BiRelKin.hh"


////    Distributions
#include "alpha0_0plus_0plus_12049_L0.h"
#include "alpha0_0plus_0plus_14032_L0.h"
#include "alpha0_0plus_0plus_15097_L0.h"
#include "alpha0_1minus_0plus_15097_L1.h"
#include "alpha0_2plus_0plus_11520_L2.h"
#include "alpha0_2plus_0plus_15097_L2.h"
#include "alpha0_3minus_0plus_15097_L3.h"
#include "alpha0_4plus_0plus_15097_L4.h"
#include "alpha0_5minus_0plus_15097_L5.h"
#include "alpha1_2plus_2plus_15097_L0.h"
#include "alpha1_3minus_2plus_15097_L1_L3_L5.h"
#include "alpha1_3minus_2plus_15097_L1.h"
#include "alpha1_3minus_2plus_15097_L3.h"
#include "alpha1_3minus_2plus_15097_L5.h"
#include "p0_2plus_12minus_14926_L1.h"
#include "p0_2plus_12minus_14926_L3.h"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0),
fEventAction(eventAction)
{
    ///////////////////////////////////////////////////////////////
    //          To generate radioactive decay - enabled particles
    ///////////////////////////////////////////////////////////////
    
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
    G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particleDefinition);
    //fParticleGun->SetParticleEnergy(0.511*MeV);
    fParticleGun->SetParticleEnergy(1.332*MeV);
    
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,25.*cm));

    ////========================================================////
    ////    Initialise a pre-calculated angular distribution
    G4String angDist_name = "a0_15097_16O";
    
    initialiseAngDist_interpolated(angDist_name);
    
    ////    Distributions - From AngCor
    initialiseAngCor_alpha0_0plus_0plus_12049_L0();
    initialiseAngCor_alpha0_0plus_0plus_14032_L0();
    initialiseAngCor_alpha0_0plus_0plus_15097_L0();
    initialiseAngCor_alpha0_1minus_0plus_15097_L1();
    initialiseAngCor_alpha0_2plus_0plus_11520_L2();
    initialiseAngCor_alpha0_2plus_0plus_15097_L2();
    initialiseAngCor_alpha0_3minus_0plus_15097_L3();
    initialiseAngCor_alpha0_4plus_0plus_15097_L4();
    initialiseAngCor_alpha0_5minus_0plus_15097_L5();
    initialiseAngCor_alpha1_2plus_2plus_15097_L0();
    initialiseAngCor_alpha1_3minus_2plus_15097_L1_L3_L5();
    initialiseAngCor_alpha1_3minus_2plus_15097_L1();
    initialiseAngCor_alpha1_3minus_2plus_15097_L3();
    initialiseAngCor_alpha1_3minus_2plus_15097_L5();
    initialiseAngCor_p0_2plus_12minus_14926_L1();
    initialiseAngCor_p0_2plus_12minus_14926_L3();
    
    
    
    ////    Decay particle energy
    //G4double decayParticleEnergy = (11.520-7.16192)*(12.0/16.0);
    //fParticleGun->SetParticleEnergy(decayParticleEnergy*MeV);

    
    //----------------------------------------------------
    particleN = 0;
    nEnergies = 0;
    
    //----------------------------------------------------
    /*
    nEnergies += 250;
    nParticlesPerEnergy = 2000000; // 502000000 particles to be simulated

    double energyMin = 0.0;
    double energyMax = 25.0;
    double energyDivision = ((energyMax-energyMin)/nEnergies);
    
    for(int i=0; i<nEnergies; i++)
    {
        initialKineticEnergies.push_back(((i+1)*energyDivision)*MeV);
    }
    
    //----------------------------------------------------
    //  The "+1" is for the specific 1.332 MeV gamma ray simulation
    initialKineticEnergies.push_back(1.332*MeV);
    nEnergies++;
    */
    
    
    //----------------------------------------------------
    nParticlesPerEnergy = 5000000;
    
    initialKineticEnergies.push_back(100.0*keV);
    initialKineticEnergies.push_back(133.2*keV);
    initialKineticEnergies.push_back(200.0*keV);
    initialKineticEnergies.push_back(400.0*keV);
    initialKineticEnergies.push_back(700.0*keV);
    initialKineticEnergies.push_back(1.0*MeV);
    initialKineticEnergies.push_back(1.332*MeV);
    initialKineticEnergies.push_back(2.0*MeV);
    initialKineticEnergies.push_back(4.0*MeV);
    initialKineticEnergies.push_back(7.0*MeV);
    initialKineticEnergies.push_back(10.0*MeV);
    initialKineticEnergies.push_back(13.32*MeV);
    initialKineticEnergies.push_back(20.0*MeV);

    //- gamma energies are 100, 133.2, 200, 400, 700, 1000, 1332, 2000, 4000, 7000, 10000, 13320, 20000 keV

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    
    ////////////////////////////////////////////////////////////////
    ////    Definition of states/resonances for 9Be(3He, t)9B
    /*
    G4int noOfStates_9Be = 5;
    G4double resonanceEnergy_9B[noOfStates_9Be]; // MeV
    G4double resonanceWidth_9B[noOfStates_9Be]; // MeV
    
    ////
    resonanceEnergy_9B[0] = 1.5;
    resonanceWidth_9B[0] = 1.2;
    
    ////    5/2-
    resonanceEnergy_9B[1] = 2.345;
    resonanceWidth_9B[1] = 0.081;
    
    ////    5/2+
    resonanceEnergy_9B[2] = 2.751;
    resonanceWidth_9B[2] = 0.614;
    
    ////    1/2-
    resonanceEnergy_9B[3] = 2.78;
    resonanceWidth_9B[3] = 3.13;
    
    ////
    resonanceEnergy_9B[4] = 4.8;
    resonanceWidth_9B[4] = 1.2;
    */
    
    
    ////////////////////////////////////////////////////////////////
    ////    Definition of states/resonances for 9Be(3He, t)9B
    /*
    G4int noOfStates_9Be = 5;
    G4double resonanceEnergy_9B[noOfStates_9Be]; // MeV
    G4double resonanceWidth_9B[noOfStates_9Be]; // MeV
    
    ////
    resonanceEnergy_9B[0] = 1.5;
    resonanceWidth_9B[0] = 1.2;
    
    ////    5/2-
    resonanceEnergy_9B[1] = 2.345;
    resonanceWidth_9B[1] = 0.081;
    
    ////    5/2+
    resonanceEnergy_9B[2] = 2.751;
    resonanceWidth_9B[2] = 0.614;
    
    ////    1/2-
    resonanceEnergy_9B[3] = 2.78;
    resonanceWidth_9B[3] = 3.13;
    
    ////
    resonanceEnergy_9B[4] = 4.8;
    resonanceWidth_9B[4] = 1.2;
    */
    
    
    
    ////    Different particles for Visualisation
    /*
    double test = G4UniformRand();
    
    G4ParticleDefinition* particleDefinitionNew;
    if(test<0.3) particleDefinitionNew = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    if(test>=0.3 && test<0.6) particleDefinitionNew = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    if(test>=0.6 && test<=1.0) particleDefinitionNew = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    fParticleGun->SetParticleDefinition(particleDefinitionNew);
    */
    
    ///////////////////////////////////////////////////////////////
    //          To generate radioactive decay - enabled particles
    ///////////////////////////////////////////////////////////////
    
    /*
    //G4int Z = 3, A = 11;
    //G4int Z = 6, A = 21;
    //G4int Z = 10, A = 18;
    //G4int Z = 11, A = 22;
    //G4int Z = 27, A = 60;
    G4int Z = 8, A = 16;
    //G4int Z = 63, A = 152;
     
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 12.049*MeV;
     
    //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    //ion->SetPDGLifeTime(1*ns);
    ion->SetPDGLifeTime(0.5485101634774202e-10*ns);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
    */

    
    ///////////////////////////////////////////////////////////
    //       Initial Position Distribution of the Particle
    ///////////////////////////////////////////////////////////
    /*
    G4double targetZoffset = -1000.0; // um
    G4double targetThickness = 2.42; // um
    G4double InitialPosition = G4RandFlat::shoot(-(targetThickness/2) + targetZoffset, (targetThickness/2) + targetZoffset); // um
    */
    
    /*
    G4double targetZoffset = 0.0; // um
    G4double targetThickness = 1.0; // um
    G4double InitialPosition = G4RandFlat::shoot(-(targetThickness/2) + targetZoffset, (targetThickness/2) + targetZoffset); // um
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,InitialPosition*um));
    */
    
    ///////////////////////////////////////////////////////////
    //       Initial Direction of the Particle
    ///////////////////////////////////////////////////////////
    /*
    G4double theta = 2*M_PI*G4UniformRand();
    
    mz = -1.0 + 2*G4UniformRand();
    
    ////    Backward Angles
    //mz = -1.0 + G4UniformRand();

    ////    Forward Angles
    //mz = 1.0 - G4UniformRand();

    ////    Very Forward Angles
    //mz = 0.001 - 0.001*G4UniformRand() + (1.0-0.001);

    G4double a = sqrt(1-(mz*mz));
    
    mx = a*cos(theta);
    my = a*sin(theta);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    */

    ////////////////////////////////////////////////////
    //      Particle Energy
    
    /*
    //--------------------------------------------------------------------
    G4double initialParticleKineticEnergy = 0.0*MeV;
    
    int energyN = (int) GetParticleN()/nParticlesPerEnergy;

    if(energyN>=(int initialKineticEnergies.size()))
    {
        initialParticleKineticEnergy = 30.0*MeV;
    }
    else
    {
        initialParticleKineticEnergy = initialKineticEnergies[energyN];
    }
        
    fParticleGun->SetParticleEnergy(initialParticleKineticEnergy);
    fEventAction->SetInitialParticleKineticEnergy(initialParticleKineticEnergy);
    */
    
    //--------------------------------------------------------------------
    G4double initialParticleKineticEnergy = 0.0*MeV;
    
    int energyN = (int) GetParticleN()/nParticlesPerEnergy;
    
    /*
    if(energyN>=(int initialKineticEnergies.size()))
    {
        initialParticleKineticEnergy = 30.0*MeV;
    }
    else
    {
        initialParticleKineticEnergy = initialKineticEnergies[energyN];
    }
    */
    
    if(energyN<((int) initialKineticEnergies.size()))
    {
        initialParticleKineticEnergy = initialKineticEnergies[energyN];
    }
    else
    {
        initialParticleKineticEnergy = 30.0*MeV;
    }
    
    fParticleGun->SetParticleEnergy(initialParticleKineticEnergy);
    fEventAction->SetInitialParticleKineticEnergy(initialParticleKineticEnergy);
    
    //--------------------------------------------------------------------
    /*
    fParticleGun->SetParticleEnergy(1.0*MeV);
    fEventAction->SetInitialParticleKineticEnergy(1.0*MeV);
    */
    

    ////////////////////////////////////////////////////
    ////    ISOTROPIC - Inverse Transform Method
    
    G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
    //G4double theta = acos(1 - (1.0*G4UniformRand()))/deg; // 0.0->90.0 deg
    //G4double theta = acos(1 - (1.0*G4UniformRand() + 1.0))/deg; // 90.0->180.0 deg

    ////    Collimator - 2.0 deg range
    //G4double theta = acos(1 - (6.09173503004822869e-04*G4UniformRand()))/deg; // deg
    
    G4double phi = 360.0*G4UniformRand(); // deg
    //G4double phi = 180.0*G4UniformRand(); // deg
    
    mx = sin(theta*deg)*cos(phi*deg);
    my = sin(theta*deg)*sin(phi*deg);
    mz = cos(theta*deg);
    
    G4ThreeVector direction_gamma0(mx, my, mz);
    //G4ThreeVector direction_gamma1 = -direction_gamma0;
    
    fParticleGun->SetParticleMomentumDirection(direction_gamma0);
    fParticleGun->GeneratePrimaryVertex(anEvent); // This generates a particle vertex (essentially produces the particle with all the previous definitons given to fParticleGun)
    
    
    //fParticleGun->SetParticleMomentumDirection(direction_gamma1);
    //fParticleGun->GeneratePrimaryVertex(anEvent); // This generates a particle vertex (essentially produces the particle with all the previous definitons given to fParticleGun)
    
    
    //--------------------------------------------------------
    //      TEST (Vertical Alignment)
    /*
    G4double theta = 135.0; // 0.0->180.0
    G4double phi = 60.0; // deg

    mx = sin(theta*deg)*cos(phi*deg);
    my = sin(theta*deg)*sin(phi*deg);
    mz = cos(theta*deg);

    my = 2.0*G4UniformRand()-1.0;
    
    G4ThreeVector direction_gamma0(mx, my, mz);
    fParticleGun->SetParticleMomentumDirection(direction_gamma0);
    fParticleGun->GeneratePrimaryVertex(anEvent); // This generates a particle vertex (essentially produces the particle with all the previous definitons given to fParticleGun)
    */
    
    //--------------------------------------------------------
    //      TEST (Along constant phi planes)
    /*
    static int particleN = 0;

    //----------------------------
    G4double phi = 0.0; // deg
    G4double theta = 0.0; // deg
    
    G4double CLOVER_phi[16];
    G4double CLOVER_theta[16];
    
    CLOVER_phi[0] = 0*deg;
    CLOVER_theta[0] = 44.9*deg;
    
    CLOVER_phi[1] = 144*deg;
    CLOVER_theta[1] = 44.9*deg;
    
    CLOVER_phi[2] = 72*deg;
    CLOVER_theta[2] = 71.7*deg;
    
    CLOVER_phi[3] = 144*deg;
    CLOVER_theta[3] = 86.7*deg;
    
    CLOVER_phi[4] = 216*deg;
    CLOVER_theta[4] = 71.7*deg;
    
    CLOVER_phi[5] = 288*deg;
    CLOVER_theta[5] = 71.7*deg;
    
    CLOVER_phi[6] = 36*deg;
    CLOVER_theta[6] = 93.3*deg;
    
    CLOVER_phi[7] = 108*deg;
    CLOVER_theta[7] = 93.3*deg;
    
    CLOVER_phi[8] = 180*deg;
    CLOVER_theta[8] = 93.3*deg;
    
    CLOVER_phi[9] = 252*deg;
    CLOVER_theta[9] = 93.3*deg;
    
    CLOVER_phi[10] = 324*deg;
    CLOVER_theta[10] = 108.3*deg;
    
    CLOVER_phi[11] = 36*deg;
    CLOVER_theta[11] = 135.1*deg;
    
    CLOVER_phi[12] = 108*deg;
    CLOVER_theta[12] = 135.1*deg;
    
    CLOVER_phi[13] = 180.*deg;
    CLOVER_theta[13] = 150.1*deg;
    
    CLOVER_phi[14] = 252*deg;
    CLOVER_theta[14] = 135.1*deg;
    
    CLOVER_phi[15] = 324*deg;
    CLOVER_theta[15] = 150.1*deg;


    //----------------------------
    for(int i=0; i<16; i++)
    {
        if(particleN>=i*100 && particleN<(i+1)*100)
        {
            phi = CLOVER_phi[i];
            theta = CLOVER_theta[i] + (20.0*G4UniformRand() - 10.0)*deg;
            //theta = CLOVER_theta[i];
        }
    }
    
    particleN++;
    
    mx = sin(theta)*cos(phi);
    my = sin(theta)*sin(phi);
    mz = cos(theta);
    
    //my += 0.2*G4UniformRand()-0.1;

    G4ThreeVector direction_gamma0(mx, my, mz);
    fParticleGun->SetParticleMomentumDirection(direction_gamma0);
    fParticleGun->GeneratePrimaryVertex(anEvent); // This generates a particle vertex (essentially produces the particle with all the previous definitons given to fParticleGun)
    */
    
    
    /*
    int nParticles = 3;
    
    for(int i=0; i<nParticles; i++)
    {
        G4double theta = acos(1 - (2.0*G4UniformRand()))/deg; // 0.0->180.0
        G4double phi = 360.0*G4UniformRand(); // deg
        
        mx = sin(theta*deg)*cos(phi*deg);
        my = sin(theta*deg)*sin(phi*deg);
        mz = cos(theta*deg);
        
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    */
    
    
    
    ////////////////////////////////////////////////////////////
    ////    Distributions - Acceptance/Rejection Method
    
    /*
    G4double range;
    G4double mx, my, mz;
    G4bool acceptance = false;
    
    G4double angDist_theta = 0.0;
    G4double function = 0.0;
    
    while(!acceptance)
    {
        angDist_theta = 180.0*G4UniformRand();
        //angDist_theta = 90.0*G4UniformRand() + 90.0;
        range = 1.001*G4UniformRand();

        //  EvalAngCor_alpha0_2plus_0plus_11520_L2: 0.25
        //  EvalAngCor_alpha0_0plus_0plus_12049_L0: 0.95
        //  EvalAngCor_alpha0_0plus_0plus_15097_L0: 0.95
        //  EvalAngCor_alpha1_2plus_2plus_15097_L0: 0.92

        
        //function = sin(angDist_theta*deg);
        
        //function = EvalAngCor_alpha0_2plus_0plus_11520_L2(angDist_theta)*sin(angDist_theta*deg);
        //function = EvalAngCor_alpha0_0plus_0plus_12049_L0(angDist_theta)*sin(angDist_theta*deg);
        //function = EvalAngCor_alpha0_0plus_0plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);
        //function = EvalAngCor_alpha1_2plus_2plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);

        function = EvalAngCor_p0_2plus_12minus_14926_L1(angDist_theta)*sin(angDist_theta*deg);

        if(range<function)
        {
            acceptance = true;
        }
        
        //G4cout << "Testing" << G4endl;
    }
    
    phi = 360.*G4UniformRand();
    
    mx = sin(angDist_theta*deg)*cos(phi*deg);
    my = sin(angDist_theta*deg)*sin(phi*deg);
    mz = cos(angDist_theta*deg);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    */
    
    
    

    ////    Rangefinder
    /*
    G4double highestValue = 0.0;
    G4double angDist_theta = 0.0;
    G4double function;
    
    for(G4double i=0; i<180.0; i+=0.1)
    {
        angDist_theta = i;
        
        //function = EvalAngCor_alpha0_2plus_0plus_11520_L2(angDist_theta)*sin(angDist_theta*deg);
        //function = EvalAngCor_alpha0_0plus_0plus_12049_L0(angDist_theta)*sin(angDist_theta*deg);
        function = EvalAngCor_alpha0_0plus_0plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);
        //function = EvalAngCor_alpha1_2plus_2plus_15097_L0(angDist_theta)*sin(angDist_theta*deg);

        
        if(function>highestValue) highestValue = function;
    
        //G4cout << "EvalAngCor_alpha0_2plus_0plus_11520_L2:    " << function << G4endl;
    }
    
    G4cout << "Here is the highestValue:    " << highestValue << G4endl;
    */

    
    
    ////////////////////////////////////
    ////    Alpha Decay
    ////    Associated Legendre Polynomial
    
    ////    Initialise a pre-calculated angular distribution
    /*
    G4String angDist_name = "a0_15097_16O";
    G4int noOfPoints;
    
    initialiseAngDist_interpolated(angDist_name, noOfPoints);
    G4cout << "noOfPoints:    " << noOfPoints << G4endl;
    */
    
    /*
    G4double range;
    G4double mx, my, mz;
    G4bool acceptance = false;
    
    G4double angDist_theta = 0.0;
    G4double function = 0.0;
    
    while(!acceptance)
    {
        
        //angDist_theta = 180.0*G4UniformRand();
        angDist_theta = 90.0*G4UniformRand() + 90.0;
        range = 1.95*G4UniformRand();
        
        //range = 180.*G4UniformRand();

        ////====================================////
        ////        NEW (Bohr's Theorem)

        ////    alpha0 decay: 2+ -> 0+
        //function += pow( (2.0*0.203684)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
        //function += pow( (0.592631)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
        //function *= sin(angDist_theta);
        
        //function += pow( (2.0*0.204898)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
        //function += pow( (0.590203)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
        //function *= sin(angDist_theta);

        
        ////    alpha0 decay: 2+ -> 0+
        //function = 0.0;
        //function += pow( (2.0*0.204898)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
        //function += pow( (0.590203)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
        //function *= sin(angDist_theta);

        function = EvaluateAngDist_interpolated(angDist_theta)*sin(angDist_theta*deg);
        //function = 0.01*angDist_theta;
        
        //if(angDist_theta>90.0) G4cout << "PROBLEM: angDist_theta>90.0" << G4endl;
            
        //function = sin(angDist_theta);
        //function = sin(angDist_theta) + 0.1*angDist_theta;
        
        //function = 0.2*angDist_theta;

        //function = angDist_theta;
        
        //G4cout << "function:    " << function << G4endl;

        
        if(range<function)
        {
            acceptance = true;
        }
    }
    
    phi = 2.0*M_PI*G4UniformRand();

    mx = sin(angDist_theta*deg)*cos(phi*deg);
    my = sin(angDist_theta*deg)*sin(phi*deg);
    mz = cos(angDist_theta*deg);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    */
    

    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////
    ////    Range-finder for the acceptance-rejection method
    
    /*
    double highestValue = 0.0;
    
    for(G4double i=0; i<180.0; i+=0.1)
    {
        
        //angDist_theta = i/10000.0;
        //angDist_theta*=0.0174533;
        
        //function = 0.0;
        
        //function += pow( (2.0*0.203684)*(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*pow(sin(angDist_theta), 2.0), 2.0);
        //function += pow( (0.592631)*(1.0/4.0)*sqrt(5.0/M_PI)*(3*cos(angDist_theta)*cos(angDist_theta) - 1), 2.0);
        //function *= sin(angDist_theta);
        
        
        angDist_theta = i;
        function = EvaluateAngDist_interpolated(angDist_theta);

        //G4cout << "function, " << angDist_theta << ":    " << function << G4endl;
        
        
        if(function>highestValue) highestValue = function;
    }

    G4cout << "Here is the highestValue:    " << highestValue << G4endl;
    */
    
    
    
    ///////////////////////////////////////////////////
    //       Initial Energy Distribution of Particle
    ///////////////////////////////////////////////////

    //G4double InitialEnergy = G4RandGauss::shoot(((12.049-7.16192)*(0.75)), 0.012);
    //G4double InitialEnergy = ((12.049-7.16192)*(0.75)); // MeV
    
    //  Empirical
    //G4double InitialEnergy = G4RandGauss::shoot(2.96058e+00, ((1.5e-3)/2.3548));
    
    //fParticleGun->SetParticleEnergy(InitialEnergy*MeV);

    
    ////////////////////////
    ////    BiRelKin
    /*
    double m[4], T[4], E[4], p[4];
    double ThSCAT, ThSCAT_Recoil, Ex, Phi;
    
    for(G4int i=0; i<4; i++)
    {
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
    ////	PR226: 16O(a, a')
    m[0] = 3.0160293191; // u
    m[1] = 9.0121822; // u
    m[2] = 3.0160492777; // u
    m[3] = 9.0133288; // u
    
    ////    Projectile Energy
    T[0] = 50.0; // MeV
    
    ////    Ejectile Energy
    T[1] = 0.0; // MeV
    
    
    ////    Uniformly distributed ThSCAT
    ThSCAT = G4RandFlat::shoot( -2., 2.); // deg
    */
    
    ////    Gaussian distributed ThSCAT
    //ThSCAT = 2.1;
    //while(pow(ThSCAT*ThSCAT, 0.5)>2.)
    //{
    //    ThSCAT = G4RandGauss::shoot(0., 1.0);
    //}
    //if(abs(ThSCAT)>2.) "while condition has been breached";
    
    
    ////    12.049 MeV 0+, 16O
    //Ex = G4RandGauss::shoot(12.049, (0.012/2.3548));
    
    ////    12.049 MeV 0+, 16O
    //Ex = G4RandGauss::shoot(15.097, (0.166/2.3548));

    
    
    ////////////////////////////////////////////
    ////    Selection of various states
    /*
    int decayPath = G4RandBit::shootBit();
    int decayProductNumber;
    int resonanceIndex;
    G4double Ex, decayProductEnergy;

    //G4ParticleDefinition* decayProduct_alpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    //G4ParticleDefinition* decayProduct_alpha = G4IonTable::GetIonTable()->GetIon(2,4,0.0);

    //G4ParticleDefinition* decayProduct_alpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    //G4ParticleDefinition* decayProduct_alpha = G4ParticleTable::GetParticleTable()->GetIon(2,4,0.0);

    G4ParticleDefinition* decayProduct;
    G4double usedEnergy;
    G4String decayModeName;
    
    if(decayPath==0)
    {
        ////////////////////////////
        ////    9B -> 4He + 5Li

        ////    Selecting the Decay Product of the Binary Decay
        decayProductNumber = G4RandBit::shootBit();
        
        if(decayProductNumber==0)
        {
            decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
            decayModeName = "9B->4He+5Li, particle=5Li";
        }
        if(decayProductNumber==1)
        {
            decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
            decayModeName = "9B->4He+5Li, particle=alpha";
        }
        
        
        ////    Override to use an alpha particle to approximate the energy loss for the 5Li particle
        //decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
        //decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
        fParticleGun->SetParticleDefinition(decayProduct);
        
        ////    Selecting the Resonance Energy
        resonanceIndex = (int) G4RandFlat::shoot( 0., noOfStates_9Be + 0.99999);
        
        
        Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
        while(Ex>20.0)
        {
            Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
        }
        
        
        //Ex = G4RandGauss::shoot(12.049, (0.012/2.3548));
        
        fEventAction->SetRecoilExcitationEnergy(Ex);
        
        ////    Setting the Decay Product Energy
        if(decayProductNumber==0) decayProductEnergy = (4.0/9.0)*Ex; // MeV
        if(decayProductNumber==1) decayProductEnergy = (5.0/9.0)*Ex; // MeV

        fParticleGun->SetParticleEnergy(decayProductEnergy*MeV);
    }

    
    if(decayPath==1)
    {
        ////////////////////////////
        ////    9B -> 8Be + p
        
        ////    Selecting the Decay Product of the Binary Decay
        decayProductNumber = G4RandBit::shootBit();
        
        if(decayProductNumber==0)
        {
            //decayProduct = G4IonTable::GetIonTable()->GetIon(1,1,0.0);
            decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("proton");
            decayModeName = "9B->8Be+p, particle=proton";
        }
        if(decayProductNumber==1)
        {
            //decayProduct = G4IonTable::GetIonTable()->GetIon(4.0,8.0,0.0);
            decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
            decayModeName = "8Be->4He+4He, particle=proton";
        }
        
        
        ////    Override to use an alpha particle to approximate the energy loss for the 5Li particle
        //decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
        //decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
        fParticleGun->SetParticleDefinition(decayProduct);
        
        ////    Selecting the Resonance Energy
        resonanceIndex = (int) G4RandFlat::shoot( 0., noOfStates_9Be + 0.99999);
        
        
        Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
        while(Ex>20.0)
        {
            Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
        }

        fEventAction->SetRecoilExcitationEnergy(Ex);
        
        
        ////    Setting the Decay Product Energy
        if(decayProductNumber==0) decayProductEnergy = (8.0/9.0)*Ex; // MeV
        
        if(decayProductNumber==1)
        {
            decayProductEnergy = (1.0/9.0)*Ex; // MeV
            
            ////    Selecting a direction for the 8Be -> 4He + 4He decay
            mz = -1.0 + 2.0*G4UniformRand();
            a = sqrt(1-(mz*mz));
            mx = a*cos(theta);
            my = a*sin(theta);
            
            G4ThreeVector energy_COM = G4ThreeVector(mx, my, mz);
            energy_COM.setMag(0.120);
            
            ////    Selecting a direction for the kineitc energy of the 8Be decay product
            mz = -1.0 + 2.0*G4UniformRand();
            a = sqrt(1-(mz*mz));
            mx = a*cos(theta);
            my = a*sin(theta);
            
            G4ThreeVector energy_8Be = G4ThreeVector(mx, my, mz);
            energy_8Be.setMag(decayProductEnergy);

            ////    Calculating the possible energy of the
            G4ThreeVector energy_8Be_alpha = 0.5*energy_COM + energy_8Be;
            usedEnergy = energy_8Be_alpha.mag();
            
            fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));

        }
        
        
        if(decayProductNumber==0) fParticleGun->SetParticleEnergy(decayProductEnergy*MeV);
        if(decayProductNumber==1) fParticleGun->SetParticleEnergy(usedEnergy*MeV);
        
    }
    
    
    
    fEventAction->SetDecayMode(decayModeName);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    
    ////    Inelastically scattered alpha
    /*
    decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    //decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
    //decayProduct = G4IonTable::GetIonTable()->GetIon(2,4,0.0);
    fParticleGun->SetParticleDefinition(decayProduct);
    //fParticleGun->SetParticleEnergy(200.*MeV);
    */
    
    /*
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    
    
    
    ////////////////////////
    ////    BiRelKin
    /*
    bool produceEjectile = true;
    bool produceRecoil = false;

    double m[4], T[4], E[4], p[4];
    double ThSCAT, ThSCAT_Recoil, Phi, Theta;
    //double Ex;
    
    for(G4int i=0; i<4; i++)
    {
        T[i] = 0.0;
        E[i] = 0.0;
        p[i] = 0.0;
    }
    
    ////	PR226: 16O(a, a')
    m[0] = 3.0160293191; // u
    m[1] = 9.0121822; // u
    m[2] = 3.0160492777; // u
    m[3] = 9.0133288; // u
    
    ////    Projectile Energy
    T[0] = 50.0; // MeV
    
    ////    Ejectile Energy
    T[1] = 0.0; // MeV
    
    
    ////    Uniformly distributed ThSCAT
    ThSCAT = G4RandFlat::shoot( 0., 2.); // deg
    
    
    ////    Gaussian distributed ThSCAT
    //ThSCAT = 2.1;
    //while(pow(ThSCAT*ThSCAT, 0.5)>2.)
    //{
    //    ThSCAT = G4RandGauss::shoot(0., 1.0);
    //}
    //if(abs(ThSCAT)>2.) "while condition has been breached";
    
    
    ////    12.049 MeV 0+, 16O
    //Ex = G4RandGauss::shoot(12.049, (0.012/2.3548));
    
    
    ////    Executing BiRelKin
    BiRelKin(m, T, E, p, ThSCAT, ThSCAT_Recoil, Ex);
        
    
    ////    Setting Recoil Energy
    if(produceEjectile)
    {
        fParticleGun->SetParticleEnergy(T[2]*MeV);
        Theta = ThSCAT;
    }
    
    if(produceRecoil)
    {
        fParticleGun->SetParticleEnergy(T[3]*MeV);
        Theta = ThSCAT_Recoil;
    }

    
    
    ////    Setting Recoil Momentum Direction
    Phi = G4RandFlat::shoot( 0., 360.);
    
    mx = sin(Theta*deg)*cos(Phi*deg);
    my = sin(Theta*deg)*sin(Phi*deg);
    mz = abs(cos(Theta*deg));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    */
    
    
    
    
    //if(decayProductNumber==1) decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    
    /*
     if(decayPath==0)
     {
     ////    Selecting the Decay Product of the Binary Decay
     decayProductNumber = (int) G4RandFlat::shoot( 0., 1 + 0.99999);
     
     ////    9B -> 4He +5Li
     if(decayProductNumber==0)
     {
     G4ParticleDefinition* decayProduct = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
     fParticleGun->SetParticleDefinition(decayProduct);
     }
     
     if(decayProductNumber==1)
     {
     G4ParticleDefinition* decayProduct = G4IonTable::GetIonTable()->GetIon(3,5,0.0);
     fParticleGun->SetParticleDefinition(decayProduct);
     }
     
     ////    Selecting the Resonance Energy
     resonanceIndex = (int) G4RandFlat::shoot( 0., noOfStates_9Be + 0.99999);
     Ex = CLHEP::RandBreitWigner::shoot(resonanceEnergy_9B[resonanceIndex], resonanceWidth_9B[resonanceIndex]);
     fEventAction->SetRecoilExcitationEnergy(Ex);
     
     ////    Setting the Decay Product Energy
     if(decayProductNumber==0) decayProductEnergy = (5.0/4.0)*Ex; // MeV
     if(decayProductNumber==1) decayProductEnergy = (4.0/5.0)*Ex; // MeV
     
     fParticleGun->SetParticleEnergy(decayProductEnergy*MeV);
     }
     */
    



    
    
    
    /*
    G4double theta_cm, phi_cm;
    
    ////    Isotropic in Centre-Of-Mass Frame
    bool isotropy = true;

    
    theta_cm = G4RandFlat::shoot( 0., 360.); // deg
    //if()
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector( mx, my, mz));
    */
    
    
    
    
    /*
     ////////    4He, +1 charge
     G4int Z = 2, A = 4;
     
     G4double ionCharge   = 1.*eplus;
     G4double excitEnergy = 0.*MeV;
     
     //G4ParticleDefinition* ion = G4ParticleTable::GetParticleTable()->GetIon(Z,A,excitEnergy);
     G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
     //ion->SetPDGLifeTime(1*ns);
     fParticleGun->SetParticleDefinition(ion);
     fParticleGun->SetParticleCharge(ionCharge);
     */
    
    
    ///////////////////////////////////////////////////
    //              TIME DISTRIBUTION                //
    ///////////////////////////////////////////////////
    
    /*
     initialclocktime = G4RandGauss::shoot( 0, 1.5);
     fParticleGun->SetParticleTime(initialclocktime*ns);
     */
    
    
    
    
    ///////////////////////////////////////////////////
    //       Initial Direction of Particle
    ///////////////////////////////////////////////////
    
    //      Flat Distribution
    //mx = G4RandFlat::shoot( -.15, .15);
    //my = G4RandFlat::shoot( -.15, .15);
    //mz = G4RandFlat::shoot( -1, 1);
    
    ////////////////////////////////
    ////    ISOTROPIC
    //      Exploiting the spherical symmetry of normal distributions
    //      Gaussian Distribution
    /*
     mx = G4RandGauss::shoot( 0, 1.);
     my = G4RandGauss::shoot( 0, 1.);
     mz = G4RandGauss::shoot( 0, 1.);
     
     if(mz<0) mz = -mz;
     //mz = abs(mz);
     //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
     //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, 10.));
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, -mz));
     */
    
    ////    Alternative Method for ISOTROPY
    /*
    G4double theta = 2*M_PI*G4UniformRand();
    //G4double mz = -1.0 + 2*G4UniformRand();
    G4double mz = -1.0 + G4UniformRand();
    
    G4double a = sqrt(1-(mz*mz));
    
    mx = a*cos(theta);
    my = a*sin(theta);
    
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
    */
    
    
    /*
     theta = M_PI*G4UniformRand();
     phi = 2*M_PI*G4UniformRand();
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta*rad)*cos(phi*rad), sin(theta*rad)*sin(phi*rad), cos(theta*rad)));
     */
    
    
    
    ///////////////////////////////////////////////////
    //                  16O(a,a')
    ///////////////////////////////////////////////////
    /*
     G4double test = G4UniformRand();
     
     //if(test<=0.5) beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     //if(test>0.5) beamEnergy = G4RandGauss::shoot( 180, 0.5); // MeV
     
     if(test>0.0 && test<=0.1) beamEnergy = G4RandGauss::shoot( 215, 0.5); // MeV
     if(test>0.1 && test<=0.2) beamEnergy = G4RandGauss::shoot( 210, 0.5); // MeV
     if(test>0.2 && test<=0.3) beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     if(test>0.3 && test<=0.4) beamEnergy = G4RandGauss::shoot( 190, 0.5); // MeV
     if(test>0.4 && test<=0.5) beamEnergy = G4RandGauss::shoot( 180, 0.5); // MeV
     if(test>0.5 && test<=0.6) beamEnergy = G4RandGauss::shoot( 185, 0.5); // MeV
     if(test>0.6 && test<=0.7) beamEnergy = G4RandGauss::shoot( 170, 0.5); // MeV
     if(test>0.7 && test<=0.8) beamEnergy = G4RandGauss::shoot( 160, 0.5); // MeV
     if(test>0.8 && test<=0.9) beamEnergy = G4RandGauss::shoot( 150, 0.5); // MeV
     if(test>0.9 && test<=1.0) beamEnergy = G4RandGauss::shoot( 140, 0.5); // MeV
     
     //beamEnergy = G4RandGauss::shoot( 200, 0.5); // MeV
     //beamEnergy = 200; // MeV
     
     recoilExcitationEnergy = 15.097;  // MeV
     alphaSeperationEnergy = 7.16192;    // MeV
     protonSeperationEnergy = 12.1274;  // MeV
     //  Elastically scattered outgoing alpha
     mx = G4RandGauss::shoot( 0, 0.00001);
     my = G4RandGauss::shoot( 0, .05);
     //mx = G4RandGauss::shoot( 0, 0.00001);
     //my = G4RandGauss::shoot( 0, 0.00001);
     //mx = 0;
     //my = 0;
     
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, 10.));
     energy = beamEnergy - recoilExcitationEnergy; // MeV
     fParticleGun->SetParticleEnergy(energy*MeV);
     fParticleGun->GeneratePrimaryVertex(anEvent);
     
     
     //////////////////////////////////////////////////////////////
     //                      16O -> 4He + 12C
     //////////////////////////////////////////////////////////////
     //  Alpha 0 isotropic decay:    15.097 MeV state of 16O -> to Ground state of 12C
     //  Alpha 1 anisotropic decay:  15.097 MeV state of 16O -> to 4.43891 MeV state of 12C
     //  Assuming recoil nucleus takes no kinetic energy from inelastic scatter
     
     G4double pAlphaDecayMode = G4UniformRand();
     
     if(pAlphaDecayMode<=0.5) alphaDecayMode = 0;
     if(pAlphaDecayMode>0.5) alphaDecayMode = 1;
     
     //alphaDecayMode = 1;
     
     ////    ALPHA 1
     if(alphaDecayMode == 0) daughterExcitationEnergy = 0; //    MeV
     
     ////    Alpha 0
     if(alphaDecayMode == 1) daughterExcitationEnergy = 4.43891; //    MeV
     
     mx = G4RandFlat::shoot( -1., 1.);
     my = G4RandFlat::shoot( -1., 1.);
     mz = G4RandFlat::shoot( -1., 1.);
     
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(3/4)*MeV);
     
     fParticleGun->SetParticleMomentumDirection(ejectileDirection);
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
    
    ////    Recoil - 12C
    /*
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(1/4)*MeV);
     
     recoilDirection = -ejectileDirection;
     fParticleGun->SetParticleMomentumDirection(recoilDirection);
     */
    
    /*
     Z = 6, A = 12;
     
     G4ParticleDefinition* recoil = G4IonTable::GetIonTable()->GetIon(Z, A, 4.43891*MeV);
     //G4ParticleDefinition* recoil = G4ParticleTable::GetParticleTable()->GetIon(Z, A, 0.*keV);
     fParticleGun->SetParticleCharge(0.*eplus);
     //recoil->SetPDGLifeTime(1*ns);
     //recoil->SetPDGStable(true);
     fParticleGun->SetParticleDefinition(recoil);
     
     fParticleGun->SetParticleEnergy( (recoilExcitationEnergy - alphaSeperationEnergy - daughterExcitationEnergy)*(1/4)*MeV);
     
     recoilDirection = -ejectileDirection;
     fParticleGun->SetParticleMomentumDirection(recoilDirection);
     */
    
    
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    //G4double Lifetime = fParticleGun->GetParticleDefinition()->GetPDGLifeTime();
    //G4double Lifetime = recoil->GetPDGLifeTime();
    
    //G4cout << "Here is the alphaSeperationEnergy    "<< alphaSeperationEnergy << G4endl;
    
    
    
    ////     Gamma Decay - from Recoil Nucleus
    /*
     if(alphaDecayMode == 1)
     {
     G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
     fParticleGun->SetParticleDefinition(particleDefinition);
     
     fParticleGun->SetParticleEnergy(daughterExcitationEnergy*MeV);
     
     mx = G4RandFlat::shoot( -1., 1.);
     my = G4RandFlat::shoot( -1., 1.);
     mz = G4RandFlat::shoot( -1., 1.);
     fParticleGun->SetParticleMomentumDirection(G4ThreeVector(mx, my, mz));
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     }
     */
    
    

    
    
    
    
    ////////////////////////////////////////////////////////
    //          PARTICLE INFORMATION OUTPUT
    ////////////////////////////////////////////////////////
    /*
     G4double Lifetime = fParticleGun->GetParticleDefinition()->GetPDGLifeTime();
     G4double Spin = fParticleGun->GetParticleDefinition()->GetPDGSpin();
     G4double Isospin = fParticleGun->GetParticleDefinition()->GetPDGIsospin();
     G4double Parity = fParticleGun->GetParticleDefinition()->GetPDGiGParity();
     G4String ParticleName = fParticleGun->GetParticleDefinition()->GetParticleName();
     
     
     
     G4cout << "Here is the Lifetime    "<< Lifetime << G4endl;
     G4cout << "Here is the Spin    "<< Spin << G4endl;
     G4cout << "Here is the Isospin    "<< Isospin << G4endl;
     G4cout << "Here is the Parity    "<< Parity << G4endl;
     G4cout << "Here is the ParticleName    "<< ParticleName << G4endl;
     */
    

    //create vertex
    //
    //fParticleGun->GeneratePrimaryVertex(anEvent);
    
    
    
    /*
     // This function is called at the begining of event
     
     // In order to avoid dependence of PrimaryGeneratorAction
     // on DetectorConstruction class we get world volume
     // from G4LogicalVolumeStore
     //
     G4double worldZHalfLength = 0;
     G4LogicalVolume* worlLV
     = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
     G4Box* worldBox = 0;
     if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid());
     if ( worldBox ) {
     worldZHalfLength = worldBox->GetZHalfLength();
     }
     else  {
     G4ExceptionDescription msg;
     msg << "World volume of box not found." << G4endl;
     msg << "Perhaps you have changed geometry." << G4endl;
     msg << "The gun will be place in the center.";
     G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002", JustWarning, msg);
     }
     
     // Set gun position
     fParticleGun
     ->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
     
     fParticleGun->GeneratePrimaryVertex(anEvent);
     */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

