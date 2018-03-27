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
// $Id$
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

//#include "SteppingAction.hh"
// use of stepping action to set the accounting volume

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4VFacet.hh"
#include "G4Isotope.hh"


G4LogicalVolume *Logic_HPGeCrystal_Walid_2;
G4LogicalVolume *Logic_CLOVER_HPGeAluminiumBackingPlate;

void DefineHPGeCrystal_Walid_2()
{
    
    G4cout << "Begin DefineHPGeCrystal_Walid()" << G4endl;
    
    G4double z , a , density;
    G4String name , symbol;
    
    G4double ddxx=(0)*mm;
    
    a = 5.323*g/cm3;
    //density = 5.323*g/cm3;
    //G4Element* elGe = new G4Element(name="Germanium",symbol="Ge",z=32., a);
    
    G4Material* world_ger = new G4Material(name="gerger", z=32.,72.64*g/mole,a );
    // G4Material*  world_ger = man->FindOrBuildMaterial("G4_Ge" , isotopes);
    
    
    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();
    G4NistManager* man = G4NistManager::Instance();
    G4bool isotopes = false;
    // Envelope parameters
    //
    // G4double env_sizeXY = 50*cm, env_sizeZ = 50*cm;
    // G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
    
    //cs source mat
    /*
    //G4Element*  Cs = man->FindOrBuildElement("Cs" , true);
    //G4Element*  Ba = man->FindOrBuildElement("Ba" , true);
    G4Isotope* Cs = new G4Isotope(name="137Cs", 55, 137, 136.907083*g/mole);
    G4Isotope* Ba = new G4Isotope(name="137Ba", 56, 137, 136.9058274*g/mole);
    G4Element* Cs137 = new G4Element("137Cs iso", "Cs", 1);   // pour un age de 29.5 ans
    Cs137->AddIsotope(Cs,  100*perCent);
    //Cs137->AddIsotope(Ba, (100-50.6389)*perCent);
    G4Element* Ba137 = new G4Element("137Cs iso", "Ba", 1);   // pour un age de 29.5 ans
    Ba137->AddIsotope(Ba,  100*perCent);
    
    G4double mde = ((0.506389*136.907083)+((1-0.506389)*136.9058274));
    
    G4Material* Cssourcemat = new G4Material("Cssource", 2.72*g/cm3, 2);
    Cssourcemat ->AddElement(Cs137, 50.6389*perCent);
    Cssourcemat ->AddElement(Ba137, (100-50.6389)*perCent);
    //-------------------------
    //60Co source mat
    
    //G4Element*  Cs = man->FindOrBuildElement("Cs" , true);
    //G4Element*  Ba = man->FindOrBuildElement("Ba" , true);
    G4Isotope* Coo60 = new G4Isotope(name="60Co", 27, 60, 59.9338171*g/mole);
    G4Isotope* Nii60 = new G4Isotope(name="60Ni", 28, 60, 59.9307864*g/mole);
    G4Element* Co60 = new G4Element("60Co iso", "Co", 1);   // pour un age de 29.5 ans
    Co60->AddIsotope(Coo60,  100*perCent);
    //Cs137->AddIsotope(Ba, (100-50.6389)*perCent);
    G4Element* Ni60 = new G4Element("60Ni iso", "Ni", 1);   // pour un age de 29.5 ans
    Ni60->AddIsotope(Nii60,  100*perCent);
    
    //56Fe----------------------------
    G4Isotope* Fee56 = new G4Isotope(name="Fe56", 26, 56, 55.9349375*g/mole);
    G4Element* Fe56 = new G4Element("56Fe iso", "Fe", 1);   // pour un age de 29.5 ans
    Fe56->AddIsotope(Fee56,  100*perCent);
    G4Material* Fe56mat = new G4Material("56Fe mat", 7.8*g/cm3, 1);
    Fe56mat ->AddElement(Fe56, 100*perCent);
    
    //--------------------------------
    
    G4double mde2 = ((0.292357*7.75)+((1-0.292357)*8.908));
    
    G4Material* Cosourcemat = new G4Material("Cosource", mde2*g/cm3, 2);
    Cosourcemat ->AddElement(Ni60, 100*(1-0.292357)*perCent);
    Cosourcemat ->AddElement(Co60, 29.2357*perCent);
    //-------------------------
    
    
    G4Element*  Ge = man->FindOrBuildElement("Ge" , isotopes);
    */
    
    G4Element*  Na = man->FindOrBuildElement("Na" , isotopes);
    G4Element* I = man->FindOrBuildElement("I", isotopes);
    G4Element* Tl = man->FindOrBuildElement("Tl", isotopes);
    G4Element* Bi = man->FindOrBuildElement("Bi", isotopes);
    
    G4Material* NAITL = new G4Material("NaITl", 7.4*g/cm3, 3);
    NAITL ->AddElement(Na, 1);
    NAITL ->AddElement(I, 1);
    NAITL ->AddElement(Tl , 1);
    
    G4Material* NaITl_mat   = nist->FindOrBuildMaterial("NaITl");
    
    
    G4Element*  La = man->FindOrBuildElement("La" , isotopes);
    G4Element* Br = man->FindOrBuildElement("Br", isotopes);
    G4Element* Ce = man->FindOrBuildElement("Ce", isotopes);
    
    G4Material* LaBr3Ce = new G4Material("LaBr3Ce", 7.4*g/cm3, 3);
    LaBr3Ce ->AddElement(La, 1);
    LaBr3Ce ->AddElement(Br, 3);
    LaBr3Ce ->AddElement(Ce , 1);
    
    G4Material* LaBr3Ce_mat   = nist->FindOrBuildMaterial("LaBr3Ce");
    
    //                                              % C 	 % Cr 	       % Ni 	 % Mo 	 % Si % Mn 	 % P 	 % S 	Autres
    // X2CrNiMo17-12-02 1.4404 	Z2CND17-12 	316 L 	0,02 	16 à 18 	10,5 à 13 	2 à 2,5 	1 	2 	0,04 	0,03 	—
    
    G4Element*  C = man->FindOrBuildElement("C" , isotopes);
    G4Element* Cr = man->FindOrBuildElement("Cr", isotopes);
    G4Element* Ni = man->FindOrBuildElement("Ni", isotopes);
    G4Element* Mo = man->FindOrBuildElement("Mo", isotopes);
    G4Element* Fe = man->FindOrBuildElement("Fe", isotopes);
    G4Element* Mn = man->FindOrBuildElement("Mn", isotopes);
    G4Element* Si = man->FindOrBuildElement("Si", isotopes);
    G4Element* P = man->FindOrBuildElement("P", isotopes);
    G4Element* S = man->FindOrBuildElement("S", isotopes);
    G4Element* O = man->FindOrBuildElement("O", isotopes);
    G4Element* Rb = man->FindOrBuildElement("Rb", isotopes);
    G4Element* Cl = man->FindOrBuildElement("Cl", isotopes);
    
    //7.6433936
    G4Material* acierinox = new G4Material("acierinox", 7.85*g/cm3, 9);
    acierinox ->AddElement(C, 0.02*perCent);
    acierinox ->AddElement(Cr, 17*perCent);
    acierinox ->AddElement(Ni , 12*perCent);
    acierinox ->AddElement(Mo , 2.25*perCent);
    acierinox ->AddElement(Fe , 65.66*perCent);
    acierinox ->AddElement(Mn , 2*perCent);
    acierinox ->AddElement(Si , 1*perCent);
    acierinox ->AddElement(P , 0.04*perCent);
    acierinox ->AddElement(S , 0.03*perCent);
    
    /*
    G4Material* BGO = new G4Material("BGO", 7.13*g/cm3, 3);
    BGO ->AddElement(Bi, 4);
    BGO ->AddElement(Ge, 3);
    BGO ->AddElement(O , 12);
    */
    
    G4Material* alu   = nist->FindOrBuildMaterial("G4_Al");
    G4Material* lithium   = nist->FindOrBuildMaterial("G4_Li");
    
    G4Material* tungsten   = nist->FindOrBuildMaterial("G4_W");
    
    G4Material* emptyspace   = nist->FindOrBuildMaterial("G4_Galactic");
    
    G4Material* cuivre   = nist->FindOrBuildMaterial("G4_Cu");
    
    G4Material* natSimat   = nist->FindOrBuildMaterial("G4_Si");
    
    G4Material* natMgmat   = nist->FindOrBuildMaterial("G4_Mg");
    
    // caoutchouc etanchiété C40H44
    G4Element* H = man->FindOrBuildElement("H", isotopes);
    
    G4Material* rubber = new G4Material("caoutchouc etanchiété", 1.0796*g/cm3, 2);
    rubber ->AddElement(C, 40);
    rubber ->AddElement(H, 44);
    
    G4Material* plexiglas = new G4Material("plexiglas", 1.18*g/cm3, 3);
    plexiglas ->AddElement(C, 5);
    plexiglas ->AddElement(H, 8);
    plexiglas ->AddElement(O, 2);
    
    // ruby
    G4Material* matruby = new G4Material("ruby", 4*g/cm3, 2);
    matruby ->AddElement(Rb, 2);
    matruby ->AddElement(O, 1);
    
    //pvc
    G4Material* matpvc = new G4Material("pvc", 1.375*g/cm3, 3);
    matpvc ->AddElement(C, 2);
    matpvc ->AddElement(H, 3);
    matpvc ->AddElement(Cl, 1);
    
    // pet
    G4Material* matpet = new G4Material("pet", 1.38*g/cm3, 3);
    matpet ->AddElement(C, 10);
    matpet ->AddElement(H, 8);
    matpet ->AddElement(O, 4);
    
    // 137Cs
    G4Isotope *s137Cs = new G4Isotope("source 137Cs",55,137,137.0058*g/mole);
    G4Isotope *s137Ba = new G4Isotope("source 137Cs",56,137,137.0058*g/mole);
    
    G4Element* source137Csa = new G4Element("source 137Csa","137Cssource rad",1);
    source137Csa->AddIsotope(s137Cs,100*perCent);
    G4Element* source137Csb = new G4Element("source 137Csb","137Cssource resid",1);
    source137Csb->AddIsotope(s137Ba,100*perCent);
    
    G4Material* mat137Cs = new G4Material("137Cssource at 29y", 137.0058*g/mole, 2);
    mat137Cs->AddElement(source137Csa,50.638915928*perCent);
    mat137Cs->AddElement(source137Csb,(100-50.638915928)*perCent);
    
    
    // 152Eu
    G4Isotope *s152Eu = new G4Isotope("source 152Eu",63,152,151.9217445*g/mole);
    G4Isotope *s152Gd = new G4Isotope("source 152Gd",64,152,151.9197910*g/mole);
    G4Isotope *s152Sm = new G4Isotope("source 152Sm",62,152,151.9197324*g/mole);
    
    G4Element* source152Eua = new G4Element("source 152Eua","152Eusource rad",1);
    source152Eua->AddIsotope(s152Eu,100*perCent);
    G4Element* source152Gda = new G4Element("source 152Gda","152Gdsource rad",1);
    source152Gda->AddIsotope(s152Gd,100*perCent);
    G4Element* source152Sma = new G4Element("source 152Sma","152Smsource rad",1);
    source152Sma->AddIsotope(s152Sm,100*perCent);
    
    G4double masmolEu = (0.6588628*5.244 )+(0.2792*(1-0.6588628)*7.6)+(0.7208*(1-0.6588628)*7.6);
    
    G4Material* mat152Eu = new G4Material("152Eusource at 9.4y", masmolEu*g/cm3, 3);
    mat152Eu->AddElement(source152Eua,(0.6588628*100)*perCent);
    mat152Eu->AddElement(source152Gda,(0.2792*(1-0.6588628)*100)*perCent);
    mat152Eu->AddElement(source152Sma,(0.7208*(1-0.6588628)*100)*perCent);
    
    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = false;
    
    //
    ///////////////////////////////////////// World//////////////////////////////////////
    //
    // G4double world_sizeXY = 4*m;
    // G4double world_sizeZ  = 1.2*env_sizeZ;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
    
    G4Box* solidWorld =
    new G4Box("World",                       //its name
              1*m, 1*m, 1*m);     //its size
    
    G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
    
    G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////      cristal NaI(Tl)   ////////////////////////////////
    
    //G4Tubs * crist_NaI1 = new G4Tubs("cristNaI1",0.*cm,6.35*cm,6.35*cm,0.*deg,360.*deg);
    //G4LogicalVolume * NaIlog1 = new G4LogicalVolume(crist_NaI1, NaITl_mat , "cristNaI1");
    //G4RotationMatrix* yRotNaI1 = new G4RotationMatrix;
    //yRotNaI1->rotateY(-45.*deg);
    //new G4PVPlacement (yRotNaI1,G4ThreeVector(25.*cm,0,25.*cm),NaIlog1,"cristNaI1",logicWorld,false,0,checkOverlaps);
    
    
    //G4Tubs * crist_NaI2 = new G4Tubs("cristNaI2",0.*cm,6.35*cm,6.35*cm,0.*deg,360.*deg);
    //G4LogicalVolume * NaIlog2 = new G4LogicalVolume(crist_NaI2, NaITl_mat , "cristNaI2");
    //G4RotationMatrix* yRotNaI2 = new G4RotationMatrix;
    //yRotNaI2->rotateY(45.*deg);
    //new G4PVPlacement (yRotNaI2,G4ThreeVector(-25.*cm,0,25.*cm),NaIlog2,"cristNaI2",logicWorld,false,0,checkOverlaps);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// cristal LaBr3Cr  //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //G4Tubs * crist_LaBr1 = new G4Tubs("cristLaBr",0.*cm,2.55*cm,2.55*cm,0.*deg,360.*deg);
    //G4LogicalVolume * LaBrlog1 = new G4LogicalVolume(crist_LaBr1, LaBr3Cr_mat , "cristLaBr1");
    //G4RotationMatrix* yRotLaBr1 = new G4RotationMatrix;
    //yRotLaBr1->rotateY(-90.*deg);
    //new G4PVPlacement (yRotLaBr1,G4ThreeVector(35.*cm,0,0.*cm),LaBrlog1,"cristLaBr1",logicWorld,false,0,checkOverlaps);
    
    //G4Tubs * crist_LaBr2 = new G4Tubs("cristLaBr2",0.*cm,1.9*cm,1.9*cm,0.*deg,360.*deg);
    //G4LogicalVolume * LaBrlog2 = new G4LogicalVolume(crist_LaBr2, LaBr3Cr_mat , "cristLaBr2");
    //G4RotationMatrix* yRotLaBr2 = new G4RotationMatrix;
    //yRotLaBr2->rotateY(90.*deg);
    ////yRotLaBr2->rotateY(90.*deg);
    //new G4PVPlacement (yRotLaBr2,G4ThreeVector(-35.*cm,0,0.*cm),LaBrlog2,"cristLaBr2",logicWorld,false,0,checkOverlaps);
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  tete du cristal  Ge  //////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    G4double length_head = (5./tan(7.1*deg)); // 5 = total enlevé vers lavant de la face tapered
    // cout << "length head ==   " << length_head  << endl;
    G4double length_head2 = (4.5/tan(7.1*deg)); // 4.5 = total enlevé vers lavant de la face tapered
    G4double rinner[6] = {0,0,0,0,0,0};
    G4double router[6] = {20.5*mm,22.4134*mm,24.0355*mm,25.1194*mm,25.5*mm,25.5*mm};
    G4double router2[6] = {20.5*mm,23.*mm,24.83012702*mm,25.0*mm,25.0*mm,25.0*mm};
    //G4double posplanZ[5] = {35.*mm,33.75*mm,32.5*mm,31.25*mm,-length_head*mm};
    G4double posplanZ[6] = {(length_head)*mm,(length_head-1.25)*mm,(length_head-2.5)*mm,
								(length_head-3.75)*mm,(length_head-5)*mm,-length_head*mm};
    
    G4Polycone* solidCone = new G4Polycone("cloverCone", 0.0*degree, 360.0*degree,
                                           6,
                                           posplanZ,
                                           rinner,
                                           router);
    
    G4Polycone* solidCone2 = new G4Polycone("cloverCone2", 0.0*degree, 360.0*degree,
                                            6,
                                            posplanZ,
                                            rinner,
                                            router2);
    
    
    G4Tubs* tubtt =
    new G4Tubs("tube",
               0.*mm, 25.5*mm,length_head*mm,0.*deg, 360.*deg);
    
    
    G4Box* boxtt =
    new G4Box("prllppd",                      //its name
              25.5*mm, 25.5*mm, length_head*mm); //its size
    
    G4Box* boxtt2 =
    new G4Box("prllppd",                      //its name
              25.*mm, 25.*mm, length_head2*mm); //its size
    
    G4ThreeVector xytrans(0.5*cm, 0.5*cm,0);
    G4ThreeVector xytrans22(0.45*cm, 0.45*cm,0);
    // G4ThreeVector xytrans2(0.435948045*cm, 0.435948045*cm,0);
    G4RotationMatrix* rot = new G4RotationMatrix;
    rot->rotateY(0.*rad);
    
    G4IntersectionSolid* subtraction = new G4IntersectionSolid("forme1",tubtt , boxtt, rot,xytrans);
    G4IntersectionSolid* subtraction2 = new G4IntersectionSolid("forme1", solidCone,boxtt, rot,xytrans);
    G4IntersectionSolid* subtraction3 = new G4IntersectionSolid("forme1", solidCone2,boxtt2, rot,xytrans22);
    
    //G4LogicalVolume *logicshape1 = new G4LogicalVolume(subtraction1, world_ger, "tube*prllppd");
    // G4LogicalVolume *logicshape2 = new G4LogicalVolume(subtraction2, world_ger, "tube*prllppd");
    
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(),         //at (0,0,0)
    //logicshape1,                //its logical volume
    //"tube-prllppd1",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(),         //at (0,0,0)
    //logicshape2,                //its logical volume
    //"tube-prllppd2",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    
    //G4Trd* tete = new G4Trd("toto",3.55*cm,3.064051955*cm,3.55*cm,3.064051955*cm,1.75*cm);
    G4Trd* tete = new G4Trd("toto",51*0.5*mm,(41*0.5)*mm,51*0.5*mm,(41*0.5)*mm,(length_head*0.5)*mm);
    
    G4Trd* tete2 = new G4Trd("toto",50*0.5*mm,(41*0.5)*mm,50*0.5*mm,(41*0.5)*mm,(length_head2*0.5)*mm);
    
    
    G4ThreeVector xytrans2(0*cm, 0*cm,(length_head*0.5)*mm);
    G4ThreeVector xytrans3(0*cm, 0*cm,(length_head2*0.5)*mm);
    G4RotationMatrix* rott = new G4RotationMatrix;
    rott->rotateY(0.*rad);
    
    G4IntersectionSolid* crist = new G4IntersectionSolid("forme2",subtraction,tete,rott,xytrans2);
    G4IntersectionSolid* crist2 = new G4IntersectionSolid("forme2",subtraction2,tete,rott,xytrans2);
    G4IntersectionSolid* crist3 = new G4IntersectionSolid("forme2",subtraction3,tete2,rott,xytrans3);
    
    //G4LogicalVolume *logcrist = new G4LogicalVolume(crist2, world_ger, "forme2");
    //// G4LogicalVolume *logcrist2 = new G4LogicalVolume(crist2, world_ger, "forme2");
    
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(0,0,0),         //at (0,0,0)
    //logcrist,                //its logical volume
    //"forme2",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(0,0,0),         //at (0,0,0)
    //logcrist2,                //its logical volume
    //"forme2",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////     reste du cor du cristal    ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    G4Tubs* moitL701 = new G4Tubs("demicyl", 0.*mm, 25.5*mm,(70.1-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc701 = new G4Box("cucub", 25.5*mm, 25.5*mm, (70.1-length_head)*0.5*mm);
    
    G4Tubs* moitL701b = new G4Tubs("demicyl", 0.*mm, 25.*mm,(70.1-length_head2)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc701b = new G4Box("cucub", 25.*mm, 25.*mm, (70.1-length_head2)*0.5*mm);
    
    G4ThreeVector xytrans701(5.*mm, 5.*mm,0.*mm); G4ThreeVector xytrans701b(4.5*mm, 4.5*mm,0.*mm);
    G4RotationMatrix* rot701 = new G4RotationMatrix;
    rot701->rotateY(0.*rad);
    
    G4IntersectionSolid* dem701 = new G4IntersectionSolid("demicyl*cucub",moitL701 ,cuc701, rot701,xytrans701);
    G4IntersectionSolid* dem701b = new G4IntersectionSolid("demicyl*cucub",moitL701b ,cuc701b, rot701,xytrans701b);
    
    
    
    G4Tubs* moitL700 = new G4Tubs("demicyl", 0.*cm, 2.55*cm,(70.0-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc700 = new G4Box("cucub", 2.55*cm, 2.55*cm, (70.0-length_head)*0.5*mm);
    
    G4ThreeVector xytrans700(0.5*cm, 0.5*cm,0.*cm);
    G4RotationMatrix* rot700 = new G4RotationMatrix;
    //rot700->rotateY(0.*rad);
    
    G4IntersectionSolid* dem700 = new G4IntersectionSolid("demicyl*cucub",moitL700 ,cuc700, rot701,xytrans700);
    
    
    
    G4Tubs* moitL699 = new G4Tubs("demicyl", 0.*cm, 2.55*cm,(69.9-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc699 = new G4Box("cucub", 2.55*cm, 2.55*cm, (69.9-length_head)*0.5*mm);
    
    G4ThreeVector xytrans699(0.5*cm, 0.5*cm,0.*cm);
    G4RotationMatrix* rot699 = new G4RotationMatrix;
    // rot699->rotateY(0.*rad);
    
    G4IntersectionSolid* dem699 = new G4IntersectionSolid("demicyl*cucub",moitL699 ,cuc699, rot701,xytrans699);
    
    
    
    G4Tubs* moitL698 = new G4Tubs("demicyl", 0.*cm, 2.55*cm,(69.8-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc698 = new G4Box("cucub", 2.55*cm, 2.55*cm, (69.8-length_head)*0.5*mm);
    
    G4ThreeVector xytrans698(0.5*cm, 0.5*cm,0.*cm);
    G4RotationMatrix* rot698 = new G4RotationMatrix;
    // rot698->rotateY(0.*rad);
    
    G4IntersectionSolid* dem698 = new G4IntersectionSolid("demicyl*cucub",moitL698 ,cuc698, rot701,xytrans698);
    
    
    G4Tubs* moitL695 = new G4Tubs("demicyl", 0.*cm, 2.55*cm,(69.5-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc695 = new G4Box("cucub", 2.55*cm, 2.55*cm, (69.5-length_head)*0.5*mm);
    
    G4ThreeVector xytrans695(0.5*cm, 0.5*cm,0.*cm);
    G4RotationMatrix* rot695 = new G4RotationMatrix;
    //rot695->rotateY(0.*rad);
    
    G4IntersectionSolid* dem695 = new G4IntersectionSolid("demicyl*cucub",moitL695 ,cuc695, rot701,xytrans695);
    
    
    G4Tubs* moitL6994 = new G4Tubs("demicyl", 0.*cm, 2.55*cm,(69.94-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc6994 = new G4Box("cucub", 2.55*cm, 2.55*cm, (69.94-length_head)*0.5*mm);
    
    G4ThreeVector xytrans6994(0.5*cm, 0.5*cm,0.*cm);
    G4RotationMatrix* rot6994 = new G4RotationMatrix;
    //rot695->rotateY(0.*rad);
    
    G4IntersectionSolid* dem6994 = new G4IntersectionSolid("demicyl*cucub",moitL6994 ,cuc6994, rot701,xytrans6994);
    
    
    G4Tubs* moitL7015 = new G4Tubs("demicyl", 0.*cm, 2.55*cm,(70.15-length_head)*0.5*mm,0.*deg, 360.*deg);
    G4Box* cuc7015 = new G4Box("cucub", 2.55*cm, 2.55*cm, (70.15-length_head)*0.5*mm);
    
    G4ThreeVector xytrans7015(0.5*cm, 0.5*cm,0.*cm);
    G4RotationMatrix* rot7015 = new G4RotationMatrix;
    //rot695->rotateY(0.*rad);
    
    G4IntersectionSolid* dem7015 = new G4IntersectionSolid("demicyl*cucub",moitL7015 ,cuc7015, rot701,xytrans7015);
    
    
    
    
    
    
    //G4LogicalVolume *logicdemicyl = new G4LogicalVolume(newdemicyl, world_ger, "demicyl*cucub");
    
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(0.*cm,0.*cm,-1.775*cm),         //at (0,0,0)
    //logicdemicyl,                //its logical volume
    //"demicyl*cucub",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////////// partie vide des diodes ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////partie cylindrik vid//////////////////////////////////////////////
    
    G4Tubs* cylduvid = new G4Tubs("cylvid", 0.*cm, 5.*mm,4*cm,0.*deg, 360.*deg);
    //  G4Sphere* teteduvid = new G4Sphere("tetvid",0.*cm, 0.55*cm, 0.*deg,360*deg,0*deg,90*deg);
    G4Tubs* teteduvid = new G4Tubs("tetevid", 0.*cm, 5.*mm,0.5*0.5*mm,0.*deg, 360.*deg);
    G4ThreeVector videtrans(0*cm, 0.*cm,4*cm);
    G4RotationMatrix* rotvid = new G4RotationMatrix;
    
    G4UnionSolid* partievide = new G4UnionSolid("cylvid+tetvid",cylduvid,teteduvid,rotvid,videtrans);
    
    //G4LogicalVolume* emptyyarea = new G4LogicalVolume(partievide,world_ger,"cylvid+tetvid");
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(0,0,0),         //at (0,0,0)
    //emptyyarea,                //its logical volume
    //"cylvid+tetvid",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    G4ThreeVector videtransGe701(0*cm, 0.*cm,-(40.  +(70.1-length_head) - (55 + 0.1) )*mm);
    G4ThreeVector videtransGe700(0*cm, 0.*cm,-(40.  +(70.0-length_head) - (55 + 0.0) )*mm);
    G4ThreeVector videtransGe699(0*cm, 0.*cm,-(40.  +(69.9-length_head) - (55 - 0.1) )*mm);
    G4ThreeVector videtransGe698(0*cm, 0.*cm,-(40.  +(69.8-length_head) - (55 - 0.2) )*mm);
    G4ThreeVector videtransGe695(0*cm, 0.*cm,-(40.  +(69.5-length_head) - (55 - 0.5) )*mm);
    G4ThreeVector videtransGe6994(0*cm, 0.*cm,-(40.  +(69.94-length_head) - (55 - 0.06) )*mm);
    G4ThreeVector videtransGe7015(0*cm, 0.*cm,-(40.  +(70.15-length_head) - (55 + 0.15) )*mm);
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  partie lithium contact  //////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    
    G4Tubs* cylduvid7015 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5+.15)*0.5*mm,0.*deg, 360.*deg);
    G4Tubs* cylduvid701 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5+.1)*0.5*mm,0.*deg, 360.*deg);
    G4Tubs* cylduvid700 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5+.0)*0.5*mm,0.*deg, 360.*deg);
    G4Tubs* cylduvid699 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5-0.1)*0.5*mm,0.*deg, 360.*deg);
    G4Tubs* cylduvid698 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5-0.2)*0.5*mm,0.*deg, 360.*deg);
    G4Tubs* cylduvid695 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5-0.5)*0.5*mm,0.*deg, 360.*deg);
    G4Tubs* cylduvid6994 = new G4Tubs("cylvid", 0.*cm, 5.5*mm,(55.5-0.06)*0.5*mm,0.*deg, 360.*deg);
    //G4Tubs* cylduvid7015 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.575*mm,0.*deg, 360.*deg);   // anciennes valeurs
    //G4Tubs* cylduvid701 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.55*mm,0.*deg, 360.*deg);
    //G4Tubs* cylduvid700 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.5*mm,0.*deg, 360.*deg);
    //G4Tubs* cylduvid699 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.45*mm,0.*deg, 360.*deg);
    //G4Tubs* cylduvid698 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.4*mm,0.*deg, 360.*deg);
    //G4Tubs* cylduvid695 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.25*mm,0.*deg, 360.*deg);
    //G4Tubs* cylduvid6994 = new G4Tubs("cylvid", 0.5*cm, 5.5*mm,22.47*mm,0.*deg, 360.*deg);
    //  G4Sphere* teteduvid2 = new G4Sphere("tetvid",0.55*cm, 0.6*cm, 0.*deg,360*deg,0*deg,90*deg);  // ancienne tete supposée ronde
    G4Tubs* teteduvid2 = new G4Tubs("tetevid", 0.*cm, 5.*mm,60.*0.5*mm,0.*deg, 360.*deg); // nouvelle tete plate
    
    G4ThreeVector videtrans7015(0*cm, 0.*cm,(0.15-2.75)*mm);
    G4ThreeVector videtrans701(0*cm, 0.*cm,(0.1-2.75)*mm);
    G4ThreeVector videtrans700(0*cm, 0.*cm,-2.75*mm);
    G4ThreeVector videtrans699(0*cm, 0.*cm,(-0.1-2.75)*mm);
    G4ThreeVector videtrans698(0*cm, 0.*cm,(-0.2-2.75)*mm);
    G4ThreeVector videtrans695(0*cm, 0.*cm,(-0.5-2.75)*mm);
    G4ThreeVector videtrans6994(0*cm, 0.*cm,(-0.06-2.75)*mm);
    
    G4RotationMatrix* rotvid2 = new G4RotationMatrix;
    
    G4SubtractionSolid* partievideLi7015 = new G4SubtractionSolid("cylvid+tetvid",cylduvid7015,teteduvid2,rotvid2,videtrans7015);
    G4SubtractionSolid* partievideLi701 = new G4SubtractionSolid("cylvid+tetvid",cylduvid701,teteduvid2,rotvid2,videtrans701);
    G4SubtractionSolid* partievideLi700 = new G4SubtractionSolid("cylvid+tetvid",cylduvid700,teteduvid2,rotvid2,videtrans700);
    G4SubtractionSolid* partievideLi699 = new G4SubtractionSolid("cylvid+tetvid",cylduvid699,teteduvid2,rotvid2,videtrans699);
    G4SubtractionSolid* partievideLi698 = new G4SubtractionSolid("cylvid+tetvid",cylduvid698,teteduvid2,rotvid2,videtrans698);
    G4SubtractionSolid* partievideLi695 = new G4SubtractionSolid("cylvid+tetvid",cylduvid695,teteduvid2,rotvid2,videtrans695);
    G4SubtractionSolid* partievideLi6994 = new G4SubtractionSolid("cylvid+tetvid",cylduvid6994,teteduvid2,rotvid2,videtrans6994);
    
    G4LogicalVolume* emptyyareaLi7015 = new G4LogicalVolume(partievideLi7015,lithium,"partievideLi7015");
    G4LogicalVolume* emptyyareaLi701 = new G4LogicalVolume(partievideLi701,lithium,"partievideLi701");
    G4LogicalVolume* emptyyareaLi700 = new G4LogicalVolume(partievideLi700,lithium,"partievideLi700");
    G4LogicalVolume* emptyyareaLi699 = new G4LogicalVolume(partievideLi699,lithium,"partievideLi699");
    G4LogicalVolume* emptyyareaLi698 = new G4LogicalVolume(partievideLi698,lithium,"partievideLi698");
    G4LogicalVolume* emptyyareaLi695 = new G4LogicalVolume(partievideLi695,lithium,"partievideLi695");
    G4LogicalVolume* emptyyareaLi6994 = new G4LogicalVolume(partievideLi6994,lithium,"partievideLi6994");
    
    
    G4ThreeVector LiTrans7015(0*cm, 0.*cm,-((55.5+0.15)*0.5  +(70.15-length_head) -55.5 -.15) *mm);
    G4ThreeVector LiTrans701(0*cm, 0.*cm,-((55.5+0.1)*0.5  +(70.1-length_head) -55.5 -0.1) *mm);
    G4ThreeVector LiTrans700(0*cm, 0.*cm,-((55.5)*0.5  +(70.0-length_head) -55.5 -0.0) *mm);
    G4ThreeVector LiTrans699(0*cm, 0.*cm,-((55.5-0.1)*0.5  +(69.9-length_head) -55.5 +0.1) *mm);
    G4ThreeVector LiTrans698(0*cm, 0.*cm,-((55.5-0.2)*0.5  +(69.8-length_head) -55.5 +0.2) *mm);
    G4ThreeVector LiTrans695(0*cm, 0.*cm,-((55.5-0.5)*0.5  +(69.5-length_head) -55.5 +0.5) *mm);
    G4ThreeVector LiTrans6994(0*cm, 0.*cm,-((55.5-0.06)*0.5  +(69.94-length_head) -55.5 +0.06) *mm);
    //G4ThreeVector LiTrans7015(0*cm, 0.*cm,-12.575*mm);     // ancienne valeurs
    //G4ThreeVector LiTrans701(0*cm, 0.*cm,-12.55*mm);
    //G4ThreeVector LiTrans700(0*cm, 0.*cm,-12.5*mm);
    //G4ThreeVector LiTrans699(0*cm, 0.*cm,-12.45*mm);
    //G4ThreeVector LiTrans698(0*cm, 0.*cm,-12.4*mm);
    //G4ThreeVector LiTrans695(0*cm, 0.*cm,-12.25*mm);
    //G4ThreeVector LiTrans6994(0*cm, 0.*cm,-12.47*mm);
    
    //new G4PVPlacement(0,                       //no rotation
    //G4ThreeVector(0,0,0),         //at (0,0,0)
    //emptyyarea,                //its logical volume
    //"cylvid+tetvid",              //its name
    //logicWorld,              //its mother  volume
    //false,                   //no boolean operation
    //0,                       //copy number
    //checkOverlaps);          //overlaps checking
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////   somme tete + crist  ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    G4ThreeVector ztrans701(0.*mm, 0.*mm,-(70.1-length_head)*0.5*mm);
    G4ThreeVector ztrans701b(0.*mm, 0.*mm,-(70.1-length_head2)*0.5*mm);
    G4ThreeVector ztrans700(0.*cm, 0.*cm,-(70.-length_head)*0.5*mm);
    G4ThreeVector ztrans699(0.*cm, 0.*cm,-(69.9-length_head)*0.5*mm);
    G4ThreeVector ztrans698(0.*cm, 0.*cm,-(69.8-length_head)*0.5*mm);
    G4ThreeVector ztrans695(0.*cm, 0.*cm,-(69.5-length_head)*0.5*mm);
    G4ThreeVector ztrans6994(0.*cm, 0.*cm,-(69.94-length_head)*0.5*mm);
    G4ThreeVector ztrans7015(0.*cm, 0.*cm,-(70.15-length_head)*0.5*mm);
    
    
    G4RotationMatrix* rotz701 = new G4RotationMatrix;
    rotz701->rotateZ(0.*rad);
    
    G4UnionSolid *cristal701 = new G4UnionSolid("newdemicyl+crist",crist2,dem701,rotz701,ztrans701 );
    G4UnionSolid *cristal701b = new G4UnionSolid("newdemicyl+crist",crist3,dem701b,rotz701,ztrans701b );
    
    G4SubtractionSolid * cristalGe701 = new G4SubtractionSolid("cristal 3g",cristal701, cylduvid, rotz701, videtransGe701);
    
    //G4LogicalVolume *crystal3a = new G4LogicalVolume(cristalGe701, world_ger,"cristal 3 green");
    // G4RotationMatrix* rot3green = new G4RotationMatrix;
    //rot3green->rotateY(-90*deg);
    //rot3green->rotateZ(-90*deg);
    
    //new G4PVPlacement(rot3green, G4ThreeVector(0,0,0), crystal3a, "cristal 3 green", logicWorld, false, 0, checkOverlaps);
    //new G4PVPlacement(rot3green, LiTrans701, emptyyareaLi701, "Li 3 green", crystal3a, false, 31, checkOverlaps);
    
    G4double volumm;
    
    
    
    G4UnionSolid *cristal700 = new G4UnionSolid("newdemicyl+crist",crist2,dem700,rotz701,ztrans700 );
    
    G4SubtractionSolid * cristalGe700 = new G4SubtractionSolid("cristal 2 g",cristal700, cylduvid, rotz701, videtransGe700);
    //volumm = cristalGe700->GetCubicVolume()/cm3;
    //G4cout << "volume =      " << volumm << "  cm3" <<"\n";
    
    //G4UnionSolid *cristal7002 = new G4UnionSolid("newdemicyl+crist",crist2,dem700,rotz701,ztrans700 );
    
    //G4SubtractionSolid * cristalGe7002 = new G4SubtractionSolid("cristal 2 g",cristal7002, partievide, rotz701, videtransGe);
    
    
    //G4LogicalVolume *logcristal7002 = new G4LogicalVolume(cristalGe7002, lithium,"cristal 3 green");
    
    //new G4PVPlacement(0, G4ThreeVector(0,0,0), logcristal7002, "cristal 3 green", logicWorld, false, 0, checkOverlaps);
    
    //G4LogicalVolume *logcristal700 = new G4LogicalVolume(cristalGe700, world_ger,"cristal 3 green");
    
    //new G4PVPlacement(0, G4ThreeVector(0,0,0), logcristal700, "cristal 3 green", logicWorld, false, 0, checkOverlaps);
    
    
    
    
    
    G4UnionSolid *cristal699 = new G4UnionSolid("newdemicyl+crist",crist2,dem699,rotz701,ztrans699 );
    
    G4SubtractionSolid * cristalGe699 = new G4SubtractionSolid("cristal 6 g",cristal699, cylduvid, rotz701, videtransGe699);
    
    G4UnionSolid *cristal698 = new G4UnionSolid("newdemicyl+crist",crist2,dem698,rotz701,ztrans698 );
    
    G4SubtractionSolid * cristalGe698 = new G4SubtractionSolid("cristal 6 g",cristal698, cylduvid, rotz701, videtransGe698);
    
    G4UnionSolid *cristal695 = new G4UnionSolid("newdemicyl+crist",crist2,dem695,rotz701,ztrans695 );
    
    G4SubtractionSolid * cristalGe695 = new G4SubtractionSolid("cristal 6 g",cristal695, cylduvid, rotz701, videtransGe695);
    
    G4UnionSolid *cristal6994 = new G4UnionSolid("newdemicyl+crist",crist2,dem6994,rotz701,ztrans6994 );
    
    G4SubtractionSolid * cristalGe6994 = new G4SubtractionSolid("cristal 6 g",cristal6994, cylduvid, rotz701, videtransGe6994);
    
    G4UnionSolid *cristal7015 = new G4UnionSolid("newdemicyl+crist",crist2,dem7015,rotz701,ztrans7015 );
    
    G4SubtractionSolid * cristalGe7015 = new G4SubtractionSolid("cristal 6 g",cristal7015, cylduvid, rotz701, videtransGe7015);
    
    
    //------------------------------------------------
    //      Standard
    Logic_HPGeCrystal_Walid_2 = new G4LogicalVolume(cristalGe700, world_ger,"cristal2green");
    
    //Logic_HPGeCrystal_Walid = new G4LogicalVolume(cristal699, world_ger,"cristal2green");
    //Logic_HPGeCrystal_Walid = new G4LogicalVolume(crist2, world_ger,"cristal2green");
    //Logic_HPGeCrystal_Walid = new G4LogicalVolume(partievide, world_ger,"cristal2green");
    
    //------------------------------------------------
    //      Test cylinder
    //G4Tubs* test = new G4Tubs("test", 0.0*mm, (50.0/2.0)*mm, (70.0/2.0)*mm, 0.*deg, 360.*deg);
    //Logic_HPGeCrystal_Walid = new G4LogicalVolume(test, world_ger,"cristal2green");
    
    
    //volumm = cristalGe700->GetCubicVolume()/cm3;
    //G4cout << "volume =      " << volumm << "  cm3" <<"\n";
    
    
    //------------------------------------------------
    //      Walid's aluminium backplate
    
    G4Tubs*corners1 = new G4Tubs("rounded corners 3",15.5*mm,30.5*mm,82*mm,0,90*deg);
    G4Box* longpartcap = new G4Box("reste du cap 1",50.5*mm,50.5*mm,81.75*mm);
    
    G4RotationMatrix* rotcorners1 = new G4RotationMatrix;
    rotcorners1->rotateZ(-90*deg);
    
    G4RotationMatrix* rotdoublcorners1 = new G4RotationMatrix;
    rotdoublcorners1->rotateZ(180*deg);
    
    G4UnionSolid* doublcorners1a = new G4UnionSolid("2 quarts 3",corners1,corners1,rotcorners1,G4ThreeVector(-70*mm,0*mm,0));
    G4UnionSolid* doublecorners1b = new G4UnionSolid("4 quarts 3",doublcorners1a,doublcorners1a,rotdoublcorners1,G4ThreeVector(-70*mm,-70*mm,0*mm));
    G4SubtractionSolid* longcap1 = new G4SubtractionSolid("long cap 1",longpartcap,doublecorners1b,0,G4ThreeVector(35*mm,35*mm,0*mm));

    //------------------------------------------------
    G4Tubs*corners2 = new G4Tubs("rounded corners 4",14*mm,30.5*mm,84*mm,0,90*deg);
    G4Box* longpartcap2 = new G4Box("reste du cap 2",49*mm,49*mm,83*mm);
    
    G4RotationMatrix* rotcorners2 = new G4RotationMatrix;
    rotcorners2->rotateZ(-90*deg);
    
    G4RotationMatrix* rotdoublcorners2 = new G4RotationMatrix;
    rotdoublcorners2->rotateZ(180*deg);
    
    G4UnionSolid* doublcorners2a = new G4UnionSolid("2 quarts 4",corners2,corners2,rotcorners2,G4ThreeVector(-70*mm,0*mm,0));
    G4UnionSolid* doublecorners2b = new G4UnionSolid("4 quarts 4",doublcorners2a,doublcorners2a,rotdoublcorners2,G4ThreeVector(-70*mm,-70*mm,0*mm));
    G4SubtractionSolid* longcap2 = new G4SubtractionSolid("long cap 2",longpartcap2,doublecorners2b,0,G4ThreeVector(35*mm,35*mm,0*mm));
    
    
    G4SubtractionSolid * longrestcap = new G4SubtractionSolid("long part cap",longcap1,longcap2,0,G4ThreeVector(0,0,0*mm));

    //------------------------------------------------

    G4Box *alplat = new G4Box("al plate", 0.5*92*mm,0.5*92*mm,0.5*10*mm);
    G4IntersectionSolid *alplat1 = new G4IntersectionSolid("plat taillée",alplat,longcap2,0,G4ThreeVector(0,0,0*mm));
    
    Logic_CLOVER_HPGeAluminiumBackingPlate = new G4LogicalVolume(alplat1,alu,"Logic_CLOVER_HPGeAluminiumBackingPlate");
    
    G4cout << "End DefineHPGeCrystal_Walid()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
