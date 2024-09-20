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
// $Id: DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class


#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include <math.h>
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"
#include "MagneticField.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4PropagatorInField.hh"
#include "G4Transform3D.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fETT_FrontPV(nullptr),
   fETT_BackPV(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density;
  G4double fractionalmass;
  G4int nAtoms;

  G4String name;
  G4int ncomponents;

  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();

  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Be");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Sn");
  nistManager->FindOrBuildMaterial("G4_Ta");
  nistManager->FindOrBuildMaterial("G4_W");

  G4Material* O = nistManager->FindOrBuildMaterial("G4_O");
  G4Material* N = nistManager->FindOrBuildMaterial("G4_N");

  //G4Material* H = nistManager->FindOrBuildMaterial("G4_H");
  //G4Material* C = nistManager->FindOrBuildMaterial("G4_C");

  G4Element* em_C = nistManager->FindOrBuildElement(name="C");
  G4Element* em_O = nistManager->FindOrBuildElement(name="O");
  G4Element* em_N = nistManager->FindOrBuildElement(name="N");
  G4Element* em_Ba = nistManager->FindOrBuildElement(name="Ba");
  G4Element* em_F = nistManager->FindOrBuildElement(name="F");
  G4Element* em_Fe = nistManager->FindOrBuildElement(name="Fe");
  G4Element* em_H = nistManager->FindOrBuildElement(name="H");
  G4Element* em_Br = nistManager->FindOrBuildElement(name="Br");
  G4Element* em_Zn = nistManager->FindOrBuildElement(name="Zn");
  G4Element* em_Mn = nistManager->FindOrBuildElement(name="Mn");
  G4Element* em_I  = nistManager->FindOrBuildElement(name="I");

  nistManager->FindOrBuildMaterial("G4_KAPTON");
  nistManager->FindOrBuildMaterial("G4_POLYCARBONATE");
  nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
  nistManager->FindOrBuildMaterial("G4_BRASS");
  nistManager->FindOrBuildMaterial("G4_POLYPROPYLENE");

  G4Material* Air = new G4Material(name="Air",density = 1.290*mg/cm3,ncomponents=2);
  Air->AddMaterial(N, fractionalmass=70*perCent);
  Air->AddMaterial(O, fractionalmass=30*perCent);

  // Image Plates: MS

  G4Material* IP_ProtectiveLayer_Mat = new G4Material(name="IP_ProtectiveLayer_Mat", density = 1.66*g/cm3, ncomponents=3);
  //C2H2O
  IP_ProtectiveLayer_Mat->AddElement(em_C,nAtoms = 2);
  IP_ProtectiveLayer_Mat->AddElement(em_H,nAtoms = 2);
  IP_ProtectiveLayer_Mat->AddElement(em_O,nAtoms = 1);

  G4Material* IP_SensitiveLayer_Mat = new G4Material(name="IP_SensitiveLayer_Mat", density=3.31*g/cm3, ncomponents= 4);

  IP_SensitiveLayer_Mat->AddElement(em_Ba, nAtoms = 20);
  IP_SensitiveLayer_Mat->AddElement(em_F, nAtoms= 20);
  IP_SensitiveLayer_Mat->AddElement(em_Br, nAtoms = 17); 
  IP_SensitiveLayer_Mat->AddElement(em_I, nAtoms = 3); //Need to add I for MS, but fractional 

  G4Material* IP_MagneticLayer_Mat = new G4Material(name="IP_MagneticLayer_Mat", density = 2.77*g/cm3,ncomponents=7);

  IP_MagneticLayer_Mat->AddElement(em_Zn,nAtoms = 1);
  IP_MagneticLayer_Mat->AddElement(em_Mn,nAtoms = 2);
  IP_MagneticLayer_Mat->AddElement(em_Fe,nAtoms = 5);
  IP_MagneticLayer_Mat->AddElement(em_N,nAtoms = 1);
  IP_MagneticLayer_Mat->AddElement(em_O,nAtoms = 40);
  IP_MagneticLayer_Mat->AddElement(em_H,nAtoms = 15);
  IP_MagneticLayer_Mat->AddElement(em_C,nAtoms = 10);


  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                kStateGas, 2.73*kelvin, 3.e-18*pascal);
  
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  // PM placement locations
  G4double PMThickness = 800.*um;
  G4double PMPos_Z = -20*cm-PMThickness;
  G4double PMSize_XY = 5.*cm;
  // Aluminum beam block placement
  G4double aluPos_Z = 60*cm;
  G4double aluThickness = 25*um;
  G4double aluSize_XY = 30*cm;
  // Kapton window placement
  G4double WinPos_Z = 100*cm;
  G4double WinThickness = 125*um;
  G4double WinSize_XY = 30*cm;

  // World Sizes
  G4double worldSize_XY = 100*cm;
  G4double worldSize_Z  = 40*m; 

  // Detector Planes
  G4double ExitDetectorThickness = 0.0025*mm/(2.*100000.);
  G4double FrontDetectorPos = PMPos_Z + PMThickness/2 + ExitDetectorThickness/2;
  //G4double BackDetectorPos = aluPos_Z+aluThickness/2+ExitDetectorThickness/2;
  G4double BackDetectorPos = 105*cm;

  // Get materials

  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  //auto PMMaterial = G4Material::GetMaterial("G4_Sn");
  //auto PMMaterial = G4Material::GetMaterial("G4_Sn");
  //auto PMMaterial = G4Material::GetMaterial("G4_Ta");
  auto PMMaterial = G4Material::GetMaterial("G4_KAPTON");
  auto FlangeMaterial = G4Material::GetMaterial("G4_Al");
  //auto FlangeMaterial = G4Material::GetMaterial("G4_Al");
  auto WindowMaterial = G4Material::GetMaterial("G4_KAPTON");
  auto PacketMaterial = G4Material::GetMaterial("G4_POLYPROPYLENE");
  
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSize_XY/2, worldSize_XY/2, worldSize_Z/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Plasma Mirror
  //
  
  /*  
  auto PlasmaMirrorS
    = new G4Box("PM",            // its name
                 PMSize_XY/2, PMSize_XY/2, PMThickness/2); // its size
                         
  auto PlasmaMirrorLV
    = new G4LogicalVolume(
                 PlasmaMirrorS,        // its solid
                 PMMaterial, // its material
                 "PM");          // its name
                                   
  fAbsorberPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., PMPos_Z), // its position
                 PlasmaMirrorLV,       // its logical volume                         
                 "PM",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  */
  
  //
  // Aluminum Beam block
  // 

  auto FlangeS
    = new G4Box("Flange",            // its name
                 aluSize_XY/2, aluSize_XY/2, aluThickness/2); // its size
                         
  auto FlangeLV
    = new G4LogicalVolume(
                 FlangeS,        // its solid
                 FlangeMaterial, // its material
                 "Flange");          // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., aluPos_Z), // its position
                 FlangeLV,       // its logical volume                         
                 "Flange",           // its name
                worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
 
   //
  // Kapton window
  // 

  auto WindowS
    = new G4Box("Window",            // its name
                 WinSize_XY/2, WinSize_XY/2, WinThickness/2); // its size
                         
  auto WindowLV
    = new G4LogicalVolume(
                 WindowS,        // its solid
                 WindowMaterial, // its material
                 "Window");          // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., WinPos_Z), // its position
                 WindowLV,       // its logical volume                         
                 "Window",           // its name
                worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  //
  // Exit Truth Tracker
  //

  auto ettS 
    = new G4Box("ETT",
                  5.0*cm/2, 5.0*cm/2, ExitDetectorThickness/2);

  auto ettLV
    = new G4LogicalVolume(
                  ettS,
                  defaultMaterial,
                  "ETT");

  fETT_FrontPV
    = new G4PVPlacement(
                  0,
                  G4ThreeVector(0., 0., FrontDetectorPos),
                  ettLV,
                  "ETT_F",
                  worldLV,
                  false,
                  0,
                  fCheckOverlaps);

  fETT_BackPV
    = new G4PVPlacement(
                  0,
                  G4ThreeVector(0., 0., BackDetectorPos),
                  ettLV,
                  "ETT_B",
                  worldLV,
                  false,
                  0,
                  fCheckOverlaps);

  //
  // Magnet Box
  //

  auto MagnetS
    = new G4Box("Mag",            // its name
                 10*cm/2, 10*cm/2, 40*cm/2); // its size
                         
  auto MagnetLV
    = new G4LogicalVolume(
                 MagnetS,        // its solid
                 defaultMaterial, // its material
                 "Mag");          // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., 0.), // its position
                 MagnetLV,       // its logical volume                         
                 "Mag",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // Initialize Absorbers
  //

  auto Brass_M= G4Material::GetMaterial("G4_BRASS");
  auto Aluminum_M = G4Material::GetMaterial("G4_Al");
  auto Steel_M = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto PMMA_M = G4Material::GetMaterial("G4_PLEXIGLASS");

  G4int N_Layers = 25;

  G4double AbsorberThickness_x = 7.*cm;
  G4double AbosrberThickness_y = 7.*cm;

  G4String AbsorberName[] = {"Absorber_Layer_1", //1
                           "Absorber_Layer_2", //2
                           "Absorber_Layer_3", //3
                           "Absorber_Layer_4", //4
                           "Absorber_Layer_5", //5
                           "Absorber_Layer_6", //6
                           "Absorber_Layer_7", //7
                           "Absorber_Layer_8", //8
                           "Absorber_Layer_9", //9
                           "Absorber_Layer_10", //10
                           "Absorber_Layer_11", //11
                           "Absorber_Layer_12", //12
                           "Absorber_Layer_13", //13
                           "Absorber_Layer_14", //14
                           "Absorber_Layer_15", //15
                           "Absorber_Layer_16", //16
                           "Absorber_Layer_17", //17
                           "Absorber_Layer_18", //18
                           "Absorber_Layer_19", //19
                           "Absorber_Layer_20", //20
                           "Absorber_Layer_21", //21
                           "Absorber_Layer_22", //22
                           "Absorber_Layer_23",  //23
                           "Absorber_Layer_24", //24
                           "Absorber_Layer_25", //25
                           "Null" };

  G4String IPNames[] = {"IP_Layer_1",  //1
                        "IP_Layer_2",  //2
                        "IP_Layer_3",  //3
                        "IP_Layer_4",  //4
                        "IP_Layer_5",  //5
                        "IP_Layer_6",  //6
                        "IP_Layer_7",  //7
                        "IP_Layer_8",  //8
                        "IP_Layer_9",  //9
                        "IP_Layer_10", //10
                        "IP_Layer_11", //11
                        "IP_Layer_12", //12
                        "IP_Layer_13", //13
                        "IP_Layer_14", //14
                        "IP_Layer_15", //15
                        "IP_Layer_16", //16
                        "IP_Layer_17", //17
                        "IP_Layer_18", //18
                        "IP_Layer_19", //19
                        "IP_Layer_20", //20
                        "IP_Layer_21", //21
                        "IP_Layer_22", //22
                        "IP_Layer_23", //23
                        "IP_Layer_24", //24
                        "IP_Layer_25", //25
                        "Null" };

  G4double AbsoLayerThickness[] = { 2*mm,  //1
                                    2*mm,  //2
                                    3*mm,  //3
                                    3*mm,  //4
                                    5*mm,  //5
                                    5*mm,  //6
                                    5*mm,  //7
                                    3*mm,  //8
                                    3*mm,  //9
                                    4*mm,  //10
                                    4*mm,  //11
                                    2*mm,  //12
                                    2*mm,  //13
                                    3*mm,  //14
                                    3*mm,  //15
                                    3*mm,  //16
                                    3*mm,  //17
                                    3*mm,  //18
                                    3*mm,  //19
                                    3*mm,  //20
                                    3*mm,  //21
                                    4*mm,  //22
                                    13*mm, //23
                                    10*mm, //24
                                    5*mm,  //25
                                    0.0*mm};


  G4Material* AbsoLayerMat[] = {PMMA_M,     //1
                                PMMA_M,     //2
                                PMMA_M,     //3
                                PMMA_M,     //4
                                PMMA_M,     //5
                                PMMA_M,     //6
                                PMMA_M,     //7
                                Aluminum_M, //8
                                Aluminum_M, //9
                                Aluminum_M, //10
                                Aluminum_M, //11
                                Brass_M   , //12
                                Brass_M   , //13
                                Brass_M   , //14
                                Brass_M   , //15
                                Brass_M   , //16
                                Brass_M   , //17
                                Brass_M   , //18
                                Steel_M   , //19
                                Steel_M   , //20
                                Steel_M   , //21
                                Steel_M   , //22
                                Steel_M   , //23
                                Steel_M   , //24
                                Steel_M   , //25
                                defaultMaterial};


  //
  // Image plate stuff
  //

  auto IP_ProtectiveLayer_M = G4Material::GetMaterial("IP_ProtectiveLayer_Mat");
  auto IP_SensitiveLayer_M = G4Material::GetMaterial("IP_SensitiveLayer_Mat");
  auto IP_MagneticLayer_M = G4Material::GetMaterial("IP_MagneticLayer_Mat");

  G4double ImagePlatePos_Z;
  G4double ImagePlateHeight = 5*cm;
  G4double ImagePlateLength = 5*cm;
  G4double LayerThickness_1 = 9.0*um; // Protective Layer
  G4double LayerThickness_2 = 115*um; // Sensitive Layer
  G4double LayerThickness_3 = 190*um; // "Support" Layer (Protective Material)
  G4double LayerThickness_4 = 160*um; // Magnetic Layer

  G4double LayerThickness_0 = 75*um;
  G4double LayerThickness_5 = 75*um;


  // Geometry
  auto IP_Layer0_S
    = new G4Box("IP_Layer0",
                 ImagePlateLength/2,ImagePlateHeight/2,LayerThickness_0/2);

  auto IP_Layer1_S
    = new G4Box("IP_Layer1",
                 ImagePlateLength/2,ImagePlateHeight/2,LayerThickness_1/2);

  auto IP_Layer2_S
    = new G4Box("IP_Sensitive",
                 ImagePlateLength/2,ImagePlateHeight/2,LayerThickness_2/2);

  auto IP_Layer3_S
    = new G4Box("IP_Layer3",
                 ImagePlateLength/2,ImagePlateHeight/2,LayerThickness_3/2);

  auto IP_Layer4_S
    = new G4Box("IP_Layer4",
                 ImagePlateLength/2,ImagePlateHeight/2,LayerThickness_4/2);


  // Logical Volumes
  auto IP_Layer0_LV
    = new G4LogicalVolume(
                  IP_Layer0_S,
                  PacketMaterial,
                  "IP_Layer0");

  auto IP_Layer1_LV
    = new G4LogicalVolume(
                  IP_Layer1_S,
                  IP_ProtectiveLayer_M,
                  "IP_Layer1");

  auto IP_Layer2_LV
    = new G4LogicalVolume(
                  IP_Layer2_S,
                  IP_SensitiveLayer_M,
                  "IP_Sensitive");

  auto IP_Layer3_LV
    = new G4LogicalVolume(
                  IP_Layer3_S,
                  IP_ProtectiveLayer_M,
                  "IP_Layer3");

  auto IP_Layer4_LV
    = new G4LogicalVolume(
                  IP_Layer4_S,
                  IP_MagneticLayer_M,
                  "IP_Layer4");

  auto IP_Layer5_LV
    = new G4LogicalVolume(
                  IP_Layer0_S,
                  PacketMaterial,
                  "IP_Layer5");

  //
  // Placements 
  //

  G4VSolid* Absorber_Ss[26];
  G4LogicalVolume* Aborbers_LVs[26];

  G4VSolid* StackLayer_S;
  G4LogicalVolume* StackLayer_LV;

  
  for (G4int LayerNo = 0; LayerNo < N_Layers; LayerNo++) {

    StackLayer_S = new G4Box(AbsorberName[LayerNo],
                 AbsorberThickness_x/2,AbosrberThickness_y/2,AbsoLayerThickness[LayerNo]/2);

    Absorber_Ss[LayerNo] = StackLayer_S;

    StackLayer_LV = new G4LogicalVolume(
                  Absorber_Ss[LayerNo],
                  AbsoLayerMat[LayerNo],
                  AbsorberName[LayerNo]);

    Aborbers_LVs[LayerNo] = StackLayer_LV;
  

  }

  for (G4int LayerNo = 0; LayerNo < N_Layers; LayerNo++) {
  G4cout << Absorber_Ss[LayerNo] << G4endl;
  G4cout << Aborbers_LVs[LayerNo] << G4endl;
  G4cout << StackLayer_LV << G4endl;
  }

  

  G4double AbsorberPos_x = 0.0*cm;
  G4double AbsorberPos_y = 0.0*cm;
  G4double AbsorberPos_z;

  G4double IPPos_x = 0.0*cm;
  G4double IPPos_y = 0.0*cm;
  G4double IPPos_z;

  G4double StackStartPos_z = BackDetectorPos + 1.0*cm;
  G4double LayerSpacing = 3*mm;

  G4double offset = 0*cm; 


  for (G4int LayerNo = 0; LayerNo < N_Layers; LayerNo++) {

    G4cout<< LayerNo <<G4endl;

    offset = offset + AbsoLayerThickness[LayerNo]/2;
    AbsorberPos_z = StackStartPos_z + offset;


    new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(AbsorberPos_x, AbsorberPos_y, AbsorberPos_z), // its position
                 Aborbers_LVs[LayerNo],       // its logical volume                         
                 AbsorberName[LayerNo],           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 LayerNo,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    
  offset = offset + AbsoLayerThickness[LayerNo]/2 + LayerSpacing/2;

  ImagePlatePos_Z = StackStartPos_z + offset;

  new G4PVPlacement(
                    0,
                    G4ThreeVector(0.,0.,ImagePlatePos_Z-LayerThickness_0/2),
                    IP_Layer0_LV,
                    "IP_Layer0",
                    worldLV,
                    false,
                    LayerNo,
                    fCheckOverlaps);

  new G4PVPlacement(
                    0,
                    G4ThreeVector(0.,0.,ImagePlatePos_Z+LayerThickness_1/2),
                    IP_Layer1_LV,
                    "IP_Layer1",
                    worldLV,
                    false,
                    LayerNo,
                    fCheckOverlaps);

  new G4PVPlacement(
                    0,
                    G4ThreeVector(0.,0.,ImagePlatePos_Z+LayerThickness_1+LayerThickness_2/2),
                    IP_Layer2_LV,
                    "IP_Sensitive",
                    worldLV,
                    false,
                    LayerNo,
                    fCheckOverlaps);

  new G4PVPlacement(
                    0,
                    G4ThreeVector(0.,0.,ImagePlatePos_Z+LayerThickness_1+LayerThickness_2+LayerThickness_3/2),
                    IP_Layer3_LV,
                    "IP_Layer3",
                    worldLV,
                    false,
                    LayerNo,
                    fCheckOverlaps);

  new G4PVPlacement(
                    0,
                    G4ThreeVector(0.,0.,ImagePlatePos_Z+LayerThickness_1+LayerThickness_2+LayerThickness_3+LayerThickness_4/2),
                    IP_Layer4_LV,
                    "IP_Layer4",
                    worldLV,
                    false,
                    LayerNo,
                    fCheckOverlaps);

  new G4PVPlacement(
                    0,
                    G4ThreeVector(0.,0.,ImagePlatePos_Z+LayerThickness_1+LayerThickness_2+LayerThickness_3+LayerThickness_4+LayerThickness_5/2),
                    IP_Layer5_LV,
                    "IP_Layer5",
                    worldLV,
                    false,
                    LayerNo,
                    fCheckOverlaps);


   offset = offset + LayerSpacing/2;
    


  }

  G4cout<< "test" <<G4endl;
  


  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  //ettLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);

  //
  // Always return the physical World
  //
  return worldPV;
}

G4ThreadLocal MagneticField* DetectorConstruction::fMagneticField;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()

{


  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  fMagneticField = new MagneticField();
  fieldMgr->SetDetectorField(fMagneticField);
  fieldMgr->CreateChordFinder(fMagneticField);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
