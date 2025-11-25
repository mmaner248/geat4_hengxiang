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
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters peek
  //
  G4double env_sizeX = 1*cm, env_sizeY = 1*cm, env_sizeZ = 1*cm;
  G4double z,a;
  G4double density;
  G4int ncomponents, natoms;
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"   ,"C" , z= 6., a=  12.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  G4Material* peek =
  new G4Material("peek", density= 1.35*g/cm3, ncomponents=3);
  peek->AddElement(H, natoms = 12);
  peek->AddElement(C, natoms = 19);
  peek->AddElement(O, natoms = 3);
  G4Material* shimo = new G4Material("shimo", density = 1.80 * g / cm3, ncomponents = 1);
  shimo->AddElement(C, natoms = 1);
  //G4Material* fDefaultMaterial = nist->FindOrBuildMaterial("G4_Galactic");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeX = 2.*env_sizeX;
  G4double world_sizeY = 2.*env_sizeY;
  G4double world_sizeZ = 2.*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  //
  // Envelope
  //
  // divide up envlope
  G4double cell_sizeX = env_sizeX/nx;
  G4double cell_sizeY = env_sizeY/ny;
  G4double cell_sizeZ = env_sizeZ/nz;
  // its solid volume
  auto solidcell = new G4Box("solidcell",                    // its name
      0.5 * cell_sizeX, 0.5 * cell_sizeY, 0.5 * cell_sizeZ);  // its size
  /*auto solidcell = new G4Tubs("solidcell",                    // its name
    0., 0.5 * env_sizeX, 0.5 * env_sizeZ, 0., CLHEP::twopi);  // its size*/
  // its logical volume(array will cause core dumped)
  auto logicalcell =
    new G4LogicalVolume(solidcell, shimo, "logicalcell", nullptr, nullptr, nullptr);

  for (G4int iz = 0; iz < nz; iz++) {
    for (G4int copyNo = iz*nxy; copyNo < (iz+1)*nxy; copyNo++) {
      G4int col = (copyNo-iz*nxy)/ny;
      G4int row = (copyNo-iz*nxy)%ny;
      fXCell[copyNo] = (col+1-G4double(nx)/2)*cell_sizeX - cell_sizeX/2;
      fYCell[copyNo] = (row+1-G4double(ny)/2)*cell_sizeY - cell_sizeY/2;
      fZCell[copyNo] = (iz+1-G4double(nz)/2) *cell_sizeZ - cell_sizeZ/2;

      new G4PVPlacement(nullptr,  // no rotation
        G4ThreeVector(fXCell[copyNo], fYCell[copyNo], fZCell[copyNo]),   // at (x,y,z)
        logicalcell,          // its logical volume
        "physicalcell",  // its name
        logicWorld,       // its mother  volume
        false,           // no boolean operations
        copyNo,              // copy number
        checkOverlaps);  // checking overlaps
    }
  }
  G4double maxStep = 0.5*cell_sizeZ;
  // G4double maxStep = 0.5*env_sizeZ/nz;
  fStepLimit = new G4UserLimits(maxStep);
  logicalcell->SetUserLimits(fStepLimit);
  // scoring volume
  fScoringVolume = logicalcell;
  // print cell's location
  //for (G4int i0 = 0; i0 < Cells; i0++){
  //  G4cout<<" copyNo "<<" x "<<fXCell[i0]/cm<<" y "<<fYCell[i0]/cm<<" z "<<fZCell[i0]/cm<<G4endl;
  //}
  //always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
