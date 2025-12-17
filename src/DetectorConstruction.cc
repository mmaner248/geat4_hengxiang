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
#include "G4PhysicalConstants.hh"
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
#include "G4Polyhedra.hh"
#include "PkaRecorder.hh"   
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
  G4bool checkOverlaps = false;

  // ===== World: air, 6 cm cube =====
  auto world_mat = nist->FindOrBuildMaterial("G4_AIR");
  const G4double worldSize = 6.0 * cm;

  auto solidWorld = new G4Box("World", 0.5 * worldSize, 0.5 * worldSize, 0.5 * worldSize);
  auto logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  auto physWorld = new G4PVPlacement(nullptr, {}, logicWorld, "World",
      nullptr, false, 0, checkOverlaps);

  // ===== Graphite block: 1 cm cube (scoring volume) =====
  auto graphite = nist->FindOrBuildMaterial("G4_GRAPHITE");
  const G4double boxSize = 1.0 * cm;

  auto solidBox = new G4Box("GraphiteBox", 0.5 * boxSize, 0.5 * boxSize, 0.5 * boxSize);
  fScoringVolume = new G4LogicalVolume(solidBox, shimo, "GraphiteBoxLV");
  new G4PVPlacement(nullptr, {}, fScoringVolume, "GraphiteBoxPV",
      logicWorld, false, 0, checkOverlaps);

  // ===== Step limit: reduce cross-cell smearing =====
  const G4double dx = boxSize / nx;
  const G4double dy = boxSize / ny;
  const G4double dz = boxSize / nz;
  const G4double maxStep = 0.5 * std::min(dx, std::min(dy, dz));

  fStepLimit = new G4UserLimits(maxStep);
  fScoringVolume->SetUserLimits(fStepLimit);

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
