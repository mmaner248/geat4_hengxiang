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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in peek volume(not in world)
  if (volume != fScoringVolume) return;
  auto phyvolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  auto copyNo = phyvolume->GetCopyNo();
  /*auto vprocess = step->GetPreStepPoint()->GetProcessDefinedStep();
  if (vprocess) processN = vprocess->GetProcessName();*/
  /*fcellTubs = dynamic_cast<G4Tubs*>(volume->GetSolid());
  G4double cyz = fcellTubs->GetZHalfLength(); // harf length of the tub
  G4double cyr = fcellTubs->GetOuterRadius(); // radius of the tub
  G4double dcyz = 2.*cyz/nz; G4double dcyr = 2.*cyr/nx; // deltaz,deltax*/
  G4double edepStep = step->GetTotalEnergyDeposit();
  if (edepStep <= 0.) return;
  G4ThreeVector prePoint  = step->GetPreStepPoint() ->GetPosition();
  G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector point = prePoint + G4UniformRand()*(postPoint - prePoint);
  if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0) point = postPoint;
  G4double rx = point.getX(); G4double ry = point.getY(); G4double rz = point.getZ();
  fEventAction->AddEdep(edepStep,copyNo);

  /*for (G4int j=0;j<nz;j++){
    for (G4int i=0;i<nx;i++){
      //if (j*dcyr-cyr <= rx && rx <= (j+1)*dcyr-cyr){ // xy plane

        //if (i*dcyr-cyr <= ry && ry <= (i+1)*dcyr-cyr){ // xy plane

      if (j*dcyz-cyz <= rz && rz <= (j+1)*dcyz-cyz){ // xz plane
        if (i*dcyr-cyr <= rx && rx <= (i+1)*dcyr-cyr){ // xz plane
          //G4int k = ny*j+i; // xy plane
          G4int k = nx*j+i; // xz plane
          // collect energy deposited in this step
          fEventAction->AddEdep(edepStep,k);
          if (processN == "eIoni") fEventAction->AddEdepeIoni(edepStep,k);
          if (processN == "eBrem") fEventAction->AddEdepeBrem(edepStep,k);
          if (processN == "msc") fEventAction->AddEdepmsc(edepStep,k);
          if (processN == "compt") fEventAction->AddEdepcompt(edepStep,k);
          if (processN == "phot") fEventAction->AddEdepphot(edepStep,k);

        }

      }

    }
  }*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
