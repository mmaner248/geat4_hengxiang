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
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 10;
  fParticleGun  = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  fBeta = particleTable->FindParticle("e-");
  fGamma = particleTable->FindParticle("gamma");
  fAlpha = particleTable->FindParticle("alpha");
  fNeutron = particleTable->FindParticle("neutron");
  fProton = particleTable->FindParticle("proton");


  // default particle kinematic
  fParticleGun->SetParticleDefinition(fGamma);
  fParticleGun->SetParticleEnergy(25.*MeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4double x0 = 0., y0 = 0., z0 = 0.; // half size of worldbox
  G4double x1 = 0., y1 = 0., z1 = 0.; // half size of cellbox
  if (!fworldBox) {
    G4LogicalVolume* worldLV
        = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    if ( worldLV ) fworldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }
  if (fworldBox) {
    x0 = fworldBox->GetXHalfLength();
    y0 = fworldBox->GetYHalfLength();
    z0 = fworldBox->GetZHalfLength();
  }


  x1 = nx*x1,y1 = ny*y1,z1 = nz*z1;
  ftheta = 6.28* CLHEP::rad;
  fphi = 6.28* CLHEP::rad;
  auto theta = (2*G4UniformRand()-1)*ftheta; //(-1,1)*ftheta
  auto phi = (2*G4UniformRand()-1)*fphi; // (-1,1)*fphi
  // G4double xp = x0*(2.*G4UniformRand()-1.); // -10~10
  // G4double yp = y0*(2.*G4UniformRand()-1.); // -10~10
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-0.8*z0));// point source

  G4ThreeVector dir(std::sin(theta) * std::cos(phi),
      std::sin(theta) * std::sin(phi),
      std::cos(theta));
  fParticleGun->SetParticleMomentumDirection(dir);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}


