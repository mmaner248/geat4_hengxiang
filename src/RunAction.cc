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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Box.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // add new units for dose
  //
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  for (G4int i = 0; i < Cells;i++){
    accumulableManager->RegisterAccumulable(fEdep[i]);
  	// accumulableManager->RegisterAccumulable(fEdep2[i]);
    /*accumulableManager->RegisterAccumulable(fEdepeIoni[i]);
    accumulableManager->RegisterAccumulable(fEdepeBrem[i]);
    accumulableManager->RegisterAccumulable(fEdepmsc[i]);
    accumulableManager->RegisterAccumulable(fEdepcompt[i]);
    accumulableManager->RegisterAccumulable(fEdepphot[i]);*/
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  G4double xcell = 0., ycell = 0., zcell = 0., volcell = 0.; // half size of cellbox
  if (!fcellBox) {
    G4LogicalVolume* cellLV
        = G4LogicalVolumeStore::GetInstance()->GetVolume("logicalcell");
    if ( cellLV ) fcellBox = dynamic_cast<G4Box*>(cellLV->GetSolid());
  }
  if (fcellBox) {
    xcell = fcellBox->GetXHalfLength(); //default mm
    ycell = fcellBox->GetYHalfLength(); //default mm
    zcell = fcellBox->GetZHalfLength(); //default mm
	volcell = fcellBox->GetCubicVolume(); // mm^3
  }

  // Compute dose = total energy deposit in a run and its variance
  // dose and its variance in every element
  for (G4int i = 0; i < Cells; i++)
  {
	  G4double edep = fEdep[i].GetValue(); // default MeV
	  // G4double edep2 = fEdep2[i].GetValue();

	  /*G4double rms = edep2 - edep * edep / nofEvents;
	  if (rms > 0.) rms = std::sqrt(rms);
	  else rms = 0.;*/

	  const auto detConstruction = static_cast<const DetectorConstruction*>(
		  G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	  G4double mass = detConstruction->GetScoringVolume()->GetMass(); // total mass

	  G4double eden = (edep/eV) / (volcell*1e18); // energy density (ev/nm^3)
	  dose[i] = edep / mass; // default gray
	  // rmsdose[i] = rms / mass;
	  pbond1[i] = pow(1. - exp(-eden / epsilonbond1), m); // assume there is only one bond in a cubic nanometer
	  pbond2[i] = pow(1. - exp(-eden / epsilonbond2), m); // hit model
	  pbond3[i] = pow(1. - exp(-eden / epsilonbond3), m);
	  pbond4[i] = pow(1. - exp(-eden / epsilonbond4), m);
	  pbond5[i] = pow(1. - exp(-eden / epsilonbond5), m);
  }
  // Run conditions
  // note: There is no primary generator action object for "master"
  //       run manager for multi-threaded mode.
  /*const auto generatorAction = static_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }*/

  // Print
  // output to an external file
	std::ofstream outfile;
  if (IsMaster()) {
	G4cout<<" zzw " << nofEvents <<G4endl;
  	outfile.open("build\\shimo.txt");
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  	for (G4int copyNo=0;copyNo<Cells;copyNo++) {
    	outfile << std::setprecision(8); // eight significant digits are resevred
       	//<< G4endl
       	//<< " The run consists of " << nofEvents << " "<< runCondition
       	//<< G4endl
		outfile
       	<< " Energy deposited "
       	//<< G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
       	<< fEdep[copyNo].GetValue()/MeV << " Absorbed dose "
       	<< dose[copyNo]   //<< " Probabilityqj "
  		//<< pbond1[copyNo] << " Probabilityty1 "
  		//<< pbond2[copyNo] << " Probabilitybm "
  		//<< pbond3[copyNo] << " Probabilitycc "
  		//<< pbond4[copyNo] << " Probabilityty2 "
  		//<< pbond5[copyNo]
        /*<< fEdepeIoni[copyNo].GetValue()/MeV << " eBrem "
        << fEdepeBrem[copyNo].GetValue()/MeV << " msc "
        << fEdepmsc[copyNo].GetValue()/MeV << " compt "
        << fEdepcompt[copyNo].GetValue()/keV << " phot "
        << fEdepphot[copyNo].GetValue()/keV*/
       	<< G4endl;
  	}
  }
	outfile.close();

  //else {
  //  G4cout
  //   << G4endl
  //   << "--------------------End of Local Run------------------------";
  //}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double *edep)
{
  for (G4int copyNo=0;copyNo<Cells;copyNo++) {
  	fEdep[copyNo]  += edep[copyNo];
  	// fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}/*
void RunAction::AddEdepeIoni(G4double *edep)
{
  for (G4int copyNo=0;copyNo<nxz;copyNo++) {
  	fEdepeIoni[copyNo]  += edep[copyNo];
  	//fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}
void RunAction::AddEdepeBrem(G4double *edep)
{
  for (G4int copyNo=0;copyNo<nxz;copyNo++) {
  	fEdepeBrem[copyNo]  += edep[copyNo];
  	//fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}
void RunAction::AddEdepmsc(G4double *edep)
{
  for (G4int copyNo=0;copyNo<nxz;copyNo++) {
  	fEdepmsc[copyNo]  += edep[copyNo];
  	//fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}
void RunAction::AddEdepcompt(G4double *edep)
{
  for (G4int copyNo=0;copyNo<nxz;copyNo++) {
  	fEdepcompt[copyNo]  += edep[copyNo];
  	//fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}
void RunAction::AddEdepphot(G4double *edep)
{
  for (G4int copyNo=0;copyNo<nxz;copyNo++) {
  	fEdepphot[copyNo]  += edep[copyNo];
  	//fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
