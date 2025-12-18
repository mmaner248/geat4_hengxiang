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
#include "PkaRecorder.hh"


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
    const G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;

    auto* accum = G4AccumulableManager::Instance();
    accum->Merge();

    // grid geometry (virtual)
    const G4double boxSize = 1.0 * cm;
    const G4double dx = boxSize / nx;
    const G4double dy = boxSize / ny;
    const G4double dz = boxSize / nz;
    const G4double volcell = dx * dy * dz; // (internal unit) ~ mm^3

    const auto det = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    const G4double massTotal = det->GetScoringVolume()->GetMass();
    const G4double massCell = massTotal / Cells;

    for (G4int i = 0; i < Cells; i++) {
        const G4double edep = fEdep[i].GetValue(); // MeV
        dose[i] = edep / massCell;                // Gy
        // 如果你还要能量密度 eV/nm^3：
        // const G4double vol_nm3 = (volcell/mm3) * 1e18;
        // eden[i] = (edep/eV) / vol_nm3;
    }

    if (IsMaster()) {
        std::ofstream outfile("shimo.txt");
        outfile << std::setprecision(10);

        for (G4int cellId = 0; cellId < Cells; cellId++) {
            const G4double edepMeV = fEdep[cellId].GetValue() / MeV;

            // 还原 ix/iy/iz（便于后处理画 3D）
            const int iz = cellId / (nx * ny);
            const int rem = cellId - iz * (nx * ny);
            const int iy = rem / nx;
            const int ix = rem - iy * nx;

            outfile
                << ix << " " << iy << " " << iz << " "
                << "Edep(MeV) " << edepMeV << " "
                << "Dose(Gy) " << dose[cellId]
                << "\n";

        }
        outfile.close();
    }
	

    //run完在各个工作线程单例里面把结果写了
    if (G4Threading::IsWorkerThread()) PkaRecorder::Instance()->WriteToFile("pka");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double *edep)
{
  for (G4int copyNo=0;copyNo<Cells;copyNo++) {
  	fEdep[copyNo]  += edep[copyNo];
  	// fEdep2[copyNo] += edep[copyNo]*edep[copyNo];
  }
}

void RunAction::AddEdep(G4double edep, G4int cellId)
{
    if (cellId < 0 || cellId >= Cells) return;
    fEdep[cellId] += edep;   // 直接累加到 accumulable
}

}
