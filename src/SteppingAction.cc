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
#include "G4SystemOfUnits.hh" //用 eV 这种单位
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "PkaRecorder.hh"
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    
    const G4Track* tr = step->GetTrack();
    //如果按照严格的pka定义，那只考虑中子第一次碰撞的C原子，那就需要加下面一句把track直接限定在入射
    // if (tr->GetParentID() != 0) return;
    // 但是我这里记录的是中子碰撞导致的所有可能引发后续级联的种子
    // 所以先注释掉
    
    // 只处理中子这条track的每一步
    if (tr->GetDefinition()->GetParticleName() != "neutron") return;


    // 只在“有物理过程发生”时看本步产生的二次粒子
    if (!step->GetPostStepPoint()->GetProcessDefinedStep()) return;
    auto* p = step->GetPostStepPoint()->GetProcessDefinedStep();
    auto name = p->GetProcessName();
    if (name == "Transportation" || name == "StepLimiter") return;

    const auto* secs = step->GetSecondaryInCurrentStep();
    if (!secs) return;

    for (const auto* secTr : *secs)
    {
        const auto* pd = secTr->GetDefinition();
        G4int Z = pd->GetAtomicNumber();
        if (Z !=6) continue;                // 只要反冲核
        if (secTr->GetKineticEnergy() < 100 * eV) continue;  // 阈值

        const auto* cproc = secTr->GetCreatorProcess();
        if (!cproc) continue;

        const auto& cname = cproc->GetProcessName();
        if (cname != "hadElastic") continue;

        // 这里的 secTr 就是“由中子这一步产生的第一代反冲”
        PkaRecorder::Instance()->RecordPka(secTr);
    }
    if (!fScoringVolume) {
        const auto detConstruction = static_cast<const DetectorConstruction*>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        fScoringVolume = detConstruction->GetScoringVolume();
    }

    auto volume = step->GetPreStepPoint()->GetTouchableHandle()
        ->GetVolume()->GetLogicalVolume();

    if (volume != fScoringVolume) return;

    const G4double edepStep = step->GetTotalEnergyDeposit();
    if (edepStep <= 0.) return;

    const G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
    const G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();

    G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);
    if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0) point = postPoint;

    // ===== Map world position -> virtual grid cellId inside 1cm graphite cube at origin =====
    const G4double boxSize = 1.0 * cm;
    const G4double half = 0.5 * boxSize;

    const G4double dx = boxSize / nx;
    const G4double dy = boxSize / ny;
    const G4double dz = boxSize / nz;

    const G4double x = point.x() + half;
    const G4double y = point.y() + half;
    const G4double z = point.z() + half;

    if (x < 0 || x >= boxSize || y < 0 || y >= boxSize || z < 0 || z >= boxSize) return;

    int ix = (int)(x / dx);
    int iy = (int)(y / dy);
    int iz = (int)(z / dz);

    // clamp (defensive)
    if (ix < 0) ix = 0; else if (ix >= nx) ix = nx - 1;
    if (iy < 0) iy = 0; else if (iy >= ny) iy = ny - 1;
    if (iz < 0) iz = 0; else if (iz >= nz) iz = nz - 1;

    const int cellId = ix + nx * (iy + ny * iz);

    fEventAction->AddEdep(edepStep, cellId);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
