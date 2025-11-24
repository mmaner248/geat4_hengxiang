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
/// \file B1/include/EventAction.hh
/// \brief Definition of the B1::EventAction class

#ifndef B1EventAction_h
#define B1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "Constants.hh"

namespace B1
{

class RunAction;

/// Event action class

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void AddEdep(G4double edep,G4int copyNo) { fEdep[copyNo] += edep; }
    /*
    void AddEdepeIoni(G4double edep,G4int copyNo) { fEdepeIoni[copyNo] += edep; }
    void AddEdepeBrem(G4double edep,G4int copyNo) { fEdepeBrem[copyNo] += edep; }
    void AddEdepmsc(G4double edep,G4int copyNo) { fEdepmsc[copyNo] += edep; }
    void AddEdepcompt(G4double edep,G4int copyNo) { fEdepcompt[copyNo] += edep; }
    void AddEdepphot(G4double edep,G4int copyNo) { fEdepphot[copyNo] += edep; }
    */

    //void AddEdepold(G4double edep) { fEdepold += edep; }

  private:
    RunAction* fRunAction = nullptr;
    //G4double   fEdepold = 0.; // for test multi thread
	  G4double   fEdep[Cells] = {0.}; // size of array depends on the mesh (plane/body)
    /*G4double   fEdepeIoni[nxz] = {0.};
    G4double   fEdepeBrem[nxz] = {0.};
    G4double   fEdepmsc[nxz] = {0.};
    G4double   fEdepcompt[nxz] = {0.};
    G4double   fEdepphot[nxz] = {0.};*/ // store energy deposit in different physics process
    //std::array<G4double, Cells> fEdep = {0.}; // doesn't work
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


