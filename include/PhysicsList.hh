//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:
    PhysicsList();
    ~PhysicsList() override;

    void ConstructParticle() override;
    void ConstructProcess() override;

private:

    G4VPhysicsConstructor*  fEmPhysics  = nullptr;
    G4VPhysicsConstructor*  fSyGnPhysics  = nullptr;
    G4VPhysicsConstructor*  fDecayPhysics = nullptr;
    G4VPhysicsConstructor*  fHadPhysics1 = nullptr;
    G4VPhysicsConstructor*  fHadPhysics2 = nullptr;
    G4VPhysicsConstructor*  fHadPhysics3 = nullptr;
    G4VPhysicsConstructor*  fHadPhysics4 = nullptr;
    G4VPhysicsConstructor*  fHadPhysics5 = nullptr;
    G4VPhysicsConstructor*  ftrackingout = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
