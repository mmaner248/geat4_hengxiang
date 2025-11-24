#include "PhysicsList.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"


#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4NeutronTrackingCut.hh"

#include "G4DecayPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4EmExtraPhysics.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
{

  SetVerboseLevel(1);
  SetDefaultCutValue(0.7*mm);

  // EM physics
  fEmPhysics = new G4EmStandardPhysics_option4();

  // Decay physics
  fDecayPhysics = new G4DecayPhysics(1);

  // Synchroton Radiation & GN Physics
  fSyGnPhysics = new G4EmExtraPhysics();

  // Hadron Physics
  fHadPhysics1 = new G4HadronElasticPhysicsXS();

  fHadPhysics2 = new G4StoppingPhysics();

  fHadPhysics3 = new G4IonPhysicsXS();

  fHadPhysics4 = new G4IonElasticPhysics();

  fHadPhysics5 = new G4HadronInelasticQBBC();

  // Neutron tracking cut
  ftrackingout = new G4NeutronTrackingCut();

}
PhysicsList::~PhysicsList()
{
  delete fSyGnPhysics;
  delete ftrackingout;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor sLivedConstructor;
  sLivedConstructor.ConstructParticle();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysics->ConstructProcess();
  fDecayPhysics->ConstructProcess();
  fSyGnPhysics->ConstructProcess();
  fHadPhysics1->ConstructProcess();
  fHadPhysics2->ConstructProcess();
  fHadPhysics3->ConstructProcess();
  fHadPhysics4->ConstructProcess();
  fHadPhysics5->ConstructProcess();
  ftrackingout->ConstructProcess();
}