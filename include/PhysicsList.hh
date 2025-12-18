#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsList : public G4VModularPhysicsList
{
public:
    // 你可以把默认改成 "FTFP_BERT_HP" 也行
    explicit PhysicsList(const G4String& refName = "QGSP_BIC_HP");
    ~PhysicsList() override;

    void ConstructParticle() override;
    void ConstructProcess() override;
    void SetCuts() override;

private:
    G4VModularPhysicsList* fRefPL = nullptr; // 参考物理表实例（我们负责 delete）
    G4String fRefName;
};

#endif

