#include "PhysicsList.hh"

#include "G4PhysListFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

PhysicsList::PhysicsList(const G4String& refName)
    : G4VModularPhysicsList(), fRefName(refName)
{
    // 1) 先把你原来的“工程参数”保留：verbosity + 默认 cut
    SetVerboseLevel(1);
    SetDefaultCutValue(0.7 * mm);

    // 2) 用工厂拿 Geant4 预定义参考物理表（HP 版）
    G4PhysListFactory factory;

    if (!factory.IsReferencePhysList(fRefName)) {
        G4cerr << "[PhysicsList] Reference physics list <" << fRefName
            << "> not found. Fallback to QGSP_BIC_HP.\n";
        fRefName = "QGSP_BIC_HP";
    }

    fRefPL = factory.GetReferencePhysList(fRefName);
    if (!fRefPL) {
        G4Exception("PhysicsList::PhysicsList", "PL001", FatalException,
            "Failed to create reference physics list.");
    }

    // 3) 让参考物理表也跟随你的 verbosity / cut（cut 我们会在 SetCuts() 再强制一次）
    fRefPL->SetVerboseLevel(GetVerboseLevel());
    fRefPL->SetDefaultCutValue(GetDefaultCutValue());

    G4cout << "[PhysicsList] Using reference physics list: " << fRefName << G4endl;
}

PhysicsList::~PhysicsList()
{
    delete fRefPL;
    fRefPL = nullptr;
}

void PhysicsList::ConstructParticle()
{
    // 参考物理表负责注册所有粒子
    fRefPL->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
    // 参考物理表负责 AddTransportation + EM/Hadronic/Decay 等过程注册
    fRefPL->ConstructProcess();
}

void PhysicsList::SetCuts()
{
    // 关键：保留你现在的 0.7 mm，让新旧结果可比
    fRefPL->SetDefaultCutValue(0.7 * mm);

    // 让参考物理表按它的方式给不同粒子/区域设置 cut
    fRefPL->SetCuts();

    // 如果你想更“硬核”一点，也可以在这里 print 一下 cut 值（可选）
    // DumpCutValuesTable();
}
