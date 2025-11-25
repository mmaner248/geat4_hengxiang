// TrackingAction.cc
#include "TrackingAction.hh"

#include "G4Track.hh"
#include "PkaRecorder.hh"     // 用来调用 PkaRecorder
#include "G4SystemOfUnits.hh" //用 eV 这种单位
namespace B1 {

    void TrackingAction::PreUserTrackingAction(const G4Track* track)
    {
        // 1. 跳过主粒子（入射中子 / gamma 等）
        if (track->GetParentID() == 0) {
            return;
        }

        // 2. 只关心有原子序的重粒子（排除中子、电子、光子等）
        const auto* pd = track->GetDefinition();
        G4int Z = pd->GetAtomicNumber();
        // G4int A = pd->GetAtomicMass();  // 如需按 A 过滤可以再用，现在先不用

        if (Z <= 0) {
            return;  // Z<=0 一定不是我们要的反冲原子
        }

        // 3. 动能阈值：比如 > 100 eV 才当 PKA（你可以以后再调这个阈值）
        G4double Ek = track->GetKineticEnergy();  // Geant4 内部单位（MeV）
        if (Ek < 100.0 * eV) {
            return;  // 太低的就先不算 PKA
        }

        // 4. 通过以上筛选，这条 track 就被认为是 PKA
        PkaRecorder::Instance()->RecordPka(track);
    }

} // namespace B1
