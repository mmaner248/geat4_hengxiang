// TrackingAction.hh
#ifndef B1_TRACKING_ACTION_HH
#define B1_TRACKING_ACTION_HH

#include "G4UserTrackingAction.hh"

class G4Track;   // 前向声明 Geant4 的 G4Track（全局命名空间）

namespace B1 {

class TrackingAction : public G4UserTrackingAction {
public:
  TrackingAction() = default;
  virtual ~TrackingAction() = default;

  // 每当一条新的 track 诞生时，Geant4 会调用这里
  virtual void PreUserTrackingAction(const G4Track* track) override;
};

} // namespace B1

#endif
