// TrackingAction.cc
#include "TrackingAction.hh"

#include "G4Track.hh"
#include "PkaRecorder.hh"     // 用来调用 PkaRecorder
#include "G4SystemOfUnits.hh" //用 eV 这种单位
namespace B1 {

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // 暂时什么都不干，先把钩子挂上
  // 后面我们会在这里识别 PKA，并调用 PkaRecorder::RecordPka(track)

  (void) track;  // 防止“未使用参数”的编译告警
}

} // namespace B1
