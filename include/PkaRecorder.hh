// PkaRecorder.hh
#ifndef PKA_RECORDER_HH
#define PKA_RECORDER_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <fstream>

class G4Track;
namespace B1 {
    // 简单单例，用来记录 PKA 事件
    class PkaRecorder {
    public:
        static PkaRecorder* Instance();

        // 在 DetectorConstruction 里告诉它体素网格信息
        void InitializeGrid(G4int nx, G4int ny, G4int nz,
            G4double envX, G4double envY, G4double envZ);

        // 在 RunAction::BeginOfRunAction 里打开文件
        void OpenEventFile(const G4String& filename);

        // 在 RunAction::EndOfRunAction 里关闭文件
        void CloseEventFile();

        // 在 TrackingAction 里发现一个 PKA 时调用
        void RecordPka(const G4Track* track);

    private:
        PkaRecorder();
        ~PkaRecorder();

        // 禁用拷贝
        PkaRecorder(const PkaRecorder&) = delete;
        PkaRecorder& operator=(const PkaRecorder&) = delete;

    private:
        // 网格参数（体素个数和整体尺寸）
        G4int    fNx{ 0 }, fNy{ 0 }, fNz{ 0 };
        G4double fEnvX{ 0. }, fEnvY{ 0. }, fEnvZ{ 0. };
        G4double fDx{ 0. }, fDy{ 0. }, fDz{ 0. };

        // 写 PKA 事件的输出文件
        std::ofstream fEventOut;
    };
}
#endif
