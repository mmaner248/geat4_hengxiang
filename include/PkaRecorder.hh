// PkaRecorder.hh
#ifndef PKA_RECORDER_HH
#define PKA_RECORDER_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <vector>
class G4Track;
namespace B1 {
    // 简单单例，用来记录 PKA 事件
    //改动，不要用单例模式，每个线程单独一个pkarecorder，这样才能真正并行，不然每个线程都在抢pkarecorder
    class PkaRecorder {
    public:
        static PkaRecorder* Instance();

        // 在 DetectorConstruction 里告诉它体素网格信息
        void InitializeGrid(G4int nx, G4int ny, G4int nz,
            G4double envX, G4double envY, G4double envZ);

   

        // 在 TrackingAction 里发现一个 PKA 时调用
        void RecordPka(const G4Track* track);
        //最后再写的成员函数
        void WriteToFile(const G4String& baseName = "pka");
         int Sizeofpka() const;
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

        // 线程局部的实例指针声明
        static G4ThreadLocal PkaRecorder* fgInstance;

        struct PkaHit {
            G4int eventID = -1;
            G4int trackID = -1;
            G4int Z = 0;
            G4int A = 0;
            G4double Ek_eV = 0.0;
            G4ThreeVector pos;  // mm
            G4ThreeVector dir;  // 方向 
        };
        // 写 PKA 事件的输出文件
        //std::ofstream fEventOut;
        //用vector来缓存pka计算结果，最后再输出
        std::vector<PkaHit> fHits;
    };
}
#endif
