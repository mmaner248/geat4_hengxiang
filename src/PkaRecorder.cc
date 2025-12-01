// PkaRecorder.cc
#include "PkaRecorder.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"   // Geant4 的自动加锁工具
#include <cmath>           // 用 std::floor 计算体素索引
/*
namespace {
    // 用于保护文件写入的全局互斥锁（多线程用）
    G4Mutex pkaMutex = G4MUTEX_INITIALIZER; 不需要了，我们现在为每个线程创建一个新的ostream，免得线程
    一直在等待IO
}*/
namespace B1 {
    //静态变量在类外初始化，这里为什么要用静态呢，和instance是一个道理，这2个东西是用来产生实例的
    //那肯定不能靠实例来初始化，所以只能是静态
    G4ThreadLocal PkaRecorder* PkaRecorder::fgInstance = nullptr;

    PkaRecorder* PkaRecorder::Instance()
    {
        if (!fgInstance) {
            fgInstance = new PkaRecorder();
        }
        return fgInstance;
    }

    PkaRecorder::PkaRecorder() {}

    PkaRecorder::~PkaRecorder()
    {
        if (fEventOut.is_open()) {
            fEventOut.close();
        }
    }

    void PkaRecorder::InitializeGrid(G4int nx, G4int ny, G4int nz,
        G4double envX, G4double envY, G4double envZ)
    {
        fNx = nx;  fNy = ny;  fNz = nz;
        fEnvX = envX;  fEnvY = envY;  fEnvZ = envZ;

        if (nx > 0) fDx = envX / nx;
        if (ny > 0) fDy = envY / ny;
        if (nz > 0) fDz = envZ / nz;
    }

    void PkaRecorder::OpenEventFile(const G4String& baseName)
    {
        if (fEventOut.is_open())
            fEventOut.close();

        // 根据线程 ID 生成文件名
        std::ostringstream ss;
        G4int tid = G4Threading::G4GetThreadId();  // worker 一般是 0,1,2...，master 可能是 -1

        if (G4Threading::IsWorkerThread()) {
            ss << baseName << "_T" << tid << ".dat";
        }
        else {
            ss << baseName << "_master.dat"; // 万一 master 用到了，也有个名
        }

        G4String filename = ss.str();

        fEventOut.open(filename, std::ios::out);
        if (fEventOut) {
            fEventOut << "# eventID trackID  Z  A  Ek_eV   x_mm  y_mm  z_mm  "
                << "ux  uy  uz\n";
        }
    }

    void PkaRecorder::CloseEventFile()
    {
        if (fEventOut.is_open()) {
            fEventOut.flush();
            fEventOut.close();
        }
    }

    void PkaRecorder::RecordPka(const G4Track* track)
    {
        // 如果这个线程还没打开文件，就打开一个
        if (!fEventOut.is_open()) {
            OpenEventFile("pka");   // 每个线程会打开 pka_T<tid>.dat
        }

        // 下面和你原来的一样，只是去掉了锁
        const auto* pd = track->GetDefinition();
        G4int Z = pd->GetAtomicNumber();
        G4int A = pd->GetAtomicMass();

        G4double Ek = track->GetKineticEnergy();   // MeV
        G4ThreeVector pos = track->GetPosition();  // mm
        G4ThreeVector dir = track->GetMomentumDirection();

        auto runManager = G4RunManager::GetRunManager();
        G4int eventID = -1;
        if (runManager && runManager->GetCurrentEvent()) {
            eventID = runManager->GetCurrentEvent()->GetEventID();
        }

        G4int trackID = track->GetTrackID();
        G4double Ek_eV = Ek / eV;

        // ... 如果你还要算 ix,iy,iz，这里照旧 ...
        // 略去体素代码，只保留简单输出：

        fEventOut << eventID << " "
            << trackID << " "
            << Z << " " << A << " "
            << Ek_eV << " "
            << pos.x() / mm << " "
            << pos.y() / mm << " "
            << pos.z() / mm << " "
            << dir.x() << " "
            << dir.y() << " "
            << dir.z() << "\n";
    }

}