// PkaRecorder.cc
#include "PkaRecorder.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"   // Geant4 的自动加锁工具
#include <cmath>           // 用 std::floor 计算体素索引
namespace {
    // 用于保护文件写入的全局互斥锁（多线程用）
    G4Mutex pkaMutex = G4MUTEX_INITIALIZER;
}
namespace B1 {


    PkaRecorder* PkaRecorder::Instance()
    {
        static PkaRecorder instance;
        return &instance;
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

    void PkaRecorder::OpenEventFile(const G4String& filename)
    {
        if (fEventOut.is_open()) fEventOut.close();

        fEventOut.open(filename, std::ios::out);
        if (fEventOut) {
            fEventOut << "# eventID trackID  Z  A  Ek_eV   x_mm  y_mm  z_mm  "
                "ux  uy  uz   ix  iy  iz\n";
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
        // 如果 master 没有打开文件，这里直接返回，不做任何事
        if (!fEventOut.is_open()) return;

        // 先从 track 和 runManager 把所有需要的物理量都取出来
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

        // 能量转成 eV，方便后处理（MD 一般习惯 eV）
        G4double Ek_eV = Ek / eV;

        // 计算体素索引 ix, iy, iz
        G4int ix = -1, iy = -1, iz = -1;

        if (fNx > 0 && fNy > 0 && fNz > 0 &&
            fDx > 0. && fDy > 0. && fDz > 0.) {

            // 假设体素盒子中心在 (0,0,0)，范围是 [-fEnvX/2, +fEnvX/2] 等
            G4double x = pos.x();  // mm
            G4double y = pos.y();
            G4double z = pos.z();

            // 平移到 [0, fEnvX] 区间，再除以 cell 尺寸
            G4double xLocal = x + 0.5 * fEnvX;
            G4double yLocal = y + 0.5 * fEnvY;
            G4double zLocal = z + 0.5 * fEnvZ;

            ix = static_cast<G4int>(std::floor(xLocal / fDx));
            iy = static_cast<G4int>(std::floor(yLocal / fDy));
            iz = static_cast<G4int>(std::floor(zLocal / fDz));

            // 越界保护：如果不在 [0, nx) 里，就标成 -1
            if (ix < 0 || ix >= fNx ||
                iy < 0 || iy >= fNy ||
                iz < 0 || iz >= fNz) {
                ix = iy = iz = -1;
            }
        }

        // ★ 下面开始是对共享 ofstream 的写操作，要加锁保护（多线程）
        {
            G4AutoLock lock(&pkaMutex);   // 作用域锁：离开这个花括号自动解锁

            fEventOut << eventID << " "
                << trackID << " "
                << Z << " " << A << " "
                << Ek_eV << " "
                << pos.x() / mm << " "
                << pos.y() / mm << " "
                << pos.z() / mm << " "
                << dir.x() << " "
                << dir.y() << " "
                << dir.z() << " "
                << ix << " "
                << iy << " "
                << iz << "\n";
        }
    }

}