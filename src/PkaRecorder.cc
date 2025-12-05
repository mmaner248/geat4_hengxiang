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





    void PkaRecorder::RecordPka(const G4Track* track)
    {
        PkaHit hit;

        // 获取粒子信息
        const auto* pd = track->GetDefinition();
        hit.Z = pd->GetAtomicNumber();
        hit.A = pd->GetAtomicMass();

        G4double Ek = track->GetKineticEnergy();   // MeV
        auto runManager = G4RunManager::GetRunManager();
        if (runManager && runManager->GetCurrentEvent()) {
            hit.eventID = runManager->GetCurrentEvent()->GetEventID();
        }

        hit.trackID = track->GetTrackID();
        hit.Ek_eV = Ek / eV;


        
        hit.pos = track->GetPosition();            // mm
        hit.dir = track->GetMomentumDirection();   // 无量纲
        fHits.push_back(std::move(hit));
    }
    // 在 run 结束时调用，一次性写文件
    void PkaRecorder::WriteToFile(const G4String& baseName)
    {
        if (fHits.empty()) {
            return; // 本线程没记录到任何 PKA，直接跳过
        }

        // 根据线程 ID 生成文件名
        std::ostringstream ss;
        G4int tid = G4Threading::G4GetThreadId();

        if (G4Threading::IsWorkerThread()) {
            ss << baseName << "_T" << tid << ".dat";
        }
        else {
            ss << baseName << "_master.dat";
        }

        G4String filename = ss.str();

        std::ofstream out(filename, std::ios::out);
        if (!out) {
            G4Exception(
                "PkaRecorder::WriteToFile",
                "PKA_FILE_OPEN_FAIL",
                FatalException,
                ("Cannot open file " + filename).c_str()
            );
            return;
        }

        // header
        out << "# eventID trackID  Z  A  Ek_eV   x_mm  y_mm  z_mm  "
            << "ux  uy  uz   \n";

        // 写所有记录
        for (const auto& h : fHits) {
            out << h.eventID << " "
                << h.trackID << " "
                << h.Z << " " << h.A << " "
                << h.Ek_eV << " "
                << h.pos.x() / mm << " "
                << h.pos.y() / mm << " "
                << h.pos.z() / mm << " "
                << h.dir.x() << " "
                << h.dir.y() << " "
                << h.dir.z() << " "
                 << "\n";
        }

        out.close();

        // 如果一个 run 结束后不再用这些数据，可以清空释放内存
        fHits.clear();
        fHits.shrink_to_fit(); // 真的想省内存的话可以加，不加也行
    }
}