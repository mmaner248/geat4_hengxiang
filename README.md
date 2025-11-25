 Geant4 模拟过程中产生的 **PKA（Primary Knock-on Atom，主击出原子 / 重反冲原子）** 事件，用于后续分子动力学分析。

- 文件名：`pka_events.txt`
- 第一行为表头（列名）：

  ```text
  # eventID trackID  Z  A  Ek_eV   x_mm  y_mm  z_mm  ux  uy  uz   ix  iy  iz
 
 在 Geant4 中，每产生一条新轨迹（主粒子或次级粒子），都会调用一次用户的 TrackingAction
 // 1）跳过入射主粒子（primary）
if (track->GetParentID() == 0)
    return;

// 2）只考虑有原子序的重粒子
Z = track->GetDefinition()->GetAtomicNumber();
if (Z <= 0)
    return; // 排除中子、电子、光子等

// 3）动能阈值（当前设为 100 eV）
Ek = track->GetKineticEnergy(); // Geant4 内部单位
if (Ek < 100 * eV)
    return;

// 4）通过以上筛选 → 视为 PKA，记录
PkaRecorder::Instance()->RecordPka(track);

表头是啥意思：Z是原子序数，A是质量数，EK是能量，x,y,z是真实坐标，ux,uy,uz是输出方向向量的分量，ix,iy,iz是网格序号，当出现负值就表示不在网格（石墨立方体）里面了

